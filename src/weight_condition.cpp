// Conditional weights: the score weights of panel SNPs that the target lacks are predicted onto the target SNPs
// that have no weight; weights on SNPs in both are kept as they are.
// The whole computation is this one entry (R only reads the tables and writes the result):
//   1. quad_match (quad_match.h): target -> panel and weights -> panel (each on the panel build that places more
//      rows), and the weights outside the panel -> target directly;
//   2. per LD block with weights to move (a weighted panel row absent from the target): O = the block rows in the
//      target without a weight (the receivers), M = every other row of the block; B_M = the weights of the rows
//      absent from the target, Bk = the kept weights (rows in both), each in panel orientation on the standardized
//      scale (times sd = sqrt(2 p (1 - p)), p = panel A1Freq); cond_weight_block (la_kernel.cpp) gives
//      Delta_O = (R_OO + delta I)^{-1} R_OM B_M and the score variance lost (with and without the move) and total;
//   3. Delta back to the allele scale of the target (/ sd, times the target's sign).
// Blocks run longest first, OpenMP over blocks; the dense algebra uses the instruction set chosen at run time.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"
#include "la_dispatch.h"
#include "quad_match.h"
#include "resources.h"

namespace {

const char* cw_method_name(int code) {
  switch (code) { case 1: return "missing"; case 2: return "eigen"; case 3: return "typed"; default: return "none"; }
}

void cw_header(const std::string& path, int& m, int& k) {
  m = 0; k = 0;
  std::ifstream in(path.c_str(), std::ios::binary);
  int32_t h[2] = {0, 0};
  if (in.read(reinterpret_cast<char*>(h), sizeof(h))) { m = h[0]; k = h[1]; }
}

std::vector<int> iv(const Rcpp::IntegerVector& x) { return std::vector<int>(x.begin(), x.end()); }
std::vector<std::string> sv(const Rcpp::CharacterVector& x) { return Rcpp::as<std::vector<std::string>>(x); }

struct BlockJob {
  int block = 0;
  std::string file;
  std::vector<int> prow;  // panel rows of the block
  double cost = 0;
};

}  // namespace

// Panel columns p_* (CHR code, BP37, BP38 with NA where absent, A1, A2, A1Freq, Block, Index); blocks/files: the
// panel's blocks and their eigen.bin; target t_* (CHR code, BP, A1, A2); weights w_* and W (one column per score).
// Returns list(w = beta + delta, beta = the weights on the target rows in target orientation (n_target x q; 0 where none),
// delta = the moved part (n_target x q; 0 on rows with a weight), target_panel / target_weight (logical), weight_status (per weight row:
// "target" = on a target row, "moved" = moved through the LD, "no_target_in_block" = in a block without target
// rows that lack a weight, "lost" = neither in the panel nor in the target), build_target, build_weight (0 = BP37, 1 = BP38),
// block_log (block, m, k, n_obs, n_miss, n_miss_weighted, method, seconds), kept, kept_drop (blocks x q:
// 1 - lost / total with and without the move, NA when total = 0), isa).
// [[Rcpp::export]]
Rcpp::List weight_condition_cpp(Rcpp::IntegerVector p_chr, Rcpp::IntegerVector p_bp37, Rcpp::IntegerVector p_bp38,
                                Rcpp::CharacterVector p_a1, Rcpp::CharacterVector p_a2, Rcpp::NumericVector p_freq,
                                Rcpp::IntegerVector p_block, Rcpp::IntegerVector p_index, Rcpp::IntegerVector blocks,
                                Rcpp::CharacterVector files, Rcpp::IntegerVector t_chr, Rcpp::IntegerVector t_bp,
                                Rcpp::CharacterVector t_a1, Rcpp::CharacterVector t_a2, Rcpp::IntegerVector w_chr,
                                Rcpp::IntegerVector w_bp, Rcpp::CharacterVector w_a1, Rcpp::CharacterVector w_a2,
                                Rcpp::NumericMatrix W, double thresh, double delta, int slice, int threads) {
  const int T = plink4r::resolve_threads(threads);
  if (!(delta > 0)) Rcpp::stop("delta must be > 0");
  const int np = p_chr.size(), nt = t_chr.size(), nw = w_chr.size(), q = W.ncol();
  if (W.nrow() != nw) Rcpp::stop("W must have one row per weight");
  const int na = NA_INTEGER;
  // plain copies: read inside the OpenMP loop
  const std::vector<double> Wv(W.begin(), W.end()), pf(p_freq.begin(), p_freq.end());
  const std::vector<int> pidx = iv(p_index);
  auto w_at = [&](int row, int c) { return Wv[static_cast<std::size_t>(c) * nw + row]; };
  // 1. matching
  const std::vector<int> pc = iv(p_chr), p37 = iv(p_bp37), p38 = iv(p_bp38), tc = iv(t_chr), tb = iv(t_bp), wc = iv(w_chr),
                         wb = iv(w_bp);
  const std::vector<std::string> pa1 = sv(p_a1), pa2 = sv(p_a2), ta1 = sv(t_a1), ta2 = sv(t_a2), wa1 = sv(w_a1),
                                 wa2 = sv(w_a2);
  const plink4r::Quad P0{&pc, &p37, &pa1, &pa2}, P1{&pc, &p38, &pa1, &pa2}, Tq{&tc, &tb, &ta1, &ta2},
      Wq{&wc, &wb, &wa1, &wa2};
  int bt = 0, bw = 0;
  const plink4r::QuadHit th = plink4r::quad_match_best(P0, P1, Tq, na, bt);
  const plink4r::QuadHit wh = plink4r::quad_match_best(P0, P1, Wq, na, bw);
  for (int j = 0; j < nt; ++j)
    if (th.row[j] == -2) Rcpp::stop("target row %d repeats a variant (CHR, BP, A1, A2) of an earlier row", j + 1);
  for (int j = 0; j < nw; ++j)
    if (wh.row[j] == -2) Rcpp::stop("weight row %d repeats a variant (CHR, BP, A1, A2) of an earlier row", j + 1);
  // weights outside the panel -> target directly (same build as the weights)
  std::vector<int> wc_out, wb_out;
  std::vector<std::string> wa1_out, wa2_out;
  std::vector<int> out_idx;
  for (int j = 0; j < nw; ++j)
    if (wh.row[j] < 0) { out_idx.push_back(j); wc_out.push_back(wc[j]); wb_out.push_back(wb[j]); wa1_out.push_back(wa1[j]); wa2_out.push_back(wa2[j]); }
  const plink4r::Quad Wo{&wc_out, &wb_out, &wa1_out, &wa2_out};
  const plink4r::QuadHit dh = plink4r::quad_match(Tq, Wo, na);
  for (std::size_t j = 0; j < out_idx.size(); ++j)
    if (dh.row[j] == -2) Rcpp::stop("weight row %d repeats a variant (CHR, BP, A1, A2) of an earlier row", out_idx[j] + 1);
  // panel row -> target row (and sign), panel row -> weight row (and sign)
  std::vector<int> p2t(np, -1), p2w(np, -1), p2ts(np, 0), p2ws(np, 0);
  for (int j = 0; j < nt; ++j) if (th.row[j] >= 0) { p2t[th.row[j]] = j; p2ts[th.row[j]] = th.sign[j]; }
  for (int j = 0; j < nw; ++j) if (wh.row[j] >= 0) { p2w[wh.row[j]] = j; p2ws[wh.row[j]] = wh.sign[j]; }
  Rcpp::NumericMatrix beta(nt, q), dlt(nt, q);
  Rcpp::LogicalVector t_in_panel(nt), t_in_w(nt);
  std::vector<std::string> wstat(nw, "lost");
  for (int j = 0; j < nt; ++j) t_in_panel[j] = th.row[j] >= 0;
  for (int r = 0; r < np; ++r)
    if (p2t[r] >= 0 && p2w[r] >= 0) {
      const int t = p2t[r], w = p2w[r];
      const double s = static_cast<double>(p2ts[r]) * p2ws[r];  // weight A1 vs target A1
      for (int c = 0; c < q; ++c) beta(t, c) = s * W(w, c);
      t_in_w[t] = true;
      wstat[w] = "target";
    }
  for (std::size_t j = 0; j < out_idx.size(); ++j)
    if (dh.row[j] >= 0) {
      const int t = dh.row[j], w = out_idx[j];
      for (int c = 0; c < q; ++c) beta(t, c) = dh.sign[j] * W(w, c);
      t_in_w[t] = true;
      wstat[w] = "target";
    }
  // 2. blocks with weights to move
  std::map<int, std::string> bfile;
  for (int i = 0; i < blocks.size(); ++i) bfile[blocks[i]] = Rcpp::as<std::string>(files[i]);
  std::map<int, std::vector<int>> rows_of;
  for (int r = 0; r < np; ++r) rows_of[p_block[r]].push_back(r);
  std::vector<BlockJob> jobs;
  for (auto& kv : rows_of) {
    bool move = false;
    for (int r : kv.second) if (p2w[r] >= 0 && p2t[r] < 0) { move = true; break; }
    if (!move) continue;
    const auto f = bfile.find(kv.first);
    if (f == bfile.end()) Rcpp::stop("no eigen.bin given for block %d", kv.first);
    BlockJob b;
    b.block = kv.first; b.file = f->second; b.prow = kv.second;
    int m, k;
    cw_header(b.file, m, k);
    const double M = m, K = k;
    b.cost = M * K * K;
    jobs.push_back(std::move(b));
  }
  const int nb = static_cast<int>(jobs.size());
  std::vector<int> order(nb);
  for (int i = 0; i < nb; ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return jobs[a].cost > jobs[b].cost; });
  const plink4r::LaFns& LA = plink4r::la();
  std::vector<std::string> err(nb), meth(nb);
  std::vector<int> mm(nb), kk(nb), nobs(nb), nmiss(nb), nmw(nb);
  std::vector<double> sec(nb);
  std::vector<std::vector<double>> lost(nb, std::vector<double>(q)), lost0(nb, std::vector<double>(q)),
      total(nb, std::vector<double>(q));
  std::vector<std::vector<std::pair<int, std::vector<double>>>> add(nb);  // (target row, delta on its allele scale)
  std::vector<char> has_obs(nb, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
  for (int o = 0; o < nb; ++o) {
    const int i = order[o];
    const BlockJob& b = jobs[i];
    const auto t0 = std::chrono::steady_clock::now();
    try {
      const plink4r::EigenBlock ld = plink4r::read_eigen_block(b.file, thresh);
      mm[i] = ld.m; kk[i] = ld.k;
      std::vector<int> row_at(ld.m, -1);  // eigen row -> panel row
      for (int r : b.prow) {
        const int x = pidx[r];
        if (x < 0 || x >= ld.m) throw std::runtime_error("panel Index outside the block");
        row_at[x] = r;
      }
      std::vector<int> obs, miss;
      for (int x = 0; x < ld.m; ++x) (row_at[x] >= 0 && p2t[row_at[x]] >= 0 && p2w[row_at[x]] < 0 ? obs : miss).push_back(x);
      const int no = static_cast<int>(obs.size()), nm = static_cast<int>(miss.size());
      auto sd_of = [&](int r) { const double p = pf[r]; return (std::isfinite(p) && p > 0 && p < 1) ? std::sqrt(2 * p * (1 - p)) : 0.0; };
      std::vector<double> B(static_cast<std::size_t>(nm) * q, 0.0), Bk(static_cast<std::size_t>(nm) * q, 0.0);
      int nw_m = 0;
      for (int a = 0; a < nm; ++a) {
        const int r = row_at[miss[a]];
        if (r < 0 || p2w[r] < 0) continue;
        const bool kept = p2t[r] >= 0;
        if (!kept) ++nw_m;
        std::vector<double>& dst = kept ? Bk : B;
        const double f = p2ws[r] * sd_of(r);
        for (int c = 0; c < q; ++c) dst[static_cast<std::size_t>(c) * nm + a] = f * w_at(p2w[r], c);
      }
      nobs[i] = no; nmiss[i] = nm; nmw[i] = nw_m; has_obs[i] = no > 0;
      std::vector<double> D(static_cast<std::size_t>(no) * q, 0.0);
      char e[256] = {0};
      const int code = LA.cond_weight_block(ld.U.data(), ld.m, ld.k, ld.lambda.data(), obs.data(), no, miss.data(), nm,
                                            B.data(), Bk.data(), q, delta, slice, D.data(), lost[i].data(), lost0[i].data(),
                                            total[i].data(), e);
      if (code < 0) throw std::runtime_error(e);
      meth[i] = cw_method_name(code);
      for (int a = 0; a < no; ++a) {
        const int r = row_at[obs[a]];
        const double sd = sd_of(r);
        if (sd <= 0) continue;
        std::vector<double> v(q);
        for (int c = 0; c < q; ++c) v[c] = p2ts[r] * D[static_cast<std::size_t>(c) * no + a] / sd;
        add[i].emplace_back(p2t[r], std::move(v));
      }
    } catch (const std::exception& ex) {
      err[i] = ex.what();
    }
    sec[i] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  for (int i = 0; i < nb; ++i)
    if (!err[i].empty()) Rcpp::stop("block %d: %s", jobs[i].block, err[i].c_str());
  // 3. delta on the target rows; status of the moved weights
  for (int i = 0; i < nb; ++i) {
    for (const auto& tv : add[i])
      for (int c = 0; c < q; ++c) dlt(tv.first, c) += tv.second[c];
    for (int r : jobs[i].prow)
      if (p2w[r] >= 0 && p2t[r] < 0) wstat[p2w[r]] = has_obs[i] ? "moved" : "no_target_in_block";
  }
  Rcpp::NumericMatrix wout(nt, q);
  for (int c = 0; c < q; ++c)
    for (int j = 0; j < nt; ++j) wout(j, c) = beta(j, c) + dlt(j, c);
  Rcpp::IntegerVector lb(nb);
  Rcpp::NumericMatrix KP(nb, q), KD(nb, q);
  for (int i = 0; i < nb; ++i) {
    lb[i] = jobs[i].block;
    for (int c = 0; c < q; ++c) {
      KP(i, c) = total[i][c] > 0 ? 1 - lost[i][c] / total[i][c] : NA_REAL;
      KD(i, c) = total[i][c] > 0 ? 1 - lost0[i][c] / total[i][c] : NA_REAL;
    }
  }
  Rcpp::DataFrame log = Rcpp::DataFrame::create(
      Rcpp::_["block"] = lb, Rcpp::_["m"] = mm, Rcpp::_["k"] = kk, Rcpp::_["n_obs"] = nobs, Rcpp::_["n_miss"] = nmiss,
      Rcpp::_["n_miss_weighted"] = nmw, Rcpp::_["method"] = meth, Rcpp::_["seconds"] = sec);
  return Rcpp::List::create(Rcpp::_["w"] = wout, Rcpp::_["beta"] = beta, Rcpp::_["delta"] = dlt, Rcpp::_["target_panel"] = t_in_panel,
                            Rcpp::_["target_weight"] = t_in_w, Rcpp::_["weight_status"] = Rcpp::wrap(wstat),
                            Rcpp::_["build_target"] = bt, Rcpp::_["build_weight"] = bw, Rcpp::_["block_log"] = log,
                            Rcpp::_["kept"] = KP, Rcpp::_["kept_drop"] = KD, Rcpp::_["isa"] = LA.isa);
}
