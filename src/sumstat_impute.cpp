// Summary-statistics imputation from eigen.bin, ported from SBayesRC.QC2 0.1.4 (Yihe Yang).
// Per block: z_M = R_MO (R_OO + delta I)^{-1} z_O with R = U diag(lam) U', solved in the smallest exact space
// (missing-SNP / Woodbury eigen-space / typed-SNP), float32, OpenMP over blocks; the per-block algebra lives in
// la_kernel.cpp (sumstat_block), compiled per instruction set.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
#include <string>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"
#include "la_dispatch.h"
#include "resources.h"

using plink4r::EigenBlock;

namespace {

std::vector<int> missing_indices(int m, const std::vector<int>& typed) {
  std::vector<int> miss;
  miss.reserve(m - typed.size());
  std::size_t j = 0;
  for (int i = 0; i < m; ++i) {
    if (j < typed.size() && typed[j] == i) ++j;
    else miss.push_back(i);
  }
  return miss;
}

// m and the stored k from the first 8 bytes of block<N>.eigen.bin (0, 0 if unreadable).
void block_header(const std::string& path, int& m, int& k) {
  m = 0; k = 0;
  std::ifstream in(path.c_str(), std::ios::binary);
  int32_t h[2] = {0, 0};
  if (in.read(reinterpret_cast<char*>(h), sizeof(h))) { m = h[0]; k = h[1]; }
}

// Flops of the chosen solver (smallest of missing / eigen / typed space), for longest-first scheduling.
double block_cost(int m, int k, int nt) {
  const double nm = std::max(0, m - nt), K = k, T = nt;
  if (nm == 0) return 0;
  if (nm <= T && nm <= K) return nm * K * K + nm * nm * K + nm * nm * nm + m * K;
  if (K < T) return T * K * K + K * K * K + nm * K * K;
  return T * T * K + T * T * T + T * T * nm;
}

const char* method_name(int code) {
  switch (code) { case 1: return "missing"; case 2: return "eigen"; case 3: return "typed"; default: return "none"; }
}

}  // namespace

// typed_index: list of 0-based strictly increasing eigen-row indices per block; z: matching z-scores.
// Returns, per block, list(z, r2) for the missing rows (all rows not in typed_index, increasing): z the imputed
// z-score and r2 = R_jO (R_OO + delta I)^{-1} R_Oj / R_jj its imputation quality; attr block_log.
// Blocks run longest first (cost from the eigen.bin headers), OpenMP over blocks; the dense algebra uses the
// instruction set chosen at run time (la_dispatch.h).
// [[Rcpp::export]]
Rcpp::List sumstat_impute_cpp(const std::string& ld_dir, Rcpp::IntegerVector blocks, Rcpp::List typed_index, Rcpp::List z,
                              double thresh, double delta, int threads) {
  const int T = plink4r::resolve_threads(threads);
  if (!(delta > 0)) Rcpp::stop("delta must be > 0");
  const int nb = blocks.size();
  if (typed_index.size() != nb || z.size() != nb) Rcpp::stop("blocks, typed_index and z differ in length");
  std::vector<int> bv(nb);
  std::vector<std::vector<int>> tv(nb);
  std::vector<std::vector<float>> zv(nb);
  for (int i = 0; i < nb; ++i) {
    bv[i] = blocks[i];
    tv[i] = Rcpp::as<std::vector<int>>(typed_index[i]);
    zv[i] = Rcpp::as<std::vector<float>>(z[i]);
  }
  std::vector<double> cost(nb);
  for (int i = 0; i < nb; ++i) {
    int m, k;
    block_header(plink4r::block_path(ld_dir, bv[i]), m, k);
    cost[i] = block_cost(m, k, static_cast<int>(tv[i].size()));
  }
  std::vector<int> order(nb);
  for (int i = 0; i < nb; ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return cost[a] > cost[b]; });
  const plink4r::LaFns& LA = plink4r::la();
  std::vector<std::vector<float>> res(nb), rq(nb);
  std::vector<std::string> err(nb), meth(nb);
  std::vector<int> mm(nb), kk(nb);
  std::vector<double> sec(nb);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
  for (int o = 0; o < nb; ++o) {
    const int i = order[o];
    const auto t0 = std::chrono::steady_clock::now();
    try {
      const EigenBlock ld = plink4r::read_eigen_block(plink4r::block_path(ld_dir, bv[i]), thresh);
      mm[i] = ld.m; kk[i] = ld.k;
      for (std::size_t j = 0; j < tv[i].size(); ++j)
        if (tv[i][j] < 0 || tv[i][j] >= ld.m || (j && tv[i][j] <= tv[i][j - 1])) throw std::runtime_error("bad typed_index");
      const std::vector<int> miss = missing_indices(ld.m, tv[i]);
      res[i].assign(miss.size(), 0.0f);
      rq[i].assign(miss.size(), 0.0f);
      char e[256] = {0};
      const int code = LA.sumstat_block(ld.U.data(), ld.m, ld.k, ld.lambda.data(), tv[i].data(), static_cast<int>(tv[i].size()),
                                        zv[i].data(), miss.data(), static_cast<int>(miss.size()), static_cast<float>(delta),
                                        res[i].data(), rq[i].data(), e);
      if (code < 0) throw std::runtime_error(e);
      meth[i] = method_name(code);
    } catch (const std::exception& e) {
      err[i] = e.what();
    }
    sec[i] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  for (int i = 0; i < nb; ++i)
    if (!err[i].empty()) Rcpp::stop("block %d: %s", bv[i], err[i].c_str());
  // per block, length m by eigen row: typed rows keep their z with r2 = 1, the others get the imputed z and r2
  Rcpp::List out(nb);
  for (int i = 0; i < nb; ++i) {
    Rcpp::NumericVector zf(mm[i]), qf(mm[i]);
    for (std::size_t j = 0; j < tv[i].size(); ++j) { zf[tv[i][j]] = zv[i][j]; qf[tv[i][j]] = 1; }
    const std::vector<int> miss = missing_indices(mm[i], tv[i]);
    for (std::size_t j = 0; j < miss.size(); ++j) { zf[miss[j]] = res[i][j]; qf[miss[j]] = rq[i][j]; }
    out[i] = Rcpp::List::create(Rcpp::_["z"] = zf, Rcpp::_["r2"] = qf);
  }
  Rcpp::IntegerVector ntyped(nb);
  for (int i = 0; i < nb; ++i) ntyped[i] = static_cast<int>(tv[i].size());
  out.attr("block_log") = Rcpp::DataFrame::create(Rcpp::_["block"] = bv, Rcpp::_["m"] = mm, Rcpp::_["k"] = kk,
                                                  Rcpp::_["n_typed"] = ntyped, Rcpp::_["method"] = meth,
                                                  Rcpp::_["seconds"] = sec, Rcpp::_["isa"] = LA.isa);
  return out;
}
