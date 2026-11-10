// Polygenic scoring kernels behind geno_score().
//
// score_stream_cpp: plink2 --score style streaming. Every variant carries a 4-entry lookup table per trait
//   (value added to the score for PGEN code 0/1/2 = ALT count, 3 = missing) and a 4-entry allele-count table.
//   R builds the tables (A1 alignment, weights, transforms, the fill of a missing call: plink2 mean, scheme 3
//   anchor, none, or the constant part of an LD-filled cell), so the kernel holds no policy. Variants are read
//   in chunks (parallel over variants), then scored in parallel over sample ranges, so no reduction and no
//   genotype matrix: memory is n x q scores plus one chunk of genovecs.
// score_ld_block_cpp: one eigen.bin LD block through impute_kernel.h (the same stage 1 / stage 2 as
//   genotype_impute), adding B[j, t] * g for every filled cell instead of returning the cells.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cstring>
#include <vector>
#include "impute_kernel.h"

namespace {
// Byte decoding: a PGEN byte holds 4 samples (2 bits each, low bits first).
struct ByteTables {
  float x[256][4];        // ALT count, missing -> 0
  unsigned char miss[256];  // 4-bit missing mask
  ByteTables() {
    for (int b = 0; b < 256; ++b) {
      miss[b] = 0;
      for (int k = 0; k < 4; ++k) {
        const int c = (b >> (2 * k)) & 3;
        x[b][k] = c == 3 ? 0.0f : static_cast<float>(c);
        if (c == 3) miss[b] |= static_cast<unsigned char>(1u << k);
      }
    }
  }
};
const ByteTables& byte_tables() { static const ByteTables T; return T; }
}  // namespace

// var_idx: 0-based PGEN variants; sample_idx: 0-based PSAM rows (row order of the output);
// tab: 4 * q values per variant, layout [variant][code][trait]; ct_tab: 4 allele counts per variant;
// gemm: 0 = per-variant lookup (one trait), 2 = float32 product per chunk (several traits), 1 = double product.
// Each chunk's genovecs are first compacted to the selected samples (a straight copy when all samples are
// taken in file order) and then decoded a byte (4 samples) at a time.
// [[Rcpp::export]]
Rcpp::List score_stream_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct,
                            Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx, Rcpp::NumericVector tab,
                            Rcpp::IntegerVector ct_tab, int q, int threads, int chunk, int gemm) {
  const int T = std::max(1, threads);
  const int p = var_idx.size(), ns = sample_idx.size();
  if (tab.size() != static_cast<R_xlen_t>(4) * q * p || ct_tab.size() != static_cast<R_xlen_t>(4) * p)
    Rcpp::stop("score tables do not match the variants");
  const std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  const std::vector<int> act(allele_ct.begin(), allele_ct.end());
  if (!gemm && q != 1) Rcpp::stop("the lookup path scores one trait");
  bool ident = ns == n_raw;
  for (int s = 0; ident && s < ns; ++s) ident = si[s] == s;
  const double* tb = tab.begin();
  const int* ctb = ct_tab.begin();
  const ByteTables& BT = byte_tables();
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  const std::size_t sb = (static_cast<std::size_t>(ns) + 3) / 4;  // bytes per compacted variant
  const int V = std::max(1, chunk);
  std::vector<unsigned char> buf(static_cast<std::size_t>(V) * sb);
  std::vector<double> S(static_cast<std::size_t>(ns) * q, 0.0);  // row-major: sample, trait
  std::vector<int> ct(ns, 0);
  std::vector<int> bad(T, -1);
  const int sblk = gemm ? 1024 : 4096;  // multiple of 4: a sample block starts on a byte
  const int nsb = (ns + sblk - 1) / sblk;
  std::vector<double> base_sum(q, 0.0);
  if (gemm) for (int j = 0; j < p; ++j) for (int t = 0; t < q; ++t) base_sum[t] += tb[static_cast<std::size_t>(j) * 4 * q + t];
  std::vector<double> L1;  // lookup path: per chunk variant, 256 x 4 values (q = 1)
  if (!gemm) L1.resize(static_cast<std::size_t>(V) * 256 * 4);
  for (int c0 = 0; c0 < p; c0 += V) {
    const int nv = std::min(V, p - c0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(T)
#endif
    for (int v = 0; v < nv; ++v) {
      const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[c0 + v]), plink4r::ik::tid());
      if (!gv) { bad[plink4r::ik::tid()] = c0 + v; continue; }
      unsigned char* dst = &buf[static_cast<std::size_t>(v) * sb];
      if (ident) { std::memcpy(dst, gv, sb); continue; }
      std::memset(dst, 0, sb);
      for (int s = 0; s < ns; ++s) {
        const int r = si[s];
        dst[s >> 2] |= static_cast<unsigned char>(((gv[r >> 2] >> (2 * (r & 3))) & 3u) << (2 * (s & 3)));
      }
    }
    for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop("PGEN read failed at variant %d", bad[t] + 1);
    if (gemm) {
      Eigen::MatrixXd SlopeD(nv, q), MissD(nv, q);
      std::vector<int> ct_obs(nv), ct_mis(nv);
      for (int v = 0; v < nv; ++v) {
        const double* tv = tb + static_cast<std::size_t>(c0 + v) * 4 * q;
        for (int t = 0; t < q; ++t) { SlopeD(v, t) = tv[q + t] - tv[t]; MissD(v, t) = tv[3 * q + t] - tv[t]; }
        ct_obs[v] = ctb[static_cast<std::size_t>(c0 + v) * 4];
        ct_mis[v] = ctb[static_cast<std::size_t>(c0 + v) * 4 + 3];
      }
      const Eigen::MatrixXf SlopeF = SlopeD.cast<float>();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
      for (int b = 0; b < nsb; ++b) {
        const int s0 = b * sblk, s1 = std::min(ns, s0 + sblk), nr = s1 - s0;
        const int nr4 = (nr + 3) / 4 * 4;
        Eigen::MatrixXf XF(nr4, nv);  // padded to whole bytes; padding rows are dropped below
        std::vector<int> mi, mv;
        for (int v = 0; v < nv; ++v) {
          const unsigned char* g = &buf[static_cast<std::size_t>(v) * sb + (s0 >> 2)];
          float* col = XF.col(v).data();
          for (int y = 0; y < nr4 / 4; ++y) {
            const unsigned byte = g[y];
            std::memcpy(col + 4 * y, BT.x[byte], 4 * sizeof(float));
            unsigned mm = BT.miss[byte];
            while (mm) {
              const int k = __builtin_ctz(mm); mm &= mm - 1;
              const int s = s0 + 4 * y + k;
              if (s < s1) { mi.push_back(s); mv.push_back(v); }
            }
          }
        }
        Eigen::MatrixXd P;
        if (gemm == 2) P = (XF.topRows(nr) * SlopeF).cast<double>();
        else P = XF.topRows(nr).cast<double>() * SlopeD;
        int ctsum = 0;
        for (int v = 0; v < nv; ++v) ctsum += ct_obs[v];
        for (int s = s0; s < s1; ++s) {
          double* Ss = &S[static_cast<std::size_t>(s) * q];
          for (int t = 0; t < q; ++t) Ss[t] += P(s - s0, t);
          ct[s] += ctsum;
        }
        for (std::size_t e = 0; e < mi.size(); ++e) {
          double* Ss = &S[static_cast<std::size_t>(mi[e]) * q];
          for (int t = 0; t < q; ++t) Ss[t] += MissD(mv[e], t);
          ct[mi[e]] += ct_mis[mv[e]] - ct_obs[mv[e]];
        }
      }
      continue;
    }
    // lookup path (q = 1): per variant a 256 x 4 table byte -> the values of its 4 samples
    for (int v = 0; v < nv; ++v) {
      const double* tv = tb + static_cast<std::size_t>(c0 + v) * 4;
      double* Lv = &L1[static_cast<std::size_t>(v) * 1024];
      for (int byte = 0; byte < 256; ++byte)
        for (int k = 0; k < 4; ++k) Lv[4 * byte + k] = tv[(byte >> (2 * k)) & 3];
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
    for (int b = 0; b < nsb; ++b) {
      const int s0 = b * sblk, s1 = std::min(ns, s0 + sblk);
      const int full = s0 + (s1 - s0) / 4 * 4;  // samples covered by whole bytes
      for (int v = 0; v < nv; ++v) {
        const unsigned char* g = &buf[static_cast<std::size_t>(v) * sb];
        const double* Lv = &L1[static_cast<std::size_t>(v) * 1024];
        const int* cv = ctb + static_cast<std::size_t>(c0 + v) * 4;
        double* Sp = S.data();
        for (int s = s0; s < full; s += 4) {
          const double* e = Lv + 4 * g[s >> 2];
          Sp[s] += e[0]; Sp[s + 1] += e[1]; Sp[s + 2] += e[2]; Sp[s + 3] += e[3];
        }
        for (int s = full; s < s1; ++s) Sp[s] += Lv[4 * g[s >> 2] + (s & 3)];
        if (cv[3] == cv[0]) {
          for (int s = s0; s < s1; ++s) ct[s] += cv[0];
        } else {
          for (int s = s0; s < s1; ++s) ct[s] += ((g[s >> 2] >> (2 * (s & 3))) & 3u) == 3u ? cv[3] : cv[0];
        }
      }
    }
  }
  if (gemm) for (int s = 0; s < ns; ++s) for (int t = 0; t < q; ++t) S[static_cast<std::size_t>(s) * q + t] += base_sum[t];
  Rcpp::NumericMatrix out(ns, q);
  for (int s = 0; s < ns; ++s) for (int t = 0; t < q; ++t) out(s, t) = S[static_cast<std::size_t>(s) * q + t];
  return Rcpp::List::create(Rcpp::_["score"] = out, Rcpp::_["allele_ct"] = Rcpp::wrap(ct));
}

// One LD block (whole block, or the eigen rows `rows`): returns the n x q matrix sum_j B[j, t] * g_ij over the
// filled cells (g = final ALT dosage of a missing call). var_idx / af per working row as in genotype_impute_cpp;
// rows of B that are all zero are imputed (they condition the others) but add nothing.
// [[Rcpp::export]]
Rcpp::List score_ld_block_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct,
                              Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                              const std::string& eigen_file, Rcpp::NumericVector af, double thresh, double delta,
                              double anchor_min_miss, Rcpp::NumericMatrix B, int threads, int nb_people,
                              double mem_gb, Rcpp::IntegerVector rows) {
  const int T = std::max(1, threads);
  const plink4r::EigenBlock eb = plink4r::read_eigen_block(eigen_file, thresh);
  const std::vector<int> rw(rows.begin(), rows.end());
  const int m = rw.empty() ? eb.m : static_cast<int>(rw.size());
  if (var_idx.size() != m || af.size() != m || B.nrow() != m) Rcpp::stop("var_idx/af/B rows != eigen m");
  const int q = B.ncol();
  const std::vector<int> vi(var_idx.begin(), var_idx.end());
  const std::vector<double> afv(af.begin(), af.end());
  std::vector<char> keep(m, 0);
  for (int j = 0; j < m; ++j) for (int t = 0; t < q; ++t) if (B(j, t) != 0) { keep[j] = 1; break; }
  const plink4r::ImputeLD L = plink4r::build_impute_ld(eb, rw, vi, afv, delta, 0, 0, T);
  const std::vector<int> sidx(sample_idx.begin(), sample_idx.end());
  const std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  plink4r::BlockGenotypes G = plink4r::read_block_genotypes(pf, vi, sidx, T);
  const int ns = G.ns;
  const plink4r::PersonPlanes P = plink4r::build_planes(L, G, T);
  const plink4r::Stage1Result R = plink4r::impute_stage1(L, P, delta, nb_people, mem_gb, T);
  const plink4r::Stage2Result A = plink4r::impute_stage2(R, vi, G.sum_obs, afv, ns, anchor_min_miss);
  Rcpp::NumericMatrix S(ns, q);
  double cells = 0;
  plink4r::for_each_cell(R, A, keep, [&](int i, int j, float, double g) {
    for (int t = 0; t < q; ++t) S(i, t) += B(j, t) * g;
    cells += 1;
  });
  return Rcpp::List::create(Rcpp::_["score"] = S, Rcpp::_["cells"] = cells,
                            Rcpp::_["persons_solved"] = static_cast<double>(R.persons_solved));
}

// Instruction set the LD imputation kernel selected at run time (avx512, avx2 or scalar).
// [[Rcpp::export]]
std::string simd_isa_cpp() { return plink4r::ik::plane_dot_impl().isa; }
