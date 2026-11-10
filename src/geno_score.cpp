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

// var_idx: 0-based PGEN variants; sample_idx: 0-based PSAM rows (row order of the output);
// tab: 4 * q values per variant, layout [variant][code][trait]; ct_tab: 4 allele counts per variant.
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
  const double* tb = tab.begin();
  const int* ctb = ct_tab.begin();
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  const std::size_t rb = (static_cast<std::size_t>(n_raw) + 3) / 4;
  const int V = std::max(1, chunk);
  std::vector<unsigned char> buf(static_cast<std::size_t>(V) * rb);
  std::vector<double> S(static_cast<std::size_t>(ns) * q, 0.0);  // row-major: sample, trait
  std::vector<int> ct(ns, 0);
  std::vector<int> bad(T, -1);
  // gemm: codes 0/1/2 add base + x * slope (affine in the ALT count x), so a chunk is scored as one dense
  // product X (samples x variants, x with missing = 0) * Slope (variants x traits), in float (gemm = 2) or
  // double (gemm = 1), plus the missing-cell corrections tab[3] - tab[0] and the constant sum of base.
  const int sblk = gemm ? 1024 : 2048;
  const int nsb = (ns + sblk - 1) / sblk;
  std::vector<double> base_sum(q, 0.0);
  if (gemm) for (int j = 0; j < p; ++j) for (int t = 0; t < q; ++t) base_sum[t] += tb[static_cast<std::size_t>(j) * 4 * q + t];
  for (int c0 = 0; c0 < p; c0 += V) {
    const int nv = std::min(V, p - c0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(T)
#endif
    for (int v = 0; v < nv; ++v) {
      const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[c0 + v]), plink4r::ik::tid());
      if (!gv) { bad[plink4r::ik::tid()] = c0 + v; continue; }
      std::memcpy(&buf[static_cast<std::size_t>(v) * rb], gv, rb);
    }
    for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop("PGEN read failed at variant %d", bad[t] + 1);
    if (gemm) {
      Eigen::MatrixXd SlopeD(nv, q), MissD(nv, q);
      for (int v = 0; v < nv; ++v) {
        const double* tv = tb + static_cast<std::size_t>(c0 + v) * 4 * q;
        for (int t = 0; t < q; ++t) { SlopeD(v, t) = tv[q + t] - tv[t]; MissD(v, t) = tv[3 * q + t] - tv[t]; }
      }
      const Eigen::MatrixXf SlopeF = SlopeD.cast<float>();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
      for (int b = 0; b < nsb; ++b) {
        const int s0 = b * sblk, s1 = std::min(ns, s0 + sblk), nr = s1 - s0;
        Eigen::MatrixXd XD;
        Eigen::MatrixXf XF;
        if (gemm == 2) XF.resize(nr, nv); else XD.resize(nr, nv);
        std::vector<int> mi, mv;
        for (int v = 0; v < nv; ++v) {
          const unsigned char* g = &buf[static_cast<std::size_t>(v) * rb];
          const int* cv = ctb + static_cast<std::size_t>(c0 + v) * 4;
          for (int s = s0; s < s1; ++s) {
            const int r = si[s];
            const unsigned c = (g[r >> 2] >> (2 * (r & 3))) & 3u;
            const unsigned x = c == 3u ? 0u : c;
            if (gemm == 2) XF(s - s0, v) = static_cast<float>(x); else XD(s - s0, v) = static_cast<double>(x);
            ct[s] += cv[c];
            if (c == 3u) { mi.push_back(s); mv.push_back(v); }
          }
        }
        Eigen::MatrixXd P = gemm == 2 ? Eigen::MatrixXd((XF * SlopeF).cast<double>()) : Eigen::MatrixXd(XD * SlopeD);
        for (int s = s0; s < s1; ++s) {
          double* Ss = &S[static_cast<std::size_t>(s) * q];
          for (int t = 0; t < q; ++t) Ss[t] += P(s - s0, t);
        }
        for (std::size_t e = 0; e < mi.size(); ++e) {
          double* Ss = &S[static_cast<std::size_t>(mi[e]) * q];
          for (int t = 0; t < q; ++t) Ss[t] += MissD(mv[e], t);
        }
      }
      continue;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
    for (int b = 0; b < nsb; ++b) {
      const int s0 = b * sblk, s1 = std::min(ns, s0 + sblk);
      for (int v = 0; v < nv; ++v) {
        const unsigned char* g = &buf[static_cast<std::size_t>(v) * rb];
        const double* tv = tb + static_cast<std::size_t>(c0 + v) * 4 * q;
        const int* cv = ctb + static_cast<std::size_t>(c0 + v) * 4;
        if (q == 1) {
          for (int s = s0; s < s1; ++s) {
            const int r = si[s];
            const unsigned c = (g[r >> 2] >> (2 * (r & 3))) & 3u;
            S[s] += tv[c];
            ct[s] += cv[c];
          }
        } else {
          for (int s = s0; s < s1; ++s) {
            const int r = si[s];
            const unsigned c = (g[r >> 2] >> (2 * (r & 3))) & 3u;
            const double* tc = tv + c * q;
            double* Ss = &S[static_cast<std::size_t>(s) * q];
            for (int t = 0; t < q; ++t) Ss[t] += tc[t];
            ct[s] += cv[c];
          }
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
