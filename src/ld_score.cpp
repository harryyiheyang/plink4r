// LD scores of one eigen LD panel, optionally stratified by annotations: for every SNP j of a block and every
// annotation column c, l_jc = sum_k R_jk^2 A_kc over the SNPs k of the same block, R the stored (truncated)
// eigen factor of the block (ld_score_kernel.h). The all-ones column gives the usual LD score. OpenMP over
// blocks, largest first.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"
#include "ld_score_kernel.h"
#include "resources.h"

// files: the eigen.bin of each block; offset: B + 1 row offsets, block i owning rows offset[i]..offset[i+1]-1 of
// A (its SNPs by eigen row, Index 0..m-1); A: annotations (rows x C). Returns rows x C LD scores.
// [[Rcpp::export]]
Rcpp::NumericMatrix ld_score_cpp(Rcpp::CharacterVector files, Rcpp::IntegerVector offset,
                                 const Rcpp::NumericMatrix& A, double thresh, int chunk, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const int nb = files.size(), n = A.nrow(), C = A.ncol();
  if (offset.size() != nb + 1 || offset[0] != 0 || offset[nb] != n) Rcpp::stop("offset does not match A");
  const std::vector<std::string> fv = Rcpp::as<std::vector<std::string>>(files);
  const std::vector<int> off(offset.begin(), offset.end());
  const double* Ap = A.begin();
  std::vector<double> out(static_cast<std::size_t>(n) * C, 0.0);
  std::vector<std::string> err(nb);
  std::vector<int> order(nb);
  for (int i = 0; i < nb; ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return off[a + 1] - off[a] > off[b + 1] - off[b]; });
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
  for (int o = 0; o < nb; ++o) {
    const int i = order[o];
    try {
      const plink4r::EigenBlock ld = plink4r::read_eigen_block(fv[i], thresh);
      const int m = off[i + 1] - off[i];
      if (ld.m != m) throw std::runtime_error("the block has " + std::to_string(ld.m) + " eigen rows, the panel table " +
                                              std::to_string(m));
      std::vector<int> all(m);
      for (int j = 0; j < m; ++j) all[j] = j;
      const std::vector<double> plus(m, 1.0);
      std::vector<double> res(m);
      for (int c = 0; c < C; ++c) {
        std::vector<int> rows;
        std::vector<double> a;
        for (int j = 0; j < m; ++j) {
          const double v = Ap[static_cast<std::size_t>(c) * n + off[i] + j];
          if (v != 0) { rows.push_back(j); a.push_back(v); }
        }
        if (rows.empty()) continue;
        plink4r::weighted_ld_score(ld.U, ld.lambda, ld.U, ld.lambda, rows, rows, a, all, all, plus, chunk, res.data());
        std::copy(res.begin(), res.end(), out.begin() + static_cast<std::size_t>(c) * n + off[i]);
      }
    } catch (const std::exception& e) {
      err[i] = e.what();
    }
  }
  for (int i = 0; i < nb; ++i)
    if (!err[i].empty()) Rcpp::stop("%s: %s", fv[i].c_str(), err[i].c_str());
  Rcpp::NumericMatrix L(n, C);
  std::copy(out.begin(), out.end(), L.begin());
  return L;
}
