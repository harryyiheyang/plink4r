// Cross-ancestry LD scores of two eigen LD panels: for a SNP j present in both, l_j = sum_k R1_jk R2_jk over the
// SNPs k present in both panels and in the same block of each (ld_score_kernel.h), R2 turned to count the alleles
// of panel 1. The SNPs are grouped by their (block of panel 1, block of panel 2) pair; OpenMP over groups,
// largest first.
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

// files1, files2: eigen.bin of the two blocks of each group; offset: G + 1 offsets into i1, i2, sign (the 0-based
// eigen rows of the group's SNPs in each block, and +1 / -1 when panel 2 counts the other allele).
// Returns the cross LD score of every SNP, in the order of i1.
// [[Rcpp::export]]
Rcpp::NumericVector ld_score_cross_cpp(Rcpp::CharacterVector files1, Rcpp::CharacterVector files2,
                                       Rcpp::IntegerVector offset, Rcpp::IntegerVector i1, Rcpp::IntegerVector i2,
                                       Rcpp::NumericVector sign, double thresh1, double thresh2, int chunk,
                                       int threads) {
  const int T = plink4r::resolve_threads(threads);
  const int ng = files1.size(), n = i1.size();
  if (files2.size() != ng || offset.size() != ng + 1 || offset[0] != 0 || offset[ng] != n || i2.size() != n ||
      sign.size() != n)
    Rcpp::stop("files, offset, i1, i2 and sign do not match");
  const std::vector<std::string> f1 = Rcpp::as<std::vector<std::string>>(files1),
                                 f2 = Rcpp::as<std::vector<std::string>>(files2);
  const std::vector<int> off(offset.begin(), offset.end()), r1(i1.begin(), i1.end()), r2(i2.begin(), i2.end());
  const std::vector<double> s(sign.begin(), sign.end());
  std::vector<double> out(n, 0.0);
  std::vector<std::string> err(ng);
  std::vector<int> order(ng);
  for (int g = 0; g < ng; ++g) order[g] = g;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return off[a + 1] - off[a] > off[b + 1] - off[b]; });
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
  for (int o = 0; o < ng; ++o) {
    const int g = order[o];
    try {
      const plink4r::EigenBlock b1 = plink4r::read_eigen_block(f1[g], thresh1);
      const plink4r::EigenBlock b2 = plink4r::read_eigen_block(f2[g], thresh2);
      const std::vector<int> a1(r1.begin() + off[g], r1.begin() + off[g + 1]), a2(r2.begin() + off[g], r2.begin() + off[g + 1]);
      const std::vector<double> sg(s.begin() + off[g], s.begin() + off[g + 1]);
      for (int j = 0; j < static_cast<int>(a1.size()); ++j)
        if (a1[j] < 0 || a1[j] >= b1.m || a2[j] < 0 || a2[j] >= b2.m) throw std::runtime_error("eigen row out of range");
      plink4r::weighted_ld_score(b1.U, b1.lambda, b2.U, b2.lambda, a1, a2, sg, a1, a2, sg, chunk, &out[off[g]]);
    } catch (const std::exception& e) {
      err[g] = e.what();
    }
  }
  for (int g = 0; g < ng; ++g)
    if (!err[g].empty()) Rcpp::stop("%s / %s: %s", f1[g].c_str(), f2[g].c_str(), err[g].c_str());
  return Rcpp::NumericVector(out.begin(), out.end());
}
