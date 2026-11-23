// Number of contiguous bootstrap groups LDSC forms from M SNPs and nblock (ldsc_kernel.h, equal_groups); R
// draws the bootstrap counts over that many groups before a fit.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "ldsc_kernel.h"

// [[Rcpp::export]]
int ldsc_groups_cpp(int M, int nblock) {
  std::vector<int> grp;
  try { return plink4r::ldsc::equal_groups(M, nblock, grp); } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}
