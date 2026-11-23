// Univariate LD score regression (SNP heritability), ported from ldscR::ldsc.univ (Yihe Yang).
// With l_j = L2_j N_j / M, Z_j^2 is regressed on (1, l_j) by weighted least squares: first with
// w_j = 1 / max(1, L2_j) / (1 + w_init l_j)^2, then n_reweight times with w_j = 1 / max(1, L2_j) / (1 + l_j h2)^2
// at the current h2. Standard errors: block bootstrap over contiguous groups of SNPs (ldsc_kernel.h), with the
// final weights; the group counts of every replicate come from R (`counts`, R x B).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include "ldsc_kernel.h"

// z, N, L2: M SNPs (in genomic order); counts: R x B bootstrap group counts (R = 0: no bootstrap), B groups from
// nblock. Returns intercept, h2, their bootstrap SEs (NA without bootstrap) and the R x 2 replicate estimates.
// [[Rcpp::export]]
Rcpp::List ldsc_univ_cpp(const Eigen::Map<Eigen::VectorXd> z, const Eigen::Map<Eigen::VectorXd> N,
                         const Eigen::Map<Eigen::VectorXd> L2, int nblock, const Rcpp::NumericMatrix& counts,
                         double w_init, int n_reweight, double tol) {
  using namespace plink4r::ldsc;
  const int M = static_cast<int>(z.size());
  if (N.size() != M || L2.size() != M) Rcpp::stop("z, N and L2 differ in length");
  if (M < 3) Rcpp::stop("fewer than 3 SNPs");
  Eigen::MatrixXd X(M, 2);
  Eigen::VectorXd y(M), w(M), lw(M);
  for (int j = 0; j < M; ++j) {
    X(j, 0) = 1.0;
    X(j, 1) = L2[j] * N[j] / M;
    y[j] = z[j] * z[j];
    lw[j] = 1.0 / std::max(1.0, L2[j]);
  }
  const Intercept ic;
  try {
    w = lw.cwiseQuotient((1.0 + w_init * X.col(1).array()).square().matrix());
    Eigen::VectorXd b = fit_all(X, y, w, ic, tol);
    for (int it = 0; it < n_reweight; ++it) {
      w = lw.cwiseQuotient((1.0 + b[1] * X.col(1).array()).square().matrix());
      b = fit_all(X, y, w, ic, tol);
    }
    Eigen::MatrixXd est(0, 2);
    Rcpp::NumericVector se(2, NA_REAL);
    if (counts.nrow() > 0) {
      std::vector<int> grp;
      const int B = equal_groups(M, nblock, grp);
      const Sums s = block_sums(X, y, w, grp, B);
      est = refit(s, counts, Eigen::MatrixXd::Zero(2, 2), ic, tol);
      se = col_sd(est);
    }
    return Rcpp::List::create(Rcpp::Named("intercept") = b[0], Rcpp::Named("intercept.se") = se[0],
                              Rcpp::Named("h2") = b[1], Rcpp::Named("h2.se") = se[1], Rcpp::Named("M") = M,
                              Rcpp::Named("boot") = Rcpp::wrap(est));
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}
