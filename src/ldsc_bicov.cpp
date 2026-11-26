// Bivariate LD score regression (genetic covariance of two traits on one LD score), ported from
// ldscR::ldsc.bicov (Yihe Yang). With l_j = L2_j sqrt(N1_j N2_j) / M, z1_j z2_j is regressed on (1, l_j):
// first with w_j = 1 / max(1, L2_j) / v_j, v_j = (1 + L2_j N1_j h21 / M)(1 + L2_j N2_j h22 / M), then n_reweight
// times with w_j = 1 / max(1, L2_j) / (v_j + 2 mu_j^2), mu_j = ecov + gcov l_j at the current fit. Standard
// errors: block bootstrap over contiguous groups of SNPs with the final weights (counts from R).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include "ldsc_kernel.h"

// Returns ecov, gcov, their bootstrap SEs (NA without bootstrap) and the R x 2 replicate estimates.
// [[Rcpp::export]]
Rcpp::List ldsc_bicov_cpp(const Eigen::Map<Eigen::VectorXd> z1, const Eigen::Map<Eigen::VectorXd> z2,
                          const Eigen::Map<Eigen::VectorXd> N1, const Eigen::Map<Eigen::VectorXd> N2,
                          const Eigen::Map<Eigen::VectorXd> L2, double h21, double h22, int nblock,
                          const Rcpp::NumericMatrix& counts, int n_reweight, double tol) {
  using namespace plink4r::ldsc;
  const int M = static_cast<int>(z1.size());
  if (z2.size() != M || N1.size() != M || N2.size() != M || L2.size() != M) Rcpp::stop("inputs differ in length");
  if (M < 3) Rcpp::stop("fewer than 3 SNPs");
  Eigen::MatrixXd X(M, 2);
  Eigen::VectorXd y(M), v(M), lw(M), w(M);
  for (int j = 0; j < M; ++j) {
    X(j, 0) = 1.0;
    X(j, 1) = L2[j] * std::sqrt(N1[j]) * std::sqrt(N2[j]) / M;
    y[j] = z1[j] * z2[j];
    v[j] = (1.0 + L2[j] * N1[j] / M * h21) * (1.0 + L2[j] * N2[j] / M * h22);
    lw[j] = 1.0 / std::max(1.0, L2[j]);
  }
  const Intercept ic;
  try {
    w = lw.cwiseQuotient(v);
    Eigen::VectorXd b = fit_all(X, y, w, ic, tol);
    for (int it = 0; it < n_reweight; ++it) {
      const Eigen::ArrayXd mu = b[0] + b[1] * X.col(1).array();
      w = lw.array() / (v.array() + 2.0 * mu.square());
      b = fit_all(X, y, w, ic, tol);
    }
    Eigen::MatrixXd est(0, 2);
    Rcpp::NumericVector se(2, NA_REAL);
    if (counts.nrow() > 0) {
      std::vector<int> grp;
      const int B = equal_groups(M, nblock, grp);
      est = refit(block_sums(X, y, w, grp, B), counts, Eigen::MatrixXd::Zero(2, 2), ic, tol);
      se = col_sd(est);
    }
    return Rcpp::List::create(Rcpp::Named("ecov") = b[0], Rcpp::Named("ecov.se") = se[0],
                              Rcpp::Named("gcov") = b[1], Rcpp::Named("gcov.se") = se[1], Rcpp::Named("M") = M,
                              Rcpp::Named("boot") = Rcpp::wrap(est));
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}
