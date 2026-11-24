// Stratified LD score regression on the squared-effect scale, ported from ldscR::ldsc.strat (Yihe Yang):
//   E[Z_j^2 / N_j] = a / N_j + sum_c tau_c L_jc,
// Z truncated to [-z_max, z_max], regressors (1 / N_j, L_j1, ..., L_jC), L_j1 the overall LD score. The
// coefficient a (the LDSC intercept) is estimated inside [lo, hi] or fixed at 1. Weights
// w_j = 1 / max(1, L_j1) / (1 / N_j + h2 L_j1 / M)^2, first with h2 = w_init, then n_reweight times at the current
// h2 = sum_c tau_c sum_j A_jc (A the annotations of the regression SNPs), then fixed. With ridge, the tau_c get
// the penalty lambda * u * sum_c v_c tau_c^2 (v_c the variance of annotation c over the SNPs, u the median of
// (X'WX)_cc / v_c), lambda given or chosen by leave-one-block-out cross-validation over lambda_grid.
// Standard errors: leave-one-block-out jackknife over the LD blocks, the bound on a applied in every refit.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "ldsc_kernel.h"

namespace {
double median(std::vector<double> v) {
  const std::size_t n = v.size();
  std::sort(v.begin(), v.end());
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
}  // namespace

// z, N: M SNPs; L, A: M x C LD scores and annotations (column 1 of L the overall LD score, column 1 of A the
// all-ones column); block: 0-based LD block of every SNP (0..B-1); lambda NaN: chosen by cross-validation.
// [[Rcpp::export]]
Rcpp::List ldsc_strat_cpp(const Eigen::Map<Eigen::VectorXd> z, const Eigen::Map<Eigen::VectorXd> N,
                          const Eigen::Map<Eigen::MatrixXd> L, const Eigen::Map<Eigen::MatrixXd> A,
                          const Rcpp::IntegerVector& block, bool intercept, double lo, double hi, double z_max,
                          double w_init, int n_reweight, bool ridge, double lambda,
                          const Rcpp::NumericVector& lambda_grid, double tol) {
  using namespace plink4r::ldsc;
  const int M = static_cast<int>(z.size()), C = static_cast<int>(L.cols()), K = C + 1;
  if (N.size() != M || L.rows() != M || A.rows() != M || block.size() != M) Rcpp::stop("inputs differ in length");
  if (A.cols() != C) Rcpp::stop("L and A differ in their number of columns");
  if (M <= K) Rcpp::stop("fewer SNPs than regressors");
  std::vector<int> grp(block.begin(), block.end());
  const int B = *std::max_element(grp.begin(), grp.end()) + 1;
  if (*std::min_element(grp.begin(), grp.end()) < 0) Rcpp::stop("block must be 0-based");

  const Eigen::VectorXd acov = A.colwise().sum().transpose();
  const Eigen::VectorXd amean = A.colwise().mean().transpose();
  const Eigen::VectorXd avar = (A.array().square().colwise().mean().transpose() - amean.array().square()).max(0.0).matrix();
  Eigen::MatrixXd X(M, K);
  Eigen::VectorXd y(M), lw(M), invN(M);
  int n_trunc = 0;
  for (int j = 0; j < M; ++j) {
    const double az = std::abs(z[j]);
    if (az > z_max) ++n_trunc;
    const double zt = std::min(az, z_max);
    invN[j] = 1.0 / N[j];
    y[j] = zt * zt * invN[j];
    lw[j] = 1.0 / std::max(1.0, L(j, 0));
    X(j, 0) = invN[j];
  }
  X.rightCols(C) = L;
  Intercept ic;
  ic.mode = intercept ? 2 : 1;
  ic.fixed = 1.0; ic.lo = lo; ic.hi = hi;

  try {
    auto weights = [&](double h2) {
      return Eigen::VectorXd(lw.array() / (invN.array() + h2 * L.col(0).array() / M).square());
    };
    Eigen::VectorXd w = weights(w_init);
    Eigen::VectorXd b = fit_all(X, y, w, ic, tol);
    double h2 = b.tail(C).dot(acov);
    for (int it = 0; it < n_reweight; ++it) {
      w = weights(h2);
      b = fit_all(X, y, w, ic, tol);
      h2 = b.tail(C).dot(acov);
    }

    const Sums s = block_sums(X, y, w, grp, B);
    Eigen::MatrixXd XtX_all;
    Eigen::VectorXd Xty_all;
    const std::vector<double> ones(B, 1.0);
    combine(s, ones.data(), XtX_all, Xty_all);
    Eigen::MatrixXd Pen = Eigen::MatrixXd::Zero(K, K);
    double lambda_used = NA_REAL;
    Rcpp::NumericVector cv;
    if (ridge) {
      Eigen::VectorXd pen_d(K);
      pen_d[0] = 0;
      pen_d.tail(C) = avar;
      std::vector<double> ratio;
      for (int c = 0; c < K; ++c) if (pen_d[c] > 0) ratio.push_back(XtX_all(c, c) / pen_d[c]);
      if (ratio.empty()) Rcpp::stop("ridge: no annotation column varies over the SNPs");
      const double unit = median(ratio);
      if (std::isnan(lambda)) {
        const int G = lambda_grid.size();
        if (G < 1) Rcpp::stop("lambda_grid is empty");
        cv = Rcpp::NumericVector(G);
        for (int g = 0; g < G; ++g) {
          const Eigen::MatrixXd Pg = (lambda_grid[g] * unit * pen_d).asDiagonal();
          double e = 0;
          for (int bb = 0; bb < B; ++bb) {
            const Eigen::VectorXd f = fit_intercept(XtX_all + Pg - s.xtx[bb], Xty_all - s.xty[bb], ic, tol);
            e += -2.0 * f.dot(s.xty[bb]) + f.dot(s.xtx[bb] * f);
          }
          cv[g] = e;
        }
        const int best = static_cast<int>(std::min_element(cv.begin(), cv.end()) - cv.begin());
        lambda_used = lambda_grid[best];
        cv = cv - cv[best];
      } else {
        lambda_used = lambda;
      }
      Pen = (lambda_used * unit * pen_d).asDiagonal();
      b = fit_intercept(XtX_all + Pen, Xty_all, ic, tol);
      h2 = b.tail(C).dot(acov);
    }

    const Eigen::MatrixXd est = leave_one_out(s, Pen, ic, tol);
    const Eigen::MatrixXd h2_vec = est.rightCols(C) * acov;
    Rcpp::NumericVector se = col_jk(est);
    const double h2_se = col_jk(h2_vec)[0];
    if (!intercept) se[0] = 0;

    return Rcpp::List::create(
        Rcpp::Named("intercept") = b[0], Rcpp::Named("intercept.se") = se[0],
        Rcpp::Named("tau") = Rcpp::NumericVector(b.data() + 1, b.data() + K),
        Rcpp::Named("tau.se") = Rcpp::NumericVector(se.begin() + 1, se.end()),
        Rcpp::Named("h2") = h2, Rcpp::Named("h2.se") = h2_se, Rcpp::Named("annot_cov") = Rcpp::wrap(acov),
        Rcpp::Named("annot_var") = Rcpp::wrap(avar), Rcpp::Named("lambda") = lambda_used, Rcpp::Named("cv") = cv,
        Rcpp::Named("n_truncated") = n_trunc, Rcpp::Named("B") = B);
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
}
