// REML of the variance components of a linear mixed model with a sparse kinship,
//   y = X beta + b + e,  b ~ N(0, s2g K),  e ~ N(0, s2e I),
// i.e. V = s2e I + s2g K, by one Brent search over h = s2g / (s2g + s2e): writing V = s2 H(h),
// H(h) = (1 - h) I + h K, s2 = s2g + s2e (the same V, s2e = (1 - h) s2, s2g = h s2), s2 is profiled out
// (s2 = y' P_H y / (n - p)) and
//   -2 l_R(h) = log|H| + log|X' H^-1 X| + (n - p) log(y' P_H y) + const.
// X'X, X'y, y'y are formed once; an evaluation at h factorises the related samples' block of H and solves with
// their rows only (kin_sparse.h). Returns the estimates and the null model for the association stream.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"

// y: n; X: n x p (with the intercept); ki, kj, kv: off-diagonal kinship entries (0-based samples, each pair
// once); kdiag: the diagonal of K. h is searched on [0, h_max]. Columns of X that are a linear combination of
// earlier ones (rank_tol) are left out. proj: also the null model for the association stream.
// Returns s2g (= tau), s2e, h, the -2 REML log-likelihood at h and at 0 (same constant), beta, keep (1-based
// columns of X used) and null (external pointer to the plink4r::KinNull, fitted at V = s2e I + s2g K).
// [[Rcpp::export]]
Rcpp::List reml_tau_gaussian_cpp(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& X,
                                 const Rcpp::IntegerVector& ki, const Rcpp::IntegerVector& kj,
                                 const Rcpp::NumericVector& kv, double kdiag, double h_max, double tol,
                                 int max_eval, double rank_tol, bool proj) {
  const int n = y.size();
  int p = X.ncol();
  if (X.nrow() != n) Rcpp::stop("y / X do not match");
  if (ki.size() != kj.size() || ki.size() != kv.size()) Rcpp::stop("kinship triplets do not match");
  if (!(h_max > 0 && h_max < 1)) Rcpp::stop("h_max must lie in (0, 1)");
  Rcpp::XPtr<plink4r::KinNull> null(new plink4r::KinNull, true);
  plink4r::KinNull& N = *null;
  try {
    N.init(plink4r::build_kin(n, ki.begin(), kj.begin(), kv.begin(), ki.size(), kdiag));
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
  const double* Xp = X.begin();
  N.set_data(Xp, y.begin(), p);
  N.cross();
  const std::vector<int> keep = N.independent(rank_tol);
  std::vector<double> Xc;
  if (static_cast<int>(keep.size()) < p) {
    p = static_cast<int>(keep.size());
    Xc.resize(static_cast<std::size_t>(n) * p);
    for (int c = 0; c < p; ++c) std::copy(Xp + static_cast<std::size_t>(keep[c]) * n,
                                          Xp + static_cast<std::size_t>(keep[c] + 1) * n, &Xc[static_cast<std::size_t>(c) * n]);
    Xp = Xc.data();
    N.set_data(Xp, y.begin(), p);
    N.cross();
  }
  if (n <= p) Rcpp::stop("more covariates than samples");
  const double big = 1e300;
  const double df = static_cast<double>(n - p);
  auto m2l = [&](double h) {
    if (!N.solve_const(h, 1.0 - h) || !(N.quad > 0)) return big;
    return N.logdet + N.logdet_xox + df * std::log(N.quad);
  };
  const double f0 = m2l(0.0);
  if (f0 >= big) Rcpp::stop("the null model without kinship cannot be fitted");
  plink4r::BrentResult br = plink4r::brent_min(m2l, 0.0, h_max, tol, tol * 1e-2, max_eval);
  double h = br.x, fh = br.fx;
  if (!(fh < f0)) { h = 0.0; fh = f0; }
  if (!N.solve_const(h, 1.0 - h)) Rcpp::stop("V is not positive definite at the estimate");
  const double s2 = N.quad / df;
  if (!N.solve_const(h * s2, (1.0 - h) * s2)) Rcpp::stop("V is not positive definite at the estimate");
  if (proj) N.project();

  Rcpp::IntegerVector kp(keep.begin(), keep.end());
  kp = kp + 1;
  return Rcpp::List::create(
      Rcpp::Named("tau") = h * s2, Rcpp::Named("s2e") = (1.0 - h) * s2, Rcpp::Named("h") = h,
      Rcpp::Named("m2l") = fh, Rcpp::Named("m2l_0") = f0, Rcpp::Named("evals") = br.evals + 1,
      Rcpp::Named("beta") = Rcpp::NumericVector(N.beta.data(), N.beta.data() + p), Rcpp::Named("keep") = kp,
      Rcpp::Named("n_related") = N.K.nr, Rcpp::Named("pairs") = static_cast<double>(N.K.pairs),
      Rcpp::Named("null") = null);
}
