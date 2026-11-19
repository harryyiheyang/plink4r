// REML of the variance components of a linear mixed model with a sparse kinship,
//   y = X beta + b + e,  b ~ N(0, s2g K),  e ~ N(0, s2e I),
// by one Brent search over h = s2g / (s2g + s2e): with H(h) = (1 - h) I + h K and V = s2 H, s2 is profiled out
// (s2 = y' P_H y / (n - p)) and
//   -2 l_R(h) = log|H| + log|X' H^-1 X| + (n - p) log(y' P_H y) + const.
// The related samples' block of H is factorised by sparse Cholesky at every h (kin_sparse.h); for everyone else
// H is a scalar. Returns the estimates and the null-model quantities the association stream needs.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"

// y: n; X: n x p full column rank (with the intercept); ki, kj, kv: off-diagonal kinship entries (0-based
// samples, each pair once); kdiag: n diagonal entries. h is searched on [0, h_max].
// Returns s2g (= tau), s2e, h, the -2 REML log-likelihood at h and at 0 (same constant), beta, and for the
// association stream d (n, = s2e), r = P y, Q, Vb (see kin_sparse.h).
// [[Rcpp::export]]
Rcpp::List reml_tau_gaussian_cpp(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& X,
                                 const Rcpp::IntegerVector& ki, const Rcpp::IntegerVector& kj,
                                 const Rcpp::NumericVector& kv, const Rcpp::NumericVector& kdiag, double h_max,
                                 double tol, int max_eval) {
  const int n = y.size(), p = X.ncol();
  if (X.nrow() != n || kdiag.size() != n) Rcpp::stop("y / X / kdiag do not match");
  if (ki.size() != kj.size() || ki.size() != kv.size()) Rcpp::stop("kinship triplets do not match");
  if (n <= p) Rcpp::stop("more covariates than samples");
  if (!(h_max > 0 && h_max < 1)) Rcpp::stop("h_max must lie in (0, 1)");
  plink4r::Kin K;
  try {
    K = plink4r::build_kin(n, ki.begin(), kj.begin(), kv.begin(), ki.size(), kdiag.begin());
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
  plink4r::KinFactor F;
  std::vector<double> d(n);
  const double big = 1e300;
  const double df = static_cast<double>(n - p);
  auto m2l = [&](double h) {
    std::fill(d.begin(), d.end(), 1.0 - h);
    if (!F.factor(K, d.data(), h)) return big;
    const plink4r::GlsFit g = plink4r::gls_fit(F, X.begin(), n, p, y.begin(), false);
    if (!g.ok || !(g.quad > 0)) return big;
    return F.logdet + g.logdet_xox + df * std::log(g.quad);
  };
  const double f0 = m2l(0.0);
  if (f0 >= big) Rcpp::stop("the null model without kinship cannot be fitted (collinear covariates?)");
  plink4r::BrentResult br = plink4r::brent_min(m2l, 0.0, h_max, tol, tol * 1e-2, max_eval);
  double h = br.x, fh = br.fx;
  if (!(fh < f0)) { h = 0.0; fh = f0; }

  std::fill(d.begin(), d.end(), 1.0 - h);
  if (!F.factor(K, d.data(), h)) Rcpp::stop("V is not positive definite at the estimate");
  plink4r::GlsFit g = plink4r::gls_fit(F, X.begin(), n, p, y.begin(), false);
  const double s2 = g.quad / df, s2g = h * s2, s2e = (1.0 - h) * s2;
  std::fill(d.begin(), d.end(), s2e);
  if (!F.factor(K, d.data(), s2g)) Rcpp::stop("V is not positive definite at the estimate");
  g = plink4r::gls_fit(F, X.begin(), n, p, y.begin(), true);
  if (!g.ok) Rcpp::stop("generalised least squares failed at the estimate");

  Rcpp::NumericMatrix Q(n, p), Vb(n, p);
  std::copy(g.Q.begin(), g.Q.end(), Q.begin());
  std::copy(g.Vb.begin(), g.Vb.end(), Vb.begin());
  return Rcpp::List::create(
      Rcpp::Named("tau") = s2g, Rcpp::Named("s2e") = s2e, Rcpp::Named("h") = h,
      Rcpp::Named("m2l") = fh, Rcpp::Named("m2l_0") = f0, Rcpp::Named("evals") = br.evals + 1,
      Rcpp::Named("beta") = Rcpp::wrap(g.beta), Rcpp::Named("d") = Rcpp::wrap(d),
      Rcpp::Named("r") = Rcpp::wrap(g.r), Rcpp::Named("Q") = Q, Rcpp::Named("Vb") = Vb,
      Rcpp::Named("n_related") = K.nr, Rcpp::Named("pairs") = static_cast<double>(K.pairs));
}
