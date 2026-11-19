// REML of the kinship correlation h of a logistic model with a sparse kinship, on the IRLS working scale,
// marginal form: logit(mu) = eta = X beta, and the Pearson residuals (y - mu) / sqrt(mu (1 - mu)) have
// covariance R(h) = (1 - h) I + h K. At the current eta: w = mu (1 - mu), S = W^-1/2, working response
// yt = eta + (y - mu) / w with covariance V = S R(h) S = diag((1 - h) / w) + h S K S. With w and yt held fixed,
// h minimises
//   -2 l_R(h) = log|V| + log|X' V^-1 X| + yt' P yt
// by one Brent search on [0, h_max]; then beta is updated by generalised least squares, eta = X beta. The two
// steps alternate until eta and h settle. A sample without a relative has V_ii = (1 - h + h k_ii)
// / w_i, which is 1 / w_i whatever h when k_ii = 1: such samples carry no information on h, and at h = 0 the fit
// is ordinary logistic regression.
// The related samples' block of V is factorised by sparse Cholesky at every h (kin_sparse.h).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"

// y: n values 0 / 1; X: n x p full column rank (with the intercept); kinship as reml_tau_gaussian_cpp.
// h_max: upper end of the search; tol: Brent tolerance and the change of h at convergence; eta_tol: largest
// change of eta at convergence; max_iter: outer iterations (each of the two phases); mu_eps: fitted
// probabilities are kept in [mu_eps, 1 - mu_eps].
// Returns h, the -2 working REML log-likelihood at h and at 0 (same working data), beta, mu, and for the
// association stream (V = diag(d) + tau Ks): tau = h, d = (1 - h) / w, the scaled kinship ks = k_ij /
// sqrt(w_i w_j) and ksdiag = k_ii / w_i, r = P yt, Q, Vb.
// [[Rcpp::export]]
Rcpp::List reml_tau_binomial_cpp(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& X,
                               const Rcpp::IntegerVector& ki, const Rcpp::IntegerVector& kj,
                               const Rcpp::NumericVector& kv, const Rcpp::NumericVector& kdiag, double h_max,
                               double tol, double eta_tol, int max_iter, int max_eval, double mu_eps) {
  const int n = y.size(), p = X.ncol();
  const std::size_t nk = ki.size();
  if (X.nrow() != n || kdiag.size() != n) Rcpp::stop("y / X / kdiag do not match");
  if (kj.size() != static_cast<R_xlen_t>(nk) || kv.size() != static_cast<R_xlen_t>(nk))
    Rcpp::stop("kinship triplets do not match");
  if (n <= p) Rcpp::stop("more covariates than samples");
  if (!(h_max > 0 && h_max < 1)) Rcpp::stop("h_max must lie in (0, 1)");
  double ybar = 0;
  for (int i = 0; i < n; ++i) ybar += y[i];
  ybar /= n;
  if (!(ybar > 0 && ybar < 1)) Rcpp::stop("the phenotype has one class only");
  const double big = 1e300;
  std::vector<double> eta(n, std::log(ybar / (1 - ybar))), mu(n), iw(n), d(n), yt(n), ks(nk), ksd(n);
  plink4r::Kin K;
  plink4r::KinFactor F;
  // working data and the scaled kinship S K S at the current eta
  auto working = [&]() {
    for (int i = 0; i < n; ++i) {
      double m = 1.0 / (1.0 + std::exp(-eta[i]));
      m = std::min(1.0 - mu_eps, std::max(mu_eps, m));
      mu[i] = m;
      iw[i] = 1.0 / (m * (1.0 - m));
      yt[i] = eta[i] + (y[i] - m) * iw[i];
      ksd[i] = kdiag[i] * iw[i];
    }
    for (std::size_t e = 0; e < nk; ++e) ks[e] = kv[e] * std::sqrt(iw[ki[e]] * iw[kj[e]]);
    K = plink4r::build_kin(n, ki.begin(), kj.begin(), ks.data(), nk, ksd.data());
  };
  auto factor = [&](double h) {
    for (int i = 0; i < n; ++i) d[i] = (1.0 - h) * iw[i];
    return F.factor(K, d.data(), h);
  };
  auto m2l = [&](double h) {
    if (!factor(h)) return big;
    const plink4r::GlsFit g = plink4r::gls_fit(F, X.begin(), n, p, yt.data(), false);
    if (!g.ok) return big;
    return F.logdet + g.logdet_xox + g.quad;
  };
  // eta <- X beta at the given h; returns the largest change of eta (-1 on failure)
  auto update = [&](double h) {
    if (!factor(h)) return -1.0;
    const plink4r::GlsFit g = plink4r::gls_fit(F, X.begin(), n, p, yt.data(), false);
    if (!g.ok) return -1.0;
    double ch = 0;
    for (int i = 0; i < n; ++i) {
      double e = 0;
      for (int c = 0; c < p; ++c) e += X(i, c) * g.beta[c];
      ch = std::max(ch, std::fabs(e - eta[i]));
      eta[i] = e;
    }
    return ch;
  };

  try {
    // phase 1: h = 0 (logistic regression)
    bool conv0 = false;
    for (int it = 0; it < max_iter; ++it) {
      working();
      const double ch = update(0.0);
      if (ch < 0) Rcpp::stop("the logistic null model cannot be fitted (collinear covariates or separation?)");
      if (ch < eta_tol) { conv0 = true; break; }
    }
    // phase 2: alternate h and beta
    double h = 0;
    int iters = 0, evals = 0;
    bool conv = false;
    for (int it = 0; it < max_iter; ++it) {
      ++iters;
      working();
      const double f0 = m2l(0.0);
      plink4r::BrentResult br = plink4r::brent_min(m2l, 0.0, h_max, tol, tol * 1e-2, max_eval);
      evals += br.evals + 1;
      const double h_new = br.fx < f0 ? br.x : 0.0;
      const double ch = update(h_new);
      if (ch < 0) Rcpp::stop("the working model cannot be fitted");
      const bool done = ch < eta_tol && std::fabs(h_new - h) <= tol;
      h = h_new;
      if (done) { conv = true; break; }
    }
    working();
    const double f0 = m2l(0.0), fh = h > 0 ? m2l(h) : f0;
    if (!factor(h)) Rcpp::stop("V is not positive definite at the estimate");
    const plink4r::GlsFit g = plink4r::gls_fit(F, X.begin(), n, p, yt.data(), true);
    if (!g.ok) Rcpp::stop("generalised least squares failed at the estimate");

    Rcpp::NumericMatrix Q(n, p), Vb(n, p);
    std::copy(g.Q.begin(), g.Q.end(), Q.begin());
    std::copy(g.Vb.begin(), g.Vb.end(), Vb.begin());
    return Rcpp::List::create(
        Rcpp::Named("tau") = h, Rcpp::Named("m2l") = fh, Rcpp::Named("m2l_0") = f0,
        Rcpp::Named("iterations") = iters, Rcpp::Named("evals") = evals,
        Rcpp::Named("converged") = conv && conv0, Rcpp::Named("beta") = Rcpp::wrap(g.beta),
        Rcpp::Named("mu") = Rcpp::wrap(mu), Rcpp::Named("d") = Rcpp::wrap(d), Rcpp::Named("r") = Rcpp::wrap(g.r),
        Rcpp::Named("Q") = Q, Rcpp::Named("Vb") = Vb, Rcpp::Named("ks") = Rcpp::wrap(ks),
        Rcpp::Named("ksdiag") = Rcpp::wrap(ksd), Rcpp::Named("n_related") = K.nr,
        Rcpp::Named("pairs") = static_cast<double>(K.pairs));
  } catch (const std::runtime_error& e) { Rcpp::stop(e.what()); }
}
