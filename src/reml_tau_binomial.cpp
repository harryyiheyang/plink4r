// The kinship correlation h of a logistic model with a sparse kinship, marginal form: logit(mu) = eta = X beta
// and Var(y) = W^1/2 H(h) W^1/2, W = diag(mu (1 - mu)), H(h) = (1 - h) I + h K: the Pearson residuals
// e = (y - mu) / sqrt(mu (1 - mu)) have covariance H(h). A sample without a relative has variance mu (1 - mu)
// whatever h (k_ii = 1), so it carries no information on h, and at h = 0 the fit is ordinary logistic regression.
//
// Nested optimisation of one function of h:
//   inner  beta_hat(h): X' W^1/2 H^-1 e = 0, by Fisher scoring to convergence. One step: at the current eta the
//          working response yt = eta + (y - mu) / w has covariance W^-1/2 H W^-1/2; generalised least squares of
//          yt on X gives beta, eta = X beta.
//   outer  f(h) = log|H(h)| + sum log w + log|X' W^1/2 H^-1 W^1/2 X| + e' H^-1 e, all at beta_hat(h): the Gaussian
//          restricted quasi-likelihood (-2 log) of y with these two moments. h minimises f; f is a function of h
//          alone, so the search is Brent's, without the derivative of beta_hat with respect to h.
// Start: h = 0 (logistic regression), then start_iter steps of the alternating scheme (h by Brent on the working
// likelihood with w and yt held fixed, then one scoring step), which costs little and lands near the minimum;
// f is then minimised on a window of half-width h_window around that value (on the whole interval if the
// minimum falls on the edge of the window). h is exactly 0 when f(0) is not improved by more than f_tol (f is
// known to the accuracy of the scoring), and when no pair is related.
// Per scoring step the scaled data S X, S yt (S = W^1/2) and their cross-products are formed once; a fit at h
// factorises the related samples' block of H and solves with their rows only (kin_sparse.h).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"

// y: n values 0 / 1; X: n x p (with the intercept); kinship and rank_tol as reml_tau_gaussian_cpp.
// h_max: upper end of the search (or the largest h with H(h) positive definite, h_bound); tol: Brent tolerance;
// eta_tol: largest change of eta at convergence of the scoring; max_iter: scoring steps per value of h;
// max_eval: evaluations per Brent search; mu_eps: fitted probabilities are kept in [mu_eps, 1 - mu_eps];
// start_iter, h_window, f_tol: as above. proj: also the null model for the association stream.
// Returns h (= tau), f at h and at 0, beta, the scoring steps (iterations) and evaluations of f (evals) made,
// keep (1-based columns of X used) and null (external pointer to the plink4r::KinNull; V^-1 = W^1/2 H^-1 W^1/2).
// [[Rcpp::export]]
Rcpp::List reml_tau_binomial_cpp(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& X,
                                 const Rcpp::IntegerVector& ki, const Rcpp::IntegerVector& kj,
                                 const Rcpp::NumericVector& kv, double kdiag, double h_max, double tol,
                                 double eta_tol, int max_iter, int max_eval, double mu_eps, int start_iter,
                                 double h_window, double f_tol, double rank_tol, bool proj) {
  const int n = y.size();
  int p = X.ncol();
  const std::size_t N = n;
  if (X.nrow() != n) Rcpp::stop("y / X do not match");
  if (ki.size() != kj.size() || ki.size() != kv.size()) Rcpp::stop("kinship triplets do not match");
  if (!(h_max > 0 && h_max < 1)) Rcpp::stop("h_max must lie in (0, 1)");
  if (!(h_window > 0)) Rcpp::stop("h_window must be > 0");
  double ybar = 0;
  for (int i = 0; i < n; ++i) ybar += y[i];
  ybar /= n;
  if (!(ybar > 0 && ybar < 1)) Rcpp::stop("the phenotype has one class only");
  Rcpp::XPtr<plink4r::KinNull> null(new plink4r::KinNull, true);
  plink4r::KinNull& M = *null;
  try {
    M.init(plink4r::build_kin(n, ki.begin(), kj.begin(), kv.begin(), ki.size(), kdiag));
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }

  const double big = 1e300;
  const double h_top = M.h_bound(h_max, tol);  // H(h) is positive definite on [0, h_top]
  const double* Xp = X.begin();
  std::vector<double> Xc, eta(N, std::log(ybar / (1 - ybar))), sw(N), yt(N), ys(N), Xs(N * p);
  double slw = 0;  // sum log w of the working data
  // working data at the current eta, scaled by S
  auto working = [&]() {
    slw = 0;
    for (int i = 0; i < n; ++i) {
      double m = 1.0 / (1.0 + std::exp(-eta[i]));
      m = std::min(1.0 - mu_eps, std::max(mu_eps, m));
      const double w = m * (1.0 - m);
      sw[i] = std::sqrt(w);
      slw += std::log(w);
      yt[i] = eta[i] + (y[i] - m) / w;
      ys[i] = sw[i] * yt[i];
    }
    for (int c = 0; c < p; ++c) {
      const double* x = Xp + c * N;
      double* xs = &Xs[c * N];
      for (std::size_t i = 0; i < N; ++i) xs[i] = sw[i] * x[i];
    }
    M.set_data(Xs.data(), ys.data(), p);
  };
  auto m2l = [&](double h) {
    if (!M.solve(h)) return big;
    return M.logdet + M.logdet_xox + M.quad;
  };
  // eta <- X beta at the given h; returns the largest change of eta (-1 on failure)
  auto update = [&](double h) {
    if (!M.solve(h)) return -1.0;
    Eigen::Map<const Eigen::MatrixXd> Xm(Xp, n, p);
    const Eigen::VectorXd e = Xm * M.beta;
    double ch = 0;
    for (int i = 0; i < n; ++i) {
      ch = std::max(ch, std::fabs(e[i] - eta[i]));
      eta[i] = e[i];
    }
    return ch;
  };

  // aliased columns: the weights of the first working data are constant, so its X block is X'X up to a factor
  working();
  const std::vector<int> keep = M.independent(rank_tol);
  if (static_cast<int>(keep.size()) < p) {
    p = static_cast<int>(keep.size());
    Xc.resize(N * p);
    for (int c = 0; c < p; ++c) std::copy(Xp + keep[c] * N, Xp + (keep[c] + 1) * N, &Xc[c * N]);
    Xp = Xc.data();
    working();
  }
  if (n <= p) Rcpp::stop("more covariates than samples");

  // f(h): scoring from the current eta to beta_hat(h); the fit, eta and working data are left at it
  int steps = 0, evals = 0;
  bool conv = true;
  auto f = [&](double h) {
    ++evals;
    for (int it = 0; it < max_iter; ++it) {
      if (steps > 0) working();
      ++steps;
      const double ch = update(h);
      if (ch < 0) return big;
      if (ch < eta_tol) return M.logdet + slw + M.logdet_xox + M.quad;
    }
    conv = false;
    return M.logdet + slw + M.logdet_xox + M.quad;
  };
  const double f0 = f(0.0);  // logistic regression
  if (f0 >= big || !conv) Rcpp::stop("the logistic null model cannot be fitted (separation?)");
  // start value: steps of the alternating scheme
  double h = 0;
  for (int it = 0; it < start_iter; ++it) {
    working();
    ++steps;
    const double w0 = m2l(0.0);
    const plink4r::BrentResult br = plink4r::brent_min(m2l, 0.0, h_top, tol, tol * 1e-2, max_eval);
    h = br.fx < w0 ? br.x : 0.0;
    if (update(h) < 0) Rcpp::stop("the working model cannot be fitted");
  }
  // f minimised around the start value; without a related pair H = I whatever h, and h = 0
  double fh = f0;
  if (M.K.nr > 0) {
    const double lo = std::max(0.0, h - h_window), hi = std::min(h_top, h + h_window);
    plink4r::BrentResult br = plink4r::brent_min(f, lo, hi, tol, tol * 1e-2, max_eval);
    const double edge = 0.01 * h_window;
    if ((lo > 0 && br.x - lo < edge) || (hi < h_top && hi - br.x < edge))
      br = plink4r::brent_min(f, 0.0, h_top, tol, tol * 1e-2, max_eval);
    h = f0 - br.fx > f_tol ? br.x : 0.0;
  } else {
    h = 0.0;
  }
  fh = f(h);  // the fit at h is the one left in M
  if (h == 0.0) fh = f0;
  if (fh >= big) Rcpp::stop("the model cannot be fitted at the estimate");
  if (proj) M.project(Xp, yt.data(), sw.data(), 1.0);

  Rcpp::IntegerVector kp(keep.begin(), keep.end());
  kp = kp + 1;
  return Rcpp::List::create(
      Rcpp::Named("tau") = h, Rcpp::Named("h_bound") = h_top, Rcpp::Named("m2l") = fh, Rcpp::Named("m2l_0") = f0,
      Rcpp::Named("iterations") = steps, Rcpp::Named("evals") = evals, Rcpp::Named("converged") = conv,
      Rcpp::Named("beta") = Rcpp::NumericVector(M.beta.data(), M.beta.data() + p), Rcpp::Named("keep") = kp,
      Rcpp::Named("n_related") = M.K.nr, Rcpp::Named("pairs") = static_cast<double>(M.K.pairs),
      Rcpp::Named("null") = null);
}
