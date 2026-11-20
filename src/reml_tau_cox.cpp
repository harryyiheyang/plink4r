// The kinship correlation h of a Cox model with a sparse kinship, on the Poisson form of the Cox likelihood: with
// the cumulative baseline hazard Lambda0, the event indicator d_i of a sample followed to time t_i has mean
//   mu_i = Lambda0(t_i) exp(eta_i),  eta = X beta,
// and the Cox likelihood in beta at a fixed Lambda0 is that of a Poisson model for d with log Lambda0(t_i) as an
// offset. Marginal form, as the logistic model (reml_tau_binomial.cpp): Var(d) = W^1/2 H(h) W^1/2, W = diag(mu),
// H(h) = (1 - h) I + h K, i.e. the Pearson residuals e = (d - mu) / sqrt(mu) (the martingale residuals, scaled)
// have covariance H(h). At h = 0 the fit is that of the Cox partial likelihood (Breslow ties).
//
// Nested optimisation of one function of h:
//   inner  beta_hat(h) and Lambda0_hat(h), by scoring until f settles. One step: at the current eta, Lambda0 is
//          the Breslow estimate (jumps d_k / sum of exp(eta) over the risk set), which gives mu; the working
//          response yt = eta + (d - mu) / mu has covariance W^-1/2 H W^-1/2, and generalised least squares of yt
//          on [1 | X] gives beta, eta = X beta. The intercept only rescales the baseline (mu does not change with
//          it) and is left out of eta; at the solution X' W^1/2 H^-1 e = 0 with X centred by the intercept, which
//          at h = 0 is the partial-likelihood score equation.
//   outer  f(h) = log|H(h)| + sum log mu + log|X1' W^1/2 H^-1 W^1/2 X1| + e' H^-1 e (X1 = [1 | X]), all at
//          beta_hat(h), Lambda0_hat(h); h minimises f, by Brent's search.
// Start, window, convergence and h = 0 as in reml_tau_binomial.cpp. The association stream then tests every
// variant in the Poisson model with the offset log Lambda0_hat and the intercept (gwas_kin_stream.cpp).
// A sample that leaves before the first event has Lambda0 = 0 and carries no information; it takes the smallest
// positive value of Lambda0.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"

// time, status: n follow-up times and event indicators 0 / 1; X: n x p with the intercept as its first column;
// kinship and rank_tol as reml_tau_gaussian_cpp.
// h_max: upper end of the search (or the largest h with H(h) positive definite, h_bound); tol: Brent tolerance
// on h; max_iter: scoring steps per value of h; max_eval: evaluations per Brent search; mu_eps: lower bound of
// the expected events mu; start_iter, h_window, f_tol: as reml_tau_binomial_cpp. proj: also the null model for
// the association stream.
// Returns h (= tau), f at h and at 0, beta (the intercept first: the log of the factor by which the Poisson fit
// rescales the baseline), the baseline (base_time: distinct event times, increasing; base_cumhaz: Lambda0 there),
// the scoring steps (iterations) and evaluations of f (evals) made, keep (1-based columns of X used) and null
// (external pointer to the plink4r::KinNull; V^-1 = W^1/2 H^-1 W^1/2).
// [[Rcpp::export]]
Rcpp::List reml_tau_cox_cpp(const Rcpp::NumericVector& time, const Rcpp::NumericVector& status,
                            const Rcpp::NumericMatrix& X, const Rcpp::IntegerVector& ki,
                            const Rcpp::IntegerVector& kj, const Rcpp::NumericVector& kv, double kdiag,
                            double h_max, double tol, int max_iter, int max_eval, double mu_eps, int start_iter,
                            double h_window, double f_tol, double rank_tol, bool proj) {
  const int n = time.size();
  int p = X.ncol();
  const std::size_t N = n;
  if (X.nrow() != n || status.size() != n || p < 1) Rcpp::stop("time / status / X do not match");
  if (ki.size() != kj.size() || ki.size() != kv.size()) Rcpp::stop("kinship triplets do not match");
  if (!(h_max > 0 && h_max < 1)) Rcpp::stop("h_max must lie in (0, 1)");
  if (!(h_window > 0)) Rcpp::stop("h_window must be > 0");
  for (int i = 0; i < n; ++i)
    if (X(i, 0) != 1.0) Rcpp::stop("the first column of X must be the intercept");
  // samples in decreasing time; an event time is a run of equal times with at least one event
  std::vector<int> ord(n), gend;  // gend: last position of every event time
  std::vector<double> gd;         // its number of events
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return time[a] > time[b]; });
  for (int pos = 0; pos < n;) {
    int e = pos;
    double dk = 0;
    while (e < n && time[ord[e]] == time[ord[pos]]) dk += status[ord[e++]];
    if (dk > 0) { gend.push_back(e - 1); gd.push_back(dk); }
    pos = e;
  }
  const int ng = static_cast<int>(gend.size());
  if (ng < 1) Rcpp::stop("no event");
  Rcpp::XPtr<plink4r::KinNull> null(new plink4r::KinNull, true);
  plink4r::KinNull& M = *null;
  try {
    M.init(plink4r::build_kin(n, ki.begin(), kj.begin(), kv.begin(), ki.size(), kdiag));
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }

  const double big = 1e300;
  const double h_top = M.h_bound(h_max, tol);  // H(h) is positive definite on [0, h_top]
  const double* Xp = X.begin();
  std::vector<double> Xc, eta(N, 0.0), we(N), jump(ng), sw(N), yt(N), ys(N), Xs(N * p);
  double emax = 0;
  // Breslow jumps of Lambda0 (times exp(emax)) at the current eta
  auto breslow = [&]() {
    emax = *std::max_element(eta.begin(), eta.end());
    for (int i = 0; i < n; ++i) we[i] = std::exp(eta[i] - emax);
    double s0 = 0;
    for (int pos = 0, g = 0; g < ng; ++g) {
      for (; pos <= gend[g]; ++pos) s0 += we[ord[pos]];
      jump[g] = gd[g] / s0;
    }
  };
  double slw = 0;  // sum log w of the working data
  // working data at the current eta, scaled by S: Lambda0, mu = w, yt
  auto working = [&]() {
    breslow();
    slw = 0;
    double cum = 0;  // Lambda0 at the time of the sample (times exp(emax)), from the earliest time on
    int g = ng - 1;
    for (int pos = n - 1; pos >= 0; --pos) {
      while (g >= 0 && gend[g] >= pos) cum += jump[g--];
      const int i = ord[pos];
      const double m = std::max(mu_eps, we[i] * std::max(cum, jump[ng - 1]));
      sw[i] = std::sqrt(m);
      slw += std::log(m);
      yt[i] = eta[i] + (status[i] - m) / m;
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
  // eta <- X beta without the intercept at the given h; -1 on failure
  auto update = [&](double h) {
    if (!M.solve(h)) return -1.0;
    Eigen::Map<const Eigen::MatrixXd> Xm(Xp, n, p);
    Eigen::Map<Eigen::VectorXd> ev(eta.data(), n);
    ev.noalias() = Xm.rightCols(p - 1) * M.beta.tail(p - 1);
    return 0.0;
  };

  // aliased columns (the weights are positive, so the rank is that of X)
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
    double prev = big, cur = big;
    for (int it = 0; it < max_iter; ++it) {
      if (steps > 0) working();
      ++steps;
      if (update(h) < 0) return big;
      cur = M.logdet + slw + M.logdet_xox + M.quad;
      if (std::fabs(cur - prev) <= f_tol * (1.0 + std::fabs(cur))) return cur;
      prev = cur;
    }
    conv = false;
    return cur;
  };
  const double f0 = f(0.0);  // the Cox model without kinship
  if (f0 >= big || !conv) Rcpp::stop("the Cox null model cannot be fitted");
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
  const double df = f_tol * (1.0 + std::fabs(f0));  // changes of f below this are not told apart
  double fh = f0;
  if (M.K.nr > 0) {
    const double lo = std::max(0.0, h - h_window), hi = std::min(h_top, h + h_window);
    plink4r::BrentResult br = plink4r::brent_min(f, lo, hi, tol, tol * 1e-2, max_eval, df);
    const double edge = 0.01 * h_window;
    if ((lo > 0 && br.x - lo < edge) || (hi < h_top && hi - br.x < edge))
      br = plink4r::brent_min(f, 0.0, h_top, tol, tol * 1e-2, max_eval, df);
    h = f0 - br.fx > df ? br.x : 0.0;
  } else {
    h = 0.0;
  }
  fh = f(h);  // the fit at h is the one left in M
  if (h == 0.0) fh = f0;
  if (fh >= big) Rcpp::stop("the model cannot be fitted at the estimate");
  if (proj) M.project(Xp, yt.data(), sw.data(), 1.0);

  // the baseline at the fit, in increasing time
  breslow();
  Rcpp::NumericVector bt(ng), bh(ng);
  double cum = 0;
  for (int g = ng - 1; g >= 0; --g) {
    cum += jump[g];
    bt[ng - 1 - g] = time[ord[gend[g]]];
    bh[ng - 1 - g] = cum * std::exp(-emax);
  }

  Rcpp::IntegerVector kp(keep.begin(), keep.end());
  kp = kp + 1;
  return Rcpp::List::create(
      Rcpp::Named("tau") = h, Rcpp::Named("h_bound") = h_top, Rcpp::Named("m2l") = fh, Rcpp::Named("m2l_0") = f0,
      Rcpp::Named("iterations") = steps, Rcpp::Named("evals") = evals, Rcpp::Named("converged") = conv,
      Rcpp::Named("beta") = Rcpp::NumericVector(M.beta.data(), M.beta.data() + p), Rcpp::Named("keep") = kp,
      Rcpp::Named("base_time") = bt, Rcpp::Named("base_cumhaz") = bh,
      Rcpp::Named("n_related") = M.K.nr, Rcpp::Named("pairs") = static_cast<double>(M.K.pairs),
      Rcpp::Named("null") = null);
}
