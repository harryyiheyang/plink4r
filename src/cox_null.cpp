// The null model of a Cox model for unrelated samples: the partial-likelihood fit (Breslow ties) on the Poisson
// form of the Cox likelihood. With the cumulative baseline hazard Lambda0, the event indicator d_i of a sample
// followed to time t_i has mean
//   mu_i = Lambda0(t_i) exp(eta_i),  eta = X beta,
// and the Cox likelihood in beta at a fixed Lambda0 is that of a Poisson model for d with log Lambda0(t_i) as an
// offset. One scoring step: at the current eta, Lambda0 is the Breslow estimate (jumps d_k / sum of exp(eta) over
// the risk set), which gives mu; the working response yt = eta + (d - mu) / mu has variance 1 / mu, and weighted
// least squares of yt on [1 | X] gives beta, eta = X beta. The intercept only rescales the baseline (mu does not
// change with it) and is left out of eta; at the solution X' (d - mu) = 0 with X centred by the intercept, the
// partial-likelihood score equation. Steps are repeated until the criterion
//   f = sum log mu + log|X1' W X1| + sum (d - mu)^2 / mu   (X1 = [1 | X], W = diag(mu))
// changes by less than f_tol (1 + |f|).
// The association stream then tests every variant in the Poisson model with the offset log Lambda0_hat and the
// intercept (gwas_kin_stream.cpp, with a null-model object that has no related sample).
// A sample that leaves before the first event has Lambda0 = 0, hence mu = 0: it is in no risk set and carries no
// information (its mu is held at the lower bound mu_eps, as every mu).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "kin_sparse.h"

// time, status: n follow-up times and event indicators 0 / 1; X: n x p with the intercept as its first column.
// max_iter: scoring steps; mu_eps: lower bound of the expected events mu; f_tol: as above; rank_tol: as
// reml_tau_gaussian_cpp. proj: also the null model for the association stream.
// Returns beta (the intercept first: the log of the factor by which the Poisson fit rescales the baseline), the
// baseline (base_time: distinct event times, increasing; base_cumhaz: Lambda0 there), the criterion (m2l), the
// scoring steps made (iterations), keep (1-based columns of X used) and null (external pointer to the
// plink4r::KinNull; V^-1 = W).
// [[Rcpp::export]]
Rcpp::List cox_null_cpp(const Rcpp::NumericVector& time, const Rcpp::NumericVector& status,
                        const Rcpp::NumericMatrix& X, int max_iter, double mu_eps, double f_tol, double rank_tol,
                        bool proj) {
  const int n = time.size();
  int p = X.ncol();
  const std::size_t N = n;
  if (X.nrow() != n || status.size() != n || p < 1) Rcpp::stop("time / status / X do not match");
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
  M.init(plink4r::build_kin(n, nullptr, nullptr, nullptr, 0, 1.0));  // no related sample

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
      const double m = std::max(mu_eps, we[i] * cum);
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

  // scoring to the partial-likelihood fit; the fit, eta and working data are left at it
  const double big = 1e300;
  int steps = 0;
  bool conv = false;
  double prev = big, cur = big;
  for (int it = 0; it < max_iter; ++it) {
    if (it > 0) working();
    ++steps;
    if (!M.solve(0.0)) Rcpp::stop("the Cox null model cannot be fitted");
    Eigen::Map<const Eigen::MatrixXd> Xm(Xp, n, p);
    Eigen::Map<Eigen::VectorXd> ev(eta.data(), n);
    ev.noalias() = Xm.rightCols(p - 1) * M.beta.tail(p - 1);  // without the intercept
    cur = slw + M.logdet_xox + M.quad;
    if (std::fabs(cur - prev) <= f_tol * (1.0 + std::fabs(cur))) { conv = true; break; }
    prev = cur;
  }
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
      Rcpp::Named("m2l") = cur, Rcpp::Named("iterations") = steps, Rcpp::Named("converged") = conv,
      Rcpp::Named("beta") = Rcpp::NumericVector(M.beta.data(), M.beta.data() + p), Rcpp::Named("keep") = kp,
      Rcpp::Named("base_time") = bt, Rcpp::Named("base_cumhaz") = bh, Rcpp::Named("null") = null);
}
