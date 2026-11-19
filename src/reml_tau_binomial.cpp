// REML of the kinship correlation h of a logistic model with a sparse kinship, on the IRLS working scale,
// marginal form: logit(mu) = eta = X beta, and the Pearson residuals (y - mu) / sqrt(mu (1 - mu)) have
// covariance H(h) = (1 - h) I + h K. At the current eta: w = mu (1 - mu), S = W^1/2, working response
// yt = eta + (y - mu) / w with covariance V = S^-1 H(h) S^-1, so V^-1 = S H^-1 S. With w and yt held fixed,
// h minimises
//   -2 l_R(h) = log|V| + log|X' V^-1 X| + yt' P yt
// by one Brent search on [0, h_max]; then beta is updated by generalised least squares, eta = X beta. The two
// steps alternate until eta and h settle. A sample without a relative has V_ii = (1 - h + h k_ii) / w_i, which is
// 1 / w_i whatever h when k_ii = 1: such samples carry no information on h, and at h = 0 the fit is ordinary
// logistic regression.
// Per outer iteration the scaled data S X, S yt and their cross-products are formed once; an evaluation at h
// factorises the related samples' block of H and solves with their rows only (kin_sparse.h). K itself is never
// rescaled. (log|V| = log|H| - sum log w; the sum does not depend on h and is left out of both likelihoods.)
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"

// y: n values 0 / 1; X: n x p (with the intercept); kinship and rank_tol as reml_tau_gaussian_cpp.
// h_max: upper end of the search; tol: Brent tolerance and the change of h at convergence; eta_tol: largest
// change of eta at convergence; max_iter: outer iterations (each of the two phases); mu_eps: fitted
// probabilities are kept in [mu_eps, 1 - mu_eps]. proj: also the null model for the association stream.
// Returns h (= tau), the -2 working REML log-likelihood at h and at 0 (same working data), beta, keep (1-based
// columns of X used) and null (external pointer to the plink4r::KinNull; V^-1 = W^1/2 H^-1 W^1/2).
// [[Rcpp::export]]
Rcpp::List reml_tau_binomial_cpp(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& X,
                                 const Rcpp::IntegerVector& ki, const Rcpp::IntegerVector& kj,
                                 const Rcpp::NumericVector& kv, double kdiag, double h_max, double tol,
                                 double eta_tol, int max_iter, int max_eval, double mu_eps, double rank_tol,
                                 bool proj) {
  const int n = y.size();
  int p = X.ncol();
  const std::size_t N = n;
  if (X.nrow() != n) Rcpp::stop("y / X do not match");
  if (ki.size() != kj.size() || ki.size() != kv.size()) Rcpp::stop("kinship triplets do not match");
  if (!(h_max > 0 && h_max < 1)) Rcpp::stop("h_max must lie in (0, 1)");
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
  const double* Xp = X.begin();
  std::vector<double> Xc, eta(N, std::log(ybar / (1 - ybar))), sw(N), yt(N), ys(N), Xs(N * p);
  // working data at the current eta, scaled by S
  auto working = [&]() {
    for (int i = 0; i < n; ++i) {
      double m = 1.0 / (1.0 + std::exp(-eta[i]));
      m = std::min(1.0 - mu_eps, std::max(mu_eps, m));
      const double w = m * (1.0 - m);
      sw[i] = std::sqrt(w);
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

  // phase 1: h = 0 (logistic regression)
  bool conv0 = false;
  for (int it = 0; it < max_iter; ++it) {
    if (it > 0) working();
    const double ch = update(0.0);
    if (ch < 0) Rcpp::stop("the logistic null model cannot be fitted (separation?)");
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
  const double f0 = m2l(0.0), fh = m2l(h);  // the fit at h is the one left in M
  if (fh >= big) Rcpp::stop("V is not positive definite at the estimate");
  if (proj) M.project(Xp, yt.data(), sw.data(), 1.0);

  Rcpp::IntegerVector kp(keep.begin(), keep.end());
  kp = kp + 1;
  return Rcpp::List::create(
      Rcpp::Named("tau") = h, Rcpp::Named("m2l") = fh, Rcpp::Named("m2l_0") = f0,
      Rcpp::Named("iterations") = iters, Rcpp::Named("evals") = evals, Rcpp::Named("converged") = conv && conv0,
      Rcpp::Named("beta") = Rcpp::NumericVector(M.beta.data(), M.beta.data() + p), Rcpp::Named("keep") = kp,
      Rcpp::Named("n_related") = M.K.nr, Rcpp::Named("pairs") = static_cast<double>(M.K.pairs),
      Rcpp::Named("null") = null);
}
