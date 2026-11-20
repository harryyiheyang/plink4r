// REML of the kinship variance tau of a logistic mixed model with a sparse kinship, on the IRLS working scale
// (penalised quasi-likelihood, as GMMAT / fastGWA-GLMM):
//   logit(mu) = eta = X beta + b,  b ~ N(0, tau K).
// At the current eta: w = mu (1 - mu), working response yt = eta + (y - mu) / w, and the working linear mixed
// model yt = X beta + b + e, e ~ N(0, W^-1), i.e. V = W^-1 + tau K. With w and yt held fixed, tau minimises
//   -2 l_R(tau) = log|V| + log|X' V^-1 X| + yt' P yt
// by one Brent search on [0, tau_max]; then beta and b are updated: eta = X beta + b = yt - W^-1 P yt.
// The two steps alternate until eta and tau settle. tau = 0 throughout is ordinary logistic IRLS, which is how
// the iteration is started. At convergence P yt = y - mu, the score residual of the association test.
// The weights differ between samples, so every evaluation forms the weighted cross-products over all samples
// and factorises the related samples' block of V (kin_sparse.h, KinNull::solve).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"
#include "resources.h"

// y: n values 0 / 1; X: n x p (with the intercept); kinship and rank_tol as reml_tau_gaussian_cpp.
// tau_max: upper end of the search; tol: Brent tolerance and the relative change of tau at convergence;
// eta_tol: largest change of eta at convergence; max_iter: outer iterations (each of the two phases);
// mu_eps: fitted probabilities are kept in [mu_eps, 1 - mu_eps]. proj: also the null model for the association
// stream. threads: OpenMP threads of the weighted cross-products (0 = from the resources).
// Returns tau, the -2 working REML log-likelihood at tau and at 0 (same working data), beta, keep (1-based
// columns of X used) and null (external pointer to the plink4r::KinNull, fitted at V = W^-1 + tau K).
// [[Rcpp::export]]
Rcpp::List reml_tau_binomial_cpp(const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& X,
                                 const Rcpp::IntegerVector& ki, const Rcpp::IntegerVector& kj,
                                 const Rcpp::NumericVector& kv, double kdiag, double tau_max, double tol,
                                 double eta_tol, int max_iter, int max_eval, double mu_eps, double rank_tol,
                                 bool proj, int threads) {
  const int n = y.size();
  int p = X.ncol();
  const std::size_t N = n;
  if (X.nrow() != n) Rcpp::stop("y / X do not match");
  if (ki.size() != kj.size() || ki.size() != kv.size()) Rcpp::stop("kinship triplets do not match");
  if (!(tau_max > 0)) Rcpp::stop("tau_max must be > 0");
  double ybar = 0;
  for (int i = 0; i < n; ++i) ybar += y[i];
  ybar /= n;
  if (!(ybar > 0 && ybar < 1)) Rcpp::stop("the phenotype has one class only");
  Rcpp::XPtr<plink4r::KinNull> null(new plink4r::KinNull, true);
  plink4r::KinNull& M = *null;
  try {
    M.init(plink4r::build_kin(n, ki.begin(), kj.begin(), kv.begin(), ki.size(), kdiag));
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
  M.threads = plink4r::resolve_threads(threads);

  const double big = 1e300;
  const double* Xp = X.begin();
  std::vector<double> Xc, eta(N, std::log(ybar / (1 - ybar))), d(N), yt(N), r(N);
  // working data at the current eta
  auto working = [&]() {
    for (int i = 0; i < n; ++i) {
      double m = 1.0 / (1.0 + std::exp(-eta[i]));
      m = std::min(1.0 - mu_eps, std::max(mu_eps, m));
      d[i] = 1.0 / (m * (1.0 - m));
      yt[i] = eta[i] + (y[i] - m) * d[i];
    }
    M.set_data(Xp, yt.data(), p);
  };
  auto m2l = [&](double tau) {
    if (!M.solve(tau, d.data())) return big;
    return M.logdet + M.logdet_xox + M.quad;
  };
  // eta <- yt - W^-1 P yt at the given tau; returns the largest change of eta (-1 on failure)
  auto update = [&](double tau) {
    if (!M.solve(tau, d.data())) return -1.0;
    M.residual(r.data());
    double ch = 0;
    for (int i = 0; i < n; ++i) {
      const double e = yt[i] - d[i] * r[i];
      ch = std::max(ch, std::fabs(e - eta[i]));
      eta[i] = e;
    }
    return ch;
  };

  // aliased columns of X
  working();
  M.cross();
  const std::vector<int> keep = M.independent(rank_tol);
  if (static_cast<int>(keep.size()) < p) {
    p = static_cast<int>(keep.size());
    Xc.resize(N * p);
    for (int c = 0; c < p; ++c) std::copy(Xp + keep[c] * N, Xp + (keep[c] + 1) * N, &Xc[c * N]);
    Xp = Xc.data();
    working();
  }
  if (n <= p) Rcpp::stop("more covariates than samples");

  // phase 1: tau = 0 (logistic regression)
  bool conv0 = false;
  for (int it = 0; it < max_iter; ++it) {
    if (it > 0) working();
    const double ch = update(0.0);
    if (ch < 0) Rcpp::stop("the logistic null model cannot be fitted (separation?)");
    if (ch < eta_tol) { conv0 = true; break; }
  }
  // phase 2: alternate tau and (beta, b)
  double tau = 0;
  int iters = 0, evals = 0;
  bool conv = false;
  for (int it = 0; it < max_iter; ++it) {
    ++iters;
    working();
    const double f0 = m2l(0.0);
    plink4r::BrentResult br = plink4r::brent_min(m2l, 0.0, tau_max, tol, tol * 1e-2, max_eval);
    evals += br.evals + 1;
    const double tau_new = br.fx < f0 ? br.x : 0.0;
    const double ch = update(tau_new);
    if (ch < 0) Rcpp::stop("the working mixed model cannot be fitted");
    const bool done = ch < eta_tol && std::fabs(tau_new - tau) <= tol * (1.0 + tau);
    tau = tau_new;
    if (done) { conv = true; break; }
  }
  working();
  const double f0 = m2l(0.0), ft = m2l(tau);  // the fit at tau is the one left in M
  if (ft >= big) Rcpp::stop("V is not positive definite at the estimate");
  if (proj) M.project();

  Rcpp::IntegerVector kp(keep.begin(), keep.end());
  kp = kp + 1;
  return Rcpp::List::create(
      Rcpp::Named("tau") = tau, Rcpp::Named("m2l") = ft, Rcpp::Named("m2l_0") = f0,
      Rcpp::Named("iterations") = iters, Rcpp::Named("evals") = evals, Rcpp::Named("converged") = conv && conv0,
      Rcpp::Named("beta") = Rcpp::NumericVector(M.beta.data(), M.beta.data() + p), Rcpp::Named("keep") = kp,
      Rcpp::Named("n_related") = M.K.nr, Rcpp::Named("pairs") = static_cast<double>(M.K.pairs),
      Rcpp::Named("null") = null);
}
