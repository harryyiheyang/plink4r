// The kinship correlation h of a Cox model with a sparse kinship, on the life table of the data.
//
// Records. Follow-up is cut into J intervals at quantiles of the event times, and a sample contributes one binary
// record to every interval at whose start it is still at risk: d_ij = 1 when its event falls in interval j, 0
// otherwise (a sample censored within an interval counts there as a non-event). Under proportional hazards
//   P(d_ij = 1) = p_ij = 1 - exp(-exp(alpha_j + x_i' beta))      (complementary log-log link),
// the alpha_j being the baseline (exp(alpha_j): the baseline hazard cumulated over interval j). The weight of a
// record, w_ij = (dp/deta)^2 / (p (1 - p)), depends on the interval and the covariates only, not on the outcome;
// with the Poisson form of the Cox likelihood it is the expected events up to the observed time, which the outcome
// determines, and a correlation of Pearson residuals is then not defined.
//
// Covariance. The records of one sample are uncorrelated (increments of a martingale); the Pearson residuals of
// the records of two relatives a, b have correlation h K_ab, whatever the two intervals:
//   Var(d) = W^1/2 R W^1/2,   R = I + h E Koff E',   E: records x samples incidence, Koff: K without its diagonal.
// With n_a the number of records of sample a, D = diag(n_a) and Et = E D^-1/2 (orthonormal columns),
//   R = (I - Et Et') + Et H Et',   H(h) = I + h D^1/2 Koff D^1/2:
// R is the identity on the contrasts within a sample and H on the (normalised) sums over a sample. So log|R| =
// log|H|, R^-1 = I + Et (H^-1 - I) Et', and the model is that of plink4r::KinNull (kin_sparse.h) for the kinship
// with entries sqrt(n_a n_b) K_ab and unit diagonal, given the cross-products of the scaled records and, for the
// related samples, their normalised sums: a search over h costs one sparse factorisation of H per value, as in
// the logistic model, whatever J.
//
// Estimation alternates two steps until the criterion settles (iteratively reweighted least squares). At the
// current alpha, beta, with the weights and the working response zt = eta + (d - p) / (dp/deta) held fixed, h
// minimises the working restricted likelihood
//   -2 l_R(h) = log|R| + log|X' W^1/2 R^-1 W^1/2 X| + zt' P zt        (X: interval indicators and covariates)
// by one Brent search on [0, h_bound]; then alpha, beta are updated by generalised least squares of zt on X at
// that h. First h = 0 until the fit settles (the life-table model without kinship), then h is searched.
// Convergence: the criterion at the chosen h changes by less than f_tol (1 + |criterion|) between two iterations.
// h is exactly 0 when the search does not improve the criterion at 0 by more than that, and with no related pair.
//
// Association (gwas_kin_stream.cpp). Given h the test of a variant g is the score test of this model with g added
// to every record of its sample (a GEE with the working covariance above). The records collapse to the samples:
// with W_a = sum_j w_aj, t_a = sum_j sqrt(w_aj) / sqrt(n_a),
//   g' Omega g = sum_a (W_a - t_a^2 [a related]) g_a^2 + |L^-1 (t_r o g_r)|^2,   H = L L',
// and Q = Omega X La^-T, r = Omega (zt - X (alpha, beta)) are summed over the records of each sample, so the
// stream reads om, L1 = [Q | r], ss = t and the factor of H as for the other models. The columns of X are not
// constant within a sample, so there is no sample-level Vb.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "brent.h"
#include "kin_sparse.h"
#include "la_dispatch.h"

// time, status: n follow-up times and event indicators 0 / 1; X: n x p with the intercept as its first column
// (the interval indicators take its place); ki, kj, kv: off-diagonal kinship entries (0-based samples, each pair
// once), unit diagonal. n_interval: intervals of the life table (fewer when there are fewer distinct event times).
// h_max: upper end of the search (or the largest h with H(h) positive definite, h_bound); tol: Brent tolerance
// on h; max_iter: iterations; max_eval: evaluations per Brent search; mu_eps: the record probabilities are kept in
// [mu_eps, 1 - mu_eps]; f_tol: as above; rank_tol: as reml_tau_gaussian_cpp. proj: also the null model for the
// association stream.
// Returns h (= tau), h_bound, the criterion at h and at 0 (same working data), beta (the J interval intercepts,
// then the covariates), the baseline (base_time: the ends of the intervals; base_cumhaz: the cumulative baseline
// hazard there), the iterations and evaluations made, keep (1-based columns of X used, the intercept among them)
// and null (external pointer to the plink4r::KinNull).
// [[Rcpp::export]]
Rcpp::List reml_tau_cox_cpp(const Rcpp::NumericVector& time, const Rcpp::NumericVector& status,
                            const Rcpp::NumericMatrix& X, const Rcpp::IntegerVector& ki,
                            const Rcpp::IntegerVector& kj, const Rcpp::NumericVector& kv, int n_interval,
                            double h_max, double tol, int max_iter, int max_eval, double mu_eps, double f_tol,
                            double rank_tol, bool proj) {
  const int n = time.size();
  int pc = X.ncol() - 1;  // covariates
  const std::size_t N = n;
  if (X.nrow() != n || status.size() != n || pc < 0) Rcpp::stop("time / status / X do not match");
  if (ki.size() != kj.size() || ki.size() != kv.size()) Rcpp::stop("kinship triplets do not match");
  if (!(h_max > 0 && h_max < 1)) Rcpp::stop("h_max must lie in (0, 1)");
  if (n_interval < 1) Rcpp::stop("n_interval must be >= 1");
  for (int i = 0; i < n; ++i)
    if (X(i, 0) != 1.0) Rcpp::stop("the first column of X must be the intercept");

  // ---- intervals: cuts at quantiles of the event times (linear interpolation of the order statistics); interval
  // j is (cuts[j - 1], cuts[j]], the first from -inf, the last to +inf; an interval without an event is merged
  std::vector<double> ev;
  for (int i = 0; i < n; ++i)
    if (status[i] == 1) ev.push_back(time[i]);
  std::sort(ev.begin(), ev.end());
  const int ne = static_cast<int>(ev.size());
  if (ne < 1) Rcpp::stop("no event");
  std::vector<double> cuts;
  const int J0 = std::min(n_interval, ne);
  for (int k = 1; k < J0; ++k) {
    const double hp = (ne - 1) * (static_cast<double>(k) / J0);
    const int lo = static_cast<int>(std::floor(hp)), hi = std::min(lo + 1, ne - 1);
    cuts.push_back(ev[lo] + (hp - lo) * (ev[hi] - ev[lo]));
  }
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  auto interval = [&](double t) {
    return static_cast<int>(std::lower_bound(cuts.begin(), cuts.end(), t) - cuts.begin());
  };
  std::vector<double> nev;  // events per interval
  for (;;) {
    nev.assign(cuts.size() + 1, 0.0);
    for (double t : ev) nev[interval(t)] += 1;
    std::size_t k = 0;
    while (k < nev.size() && nev[k] > 0) ++k;
    if (k == nev.size()) break;
    cuts.erase(cuts.begin() + (k < cuts.size() ? k : k - 1));
  }
  const int J = static_cast<int>(cuts.size()) + 1;
  std::vector<int> nint(n);     // records of a sample: intervals 0 .. nint - 1
  std::vector<double> nrisk(J, 0.0);
  double tmax = 0;
  for (int i = 0; i < n; ++i) {
    nint[i] = interval(time[i]) + 1;
    for (int j = 0; j < nint[i]; ++j) nrisk[j] += 1;
    tmax = std::max(tmax, static_cast<double>(time[i]));
  }

  // ---- kinship between the sums over the samples: sqrt(n_a n_b) K_ab, unit diagonal
  std::vector<double> kvs(kv.size());
  for (R_xlen_t e = 0; e < kv.size(); ++e) {
    if (ki[e] < 0 || ki[e] >= n || kj[e] < 0 || kj[e] >= n) Rcpp::stop("kinship pair out of range");
    kvs[e] = kv[e] * std::sqrt(static_cast<double>(nint[ki[e]]) * nint[kj[e]]);
  }
  Rcpp::XPtr<plink4r::KinNull> null(new plink4r::KinNull, true);
  plink4r::KinNull& M = *null;
  try {
    M.init(plink4r::build_kin(n, ki.begin(), kj.begin(), kvs.data(), ki.size(), 1.0));
  } catch (const std::exception& e) { Rcpp::stop(e.what()); }
  const int nr = M.K.nr;
  const std::vector<int>& rel = M.K.rel;

  const double big = 1e300;
  const plink4r::LaFns& LA = plink4r::la();
  const double* Xp = X.begin() + N;  // the covariates, n x pc
  int P = J + pc;
  std::vector<double> Xc, par(P, 0.0), xb(N), Wm(N * J, 0.0), Zm(N * J, 0.0), Wp(N), WZp(N), Xs, tt(nr);
  for (int j = 0; j < J; ++j) {  // start: the life table without covariates
    const double rate = std::min(1.0 - 1e-6, nev[j] / nrisk[j]);
    par[j] = std::log(-std::log1p(-rate));
  }
  Eigen::MatrixXd G, Zr;
  // working data at the current alpha, beta: weights, working response, their cross-products, the related sums
  auto working = [&]() {
    const int P1 = P + 1;
    Eigen::Map<const Eigen::MatrixXd> Xm(Xp, n, pc);
    Eigen::Map<Eigen::VectorXd> xbv(xb.data(), n);
    if (pc > 0) xbv.noalias() = Xm * Eigen::Map<const Eigen::VectorXd>(par.data() + J, pc); else xbv.setZero();
    G.setZero(P1, P1);
    double zz = 0;
    for (int i = 0; i < n; ++i) {
      double Wi = 0, WZi = 0;
      const int ni = nint[i];
      for (int j = 0; j < ni; ++j) {
        const double eta = par[j] + xb[i];
        const double pr = std::min(1.0 - mu_eps, std::max(mu_eps, -std::expm1(-std::exp(eta))));
        const double dp = -(1.0 - pr) * std::log1p(-pr);  // dp / deta
        const double d = (j == ni - 1 && status[i] == 1) ? 1.0 : 0.0;
        const double w = dp * dp / (pr * (1.0 - pr)), zt = eta + (d - pr) / dp;
        Wm[j * N + i] = w;
        Zm[j * N + i] = zt;
        G(j, j) += w;
        G(j, P) += w * zt;
        zz += w * zt * zt;
        Wi += w;
        WZi += w * zt;
      }
      Wp[i] = Wi;
      WZp[i] = WZi;
    }
    G(P, P) = zz;
    if (pc > 0) {
      Xs.resize(N * pc);
      for (int c = 0; c < pc; ++c) {
        const double* x = Xp + c * N;
        double* xs = &Xs[c * N];
        for (std::size_t i = 0; i < N; ++i) xs[i] = std::sqrt(Wp[i]) * x[i];
      }
      LA.dgemm_tn(Wm.data(), n, J, n, Xp, pc, n, G.data() + static_cast<std::size_t>(J) * P1, P1, 0);
      LA.dgemm_tn(Xs.data(), n, pc, n, Xs.data(), pc, n, G.data() + static_cast<std::size_t>(J) * P1 + J, P1, 0);
      LA.dgemm_tn(Xp, n, pc, n, WZp.data(), 1, n, G.data() + static_cast<std::size_t>(P) * P1 + J, P1, 0);
    }
    for (int a = 0; a < P1; ++a)
      for (int b = 0; b < a; ++b) G(a, b) = G(b, a);
    Zr.setZero(nr, P1);
    for (int a = 0; a < nr; ++a) {
      const int i = rel[a], ni = nint[i];
      const double rn = 1.0 / std::sqrt(static_cast<double>(ni));
      double t = 0, sz = 0;
      for (int j = 0; j < ni; ++j) {
        const double s = std::sqrt(Wm[j * N + i]);
        Zr(a, j) = s * rn;
        t += s;
        sz += s * Zm[j * N + i];
      }
      t *= rn;
      tt[a] = t;
      for (int c = 0; c < pc; ++c) Zr(a, J + c) = t * Xp[c * N + i];
      Zr(a, P) = sz * rn;
    }
    M.set_cross(G, Zr, P);
  };
  // -2 l_R(h) of the working data
  auto m2l = [&](double h) {
    if (!M.solve(h)) return big;
    return M.logdet + M.logdet_xox + M.quad;
  };
  // alpha, beta by generalised least squares at h; false on failure
  auto update = [&](double h) {
    if (!M.solve(h)) return false;
    for (int a = 0; a < P; ++a) par[a] = M.beta[a];
    return true;
  };

  // aliased covariate columns (the weights are positive, so the rank is that of the design)
  working();
  std::vector<int> keepx(1, 0);  // columns of X used
  {
    const std::vector<int> keep = M.independent(rank_tol);
    std::vector<int> kc;  // covariates kept
    int nd = 0;
    for (int a : keep) {
      if (a < J) ++nd; else { kc.push_back(a - J); keepx.push_back(a - J + 1); }
    }
    if (nd < J) Rcpp::stop("an interval of the life table carries no information");
    if (static_cast<int>(kc.size()) < pc) {
      const int pk = static_cast<int>(kc.size());
      Xc.resize(N * std::max(pk, 1));
      for (int c = 0; c < pk; ++c) std::copy(Xp + kc[c] * N, Xp + (kc[c] + 1) * N, &Xc[c * N]);
      Xp = Xc.data();
      pc = pk;
      P = J + pc;
      par.resize(P);
      working();
    }
  }
  if (n <= P) Rcpp::stop("more covariates than samples");

  // h = 0 first (the life-table model without kinship), then h and alpha, beta in turn
  const double tol_h = tol / J;                 // h is a correlation of records: its scale is 1 / J
  const double top = M.h_bound(h_max, tol_h);  // H(h) is positive definite on [0, top]
  double h = 0, f0 = big, fh = big, prev = big;
  int iters = 0, evals = 0;
  bool conv = false, free_h = false;  // free_h: h is searched (after the fit at 0 has settled)
  for (int it = 0; it < max_iter; ++it) {
    if (it > 0) working();
    ++iters;
    f0 = m2l(0.0);
    if (f0 >= big) Rcpp::stop("the Cox null model cannot be fitted");
    fh = f0;
    h = 0;
    const double df = f_tol * (1.0 + std::fabs(f0));  // changes of the criterion below this are not told apart
    if (free_h) {
      const plink4r::BrentResult br = plink4r::brent_min(m2l, 0.0, top, tol, tol_h * 1e-2, max_eval, df);
      evals += br.evals;
      if (f0 - br.fx > df) { h = br.x; fh = br.fx; }
    }
    if (!update(h)) Rcpp::stop("the working model cannot be fitted");
    const bool settled = std::fabs(fh - prev) <= f_tol * (1.0 + std::fabs(fh));
    prev = fh;
    if (!settled) continue;
    if (free_h || nr == 0) { conv = true; break; }
    free_h = true;  // the fit without kinship has settled: start the search from it
    prev = big;
  }
  // the fit at h on the working data of the final alpha, beta
  working();
  f0 = m2l(0.0);
  fh = h > 0 ? m2l(h) : f0;
  if (fh >= big) Rcpp::stop("H is not positive definite at the estimate");

  if (proj) {
    const int P1 = P + 1;
    Eigen::MatrixXd Li = Eigen::MatrixXd::Identity(P, P);
    M.xox.matrixL().solveInPlace(Li);  // La^-1
    M.om.assign(Wp.begin(), Wp.end());
    M.Vb.clear();
    M.L1.assign(N * P1, 0.0);
    Eigen::Map<Eigen::MatrixXd> Lm(M.L1.data(), n, P1);
    Eigen::Map<const Eigen::MatrixXd> Xm(Xp, n, pc), Wmm(Wm.data(), n, J);
    Eigen::Map<const Eigen::VectorXd> Wv(Wp.data(), n), WZv(WZp.data(), n);
    // sums over the records of a sample of W X La^-T and of W (zt - X (alpha, beta))
    Lm.leftCols(P).noalias() = Wmm * Li.leftCols(J).transpose();
    Lm.col(P) = WZv - Wmm * M.beta.head(J);
    if (pc > 0) {
      Lm.leftCols(P).noalias() += Wv.asDiagonal() * (Xm * Li.rightCols(pc).transpose());
      Lm.col(P).array() -= Wv.array() * (Xm * M.beta.tail(pc)).array();
    }
    // the related samples: + t_a ((H^-1 - I) applied to the normalised sums)_a
    M.rs.resize(nr);
    M.ss.resize(nr);
    if (nr > 0) {
      Eigen::MatrixXd T(nr, P1), Tp(nr, P1);
      T.leftCols(P).noalias() = Zr.leftCols(P) * Li.transpose();
      T.col(P) = Zr.col(P) - Zr.leftCols(P) * M.beta;
      for (int q = 0; q < nr; ++q) Tp.row(M.llt.permutationP().indices()[q]) = T.row(q);
      M.llt.matrixL().solveInPlace(Tp);
      M.llt.matrixU().solveInPlace(Tp);
      for (int q = 0; q < nr; ++q) {
        const int pos = M.llt.permutationP().indices()[q], i = rel[q];
        Lm.row(i) += tt[q] * (Tp.row(pos) - T.row(q));
        M.om[i] = std::max(0.0, Wp[i] - tt[q] * tt[q]);
        M.rs[pos] = i;
        M.ss[pos] = tt[q];
      }
    }
  }

  // the baseline: cumulative hazard at the end of every interval
  Rcpp::NumericVector bt(J), bh(J);
  double cum = 0;
  for (int j = 0; j < J; ++j) {
    cum += std::exp(M.beta[j]);
    bt[j] = j < J - 1 ? cuts[j] : tmax;
    bh[j] = cum;
  }
  Rcpp::IntegerVector kp(keepx.begin(), keepx.end());
  kp = kp + 1;
  return Rcpp::List::create(
      Rcpp::Named("tau") = h, Rcpp::Named("h_bound") = top, Rcpp::Named("m2l") = fh, Rcpp::Named("m2l_0") = f0,
      Rcpp::Named("iterations") = iters, Rcpp::Named("evals") = evals, Rcpp::Named("converged") = conv,
      Rcpp::Named("beta") = Rcpp::NumericVector(M.beta.data(), M.beta.data() + P), Rcpp::Named("keep") = kp,
      Rcpp::Named("n_interval") = J, Rcpp::Named("base_time") = bt, Rcpp::Named("base_cumhaz") = bh,
      Rcpp::Named("n_related") = nr, Rcpp::Named("pairs") = static_cast<double>(M.K.pairs),
      Rcpp::Named("null") = null);
}
