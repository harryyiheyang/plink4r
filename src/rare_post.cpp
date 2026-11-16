// post_mix_mv_cpp: posterior mean and sd of the S logistic G x E coefficients of every variant under a product
// mixture prior (rare_variant_calibration.md sec. 15). Main-effect component k and interaction component l give the
// prior N(0, D C D), D = diag(d_k, d_l, ..., d_l), d^2 = scaled * base_j + abs, C = I / (tr(I) / S) the shape of
// the variant's information I (optionally its eigenvalues raised to gamma, trace kept at S). With the score
// u = I bhat, for a prior N(0, Sig): posterior mean Sig (1 + I Sig)^-1 u, covariance Sig (1 + I Sig)^-1, and log Bayes
// factor against beta = 0: -log det(1 + I Sig) / 2 + u' mean / 2 (no inverse of I needed, so a near-singular I
// of a rare variant is fine). Component weights w_k v_l exp(lbf); the mixture moments by the law of total variance.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "resources.h"

// [[Rcpp::export]]
Rcpp::List post_mix_mv_cpp(Rcpp::NumericMatrix score, Rcpp::NumericMatrix info, Rcpp::NumericVector base,
                           Rcpp::NumericVector m_scaled, Rcpp::NumericVector m_abs, Rcpp::NumericVector m_w,
                           Rcpp::NumericVector i_scaled, Rcpp::NumericVector i_abs, Rcpp::NumericVector i_w,
                           double gamma, int threads) {
  const int p = score.nrow(), S = score.ncol(), Km = m_w.size(), Ki = i_w.size();
  if (info.nrow() != p || info.ncol() != S * S || base.size() != p) Rcpp::stop("score / info / base do not match");
  if (m_scaled.size() != Km || m_abs.size() != Km || i_scaled.size() != Ki || i_abs.size() != Ki)
    Rcpp::stop("component vectors do not match");
  const int T = plink4r::resolve_threads(threads);
  Rcpp::NumericMatrix mean(p, S), sd(p, S);
  std::fill(mean.begin(), mean.end(), NA_REAL); std::fill(sd.begin(), sd.end(), NA_REAL);
  const double* U = score.begin();
  const double* IN = info.begin();
  const double* B = base.begin();
  double* MO = mean.begin();
  double* SO = sd.begin();
  const std::vector<double> ms(m_scaled.begin(), m_scaled.end()), ma(m_abs.begin(), m_abs.end()),
      mw(m_w.begin(), m_w.end()), is(i_scaled.begin(), i_scaled.end()), ia(i_abs.begin(), i_abs.end()),
      iw(i_w.begin(), i_w.end());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) num_threads(T)
#endif
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXd u(S);
    Eigen::MatrixXd I(S, S);
    bool ok = std::isfinite(B[j]);
    for (int q = 0; q < S && ok; ++q) { u[q] = U[j + (std::size_t)q * p]; ok = std::isfinite(u[q]); }
    for (int q = 0; q < S * S && ok; ++q) { I(q % S, q / S) = IN[j + (std::size_t)q * p]; ok = std::isfinite(I(q % S, q / S)); }
    if (!ok) continue;
    I = 0.5 * (I + I.transpose());
    const double tr = I.trace();
    if (!(tr > 0)) continue;
    Eigen::MatrixXd C = I * (S / tr);
    if (gamma != 1.0) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
      Eigen::VectorXd ev = es.eigenvalues().cwiseMax(0.0).array().pow(gamma).matrix();
      const double t2 = ev.sum();
      if (t2 > 0) ev *= S / t2;
      C = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
    }
    const int K = Km * Ki;
    std::vector<double> lw(K);
    std::vector<Eigen::VectorXd> mk(K);
    std::vector<Eigen::VectorXd> vk(K);
    double mx = -INFINITY;
    Eigen::VectorXd d(S);
    const Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(S, S);
    for (int a = 0; a < Km; ++a)
      for (int b = 0; b < Ki; ++b) {
        const int k = a * Ki + b;
        d[0] = std::sqrt(std::max(0.0, ms[a] * B[j] + ma[a]));
        const double di = std::sqrt(std::max(0.0, is[b] * B[j] + ia[b]));
        for (int q = 1; q < S; ++q) d[q] = di;
        const Eigen::MatrixXd Sig = d.asDiagonal() * C * d.asDiagonal();
        Eigen::PartialPivLU<Eigen::MatrixXd> lu(Id + I * Sig);
        const Eigen::MatrixXd Vp = Sig * lu.inverse();
        mk[k] = Vp * u;
        vk[k] = Vp.diagonal();
        const double det = lu.determinant();
        lw[k] = (mw[a] > 0 && iw[b] > 0 && det > 0) ? std::log(mw[a]) + std::log(iw[b]) - 0.5 * std::log(det) +
                                                         0.5 * u.dot(mk[k])
                                                   : -INFINITY;
        if (lw[k] > mx) mx = lw[k];
      }
    if (!std::isfinite(mx)) continue;
    double tot = 0.0;
    for (int k = 0; k < K; ++k) { lw[k] = std::exp(lw[k] - mx); tot += lw[k]; }
    Eigen::VectorXd m1 = Eigen::VectorXd::Zero(S), m2 = Eigen::VectorXd::Zero(S);
    for (int k = 0; k < K; ++k) {
      const double w = lw[k] / tot;
      if (w == 0) continue;
      m1 += w * mk[k];
      m2 += w * (vk[k] + mk[k].cwiseProduct(mk[k]));
    }
    for (int q = 0; q < S; ++q) {
      MO[j + (std::size_t)q * p] = m1[q];
      SO[j + (std::size_t)q * p] = std::sqrt(std::max(0.0, m2[q] - m1[q] * m1[q]));
    }
  }
  return Rcpp::List::create(Rcpp::Named("mean") = mean, Rcpp::Named("sd") = sd);
}
