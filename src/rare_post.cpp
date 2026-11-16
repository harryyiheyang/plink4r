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

// int_cauchy_cpp: multi-exposure interaction P of every variant. Single-exposure P values (p1, p x q) carry half
// of the weight, equally; the variance-component test carries the other half: Q = z'z, z_j = b_j / se_j, whose
// null law is sum_k lambda_k chi^2_1 with lambda the eigenvalues of the correlation of the q estimates; its upper
// tail by saddlepoint (Lugannani-Rice; K(t) = -sum log(1 - 2 lambda t) / 2). Cauchy combination
// T = sum w_i tan((1/2 - p_i) pi), P = 1/2 - atan(T) / pi (tan(.) ~ 1 / (p pi) for p < 1e-15).
namespace {

double wchisq_spa(double Q, const Eigen::VectorXd& lam) {
  const double mu = lam.sum(), var = 2.0 * lam.squaredNorm(), lmax = lam.maxCoeff();
  if (!(lmax > 0)) return NA_REAL;
  auto K = [&](double t, double& k0, double& k1, double& k2) {
    k0 = k1 = k2 = 0.0;
    for (int i = 0; i < lam.size(); ++i) {
      const double a = 1.0 - 2.0 * lam[i] * t;
      k0 -= 0.5 * std::log(a); k1 += lam[i] / a; k2 += 2.0 * lam[i] * lam[i] / (a * a);
    }
  };
  if (std::fabs(Q - mu) < 1e-3 * std::sqrt(var)) return 0.5 * std::erfc((Q - mu) / std::sqrt(2.0 * var));
  double lo = Q > mu ? 0.0 : -1e300, hi = Q > mu ? 0.5 / lmax : 0.0;
  if (Q <= mu) lo = -(1.0 / (2.0 * lmax)) * 1e8;
  double t = Q > mu ? 0.5 * hi : -1.0, k0, k1, k2;
  for (int it = 0; it < 200; ++it) {
    K(t, k0, k1, k2);
    const double f = k1 - Q;
    if (f > 0) hi = t; else lo = t;
    if (std::fabs(f) <= 1e-12 * std::max(1.0, Q)) break;
    double tn = t - f / k2;
    if (!(tn > lo && tn < hi)) tn = 0.5 * (lo + hi);
    t = tn;
  }
  K(t, k0, k1, k2);
  const double w2 = 2.0 * (t * Q - k0);
  if (!(w2 > 0)) return 0.5 * std::erfc((Q - mu) / std::sqrt(2.0 * var));
  const double w = (t > 0 ? 1.0 : -1.0) * std::sqrt(w2), u = t * std::sqrt(k2);
  const double p = 0.5 * std::erfc(w / std::sqrt(2.0)) + std::exp(-0.5 * w * w) * 0.3989422804014327 * (1.0 / u - 1.0 / w);
  return std::min(1.0, std::max(0.0, p));
}

double cauchy_term(double p) {
  if (p < 1e-15) return 1.0 / (p * M_PI);
  return std::tan((0.5 - p) * M_PI);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector int_cauchy_cpp(Rcpp::NumericMatrix b, Rcpp::NumericMatrix v, Rcpp::NumericMatrix p1, int threads) {
  const int p = b.nrow(), q = b.ncol();
  if (v.nrow() != p || v.ncol() != q * q || p1.nrow() != p || p1.ncol() != q) Rcpp::stop("b / v / p1 do not match");
  const int T = plink4r::resolve_threads(threads);
  Rcpp::NumericVector out(p, NA_REAL);
  const double* B = b.begin();
  const double* V = v.begin();
  const double* P1 = p1.begin();
  double* O = out.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1024) num_threads(T)
#endif
  for (int j = 0; j < p; ++j) {
    Eigen::MatrixXd C(q, q);
    Eigen::VectorXd z(q);
    bool ok = true;
    for (int a = 0; a < q * q && ok; ++a) { C(a % q, a / q) = V[j + (std::size_t)a * p]; ok = std::isfinite(C(a % q, a / q)); }
    for (int a = 0; a < q && ok; ++a) {
      ok = std::isfinite(B[j + (std::size_t)a * p]) && C(a, a) > 0 && std::isfinite(P1[j + (std::size_t)a * p]);
      if (ok) z[a] = B[j + (std::size_t)a * p] / std::sqrt(C(a, a));
    }
    if (!ok) continue;
    const Eigen::VectorXd d = C.diagonal().cwiseSqrt().cwiseInverse();
    const Eigen::MatrixXd R = d.asDiagonal() * (0.5 * (C + C.transpose())) * d.asDiagonal();
    const Eigen::VectorXd lam = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(R, Eigen::EigenvaluesOnly)
                                    .eigenvalues().cwiseMax(0.0);
    const double pvc = wchisq_spa(z.squaredNorm(), lam);
    if (ISNAN(pvc)) continue;
    double Tst = 0.5 * cauchy_term(std::max(pvc, 1e-300));
    for (int a = 0; a < q; ++a) Tst += 0.5 / q * cauchy_term(std::max(P1[j + (std::size_t)a * p], 1e-300));
    const double pc = Tst > 1e15 ? 1.0 / (Tst * M_PI) : 0.5 - std::atan(Tst) / M_PI;
    O[j] = std::min(1.0, std::max(pc, 0.0));
  }
  return out;
}
