// Annotation mixture partition of polygenic score weights (top layer: functional unit built on vi_panel.h /
// vi_store.h). Model (b_hat ~ N(R beta, R / n), R = V V' on the eigen panel, never formed):
//   classes of SNP j: null, and (g, c) for each group g with A_jg > 0 (groups = annotations + residual group) and
//   size c (gamma_c sigma2, shared across groups); P(z_j = (g, c)) = A_jg w_gc / (1 + s_j), P(null) = 1 / (1 + s_j),
//   s_j = sum_gc A_jg w_gc; beta_j | (g, c) ~ N(0, gamma_c sigma2).
// Fit: mean-field coordinate ascent over SNPs (per-SNP posterior = null + one normal per size class), sigma2 by
// closed-form EM, w by an MM update of the bound on -log(1 + s_j), with a spike-and-slab lasso type penalty on each
// annotation's total weight w_g = sum_c w_gc (adaptive lambda between Exp(lam0) spike and Exp(lam1) slab), residual
// group unpenalised. Annotation-specific weights: part_jg = sum_c q(z_j = (g, c)) mu_jc, summing to the posterior
// mean over groups. Ported from a numpy prototype (mixpart.py) to the eigen panel.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "resources.h"
#include "vi_panel.h"

using plink4r::EigenPanel;
using plink4r::RowStore;

namespace {

constexpr int NC = 3;  // size classes per group

struct Ctl {
  std::vector<double> gam;  // size-class variances relative to sigma2 (NC)
  double w_init, h2_init, lam0, lam1, theta, tol;
  int warm, max_iter, mm_iter;
  bool penalty;
};

template <typename T>
Rcpp::List fit(const EigenPanel<T>& P, const std::vector<float>& n, const std::vector<double>& Ae, int G, int Kann,
               const Ctl& c, int threads) {
  const int64_t M = P.M;
  std::vector<double> w(G * NC, c.w_init / NC);
  double gsum = 0;
  for (double x : c.gam) gsum += x;
  double sigma2 = c.h2_init / (M * c.w_init * gsum / NC);
  std::vector<float> beta(M, 0.f), s(P.S, 0.f);
  std::vector<double> gm(M * NC, 0.), mu(M * NC, 0.), sv(M * NC, 0.), Sc(M * NC, 0.);
  std::vector<float> bmax(P.nb, 0.f);
  int it = 0;
  bool conv = false;
  for (it = 0; it < c.max_iter; ++it) {
    // prior odds of each size class against null: Sc_jc = sum_g A_jg w_gc
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (int64_t j = 0; j < M; ++j)
      for (int cc = 0; cc < NC; ++cc) {
        double x = 0;
        for (int g = 0; g < G; ++g) x += Ae[j * G + g] * w[g * NC + cc];
        Sc[j * NC + cc] = x;
      }
    double v[NC];
    for (int cc = 0; cc < NC; ++cc) v[cc] = c.gam[cc] * sigma2;
    double bnorm = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(max : bnorm) num_threads(threads)
#endif
    for (int oi = 0; oi < P.nb; ++oi) {
      const int b = P.order[oi];
      const int k = P.k[b];
      float* sb = s.data() + P.s_off[b];
      float mx = 0.f;
      for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j) {
        const T* vj = P.row(b, j);
        const double dj = P.d[j], nj = n[j];
        const double rj = P.bt[j] - (RowStore<T>::dot(vj, sb, k) - dj * beta[j]);
        double lo[NC], m_[NC], s_[NC], top = 0.;
        for (int cc = 0; cc < NC; ++cc) {
          s_[cc] = 1. / (nj * dj + 1. / v[cc]);
          m_[cc] = nj * rj * s_[cc];
          const double sc = Sc[j * NC + cc];
          lo[cc] = sc > 0 ? std::log(sc) + .5 * std::log(s_[cc] / v[cc]) + m_[cc] * m_[cc] / (2. * s_[cc]) : -INFINITY;
          top = std::max(top, lo[cc]);
        }
        double den = std::exp(-top), bn = 0, e[NC];
        for (int cc = 0; cc < NC; ++cc) { e[cc] = std::exp(lo[cc] - top); den += e[cc]; }
        for (int cc = 0; cc < NC; ++cc) {
          const double g = e[cc] / den;
          gm[j * NC + cc] = g; mu[j * NC + cc] = m_[cc]; sv[j * NC + cc] = s_[cc];
          bn += g * m_[cc];
        }
        const float de = (float)bn - beta[j];
        if (de != 0.f) { RowStore<T>::axpy(sb, vj, de, k); beta[j] += de; }
        mx = std::max(mx, std::fabs(de));
        bnorm = std::max(bnorm, std::fabs((double)beta[j]));
      }
      bmax[b] = mx;
    }
    // sigma2 and expected counts Rgc = w_gc sum_j A_jg gm_jc / Sc_jc
    double num = 0, gs = 0;
    std::vector<double> Rgc(G * NC, 0.);
    for (int64_t j = 0; j < M; ++j)
      for (int cc = 0; cc < NC; ++cc) {
        const double g = gm[j * NC + cc];
        num += g * (mu[j * NC + cc] * mu[j * NC + cc] + sv[j * NC + cc]) / c.gam[cc];
        gs += g;
        const double t = Sc[j * NC + cc] > 0 ? g / Sc[j * NC + cc] : 0.;
        for (int gg = 0; gg < G; ++gg) Rgc[gg * NC + cc] += Ae[j * G + gg] * t;
      }
    for (int i = 0; i < G * NC; ++i) Rgc[i] *= w[i];
    if (gs > 0) sigma2 = num / gs;
    for (int r = 0; r < c.mm_iter; ++r) {
      std::vector<double> D(G, 0.);
      for (int64_t j = 0; j < M; ++j) {
        double sj = 0;
        for (int gg = 0; gg < G; ++gg) for (int cc = 0; cc < NC; ++cc) sj += Ae[j * G + gg] * w[gg * NC + cc];
        const double q = 1. / (1. + sj);
        for (int gg = 0; gg < G; ++gg) D[gg] += Ae[j * G + gg] * q;
      }
      std::vector<double> lam(G, 0.);
      if (c.penalty && it >= c.warm)
        for (int gg = 0; gg < Kann; ++gg) {
          double wg = 0;
          for (int cc = 0; cc < NC; ++cc) wg += w[gg * NC + cc];
          const double ps = 1. / (1. + (1. - c.theta) / c.theta * c.lam0 / c.lam1 *
                                           std::exp(-std::min((c.lam0 - c.lam1) * wg, 700.)));
          lam[gg] = c.lam1 * ps + c.lam0 * (1. - ps);
        }
      for (int gg = 0; gg < G; ++gg)
        for (int cc = 0; cc < NC; ++cc) w[gg * NC + cc] = std::max(Rgc[gg * NC + cc] / (D[gg] + lam[gg]), 1e-12);
    }
    const double md = *std::max_element(bmax.begin(), bmax.end());
    if (it > c.warm + 5 && md < c.tol * std::max(bnorm, 1e-12)) { conv = true; ++it; break; }
  }
  // annotation-specific weights and shares
  Rcpp::NumericMatrix part(M, G);
  std::vector<double> h2post(G, 0.), h2prior(G, 0.);
  for (int64_t j = 0; j < M; ++j) {
    double sj = 0;
    for (int gg = 0; gg < G; ++gg) for (int cc = 0; cc < NC; ++cc) sj += Ae[j * G + gg] * w[gg * NC + cc];
    for (int gg = 0; gg < G; ++gg) {
      const double a = Ae[j * G + gg];
      if (a == 0.) continue;
      double p = 0;
      for (int cc = 0; cc < NC; ++cc) {
        const double sc = Sc[j * NC + cc];
        if (sc <= 0) continue;
        const double resp = gm[j * NC + cc] * a * w[gg * NC + cc] / sc;  // q(z_j = (g, c))
        p += resp * mu[j * NC + cc];
        h2post[gg] += resp * (mu[j * NC + cc] * mu[j * NC + cc] + sv[j * NC + cc]);
        h2prior[gg] += sigma2 * c.gam[cc] * a * w[gg * NC + cc] / (1. + sj);
      }
      part(j, gg) = p;
    }
  }
  Rcpp::NumericVector bt(M);
  for (int64_t j = 0; j < M; ++j) bt[j] = beta[j];
  Rcpp::NumericMatrix W(G, NC);
  for (int gg = 0; gg < G; ++gg) for (int cc = 0; cc < NC; ++cc) W(gg, cc) = w[gg * NC + cc];
  return Rcpp::List::create(Rcpp::_["BETA_STD"] = bt, Rcpp::_["PART_STD"] = part,
                            Rcpp::_["h2_post"] = Rcpp::NumericVector(h2post.begin(), h2post.end()),
                            Rcpp::_["h2_prior"] = Rcpp::NumericVector(h2prior.begin(), h2prior.end()),
                            Rcpp::_["w"] = W, Rcpp::_["sigma2"] = sigma2, Rcpp::_["iters"] = it,
                            Rcpp::_["converged"] = conv);
}

}  // namespace

// Entry: blocks, keep (0-based eigen rows with data per block), z / n per kept SNP (blocks concatenated), A (kept SNPs
// x annotations, >= 0), resid_all (residual group on every SNP, else only on SNPs without annotation), thresh, fp16,
// control (mixpart_control()).
// [[Rcpp::export]]
Rcpp::List mixpart_cpp(const std::vector<std::string>& files, Rcpp::IntegerVector blocks, Rcpp::List keep, Rcpp::NumericVector z,
                       Rcpp::NumericVector n, Rcpp::NumericMatrix A, bool resid_all, double thresh, bool fp16,
                       Rcpp::List control, int threads) {
  const int T = plink4r::resolve_threads(threads);
  std::vector<int> bl = Rcpp::as<std::vector<int>>(blocks);
  std::vector<std::vector<int>> kp;
  std::size_t M = 0;
  for (int b = 0; b < keep.size(); ++b) { kp.push_back(Rcpp::as<std::vector<int>>(keep[b])); M += kp.back().size(); }
  if ((std::size_t)z.size() != M || (std::size_t)n.size() != M || (std::size_t)A.nrow() != M)
    Rcpp::stop("z, n and A need one entry / row per kept SNP");
  std::vector<float> bstd(M), nf(M);
  for (std::size_t j = 0; j < M; ++j) {
    bstd[j] = (float)(z[j] / std::sqrt(n[j] - 1. + z[j] * z[j]));
    nf[j] = (float)n[j];
  }
  const int Kann = A.ncol(), G = Kann + 1;
  std::vector<double> Ae(M * G, 0.);
  for (std::size_t j = 0; j < M; ++j) {
    bool any = false;
    for (int k = 0; k < Kann; ++k) {
      const double x = A(j, k);
      if (!(x >= 0) || !std::isfinite(x)) Rcpp::stop("annotations must be finite and >= 0");
      Ae[j * G + k] = x;
      any = any || x != 0;
    }
    Ae[j * G + Kann] = (resid_all || !any) ? 1. : 0.;
  }
  Ctl c;
  c.gam = Rcpp::as<std::vector<double>>(control["gamma"]);
  if ((int)c.gam.size() != NC) Rcpp::stop("gamma needs 3 values");
  c.w_init = Rcpp::as<double>(control["w_init"]);
  c.h2_init = Rcpp::as<double>(control["h2_init"]);
  c.lam0 = Rcpp::as<double>(control["lam0_per_snp"]) * (double)M;
  c.lam1 = Rcpp::as<double>(control["lam1"]);
  c.theta = Rcpp::as<double>(control["theta"]);
  c.tol = Rcpp::as<double>(control["tol"]);
  c.warm = Rcpp::as<int>(control["warm"]);
  c.max_iter = Rcpp::as<int>(control["max_iter"]);
  c.mm_iter = Rcpp::as<int>(control["mm_iter"]);
  c.penalty = Rcpp::as<bool>(control["penalty"]);
  try {
    if (fp16) {
      const EigenPanel<uint16_t> P = plink4r::load_eigen_panel<uint16_t>(files, bl, kp, bstd, thresh, T);
      return fit(P, nf, Ae, G, Kann, c, T);
    }
    const EigenPanel<float> P = plink4r::load_eigen_panel<float>(files, bl, kp, bstd, thresh, T);
    return fit(P, nf, Ae, G, Kann, c, T);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}
