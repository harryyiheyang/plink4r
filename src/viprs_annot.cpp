// VIPRS with an additive annotation model of the per-SNP effect variance, and the partition of the effects into
// annotation-specific components (PartitionPRS, Research Design D.1.2). Top layer: the package's annotation VIPRS
// functional unit, built on vi_panel.h / vi_store.h / ld_bins.h; copied from viprs.cpp and edited.
// Model: beta_j = sum_k beta_j^(k), beta_j^(k) independent with variance A_jk sigma2_k (A_jk >= 0; k = 0 is the
// residual annotation A_R: 1 for every SNP, or only for SNPs without annotation), so h2_j = Var(beta_j) = sum_k A_jk sigma2_k. Prior on beta_j: VIPRS
// spike and slab with fixed pi; the slab is a mixture of C classes with fixed shares p_c (mix_prop, default small /
// medium / large = 5 : 3 : 2) and variances kap_c h2_j / pi (mix_var relative, normalized to sum_c p_c kap_c = 1),
// so Var(beta_j) = h2_j. Coordinate update (D_j = n_j r_j / sigma; tau_j = pi / h2_j):
//   vt_jc = n_j d_j / sigma + tau_j / kap_c,  u_c = log p_c + .5 D_j^2 / vt_jc + .5 (log tau_j - log kap_c - log vt_jc)
//   gamma_j = sigmoid(logit(pi) + logsumexp_c u_c), psi_jc = softmax(u)_c, eta_j = gamma_j sum_c psi_jc D_j / vt_jc
// (one class: the VIPRS update with tau_j in place of tau).
// M-step (D.1.2): w_jk = A_jk sigma2_k / h2_j, E[(beta_j^(k))^2] = w_jk^2 zeta_j + w_jk (1 - w_jk) h2_j with
// zeta_j = gamma_j sum_c psi_jc (mu_jc^2 + 1 / vt_jc) / kap_c (the second moment on the scale of h2_j: the EM
// statistic of the scale of a Gaussian scale mixture), sigma2_k <- (S_k + 2 b0) / (sum_j A_jk + 2 (a0 + 1)),
// S_k = sum_j E[(beta_j^(k))^2]. With one slab class, the residual annotation only and a0 = -1, b0 = 0 this is the
// VIPRS 0.1.5 tau update (tau = pi M / sum zeta), i.e. viprs(). Start: sigma2_init (e.g. ldsc_strat() tau, cut at 0)
// or all sigma2_k equal with sum_j h2_j = the naive LDSC h2.
// Output: posterior mean eta_j and its MMSE allocation eta_j^(k) = w_jk eta_j (sum_k eta_j^(k) = eta_j), the
// weights of the annotation-specific PRS.
// LD-bin suppression (as sbayesrc()): inside an LD bin (ld_bins.h) with at least one annotated member (any nonzero
// binary 0/1 annotation column), the members without annotation are fixed at gamma_j = 0 (beta_j = 0, never updated,
// no residual update) and left out of the M-step.
// Hyperparameter search as viprs.cpp (pi grid ascending with warm start, PUMAS, best converged point refitted); the
// sigma2_k (variance of beta) are carried over when pi moves (the start sets those of the first grid point).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "resources.h"
#include "ld_bins.h"
#include "vi_panel.h"

using plink4r::EigenPanel;
using plink4r::RowStore;

namespace {

struct ViprsControl {
  bool use_active;
  int full_every, refresh_every, min_iter, patience, max_iter, grid_max_iter, stop_after, pi_steps, grid_extend;
  double act_tol, x_tol, f_tol, prop_train, sigma_fixed, h2_min, h2_max, pi_lo_count, pi_lo_min, pi_hi_count,
      pi_hi_max, pi_cap, a0, b0, init_floor;
  unsigned long long seed;
  std::vector<double> pi_grid_user;  // empty: official grid
  // slab mixture: class c has prior share mix_prop[c] and variance kap[c] h2_j / pi (sum_c mix_prop[c] kap[c] = 1)
  std::vector<float> lprop, kap, lkap;
};

ViprsControl read_control(const Rcpp::List& c) {
  ViprsControl x;
  x.use_active = Rcpp::as<bool>(c["use_active"]);
  x.full_every = Rcpp::as<int>(c["full_every"]);
  x.refresh_every = Rcpp::as<int>(c["refresh_every"]);
  x.min_iter = Rcpp::as<int>(c["min_iter"]);
  x.patience = Rcpp::as<int>(c["patience"]);
  x.max_iter = Rcpp::as<int>(c["max_iter"]);
  x.grid_max_iter = Rcpp::as<int>(c["grid_max_iter"]);
  x.stop_after = Rcpp::as<int>(c["stop_after"]);
  x.pi_steps = Rcpp::as<int>(c["pi_steps"]);
  x.grid_extend = Rcpp::as<int>(c["grid_extend"]);
  x.act_tol = Rcpp::as<double>(c["act_tol"]);
  x.x_tol = Rcpp::as<double>(c["x_tol"]);
  x.f_tol = Rcpp::as<double>(c["f_tol"]);
  x.prop_train = Rcpp::as<double>(c["prop_train"]);
  x.sigma_fixed = Rcpp::as<double>(c["sigma_fixed"]);
  x.h2_min = Rcpp::as<double>(c["h2_min"]);
  x.h2_max = Rcpp::as<double>(c["h2_max"]);
  x.pi_lo_count = Rcpp::as<double>(c["pi_lo_count"]);
  x.pi_lo_min = Rcpp::as<double>(c["pi_lo_min"]);
  x.pi_hi_count = Rcpp::as<double>(c["pi_hi_count"]);
  x.pi_hi_max = Rcpp::as<double>(c["pi_hi_max"]);
  x.pi_cap = Rcpp::as<double>(c["pi_cap"]);
  x.a0 = Rcpp::as<double>(c["a0"]);
  x.b0 = Rcpp::as<double>(c["b0"]);
  x.init_floor = Rcpp::as<double>(c["init_floor"]);
  x.seed = static_cast<unsigned long long>(Rcpp::as<double>(c["seed"]));
  if (!Rf_isNull(c["pi_grid"])) x.pi_grid_user = Rcpp::as<std::vector<double>>(c["pi_grid"]);
  std::sort(x.pi_grid_user.begin(), x.pi_grid_user.end());
  const std::vector<double> mp = Rcpp::as<std::vector<double>>(c["mix_prop"]), mv = Rcpp::as<std::vector<double>>(c["mix_var"]);
  if (mp.empty() || mp.size() != mv.size() || mp.size() > 16)
    Rcpp::stop("mix_prop and mix_var need the same length, 1 to 16");
  double sp = 0, spv = 0;
  for (std::size_t i = 0; i < mp.size(); ++i) {
    if (!(mp[i] > 0) || !(mv[i] > 0)) Rcpp::stop("mix_prop and mix_var must be > 0");
    sp += mp[i];
  }
  for (std::size_t i = 0; i < mp.size(); ++i) spv += mp[i] / sp * mv[i];
  for (std::size_t i = 0; i < mp.size(); ++i) {
    x.lprop.push_back((float)std::log(mp[i] / sp));
    x.kap.push_back((float)(mv[i] / spv));
    x.lkap.push_back((float)std::log(mv[i] / spv));
  }
  return x;
}

struct State {
  std::vector<float> gamma, mu, eta, s, r;
  std::vector<uint8_t> active;
  State(int64_t M, int64_t S, double pi)
      : gamma(M, (float)pi), mu(M, 0.f), eta(M, 0.f), s(S, 0.f), r(M, 0.f), active(M, 1) {}
};

// ---- kernels on the panel

// Posterior of one SNP under the slab mixture, from the data term D = n r / sigma and base = n d / sigma:
// class c: vt_c = base + tau / kap_c, u_c = log prop_c + .5 D^2 / vt_c + .5 (log tau - log kap_c - log vt_c);
// logit gamma = lp + logsumexp(u), psi_c = softmax(u), slab mean sum_c psi_c D / vt_c. With one class this is the
// VIPRS update. psi (C entries) and vt are written when non-null.
inline float mix_post(float D, float base, float tau, float ltau, float lp, const ViprsControl& c, float& slab_mean,
                      float* psi, float* vt) {
  const int C = (int)c.kap.size();
  float u[16], v[16], umax = -INFINITY;
  for (int i = 0; i < C; ++i) {
    v[i] = base + tau / c.kap[i];
    u[i] = c.lprop[i] + .5f * D * D / v[i] + .5f * (ltau - c.lkap[i] - std::log(v[i]));
    umax = std::max(umax, u[i]);
  }
  float se = 0.f, m = 0.f;
  for (int i = 0; i < C; ++i) { u[i] = std::exp(u[i] - umax); se += u[i]; m += u[i] * D / v[i]; }
  slab_mean = m / se;
  if (psi) for (int i = 0; i < C; ++i) { psi[i] = u[i] / se; vt[i] = v[i]; }
  const float x = lp + umax + std::log(se);
  return x < 0 ? std::exp(x) / (1.f + std::exp(x)) : 1.f / (1.f + std::exp(-x));
}

// st.mu holds the data term D_j = n_j r_j / sigma of the last update (the class means are D_j / vt_jc)
template <typename T>
void sweep(const EigenPanel<T>& P, const float* n, const float* bt, float pi, const float* tau, const float* ltau,
           const char* sup, const ViprsControl& c,
           float sigma, State& st, bool active_only, float act_tol, std::vector<float>& bmax, std::vector<int64_t>& nvis, int threads) {
  const float eps = std::max(std::numeric_limits<float>::epsilon(), 1e-8f);
  const float lp = std::log(pi) - std::log1p(-pi);
  float* gamma = st.gamma.data(); float* mu = st.mu.data(); float* eta = st.eta.data();
  float* r = st.r.data(); uint8_t* active = st.active.data();
  const float* d = P.d.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int oi = 0; oi < P.nb; ++oi) {
    const int b = P.order[oi];
    const int k = P.k[b];
    float* sb = st.s.data() + P.s_off[b];
    float mx = 0.f;
    int64_t cnt = 0;
    for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j) {
      if (sup[j]) continue;  // suppressed: gamma = eta = 0 throughout
      const T* vj = P.row(b, j);
      const bool lazy = active_only && !active[j];
      if (!lazy) {
        ++cnt;
        r[j] = bt[j] - (RowStore<T>::dot(vj, sb, k) - d[j] * eta[j]);
      }
      const float D = n[j] * r[j] / sigma;
      float m_;
      const float g = mix_post(D, n[j] * d[j] / sigma, tau[j], ltau[j], lp, c, m_, nullptr, nullptr);
      const float de = g * m_ - eta[j];
      const float ad = std::fabs(de);
      if (ad >= eps) {
        mu[j] = D;
        gamma[j] = g;
        if (!lazy || ad > act_tol) {
          RowStore<T>::axpy(sb, vj, de, k);
          eta[j] += de;
          if (lazy) active[j] = 1;
        }
      }
      if (ad > mx) mx = ad;
      if (!active_only) active[j] = ad > act_tol;
    }
    bmax[b] = mx;
    nvis[b] = cnt;
  }
}

// s = V' x per block (double accumulation)
template <typename T>
void vtx(const EigenPanel<T>& P, const float* x, float* s, int threads) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int b = 0; b < P.nb; ++b) {
    const int k = P.k[b];
    std::vector<double> acc(k, 0.);
    for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j) {
      const double xj = x[j];
      if (xj == 0.) continue;
      const T* vj = P.row(b, j);
      for (int i = 0; i < k; ++i) acc[i] += xj * RowStore<T>::get(vj, i);
    }
    for (int i = 0; i < k; ++i) s[P.s_off[b] + i] = (float)acc[i];
  }
}

// y = V z per block (z in eigen coordinates)
template <typename T>
void vz(const EigenPanel<T>& P, const float* z, float* y, int threads) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int b = 0; b < P.nb; ++b)
    for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j) y[j] = RowStore<T>::dot(P.row(b, j), z + P.s_off[b], P.k[b]);
}

// Sums for the M-step and the ELBO with (tau, sigma) of the last sweep, the class means D_j / vt_jc:
// [0] sum gamma, [1] sum zeta, [2] sum d E[beta^2], [3] eta'bt, [4] sum d eta^2, [5] s's, [6] sum gamma log gamma,
// [7] sum (1-gamma) log(1-gamma), [8] sum gamma sum_c psi_c (log vt_c + log kap_c - 2 log prop_c + 2 log psi_c),
// [9] eta'bv (bv may be null). zeta_j = gamma sum_c psi_c (mu_c^2 + 1/vt_c) / kap_c (E[beta^2] on the scale of
// h2_j, prior mean h2_j; the M-step statistic), also returned per SNP. One class: zeta = E[beta^2], [8] = sum gamma
// log vt (VIPRS).
void stats(int64_t M, const float* d, const float* n, const float* bt, const float* bv, const State& st, int64_t S,
           const float* tau, const float* ltau, const ViprsControl& c, float sigma, double* out, double* zeta,
           int threads) {
  double a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0, a8 = 0, a9 = 0;
  const double lo = 1e-12;
  const int C = (int)c.kap.size();
  const float* gamma = st.gamma.data(); const float* mu = st.mu.data(); const float* eta = st.eta.data();
  const float* s = st.s.data();
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : a0, a1, a2, a3, a4, a6, a7, a8, a9) num_threads(threads)
#endif
  for (int64_t j = 0; j < M; ++j) {
    float psi[16], vt[16], m_;
    const float D = mu[j];
    mix_post(D, n[j] * d[j] / sigma, tau[j], ltau[j], 0.f, c, m_, psi, vt);
    double e2 = 0, zs = 0, lv = 0;
    for (int i = 0; i < C; ++i) {
      const double mc = (double)D / vt[i], q = (double)mc * mc + 1. / vt[i];
      e2 += psi[i] * q;
      zs += psi[i] * q / c.kap[i];
      lv += psi[i] * (std::log((double)vt[i]) + c.lkap[i] - 2. * c.lprop[i] +
                      (psi[i] > 0 ? 2. * std::log((double)psi[i]) : 0.));
    }
    const double g = std::min(std::max((double)gamma[j], lo), 1. - lo);
    const double z = gamma[j] * zs;
    a0 += gamma[j];
    a1 += z;
    zeta[j] = z;
    a2 += d[j] * gamma[j] * e2;
    a3 += (double)eta[j] * bt[j];
    a4 += (double)d[j] * eta[j] * eta[j];
    a6 += g * std::log(g);
    a7 += (1. - g) * std::log(1. - g);
    a8 += gamma[j] * lv;
    if (bv) a9 += (double)eta[j] * bv[j];
  }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : a5) num_threads(threads)
#endif
  for (int64_t i = 0; i < S; ++i) a5 += (double)s[i] * s[i];
  out[0] = a0; out[1] = a1; out[2] = a2; out[3] = a3; out[4] = a4; out[5] = a5;
  out[6] = a6; out[7] = a7; out[8] = a8; out[9] = a9;
}

// ---- variance model: h2_j = sum_k a_jk sigma2_k, a_j = (1, A_j) (A row-major M x (K - 1), A >= 0)

struct VarModel {
  int64_t M;
  int K;                         // annotations, residual annotation first
  const std::vector<double>& A;
  const std::vector<char>& sup;  // suppressed SNPs: left out of the M-step
  std::vector<double> aR;        // residual annotation: 1 for every SNP (resid_all) or only for SNPs without annotation
  std::vector<double> par;       // sigma2_k
  std::vector<double> h2;        // h2_j
  std::vector<float> tau, ltau;  // tau_j = pi / h2_j and log tau_j, as used by the sweep
  double pi = 1;
  int threads;

  VarModel(int64_t M_, int K_, const std::vector<double>& A_, const std::vector<char>& sup_, int threads_)
      : M(M_), K(K_), A(A_), sup(sup_), aR(M_, 1.), par(K_, 0.), h2(M_, 1.), tau(M_, 1.f), ltau(M_, 0.f),
        threads(threads_) {}
  void set_residual(bool all) {
    for (int64_t j = 0; j < M; ++j) {
      bool any = false;
      const double* r = row(j);
      for (int k = 1; k < K && !any; ++k) any = r[k - 1] != 0;
      aR[j] = (all || !any) ? 1. : 0.;
    }
  }

  const double* row(int64_t j) const { return A.data() + j * (K - 1); }
  void set_all() {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (int64_t j = 0; j < M; ++j) {
      double h = aR[j] * par[0];
      const double* r = row(j);
      for (int k = 1; k < K; ++k) h += r[k - 1] * par[k];
      h2[j] = h;
      const double t = pi / h;
      tau[j] = (float)t; ltau[j] = (float)std::log(t);
    }
  }
  // start: all sigma2_k equal, sum_j h2_j = h2 (LDSC start)
  void start(double h2_0, double pi_) {
    pi = pi_;
    double s = 0;
    for (int64_t j = 0; j < M; ++j) {
      s += aR[j];
      const double* r = row(j);
      for (int k = 1; k < K; ++k) s += r[k - 1];
    }
    std::fill(par.begin(), par.end(), h2_0 / s);
    set_all();
  }
  // start from given sigma2_k (e.g. stratified LDSC tau_k): values below floor * (mean h2_j of the values cut at 0)
  // are raised to it, so no sigma2_k starts at 0 (0 is a fixed point of the M-step); returns sum_j h2_j over the
  // SNPs not suppressed, 0 when the values cut at 0 give no positive h2_j (then nothing is set)
  double start_at(const std::vector<double>& s2, double floor, double pi_) {
    std::vector<double> p(K);
    for (int k = 0; k < K; ++k) p[k] = std::max(s2[k], 0.);
    auto total = [&](const std::vector<double>& q) {
      double tot = 0;
      for (int64_t j = 0; j < M; ++j) {
        if (sup[j]) continue;
        double h = aR[j] * q[0];
        const double* r = row(j);
        for (int k = 1; k < K; ++k) h += r[k - 1] * q[k];
        tot += h;
      }
      return tot;
    };
    int64_t m = 0;
    for (int64_t j = 0; j < M; ++j) m += !sup[j];
    const double t0 = total(p);
    if (!(t0 > 0) || m == 0) return 0.;
    const double lo = floor * t0 / m;
    for (int k = 0; k < K; ++k) p[k] = std::max(p[k], lo);
    par = p; pi = pi_;
    set_all();
    return total(p);
  }
  // pi moves: sigma2_k (variance of beta) kept, the slab variance h2_j / pi follows
  void pi_move(double pi_) { pi = pi_; set_all(); }

  // one update of sigma2_k (D.1.2); returns nothing, refreshes h2, tau
  void mstep(const double* zeta, double a0, double b0) {
    const int K1 = K - 1;
    std::vector<double> SN(2 * K, 0.);  // S_k, then sum_j A_jk
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
    {
      std::vector<double> loc(2 * K, 0.);
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (int64_t j = 0; j < M; ++j) {
        if (sup[j]) continue;
        const double h = h2[j], z = zeta[j];
        const double* r = row(j);
        for (int k = 0; k < K; ++k) {
          const double a = k ? r[k - 1] : aR[j];
          if (a == 0.) continue;
          const double w = a * par[k] / h;
          loc[k] += w * w * z + w * (1. - w) * h;
          loc[K + k] += a;
        }
      }
#ifdef _OPENMP
#pragma omp critical
#endif
      for (int i = 0; i < 2 * K; ++i) SN[i] += loc[i];
    }
    for (int k = 0; k < K; ++k) {
      const double den = SN[K + k] + 2. * (a0 + 1.);
      par[k] = den > 0 ? std::max((SN[k] + 2. * b0) / den, 0.) : 0.;
    }
    (void)K1;
    set_all();
  }
};

// ---- one EM run with some hyperparameters fixed

struct RunResult {
  double pi, sigma, h2, elbo;
  int iters;
  bool converged, sigma_fixed;
  std::string message;
};

template <typename T>
struct Fitter {
  const EigenPanel<T>& P;
  const ViprsControl& c;
  const std::vector<float>& n;  // per-SNP sample size
  std::vector<float> z;         // per-SNP z (LDSC start)
  int threads;
  double n_max;
  int64_t sweeps_full = 0, sweeps_part = 0, snp_visits = 0;

  Fitter(const EigenPanel<T>& P_, const ViprsControl& c_, const std::vector<float>& n_, const std::vector<float>& z_,
         int threads_)
      : P(P_), c(c_), n(n_), z(z_), threads(threads_) {
    n_max = *std::max_element(n.begin(), n.end());
  }

  // VIPRS-style naive LDSC: h2 = M (mean chi2 - 1) / (median N mean l), clipped to [h2_min, h2_max]
  double h2_ldsc() const {
    double chi = 0, l = 0;
    for (int64_t j = 0; j < P.M; ++j) { chi += (double)z[j] * z[j]; l += P.ldscore[j]; }
    std::vector<float> nn(n);
    const std::size_t h = nn.size() / 2;
    std::nth_element(nn.begin(), nn.begin() + h, nn.end());
    double med = nn[h];
    if (nn.size() % 2 == 0) med = (med + *std::max_element(nn.begin(), nn.begin() + h)) / 2;
    const double h2 = P.M * (chi / P.M - 1.) / (med * (l / P.M));
    return std::min(std::max(h2, c.h2_min), c.h2_max);
  }

  // official VIPRS GS grid: pi_steps values, log10 spaced from max(pi_lo_count / M, pi_lo_min) to
  // min(pi_hi_count / M, pi_hi_max), extended upward by grid_extend points of the same spacing (capped at pi_cap);
  // or the user grid (sorted)
  std::vector<double> pi_grid() const {
    if (!c.pi_grid_user.empty()) return c.pi_grid_user;
    const double M = (double)P.M, steps = c.pi_steps;
    const double a = std::log10(std::max(c.pi_lo_count / M, c.pi_lo_min)), e = std::log10(std::min(c.pi_hi_count / M, c.pi_hi_max));
    std::vector<double> g(c.pi_steps);
    for (int i = 0; i < c.pi_steps; ++i) g[i] = std::pow(10., c.pi_steps > 1 ? a + i * (e - a) / (steps - 1) : a);
    if (c.pi_steps > 1) g.back() = std::pow(10., e);
    if (c.grid_extend && c.pi_steps > 1) {
      const double r = g[1] / g[0], last = g.back();
      for (int i = 1; i <= c.grid_extend; ++i) {
        const double x = last * std::pow(r, (double)i);
        if (x <= c.pi_cap) g.push_back(x);
      }
    }
    return g;
  }

  // EM at fixed pi: E-step sweeps, M-step of the variance model (vm, updated in place) and of sigma
  RunResult run(State& st, const float* nv, const float* bt, double pi, VarModel& vm, double sigma, int max_iter,
                const float* bv) {
    const int64_t M = P.M;
    std::vector<float> bmax(P.nb, 0.f);
    std::vector<int64_t> nvis(P.nb, 0);
    double out[10];
    std::vector<double> zeta(M);
    vtx(P, st.eta.data(), st.s.data(), threads);
    std::fill(st.active.begin(), st.active.end(), 1);
    double prev_elbo = -INFINITY, prev_sg = NAN, elbo = NAN, sigma_g = 0;
    int sg_count = 0, it = 0;
    bool have_prev_sg = false, force_full = true, fixed_sigma = false, ok = false;
    std::string msg = "max_iter reached";
    for (it = 1; it <= max_iter; ++it) {
      const bool full = !c.use_active || force_full || (it % c.full_every == 0);
      sweep(P, nv, bt, (float)pi, vm.tau.data(), vm.ltau.data(), vm.sup.data(), c, (float)sigma, st, !full, (float)c.act_tol, bmax, nvis, threads);
      for (int b = 0; b < P.nb; ++b) snp_visits += nvis[b];
      if (full) ++sweeps_full; else ++sweeps_part;
      force_full = false;
      if (it % c.refresh_every == 0) vtx(P, st.eta.data(), st.s.data(), threads);
      stats(M, P.d.data(), nv, bt, bv, st, P.S, vm.tau.data(), vm.ltau.data(), c, (float)sigma, out, zeta.data(), threads);
      const double sg_sum = out[0], sdz = out[2], eb = out[3], sde2 = out[4], ss = out[5], glg = out[6],
                   g1lg1 = out[7], glvt = out[8];
      vm.mstep(zeta.data(), c.a0, c.b0);
      // prior part of the ELBO with the updated variance model (as VIPRS 0.1.5 with the updated tau)
      double glt = 0, tz = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : glt, tz) num_threads(threads)
#endif
      for (int64_t j = 0; j < M; ++j) {
        glt += st.gamma[j] * (double)vm.ltau[j];
        tz += (double)vm.tau[j] * zeta[j];
      }
      const double prior = .5 * (sg_sum - glvt + glt) - .5 * tz;
      sigma_g = sdz + (ss - sde2);
      // sigma M-step and the likelihood part of the ELBO (VIPRS.elbo: y'y / N = 1)
      if (!fixed_sigma) {
        const double ns = 1. - 2. * eb + sigma_g;
        if (ns <= 0. || !std::isfinite(ns)) { fixed_sigma = true; sigma = c.sigma_fixed; }
        else sigma = ns;
      }
      const double lik = !fixed_sigma ? -std::log(2 * M_PI * sigma) - 1.
                                      : -std::log(2 * M_PI * sigma) - (1. / sigma) * (1. - 2. * eb + sigma_g);
      elbo = .5 * n_max * lik;
      elbo += -glg + sg_sum * std::log(pi) - g1lg1 + (M - sg_sum) * std::log1p(-pi);
      elbo += prior;
      const double max_diff = *std::max_element(bmax.begin(), bmax.end());
      if (!std::isfinite(elbo)) { msg = "ELBO undefined"; break; }
      if (!full) {  // never declare convergence on a partial sweep; check with a full sweep next
        force_full = max_diff < c.x_tol;
        prev_elbo = elbo;
        continue;
      }
      if (it > c.min_iter && std::fabs(elbo - prev_elbo) <= c.f_tol) { msg = "ELBO converged"; ok = true; break; }
      if (it > c.min_iter && max_diff < c.x_tol) { msg = "parameters converged"; ok = true; break; }
      if (have_prev_sg && it > c.min_iter && std::fabs(sigma_g - prev_sg) < c.x_tol && max_diff < 10 * c.x_tol) {
        ++sg_count;
        if (sg_count > c.patience) { msg = "LD-weighted effects converged"; ok = true; break; }
      }
      prev_sg = sigma_g; have_prev_sg = true; prev_elbo = elbo;
    }
    if (it > max_iter) it = max_iter;
    vtx(P, st.eta.data(), st.s.data(), threads);
    return RunResult{pi, sigma, sigma_g / (sigma_g + sigma), elbo, it, ok, fixed_sigma, msg};
  }

  // PUMAS split (magenpy / VIPRS): training half bt + sqrt(c) V xi, xi ~ N(0, I_S), c = 1/(prop_train N) - 1/N
  void pumas(std::vector<float>& train, std::vector<float>& valid, std::vector<float>& n_tr, std::vector<double>& xi) {
    std::mt19937_64 eng(c.seed);
    std::normal_distribution<double> N01;
    xi.resize(P.S);
    for (auto& x : xi) x = N01(eng);
    const double cc = 1. / (c.prop_train * n_max) - 1. / n_max, sc = std::sqrt(cc);
    std::vector<float> xf(xi.begin(), xi.end()), y(P.M);
    vz(P, xf.data(), y.data(), threads);
    train.resize(P.M); valid.resize(P.M); n_tr.resize(P.M);
    const float pt = (float)(1. - c.prop_train), ratio = (float)(c.prop_train / (1. - c.prop_train));
    for (int64_t j = 0; j < P.M; ++j) {
      train[j] = (float)(P.bt[j] + sc * y[j]);
      valid[j] = P.bt[j] / pt - train[j] * ratio;
      n_tr[j] = n[j] * (float)c.prop_train;
    }
  }

  static void pseudo_r2(const State& st, const std::vector<float>& bv, double& r2, double& rb) {
    rb = 0;
    double bsb = 0;
    for (std::size_t j = 0; j < st.eta.size(); ++j) rb += (double)st.eta[j] * bv[j];
    for (float x : st.s) bsb += (double)x * x;
    r2 = bsb > 0 ? rb * rb / bsb : 0.;
  }
};

// ---- output

// annotation heritability h2_k = (sum_j A_jk) sigma2_k over the SNPs not suppressed (k = 0: residual annotation)
Rcpp::NumericVector h2_annot(const VarModel& vm) {
  Rcpp::NumericVector h(vm.K);
  for (int64_t j = 0; j < vm.M; ++j) {
    if (vm.sup[j]) continue;
    const double* r = vm.row(j);
    h[0] += vm.aR[j] * vm.par[0];
    for (int k = 1; k < vm.K; ++k) h[k] += r[k - 1] * vm.par[k];
  }
  return h;
}

template <typename T>
Rcpp::List pack(const Fitter<T>& F, const State& st, const RunResult& res, const VarModel& vm,
                const std::vector<double>& af, double h2_0, Rcpp::List extra) {
  const int64_t M = F.P.M;
  const int K = vm.K;
  Rcpp::NumericVector beta(M), bstd(M), pip(M), h2j(M);
  Rcpp::NumericMatrix part(M, K);  // eta_j^(k) = w_jk eta_j (standardized scale), k = 0 residual annotation
  for (int64_t j = 0; j < M; ++j) {
    const double p = std::min(std::max(af[j], 1e-6), 1 - 1e-6);
    bstd[j] = st.eta[j];
    beta[j] = st.eta[j] / std::sqrt(2. * p * (1. - p));
    pip[j] = st.gamma[j];
    h2j[j] = vm.h2[j];
    const double* r = vm.row(j);
    for (int k = 0; k < K; ++k) part(j, k) = (k ? r[k - 1] : vm.aR[j]) * vm.par[k] / vm.h2[j] * st.eta[j];
  }
  Rcpp::List hyp = Rcpp::List::create(
      Rcpp::_["pi"] = res.pi, Rcpp::_["sigma"] = res.sigma, Rcpp::_["h2"] = res.h2,
      Rcpp::_["elbo"] = res.elbo, Rcpp::_["iters"] = res.iters,
      Rcpp::_["converged"] = res.converged, Rcpp::_["message"] = res.message, Rcpp::_["sigma_fixed"] = res.sigma_fixed,
      Rcpp::_["h2_ldsc_init"] = h2_0, Rcpp::_["n_snps"] = (double)M, Rcpp::_["sum_k"] = (double)F.P.S,
      Rcpp::_["sweeps_full"] = (double)F.sweeps_full, Rcpp::_["sweeps_partial"] = (double)F.sweeps_part,
      Rcpp::_["snp_visits_per_M"] = (double)F.snp_visits / M);
  return Rcpp::List::create(Rcpp::_["BETA"] = beta, Rcpp::_["BETA_STD"] = bstd, Rcpp::_["PIP"] = pip,
                            Rcpp::_["H2_SNP"] = h2j, Rcpp::_["PART_STD"] = part,
                            Rcpp::_["sigma2"] = Rcpp::NumericVector(vm.par.begin(), vm.par.end()),
                            Rcpp::_["h2_annot"] = h2_annot(vm), Rcpp::_["hyp"] = hyp, Rcpp::_["extra"] = extra);
}

// shared set-up of the three entries: inputs -> standardized b, panel
struct Inputs {
  std::vector<int> blocks;
  std::vector<std::vector<int>> keep;
  std::vector<float> z, n, bstd;
  std::vector<double> af;
};

Inputs read_inputs(Rcpp::IntegerVector blocks, Rcpp::List keep, Rcpp::NumericVector z, Rcpp::NumericVector n,
                   Rcpp::NumericVector af) {
  Inputs in;
  in.blocks = Rcpp::as<std::vector<int>>(blocks);
  if (keep.size() != blocks.size()) Rcpp::stop("blocks and keep differ in length");
  std::size_t M = 0;
  for (int b = 0; b < keep.size(); ++b) {
    in.keep.push_back(Rcpp::as<std::vector<int>>(keep[b]));
    M += in.keep.back().size();
  }
  if ((std::size_t)z.size() != M || (std::size_t)n.size() != M || (std::size_t)af.size() != M)
    Rcpp::stop("z, n and af must have one entry per kept SNP");
  in.z.resize(M); in.n.resize(M); in.bstd.resize(M); in.af.resize(M);
  for (std::size_t j = 0; j < M; ++j) {
    if (!std::isfinite(z[j]) || !(n[j] > 1)) Rcpp::stop("z must be finite and n > 1");
    in.z[j] = (float)z[j];
    in.n[j] = (float)n[j];
    in.af[j] = af[j];
    in.bstd[j] = (float)(z[j] / std::sqrt(n[j] - 1. + z[j] * z[j]));  // marginal standardized effect
  }
  return in;
}

// pi ascending with warm start (state and variance model carried over), PUMAS pseudo-validation, early stop after
// stop_after consecutive drops of pseudo R2, refit of the best converged point on all data (points that did not
// converge are never chosen)
template <typename T>
Rcpp::List fit_path(const Inputs& in, const std::vector<double>& A, int K, const std::vector<char>& sup, double n_bins,
                    bool resid_all, const std::vector<double>& s2_init, const std::vector<std::string>& files,
                    double thresh, const ViprsControl& c, int threads) {
  const EigenPanel<T> P = plink4r::load_eigen_panel<T>(files, in.blocks, in.keep, in.bstd, thresh, threads);
  Fitter<T> F(P, c, in.n, in.z, threads);
  const double h2_0 = F.h2_ldsc();
  const std::vector<double> grid = F.pi_grid();
  std::vector<float> train, valid, n_tr;
  std::vector<double> xi;
  F.pumas(train, valid, n_tr, xi);
  std::vector<double> gpi, gsig, gh2, gr2, grb, gsec;
  std::vector<int> git, gconv;
  Rcpp::NumericMatrix gcoef(grid.size(), K);
  double best_r2 = -INFINITY, bpi = 0, bsig = 0, best_here = -INFINITY;
  State st(P.M, P.S, grid[0]), best(0, 0, 0.);
  for (int64_t j = 0; j < P.M; ++j) if (sup[j]) st.gamma[j] = 0.f;
  VarModel vm(P.M, K, A, sup, threads), bvm(P.M, K, A, sup, threads);
  vm.set_residual(resid_all); bvm.set_residual(resid_all);
  bool have_best = false;
  double h2_start = h2_0;
  if (!s2_init.empty()) {
    const double tot = vm.start_at(s2_init, c.init_floor, grid[0]);
    if (tot > 0) h2_start = std::min(std::max(tot, c.h2_min), c.h2_max);
    else Rcpp::warning("sigma2_init has no positive value: the LDSC start is used");
    if (!(tot > 0)) vm.start(h2_0, grid[0]);
  } else vm.start(h2_0, grid[0]);
  double sigma = 1. - h2_start;
  int drops = 0;
  for (double pi : grid) {
    vm.pi_move(pi);
    const auto t1 = std::chrono::steady_clock::now();
    const RunResult r = F.run(st, n_tr.data(), train.data(), pi, vm, sigma, c.grid_max_iter, valid.data());
    sigma = r.sigma;
    double r2, rb;
    Fitter<T>::pseudo_r2(st, valid, r2, rb);
    for (int k = 0; k < K; ++k) gcoef(gpi.size(), k) = vm.par[k];
    gpi.push_back(pi); gsig.push_back(sigma); gh2.push_back(r.h2); gr2.push_back(r2);
    grb.push_back(rb); git.push_back(r.iters); gconv.push_back(r.converged);
    gsec.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count());
    if (r.converged && rb > 0 && r2 > best_r2) {
      best_r2 = r2; best = st; have_best = true;
      bpi = pi; bsig = sigma; bvm.par = vm.par; bvm.pi = pi; bvm.set_all();
    }
    if (r2 > best_here) { best_here = r2; drops = 0; }
    else if (c.stop_after && ++drops >= c.stop_after) break;
  }
  if (!have_best) Rcpp::stop("no grid point converged with a positive pseudo R2");
  RunResult res = F.run(best, in.n.data(), P.bt.data(), bpi, bvm, bsig, c.max_iter, nullptr);
  const int ng = gpi.size();
  Rcpp::DataFrame tab = Rcpp::DataFrame::create(
      Rcpp::_["pi"] = gpi, Rcpp::_["sigma"] = gsig, Rcpp::_["h2"] = gh2, Rcpp::_["pseudo_r2"] = gr2,
      Rcpp::_["rb"] = grb, Rcpp::_["iters"] = git, Rcpp::_["converged"] = gconv,
      Rcpp::_["seconds"] = gsec);
  return pack(F, best, res, bvm, in.af, h2_0,
              Rcpp::List::create(Rcpp::_["grid"] = tab, Rcpp::_["grid_sigma2"] = gcoef(Rcpp::Range(0, ng - 1), Rcpp::_),
                                 Rcpp::_["best_pseudo_r2"] = best_r2, Rcpp::_["n_bins"] = n_bins,
                                 Rcpp::_["suppressed"] = Rcpp::LogicalVector(sup.begin(), sup.end()),
                                 Rcpp::_["pumas_xi"] = Rcpp::NumericVector(xi.begin(), xi.end())));
}

}  // namespace

// Entry: files (eigen.bin of each block), blocks (ids), keep (0-based eigen rows with data per block, increasing), z / n / af per kept SNP (blocks
// concatenated in that order; af = frequency of the effect allele, for BETA = BETA_STD / sqrt(2 af (1 - af))),
// A (kept SNPs x annotations, same order; the intercept is added here), bin_file (LD-bin list, "" = no suppression),
// bin_start (first panel row of each kept block), n_panel (SNPs of the panel), sigma2_init (empty, or sigma2_k of the
// first grid point, residual first; raised to at least init_floor times the mean h2_j), thresh (eigen truncation, <= 0 = every stored
// pair), fp16 (store V as float16), control (viprs_annot_control()).

// [[Rcpp::export]]
Rcpp::List viprs_annot_cpp(const std::vector<std::string>& files, Rcpp::IntegerVector blocks, Rcpp::List keep,
                           Rcpp::NumericVector z, Rcpp::NumericVector n, Rcpp::NumericVector af, Rcpp::NumericMatrix A,
                           const std::string& bin_file, Rcpp::NumericVector bin_start, double n_panel, bool resid_all,
                           Rcpp::NumericVector sigma2_init, double thresh, bool fp16, Rcpp::List control,
                           int threads) {
  const int T = plink4r::resolve_threads(threads);
  const ViprsControl c = read_control(control);
  const Inputs in = read_inputs(blocks, keep, z, n, af);
  const int64_t M = in.z.size();
  if (A.nrow() != M) Rcpp::stop("A must have one row per kept SNP");
  const int K = A.ncol() + 1;
  const std::vector<double> s2(sigma2_init.begin(), sigma2_init.end());
  if (!s2.empty() && (int)s2.size() != K) Rcpp::stop("sigma2_init needs one value per annotation plus the residual");
  for (double x : s2) if (!std::isfinite(x)) Rcpp::stop("sigma2_init must be finite");
  std::vector<double> Ar(M * (K - 1));
  for (int64_t j = 0; j < M; ++j)
    for (int k = 0; k < K - 1; ++k) {
      const double x = A(j, k);
      if (!std::isfinite(x)) Rcpp::stop("A must be finite");
      Ar[j * (K - 1) + k] = x;
    }
  for (double x : Ar) if (x < 0) Rcpp::stop("annotations must be >= 0 (variance model sum_k A_jk sigma2_k)");
  // LD-bin suppression: annotated = any nonzero binary (0/1) annotation column
  std::vector<char> sup(M, 0);
  double n_bins = 0;
  std::vector<int> bincol;
  for (int a = 0; a < A.ncol(); ++a) {
    bool bin01 = true;
    for (int64_t j = 0; j < M && bin01; ++j) bin01 = A(j, a) == 0 || A(j, a) == 1;
    if (bin01) bincol.push_back(a);
  }
  try {
    if (!bin_file.empty() && !bincol.empty()) {
      if ((std::size_t)bin_start.size() != in.keep.size()) Rcpp::stop("bin_start needs one entry per block");
      const std::vector<int> bid = plink4r::read_bins16(
          bin_file, (int64_t)n_panel, std::vector<int64_t>(bin_start.begin(), bin_start.end()), in.keep);
      std::vector<char> ann(M, 0);
      for (int64_t j = 0; j < M; ++j)
        for (int a : bincol) if (A(j, a) != 0) { ann[j] = 1; break; }
      int64_t off = 0;
      for (const auto& kb : in.keep) {
        const int m = (int)kb.size();
        plink4r::bin_suppress(m, bid.data() + off, ann.data() + off, sup.data() + off);
        int mx = 0;
        for (int j = 0; j < m; ++j) mx = std::max(mx, bid[off + j] + 1);
        n_bins += mx;
        off += m;
      }
    }
    return fp16 ? fit_path<uint16_t>(in, Ar, K, sup, n_bins, resid_all, s2, files, thresh, c, T)
                : fit_path<float>(in, Ar, K, sup, n_bins, resid_all, s2, files, thresh, c, T);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}
