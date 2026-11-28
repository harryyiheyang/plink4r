// VIPRS with an additive annotation model of the per-SNP effect variance, and the partition of the effects into
// annotation-specific components (PartitionPRS, Research Design D.1.2). Top layer: the package's annotation VIPRS
// functional unit, built on vi_panel.h / vi_store.h / ld_bins.h; copied from viprs.cpp and edited.
// Model: beta_j = sum_k beta_j^(k), beta_j^(k) independent with variance A_jk sigma2_k (A_jk >= 0; k = 0 is the
// intercept, A_j0 = 1 for every SNP), so h2_j = E[beta_j^2] = sigma2_0 + sum_k A_jk sigma2_k. Prior on beta_j: VIPRS
// spike and slab, one pi (chosen on the grid); the slab N(0, h2_j / pi) (default; mix_prop / mix_var allow a fixed
// mixture of C classes with shares p_c and variances kap_c h2_j / pi, sum_c p_c kap_c = 1). Coordinate update
// (D_j = n_j r_j / sigma; tau_j = pi / h2_j):
//   vt_jc = n_j d_j / sigma + tau_j / kap_c,  u_c = log p_c + .5 D_j^2 / vt_jc + .5 (log tau_j - log kap_c - log vt_jc)
//   gamma_j = sigmoid(logit(pi) + logsumexp_c u_c), psi_jc = softmax(u)_c, eta_j = gamma_j sum_c psi_jc D_j / vt_jc.
// M-step (Research Design D.1.2): sigma2_k <- S_k / sum_j A_jk, S_k = sum_j [w_jk^2 E[beta_j^2] + w_jk (1 - w_jk) h2_j],
// whose fixed point sum_j A_jk (E[beta_j^2] - h2_j) / h2_j^2 = 0 is the stationary point of
// Fm = -1/2 sum_j (log h2_j + E[beta_j^2] / h2_j) (SNPs not suppressed); MASH prior on the annotation variances,
// sigma2_k = a_k sigma2_0, a_k in mash_grid (0, 1, 10, 100) with fixed weights mash_w common to every annotation, on
// the enrichment: E_k = sigma2_k / sigma2_0 of the free fit with its standard error se_k (sandwich over the LD blocks),
// one observation per annotation, q(a_k = m) proportional to w_m N(E_k; grid[m], se_k^2); VI uses the mean E_q[a_k],
// sigma2_0 by MM steps on Fm (VarModel::mstep); the Gibbs fit draws a_k from q. Start: sigma2_init (e.g.
// ldsc_strat() tau) or all sigma2_k equal with sum_j h2_j = the naive LDSC h2.
// Search by VI (below); final fit at the chosen pi by a short Gibbs sampler with Rao-Blackwell means
// (viprs_annot_gibbs.h; gibbs_draws = 0: the VI refit). Output: posterior mean b_j and its MMSE allocation
// b_j^(k) = w_jk b_j (Gibbs: the mean over the draws of w_jk m_j; sum_k b_j^(k) = b_j, k = 0 the intercept part), the
// weights of the annotation-specific PRS, and the Shapley value of each column in beta' R beta.
// LD-bin suppression (as sbayesrc()): inside an LD bin (ld_bins.h) with at least one annotated member (nonzero in
// any of the columns skip_cols), the members without annotation are fixed at gamma_j = 0 (beta_j = 0, never updated,
// no residual update) and left out of the M-step.
// Eigen cutoff, then pi, by VI. Stage 1, the cutoff (as sbayesrc()) on one PUMAS split (xi drawn once on the full
// panel, the truncated panel of each cutoff uses its leading entries), r = eta' valid / sqrt(eta' R eta): a chain over
// the cutoffs from the lowest up, each level warm-started from the one below, at the fixed pi cut_pi (default 0.01);
// the largest cutoff unless a smaller one has r / r_largest > rel_r_min. Stage 2, pi at the chosen cutoff,
// warm-started from its stage-1 fit: the pi grid ascending on all data with warm start (sigma2_k carried over when pi
// moves); the converged point of the largest ELBO + E_q[log w(a)] - E_q[log q(a)] is refitted to full convergence,
// then the Gibbs fit.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
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
  int full_every, refresh_every, min_iter, patience, max_iter, grid_max_iter, pi_steps, grid_extend,
      gibbs_burn, gibbs_draws;
  double rel_r_min, cut_pi;  // cut_pi: the fixed pi of the cutoff choice
  bool fix_sigma2;  // sigma2_k kept at the start (diagnostic: no M-step)
  std::vector<double> mash_grid, mash_w;  // MASH: sigma2_k in mash_grid * b with prior weights mash_w
  double start_scale, start_zero_share;  // given start (sigma2_init): see VarModel::start_at
  std::vector<double> cutoffs;  // eigen cutoffs compared at pi = cut_pi (one value: fixed)
  double act_tol, x_tol, f_tol, prop_train, sigma_fixed, h2_min, h2_max, pi_lo_count, pi_lo_min, pi_hi_count,
      pi_hi_max, pi_cap;
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
  x.seed = static_cast<unsigned long long>(Rcpp::as<double>(c["seed"]));
  x.gibbs_burn = Rcpp::as<int>(c["gibbs_burn"]);
  x.gibbs_draws = Rcpp::as<int>(c["gibbs_draws"]);
  x.cutoffs = Rcpp::as<std::vector<double>>(c["cutoffs"]);
  x.rel_r_min = Rcpp::as<double>(c["rel_r_min"]);
  x.cut_pi = Rcpp::as<double>(c["cut_pi"]);
  x.mash_grid = Rcpp::as<std::vector<double>>(c["mash_grid"]); x.mash_w = Rcpp::as<std::vector<double>>(c["mash_w"]);
  if (x.mash_grid.size() != x.mash_w.size() || x.mash_grid.empty()) Rcpp::stop("mash_grid and mash_w must have the same length");
  if (!(x.cut_pi > 0 && x.cut_pi < 1)) Rcpp::stop("cut_pi must be in (0, 1)");
  x.fix_sigma2 = Rcpp::as<bool>(c["fix_sigma2"]);
  x.start_scale = Rcpp::as<double>(c["start_scale"]);
  x.start_zero_share = Rcpp::as<double>(c["start_zero_share"]);
  if (x.cutoffs.empty()) Rcpp::stop("cutoffs: at least one eigen cutoff");
  for (double v : x.cutoffs) if (!(v > 0 && v <= 1)) Rcpp::stop("cutoffs must be in (0, 1]");
  std::sort(x.cutoffs.begin(), x.cutoffs.end());

  if (x.gibbs_burn < 0 || x.gibbs_draws < 0) Rcpp::stop("gibbs_burn and gibbs_draws must be >= 0");
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
// h2_j, for the ELBO) and e2v_j = E[beta_j^2 | data] (the M-step statistic) also returned per SNP. One class:
// zeta = e2v, [8] = sum gamma log vt (VIPRS).
void stats(int64_t M, const float* d, const float* n, const float* bt, const float* bv, const State& st, int64_t S,
           const float* tau, const float* ltau, const char* sup, const ViprsControl& c, float sigma, double* out,
           double* zeta, double* e2v, int threads) {
  double a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0, a8 = 0, a9 = 0;
  const double lo = 1e-12;
  const int C = (int)c.kap.size();
  const float* gamma = st.gamma.data(); const float* mu = st.mu.data(); const float* eta = st.eta.data();
  const float* s = st.s.data();
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : a0, a1, a2, a3, a4, a6, a7, a8, a9) num_threads(threads)
#endif
  for (int64_t j = 0; j < M; ++j) {
    if (sup[j]) { zeta[j] = 0; e2v[j] = 0; continue; }  // beta_j = 0 throughout: no term
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
    e2v[j] = gamma[j] * e2;  // E[beta_j^2 | data]
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
  int K;                         // annotations, intercept first
  const std::vector<double>& A;
  const std::vector<char>& sup;  // suppressed SNPs: left out of the M-step
  std::vector<double> par;       // sigma2_k
  std::vector<double> h2;        // h2_j
  std::vector<float> tau, ltau;  // tau_j = pi / h2_j and log tau_j, as used by the sweep
  double pi = 1;
  int threads;
  // MASH prior on the annotation variances: sigma2_k = a_k sigma2_0, a_k in grid (0, 1, 10, 100) with fixed prior
  // weights w (log: lw), common to every annotation. VI keeps q(a_k) (q, K x atoms, row 0 unused) and uses its mean
  // cbar_k = E_q[a_k] (par[k] = cbar_k sigma2_0); Gibbs draws a_k (cbar_k = the drawn atom). fixed: sigma2_k given
  // (no MASH, e.g. the oracle).
  std::vector<double> grid, lw, q, cbar;
  bool fixed = false;
  // M-step objective (Research Design D.1.2): Fm = -1/2 sum_j (log h2_j + e_j / h2_j) over the SNPs not suppressed,
  // e_j = E[beta_j^2 | data] the posterior second moment (spike included). Its stationarity in sigma2_k,
  // sum_j A_jk (e_j - h2_j) / h2_j^2 = 0, is the fixed point of the update sigma2_k <- S_k / sum_j A_jk,
  // S_k = sum_j [w_jk^2 e_j + w_jk (1 - w_jk) h2_j] (0/1 annotations).
  // MASH on the enrichment: the free fit fpar (every sigma2_k by MM on Fm) gives E_k = fpar_k / fpar_0 and its standard
  // error se_k (sandwich over the LD blocks blk: I^-1 (sum_b g_b g_b') I^-1, g_b the gradient of the SNPs of block b,
  // I = 1/2 sum_j a_j a_j' / h2_j^2; delta method for the ratio), one observation per annotation; then
  // q(a_k = m) proportional to w_m N(E_k; grid[m], se_k^2).
  std::vector<int> blk;  // LD block of each SNP
  int nblk = 0;
  std::vector<double> fpar, enr, enr_se;  // free fit, E_k, se_k (reported)
  void free_fit(const double* e, int steps) {
    if ((int)fpar.size() != K) fpar = par;
    for (int k = 0; k < K; ++k) if (!(fpar[k] > 0)) fpar[k] = 1e-3 * fpar[0];  // 0 is a fixed point of the MM step
    std::vector<double> hf(M);
    for (int it = 0; it < steps; ++it) {
      std::vector<double> NU(2 * K, 0.);
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
          const double* r = row(j);
          double h = fpar[0];
          for (int k = 1; k < K; ++k) h += r[k - 1] * fpar[k];
          const double u = e[j] / (h * h), v = 1. / h;
          loc[0] += u; loc[K] += v;
          for (int k = 1; k < K; ++k) { const double a = r[k - 1]; if (a != 0.) { loc[k] += a * u; loc[K + k] += a * v; } }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        for (int i = 0; i < 2 * K; ++i) NU[i] += loc[i];
      }
      for (int k = 0; k < K; ++k) if (NU[K + k] > 0 && NU[k] > 0) fpar[k] *= std::sqrt(NU[k] / NU[K + k]);
      fpar[0] = std::max(fpar[0], 1e-20);
    }
  }
  // E_k and se_k at the free fit
  void enrichment(const double* e) {
    Eigen::MatrixXd I = Eigen::MatrixXd::Zero(K, K), G = Eigen::MatrixXd::Zero(nblk, K);
    std::vector<double> a(K);
    for (int64_t j = 0; j < M; ++j) {
      if (sup[j]) continue;
      const double* r = row(j);
      a[0] = 1.;
      double h = fpar[0];
      for (int k = 1; k < K; ++k) { a[k] = r[k - 1]; h += a[k] * fpar[k]; }
      const double gi = .5 * (e[j] - h) / (h * h), ii = .5 / (h * h);
      for (int k = 0; k < K; ++k) {
        if (a[k] == 0.) continue;
        G(blk[j], k) += a[k] * gi;
        for (int l = k; l < K; ++l) if (a[l] != 0.) I(k, l) += a[k] * a[l] * ii;
      }
    }
    for (int k = 0; k < K; ++k) for (int l = 0; l < k; ++l) I(k, l) = I(l, k);
    // centre the block gradients (their sum is ~0 at the free fit; exactly 0 only at convergence)
    const Eigen::RowVectorXd gm = G.colwise().mean();
    G.rowwise() -= gm;
    const Eigen::MatrixXd Ii = I.completeOrthogonalDecomposition().pseudoInverse();
    const Eigen::MatrixXd V = Ii * (G.transpose() * G) * Ii;
    enr.assign(K, 0.); enr_se.assign(K, INFINITY);
    const double s0 = fpar[0];
    for (int k = 1; k < K; ++k) {
      enr[k] = fpar[k] / s0;
      const double d0 = -fpar[k] / (s0 * s0), dk = 1. / s0;
      const double v = dk * dk * V(k, k) + 2 * dk * d0 * V(0, k) + d0 * d0 * V(0, 0);
      enr_se[k] = v > 0 ? std::sqrt(v) : INFINITY;
    }
  }
  // q(a_k) from the enrichment; VI: cbar_k = E_q[a_k]; Gibbs (eng): a_k drawn from q. Returns the largest change of cbar.
  double update_q(std::mt19937_64* eng) {
    const size_t A_ = grid.size();
    if (q.size() != (size_t)K * A_) q.assign((size_t)K * A_, 0.);
    double dmax = 0;
    std::vector<double> u(A_);
    for (int k = 1; k < K; ++k) {
      double um = -INFINITY;
      const double se = enr_se[k];
      for (size_t m = 0; m < A_; ++m) {
        const double d = (enr[k] - grid[m]) / se;
        u[m] = lw[m] + (std::isfinite(se) ? -.5 * d * d : 0.);
        um = std::max(um, u[m]);
      }
      double z = 0;
      for (size_t m = 0; m < A_; ++m) { u[m] = std::exp(u[m] - um); z += u[m]; }
      double cb = 0;
      for (size_t m = 0; m < A_; ++m) { q[(size_t)k * A_ + m] = u[m] / z; cb += u[m] / z * grid[m]; }
      if (eng) {
        double x = std::uniform_real_distribution<double>(0., z)(*eng);
        size_t pick = 0;
        for (; pick < A_ - 1 && (x -= u[pick]) > 0; ++pick) {}
        cb = grid[pick];
      }
      dmax = std::max(dmax, std::fabs(cb - cbar[k]) / std::max(cbar[k], 1.));
      cbar[k] = cb; par[k] = cb * par[0];
    }
    set_all();
    return dmax;
  }
  // E_q[log w(a_k)] - E_q[log q(a_k)] summed over the annotations: the MASH part of the ELBO (0 when fixed)
  double mash_kl() const {
    double f = 0;
    if (fixed || (int)q.size() != K * (int)grid.size()) return 0;
    const size_t A_ = grid.size();
    for (int k = 1; k < K; ++k)
      for (size_t m = 0; m < A_; ++m) {
        const double p = q[(size_t)k * A_ + m];
        if (p > 0) f += p * (lw[m] - std::log(p));
      }
    return f;
  }
  // MM steps on sigma2_0 with cbar fixed (h2_j = sigma2_0 abar_j, abar_j = 1 + sum_k A_jk cbar_k):
  // sigma2_0 <- sigma2_0 sqrt(sum_j abar_j e_j / h2_j^2 / sum_j abar_j / h2_j) (Fm never decreases)
  void mm_s0(const double* e, int steps) {
    for (int it = 0; it < steps; ++it) {
      double Nb = 0, Db = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : Nb, Db) num_threads(threads)
#endif
      for (int64_t j = 0; j < M; ++j) {
        if (sup[j]) continue;
        const double h = h2[j], a = h / par[0];
        Nb += a * e[j] / (h * h); Db += a / h;
      }
      if (Db > 0 && Nb > 0) par[0] *= std::sqrt(Nb / Db);
      par[0] = std::max(par[0], 1e-20);
      for (int k = 1; k < K; ++k) par[k] = cbar[k] * par[0];
      set_all();
    }
  }
  VarModel(int64_t M_, int K_, const std::vector<double>& A_, const std::vector<char>& sup_, int threads_)
      : M(M_), K(K_), A(A_), sup(sup_), par(K_, 0.), h2(M_, 1.), tau(M_, 1.f), ltau(M_, 0.f), threads(threads_) {}

  const double* row(int64_t j) const { return A.data() + j * (K - 1); }
  // MMSE weight of piece k of SNP j: w_jk = A_jk sigma2_k / h2_j (A_j0 = 1)
  double w(int64_t j, int k) const { return (k ? row(j)[k - 1] : 1.) * par[k] / h2[j]; }
  void set_all() {
    if (!(par[0] > 0)) throw std::runtime_error("the intercept sigma2_0 is not positive");
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (int64_t j = 0; j < M; ++j) {
      double h = par[0];
      const double* r = row(j);
      for (int k = 1; k < K; ++k) h += r[k - 1] * par[k];
      h2[j] = h;
      const double t = pi / h;
      tau[j] = (float)t; ltau[j] = (float)std::log(t);
    }
  }
  // start: all sigma2_k equal, sum_j h2_j = h2 over the SNPs not suppressed (LDSC start)
  void start(double h2_0, double pi_) {
    pi = pi_;
    double s = 0;
    for (int64_t j = 0; j < M; ++j) {
      if (sup[j]) continue;
      s += 1.;
      const double* r = row(j);
      for (int k = 1; k < K; ++k) s += r[k - 1];
    }
    std::fill(par.begin(), par.end(), h2_0 / s);
    set_all();
  }
  // M-step of the variance model with the MASH prior on the enrichment: `inner` MM steps of the free fit, E_k and se_k,
  // q(a_k) (VI: the mean; Gibbs, eng: a draw), then `inner` MM steps on sigma2_0 with sigma2_k = cbar_k sigma2_0
  // (e: VI second moments; Gibbs: the stochastic-approximation average of the RB second moments).
  void mstep(const double* e, int inner = 10, std::mt19937_64* eng = nullptr) {
    if (fixed) return;
    if ((int)cbar.size() != K) { cbar.assign(K, 0.); for (int k = 1; k < K; ++k) cbar[k] = par[k] / par[0]; }
    free_fit(e, inner);
    if (K > 1) { enrichment(e); update_q(eng); }
    mm_s0(e, inner);
  }
  // start from given sigma2_k (e.g. stratified LDSC tau_k), h0 the LDSC h2: negative values set to 0; every sigma2_k
  // times scale (start_scale); a column at 0 (intercept included) set so that its h2 (sum_j A_jk sigma2_k) is zero_share
  // (start_zero_share) of h0.
  // Returns sum_j h2_j over the SNPs not suppressed.
  double start_at(const std::vector<double>& s2, double pi_, double h0, double scale, double zero_share) {
    std::vector<double> sa(K, 0.);  // sum_j A_jk over the SNPs not suppressed
    for (int64_t j = 0; j < M; ++j) {
      if (sup[j]) continue;
      sa[0] += 1.;
      const double* r = row(j);
      for (int k = 1; k < K; ++k) sa[k] += r[k - 1];
    }
    std::vector<double> p(K);
    for (int k = 0; k < K; ++k)
      p[k] = s2[k] > 0 ? scale * s2[k] : (sa[k] > 0 ? zero_share * h0 / sa[k] : 0.);
    double tot = 0;
    for (int k = 0; k < K; ++k) tot += sa[k] * p[k];
    par = p; pi = pi_;
    set_all();
    return tot;
  }
  // pi moves: sigma2_k (variance of beta) kept, the slab variance h2_j / pi follows
  void pi_move(double pi_) { pi = pi_; set_all(); }

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
    std::vector<double> zeta(M), e2v(M);
    int64_t n_free = 0;  // SNPs not suppressed (the suppressed ones are outside the model, as in stats())
    for (int64_t j = 0; j < M; ++j) n_free += !vm.sup[j];
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
      stats(M, P.d.data(), nv, bt, bv, st, P.S, vm.tau.data(), vm.ltau.data(), vm.sup.data(), c, (float)sigma, out, zeta.data(), e2v.data(), threads);
      const double sg_sum = out[0], sdz = out[2], eb = out[3], sde2 = out[4], ss = out[5], glg = out[6],
                   g1lg1 = out[7], glvt = out[8];
      if (!c.fix_sigma2) {
        vm.mstep(e2v.data());
      }
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
      elbo += -glg + sg_sum * std::log(pi) - g1lg1 + (n_free - sg_sum) * std::log1p(-pi);
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

  // PUMAS split (magenpy / VIPRS): training half bt + sqrt(c) V xi, xi ~ N(0, I_S), c = 1/(prop_train N) - 1/N.
  // xi of block b from its own stream (seed, b), so the panel at fewer eigen components uses the leading entries of the
  // same draws (one split across eigen cutoffs).
  void pumas(std::vector<float>& train, std::vector<float>& valid, std::vector<float>& n_tr, std::vector<double>& xi) {
    std::normal_distribution<double> N01;
    xi.resize(P.S);
    for (int b = 0; b < P.nb; ++b) {
      std::mt19937_64 eng(c.seed + 0x9E3779B97F4A7C15ull * (unsigned long long)(b + 1));
      for (int i = 0; i < P.k[b]; ++i) xi[P.s_off[b] + i] = N01(eng);
    }
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
};

// ---- output

// annotation heritability h2_k = (sum_j A_jk) sigma2_k over the SNPs not suppressed (k = 0: intercept)
Rcpp::NumericVector h2_annot(const VarModel& vm) {
  Rcpp::NumericVector h(vm.K);
  for (int64_t j = 0; j < vm.M; ++j) {
    if (vm.sup[j]) continue;
    const double* r = vm.row(j);
    h[0] += vm.par[0];
    for (int k = 1; k < vm.K; ++k) h[k] += r[k - 1] * vm.par[k];
  }
  return h;
}

#include "viprs_annot_gibbs.h"

// Output from the final fit: Gibbs + RB (G.draws > 0) or the VI state (st, vm). Columns: the mean of w_jk m_j (Gibbs)
// or w_jk eta_j (VI); per SNP they sum to BETA_STD. Heritability per column: the Shapley value in beta'R beta,
// h2_shapley = mean over the draws of beta^(k)' R beta (VI: b_k' R b).
template <typename T>
Rcpp::List pack(const Fitter<T>& F, const State& st, const RunResult& res, const VarModel& vm, const GibbsOut& G,
                const std::vector<double>& af, double h2_0, Rcpp::List extra) {
  const EigenPanel<T>& P = F.P;
  const int64_t M = P.M;
  const int K = vm.K;
  const bool gib = G.draws > 0;
  Rcpp::NumericVector beta(M), bstd(M), pip(M), h2j(M);
  Rcpp::NumericMatrix part(M, K);  // standardized scale, k = 0 intercept
  std::vector<float> bf(M);
  for (int64_t j = 0; j < M; ++j) {
    const double p = std::min(std::max(af[j], 1e-6), 1 - 1e-6);
    const double b = gib ? G.bmean[j] : st.eta[j];
    bstd[j] = b; bf[j] = (float)b;
    beta[j] = b / std::sqrt(2. * p * (1. - p));
    pip[j] = gib ? G.pip[j] : st.gamma[j];
    h2j[j] = vm.h2[j];
    for (int k = 0; k < K; ++k) part(j, k) = gib ? G.part[(size_t)j * K + k] : vm.w(j, k) * b;
  }
  std::vector<float> sb(P.S);
  vtx(P, bf.data(), sb.data(), F.threads);
  std::vector<double> q(M);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(F.threads)
#endif
  for (int b = 0; b < P.nb; ++b)
    for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j)
      q[j] = RowStore<T>::dot(P.row(b, j), sb.data() + P.s_off[b], P.k[b]);
  Rcpp::NumericVector hsh(K);
  if (gib) for (int k = 0; k < K; ++k) hsh[k] = G.phi[k];
  else for (int64_t j = 0; j < M; ++j) for (int k = 0; k < K; ++k) hsh[k] += part(j, k) * q[j];
  Rcpp::List hyp = Rcpp::List::create(
      Rcpp::_["pi"] = res.pi, Rcpp::_["sigma"] = gib ? G.sigma : res.sigma, Rcpp::_["h2"] = gib ? G.h2 : res.h2,
      Rcpp::_["elbo"] = res.elbo, Rcpp::_["iters"] = res.iters,
      Rcpp::_["converged"] = res.converged, Rcpp::_["message"] = res.message, Rcpp::_["sigma_fixed"] = res.sigma_fixed,
      Rcpp::_["gibbs_draws"] = G.draws,
      Rcpp::_["h2_ldsc_init"] = h2_0, Rcpp::_["n_snps"] = (double)M, Rcpp::_["sum_k"] = (double)F.P.S,
      Rcpp::_["sweeps_full"] = (double)F.sweeps_full, Rcpp::_["sweeps_partial"] = (double)F.sweeps_part,
      Rcpp::_["snp_visits_per_M"] = (double)F.snp_visits / M);
  return Rcpp::List::create(Rcpp::_["BETA"] = beta, Rcpp::_["BETA_STD"] = bstd, Rcpp::_["PIP"] = pip,
                            Rcpp::_["H2_SNP"] = h2j, Rcpp::_["PART_STD"] = part,
                            Rcpp::_["sigma2"] = Rcpp::NumericVector(vm.par.begin(), vm.par.end()),
                            Rcpp::_["h2_annot"] = h2_annot(vm), Rcpp::_["h2_shapley"] = hsh,
                            Rcpp::_["hyp"] = hyp,
                            Rcpp::_["extra"] = extra);
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

// variance model start: sigma2_init (intercept > 0) or all sigma2_k equal with sum_j h2_j = the LDSC h2; returns the
// start h2 (sigma = 1 - it)
double start_model(VarModel& vm, const std::vector<double>& s2_init, double h2_0, double pi, const ViprsControl& c) {
  if (!s2_init.empty()) {
    const double tot = vm.start_at(s2_init, pi, h2_0, c.start_scale, c.start_zero_share);
    return std::min(std::max(tot, c.h2_min), c.h2_max);
  }
  vm.start(h2_0, pi);
  return h2_0;
}

// Cutoff, then pi, by VI on one PUMAS split: the panel is read once with every stored pair (P0); xi is drawn once per
// block (Fitter::pumas), the panel at a smaller cutoff uses its leading entries, so every cutoff trains on the same
// pseudo data, and every fit is scored on the validation half of P0 with R of P0, r = b' valid / sqrt(b' R b).
// Stage 1: a chain over the cutoffs in ascending order at pi = cut_pi, each level warm-started from the one below
// (max_iter iterations, scored converged or not); the largest cutoff unless a smaller one has |r / r_largest| >
// rel_r_min (SBayesRC rule). Stage 2: at the chosen cutoff, from its stage-1 fit, the pi grid ascending on all data
// with warm start (state and variance model carried over); the converged point of the largest ELBO (plus the MASH part),
// refitted to full convergence, then the final Gibbs fit.
template <typename T>
Rcpp::List fit_path(const Inputs& in, const std::vector<double>& A, int K, const std::vector<char>& sup, double n_bins,
                    const std::vector<double>& s2_init, const std::vector<std::string>& files,
                    const ViprsControl& c, int threads) {
  std::unique_ptr<EigenPanel<T>> P0(
      new EigenPanel<T>(plink4r::load_eigen_panel<T>(files, in.blocks, in.keep, in.bstd, 0., threads)));
  const int nc = (int)c.cutoffs.size();
  std::vector<int> blk0(P0->M);  // LD block of each SNP (the jackknife units of the enrichment se)
  for (int b = 0; b < P0->nb; ++b) for (int64_t j = P0->snp_off[b]; j < P0->snp_off[b + 1]; ++j) blk0[j] = b;
  const int nblk0 = P0->nb;
  std::vector<std::vector<int>> ks(nc, std::vector<int>(P0->nb));
  for (int i = 0; i < nc; ++i)
    for (int b = 0; b < P0->nb; ++b) {
      int m = 0;
      plink4r::eigen_block_dims(files[b], c.cutoffs[i], m, ks[i][b], true);
      ks[i][b] = std::min(ks[i][b], P0->k[b]);
    }
  std::vector<float> valid0;
  std::vector<double> xi0;
  {
    Fitter<T> F0(*P0, c, in.n, in.z, threads);
    std::vector<float> tr0, ntr0;
    F0.pumas(tr0, valid0, ntr0, xi0);
  }
  std::vector<float> s0(P0->S);
  auto score = [&](const State& st, double& rb, double& bsb) {
    rb = 0; bsb = 0;
    for (int64_t j = 0; j < P0->M; ++j) rb += (double)st.eta[j] * valid0[j];
    vtx(*P0, st.eta.data(), s0.data(), threads);
    for (float x : s0) bsb += (double)x * x;
  };
  std::vector<double> gcut, gpi, gsig, gh2, gr, grb, gelbo, gsec;
  std::vector<int> gstage, git, gconv;
  std::vector<std::vector<double>> gpar;
  struct Best { bool have = false; double r = -INFINITY, pi = 0, sigma = 0; State st{0, 0, 0.}; std::vector<double> par, cbar, q, fpar; };
  double h2_0 = 0;
  // MASH settings of a variance model (fix_sigma2: sigma2_k kept at the start, no MASH)
  auto setsel = [&](VarModel& v) {
    v.grid = c.mash_grid; v.lw.clear(); v.fixed = c.fix_sigma2; v.blk = blk0; v.nblk = nblk0;
    for (double w : c.mash_w) v.lw.push_back(w > 0 ? std::log(w) : -INFINITY);
  };
  auto keep = [&](Best& o, const VarModel& v, double r, double pi, double sigma, const State& st) {
    o.have = true; o.r = r; o.pi = pi; o.sigma = sigma; o.st = st; o.par = v.par; o.cbar = v.cbar; o.q = v.q; o.fpar = v.fpar;
  };
  // pi path on panel P at cutoff ic (ascending, warm start), each point recorded in the grid table. Stage 1: fitted on
  // the PUMAS training half with max_iter iterations, scored converged or not (as the fixed chain of SBayesRC). Stage 2: fitted on all data, the converged
  // point of the largest ELBO
  // init: warm start (effects, sigma2_k, sigma of an earlier fit; s = V' eta is recomputed on this panel)
  auto path = [&](const EigenPanel<T>& P, int ic, int stage, const std::vector<double>& grid, const Best* init) {
    Fitter<T> F(P, c, in.n, in.z, threads);
    const double h0 = F.h2_ldsc();
    std::vector<float> train, valid, n_tr;
    std::vector<double> xi;
    const bool pseudo = stage == 1;
    if (pseudo) F.pumas(train, valid, n_tr, xi);
    const float* nv = pseudo ? n_tr.data() : in.n.data();
    const float* bv = pseudo ? train.data() : P.bt.data();
    State st(P.M, P.S, grid[0]);
    VarModel vm(P.M, K, A, sup, threads);
    setsel(vm);
    double sigma;
    if (init) {
      st = init->st; st.s.assign(P.S, 0.f);
      vm.par = init->par; vm.cbar = init->cbar; vm.q = init->q; vm.fpar = init->fpar; vm.pi = grid[0]; vm.set_all();
      sigma = init->sigma;
    } else {
      for (int64_t j = 0; j < P.M; ++j) if (sup[j]) st.gamma[j] = 0.f;
      sigma = 1. - start_model(vm, s2_init, h0, grid[0], c);
    }
    Best bb, bn;  // bn: the best point when none converged (stage 2 falls back on it with a warning)
    for (double pi : grid) {
      vm.pi_move(pi);
      const auto t1 = std::chrono::steady_clock::now();
      const RunResult r = F.run(st, nv, bv, pi, vm, sigma, stage == 1 ? c.max_iter : c.grid_max_iter, nullptr);
      sigma = r.sigma;
      double rb = NAN, rr = NAN;
      if (pseudo) {
        double bsb;
        score(st, rb, bsb);
        rr = bsb > 0 ? rb / std::sqrt(bsb) : 0.;
      }
      const double elbo = r.elbo + vm.mash_kl();  // ELBO with the MASH part E_q[log w(a)] - E_q[log q(a)]
      const double crit = pseudo ? rr : elbo;  // pseudo r: only r > 0 can be chosen
      gstage.push_back(stage); gcut.push_back(c.cutoffs[ic]); gpi.push_back(pi); gsig.push_back(sigma);
      gh2.push_back(r.h2); gr.push_back(rr); grb.push_back(rb); gelbo.push_back(elbo); git.push_back(r.iters);
      gconv.push_back(r.converged); gpar.push_back(vm.par);
      gsec.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count());
      const double rr_ = crit;
      if ((r.converged || stage == 1) && std::isfinite(rr_) && (!pseudo || rr_ > 0) && rr_ > bb.r) {
        keep(bb, vm, rr_, pi, sigma, st);
      } else if (!bb.have && std::isfinite(rr_) && (!pseudo || rr_ > 0) && rr_ > bn.r) {
        keep(bn, vm, rr_, pi, sigma, st);
      }
      Rcpp::checkUserInterrupt();
    }
    if (!bb.have && bn.have) {
      Rcpp::warning("no pi of the grid converged in grid_max_iter iterations; the best unconverged point (pi = %g) is "
                    "refitted", bn.pi);
      return bn;
    }
    return bb;
  };
  // stage 1: the cutoff at pi = cut_pi, a chain from the lowest cutoff (most shrinkage) up, each level warm-started from
  // the one below (as sbayesrc())
  std::vector<Best> bc(nc);
  std::vector<int> ord(nc);
  for (int i = 0; i < nc; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(), [&](int x, int y) { return c.cutoffs[x] < c.cutoffs[y]; });
  int ibest = -1;
  const Best* prev = nullptr;
  for (int ic : ord) {
    std::unique_ptr<EigenPanel<T>> Pi;
    if (ks[ic] != P0->k) Pi.reset(new EigenPanel<T>(plink4r::truncate_eigen_panel(*P0, ks[ic], threads)));
    bc[ic] = path(Pi ? *Pi : *P0, ic, 1, std::vector<double>{c.cut_pi}, prev);
    if (bc[ic].have) prev = &bc[ic];
    if (bc[ic].have && (ibest < 0 || bc[ic].r > bc[ibest].r || (bc[ic].r == bc[ibest].r && c.cutoffs[ic] < c.cutoffs[ibest])))
      ibest = ic;
  }
  int ib = -1;  // base: the largest cutoff with a usable fit
  for (int i = nc - 1; i >= 0 && ib < 0; --i) if (bc[ord[i]].have) ib = ord[i];
  if (ib < 0) Rcpp::stop("no cutoff with a positive pseudo r at pi = cut_pi");
  int ich = ib;
  if (ibest >= 0 && ibest != ib && bc[ibest].r / bc[ib].r > c.rel_r_min) ich = ibest;
  // stage 2: pi at the chosen cutoff
  std::unique_ptr<EigenPanel<T>> Pc;
  if (ks[ich] != P0->k) { Pc.reset(new EigenPanel<T>(plink4r::truncate_eigen_panel(*P0, ks[ich], threads))); }
  const EigenPanel<T>& P = Pc ? *Pc : *P0;
  Best best;
  {
    Fitter<T> Fg(P, c, in.n, in.z, threads);
    best = path(P, ich, 2, Fg.pi_grid(), &bc[ich]);  // warm from the stage-1 fit of the chosen cutoff
  }
  if (!best.have)
    Rcpp::stop("no pi of the grid converged at the chosen cutoff");
  if (Pc) P0.reset();
  Fitter<T> F(P, c, in.n, in.z, threads);
  h2_0 = F.h2_ldsc();
  VarModel bvm(P.M, K, A, sup, threads);
  setsel(bvm);
  bvm.par = best.par; bvm.cbar = best.cbar; bvm.q = best.q; bvm.fpar = best.fpar; bvm.pi = best.pi; bvm.set_all();
  State bst = best.st;
  RunResult res = F.run(bst, in.n.data(), P.bt.data(), best.pi, bvm, best.sigma, c.max_iter, nullptr);
  const std::vector<double> qvi = bvm.q, cbvi = bvm.cbar, envi = bvm.enr, sevi = bvm.enr_se;
  GibbsOut G;
  if (c.gibbs_draws > 0) G = gibbs_fit(P, in.n, c, bst, bvm, res.sigma, c.gibbs_burn, c.gibbs_draws, threads);
  // q(a_k) of the VI refit (rows: annotations, intercept row 0 unused; columns: the atoms)
  const int na = (int)c.mash_grid.size();
  Rcpp::NumericMatrix qmat(K, na);
  if (qvi.size() == (size_t)K * na) for (int k = 0; k < K; ++k) for (int m = 0; m < na; ++m) qmat(k, m) = qvi[(size_t)k * na + m];
  const int ng = gpi.size();
  Rcpp::NumericMatrix gcoef(ng, K);
  for (int g = 0; g < ng; ++g) for (int k = 0; k < K; ++k) gcoef(g, k) = gpar[g][k];
  Rcpp::DataFrame tab = Rcpp::DataFrame::create(
      Rcpp::_["stage"] = gstage, Rcpp::_["cutoff"] = gcut, Rcpp::_["pi"] = gpi, Rcpp::_["sigma"] = gsig, Rcpp::_["h2"] = gh2,
      Rcpp::_["pseudo_r"] = gr, Rcpp::_["rb"] = grb, Rcpp::_["elbo"] = gelbo, Rcpp::_["iters"] = git, Rcpp::_["converged"] = gconv,
      Rcpp::_["seconds"] = gsec);
  return pack(F, bst, res, bvm, G, in.af, h2_0,
              Rcpp::List::create(Rcpp::_["grid"] = tab, Rcpp::_["grid_sigma2"] = gcoef,
                                 Rcpp::_["mash"] = Rcpp::List::create(
                                     Rcpp::_["sigma2_0"] = bvm.par[0],  // Gibbs: the mean over the draws
                                     Rcpp::_["q"] = qmat, Rcpp::_["mean_atom"] = Rcpp::NumericVector(cbvi.begin(), cbvi.end()),
                                     Rcpp::_["enrichment"] = Rcpp::NumericVector(envi.begin(), envi.end()),
                                     Rcpp::_["enrichment_se"] = Rcpp::NumericVector(sevi.begin(), sevi.end()),
                                     Rcpp::_["pzero_gibbs"] = Rcpp::NumericVector(G.pzero.begin(), G.pzero.end()),
                                     Rcpp::_["draws"] = G.draws),
                                 Rcpp::_["best_pseudo_r"] = bc[ich].r, Rcpp::_["best_pi_crit"] = best.r, Rcpp::_["n_bins"] = n_bins,
                                 Rcpp::_["cutoff"] = c.cutoffs[ich],
                                 Rcpp::_["suppressed"] = Rcpp::LogicalVector(sup.begin(), sup.end())));
}

}  // namespace

// Entry: files (eigen.bin of each block), blocks (ids), keep (0-based eigen rows with data per block, increasing), z / n / af per kept SNP (blocks
// concatenated in that order; af = frequency of the effect allele, for BETA = BETA_STD / sqrt(2 af (1 - af))),
// A (kept SNPs x annotations, same order; the intercept is added here), bin_file (LD-bin list, "" = no
// suppression), skip_cols (0-based columns of A that make a SNP annotated for the suppression),
// bin_start (first panel row of each kept block), n_panel (SNPs of the panel), sigma2_init (empty, or sigma2_k of the
// first grid point, intercept first; annotation values cut at 0), fp16 (store V as float16), control (viprs_annot_control()).

// [[Rcpp::export]]
Rcpp::List viprs_annot_cpp(const std::vector<std::string>& files, Rcpp::IntegerVector blocks, Rcpp::List keep,
                           Rcpp::NumericVector z, Rcpp::NumericVector n, Rcpp::NumericVector af, Rcpp::NumericMatrix A,
                           const std::string& bin_file, Rcpp::NumericVector bin_start, double n_panel,
                           Rcpp::NumericVector sigma2_init, Rcpp::IntegerVector skip_cols, bool fp16,
                           Rcpp::List control, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const ViprsControl c = read_control(control);
  const Inputs in = read_inputs(blocks, keep, z, n, af);
  const int64_t M = in.z.size();
  if (A.nrow() != M) Rcpp::stop("A must have one row per kept SNP");
  const int K = A.ncol() + 1;
  const std::vector<double> s2(sigma2_init.begin(), sigma2_init.end());
  if (!s2.empty() && (int)s2.size() != K) Rcpp::stop("sigma2_init needs one value per annotation plus the intercept");
  for (double x : s2) if (!std::isfinite(x)) Rcpp::stop("sigma2_init must be finite");
  std::vector<double> Ar(M * (K - 1));
  for (int64_t j = 0; j < M; ++j)
    for (int k = 0; k < K - 1; ++k) {
      const double x = A(j, k);
      if (!std::isfinite(x)) Rcpp::stop("A must be finite");
      Ar[j * (K - 1) + k] = x;
    }
  for (double x : Ar) if (x < 0) Rcpp::stop("annotations must be >= 0 (variance model sum_k A_jk sigma2_k)");
  // LD-bin suppression: annotated = any nonzero value in the columns skip_cols (0-based)
  std::vector<char> sup(M, 0);
  double n_bins = 0;
  const std::vector<int> bincol(skip_cols.begin(), skip_cols.end());
  for (int a : bincol) if (a < 0 || a >= A.ncol()) Rcpp::stop("skip_cols outside A");
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
    return fp16 ? fit_path<uint16_t>(in, Ar, K, sup, n_bins, s2, files, c, T)
                : fit_path<float>(in, Ar, K, sup, n_bins, s2, files, c, T);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}
