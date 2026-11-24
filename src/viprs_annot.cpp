// VIPRS with an annotation model of the per-SNP slab variance, log link (top layer: the package's annotation VIPRS
// functional unit, built on vi_panel.h / vi_store.h; copied from viprs.cpp and edited).
// Slab of SNP j: beta_j | gamma_j = 1 ~ N(0, v_j), v_j = exp(a_j' alpha), a_j = (1, annotations of j) (LDSC-type
// variance model on the log scale). Coordinate update as viprs.cpp with tau_j = 1 / v_j in place of tau:
//   vt_j = n_j d_j / sigma + tau_j,  gamma_j = sigmoid(.5 vt_j mu_j^2 + logit(pi) + .5 (log tau_j - log vt_j))
// M-step of alpha (pi fixed on the grid): maximize Q(alpha) = sum_j [-.5 pi a_j' alpha - .5 zeta_j exp(-a_j' alpha)],
// zeta_j = gamma_j (mu_j^2 + 1 / vt_j), weight pi per SNP as the VIPRS 0.1.5 tau update (tau = pi M / sum zeta), so
// the intercept-only model is viprs(). Q is concave in alpha: Newton steps with step halving, to mstep_tol.
// LD-bin suppression (as sbayesrc()): inside an LD bin (ld_bins.h) with at least one annotated member (any nonzero
// binary 0/1 annotation column), the members without annotation are fixed at gamma_j = 0 (beta_j = 0, never updated,
// no residual update) and left out of the M-step of alpha.
// Hyperparameter search as viprs.cpp (pi grid ascending with warm start, PUMAS, best converged point refitted);
// a move of pi shifts the intercept by log(prev_pi / pi) (keeps the implied h2).
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
  int mstep_max_iter, mstep_max_halving;
  double act_tol, x_tol, f_tol, prop_train, sigma_fixed, h2_min, h2_max, pi_lo_count, pi_lo_min, pi_hi_count,
      pi_hi_max, pi_cap, mstep_tol;
  unsigned long long seed;
  std::vector<double> pi_grid_user;  // empty: official grid
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
  x.mstep_max_iter = Rcpp::as<int>(c["mstep_max_iter"]);
  x.mstep_tol = Rcpp::as<double>(c["mstep_tol"]);
  x.mstep_max_halving = Rcpp::as<int>(c["mstep_max_halving"]);
  x.seed = static_cast<unsigned long long>(Rcpp::as<double>(c["seed"]));
  if (!Rf_isNull(c["pi_grid"])) x.pi_grid_user = Rcpp::as<std::vector<double>>(c["pi_grid"]);
  std::sort(x.pi_grid_user.begin(), x.pi_grid_user.end());
  return x;
}

struct State {
  std::vector<float> gamma, mu, eta, s, r;
  std::vector<uint8_t> active;
  State(int64_t M, int64_t S, double pi)
      : gamma(M, (float)pi), mu(M, 0.f), eta(M, 0.f), s(S, 0.f), r(M, 0.f), active(M, 1) {}
};

// ---- kernels on the panel

template <typename T>
void sweep(const EigenPanel<T>& P, const float* n, const float* bt, float pi, const float* tau, const char* sup,
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
      const float vt = n[j] * d[j] / sigma + tau[j];
      const float mm = n[j] / (vt * sigma);
      const bool lazy = active_only && !active[j];
      if (!lazy) {
        ++cnt;
        r[j] = bt[j] - (RowStore<T>::dot(vj, sb, k) - d[j] * eta[j]);
      }
      const float m_ = mm * r[j];
      const float u = std::sqrt(.5f * vt) * m_;
      const float x = u * u + lp + .5f * (std::log(tau[j]) - std::log(vt));
      const float g = x < 0 ? std::exp(x) / (1.f + std::exp(x)) : 1.f / (1.f + std::exp(-x));
      const float de = g * m_ - eta[j];
      const float ad = std::fabs(de);
      if (ad >= eps) {
        mu[j] = m_;
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

// Sums for the M-step and the ELBO with vt from (tau, sigma) of the last sweep:
// [0] sum gamma, [1] sum zeta, [2] sum d zeta, [3] eta'bt, [4] sum d eta^2, [5] s's, [6] sum gamma log gamma,
// [7] sum (1-gamma) log(1-gamma), [8] sum gamma log vt, [9] eta'bv (bv may be null); zeta = gamma (mu^2 + 1/vt),
// also returned per SNP
void stats(int64_t M, const float* d, const float* n, const float* bt, const float* bv, const State& st, int64_t S,
           const float* tau, float sigma, double* out, double* zeta, int threads) {
  double a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0, a7 = 0, a8 = 0, a9 = 0;
  const double lo = 1e-12;
  const float* gamma = st.gamma.data(); const float* mu = st.mu.data(); const float* eta = st.eta.data();
  const float* s = st.s.data();
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : a0, a1, a2, a3, a4, a6, a7, a8, a9) num_threads(threads)
#endif
  for (int64_t j = 0; j < M; ++j) {
    const double vt = (double)n[j] * d[j] / sigma + tau[j];
    const double g = std::min(std::max((double)gamma[j], lo), 1. - lo);
    const double z = gamma[j] * ((double)mu[j] * mu[j] + 1. / vt);
    a0 += gamma[j];
    a1 += z;
    zeta[j] = z;
    a2 += d[j] * z;
    a3 += (double)eta[j] * bt[j];
    a4 += (double)d[j] * eta[j] * eta[j];
    a6 += g * std::log(g);
    a7 += (1. - g) * std::log(1. - g);
    a8 += gamma[j] * std::log(vt);
    if (bv) a9 += (double)eta[j] * bv[j];
  }
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : a5) num_threads(threads)
#endif
  for (int64_t i = 0; i < S; ++i) a5 += (double)s[i] * s[i];
  out[0] = a0; out[1] = a1; out[2] = a2; out[3] = a3; out[4] = a4; out[5] = a5;
  out[6] = a6; out[7] = a7; out[8] = a8; out[9] = a9;
}

// ---- variance model: v_j = exp(a_j' alpha), a_j = (1, A_j) (A row-major M x (K - 1))

struct VarModel {
  int64_t M;
  int K;                      // coefficients, intercept first
  const std::vector<double>& A;
  const std::vector<char>& sup;  // suppressed SNPs: left out of Q
  std::vector<double> par;    // alpha
  std::vector<float> tau;     // 1 / v_j
  int threads;

  VarModel(int64_t M_, int K_, const std::vector<double>& A_, const std::vector<char>& sup_, int threads_)
      : M(M_), K(K_), A(A_), sup(sup_), par(K_, 0.), tau(M_, 1.f), threads(threads_) {}

  double lin(int64_t j, const double* a) const {
    double x = a[0];
    const double* r = A.data() + j * (K - 1);
    for (int k = 1; k < K; ++k) x += r[k - 1] * a[k];
    return x;
  }
  void set_tau() {
    for (int64_t j = 0; j < M; ++j) tau[j] = (float)std::exp(-lin(j, par.data()));
  }
  // start: v_j = h2 / (pi M) for all SNPs (M over all SNPs, as the sigma2 start of sbayesrc())
  void start(double h2, double pi) {
    std::fill(par.begin(), par.end(), 0.);
    par[0] = std::log(h2 / (pi * M));
    set_tau();
  }
  // pi moves from prev to pi: v_j * prev / pi (keeps the implied h2)
  void pi_move(double prev, double pi) { par[0] += std::log(prev / pi); set_tau(); }

  double q(const double* a, const double* zeta, double pi) const {
    double s = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : s) num_threads(threads)
#endif
    for (int64_t j = 0; j < M; ++j) {
      if (sup[j]) continue;
      const double x = lin(j, a); s += -.5 * pi * x - .5 * zeta[j] * std::exp(-x); }
    return s;
  }
  // Newton on the concave Q with step halving; returns the Newton steps taken
  int mstep(const double* zeta, double pi, int max_iter, double tol, int max_halving) {
    int it = 0;
    double q0 = q(par.data(), zeta, pi);
    for (; it < max_iter; ++it) {
      Eigen::MatrixXd H = Eigen::MatrixXd::Zero(K, K);
      Eigen::VectorXd g = Eigen::VectorXd::Zero(K);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
      {
        Eigen::MatrixXd Hl = Eigen::MatrixXd::Zero(K, K);
        Eigen::VectorXd gl = Eigen::VectorXd::Zero(K), a(K);
#ifdef _OPENMP
#pragma omp for nowait
#endif
        for (int64_t j = 0; j < M; ++j) {
          if (sup[j]) continue;
          const double e = zeta[j] * std::exp(-lin(j, par.data()));
          a[0] = 1.;
          const double* r = A.data() + j * (K - 1);
          for (int k = 1; k < K; ++k) a[k] = r[k - 1];
          gl.noalias() += (.5 * (e - pi)) * a;
          Hl.selfadjointView<Eigen::Lower>().rankUpdate(a, .5 * e);
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        { H += Hl; g += gl; }
      }
      const Eigen::VectorXd step = H.selfadjointView<Eigen::Lower>().ldlt().solve(g);
      if (!step.allFinite()) break;
      double t = 1.;
      std::vector<double> nw(K);
      bool moved = false;
      for (int h = 0; h < max_halving; ++h, t *= .5) {
        for (int k = 0; k < K; ++k) nw[k] = par[k] + t * step[k];
        const double q1 = q(nw.data(), zeta, pi);
        if (std::isfinite(q1) && q1 >= q0) { par = nw; q0 = q1; moved = true; break; }
      }
      if (!moved || t * step.cwiseAbs().maxCoeff() < tol) { ++it; break; }
    }
    set_tau();
    return it;
  }
};

// ---- one EM run with some hyperparameters fixed

struct RunResult {
  double pi, sigma, h2, elbo;
  int iters, mstep_steps;
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
    int msteps = 0;
    vtx(P, st.eta.data(), st.s.data(), threads);
    std::fill(st.active.begin(), st.active.end(), 1);
    double prev_elbo = -INFINITY, prev_sg = NAN, elbo = NAN, sigma_g = 0;
    int sg_count = 0, it = 0;
    bool have_prev_sg = false, force_full = true, fixed_sigma = false, ok = false;
    std::string msg = "max_iter reached";
    for (it = 1; it <= max_iter; ++it) {
      const bool full = !c.use_active || force_full || (it % c.full_every == 0);
      sweep(P, nv, bt, (float)pi, vm.tau.data(), vm.sup.data(), (float)sigma, st, !full, (float)c.act_tol, bmax, nvis, threads);
      for (int b = 0; b < P.nb; ++b) snp_visits += nvis[b];
      if (full) ++sweeps_full; else ++sweeps_part;
      force_full = false;
      if (it % c.refresh_every == 0) vtx(P, st.eta.data(), st.s.data(), threads);
      stats(M, P.d.data(), nv, bt, bv, st, P.S, vm.tau.data(), (float)sigma, out, zeta.data(), threads);
      const double sg_sum = out[0], sdz = out[2], eb = out[3], sde2 = out[4], ss = out[5], glg = out[6],
                   g1lg1 = out[7], glvt = out[8];
      msteps += vm.mstep(zeta.data(), pi, c.mstep_max_iter, c.mstep_tol, c.mstep_max_halving);
      // prior part of the ELBO with the updated variance model (as VIPRS 0.1.5 with the updated tau)
      double glt = 0, tz = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : glt, tz) num_threads(threads)
#endif
      for (int64_t j = 0; j < M; ++j) {
        glt += st.gamma[j] * std::log((double)vm.tau[j]);
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
    return RunResult{pi, sigma, sigma_g / (sigma_g + sigma), elbo, it, msteps, ok, fixed_sigma, msg};
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

// per annotation column: mean slab variance of its SNPs (weighted by the annotation) over the mean of all SNPs
// (suppressed SNPs left out)
Rcpp::NumericVector enrichment(const VarModel& vm) {
  const int K1 = vm.K - 1;
  std::vector<double> num(K1, 0.), den(K1, 0.);
  double all = 0;
  int64_t nfree = 0;
  for (int64_t j = 0; j < vm.M; ++j) {
    if (vm.sup[j]) continue;
    ++nfree;
    const double v = 1. / vm.tau[j];
    all += v;
    const double* r = vm.A.data() + j * K1;
    for (int k = 0; k < K1; ++k) { num[k] += r[k] * v; den[k] += r[k]; }
  }
  Rcpp::NumericVector e(K1);
  for (int k = 0; k < K1; ++k) e[k] = den[k] != 0 ? (num[k] / den[k]) / (all / nfree) : NA_REAL;
  return e;
}

template <typename T>
Rcpp::List pack(const Fitter<T>& F, const State& st, const RunResult& res, const VarModel& vm,
                const std::vector<double>& af, double h2_0, Rcpp::List extra) {
  const int64_t M = F.P.M;
  Rcpp::NumericVector beta(M), bstd(M), pip(M), var(M);
  for (int64_t j = 0; j < M; ++j) {
    const double p = std::min(std::max(af[j], 1e-6), 1 - 1e-6);
    bstd[j] = st.eta[j];
    beta[j] = st.eta[j] / std::sqrt(2. * p * (1. - p));
    pip[j] = st.gamma[j];
    var[j] = 1. / vm.tau[j];
  }
  Rcpp::List hyp = Rcpp::List::create(
      Rcpp::_["pi"] = res.pi, Rcpp::_["sigma"] = res.sigma, Rcpp::_["h2"] = res.h2,
      Rcpp::_["elbo"] = res.elbo, Rcpp::_["iters"] = res.iters, Rcpp::_["mstep_steps"] = res.mstep_steps,
      Rcpp::_["converged"] = res.converged, Rcpp::_["message"] = res.message, Rcpp::_["sigma_fixed"] = res.sigma_fixed,
      Rcpp::_["h2_ldsc_init"] = h2_0, Rcpp::_["n_snps"] = (double)M, Rcpp::_["sum_k"] = (double)F.P.S,
      Rcpp::_["sweeps_full"] = (double)F.sweeps_full, Rcpp::_["sweeps_partial"] = (double)F.sweeps_part,
      Rcpp::_["snp_visits_per_M"] = (double)F.snp_visits / M);
  return Rcpp::List::create(Rcpp::_["BETA"] = beta, Rcpp::_["BETA_STD"] = bstd, Rcpp::_["PIP"] = pip,
                            Rcpp::_["VAR"] = var, Rcpp::_["coef"] = Rcpp::NumericVector(vm.par.begin(), vm.par.end()),
                            Rcpp::_["enrichment"] = enrichment(vm), Rcpp::_["hyp"] = hyp, Rcpp::_["extra"] = extra);
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
                    const std::string& ld_dir, double thresh, const ViprsControl& c, int threads) {
  const EigenPanel<T> P = plink4r::load_eigen_panel<T>(ld_dir, in.blocks, in.keep, in.bstd, thresh, threads);
  Fitter<T> F(P, c, in.n, in.z, threads);
  const double h2_0 = F.h2_ldsc();
  const std::vector<double> grid = F.pi_grid();
  std::vector<float> train, valid, n_tr;
  std::vector<double> xi;
  F.pumas(train, valid, n_tr, xi);
  std::vector<double> gpi, gsig, gh2, gr2, grb, gsec;
  std::vector<int> git, gms, gconv;
  Rcpp::NumericMatrix gcoef(grid.size(), K);
  double best_r2 = -INFINITY, bpi = 0, bsig = 0, best_here = -INFINITY;
  State st(P.M, P.S, grid[0]), best(0, 0, 0.);
  for (int64_t j = 0; j < P.M; ++j) if (sup[j]) st.gamma[j] = 0.f;
  VarModel vm(P.M, K, A, sup, threads), bvm(P.M, K, A, sup, threads);
  bool have_best = false;
  vm.start(h2_0, grid[0]);
  double sigma = 1. - h2_0, prev_pi = grid[0];
  int drops = 0;
  for (double pi : grid) {
    vm.pi_move(prev_pi, pi);
    const auto t1 = std::chrono::steady_clock::now();
    const RunResult r = F.run(st, n_tr.data(), train.data(), pi, vm, sigma, c.grid_max_iter, valid.data());
    sigma = r.sigma; prev_pi = pi;
    double r2, rb;
    Fitter<T>::pseudo_r2(st, valid, r2, rb);
    for (int k = 0; k < K; ++k) gcoef(gpi.size(), k) = vm.par[k];
    gpi.push_back(pi); gsig.push_back(sigma); gh2.push_back(r.h2); gr2.push_back(r2);
    grb.push_back(rb); git.push_back(r.iters); gms.push_back(r.mstep_steps); gconv.push_back(r.converged);
    gsec.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count());
    if (r.converged && rb > 0 && r2 > best_r2) {
      best_r2 = r2; best = st; have_best = true;
      bpi = pi; bsig = sigma; bvm.par = vm.par; bvm.tau = vm.tau;
    }
    if (r2 > best_here) { best_here = r2; drops = 0; }
    else if (c.stop_after && ++drops >= c.stop_after) break;
  }
  if (!have_best) Rcpp::stop("no grid point converged with a positive pseudo R2");
  RunResult res = F.run(best, in.n.data(), P.bt.data(), bpi, bvm, bsig, c.max_iter, nullptr);
  const int ng = gpi.size();
  Rcpp::DataFrame tab = Rcpp::DataFrame::create(
      Rcpp::_["pi"] = gpi, Rcpp::_["sigma"] = gsig, Rcpp::_["h2"] = gh2, Rcpp::_["pseudo_r2"] = gr2,
      Rcpp::_["rb"] = grb, Rcpp::_["iters"] = git, Rcpp::_["mstep_steps"] = gms, Rcpp::_["converged"] = gconv,
      Rcpp::_["seconds"] = gsec);
  return pack(F, best, res, bvm, in.af, h2_0,
              Rcpp::List::create(Rcpp::_["grid"] = tab, Rcpp::_["grid_coef"] = gcoef(Rcpp::Range(0, ng - 1), Rcpp::_),
                                 Rcpp::_["best_pseudo_r2"] = best_r2, Rcpp::_["n_bins"] = n_bins,
                                 Rcpp::_["suppressed"] = Rcpp::LogicalVector(sup.begin(), sup.end()),
                                 Rcpp::_["pumas_xi"] = Rcpp::NumericVector(xi.begin(), xi.end())));
}

}  // namespace

// Entry: blocks (ids), keep (0-based eigen rows with data per block, increasing), z / n / af per kept SNP (blocks
// concatenated in that order; af = frequency of the effect allele, for BETA = BETA_STD / sqrt(2 af (1 - af))),
// A (kept SNPs x annotations, same order; the intercept is added here), bin_file (LD-bin list, "" = no suppression),
// bin_start (first panel row of each kept block), n_panel (SNPs of the panel), thresh (eigen truncation, <= 0 = every
// stored pair), fp16 (store V as float16), control (viprs_annot_control()).

// [[Rcpp::export]]
Rcpp::List viprs_annot_cpp(const std::string& ld_dir, Rcpp::IntegerVector blocks, Rcpp::List keep,
                           Rcpp::NumericVector z, Rcpp::NumericVector n, Rcpp::NumericVector af, Rcpp::NumericMatrix A,
                           const std::string& bin_file, Rcpp::NumericVector bin_start, double n_panel, double thresh,
                           bool fp16, Rcpp::List control, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const ViprsControl c = read_control(control);
  const Inputs in = read_inputs(blocks, keep, z, n, af);
  const int64_t M = in.z.size();
  if (A.nrow() != M) Rcpp::stop("A must have one row per kept SNP");
  const int K = A.ncol() + 1;
  std::vector<double> Ar(M * (K - 1));
  for (int64_t j = 0; j < M; ++j)
    for (int k = 0; k < K - 1; ++k) {
      const double x = A(j, k);
      if (!std::isfinite(x)) Rcpp::stop("A must be finite");
      Ar[j * (K - 1) + k] = x;
    }
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
    return fp16 ? fit_path<uint16_t>(in, Ar, K, sup, n_bins, ld_dir, thresh, c, T)
                : fit_path<float>(in, Ar, K, sup, n_bins, ld_dir, thresh, c, T);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}
