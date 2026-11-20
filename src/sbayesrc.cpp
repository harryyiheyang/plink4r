// SBayesRC on the eigen LD (top layer: the package's SBayesRC functional unit, built on eigen_int16.h).
// Ported from fastSBayesRC (Yihe Yang; sbrc_mc.py, sbrc_vi.py, fast_chain.py, kernel.cpp, 2026-11-10/11).
//
// Model (SBayesRC, Zheng et al. 2024): per LD block w = X' beta + e, X = U Lambda^{1/2} = Q diag(c) (int16 rows),
// w = Lambda^{-1} X' bhat, e ~ N(0, I / n); beta_j ~ sum_k pi_jk N(0, gam_k sigma^2), gam_0 = 0. Annotations enter the
// prior weights through the nested logistic structure P(delta_j >= k | delta_j >= k-1) = logistic(f_j' theta_k),
// f_j = (1, A_j), with the adaptive spike-and-slab lasso on the slopes (lam0 = ratio lam1).
// Fit (sbrc_mc.fit): per iteration one Gibbs sweep per block (blocks in parallel; Rao-Blackwell means kept), sigma^2
// from its scaled inverse chi^2, theta by stochastic-approximation EM on the soft counts (Newton / coordinate descent
// with soft thresholding on the q x q quadratic model; Hessian refreshed every hess_every iterations).
// Pipeline (fast_chain.py): SBayesRC pseudo training / validation summaries, warm-start chain over eigen cutoffs,
// SBayesRC rule for the cutoff, main fit on the full data warm from the chosen level.
// LD-bin suppression (ld_bins.h): inside a bin of r2 > bin_r2 with an annotated member, SNPs without annotation get
// delta = 0 (no prior mass on the nonzero classes) in every fit and are left out of the theta regression.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_int16.h"
#include "ld_bins.h"
#include "resources.h"

using plink4r::Block16;

namespace {

// ---------------------------------------------------------------- random numbers
// Own streams (mt19937_64), or a replay of recorded draws (validation against the original, which draws with numpy).
struct Rng {
  virtual ~Rng() {}
  virtual double normal() = 0;
  virtual unsigned long long block_seed() = 0;  // in [1, 2^62)
  virtual double chisq(double df) = 0;
  virtual void restart_fit() = 0;  // a new fit starts its own stream (sbrc_mc.fit: default_rng(seed) per call)
};

struct RngStd : Rng {
  std::mt19937_64 pseudo, fit;
  unsigned long long fit_seed;
  RngStd(unsigned long long pseudo_seed, unsigned long long fit_seed_) : pseudo(pseudo_seed), fit(fit_seed_), fit_seed(fit_seed_) {}
  double normal() override { return std::normal_distribution<double>()(pseudo); }
  unsigned long long block_seed() override {
    return std::uniform_int_distribution<unsigned long long>(1ull, (1ull << 62) - 1)(fit);
  }
  double chisq(double df) override { return std::chi_squared_distribution<double>(df)(fit); }
  void restart_fit() override { fit.seed(fit_seed); }
};

struct RngReplay : Rng {
  std::vector<double> nrm, chi;
  std::vector<unsigned long long> sd;
  std::size_t in = 0, ic = 0, is = 0;
  double normal() override {
    if (in >= nrm.size()) throw std::runtime_error("replay: normals exhausted");
    return nrm[in++];
  }
  unsigned long long block_seed() override {
    if (is >= sd.size()) throw std::runtime_error("replay: seeds exhausted");
    return sd[is++];
  }
  double chisq(double) override {
    if (ic >= chi.size()) throw std::runtime_error("replay: chi-square draws exhausted");
    return chi[ic++];
  }
  void restart_fit() override {}
};

// ---------------------------------------------------------------- int16 row kernels (kernel.cpp)

inline float dot16(const int16_t* q, const float* rc, int k) {
  float acc = 0.f;
#pragma omp simd reduction(+ : acc)
  for (int t = 0; t < k; t++) acc += (float)q[t] * rc[t];
  return acc;
}

inline void axpy16(const int16_t* q, const float* c2, float* rc, float f, int k) {
#pragma omp simd
  for (int t = 0; t < k; t++) rc[t] -= c2[t] * (float)q[t] * f;
}

// d_j = sum_t (q_jt c_t)^2 over the first k columns
void diag16(int m, int ld, int k, const int16_t* Q, const float* c2, double* d) {
  for (int j = 0; j < m; j++) {
    const int16_t* q = Q + (size_t)j * ld;
    float acc = 0.f;
#pragma omp simd reduction(+ : acc)
    for (int t = 0; t < k; t++) acc += (float)q[t] * (float)q[t] * c2[t];
    d[j] = acc;
  }
}

// One Gibbs sweep of one block: sample delta_j ~ phi_j and beta_j | delta_j, residual updated with the sampled beta
// (rc = c * r, scaled); rb gets the Rao-Blackwell mean, phi the class probabilities, cnt the class counts,
// ssb += sum beta^2 / gam.
void gibbs_block16(int m, int ld, int k, int K, const int16_t* Q, const float* c2, const double* d, float* rc, double* mb,
                   const double* logpi, double s2, const double* gam, double prec, double* rb, double* cnt, double* ssb,
                   double* phi, unsigned long long seed) {
  std::mt19937_64 eng(seed);
  std::uniform_real_distribution<double> U01(0.0, 1.0);
  std::normal_distribution<double> N01;
  std::vector<double> lw(K), mu(K), v(K);
  for (int j = 0; j < m; j++) {
    const int16_t* q = Q + (size_t)j * ld;
    const double dj = d[j];
    const double xr = (double)dot16(q, rc, k) + dj * mb[j];
    const double* lp = logpi + (size_t)j * K;
    double mx = -INFINITY;
    for (int c = 0; c < K; c++) {
      if (gam[c] == 0) { lw[c] = lp[c]; mu[c] = 0; v[c] = 0; }
      else {
        const double g = gam[c] * s2;
        v[c] = 1.0 / (prec * dj + 1.0 / g);
        mu[c] = v[c] * prec * xr;
        lw[c] = lp[c] + 0.5 * std::log(v[c] / g) + 0.5 * mu[c] * mu[c] / v[c];
      }
      mx = std::max(mx, lw[c]);
    }
    double sw = 0;
    for (int c = 0; c < K; c++) { lw[c] = std::exp(lw[c] - mx); sw += lw[c]; }
    double u = U01(eng) * sw, acc = 0, rbm = 0;
    int cs = K - 1;
    for (int c = 0; c < K; c++) { rbm += lw[c] / sw * mu[c]; phi[(size_t)j * K + c] = lw[c] / sw; }
    for (int c = 0; c < K; c++) { acc += lw[c]; if (u <= acc) { cs = c; break; } }
    const double bnew = gam[cs] == 0 ? 0.0 : mu[cs] + std::sqrt(v[cs]) * N01(eng);
    rb[j] = rbm;
    cnt[cs] += 1;
    if (gam[cs] > 0) ssb[0] += bnew * bnew / gam[cs];
    const double dm = bnew - mb[j];
    if (dm != 0) { axpy16(q, c2, rc, (float)dm, k); mb[j] = bnew; }
  }
}

// ---------------------------------------------------------------- annotation prior (logistic nested structure)

struct Ctl {
  std::vector<double> gam, pi0;
  double ve, het_floor, lam_eps, vg_floor, ridge_prec, ssl_ratio, ssl_l1, nu, eta_clip, w_floor, h_floor, pen_intercept, cd_tol, irls_tol;
  int warm, newton_steps, hess_every, slope_start, cd_maxit;
  std::vector<double> ssl_a, ssl_b;  // Beta mixture prior of the SSL slab weight
  double ssl_gs, ssl_gr, ssl_th_lo, ssl_th_hi;
};

// log pi (M x K): pi_j0 = 1 - p_1, pi_jk = p_1..p_k (1 - p_{k+1}), pi_j,K-1 = p_1..p_{K-1}, p_k = logistic(f_j' th_k)
void logit_logpi(int M, int K, int q, const double* F, const double* th, double clip, double* out, int threads) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
  for (int j = 0; j < M; j++) {
    const double* f = F + (size_t)j * q;
    double* o = out + (size_t)j * K;
    double cum = 0;
    for (int k = 0; k < K - 1; k++) {
      double eta = 0;
      const double* t = th + (size_t)k * q;
      for (int c = 0; c < q; c++) eta += f[c] * t[c];
      eta = std::min(std::max(eta, -clip), clip);
      const double lup = eta > 0 ? -std::log1p(std::exp(-eta)) : eta - std::log1p(std::exp(eta));
      o[k] = cum + lup - eta;
      cum += lup;
    }
    o[K - 1] = cum;
  }
}

// gradient g (L x q) and, with want_G, Gram G (L x q x q) of sum_j s_lj log p_lj + f_lj log(1 - p_lj) at th (L x q)
void logit_stats_multi(int M, int q, int L, const double* F, const double* th, const double* s, const double* fl,
                       bool want_G, double clip, double w_floor, double h_floor, double* g, double* G, int threads) {
  std::fill(g, g + (size_t)L * q, 0.0);
  if (want_G) std::fill(G, G + (size_t)L * q * q, 0.0);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
  {
    std::vector<double> gl((size_t)L * q, 0.0), Gl(want_G ? (size_t)L * q * q : 0, 0.0);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int j = 0; j < M; j++) {
      const double* f = F + (size_t)j * q;
      for (int l = 0; l < L; l++) {
        const double sj = s[(size_t)l * M + j], w = sj + fl[(size_t)l * M + j];
        if (w <= w_floor) continue;
        const double* t = th + (size_t)l * q;
        double eta = 0;
        for (int c = 0; c < q; c++) eta += f[c] * t[c];
        eta = std::min(std::max(eta, -clip), clip);
        const double p = 1.0 / (1.0 + std::exp(-eta)), gj = sj - w * p;
        double* gg = &gl[(size_t)l * q];
        for (int a = 0; a < q; a++) gg[a] += f[a] * gj;
        if (want_G) {
          const double hj = std::max(w * p * (1 - p), h_floor);
          double* GG = &Gl[(size_t)l * q * q];
          for (int a = 0; a < q; a++) {
            const double ha = hj * f[a];
            for (int b = a; b < q; b++) GG[(size_t)a * q + b] += ha * f[b];
          }
        }
      }
    }
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      for (size_t i = 0; i < gl.size(); i++) g[i] += gl[i];
      for (size_t i = 0; i < Gl.size(); i++) G[i] += Gl[i];
    }
  }
  if (want_G)
    for (int l = 0; l < L; l++) {
      double* GG = G + (size_t)l * q * q;
      for (int a = 0; a < q; a++)
        for (int b = 0; b < a; b++) GG[(size_t)a * q + b] = GG[(size_t)b * q + a];
    }
}

// max_b g'(b - th) - 1/2 (b - th)' G (b - th) - sum lam_c |b_c| by coordinate descent with soft thresholding
std::vector<double> cd_quad(int q, const std::vector<double>& G, const std::vector<double>& g, const std::vector<double>& th,
                            const std::vector<double>& lam, int maxit, double tol) {
  std::vector<double> r(q), b(q);
  for (int c = 0; c < q; c++) { b[c] = th[c]; r[c] = g[c]; }
  for (int it = 0; it < maxit; it++) {
    double dmax = 0;
    for (int c = 0; c < q; c++) {
      const double Gcc = G[(size_t)c * q + c];
      const double z = b[c] + r[c] / Gcc;
      const double nb = lam[c] > 0 ? (z > 0 ? 1 : -1) * std::max(std::fabs(z) - lam[c] / Gcc, 0.0) : z;
      const double dl = nb - b[c];
      if (dl != 0) {
        for (int a = 0; a < q; a++) r[a] -= G[(size_t)a * q + c] * dl;
        b[c] = nb;
        dmax = std::max(dmax, std::fabs(dl));
      }
    }
    if (dmax < tol) break;
  }
  return b;
}

// intercept-only level M-step (sbrc_vi._irls_cd with one column of ones): Newton steps with the ridge pen
double irls_intercept(int M, const double* s, const double* fl, double th, double pen, const Ctl& c, int threads) {
  const std::vector<double> ones(M, 1.0);
  for (int o = 0; o < c.newton_steps; ++o) {
    std::vector<double> g0(1), G0(1);
    logit_stats_multi(M, 1, 1, ones.data(), &th, s, fl, true, c.eta_clip, c.w_floor, c.h_floor, g0.data(), G0.data(),
                      threads);
    const std::vector<double> b = cd_quad(1, {G0[0] + pen}, {g0[0] - pen * th}, {th}, {0.0}, c.cd_maxit, c.cd_tol);
    const double step = b[0] - th;
    th = b[0];
    if (std::fabs(step) < c.irls_tol) break;
  }
  return th;
}

// intercepts of the nested regressions that reproduce a global pi
std::vector<double> theta_from_pi(const std::vector<double>& pi) {
  const int K = (int)pi.size();
  std::vector<double> th(K - 1);
  for (int k = 1; k < K; ++k) {
    double a = 0, b = 0;
    for (int l = k; l < K; ++l) a += pi[l];
    for (int l = k - 1; l < K; ++l) b += pi[l];
    const double p = std::min(std::max(a / b, 1e-10), 1 - 1e-10);
    th[k - 1] = std::log(p / (1 - p));
  }
  return th;
}

// adaptive spike-and-slab lasso on the slopes of one level (featprs): lam0 = ratio lam1, lam1 ~ Gamma MAP, slab weight
// ~ equal mixture of Beta priors
struct SSL {
  double ratio, l1, th;
  const Ctl* c;
  SSL(const Ctl& c_) : ratio(c_.ssl_ratio), l1(c_.ssl_l1), c(&c_) {
    double m = 0;
    for (std::size_t i = 0; i < c->ssl_a.size(); ++i) m += c->ssl_a[i] / (c->ssl_a[i] + c->ssl_b[i]);
    th = m / c->ssl_a.size();
  }
  double pstar(double a) const {
    const double l0 = ratio * l1;
    const double ls = std::log(th) + std::log(l1 / 2) - l1 * std::fabs(a);
    const double lp = std::log(1 - th) + std::log(l0 / 2) - l0 * std::fabs(a);
    return 1 / (1 + std::exp(lp - ls));
  }
  double weight(double a) const {
    const double p = pstar(a);
    return p * l1 + (1 - p) * ratio * l1;
  }
  void update(const std::vector<double>& a) {
    const double K = (double)a.size();
    double s = 0, pa = 0;
    for (double x : a) { const double p = pstar(x); s += p; pa += p * std::fabs(x); }
    const int nb = (int)c->ssl_a.size();
    std::vector<double> lml(nb);
    double mx = -INFINITY;
    for (int i = 0; i < nb; ++i) {
      const double x = c->ssl_a[i], y = c->ssl_b[i];
      lml[i] = std::lgamma(x + s) + std::lgamma(y + K - s) - std::lgamma(x + y + K) - std::lgamma(x) - std::lgamma(y) +
               std::lgamma(x + y);
      mx = std::max(mx, lml[i]);
    }
    double sw = 0, e = 0;
    for (int i = 0; i < nb; ++i) { lml[i] = std::exp(lml[i] - mx); sw += lml[i]; }
    for (int i = 0; i < nb; ++i) e += lml[i] / sw * ((c->ssl_a[i] + s) / (c->ssl_a[i] + c->ssl_b[i] + K));
    th = std::min(std::max(e, c->ssl_th_lo), c->ssl_th_hi);
    l1 = (s + c->ssl_gs - 1) / (pa + c->ssl_gr);
  }
};

// ---------------------------------------------------------------- one fit (sbrc_mc.fit, int16 path, qn M-step)

struct FitState {  // warm start carried between fits (cutoff chain)
  std::vector<double> beta_last, theta;  // M; (K-1) x q
  double s2 = 0;
  std::vector<double> ssl_l1, ssl_th;
};

struct FitOut {
  std::vector<double> beta;  // Rao-Blackwell posterior mean after burn-in
  FitState st;
  std::vector<double> tr_vg, tr_nnz, tr_s2;
};

FitOut fit(const std::vector<Block16>& P, const std::vector<int64_t>& off, const std::vector<int>& kact,
           const std::vector<std::vector<double>>& W, double n, const std::vector<double>& F, int q, int niter, int burn,
           const FitState* init, const std::vector<char>& sup, const Ctl& c, Rng& rng, int threads) {
  const int nb = (int)P.size();
  const int64_t M = off[nb];
  const bool any_sup = !sup.empty();
  int64_t M_free = M;
  if (any_sup) for (char x : sup) M_free -= x;
  const int K = (int)c.gam.size(), L = K - 1;
  std::vector<bool> nz(K);
  for (int k = 0; k < K; ++k) nz[k] = c.gam[k] > 0;
  std::vector<bool> slopes(q, false);
  for (int a = 1; a < q; ++a) slopes[a] = true;
  const bool use_ssl = q > 1;
  std::vector<SSL> ssls;
  if (use_ssl) for (int l = 0; l < L; ++l) ssls.emplace_back(c);
  std::vector<double> th((size_t)L * q, 0.0);
  const std::vector<double> th0 = theta_from_pi(c.pi0);
  for (int l = 0; l < L; ++l) th[(size_t)l * q] = th0[l];
  const double prec = n / c.ve;
  double vg0 = 0, pg = 0;
  for (int b = 0; b < nb; ++b) {
    double ww = 0;
    for (int t = 0; t < kact[b]; ++t) ww += W[b][t] * W[b][t];
    vg0 += std::max(ww - kact[b] / prec, 0.0);
  }
  for (int k = 0; k < K; ++k) pg += c.pi0[k] * c.gam[k];
  double s2 = std::max(vg0, c.vg_floor) / std::max<int64_t>(M_free, 1) / pg;
  const bool has_init = init != nullptr;
  if (has_init) {
    th = init->theta;
    s2 = init->s2;
    for (int l = 0; l < (int)ssls.size() && l < (int)init->ssl_l1.size(); ++l) { ssls[l].l1 = init->ssl_l1[l]; ssls[l].th = init->ssl_th[l]; }
  }
  // per-block state: scaled residual rc = c r (float), d, current beta
  std::vector<std::vector<float>> rc(nb), c2(nb), cc(nb);
  std::vector<std::vector<double>> d(nb), mb(nb);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int b = 0; b < nb; ++b) {
    const Block16& B = P[b];
    const int ka = kact[b];
    cc[b].assign(B.c.begin(), B.c.begin() + ka);
    c2[b].resize(ka);
    for (int t = 0; t < ka; ++t) c2[b][t] = cc[b][t] * cc[b][t];
    d[b].resize(B.nk);
    diag16(B.nk, B.ld, ka, B.Q.data(), c2[b].data(), d[b].data());
    mb[b].assign(B.nk, 0.0);
    bool any = false;
    if (has_init)
      for (int j = 0; j < B.nk; ++j) { mb[b][j] = init->beta_last[off[b] + j]; any = any || mb[b][j] != 0; }
    std::vector<double> r(W[b].begin(), W[b].begin() + ka);
    if (any) {  // r = w - c * (Q' m) (float32 product, as the original)
      std::vector<float> acc(ka, 0.f);
      for (int j = 0; j < B.nk; ++j) {
        const float mj = (float)mb[b][j];
        if (mj == 0.f) continue;
        const int16_t* qj = B.Q.data() + (size_t)j * B.ld;
        for (int t = 0; t < ka; ++t) acc[t] += (float)qj[t] * mj;
      }
      for (int t = 0; t < ka; ++t) r[t] -= (double)cc[b][t] * acc[t];
    }
    rc[b].resize(ka);
    for (int t = 0; t < ka; ++t) rc[b][t] = (float)(cc[b][t] * r[t]);
  }
  std::vector<double> logpi((size_t)M * K), phi((size_t)M * K, 0.0), rb(M, 0.0), bsum(M, 0.0);
  std::vector<double> Sg((size_t)L * M, 0.0), Fs((size_t)L * M, 0.0), gA((size_t)L * q), GA((size_t)L * q * q);
  std::vector<std::vector<double>> Hc;
  int nsum = 0;
  bool first = true;
  FitOut out;
  rng.restart_fit();
  std::vector<double> cnts((size_t)nb * K), ssbs(nb);
  std::vector<unsigned long long> seeds(nb);
  for (int it = 1; it <= niter; ++it) {
    logit_logpi((int)M, K, q, F.data(), th.data(), c.eta_clip, logpi.data(), threads);
    if (any_sup)  // suppressed SNPs: all prior mass on the zero class
      for (int64_t j = 0; j < M; ++j)
        if (sup[j]) {
          double* lp = logpi.data() + (size_t)j * K;
          for (int k = 0; k < K; ++k) lp[k] = c.gam[k] == 0 ? 0.0 : -INFINITY;
        }
    std::fill(cnts.begin(), cnts.end(), 0.0);
    std::fill(ssbs.begin(), ssbs.end(), 0.0);
    for (int b = 0; b < nb; ++b) seeds[b] = rng.block_seed();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
    for (int b = 0; b < nb; ++b) {
      const Block16& B = P[b];
      const int64_t i0 = off[b];
      gibbs_block16(B.nk, B.ld, kact[b], K, B.Q.data(), c2[b].data(), d[b].data(), rc[b].data(), mb[b].data(),
                    logpi.data() + (size_t)i0 * K, s2, c.gam.data(), prec, rb.data() + i0, cnts.data() + (size_t)b * K,
                    &ssbs[b], phi.data() + (size_t)i0 * K, seeds[b]);
    }
    double nnz = 0, ssb = 0;
    for (int b = 0; b < nb; ++b) {
      for (int k = 0; k < K; ++k) if (nz[k]) nnz += cnts[(size_t)b * K + k];
      ssb += ssbs[b];
    }
    // sigma^2 and SAEM M-step for theta
    s2 = (ssb + c.nu * s2) / rng.chisq(nnz + c.nu);
    const double a = (it <= burn || first) ? 1.0 : 1.0 / (it - burn + 1);
    first = false;
    const bool slopes_on = it > c.slope_start || has_init;
    const bool in_warm = use_ssl && it <= c.slope_start + c.warm && !has_init;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (int64_t j = 0; j < M; ++j) {
      if (any_sup && sup[j]) continue;  // left out of the regression (soft counts stay 0)
      const double* p = phi.data() + (size_t)j * K;
      double up = 0;
      for (int l = K - 1; l >= 1; --l) {
        up += p[l];
        Sg[(size_t)(l - 1) * M + j] = (1 - a) * Sg[(size_t)(l - 1) * M + j] + a * up;
        Fs[(size_t)(l - 1) * M + j] = (1 - a) * Fs[(size_t)(l - 1) * M + j] + a * p[l - 1];
      }
    }
    const bool fused = slopes_on || q == 1;
    if (fused) {
      const bool want = Hc.empty() || it % c.hess_every == 0;
      logit_stats_multi((int)M, q, L, F.data(), th.data(), Sg.data(), Fs.data(), want, c.eta_clip, c.w_floor, c.h_floor,
                        gA.data(), GA.data(), threads);
      if (want) {
        Hc.assign(L, std::vector<double>((size_t)q * q));
        for (int l = 0; l < L; ++l) std::copy(GA.begin() + (size_t)l * q * q, GA.begin() + (size_t)(l + 1) * q * q, Hc[l].begin());
      }
    }
    for (int k = 1; k < K; ++k) {
      SSL* sk = (use_ssl && !in_warm) ? &ssls[k - 1] : nullptr;
      std::vector<double> pk(q);
      for (int a2 = 0; a2 < q; ++a2) pk[a2] = (slopes[a2] && in_warm) ? c.ridge_prec : (a2 == 0 ? c.pen_intercept : 0.0);
      double* tk = th.data() + (size_t)(k - 1) * q;
      if (fused) {
        std::vector<double> t0(tk, tk + q), lam(q, 0.0), G = Hc[k - 1], g(q);
        if (sk) for (int a2 = 1; a2 < q; ++a2) lam[a2] = sk->weight(t0[a2]);
        for (int a2 = 0; a2 < q; ++a2) { G[(size_t)a2 * q + a2] += pk[a2]; g[a2] = gA[(size_t)(k - 1) * q + a2] - pk[a2] * t0[a2]; }
        const std::vector<double> nb2 = cd_quad(q, G, g, t0, lam, c.cd_maxit, c.cd_tol);
        std::copy(nb2.begin(), nb2.end(), tk);
      } else {
        tk[0] = irls_intercept((int)M, Sg.data() + (size_t)(k - 1) * M, Fs.data() + (size_t)(k - 1) * M, tk[0], pk[0], c,
                               threads);
      }
      if (sk && it > c.slope_start) sk->update(std::vector<double>(tk + 1, tk + q));
    }
    double vg = 0;
    for (int b = 0; b < nb; ++b)
      for (int t = 0; t < kact[b]; ++t) {
        const double x = W[b][t] - (cc[b][t] > 0 ? rc[b][t] / (double)cc[b][t] : 0.0);
        vg += x * x;
      }
    out.tr_vg.push_back(vg); out.tr_nnz.push_back(nnz); out.tr_s2.push_back(s2);
    if (it > burn) {
      for (int64_t j = 0; j < M; ++j) bsum[j] += rb[j];
      ++nsum;
    }
  }
  out.beta.resize(M);
  for (int64_t j = 0; j < M; ++j) out.beta[j] = bsum[j] / std::max(nsum, 1);
  out.st.beta_last.resize(M);
  for (int b = 0; b < nb; ++b) std::copy(mb[b].begin(), mb[b].end(), out.st.beta_last.begin() + off[b]);
  out.st.theta = th;
  out.st.s2 = s2;
  for (auto& x : ssls) { out.st.ssl_l1.push_back(x.l1); out.st.ssl_th.push_back(x.th); }
  return out;
}

// ---------------------------------------------------------------- panel-level helpers (fast_chain.py)

// w_b = c * (Q' b_b) / lambda (lambda <= 1e-8 -> 0), Q' b in float32
std::vector<std::vector<double>> wfull(const std::vector<Block16>& P, const std::vector<int64_t>& off,
                                       const std::vector<double>& bh, double lam_eps, int threads) {
  const int nb = (int)P.size();
  std::vector<std::vector<double>> W(nb);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int b = 0; b < nb; ++b) {
    const Block16& B = P[b];
    std::vector<float> acc(B.ld, 0.f);
    for (int j = 0; j < B.nk; ++j) {
      const float x = (float)bh[off[b] + j];
      const int16_t* qj = B.Q.data() + (size_t)j * B.ld;
      for (int t = 0; t < B.ld; ++t) acc[t] += (float)qj[t] * x;
    }
    W[b].resize(B.ld);
    for (int t = 0; t < B.ld; ++t) W[b][t] = B.lambda[t] > lam_eps ? (double)B.c[t] * acc[t] / B.lambda[t] : 0.0;
  }
  return W;
}

// beta' R beta on the full stored panel = sum_b ||X_b' beta_b||^2
double prs_var(const std::vector<Block16>& P, const std::vector<int64_t>& off, const std::vector<double>& beta,
               int threads) {
  const int nb = (int)P.size();
  double tot = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : tot) num_threads(threads)
#endif
  for (int b = 0; b < nb; ++b) {
    const Block16& B = P[b];
    std::vector<float> acc(B.ld, 0.f);
    for (int j = 0; j < B.nk; ++j) {
      const float x = (float)beta[off[b] + j];
      if (x == 0.f) continue;
      const int16_t* qj = B.Q.data() + (size_t)j * B.ld;
      for (int t = 0; t < B.ld; ++t) acc[t] += (float)qj[t] * x;
    }
    double s = 0;
    for (int t = 0; t < B.ld; ++t) { const double v = (double)B.c[t] * acc[t]; s += v * v; }
    tot += s;
  }
  return tot;
}

// eigen components kept at cutoff cut: all stored when cut >= file cutoff, else the fewest leading ones with
// cumsum(lambda) >= cut * sumLambda
std::vector<int> kcut(const std::vector<Block16>& P, double cut) {
  std::vector<int> k(P.size());
  for (std::size_t b = 0; b < P.size(); ++b) {
    const Block16& B = P[b];
    if (cut >= B.file_cut - 1e-9) { k[b] = B.ld; continue; }
    double cum = 0;
    int i = 0;
    for (; i < B.ld; ++i) { cum += B.lambda[i]; if (cum >= cut * B.sum_lambda) break; }
    k[b] = std::min(B.ld, i + 1);
  }
  return k;
}

Ctl read_ctl(const Rcpp::List& x) {
  Ctl c;
  c.gam = Rcpp::as<std::vector<double>>(x["gam"]);
  c.pi0 = Rcpp::as<std::vector<double>>(x["pi0"]);
  c.ve = Rcpp::as<double>(x["ve"]);
  c.het_floor = Rcpp::as<double>(x["het_floor"]);
  c.lam_eps = Rcpp::as<double>(x["lam_eps"]);
  c.vg_floor = Rcpp::as<double>(x["vg_floor"]);
  c.ridge_prec = Rcpp::as<double>(x["ridge_prec"]);
  c.ssl_ratio = Rcpp::as<double>(x["ssl_ratio"]);
  c.ssl_l1 = Rcpp::as<double>(x["ssl_l1"]);
  c.nu = Rcpp::as<double>(x["nu"]);
  c.eta_clip = Rcpp::as<double>(x["eta_clip"]);
  c.w_floor = Rcpp::as<double>(x["w_floor"]);
  c.h_floor = Rcpp::as<double>(x["h_floor"]);
  c.pen_intercept = Rcpp::as<double>(x["pen_intercept"]);
  c.cd_tol = Rcpp::as<double>(x["cd_tol"]);
  c.irls_tol = Rcpp::as<double>(x["irls_tol"]);
  c.warm = Rcpp::as<int>(x["warm"]);
  c.newton_steps = Rcpp::as<int>(x["newton_steps"]);
  c.hess_every = Rcpp::as<int>(x["hess_every"]);
  c.slope_start = Rcpp::as<int>(x["slope_start"]);
  c.cd_maxit = Rcpp::as<int>(x["cd_maxit"]);
  Rcpp::NumericMatrix sb = x["ssl_beta"];
  for (int i = 0; i < sb.nrow(); ++i) { c.ssl_a.push_back(sb(i, 0)); c.ssl_b.push_back(sb(i, 1)); }
  Rcpp::NumericVector sg = x["ssl_gamma"];
  c.ssl_gs = sg[0]; c.ssl_gr = sg[1];
  Rcpp::NumericVector sr = x["ssl_th_range"];
  c.ssl_th_lo = sr[0]; c.ssl_th_hi = sr[1];
  if (c.gam.size() < 2 || c.gam.size() != c.pi0.size()) Rcpp::stop("gam and pi0 must have the same length >= 2");
  return c;
}

}  // namespace

// SBayesRC pipeline of fast_chain.py. blocks / keep: panel blocks and their 0-based eigen rows with data; z, n: z-score
// and sample size per kept SNP (blocks concatenated): bhat = z / sqrt(n + z^2) (summarystat_impute BETA), N = median n;
// af: effect-allele frequency; A: kept SNPs x C annotation matrix (C may be 0); control: sbayesrc_control();
// replay: NULL, or list(normal, seed (character), chisq) of recorded draws (validation against the original only).
// [[Rcpp::export]]
Rcpp::List sbayesrc_cpp(const std::string& ld_dir, Rcpp::IntegerVector blocks, Rcpp::List keep, Rcpp::NumericVector z,
                        Rcpp::NumericVector n, Rcpp::NumericVector af, Rcpp::NumericMatrix A, Rcpp::List control,
                        Rcpp::Nullable<Rcpp::List> replay, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const Ctl c = read_ctl(control);
  const std::vector<double> cutoffs = Rcpp::as<std::vector<double>>(control["cutoffs"]);
  const int chain_iter = Rcpp::as<int>(control["chain_iter"]), chain_burn = Rcpp::as<int>(control["chain_burn"]);
  const int main_iter = Rcpp::as<int>(control["main_iter"]), main_burn = Rcpp::as<int>(control["main_burn"]);
  const double train_frac = Rcpp::as<double>(control["train_frac"]), rel_r_min = Rcpp::as<double>(control["rel_r_min"]);
  const double bin_r2 = Rcpp::as<double>(control["bin_r2"]);
  const int bin_tile = Rcpp::as<int>(control["bin_tile"]);
  const int nb = blocks.size();
  if (keep.size() != nb) Rcpp::stop("blocks and keep differ in length");
  std::vector<std::vector<int>> kv(nb);
  std::vector<int64_t> off(nb + 1, 0);
  for (int b = 0; b < nb; ++b) { kv[b] = Rcpp::as<std::vector<int>>(keep[b]); off[b + 1] = off[b] + (int64_t)kv[b].size(); }
  const int64_t M = off[nb];
  if (z.size() != M || n.size() != M || af.size() != M || A.nrow() != M)
    Rcpp::stop("z, n, af and A must have one row per kept SNP");
  std::vector<double> bh(M), nn(n.begin(), n.end());
  for (int64_t j = 0; j < M; ++j) bh[j] = z[j] / std::sqrt(n[j] + z[j] * z[j]);
  std::nth_element(nn.begin(), nn.begin() + M / 2, nn.end());
  double N = nn[M / 2];
  if (M % 2 == 0) N = (N + *std::max_element(nn.begin(), nn.begin() + M / 2)) / 2;
  const int q = 1 + A.ncol();
  std::vector<double> F((size_t)M * q);
  for (int64_t j = 0; j < M; ++j) {
    F[(size_t)j * q] = 1.0;
    for (int a = 1; a < q; ++a) F[(size_t)j * q + a] = A(j, a - 1);
  }
  std::unique_ptr<Rng> rng;
  if (replay.isNotNull()) {
    Rcpp::List rp(replay);
    RngReplay* r = new RngReplay();
    r->nrm = Rcpp::as<std::vector<double>>(rp["normal"]);
    r->chi = Rcpp::as<std::vector<double>>(rp["chisq"]);
    for (const std::string& s : Rcpp::as<std::vector<std::string>>(rp["seed"])) r->sd.push_back(std::stoull(s));
    rng.reset(r);
  } else {
    rng.reset(new RngStd((unsigned long long)Rcpp::as<double>(control["pseudo_seed"]),
                         (unsigned long long)Rcpp::as<double>(control["seed"])));
  }
  try {
    // panel
    std::vector<Block16> P(nb);
    std::vector<std::string> err(nb);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
    for (int b = 0; b < nb; ++b) {
      try { P[b] = plink4r::read_block16(plink4r::block_path(ld_dir, blocks[b]), kv[b]); }
      catch (const std::exception& e) { err[b] = e.what(); }
    }
    for (int b = 0; b < nb; ++b) if (!err[b].empty()) Rcpp::stop("block %d: %s", blocks[b], err[b].c_str());
    // LD-bin suppression: annotated = any nonzero binary (0/1) annotation column; bins at r2 > bin_r2
    std::vector<char> sup;
    int64_t n_bins = 0, n_sup = 0;
    std::vector<int> bincol;
    for (int a = 0; a < A.ncol(); ++a) {
      bool bin01 = true;
      for (int64_t j = 0; j < M && bin01; ++j) bin01 = A(j, a) == 0 || A(j, a) == 1;
      if (bin01) bincol.push_back(a);
    }
    if (bin_r2 > 0 && !bincol.empty()) {
      sup.assign(M, 0);
      std::vector<int64_t> nbin(nb, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
      for (int b = 0; b < nb; ++b) {
        const std::vector<int> bid = plink4r::block_bins(plink4r::block_neighbours(P[b], bin_r2, bin_tile));
        std::vector<char> ann(P[b].nk, 0);
        for (int j = 0; j < P[b].nk; ++j)
          for (int a : bincol) if (A(off[b] + j, a) != 0) { ann[j] = 1; break; }
        const std::vector<char> s = plink4r::block_suppress(bid, ann);
        std::copy(s.begin(), s.end(), sup.begin() + off[b]);
        for (int x : bid) nbin[b] = std::max<int64_t>(nbin[b], x + 1);
      }
      for (int b = 0; b < nb; ++b) n_bins += nbin[b];
      for (char x : sup) n_sup += x;
    }
    // pseudo summaries (SBayesRC pseudoSuma): bt = X z per block, z ~ N(0, I_k) on the stored panel
    const double n_t = std::nearbyint(train_frac * N), n_v = N - n_t;
    std::vector<double> bt(M);
    for (int b = 0; b < nb; ++b) {
      const Block16& B = P[b];
      std::vector<float> v(B.ld);
      for (int t = 0; t < B.ld; ++t) v[t] = B.c[t] * (float)rng->normal();
      for (int j = 0; j < B.nk; ++j) {
        const int16_t* qj = B.Q.data() + (size_t)j * B.ld;
        float s = 0.f;
        for (int t = 0; t < B.ld; ++t) s += (float)qj[t] * v[t];
        bt[off[b] + j] = s;
      }
    }
    std::vector<double> bh_t(M), bh_v(M);
    const double sc = std::sqrt(1 / n_t - 1 / N);
    for (int64_t j = 0; j < M; ++j) { bh_t[j] = bh[j] + sc * bt[j]; bh_v[j] = (bh[j] * N - bh_t[j] * n_t) / n_v; }
    const std::vector<std::vector<double>> W_t = wfull(P, off, bh_t, c.lam_eps, T), W = wfull(P, off, bh, c.lam_eps, T);
    // warm-start chain over the cutoffs on the pseudo training data
    std::vector<double> ch_cut, ch_r, ch_r2, ch_sec;
    std::vector<int> ch_k;
    std::vector<FitState> states;
    states.reserve(cutoffs.size());
    const FitState* init = nullptr;
    for (double cut : cutoffs) {
      const auto t0 = std::chrono::steady_clock::now();
      const std::vector<int> ks = kcut(P, cut);
      FitOut r = fit(P, off, ks, W_t, n_t, F, q, chain_iter, chain_burn, init, sup, c, *rng, T);
      double bv = 0, bb = 0;
      for (int64_t j = 0; j < M; ++j) { bv += r.beta[j] * bh_v[j]; bb += r.beta[j] * r.beta[j]; }
      states.push_back(r.st);
      init = &states.back();
      int ktot = 0;
      for (int k : ks) ktot += k;
      ch_cut.push_back(cut); ch_k.push_back(ktot); ch_r.push_back(bv / std::sqrt(prs_var(P, off, r.beta, T)));
      ch_r2.push_back(bv / std::sqrt(bb));
      ch_sec.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
      Rcpp::checkUserInterrupt();
    }
    // SBayesRC rule: base = the largest cutoff; take the best rel_r = |r / r_base| over r > 0 when it exceeds rel_r_min,
    // or when the base r is negative
    const int nc = (int)cutoffs.size();
    const int ib = (int)(std::max_element(cutoffs.begin(), cutoffs.end()) - cutoffs.begin());
    const double r0 = ch_r[ib];
    std::vector<double> rel(nc);
    int best = -1;
    for (int i = 0; i < nc; ++i) {
      rel[i] = std::fabs(ch_r[i] / r0);
      if (ch_r[i] > 0 && (best < 0 || rel[i] > rel[best])) best = i;
    }
    int ich = ib;
    if (best >= 0 && (rel[best] > rel_r_min || (best != ib && r0 < 0))) ich = best;
    Rcpp::DataFrame chain = Rcpp::DataFrame::create(Rcpp::_["cut"] = ch_cut, Rcpp::_["k_total"] = ch_k,
                                                    Rcpp::_["r"] = ch_r, Rcpp::_["r_sbrc_norm"] = ch_r2,
                                                    Rcpp::_["rel_r"] = rel, Rcpp::_["seconds"] = ch_sec);
    const double cmin = *std::min_element(cutoffs.begin(), cutoffs.end());
    if (best < 0) Rcpp::stop("all chain correlations are negative");
    if (cutoffs[ich] == cmin && nc > 1)
      Rcpp::stop("chosen cutoff = smallest step %g (SBayesRC rule): widen the cutoff grid", cmin);
    // main fit on the full data at the chosen cutoff, warm from that level
    const auto t0 = std::chrono::steady_clock::now();
    FitOut r = fit(P, off, kcut(P, cutoffs[ich]), W, N, F, q, main_iter, main_burn, &states[ich], sup, c, *rng, T);
    const double main_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    Rcpp::NumericVector beta(M), bstd(M);
    for (int64_t j = 0; j < M; ++j) {
      const double p = af[j];
      bstd[j] = r.beta[j];
      beta[j] = r.beta[j] / std::sqrt(std::max(2 * p * (1 - p), c.het_floor));
    }
    double nnz = 0;
    for (int i = main_burn; i < (int)r.tr_nnz.size(); ++i) nnz += r.tr_nnz[i];
    nnz /= std::max(1, (int)r.tr_nnz.size() - main_burn);
    Rcpp::NumericMatrix theta(c.gam.size() - 1, q);
    for (int l = 0; l < theta.nrow(); ++l) for (int a = 0; a < q; ++a) theta(l, a) = r.st.theta[(size_t)l * q + a];
    return Rcpp::List::create(
        Rcpp::_["BETA"] = beta, Rcpp::_["BETA_STD"] = bstd, Rcpp::_["chain"] = chain,
        Rcpp::_["cutoff"] = cutoffs[ich], Rcpp::_["theta"] = theta, Rcpp::_["sigma2"] = r.st.s2,
        Rcpp::_["nnz"] = nnz, Rcpp::_["main_seconds"] = main_sec, Rcpp::_["n_bins"] = (double)n_bins,
        Rcpp::_["suppressed"] = Rcpp::LogicalVector(sup.begin(), sup.end()),
        Rcpp::_["trace"] = Rcpp::DataFrame::create(Rcpp::_["vg"] = r.tr_vg, Rcpp::_["nnz"] = r.tr_nnz,
                                                   Rcpp::_["sigma2"] = r.tr_s2));
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}
