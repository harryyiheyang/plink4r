// Annotation-partitioned PRS by a short Gibbs sampler with Rao-Blackwellization (PartitionPRS, Research Design
// D.1.2). Top layer: the package's partitioned-PRS functional unit; copied from viprs_annot.cpp and edited (the
// variational fit of viprs_annot() chooses pi and gives the start; this file only samples).
// Model as viprs_annot(): h2_j = sigma2_0 + sum_k A_jk sigma2_k (sigma2_0 the intercept on every SNP), beta_j = 0 with
// probability 1 - pi, else class c (prior share p_c, fixed) with beta_j ~ N(0, kap_c h2_j / pi), sum_c p_c kap_c = 1.
// One sweep (blocks in parallel: the LD is block diagonal; SNPs of a block in order): D_j = n_j r_j / sigma,
// vt_jc = n_j d_j / sigma + tau_j / kap_c (tau_j = pi / h2_j); the exact conditional of beta_j given the others is the
// spike with probability 1 - gamma_j and class c with gamma_j psi_jc, beta_j | c ~ N(D_j / vt_jc, 1 / vt_jc)
// (gamma, psi as in viprs_annot). Rao-Blackwellized mean of a draw: gamma_j sum_c psi_jc D_j / vt_jc.
// sigma2_k and sigma by SAEM (step 1 during burn-in, 1 / (t - burn + 1) after): the M-step of D.1.2 (moment form,
// sigma2_k <- (S_k + 2 b0) / (sum_j A_jk + 2 (a0 + 1)), S_k = sum_j w_jk^2 zeta_j + w_jk (1 - w_jk) h2_j) on the
// stochastic average of the Rao-Blackwellized zeta_j, and sigma on the average of 1 - 2 beta'b + beta'R beta.
// Output: the Rao-Blackwellized posterior mean b_j (average over the draws after burn-in), its MMSE allocation
// b_j^(k) = w_jk b_j with w_jk = A_jk sigma2_k / h2_j at the final sigma2 (sum_k b_j^(k) = b_j: the partitioned PRS add
// up), the model heritability h2_k = (sum_j A_jk) sigma2_k (additive) and the Shapley share of the PRS variance
// phi_k = b^(k)' R b (the Shapley value of v(S) = b_S' R b_S; sum_k phi_k = b' R b).
// Suppressed SNPs (LD-bin skip, from viprs_annot()) stay at beta = 0 and are left out of the M-step.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <random>
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

struct GibbsControl {
  int iter, burn;
  double a0, b0, var_floor, sigma_fixed, af_clip;
  unsigned long long seed;
  std::vector<float> lprop, kap, lkap;
};

GibbsControl read_control(const Rcpp::List& c) {
  GibbsControl x;
  x.iter = Rcpp::as<int>(c["iter"]);
  x.burn = Rcpp::as<int>(c["burn"]);
  x.a0 = Rcpp::as<double>(c["a0"]);
  x.b0 = Rcpp::as<double>(c["b0"]);
  x.var_floor = Rcpp::as<double>(c["var_floor"]);
  x.sigma_fixed = Rcpp::as<double>(c["sigma_fixed"]);
  x.af_clip = Rcpp::as<double>(c["af_clip"]);
  x.seed = static_cast<unsigned long long>(Rcpp::as<double>(c["seed"]));
  if (x.iter <= x.burn || x.burn < 0) Rcpp::stop("iter must exceed burn >= 0");
  const std::vector<double> mp = Rcpp::as<std::vector<double>>(c["mix_prop"]), mv = Rcpp::as<std::vector<double>>(c["mix_var"]);
  if (mp.empty() || mp.size() != mv.size()) Rcpp::stop("mix_prop and mix_var need the same length >= 1");
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

// variance model h2_j = sigma2_0 + sum_k A_jk sigma2_k (A row-major M x (K - 1)), M-step of D.1.2
struct VarModel {
  int64_t M;
  int K;
  const std::vector<double>& A;
  const std::vector<char>& sup;
  std::vector<double> par, h2;
  std::vector<float> tau, ltau;
  double pi = 1, vfloor = 0;
  int threads;
  VarModel(int64_t M_, int K_, const std::vector<double>& A_, const std::vector<char>& sup_, int threads_)
      : M(M_), K(K_), A(A_), sup(sup_), par(K_, 0.), h2(M_, 1.), tau(M_, 1.f), ltau(M_, 0.f), threads(threads_) {}
  const double* row(int64_t j) const { return A.data() + j * (K - 1); }
  void set_all() {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (int64_t j = 0; j < M; ++j) {
      double h = par[0];
      const double* r = row(j);
      for (int k = 1; k < K; ++k) h += r[k - 1] * par[k];
      h = std::max(h, vfloor);
      h2[j] = h;
      const double t = pi / h;
      tau[j] = (float)t; ltau[j] = (float)std::log(t);
    }
  }
  void mstep(const double* zeta, double a0, double b0) {
    std::vector<double> SN(2 * K, 0.);
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
          const double a = k ? r[k - 1] : 1.;
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
      par[k] = den > 0 ? std::max((SN[k] + 2. * b0) / den, vfloor) : vfloor;
    }
    set_all();
  }
};

// One Gibbs sweep. beta: current draw (updated, s = V' beta kept in step); rb: Rao-Blackwellized mean of each SNP at
// its update; gam: slab probability; zeta: Rao-Blackwellized E[beta_j^2] on the scale of h2_j (the M-step statistic).
template <typename T>
void gibbs_sweep(const EigenPanel<T>& P, const float* n, const float* bt, const VarModel& vm, const GibbsControl& c,
                 float sigma, unsigned long long seed, std::vector<float>& beta, std::vector<float>& s,
                 std::vector<float>& rb, std::vector<float>& gam, std::vector<double>& zeta, int threads) {
  const int C = (int)c.kap.size();
  const float lp = (float)(std::log(vm.pi) - std::log1p(-vm.pi));
  const float* d = P.d.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int oi = 0; oi < P.nb; ++oi) {
    const int b = P.order[oi];
    const int k = P.k[b];
    float* sb = s.data() + P.s_off[b];
    std::mt19937_64 eng(seed * 1000003ull + (unsigned long long)b);
    std::uniform_real_distribution<float> U01(0.f, 1.f);
    std::normal_distribution<float> N01;
    std::vector<float> u(C), v(C);
    for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j) {
      if (vm.sup[j]) continue;
      const T* vj = P.row(b, j);
      const float r = bt[j] - (RowStore<T>::dot(vj, sb, k) - d[j] * beta[j]);
      const float D = n[j] * r / sigma, base = n[j] * d[j] / sigma;
      float umax = -INFINITY;
      for (int i = 0; i < C; ++i) {
        v[i] = base + vm.tau[j] / c.kap[i];
        u[i] = c.lprop[i] + .5f * D * D / v[i] + .5f * (vm.ltau[j] - c.lkap[i] - std::log(v[i]));
        umax = std::max(umax, u[i]);
      }
      float se = 0.f;
      for (int i = 0; i < C; ++i) { u[i] = std::exp(u[i] - umax); se += u[i]; }
      const float x = lp + umax + std::log(se);
      const float g = x < 0 ? std::exp(x) / (1.f + std::exp(x)) : 1.f / (1.f + std::exp(-x));
      float m = 0.f;
      double z = 0;
      for (int i = 0; i < C; ++i) {
        const float psi = u[i] / se, mc = D / v[i];
        m += psi * mc;
        z += psi * ((double)mc * mc + 1. / v[i]) / c.kap[i];
      }
      rb[j] = g * m;
      gam[j] = g;
      zeta[j] = g * z;
      float nb_ = 0.f;
      if (U01(eng) < g) {
        float q = U01(eng) * se, acc = 0.f;
        int ci = C - 1;
        for (int i = 0; i < C; ++i) { acc += u[i]; if (q <= acc) { ci = i; break; } }
        nb_ = D / v[ci] + N01(eng) / std::sqrt(v[ci]);
      }
      const float de = nb_ - beta[j];
      if (de != 0.f) { RowStore<T>::axpy(sb, vj, de, k); beta[j] = nb_; }
    }
  }
}

double dot_f(const std::vector<float>& a, const std::vector<float>& b) {
  double x = 0;
  for (std::size_t i = 0; i < a.size(); ++i) x += (double)a[i] * b[i];
  return x;
}

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
    in.bstd[j] = (float)(z[j] / std::sqrt(n[j] - 1. + z[j] * z[j]));
  }
  return in;
}

template <typename T>
Rcpp::List fit(const Inputs& in, const std::vector<double>& A, int K, const std::vector<char>& sup, double pi,
               const std::vector<double>& s2, double sigma0, const std::vector<float>& beta0, double h2_ref,
               const std::vector<std::string>& files, double thresh, const GibbsControl& c, int threads) {
  const EigenPanel<T> P = plink4r::load_eigen_panel<T>(files, in.blocks, in.keep, in.bstd, thresh, threads);
  const int64_t M = P.M;
  VarModel vm(M, K, A, sup, threads);
  vm.pi = pi;
  vm.vfloor = c.var_floor * h2_ref / (double)M;
  for (int k = 0; k < K; ++k) vm.par[k] = std::max(s2[k], vm.vfloor);
  vm.set_all();
  std::vector<float> beta(beta0), s(P.S), rb(M, 0.f), gam(M, 0.f);
  for (int64_t j = 0; j < M; ++j) if (sup[j]) beta[j] = 0.f;
  vtx(P, beta.data(), s.data(), threads);
  std::vector<double> zeta(M, 0.), zbar(M, 0.), rbsum(M, 0.), gsum(M, 0.);
  double sigma = sigma0, sbar = sigma0, gv_sum = 0;
  const int keep_n = c.iter - c.burn;
  Rcpp::NumericMatrix trace(c.iter, K + 2);  // sigma2_k, sigma, beta' R beta of the draw
  const auto t0 = std::chrono::steady_clock::now();
  for (int t = 1; t <= c.iter; ++t) {
    gibbs_sweep(P, in.n.data(), P.bt.data(), vm, c, (float)sigma, c.seed + (unsigned long long)t, beta, s, rb, gam, zeta,
                threads);
    vtx(P, beta.data(), s.data(), threads);  // refresh V' beta (float drift)
    const double bRb = dot_f(s, s), bb = [&] { double x = 0; for (int64_t j = 0; j < M; ++j) x += (double)beta[j] * P.bt[j]; return x; }();
    const double a = t <= c.burn ? 1. : 1. / (double)(t - c.burn + 1);
    for (int64_t j = 0; j < M; ++j) zbar[j] = (1. - a) * zbar[j] + a * zeta[j];
    const double ns = 1. - 2. * bb + bRb;
    sbar = (1. - a) * sbar + a * (ns > 0 && std::isfinite(ns) ? ns : c.sigma_fixed);
    sigma = sbar;
    vm.mstep(zbar.data(), c.a0, c.b0);
    if (t > c.burn) {
      for (int64_t j = 0; j < M; ++j) { rbsum[j] += rb[j]; gsum[j] += gam[j]; }
      gv_sum += bRb;
    }
    for (int k = 0; k < K; ++k) trace(t - 1, k) = vm.par[k];
    trace(t - 1, K) = sigma; trace(t - 1, K + 1) = bRb;
    Rcpp::checkUserInterrupt();
  }
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  // outputs at the final sigma2
  Rcpp::NumericVector bstd(M), beta_out(M), pip(M), h2j(M), h2a(K), phi(K);
  Rcpp::NumericMatrix part(M, K);
  std::vector<float> bm(M), x(M), sk(P.S), sb(P.S);
  for (int64_t j = 0; j < M; ++j) {
    bm[j] = (float)(rbsum[j] / keep_n);
    const double p = std::min(std::max(in.af[j], c.af_clip), 1 - c.af_clip);
    bstd[j] = bm[j];
    beta_out[j] = bm[j] / std::sqrt(2. * p * (1. - p));
    pip[j] = gsum[j] / keep_n;
    h2j[j] = vm.h2[j];
    const double* r = vm.row(j);
    for (int k = 0; k < K; ++k) {
      const double a = k ? r[k - 1] : 1.;
      part(j, k) = a * vm.par[k] / vm.h2[j] * bm[j];
      if (!sup[j]) h2a[k] += a * vm.par[k];
    }
  }
  vtx(P, bm.data(), sb.data(), threads);
  const double bRb_mean = dot_f(sb, sb);
  for (int k = 0; k < K; ++k) {
    for (int64_t j = 0; j < M; ++j) x[j] = (float)part(j, k);
    vtx(P, x.data(), sk.data(), threads);
    phi[k] = dot_f(sk, sb);
  }
  const double gv = gv_sum / keep_n;
  Rcpp::List hyp = Rcpp::List::create(
      Rcpp::_["pi"] = pi, Rcpp::_["sigma"] = sigma, Rcpp::_["h2"] = gv / (gv + sigma), Rcpp::_["bRb"] = bRb_mean,
      Rcpp::_["iter"] = c.iter, Rcpp::_["burn"] = c.burn, Rcpp::_["n_snps"] = (double)M, Rcpp::_["sum_k"] = (double)P.S,
      Rcpp::_["seconds"] = secs);
  return Rcpp::List::create(Rcpp::_["BETA"] = beta_out, Rcpp::_["BETA_STD"] = bstd, Rcpp::_["PIP"] = pip,
                            Rcpp::_["H2_SNP"] = h2j, Rcpp::_["PART_STD"] = part,
                            Rcpp::_["sigma2"] = Rcpp::NumericVector(vm.par.begin(), vm.par.end()),
                            Rcpp::_["h2_annot"] = h2a, Rcpp::_["h2_shapley"] = phi, Rcpp::_["hyp"] = hyp,
                            Rcpp::_["trace"] = trace);
}

}  // namespace

// Entry: files (eigen.bin of each block), blocks, keep (0-based eigen rows with data per block), z / n / af per kept
// SNP, A (kept SNPs x annotations; the intercept is added here), suppressed (LD-bin skip of viprs_annot()), pi (chosen
// by viprs_annot()), sigma2 (start, intercept first), sigma (start), beta (start, standardized), h2_ref (scale of
// var_floor: var_floor * h2_ref / M), thresh, fp16, control (partprs_control()).
// [[Rcpp::export]]
Rcpp::List partprs_cpp(const std::vector<std::string>& files, Rcpp::IntegerVector blocks, Rcpp::List keep,
                       Rcpp::NumericVector z, Rcpp::NumericVector n, Rcpp::NumericVector af, Rcpp::NumericMatrix A,
                       Rcpp::LogicalVector suppressed, double pi, Rcpp::NumericVector sigma2, double sigma,
                       Rcpp::NumericVector beta, double h2_ref, double thresh, bool fp16, Rcpp::List control,
                       int threads) {
  const int T = plink4r::resolve_threads(threads);
  const GibbsControl c = read_control(control);
  const Inputs in = read_inputs(blocks, keep, z, n, af);
  const int64_t M = in.z.size();
  if (A.nrow() != M || suppressed.size() != M || beta.size() != M) Rcpp::stop("A, suppressed and beta need one row per kept SNP");
  const int K = A.ncol() + 1;
  if (sigma2.size() != K) Rcpp::stop("sigma2 needs one value per annotation plus the intercept");
  if (!(pi > 0 && pi < 1) || !(sigma > 0)) Rcpp::stop("pi must be in (0, 1) and sigma > 0");
  std::vector<double> Ar(M * (K - 1));
  for (int64_t j = 0; j < M; ++j)
    for (int k = 0; k < K - 1; ++k) {
      const double x = A(j, k);
      if (!std::isfinite(x) || x < 0) Rcpp::stop("annotations must be finite and >= 0");
      Ar[j * (K - 1) + k] = x;
    }
  std::vector<char> sup(M);
  for (int64_t j = 0; j < M; ++j) sup[j] = suppressed[j] == TRUE;
  const std::vector<double> s2(sigma2.begin(), sigma2.end());
  std::vector<float> b0(M);
  for (int64_t j = 0; j < M; ++j) b0[j] = (float)beta[j];
  try {
    return fp16 ? fit<uint16_t>(in, Ar, K, sup, pi, s2, sigma, b0, h2_ref, files, thresh, c, T)
                : fit<float>(in, Ar, K, sup, pi, s2, sigma, b0, h2_ref, files, thresh, c, T);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
}
