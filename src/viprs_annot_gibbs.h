// Final fit of viprs_annot() by a short Gibbs sampler with Rao-Blackwell means (PartitionPRS). Included by
// viprs_annot.cpp inside its anonymous namespace, after the VI parts it reuses (EigenPanel, State, VarModel, mix_post,
// ViprsControl). The VI path chooses pi (grid, warm start) and gives the start (effects, sigma2_0, q(a_k), sigma); here,
// at that pi, with sigma2_0 kept:
//   per sweep and SNP j (blocks in parallel, own random stream per block and sweep): the conditional posterior of
//   beta_j given the others (spike + slab, mix_post), its mean m_j (Rao-Blackwell), then a draw of beta_j and the
//   update of s = V' beta; after the sweep the MASH atoms a_k drawn from q(a_k) of the M-step objective
//   (VarModel::mstep) on the stochastic-approximation average of the RB second moments E[beta_j^2 | others], and sigma = 1 - 2 beta'bt + beta'R beta by stochastic approximation (step 1 during burn-in,
//   1 / (t - burn + 1) after, as sbayesrc()).
// After burn-in, per sweep: m_j, w_jk m_j (w_jk = A_jk sigma2_k / h2_j of that sweep; the columns sum to m_j in every
// sweep, so their means sum to the mean of m_j) and the Shapley value of column k in beta' R beta,
// phi_k = beta^(k)' R beta = sum_j w_jk beta_j (V_j . s) (sum_k phi_k = beta' R beta exactly; the cross terms
// beta^(k)' R beta^(l) are split evenly between k and l).
#ifndef PLINK4R_VIPRS_ANNOT_GIBBS_H
#define PLINK4R_VIPRS_ANNOT_GIBBS_H

struct GibbsOut {
  std::vector<double> bmean, pip;  // RB mean of beta_j (standardized), mean slab probability
  std::vector<double> part;        // M x K row-major: mean of w_jk m_j
  std::vector<double> phi;         // per column: mean of phi_k over the draws
  std::vector<double> pzero;       // MASH: mean over the draws of q(a_k = 0)
  double sigma = 0, h2 = 0;        // mean of sigma and of beta'R beta / (beta'R beta + sigma) over the draws
  int draws = 0;
};

template <typename T>
GibbsOut gibbs_fit(const EigenPanel<T>& P, const std::vector<float>& n, const ViprsControl& c, const State& st0,
                   VarModel& vm, double sigma, int burn, int draws, int threads) {
  const int64_t M = P.M;
  const int K = vm.K, C = (int)c.kap.size();
  const char* sup = vm.sup.data();
  const float lp = (float)(std::log(vm.pi) - std::log1p(-vm.pi));
  std::vector<float> beta(st0.eta), m(M, 0.f), g(M, 0.f), s(P.S, 0.f);
  std::vector<double> e(M, 0.), ebar(M, 0.);  // RB second moment E[beta_j^2 | others] of the sweep, its stochastic approximation
  vtx(P, beta.data(), s.data(), threads);
  GibbsOut o;
  o.bmean.assign(M, 0.); o.pip.assign(M, 0.); o.part.assign((size_t)M * K, 0.);
  o.phi.assign(K, 0.); o.pzero.assign(K, 0.);
  o.draws = draws;
  std::vector<double> psum(K, 0.);  // sigma2_k summed over the draws (reported: the posterior mean)
  std::vector<double> q(M, 0.);
  for (int it = 1; it <= burn + draws; ++it) {
    const float sg = (float)sigma;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
    for (int oi = 0; oi < P.nb; ++oi) {
      const int b = P.order[oi];
      const int k = P.k[b];
      float* sb = s.data() + P.s_off[b];
      std::mt19937_64 eng(c.seed + 0x9E3779B97F4A7C15ull * (unsigned long long)it + 1000003ull * (unsigned long long)b);
      std::uniform_real_distribution<float> U01(0.f, 1.f);
      std::normal_distribution<float> N01;
      for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j) {
        if (sup[j]) continue;
        const T* vj = P.row(b, j);
        const float r = P.bt[j] - (RowStore<T>::dot(vj, sb, k) - P.d[j] * beta[j]);
        const float D = n[j] * r / sg;
        float psi[16], vt[16], m_;
        const float gj = mix_post(D, n[j] * P.d[j] / sg, vm.tau[j], vm.ltau[j], lp, c, m_, psi, vt);
        m[j] = gj * m_; g[j] = gj;
        double e2 = 0;
        for (int i = 0; i < C; ++i) { const double mc = (double)D / vt[i]; e2 += psi[i] * (mc * mc + 1. / vt[i]); }
        e[j] = gj * e2;
        float bn = 0.f;
        if (U01(eng) < gj) {
          float u = U01(eng), acc = 0.f;
          int cls = C - 1;
          for (int i = 0; i < C; ++i) { acc += psi[i]; if (u < acc) { cls = i; break; } }
          bn = D / vt[cls] + N01(eng) / std::sqrt(vt[cls]);
        }
        const float de = bn - beta[j];
        if (de != 0.f) { RowStore<T>::axpy(sb, vj, de, k); beta[j] = bn; }
      }
    }
    // stochastic-approximation M-steps
    const double a = it <= burn ? 1. : 1. / (it - burn + 1);
    double eb = 0, ss = 0;
    for (int64_t j = 0; j < M; ++j) eb += (double)beta[j] * P.bt[j];
    for (int64_t i = 0; i < P.S; ++i) ss += (double)s[i] * s[i];
    if (it > burn) {
      // draws: quantities with the variance model and sigma of this sweep
      const int t = it - burn;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
      for (int b = 0; b < P.nb; ++b)
        for (int64_t j = P.snp_off[b]; j < P.snp_off[b + 1]; ++j)
          q[j] = beta[j] != 0.f ? (double)RowStore<T>::dot(P.row(b, j), s.data() + P.s_off[b], P.k[b]) : 0.;
      std::vector<double> ph(K, 0.);
      for (int64_t j = 0; j < M; ++j) {
        if (sup[j]) continue;
        o.bmean[j] += m[j]; o.pip[j] += g[j];
        double* pj = o.part.data() + (size_t)j * K;
        for (int kk = 0; kk < K; ++kk) {
          const double w = vm.w(j, kk);
          if (w == 0.) continue;
          pj[kk] += w * m[j];
          ph[kk] += w * beta[j] * q[j];
        }
      }
      for (int kk = 0; kk < K; ++kk) o.phi[kk] += ph[kk];
      o.sigma += sigma; o.h2 += ss / (ss + sigma);
      (void)t;
    }
    // the MASH atoms drawn from q(a_k) of the M-step objective on the stochastic-approximation average of the RB
    // second moments (sigma2_0 kept)
    const double ae = it == 1 ? 1. : a;  // the first sweep starts the average
    for (int64_t j = 0; j < M; ++j) ebar[j] = (1. - ae) * ebar[j] + ae * e[j];
    if (!c.fix_sigma2) {
      if (it > burn) for (int kk = 0; kk < K; ++kk) psum[kk] += vm.par[kk];  // sigma2_k of this sweep (as w_jk above)
      std::mt19937_64 meng(c.seed + 0x2545F4914F6CDD1Dull * (unsigned long long)it);
      vm.mstep(ebar.data(), 10, &meng);
      if (it > burn) for (int kk = 1; kk < K; ++kk) o.pzero[kk] += vm.q[(size_t)kk * vm.grid.size()];
    }
    vm.set_all();
    const double ns = 1. - 2. * eb + ss;
    if (ns > 0. && std::isfinite(ns)) sigma = (1. - a) * sigma + a * ns;
  }
  if (draws > 0) {
    const double inv = 1. / draws;
    for (auto& x : o.bmean) x *= inv;
    for (auto& x : o.pip) x *= inv;
    for (auto& x : o.part) x *= inv;
    for (auto& x : o.pzero) x *= inv;
    for (int kk = 0; kk < K; ++kk) o.phi[kk] *= inv;
    o.sigma *= inv; o.h2 *= inv;
    if (!c.fix_sigma2) { for (int kk = 0; kk < K; ++kk) vm.par[kk] = psum[kk] * inv; vm.set_all(); }
  }
  return o;
}

#endif
