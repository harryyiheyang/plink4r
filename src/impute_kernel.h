// Core of the eigen.bin LD imputation of hard-call missing cells, free of Rcpp types so that several
// R entry points share it: genotype_impute_cpp (returns every filled cell) and the LD path of geno_score
// (adds w * filled cell into scores without keeping the cells).
//
// Working rows S: the panel rows of the block (or window) that are also in the PGEN (the intersection; rows
// absent from the PGEN are not carried). Each has a PGEN variant, the ALT frequency p in PGEN orientation and a
// sign s = +1 when PGEN ALT is the panel alt, -1 when REF/ALT are swapped (the panel is not assumed to share the
// PGEN's orientation).
// Stage 1 (per person i):
//   z~ = (g - 2p) / (s sqrt(2p(1-p))) in panel orientation, 0 at missing;  C_S = I - delta (R_SS + delta I)^{-1}
//   (la_kernel.cpp: Woodbury U_S G U_S' or direct, whichever of k and |S| is smaller)
//   z_M = (I - C_MM)^{-1} C_{M,.} z~            ( = R_MO (R_OO + delta I)^{-1} z_O )
//   delta must be > 0 (at delta = 0, C = I - delta (R + delta I)^{-1} degenerates).
// Stage 2 (per SNP j with missing rate >= anchor_min_miss; < 0 = never): shift the missing cells so that
//   sum_obs g + sum_miss g_hat = 2 n p_j, clip [0,2]. Below the cutoff stage 1 is kept: the anchor error is
//   amplified by 2/missing-rate.
//
// Layout for speed (exact, same result as the per-person dot product):
//   with a_lj = C_lj / sd_l (sd signed) and A_j = -sum_l a_lj mu_l,
//   rhs_ij = C_{.j}' z~_i = A_j + sum_{l in M_i} a_lj mu_l + sum_{l het_i} a_lj + 2 sum_{l homalt_i} a_lj.
//   People are processed in blocks; inside a block the loop runs over SNP columns j (column a_.j / C_.j hot in
//   cache) and over the people missing j, using 1-bit planes (het / homalt; la_kernel.cpp plane_dot), so C is
//   streamed from memory once per person block instead of once per (person, missing SNP). C_MM is gathered
//   from the same hot column.
//
// Pipeline: build_impute_ld -> read_block_genotypes -> build_planes -> impute_stage1 -> impute_stage2,
// then for_each_cell visits every filled cell (i, j, g1, g). Errors are thrown as std::runtime_error.
#ifndef PLINK4R_IMPUTE_KERNEL_H
#define PLINK4R_IMPUTE_KERNEL_H

#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"
#include "pgen_reader.h"
#include "la_dispatch.h"
#include "resources.h"

namespace plink4r {
namespace ik {

using Eigen::MatrixXf;
using Eigen::VectorXf;

inline double secs(std::chrono::steady_clock::time_point a) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
}
inline int tid() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

// Size of the current team (read-only; may be below num_threads(T) when the runtime limits threads).
inline int team() {
#ifdef _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}
}  // namespace ik

// ---- LD of one block: C, the scaled columns a_lj and A_j, the anchor mean / sd per row.
struct ImputeLD {
  int m = 0, k = 0, mp = 0;
  std::vector<int> vi;                  // 0-based PGEN variant per working row
  std::vector<float> mu, sd, isd, Aj;   // 2p, s sqrt(2p(1-p)) (signed), 1/sd, A_j
  Eigen::MatrixXf C, Aw;
  double t_eig = 0, t_C = 0;
};

// ld: the eigen block; rows: 0-based eigen rows of the working set S (panel rows present in the PGEN,
// increasing); vi, af, sign: per working row (PGEN variant, PGEN-ALT frequency in (0,1), +1 / -1).
inline ImputeLD build_impute_ld(const EigenBlock& ld, const std::vector<int>& rows, const std::vector<int>& vi,
                                const std::vector<double>& af, const std::vector<int>& sign, double delta,
                                int T, double t_eig = 0) {
  using Eigen::MatrixXf;
  using Eigen::VectorXf;
  auto t0 = std::chrono::steady_clock::now();
  ImputeLD L;
  L.t_eig = t_eig;
  const int m = static_cast<int>(rows.size()), k = ld.k;
  if (m == 0) throw std::runtime_error("no working rows (no panel SNP of the block is in the PGEN)");
  if (!(delta > 0)) throw std::runtime_error("delta must be > 0");
  L.m = m; L.k = k;
  if (static_cast<int>(vi.size()) != m || static_cast<int>(af.size()) != m || static_cast<int>(sign.size()) != m)
    throw std::runtime_error("rows / var_idx / af / sign differ in length");
  L.vi = vi;
  L.mu.assign(m, 0.0f); L.sd.assign(m, 1.0f); L.isd.assign(m, 0.0f);
  for (int j = 0; j < m; ++j) {
    if (vi[j] < 0) throw std::runtime_error("negative var_idx at working row " + std::to_string(j + 1));
    if (sign[j] != 1 && sign[j] != -1) throw std::runtime_error("sign must be +1 or -1");
    if (j > 0 && rows[j] <= rows[j - 1]) throw std::runtime_error("rows must be increasing");
    const double p = af[j];
    if (!(p > 0 && p < 1))
      throw std::runtime_error("af must be in (0,1) for every working row (row " + std::to_string(j + 1) + ")");
    L.mu[j] = static_cast<float>(2 * p);
    L.sd[j] = static_cast<float>(sign[j] * std::sqrt(2 * p * (1 - p)));
    L.isd[j] = 1.0f / L.sd[j];
  }
  L.t_eig += ik::secs(t0); t0 = std::chrono::steady_clock::now();
  const int mp = (m + 15) / 16 * 16;  // padded rows for the 16-wide plane dot
  L.mp = mp;
  MatrixXf& C = L.C;
  C.resize(m, m);
  for (int r = 0; r < static_cast<int>(rows.size()); ++r)
    if (rows[r] < 0 || rows[r] >= ld.m) throw std::runtime_error("rows out of range");
  {
    char err[256] = {0};
    if (la().build_C(ld.U.data(), ld.m, k, ld.lambda.data(), delta, rows.data(), static_cast<int>(rows.size()),
                     C.data(), T, err) != 0)
      throw std::runtime_error(err);
  }
  MatrixXf& Aw = L.Aw;
  Aw = MatrixXf::Zero(mp, m);  // a_lj, padded rows are 0
  L.Aj.assign(m, 0.0f);
  const std::vector<float>& isd = L.isd;
  const std::vector<float>& mu = L.mu;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(T)
#endif
  for (int j = 0; j < m; ++j) {
    float s = 0;
    for (int l = 0; l < m; ++l) { const float v = C(l, j) * isd[l]; Aw(l, j) = v; }
    for (int l = 0; l < m; ++l) s += Aw(l, j) * mu[l];
    L.Aj[j] = -s;
  }
  L.t_C = ik::secs(t0);
  return L;
}

// ---- genotypes of one block: SNP-major genovecs of the selected samples (codes 0/1/2 ALT count, 3 missing)
struct BlockGenotypes {
  int m = 0, ns = 0, sb = 0;
  std::vector<uint8_t> GS;       // m x sb bytes
  std::vector<double> sum_obs;   // observed ALT count per row
  double t_geno = 0;
};

inline BlockGenotypes read_block_genotypes(PgenFile& pf, const std::vector<int>& vi, const std::vector<int>& sidx, int T) {
  auto t0 = std::chrono::steady_clock::now();
  BlockGenotypes B;
  const int m = static_cast<int>(vi.size()), ns = static_cast<int>(sidx.size()), sb = (ns + 3) / 4;
  B.m = m; B.ns = ns; B.sb = sb;
  B.GS.assign(static_cast<std::size_t>(m) * sb, 0);
  B.sum_obs.assign(m, 0.0);
  std::vector<int> bad(T, -1);
  const bool ident = identity_samples(sidx.data(), ns, static_cast<int>(pf.samples()));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
  for (int j = 0; j < m; ++j) {
    const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[j]), ik::tid());
    if (!gv) { bad[ik::tid()] = j; continue; }
    uint8_t* dst = &B.GS[static_cast<std::size_t>(j) * sb];
    compact_genovec(gv, sidx.data(), ns, ident, dst);
    long long miss = 0, alt = 0;
    count_genovec(dst, sb, &miss, &alt);
    B.sum_obs[j] = static_cast<double>(alt);
  }
  for (int t = 0; t < T; ++t)
    if (bad[t] >= 0) throw std::runtime_error("PGEN read failed at eigen row " + std::to_string(bad[t] + 1));
  B.t_geno = ik::secs(t0);
  return B;
}

// ---- bit-planes and missing lists of every person (chunks of 256 people = one cache line per SNP row)
struct PersonPlanes {
  int ns = 0, nbyte = 0;
  std::vector<uint8_t> het, hom;          // ns x nbyte
  std::vector<std::vector<int>> Mi;       // missing rows per person, ascending
  double t_planes = 0;
};

// Consumes B.GS (released on return).
inline PersonPlanes build_planes(const ImputeLD& L, BlockGenotypes& B, int T) {
  auto t0 = std::chrono::steady_clock::now();
  PersonPlanes P;
  const int m = L.m, ns = B.ns, sb = B.sb;
  const int nbyte = L.mp / 8;
  P.ns = ns; P.nbyte = nbyte;
  P.het.assign(static_cast<std::size_t>(ns) * nbyte, 0);
  P.hom.assign(static_cast<std::size_t>(ns) * nbyte, 0);
  P.Mi.assign(ns, std::vector<int>());
  const std::vector<uint8_t>& GS = B.GS;
  const int nch = (ns + 255) / 256;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
  for (int ch = 0; ch < nch; ++ch) {
    const int i0 = ch * 256, i1 = std::min(ns, i0 + 256);
    for (int lb = 0; lb < nbyte; ++lb) {  // 8 SNPs -> one plane byte per person, branch-free
      const int l0 = 8 * lb, nl = std::min(8, m - l0);
      if (nl <= 0) break;
      const uint8_t* g[8];
      for (int b = 0; b < nl; ++b) g[b] = &GS[static_cast<std::size_t>(l0 + b) * sb];
      for (int s = i0; s < i1; ++s) {
        const int byte = s >> 2, sh = 2 * (s & 3);
        unsigned h = 0, o = 0, x = 0;
        for (int b = 0; b < nl; ++b) {
          const unsigned c = (g[b][byte] >> sh) & 3u;
          h |= static_cast<unsigned>(c == 1u) << b;
          o |= static_cast<unsigned>(c == 2u) << b;
          x |= static_cast<unsigned>(c == 3u) << b;
        }
        P.het[static_cast<std::size_t>(s) * nbyte + lb] = static_cast<uint8_t>(h);
        P.hom[static_cast<std::size_t>(s) * nbyte + lb] = static_cast<uint8_t>(o);
        while (x) { P.Mi[s].push_back(l0 + __builtin_ctz(x)); x &= x - 1; }
      }
    }
  }
  std::vector<uint8_t>().swap(B.GS);
  P.t_planes = ik::secs(t0);
  return P;
}

// ---- stage 1: the conditional mean of every missing cell, per thread lists (person, row, g1 clipped to [0,2])
struct Stage1Result {
  std::vector<std::vector<int>> oi, oj;
  std::vector<std::vector<float>> og;
  long long persons_solved = 0;
  int n_blocks = 0;
  double t_lists = 0, t_rhs = 0, t_chol = 0, t_solve = 0;
};

inline Stage1Result impute_stage1(const ImputeLD& L, const PersonPlanes& P, double delta, int nb_people,
                                  double mem_gb, int T) {
  using Eigen::MatrixXf;
  using Eigen::VectorXf;
  const auto t_start = std::chrono::steady_clock::now();
  Stage1Result R;
  const int m = L.m, ns = P.ns, nbyte = P.nbyte;
  const MatrixXf& C = L.C;
  const MatrixXf& Aw = L.Aw;
  const std::vector<float>& mu = L.mu;
  const std::vector<float>& sd = L.sd;
  const std::vector<std::vector<int>>& MiA = P.Mi;
  const auto plane_dot = la().plane_dot;
  // person blocks: fixed size, or as large as the memory budget for S_i, rhs and the SNP lists allows
  std::vector<int> bstart{0};
  {
    const double budget = mem_gb * 1e9;
    double acc = 0;
    for (int s = 0; s < ns; ++s) {
      const double q = static_cast<double>(MiA[s].size());
      const double cost = 4.0 * q * q + 12.0 * q;
      const bool full = nb_people > 0 ? (s - bstart.back() >= nb_people) : (acc + cost > budget && s > bstart.back());
      if (full) { bstart.push_back(s); acc = 0; }
      acc += cost;
    }
    bstart.push_back(ns);
  }
  const int nblk = static_cast<int>(bstart.size()) - 1;
  R.n_blocks = nblk;
  R.oi.assign(T, std::vector<int>()); R.oj.assign(T, std::vector<int>()); R.og.assign(T, std::vector<float>());
  std::vector<int> moff, jcnt(m + 1), jlist_i, jlist_a;
  std::vector<long long> soff;
  std::vector<float> rhs, S;
  for (int bk = 0; bk < nblk; ++bk) {
    const int b0 = bstart[bk], nb = bstart[bk + 1] - b0;
    auto tb = std::chrono::steady_clock::now();
    const uint8_t* het = &P.het[static_cast<std::size_t>(b0) * nbyte];
    const uint8_t* hom = &P.hom[static_cast<std::size_t>(b0) * nbyte];
    const std::vector<int>* Mi = &MiA[b0];
    // offsets and per-SNP lists of (person, position in M)
    moff.assign(nb + 1, 0);
    soff.assign(nb + 1, 0);
    for (int i = 0; i < nb; ++i) {
      const int q = static_cast<int>(Mi[i].size());
      moff[i + 1] = moff[i] + q;
      soff[i + 1] = soff[i] + static_cast<long long>(q) * q;
    }
    // per-SNP lists in parallel: thread t counts its person range, then fills its own slots (order = person order)
    {
      std::vector<int> tcnt(static_cast<std::size_t>(T) * m, 0);
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
      {
        const int th = ik::tid(), nt = ik::team();  // split by the team actually granted, not by T
        const int i0 = static_cast<int>(static_cast<long long>(nb) * th / nt), i1 = static_cast<int>(static_cast<long long>(nb) * (th + 1) / nt);
        int* c = &tcnt[static_cast<std::size_t>(th) * m];
        for (int i = i0; i < i1; ++i) for (int l : Mi[i]) ++c[l];
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
        {
          jcnt[0] = 0;
          for (int j = 0; j < m; ++j) {
            int run = jcnt[j];
            for (int t = 0; t < T; ++t) { const int v = tcnt[static_cast<std::size_t>(t) * m + j]; tcnt[static_cast<std::size_t>(t) * m + j] = run; run += v; }
            jcnt[j + 1] = run;
          }
          jlist_i.resize(moff[nb]); jlist_a.resize(moff[nb]);
        }
        for (int i = i0; i < i1; ++i)
          for (int a = 0; a < static_cast<int>(Mi[i].size()); ++a) {
            const int p = c[Mi[i][a]]++;
            jlist_i[p] = i; jlist_a[p] = a;
          }
      }
    }
    rhs.resize(moff[nb]);                          // fully overwritten below
    S.resize(static_cast<std::size_t>(soff[nb]));  // every column of every S_i is written below
    R.t_lists += ik::secs(tb); tb = std::chrono::steady_clock::now();

    // rhs and C_MM, column by column
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(T)
#endif
    for (int j = 0; j < m; ++j) {
      const int e0 = jcnt[j], e1 = jcnt[j + 1];
      if (e0 == e1) continue;
      const float* aj = Aw.col(j).data();
      const float* cj = C.col(j).data();
      for (int e = e0; e < e1; ++e) {
        const int i = jlist_i[e], a = jlist_a[e];
        const std::vector<int>& M = Mi[i];
        const int q = static_cast<int>(M.size());
        const uint8_t* hi = &het[static_cast<std::size_t>(i) * nbyte];
        const uint8_t* oi_ = &hom[static_cast<std::size_t>(i) * nbyte];
        float r = L.Aj[j] + plane_dot(aj, hi, oi_, nbyte);
        float* Sc = &S[static_cast<std::size_t>(soff[i]) + static_cast<std::size_t>(a) * q];  // column a of S_i
        for (int bq = 0; bq < q; ++bq) {
          const int l = M[bq];
          r += aj[l] * mu[l];
          Sc[bq] = -cj[l];
        }
        Sc[a] += 1.0f;
        rhs[moff[i] + a] = r;
      }
    }
    R.t_rhs += ik::secs(tb); tb = std::chrono::steady_clock::now();

    // per-person solve
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
    {
      const int th = ik::tid();
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
      for (int i = 0; i < nb; ++i) {
        const std::vector<int>& M = Mi[i];
        const int q = static_cast<int>(M.size());
        if (q == 0) continue;
        Eigen::Map<MatrixXf> Si(&S[static_cast<std::size_t>(soff[i])], q, q);
        Eigen::Map<VectorXf> ri(&rhs[moff[i]], q);
        VectorXf zm;
        Eigen::LLT<Eigen::Ref<MatrixXf>> llt(Si);
        if (llt.info() == Eigen::Success) zm = llt.solve(ri);
        else {
          MatrixXf Sf(q, q);  // LLT overwrote the lower triangle and diagonal; rebuild from the strict upper one
          for (int x = 0; x < q; ++x)
            for (int y = 0; y < q; ++y) Sf(y, x) = y < x ? Si(y, x) : (y > x ? Si(x, y) : 1.0f - C(M[x], M[x]));
          zm = Sf.ldlt().solve(ri);
        }
        for (int x = 0; x < q; ++x) {
          const int j = M[x];
          R.oi[th].push_back(b0 + i); R.oj[th].push_back(j);
          R.og[th].push_back(std::min(2.0f, std::max(0.0f, mu[j] + sd[j] * zm[x])));
        }
      }
    }
    R.t_chol += ik::secs(tb);
    for (int i = 0; i < nb; ++i) if (!Mi[i].empty()) ++R.persons_solved;
  }
  R.t_solve = ik::secs(t_start);
  return R;
}

// ---- stage 2: anchor per SNP. NaN marks rows without filled cells (nothing missing).
struct Stage2Result {
  std::vector<double> sum1, cnt1, craw, cc, shift, mean1;
  double t_anchor = 0;
};

inline Stage2Result impute_stage2(const Stage1Result& R, const std::vector<double>& sum_obs,
                                  const std::vector<double>& af, int ns, double anchor_min_miss) {
  const int m = static_cast<int>(sum_obs.size());
  const double nan = std::numeric_limits<double>::quiet_NaN();
  Stage2Result A;
  A.sum1.assign(m, 0.0); A.cnt1.assign(m, 0.0);
  const int T = static_cast<int>(R.og.size());
  for (int t = 0; t < T; ++t)
    for (std::size_t r = 0; r < R.og[t].size(); ++r) { A.sum1[R.oj[t][r]] += R.og[t][r]; A.cnt1[R.oj[t][r]] += 1; }
  A.craw.assign(m, nan); A.cc.assign(m, nan); A.shift.assign(m, 0.0); A.mean1.assign(m, nan);
  for (int j = 0; j < m; ++j) {
    if (A.cnt1[j] == 0) continue;
    A.mean1[j] = A.sum1[j] / A.cnt1[j];
    A.craw[j] = (2.0 * ns * af[j] - sum_obs[j]) / A.cnt1[j];
    A.cc[j] = std::min(2.0, std::max(0.0, A.craw[j]));
    if (anchor_min_miss >= 0 && A.cnt1[j] >= anchor_min_miss * ns) A.shift[j] = A.cc[j] - A.mean1[j];
  }
  return A;
}

// Final fill of a cell from its stage-1 value and the row's shift.
inline double final_fill(float g1, double shift) { return std::min(2.0, std::max(0.0, g1 + shift)); }

// Visits every filled cell of the kept rows: f(i = person, j = row, g1 = stage 1, g = final), thread-list order.
template <class F>
inline void for_each_cell(const Stage1Result& R, const Stage2Result& A, const std::vector<char>& keep, F f) {
  const int T = static_cast<int>(R.og.size());
  for (int t = 0; t < T; ++t)
    for (std::size_t r = 0; r < R.oj[t].size(); ++r) {
      const int j = R.oj[t][r];
      if (!keep[j]) continue;
      f(R.oi[t][r], j, R.og[t][r], final_fill(R.og[t][r], A.shift[j]));
    }
}

}  // namespace plink4r
#endif
