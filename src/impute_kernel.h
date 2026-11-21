// Core of the eigen.bin LD imputation of missing cells (hard calls and dosages), free of Rcpp types so that two
// R entry points share it: the LD fill of genotype_fill (genotype_fill_ld.cpp: writes the filled cells into a geno
// object) and the LD path of geno_score (adds w * filled cell into scores without keeping the cells).
//
// Working rows S: the panel rows of the block (or window) that are also in the PGEN (the intersection; rows
// absent from the PGEN are not carried). Each has a PGEN variant, the ALT frequency p in PGEN orientation and a
// sign s = +1 when PGEN ALT is the panel alt, -1 when REF/ALT are swapped (the panel is not assumed to share the
// PGEN's orientation).
// Stage 1 (per person i):
//   z~ = (g - 2p) / (s sqrt(2p(1-p))) in panel orientation, 0 at missing (g = observed hard call or dosage);  C_S = I - delta (R_SS + delta I)^{-1}
//   (la_kernel.cpp: Woodbury U_S G U_S' or direct, whichever of k and |S| is smaller)
//   z_M = (I - C_MM)^{-1} C_{M,.} z~            ( = R_MO (R_OO + delta I)^{-1} z_O )
//   delta must be > 0 (at delta = 0, C = I - delta (R + delta I)^{-1} degenerates).
// Stage 2 (per SNP j with missing rate >= anchor_min_miss; < 0 = never): shift the missing cells so that
//   sum_obs g + sum_miss g_hat = 2 n p_j, clip [0,2]. Below the cutoff stage 1 is kept: the anchor error is
//   amplified by 2/missing-rate.
//
// Layout for speed (exact, same result as the per-person dot product):
//   with a_lj = C_lj / sd_l (sd signed) and A_j = -sum_l a_lj mu_l,
//   rhs_ij = C_{.j}' z~_i = A_j + sum_{l in M_i} a_lj mu_l + sum_{l het_i} a_lj + 2 sum_{l homalt_i} a_lj
//            + step sum_{(l, d) in D_i} a_lj d.
//   An observed dosage g enters as H + step d: H = round(g) (halves up) in the planes, d = round((g - H) / step)
//   (int8, |d| <= 0.5 / step) in the person's sparse list D_i of observed cells with d != 0; hard calls have
//   D_i empty. step (dosage_step) quantizes only this LD input; stored observed values are never changed.
//   People are processed in blocks; inside a block the loop runs over SNP columns j (column a_.j / C_.j hot in
//   cache) and over the people missing j, using 1-bit planes (het / homalt; la_kernel.cpp plane_dot), so C is
//   streamed from memory once per person block instead of once per (person, missing SNP). C_MM is gathered
//   from the same hot column.
//
// Pipeline: build_impute_ld -> read_block_genotypes (or geno_genovec) -> build_planes -> impute_stage1 -> impute_stage2,
// then for_each_cell visits every filled cell (i, j, g1, g). Errors are thrown as std::runtime_error.
#ifndef PLINK4R_IMPUTE_KERNEL_H
#define PLINK4R_IMPUTE_KERNEL_H

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <memory>
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

inline int tid() { return thread_id(); }

inline int team() { return thread_team(); }
}  // namespace ik

// ---- LD of one block: C, the scaled columns a_lj and A_j, the anchor mean / sd per row.
struct ImputeLD {
  int m = 0, k = 0, mp = 0;
  std::vector<int> vi;                  // 0-based PGEN variant per working row
  std::vector<float> mu, sd, isd, Aj;   // 2p, s sqrt(2p(1-p)) (signed), 1/sd, A_j
  Eigen::MatrixXf C, Aw;
};

// ld: the eigen block; rows: 0-based eigen rows of the working set S (panel rows present in the PGEN,
// increasing); vi, af, sign: per working row (PGEN variant, PGEN-ALT frequency in (0,1), +1 / -1).
inline ImputeLD build_impute_ld(const EigenBlock& ld, const std::vector<int>& rows, const std::vector<int>& vi,
                                const std::vector<double>& af, const std::vector<int>& sign, double delta,
                                int T) {
  using Eigen::MatrixXf;
  using Eigen::VectorXf;
  ImputeLD L;
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
  return L;
}

// ---- genotypes of one block: SNP-major genovecs of the selected samples (codes 0/1/2 ALT count, 3 missing)
// plus the dosage part of the observed cells: d = round((g - H) / step) != 0, per row j the selected samples
// di[dptr[j] .. dptr[j + 1]) (ascending) and their d.
struct BlockGenotypes {
  int m = 0, ns = 0, sb = 0;
  std::vector<uint8_t> GS;       // m x sb bytes
  std::vector<double> sum_obs;   // observed ALT sum per row (H + step d)
  double step = 0.01;            // dosage grid of the LD input
  std::vector<long long> dptr;   // m + 1
  std::vector<int> di;
  std::vector<int8_t> dd;
};

// d of an observed ALT dosage v (1/16384 units, 0..32768) with H = round(v / 16384) (halves up).
inline int dosage_grid(std::uint16_t v, unsigned h, double step) {
  const int q = static_cast<int>(v) - static_cast<int>(h << 14);
  const long d = std::lround(q / (16384.0 * step));
  return static_cast<int>(std::max(-127L, std::min(127L, d)));
}

// Gathers per-row (sample, d) lists into the CSR of B.
inline void set_dosage_lists(BlockGenotypes& B, std::vector<std::vector<int>>& li, std::vector<std::vector<int8_t>>& ld) {
  const int m = B.m;
  B.dptr.assign(m + 1, 0);
  for (int j = 0; j < m; ++j) B.dptr[j + 1] = B.dptr[j] + static_cast<long long>(li[j].size());
  B.di.resize(B.dptr[m]); B.dd.resize(B.dptr[m]);
  for (int j = 0; j < m; ++j) {
    std::copy(li[j].begin(), li[j].end(), B.di.begin() + B.dptr[j]);
    std::copy(ld[j].begin(), ld[j].end(), B.dd.begin() + B.dptr[j]);
    std::vector<int>().swap(li[j]); std::vector<int8_t>().swap(ld[j]);
  }
}

inline BlockGenotypes read_block_genotypes(PgenFile& pf, const std::vector<int>& vi, const std::vector<int>& sidx,
                                           double step, int T) {
  BlockGenotypes B;
  const int m = static_cast<int>(vi.size()), ns = static_cast<int>(sidx.size()), sb = (ns + 3) / 4;
  B.m = m; B.ns = ns; B.sb = sb; B.step = step;
  B.GS.assign(static_cast<std::size_t>(m) * sb, 0);
  B.sum_obs.assign(m, 0.0);
  std::vector<int> bad(T, -1);
  if (pf.dosage()) {  // H = round(dosage) into the codes, the rest on the step grid into the lists
    std::vector<std::vector<int>> li(m);
    std::vector<std::vector<int8_t>> ld(m);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
    for (int j = 0; j < m; ++j) {
      const std::uint16_t* x = pf.read_dosage(static_cast<std::size_t>(vi[j]), ik::tid());
      if (!x) { bad[ik::tid()] = j; continue; }
      uint8_t* dst = &B.GS[static_cast<std::size_t>(j) * sb];
      long long hs = 0, ds = 0;
      for (int i = 0; i < ns; ++i) {
        const std::uint16_t v = x[sidx[i]];
        unsigned code = 3u;
        if (v != PgenFile::kMissing) {
          code = (static_cast<unsigned>(v) + 8192u) >> 14;
          const int d = dosage_grid(v, code, step);
          if (d) { li[j].push_back(i); ld[j].push_back(static_cast<int8_t>(d)); ds += d; }
          hs += code;
        }
        dst[i >> 2] |= static_cast<uint8_t>(code << (2 * (i & 3)));
      }
      B.sum_obs[j] = static_cast<double>(hs) + step * static_cast<double>(ds);
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) throw std::runtime_error("PGEN read failed at eigen row " + std::to_string(bad[t] + 1));
    set_dosage_lists(B, li, ld);
    return B;
  }
  B.dptr.assign(m + 1, 0);
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
  return B;
}

// ---- bit-planes and missing lists of every person (chunks of 256 people = one cache line per SNP row)
struct PersonPlanes {
  int ns = 0, nbyte = 0;
  std::vector<uint8_t> het, hom;          // ns x nbyte
  std::vector<std::vector<int>> Mi;       // missing rows per person, ascending
  double step = 0.01;                     // dosage grid
  std::vector<long long> doff;            // ns + 1: person i's observed (row, d) pairs dl / dv[doff[i] .. doff[i + 1])
  std::vector<int> dl;
  std::vector<int8_t> dv;
};

// Consumes B.GS (released on return).
inline PersonPlanes build_planes(const ImputeLD& L, BlockGenotypes& B, int T) {
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
  // dosage lists, person-major (counting sort of the per-row lists)
  P.step = B.step;
  P.doff.assign(static_cast<std::size_t>(ns) + 1, 0);
  for (int i : B.di) ++P.doff[static_cast<std::size_t>(i) + 1];
  for (int s = 0; s < ns; ++s) P.doff[s + 1] += P.doff[s];
  P.dl.resize(B.di.size()); P.dv.resize(B.di.size());
  {
    std::vector<long long> pos(P.doff.begin(), P.doff.end() - 1);
    for (int j = 0; j < m; ++j)
      for (long long e = B.dptr[j]; e < B.dptr[j + 1]; ++e) {
        const long long t = pos[B.di[e]]++;
        P.dl[t] = j; P.dv[t] = B.dd[e];
      }
  }
  std::vector<int>().swap(B.di); std::vector<int8_t>().swap(B.dd);
  return P;
}

// ---- stage 1: the conditional mean of every missing cell, per thread lists (person, row, g1 clipped to [0,2])
struct Stage1Result {
  std::vector<std::vector<int>> oi, oj;
  std::vector<std::vector<float>> og;
  long long persons_solved = 0;
};

inline Stage1Result impute_stage1(const ImputeLD& L, const PersonPlanes& P, double delta, int nb_people,
                                  double mem_gb, int T) {
  using Eigen::MatrixXf;
  using Eigen::VectorXf;
  Stage1Result R;
  const int m = L.m, ns = P.ns, nbyte = P.nbyte;
  const MatrixXf& C = L.C;
  const MatrixXf& Aw = L.Aw;
  const std::vector<float>& mu = L.mu;
  const std::vector<float>& sd = L.sd;
  const std::vector<std::vector<int>>& MiA = P.Mi;
  const auto plane_dot = la().plane_dot;
  // S_i = I - C_MM: gathered from C in the per-person solve when C is small enough to stay in cache (no S
  // buffer, no scattered writes in the rhs pass); otherwise written column by column in the rhs pass, where
  // column j of C is hot.
  const bool gatherS = static_cast<double>(m) * m * 4.0 <= 32e6;
  // person blocks: fixed size, or as large as the memory budget for S_i, rhs and the SNP lists allows
  std::vector<int> bstart{0};
  {
    const double budget = mem_gb * 1e9;
    double acc = 0;
    for (int s = 0; s < ns; ++s) {
      const double q = static_cast<double>(MiA[s].size());
      const double cost = (gatherS ? 0.0 : 4.0 * q * q) + 36.0 * q;  // S_i, rhs, flat M_i, cell refs
      const bool full = nb_people > 0 ? (s - bstart.back() >= nb_people) : (acc + cost > budget && s > bstart.back());
      if (full) { bstart.push_back(s); acc = 0; }
      acc += cost;
    }
    bstart.push_back(ns);
  }
  const int nblk = static_cast<int>(bstart.size()) - 1;
  R.oi.assign(T, std::vector<int>()); R.oj.assign(T, std::vector<int>()); R.og.assign(T, std::vector<float>());
  // one entry per missing cell in SNP order, carrying everything the rhs loop needs about its person so the
  // loop reads it sequentially (no dependent loads of the person's list headers)
  struct CellRef { int i, a, q, mo; long long so; };
  std::vector<int> moff, jcnt(m + 1);
  std::unique_ptr<CellRef[]> cref_buf;
  std::unique_ptr<int[]> mflat_buf;
  std::size_t cref_cap = 0;
  std::vector<long long> soff;
  // rhs and S_i buffers: grown without zero fill (every entry is written before it is read; the first touch
  // then happens in the parallel rhs loop, not in a serial memset)
  std::unique_ptr<float[]> rhs_buf, S_buf;
  std::size_t rhs_cap = 0, S_cap = 0;
  float* rhs = nullptr;
  float* S = nullptr;
  for (int bk = 0; bk < nblk; ++bk) {
    const int b0 = bstart[bk], nb = bstart[bk + 1] - b0;
    const uint8_t* het = &P.het[static_cast<std::size_t>(b0) * nbyte];
    const uint8_t* hom = &P.hom[static_cast<std::size_t>(b0) * nbyte];
    const std::vector<int>* Mi = &MiA[b0];
    // offsets and per-SNP lists of (person, position in M)
    moff.assign(nb + 1, 0);
    soff.assign(nb + 1, 0);
    for (int i = 0; i < nb; ++i) {
      const int q = static_cast<int>(Mi[i].size());
      moff[i + 1] = moff[i] + q;
      soff[i + 1] = soff[i] + (gatherS ? 0LL : static_cast<long long>(q) * q);
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
        if (static_cast<std::size_t>(moff[nb]) > cref_cap) {
          cref_cap = moff[nb]; cref_buf.reset(new CellRef[cref_cap]); mflat_buf.reset(new int[cref_cap]);
        }
        int* mflat = mflat_buf.get();
        for (int i = i0; i < i1; ++i) std::copy(Mi[i].begin(), Mi[i].end(), mflat + moff[i]);
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
        }
        CellRef* cref = cref_buf.get();
        for (int i = i0; i < i1; ++i) {
          const int q = static_cast<int>(Mi[i].size());
          for (int a = 0; a < q; ++a) {
            const int p = c[Mi[i][a]]++;
            cref[p] = CellRef{i, a, q, moff[i], soff[i]};
          }
        }
      }
    }
    if (static_cast<std::size_t>(moff[nb]) > rhs_cap) { rhs_cap = moff[nb]; rhs_buf.reset(new float[rhs_cap]); }
    if (static_cast<std::size_t>(soff[nb]) > S_cap) { S_cap = static_cast<std::size_t>(soff[nb]); S_buf.reset(new float[S_cap]); }
    rhs = rhs_buf.get();  // fully overwritten below
    S = S_buf.get();      // every column of every S_i is written below

    // rhs and C_MM, column by column
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(T)
#endif
    for (int j = 0; j < m; ++j) {
      const int e0 = jcnt[j], e1 = jcnt[j + 1];
      if (e0 == e1) continue;
      const float* aj = Aw.col(j).data();
      const float* cj = C.col(j).data();
      const CellRef* cref = cref_buf.get();
      const int* mflat = mflat_buf.get();
      constexpr int kAhead = 6;  // prefetch the planes, M_i and S_i column of the cell kAhead ahead
      for (int e = e0; e < e1; ++e) {
        if (e + kAhead < e1) {
          const CellRef& cp = cref[e + kAhead];
          const uint8_t* hp = &het[static_cast<std::size_t>(cp.i) * nbyte];
          const uint8_t* op = &hom[static_cast<std::size_t>(cp.i) * nbyte];
          for (int b = 0; b < nbyte; b += 64) { __builtin_prefetch(hp + b); __builtin_prefetch(op + b); }
          __builtin_prefetch(mflat + cp.mo);
          if (!gatherS) __builtin_prefetch(&S[static_cast<std::size_t>(cp.so) + static_cast<std::size_t>(cp.a) * cp.q], 1);
          __builtin_prefetch(&rhs[cp.mo + cp.a], 1);
        }
        const CellRef& ce = cref[e];
        const int i = ce.i, a = ce.a, q = ce.q;
        const int* M = mflat + ce.mo;
        const uint8_t* hi = &het[static_cast<std::size_t>(i) * nbyte];
        const uint8_t* oi_ = &hom[static_cast<std::size_t>(i) * nbyte];
        float r = L.Aj[j] + plane_dot(aj, hi, oi_, nbyte);
        {  // dosage part of the observed cells
          const long long d0 = P.doff[b0 + i], d1 = P.doff[b0 + i + 1];
          if (d0 < d1) {
            float sd_ = 0;
            for (long long e = d0; e < d1; ++e) sd_ += aj[P.dl[e]] * static_cast<float>(P.dv[e]);
            r += static_cast<float>(P.step) * sd_;
          }
        }
        float r0 = 0, r1 = 0, r2 = 0, r3 = 0;  // independent partial sums (no single add chain)
        int bq = 0;
        for (; bq + 4 <= q; bq += 4) {
          const int l0 = M[bq], l1 = M[bq + 1], l2 = M[bq + 2], l3 = M[bq + 3];
          r0 += aj[l0] * mu[l0]; r1 += aj[l1] * mu[l1]; r2 += aj[l2] * mu[l2]; r3 += aj[l3] * mu[l3];
        }
        for (; bq < q; ++bq) { const int l = M[bq]; r0 += aj[l] * mu[l]; }
        r += (r0 + r1) + (r2 + r3);
        if (!gatherS) {  // column a of S_i
          float* Sc = &S[static_cast<std::size_t>(ce.so) + static_cast<std::size_t>(a) * q];
          for (int b = 0; b < q; ++b) Sc[b] = -cj[M[b]];
          Sc[a] += 1.0f;
        }
        rhs[ce.mo + a] = r;
      }
    }

    // per-person solve
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
    {
      const int th = ik::tid();
      std::vector<float> sbuf;  // this thread's S_i when gathered
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
      for (int i = 0; i < nb; ++i) {
        const std::vector<int>& M = Mi[i];
        const int q = static_cast<int>(M.size());
        if (q == 0) continue;
        float* sp;
        if (gatherS) {  // lower triangle and diagonal of I - C_MM (what LLT reads); the strict upper part too,
          sbuf.resize(static_cast<std::size_t>(q) * q);  // for the LDLT fallback below
          sp = sbuf.data();
          for (int x = 0; x < q; ++x) {
            const float* cx = C.col(M[x]).data();
            float* col = sp + static_cast<std::size_t>(x) * q;
            for (int y = 0; y < q; ++y) col[y] = -cx[M[y]];
            col[x] += 1.0f;
          }
        } else {
          sp = &S[static_cast<std::size_t>(soff[i])];
        }
        Eigen::Map<MatrixXf> Si(sp, q, q);
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
    for (int i = 0; i < nb; ++i) if (!Mi[i].empty()) ++R.persons_solved;
  }
  return R;
}

// ---- stage 2: anchor per SNP. cnt1 = filled cells of each row, shift = what stage 2 adds to them (0: stage 1 is kept).
struct Stage2Result {
  std::vector<double> cnt1, shift;
};

inline Stage2Result impute_stage2(const Stage1Result& R, const std::vector<double>& sum_obs,
                                  const std::vector<double>& af, int ns, double anchor_min_miss) {
  const int m = static_cast<int>(sum_obs.size());
  Stage2Result A;
  std::vector<double> sum1(m, 0.0);
  A.cnt1.assign(m, 0.0); A.shift.assign(m, 0.0);
  const int T = static_cast<int>(R.og.size());
  for (int t = 0; t < T; ++t)
    for (std::size_t r = 0; r < R.og[t].size(); ++r) { sum1[R.oj[t][r]] += R.og[t][r]; A.cnt1[R.oj[t][r]] += 1; }
  for (int j = 0; j < m; ++j) {
    const bool anchored = A.cnt1[j] > 0 && anchor_min_miss >= 0 && A.cnt1[j] >= anchor_min_miss * ns;
    if (!anchored) continue;
    const double mean1 = sum1[j] / A.cnt1[j];
    const double c = std::min(2.0, std::max(0.0, (2.0 * ns * af[j] - sum_obs[j]) / A.cnt1[j]));
    A.shift[j] = c - mean1;
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
