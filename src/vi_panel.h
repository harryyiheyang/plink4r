// Second layer of the summary-statistics PRS methods: the genome-wide eigen panel in memory, restricted to the SNPs
// with data, in the form the coordinate updates read (built from eigen_bin.h and vi_store.h).
//
// Per block b with eigenpairs (lambda, U) of R = U diag(lambda) U' (block<N>.eigen.bin, either layout) and kept rows
// keep_b (0-based eigen rows with data, increasing):
//   V_b   = U[keep_b, ] diag(lambda)^{1/2}, one row of length k_b per SNP, stored as T (float or float16)
//   d_j   = V_j . V_j                 (diagonal of V V' from the stored rows)
//   l_j   = sum_i lambda_i V_ji^2     (LD score of the eigen LD)
//   w_b   = Lambda^{-1} V_b' b_b      (SBayesRC eigen-space data, = Lambda^{-1/2} U' b; 0 where lambda <= 1e-8)
//   bt_j  = V_j . w_b                 (= (U U' b)_j)
// d, l, w and bt are computed from the stored (rounded) rows, so the data term and the LD of the fit are one matrix.
// Ported from viprs_eigen/src/_engine.cpp (load_one / load_all); blocks are read in parallel, largest first.
#ifndef PLINK4R_VI_PANEL_H
#define PLINK4R_VI_PANEL_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"
#include "vi_store.h"

namespace plink4r {

template <typename T>
struct EigenPanel {
  int nb = 0;
  int64_t M = 0, S = 0;
  std::vector<int64_t> snp_off, v_off, s_off;  // nb + 1 / nb / nb + 1
  std::vector<int> k;                          // eigen dimension per block
  std::vector<int> order;                      // blocks by decreasing nk * k (sweep scheduling)
  std::vector<T> V;                            // concatenated row-major V_b
  std::vector<float> d, ldscore, lam, w, bt;   // M / M / S / S / M

  const T* row(int b, int64_t j) const { return V.data() + v_off[b] + (j - snp_off[b]) * (int64_t)k[b]; }
};

// m and the eigen dimension k kept by read_eigen_block(path, thresh), from the header and eigenvalues only.
inline void eigen_block_dims(const std::string& path, double thresh, int& m, int& k) {
  std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
  if (!in) throw std::runtime_error("cannot open " + path);
  const std::uint64_t size = static_cast<std::uint64_t>(in.tellg());
  in.seekg(0);
  int32_t mm = 0, kf = 0;
  read_n(in, &mm, 1, path);
  read_n(in, &kf, 1, path);
  if (mm <= 0 || kf <= 0 || kf > mm) throw std::runtime_error("bad dimensions: " + path);
  std::vector<double> lam(kf);
  double sum_lambda = 0;
  float file_thresh = 0;
  if (size == 20ull + 12ull * kf + 2ull * mm * kf) {
    read_n(in, &sum_lambda, 1, path);
    read_n(in, &file_thresh, 1, path);
    read_n(in, lam.data(), kf, path);
  } else if (size == 16ull + 4ull * kf + 4ull * mm * kf) {
    float s = 0;
    read_n(in, &s, 1, path);
    sum_lambda = s;
    read_n(in, &file_thresh, 1, path);
    std::vector<float> lf(kf);
    read_n(in, lf.data(), kf, path);
    for (int j = 0; j < kf; ++j) lam[j] = lf[j];
  } else {
    throw std::runtime_error("unknown eigen.bin layout (size mismatch): " + path);
  }
  m = mm;
  k = kf;
  if (thresh > 0) {
    if (thresh > file_thresh + 1e-6) throw std::runtime_error("thresh exceeds the file cutoff: " + path);
    double cum = 0;
    for (int j = 0; j < kf; ++j) {
      cum += lam[j];
      if (cum >= thresh * sum_lambda) { k = j + 1; break; }
    }
  }
}

// blocks: block ids; keep[b]: kept eigen rows of block b (increasing); bstd: standardized marginal effect of every kept
// SNP, blocks concatenated in the same order. thresh: eigen truncation of read_eigen_block (<= 0 keeps every stored
// pair). Errors name the block.
template <typename T>
EigenPanel<T> load_eigen_panel(const std::vector<std::string>& files, const std::vector<int>& blocks,
                               const std::vector<std::vector<int>>& keep, const std::vector<float>& bstd,
                               double thresh, int threads) {
  EigenPanel<T> P;
  const int nb = static_cast<int>(blocks.size());
  if (files.size() != blocks.size()) throw std::runtime_error("one eigen.bin file per block needed");
  P.nb = nb;
  P.snp_off.assign(nb + 1, 0);
  for (int b = 0; b < nb; ++b) P.snp_off[b + 1] = P.snp_off[b] + static_cast<int64_t>(keep[b].size());
  P.M = P.snp_off[nb];
  if (static_cast<int64_t>(bstd.size()) != P.M) throw std::runtime_error("bstd length differs from the kept SNPs");
  // pass 1: eigen dimension after truncation, from the header and eigenvalues only (no U is held)
  std::vector<int> kdim(nb), mdim(nb);
  std::vector<std::string> err(nb);
  for (int b = 0; b < nb; ++b) {
    try {
      eigen_block_dims(files[b], thresh, mdim[b], kdim[b]);
      for (std::size_t r = 0; r < keep[b].size(); ++r)
        if (keep[b][r] < 0 || keep[b][r] >= mdim[b] || (r && keep[b][r] <= keep[b][r - 1]))
          throw std::runtime_error("ld_index outside the block or not increasing");
    } catch (const std::exception& e) {
      err[b] = e.what();
    }
  }
  for (int b = 0; b < nb; ++b)
    if (!err[b].empty()) throw std::runtime_error("block " + std::to_string(blocks[b]) + ": " + err[b]);
  P.k.resize(nb);
  P.v_off.resize(nb);
  P.s_off.assign(nb + 1, 0);
  int64_t tot = 0;
  for (int b = 0; b < nb; ++b) {
    P.k[b] = kdim[b];
    P.v_off[b] = tot;
    tot += (P.snp_off[b + 1] - P.snp_off[b]) * (int64_t)P.k[b];
    P.s_off[b + 1] = P.s_off[b] + P.k[b];
  }
  P.S = P.s_off[nb];
  P.V.resize(tot);
  P.d.resize(P.M);
  P.ldscore.resize(P.M);
  P.lam.resize(P.S);
  P.w.resize(P.S);
  P.bt.resize(P.M);
  P.order.resize(nb);
  for (int b = 0; b < nb; ++b) P.order[b] = b;
  std::stable_sort(P.order.begin(), P.order.end(), [&](int a, int c) {
    return (P.snp_off[a + 1] - P.snp_off[a]) * (int64_t)P.k[a] > (P.snp_off[c + 1] - P.snp_off[c]) * (int64_t)P.k[c];
  });
  // pass 2: read, rotate the kept rows, data terms (each block freed once rotated)
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int o = 0; o < nb; ++o) {
    const int b = P.order[o];
    EigenBlock E;
    try {
      E = read_eigen_block(files[b], thresh);
      if (E.k != P.k[b]) throw std::runtime_error("eigen dimension changed between reads");
    } catch (const std::exception& e) {
      err[b] = e.what();
      continue;
    }
    const int kb = P.k[b];
    const int64_t o0 = P.snp_off[b], nk = P.snp_off[b + 1] - o0;
    std::vector<float> lamc(kb), sq(kb);
    for (int i = 0; i < kb; ++i) {
      lamc[i] = std::max(E.lambda[i], 0.f);
      sq[i] = std::sqrt(lamc[i]);
      P.lam[P.s_off[b] + i] = lamc[i];
    }
    std::vector<double> wacc(kb, 0.);
    T* Vb = P.V.data() + P.v_off[b];
    for (int64_t r = 0; r < nk; ++r) {
      const int j = keep[b][r];
      T* row = Vb + r * (int64_t)kb;
      const double br = bstd[o0 + r];
      double dd = 0., ll = 0.;
      for (int i = 0; i < kb; ++i) {
        const float v = RowStore<T>::put(row, i, E.U(j, i) * sq[i]);
        dd += (double)v * v;
        ll += (double)lamc[i] * v * v;
        wacc[i] += (double)v * br;
      }
      P.d[o0 + r] = (float)dd;
      P.ldscore[o0 + r] = (float)ll;
    }
    float* wb = P.w.data() + P.s_off[b];
    for (int i = 0; i < kb; ++i) wb[i] = (float)(lamc[i] > 1e-8f ? wacc[i] / lamc[i] : 0.);
    for (int64_t r = 0; r < nk; ++r) P.bt[o0 + r] = RowStore<T>::dot(Vb + r * (int64_t)kb, wb, kb);
  }
  for (int b = 0; b < nb; ++b)
    if (!err[b].empty()) throw std::runtime_error("block " + std::to_string(blocks[b]) + ": " + err[b]);
  return P;
}

}  // namespace plink4r
#endif
