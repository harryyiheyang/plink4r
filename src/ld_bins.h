// High-LD bins of one int16 eigen block and the prior suppression of sbayesrc() built on them (layer above
// eigen_int16.h).
//
// Bins (ld_bins.py, Yihe 2026-11-09): R = X X' over the kept rows, X = Q diag(c) = U Lambda^{1/2} at the stored k,
// r2_ij = R_ij^2 / (R_ii R_jj). Greedy: repeatedly take the unassigned SNP with the most unassigned neighbours at
// r2 > t (ties: the first in panel order) as the bin centre; its bin = itself + those neighbours. Singletons stay alone.
// R is formed in row tiles by the package GEMM (la_dispatch.h) and only the neighbour lists are kept; each pair is
// decided once, so the lists are symmetric.
//
// Suppression (Yihe 2026-11-20): inside a bin with at least one annotated member, members without annotation get
// no prior mass on the nonzero classes (delta = 0) and leave the annotation regression.
#ifndef PLINK4R_LD_BINS_H
#define PLINK4R_LD_BINS_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <utility>
#include <vector>
#include "eigen_int16.h"
#include "la_dispatch.h"

namespace plink4r {

// neighbour lists (r2 > t, self excluded) of the kept rows of one block; tile = rows per GEMM call
inline std::vector<std::vector<int>> block_neighbours(const Block16& B, double t, int tile) {
  const int m = B.nk, k = B.ld;
  std::vector<std::vector<int>> nb(m);
  if (m < 2) return nb;
  std::vector<float> Xr((std::size_t)m * k), Xc((std::size_t)m * k);  // row-major (= col-major k x m), col-major m x k
  for (int j = 0; j < m; ++j)
    for (int s = 0; s < k; ++s) {
      const float x = (float)B.Q[(std::size_t)j * k + s] * B.c[s];
      Xr[(std::size_t)j * k + s] = x;
      Xc[(std::size_t)s * m + j] = x;
    }
  std::vector<double> dg(m);
  for (int j = 0; j < m; ++j) {
    double a = 0;
    for (int s = 0; s < k; ++s) a += (double)Xr[(std::size_t)j * k + s] * Xr[(std::size_t)j * k + s];
    dg[j] = std::max(a, 1e-12);
  }
  const int T = std::max(1, std::min(tile, m));
  std::vector<float> R((std::size_t)m * T);
  for (int j0 = 0; j0 < m; j0 += T) {
    const int nt = std::min(T, m - j0);
    la().sgemm_nn(Xc.data(), m, k, m, Xr.data() + (std::size_t)j0 * k, nt, k, R.data(), m);  // R = X X[j0:j0+nt]'
    for (int jj = 0; jj < nt; ++jj) {
      const int j = j0 + jj;
      const float* col = R.data() + (std::size_t)jj * m;
      for (int i = j + 1; i < m; ++i) {
        const double r = col[i];
        if (r * r > t * dg[i] * dg[j]) { nb[j].push_back(i); nb[i].push_back(j); }
      }
    }
  }
  return nb;
}

// bin id (0..) per kept row, greedy as ld_bins.block_bins
inline std::vector<int> block_bins(const std::vector<std::vector<int>>& nb) {
  const int m = (int)nb.size();
  std::vector<int> bid(m, -1), cnt(m);
  for (int j = 0; j < m; ++j) cnt[j] = (int)nb[j].size() + 1;  // self included, as r2_jj = 1 > t
  // max-heap on (cnt, -index); stale entries are skipped
  std::priority_queue<std::pair<int, int>> pq;
  for (int j = 0; j < m; ++j) pq.push({cnt[j], -j});
  int nbins = 0;
  std::vector<int> mem;
  while (!pq.empty()) {
    const auto top = pq.top();
    pq.pop();
    const int c = -top.second;
    if (bid[c] >= 0 || top.first != cnt[c]) continue;
    mem.assign(1, c);
    for (int i : nb[c]) if (bid[i] < 0) mem.push_back(i);
    for (int i : mem) bid[i] = nbins;
    ++nbins;
    for (int i : mem) {
      --cnt[i];
      for (int h : nb[i]) {
        --cnt[h];
        if (bid[h] < 0) pq.push({cnt[h], -h});
      }
    }
  }
  return bid;
}

// suppressed rows of one block: no annotation (annotated[j] false) in a bin with an annotated member
inline std::vector<char> block_suppress(const std::vector<int>& bid, const std::vector<char>& annotated) {
  const int m = (int)bid.size();
  int nbins = 0;
  for (int b : bid) nbins = std::max(nbins, b + 1);
  std::vector<char> has(nbins, 0), out(m, 0);
  for (int j = 0; j < m; ++j) if (annotated[j]) has[bid[j]] = 1;
  for (int j = 0; j < m; ++j) out[j] = !annotated[j] && has[bid[j]];
  return out;
}

}  // namespace plink4r
#endif
