// High-LD bins from a pairwise LD list (GCTB GWFM .pwld: SNP1 SNP2 LDcorrelation, each pair once) and the prior
// suppression of sbayesrc() built on them (layer below sbayesrc.cpp).
//
// Bins (ld_bins.py, Yihe 2026-11-09): neighbours are the pairs with r^2 > t. Greedy over the SNPs in panel order:
// repeatedly take the unassigned SNP with the most unassigned neighbours (itself included; ties: the first in panel
// order) as the bin centre; its bin = itself + those neighbours. SNPs without a pair are singletons.
// Suppression (Yihe 2026-11-20): inside a bin with at least one annotated member, members without annotation get no
// prior mass on the nonzero classes (delta = 0) and leave the annotation regression.
#ifndef PLINK4R_LD_BINS_H
#define PLINK4R_LD_BINS_H

#include <algorithm>
#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

namespace plink4r {

// neighbour lists over M SNPs from pairs (i, j) (0-based, -1 = SNP not kept) with r^2 > t; duplicates and self pairs
// dropped
inline std::vector<std::vector<int>> pair_neighbours(int64_t M, const int* pi, const int* pj, const double* r,
                                                     int64_t npair, double t) {
  std::vector<std::vector<int>> nb(M);
  for (int64_t p = 0; p < npair; ++p) {
    const int i = pi[p], j = pj[p];
    if (i < 0 || j < 0 || i == j || !(r[p] * r[p] > t)) continue;
    nb[i].push_back(j);
    nb[j].push_back(i);
  }
  for (auto& v : nb) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }
  return nb;
}

// bin id (0..) per SNP, greedy as ld_bins.block_bins
inline std::vector<int> greedy_bins(const std::vector<std::vector<int>>& nb) {
  const int64_t m = (int64_t)nb.size();
  std::vector<int> bid(m, -1), cnt(m);
  for (int64_t j = 0; j < m; ++j) cnt[j] = (int)nb[j].size() + 1;
  std::priority_queue<std::pair<int, int64_t>> pq;  // (cnt, -index): max count, then first index; stale entries skipped
  for (int64_t j = 0; j < m; ++j) pq.push({cnt[j], -j});
  int nbins = 0;
  std::vector<int64_t> mem;
  while (!pq.empty()) {
    const auto top = pq.top();
    pq.pop();
    const int64_t c = -top.second;
    if (bid[c] >= 0 || top.first != cnt[c]) continue;
    mem.assign(1, c);
    for (int i : nb[c]) if (bid[i] < 0) mem.push_back(i);
    for (int64_t i : mem) bid[i] = nbins;
    ++nbins;
    for (int64_t i : mem) {
      --cnt[i];
      for (int h : nb[i]) {
        --cnt[h];
        if (bid[h] < 0) pq.push({cnt[h], -(int64_t)h});
      }
    }
  }
  return bid;
}

// suppressed SNPs: no annotation in a bin with an annotated member
inline std::vector<char> bin_suppress(const std::vector<int>& bid, const std::vector<char>& annotated) {
  const int64_t m = (int64_t)bid.size();
  int nbins = 0;
  for (int b : bid) nbins = std::max(nbins, b + 1);
  std::vector<char> has(nbins, 0), out(m, 0);
  for (int64_t j = 0; j < m; ++j) if (annotated[j]) has[bid[j]] = 1;
  for (int64_t j = 0; j < m; ++j) out[j] = !annotated[j] && has[bid[j]];
  return out;
}

}  // namespace plink4r
#endif
