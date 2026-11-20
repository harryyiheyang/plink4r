// LD-bin list of the panel and the prior suppression of sbayesrc() built on it (layer below sbayesrc.cpp).
//
// Bin file (ld_dir/ldbin<r2>.bin): one int16 per panel SNP in panel order (block, ld_index), the SNP's bin id inside
// its block (0-based), -1 = in no bin. The panel order is the index, so no IDs are stored (2 bytes per SNP).
// Suppression (Yihe 2026-11-20): inside a bin with at least one annotated member, members without annotation get no
// prior mass on the nonzero classes (delta = 0) and leave the annotation regression.
#ifndef PLINK4R_LD_BINS_H
#define PLINK4R_LD_BINS_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace plink4r {

// bin ids of the kept SNPs: rows start[b] + keep[b][r] of the file, blocks concatenated; n_panel = SNPs in the file
inline std::vector<int> read_bins16(const std::string& path, int64_t n_panel, const std::vector<int64_t>& start,
                                    const std::vector<std::vector<int>>& keep) {
  std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
  if (!in) throw std::runtime_error("cannot open " + path);
  if ((int64_t)in.tellg() != 2 * n_panel)
    throw std::runtime_error("LD-bin file size does not match the panel (2 bytes per SNP): " + path);
  std::vector<int> out;
  for (std::size_t b = 0; b < keep.size(); ++b) {
    if (keep[b].empty()) continue;
    const int64_t lo = start[b] + keep[b].front(), hi = start[b] + keep[b].back();
    if (lo < 0 || hi >= n_panel) throw std::runtime_error("LD-bin rows outside the file: " + path);
    std::vector<int16_t> v(hi - lo + 1);
    in.seekg(2 * lo);
    in.read(reinterpret_cast<char*>(v.data()), 2 * (std::streamsize)v.size());
    if (!in) throw std::runtime_error("read error: " + path);
    for (int r : keep[b]) out.push_back(v[start[b] + r - lo]);
  }
  return out;
}

// suppressed SNPs of one block: no annotation in a bin with an annotated member (bid < 0: no bin)
inline void bin_suppress(int m, const int* bid, const char* annotated, char* out) {
  int nbins = 0;
  for (int j = 0; j < m; ++j) nbins = std::max(nbins, bid[j] + 1);
  std::vector<char> has(nbins, 0);
  for (int j = 0; j < m; ++j) if (annotated[j] && bid[j] >= 0) has[bid[j]] = 1;
  for (int j = 0; j < m; ++j) out[j] = !annotated[j] && bid[j] >= 0 && has[bid[j]];
}

}  // namespace plink4r
#endif
