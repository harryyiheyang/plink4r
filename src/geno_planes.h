// The 2-bit planes of a genotype matrix: two bit planes per variant, lo = [x == 1] and hi = [x == 2],
// with their packed-code encoder. The cross-product kernels on them are in geno_kernel.h.
#ifndef PLINK4R_GENO_PLANES_H
#define PLINK4R_GENO_PLANES_H

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

using u64 = std::uint64_t;
constexpr std::size_t kTile = 64;          // variants per tile side
constexpr std::size_t kChunk = 512;        // 64-bit words per cache chunk

struct Planes {
  std::size_t words = 0;              // words is a multiple of 4
  std::vector<u64> bits;              // variant j: lo at 2jW, hi at (2j+1)W
  std::vector<std::int64_t> s, v;     // sum x and n*sum x^2 - (sum x)^2
  const u64* lo(std::size_t j) const { return bits.data() + 2 * j * words; }
  const u64* hi(std::size_t j) const { return bits.data() + (2 * j + 1) * words; }
};

inline std::size_t plane_words(std::size_t n) { return ((n + 63) / 64 + 3) / 4 * 4; }

inline int popcount64(u64 x) { return __builtin_popcountll(x); }

inline void fill_median(std::size_t n, std::size_t words, u64* lo, u64* hi,
                        const u64* miss, std::int64_t& s, std::int64_t& v) {
  std::int64_t c1 = 0, c2 = 0, cm = 0;
  for (std::size_t w = 0; w < words; ++w) {
    c1 += popcount64(lo[w]);
    c2 += popcount64(hi[w]);
    cm += popcount64(miss[w]);
  }
  const std::int64_t observed = static_cast<std::int64_t>(n) - cm;
  const std::int64_t c0 = observed - c1 - c2;
  if (cm && observed) {
    u64* fill = 2 * c0 >= observed ? nullptr : (2 * (c0 + c1) >= observed ? lo : hi);
    if (fill) {
      for (std::size_t w = 0; w < words; ++w) fill[w] |= miss[w];
      (fill == lo ? c1 : c2) += cm;
    }
  }
  s = c1 + 2 * c2;
  v = static_cast<std::int64_t>(n) * (c1 + 4 * c2) - s * s;
}

// Packed 2-bit genotypes, 4 samples per byte with the first sample in the low
// bits (PLINK 1 BED columns and pgenlib genovecs). dosage[code] is the
// genotype x in {0, 1, 2} for that code, or -1 for missing.
struct NibbleTable {
  std::array<unsigned char, 256> lo{}, hi{}, miss{};
  explicit NibbleTable(const std::array<int, 4>& dosage) {
    for (unsigned int b = 0; b < 256; ++b)
      for (unsigned int k = 0; k < 4; ++k) {
        const int x = dosage[(b >> (2 * k)) & 3u];
        if (x == 2) hi[b] |= 1u << k;
        else if (x == 1) lo[b] |= 1u << k;
        else if (x < 0) miss[b] |= 1u << k;
      }
  }
};

// With miss_out, the missing plane is written there and left unfilled (s and
// v then describe the observed calls only); otherwise missing calls are filled
// with the variant's lower median.
inline void encode_packed(const NibbleTable& table, const unsigned char* column,
                          std::size_t n, std::size_t stride, std::size_t words,
                          u64* lo, u64* hi, std::int64_t& s, std::int64_t& v,
                          u64* miss_out = nullptr) {
  std::vector<u64> miss_local(miss_out ? 0 : words, 0);
  u64* miss = miss_out ? miss_out : miss_local.data();
  for (std::size_t w = 0; w < words; ++w) {
    u64 l = 0, h = 0, m = 0;
    const std::size_t b0 = 16 * w, b1 = std::min(stride, b0 + 16);
    for (std::size_t b = b0; b < b1; ++b) {
      const unsigned char byte = column[b];
      const unsigned int shift = 4 * (b - b0);
      l |= static_cast<u64>(table.lo[byte]) << shift;
      h |= static_cast<u64>(table.hi[byte]) << shift;
      m |= static_cast<u64>(table.miss[byte]) << shift;
    }
    lo[w] = l;
    hi[w] = h;
    miss[w] = m;
  }
  const std::size_t full = n / 64;
  if (n % 64) {
    const u64 mask = (u64(1) << (n % 64)) - 1;
    lo[full] &= mask;
    hi[full] &= mask;
    miss[full] &= mask;
  }
  if (miss_out) {
    std::int64_t c1 = 0, c2 = 0;
    for (std::size_t w = 0; w < words; ++w) {
      c1 += popcount64(lo[w]);
      c2 += popcount64(hi[w]);
    }
    s = c1 + 2 * c2;
    v = c1 + 4 * c2;  // sum of squares of the observed calls
    return;
  }
  fill_median(n, words, lo, hi, miss, s, v);
}

inline void reset_planes(Planes& out, std::size_t n, std::size_t m) {
  out.words = plane_words(n);
  out.bits.assign(2 * m * out.words, 0);
  out.s.assign(m, 0);
  out.v.assign(m, 0);
}

#endif
