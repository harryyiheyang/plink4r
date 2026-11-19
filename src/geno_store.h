// geno: a PLINK BED / PGEN genotype matrix held in memory as X = H + D with
// the products SuSiE-style fine-mapping needs: X'X by exact popcount plus the
// D terms, X'M and X B by plane scans, and dense sub-blocks on demand.
// H is a 0/1/2 matrix in 2-bit planes (n p / 4 bytes). D is a sparse
// deviation: per variant a 1-bit nonzero mask and int16 values in units of
// that variant's `unit` (1/16384 for hard-call inputs, 0.01 for a dosage PGEN),
// allocated only for variants with a nonzero D. It holds
//   - dosage PGEN: dosage - H, the file's 1/16384 dosage taken to the nearest 0.01 (PI 2026-11-18: the
//     precision of imputed dosages) and H that value rounded;
//   - every input: the fill of a missing entry, f - H with H = round(f) and
//     f the variant's lower median (impute = "median", D = 0) or observed
//     mean (impute = "mean", and dosage PGEN; quantized to the variant's unit);
//   - appended columns (geno_append_cpp): all of x, H = 0, own unit.
// Values are A1 (BED) / ALT (PGEN) dosages as in BEDMatrix. Centering and
// scaling are applied on the R side from mu and sd of the stored X.
// The rows of the cells that were missing in the file are kept per variant
// (M, ascending); a missing cell stores its current fill as H + D like any
// other cell. geno_fill_cpp rewrites those cells (the one fill routine),
// geno_take_cpp copies rows / columns with an optional exact flip x -> 2 - x.

#ifndef PLINK4R_GENO_STORE_H
#define PLINK4R_GENO_STORE_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include "geno_planes.h"
#include "geno_access.h"
#include "pgen_reader.h"
#include "resources.h"
#include "la_dispatch.h"
#include <R_ext/Rdynload.h>

// The store itself (struct Geno) and its cell accessors, shared by geno_store.cpp (opening, editing, filling),
// geno_products.cpp (the products) and geno_api.cpp (the C-callable entry points).
namespace plink4r {
namespace gstore {


using plink4r::PgenFile;
// PGEN genovec codes as ALT dosage for the geno object: 0 hom REF, 1 het, 2 two non-REF alleles, 3 missing.
// This matches the BED table of geno (A1 dosage, as BEDMatrix), since plink2 imports BED with A1 as ALT.
inline const NibbleTable& pgen_table() {
  static const NibbleTable t({0, 1, 2, -1});
  return t;
}

constexpr std::size_t kGroup = 64;          // variants per read group
constexpr std::size_t kDBlock = 8;          // mask words per prefix count (512 rows)
constexpr double kDUnit = 1.0 / 16384.0;    // value of one int16 step of D, hard-call inputs
constexpr double kDosUnit = 0.01;           // value of one int16 step of D, dosage PGEN
constexpr int kDosSteps = 100;              // steps of kDosUnit per allele

// int16 steps per allele of a genotype column's unit (16384 or 100); 0 for an appended column.
inline int unit_steps(double unit) {
  return unit == kDUnit ? 16384 : unit == kDosUnit ? kDosSteps : 0;
}

// D of one variant: nonzero mask over the samples (H words), int16 values
// q = D / unit in ascending row order, and the nonzero count before each block
// of kDBlock mask words. All empty when the variant's D is zero. unit is
// 1/16384 (hard-call inputs) or 0.01 (dosage PGEN) for genotype columns; an
// appended column (H = 0, all of x in D) has its own unit max|x| / 32767.
struct DCol {
  std::vector<u64> mask;
  std::vector<std::int16_t> q;
  std::vector<std::uint32_t> pre;
  double unit = kDUnit;
};

struct Geno {
  std::size_t n = 0, p = 0;
  bool mean = false;            // impute = "mean"
  bool dosage = false;          // PGEN with dosages
  Planes G;                     // H: lo = [h == 1], hi = [h == 2]
  std::vector<DCol> D;
  std::vector<double> dsum, dsq, dhd;  // column sums of D, D^2 and H D
  std::vector<double> mu, sd;   // column mean and sd (n - 1) of X = H + D
  std::vector<double> nmiss;    // missing calls per variant before filling
  std::vector<std::vector<std::uint32_t>> M;  // rows of those missing calls per variant, ascending
  std::size_t dnnz() const {
    std::size_t s = 0;
    for (const DCol& c : D) s += c.q.size();
    return s;
  }
};

// Appends nonzero D entries of variant j in ascending row order and keeps its
// column sums; finish() stores the sums and builds the prefix counts.
struct DBuilder {
  Geno& g;
  std::size_t j;
  double ds = 0.0, dsq = 0.0, dh = 0.0;
  DBuilder(Geno& g_, std::size_t j_) : g(g_), j(j_) {}
  void add(std::size_t i, unsigned int h, int q) {
    DCol& c = g.D[j];
    if (c.mask.empty()) c.mask.assign(g.G.words, 0);
    c.mask[i >> 6] |= u64(1) << (i & 63u);
    c.q.push_back(static_cast<std::int16_t>(q));
    const double dd = q * g.D[j].unit;
    ds += dd;
    dsq += dd * dd;
    dh += h * dd;
  }
  void finish() {
    DCol& c = g.D[j];
    g.dsum[j] = ds;
    g.dsq[j] = dsq;
    g.dhd[j] = dh;
    if (c.mask.empty()) return;
    c.q.shrink_to_fit();
    c.pre.assign((c.mask.size() + kDBlock - 1) / kDBlock, 0);
    std::uint32_t cnt = 0;
    for (std::size_t w = 0; w < c.mask.size(); ++w) {
      if (w % kDBlock == 0) c.pre[w / kDBlock] = cnt;
      cnt += static_cast<std::uint32_t>(popcount64(c.mask[w]));
    }
  }
};

// A fill value f as H = round(f) (halves up, within 0..2) and the int16
// deviation q = round((f - H) / unit).
inline void quantize_fill(double f, double unit, unsigned int& h, int& q) {
  const double hf = std::min(2.0, std::max(0.0, std::floor(f + 0.5)));
  h = static_cast<unsigned int>(hf);
  q = static_cast<int>(std::lround((f - hf) / unit));
}

// Nonzeros of D column c in mask words before word w.
inline std::size_t d_rank(const DCol& c, std::size_t w) {
  const std::size_t blk = w / kDBlock;
  std::size_t k = blk < c.pre.size() ? c.pre[blk] : c.q.size();
  for (std::size_t x = blk * kDBlock; x < w; ++x) k += popcount64(c.mask[x]);
  return k;
}

// int16 value q of D column c at sample i (0 when not stored).
inline int q_at(const DCol& c, std::size_t i) {
  if (c.mask.empty()) return 0;
  const std::size_t w = i / 64, b = i % 64;
  if (!((c.mask[w] >> b) & 1u)) return 0;
  const std::size_t e = d_rank(c, w) + popcount64(c.mask[w] & ((u64(1) << b) - 1));
  return c.q[e];
}

// D of column c at sample i.
inline double d_at(const DCol& c, std::size_t i) {
  const int q = q_at(c, i);
  return q ? q * c.unit : 0.0;
}

// H of column j at sample i.
inline unsigned int h_at(const Planes& G, std::size_t j, std::size_t i) {
  const u64 bit = u64(1) << (i % 64);
  return (G.lo(j)[i / 64] & bit) ? 1u : (G.hi(j)[i / 64] & bit) ? 2u : 0u;
}


inline int resolve_threads(int threads) {
#ifdef _OPENMP
  return plink4r::resolve_threads(threads);
#else
  return 1;
#endif
}

// Sum of the stored values of column j over all cells, in int16 steps of a genotype column's unit (k per allele).
inline std::int64_t step_sum_col(const Geno& g, std::size_t j, int k) {
  std::int64_t tot = k * static_cast<std::int64_t>(g.G.s[j]);
  for (std::int16_t q : g.D[j].q) tot += q;
  return tot;
}

}  // namespace gstore
}  // namespace plink4r

#endif
