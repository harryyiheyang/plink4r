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

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include "geno_kernel.h"
#include "geno_access.h"
#include "pgen_reader.h"
#include "resources.h"
#include "la_dispatch.h"
#include <R_ext/Rdynload.h>

namespace {

using plink4r::PgenFile;
// PGEN genovec codes as ALT dosage for the geno object: 0 hom REF, 1 het, 2 two non-REF alleles, 3 missing.
// This matches the BED table of geno (A1 dosage, as BEDMatrix), since plink2 imports BED with A1 as ALT.
const NibbleTable& pgen_table() {
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

// Packed column of the selected samples, in their requested order.
inline void gather_packed(const unsigned char* column, const std::vector<int>& rows,
                          unsigned char* buf) {
  std::fill(buf, buf + (rows.size() + 3) / 4, 0);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const unsigned int r = static_cast<unsigned int>(rows[i]);
    const unsigned int code = (column[r >> 2] >> (2 * (r & 3u))) & 3u;
    buf[i >> 2] |= static_cast<unsigned char>(code << (2 * (i & 3u)));
  }
}

// Fills the missing entries (mask miss) of variant j with the constant f:
// round(f) into the planes, the rest into D; records their rows in M[j].
void fill_missing(Geno& g, std::size_t j, const u64* miss, double f, bool fill, DBuilder& db) {
  const std::size_t W = g.G.words;
  u64* lo = g.G.bits.data() + 2 * j * W;
  u64* hi = lo + W;
  unsigned int h;
  int q;
  quantize_fill(f, g.D[j].unit, h, q);
  std::vector<std::uint32_t>& rows = g.M[j];
  for (std::size_t w = 0; w < W; ++w)
    for (u64 bits = miss[w]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * w + __builtin_ctzll(bits);
      rows.push_back(static_cast<std::uint32_t>(i));
      if (!fill) continue;
      if (h == 1) lo[w] |= u64(1) << (i & 63u);
      else if (h == 2) hi[w] |= u64(1) << (i & 63u);
      if (q) db.add(i, h, q);
    }
}

// G.s[j] = sum h and G.v[j] = n sum h^2 - (sum h)^2 from the planes.
void plane_stats(Geno& g, std::size_t j) {
  std::int64_t c1 = 0, c2 = 0;
  const u64* lo = g.G.lo(j);
  const u64* hi = g.G.hi(j);
  for (std::size_t w = 0; w < g.G.words; ++w) {
    c1 += popcount64(lo[w]);
    c2 += popcount64(hi[w]);
  }
  const std::int64_t s = c1 + 2 * c2;
  g.G.s[j] = s;
  g.G.v[j] = static_cast<std::int64_t>(g.n) * (c1 + 4 * c2) - s * s;
}

// Hardcall column: planes of the observed calls, then the missing calls filled
// with the variant's lower median or observed mean (left 0 when every call is
// missing); their rows go to M[j]. Other fills: geno_fill_cpp afterwards.
void encode_column(Geno& g, const NibbleTable& table, const unsigned char* buf,
                   std::size_t j, u64* miss) {
  const std::size_t W = g.G.words, n = g.n;
  u64* lo = g.G.bits.data() + 2 * j * W;
  std::int64_t s, v;
  encode_packed(table, buf, n, (n + 3) / 4, W, lo, lo + W, s, v, miss);
  std::int64_t cm = 0;
  for (std::size_t w = 0; w < W; ++w) cm += popcount64(miss[w]);
  const std::int64_t obs = static_cast<std::int64_t>(n) - cm;
  g.nmiss[j] = static_cast<double>(cm);
  DBuilder db(g, j);
  if (cm) {
    double f = 0.0;
    if (obs) {
      std::int64_t c1 = 0, c2 = 0;
      for (std::size_t w = 0; w < W; ++w) {
        c1 += popcount64(lo[w]);
        c2 += popcount64(lo[W + w]);
      }
      const std::int64_t c0 = obs - c1 - c2;
      f = g.mean ? static_cast<double>(c1 + 2 * c2) / static_cast<double>(obs)
        : (2 * c0 >= obs ? 0.0 : (2 * (c0 + c1) >= obs ? 1.0 : 2.0));
    }
    g.M[j].reserve(static_cast<std::size_t>(cm));
    fill_missing(g, j, miss, f, obs > 0, db);
  }
  db.finish();
  plane_stats(g, j);
}

void init_geno(Geno& g, std::size_t n, std::size_t p, bool mean) {
  if (n < 2) Rcpp::stop("The genotype matrix needs at least two samples.");
  if (!p) Rcpp::stop("No variants selected.");
  g.n = n;
  g.p = p;
  g.mean = mean;
  reset_planes(g.G, n, p);
  g.D.assign(p, DCol());
  g.dsum.assign(p, 0.0);
  g.dsq.assign(p, 0.0);
  g.dhd.assign(p, 0.0);
  g.nmiss.assign(p, 0.0);
  g.M.assign(p, std::vector<std::uint32_t>());
}

// mu[j] and sd[j] (n - 1) of X = H + D from the plane and D sums.
void col_stats(Geno& g, std::size_t j) {
  const double n = static_cast<double>(g.n);
  const double s = static_cast<double>(g.G.s[j]);
  const double v = static_cast<double>(g.G.v[j]);
  const double mean = (s + g.dsum[j]) / n;
  const double sumsq = (v + s * s) / n + 2.0 * g.dhd[j] + g.dsq[j];
  g.mu[j] = mean;
  const double var = (sumsq - n * mean * mean) / (n - 1);
  g.sd[j] = var > 0 ? std::sqrt(var) : 0.0;
}

void finish_stats(Geno& g) {
  g.mu.assign(g.p, 0.0);
  g.sd.assign(g.p, 0.0);
  for (std::size_t j = 0; j < g.p; ++j) col_stats(g, j);
}

std::vector<int> to_rows(const Rcpp::IntegerVector& x) {
  return std::vector<int>(x.begin(), x.end());
}

// Dosage PGEN: the file's dosage (1/16384 steps) is taken to the nearest 0.01 (plink4r::dosage_grid); H = that
// value rounded to the nearest integer (halves up) into the planes, D = value - H as int16 in [-50, 49] steps
// of 0.01. Missing dosages take the variant's observed mean of those values, through D; their rows go to M[j].
void load_dosage(Geno& g, PgenFile& file, const Rcpp::IntegerVector& variants,
                 const std::vector<int>& rows, int threads, const std::string& path) {
  const std::size_t n = g.n, p = g.p, sel_stride = (n + 3) / 4, W = g.G.words;
  g.dosage = true;
  const NibbleTable& table = pgen_table();
  const std::uint16_t kMiss = PgenFile::kMissing;
  std::vector<unsigned char> packed(static_cast<std::size_t>(threads) * sel_stride);
  for (std::size_t g0 = 0; g0 < p; g0 += kGroup) {
    const std::size_t count = std::min(kGroup, p - g0);
    int failed = 0, range = 0;
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const int t = plink4r::thread_id();
      const std::size_t j = g0 + static_cast<std::size_t>(k);
      const std::uint16_t* x = file.read_dosage(static_cast<std::size_t>(variants[j]), t);
      if (!x) {
        #pragma omp atomic write
        failed = 1;
        continue;
      }
      std::uint64_t obs_sum = 0;
      std::size_t miss = 0;
      bool bad_range = false;
      for (std::size_t i = 0; i < n; ++i) {
        const std::uint16_t v = x[rows[i]];
        if (v == kMiss) ++miss;
        else if (v > 32768u) bad_range = true;
        else obs_sum += plink4r::dosage_grid(v);
      }
      if (bad_range) {
        #pragma omp atomic write
        range = 1;
        continue;
      }
      g.nmiss[j] = static_cast<double>(miss);
      unsigned int hfill = 0;
      int qfill = 0;
      const double fmean = miss < n ? static_cast<double>(obs_sum) / kDosSteps / static_cast<double>(n - miss) : 0.0;
      if (miss) quantize_fill(fmean, kDosUnit, hfill, qfill);
      g.M[j].reserve(miss);
      unsigned char* buf = packed.data() + static_cast<std::size_t>(t) * sel_stride;
      std::fill(buf, buf + sel_stride, 0);
      g.D[j].unit = kDosUnit;
      DBuilder db(g, j);
      for (std::size_t i = 0; i < n; ++i) {
        const std::uint16_t v = x[rows[i]];
        unsigned int h;
        int q;
        if (v == kMiss) {
          h = hfill;
          q = qfill;
          g.M[j].push_back(static_cast<std::uint32_t>(i));
        } else {
          const unsigned int c = plink4r::dosage_grid(v);
          h = (c + kDosSteps / 2) / kDosSteps;
          q = static_cast<int>(c) - static_cast<int>(h) * kDosSteps;
        }
        buf[i >> 2] |= static_cast<unsigned char>(h << (2 * (i & 3u)));
        if (q) db.add(i, h, q);
      }
      db.finish();
      u64* lo = g.G.bits.data() + 2 * j * W;
      encode_packed(table, buf, n, sel_stride, W, lo, lo + W, g.G.s[j], g.G.v[j]);
    }
    if (failed) Rcpp::stop("Could not read PGEN file: " + path);
    if (range) Rcpp::stop("Dosage above 2 in PGEN file: " + path);
    Rcpp::checkUserInterrupt();
  }
}

// Bytes held by H (2-bit planes) and by D (masks, int16 values, prefix
// counts, per-variant headers and the three column sums).
double bytes_h(const Geno& g) { return 8.0 * g.G.bits.size(); }
double bytes_d(const Geno& g) {
  double b = (static_cast<double>(sizeof(DCol)) + 3.0 * 8.0) * g.p;
  for (const DCol& c : g.D) b += 8.0 * c.mask.size() + 2.0 * c.q.size() + 4.0 * c.pre.size();
  return b;
}

constexpr std::size_t kDRows = 2048;  // rows per block of the D accumulation (a multiple of 64)
constexpr std::size_t kDCols = 256;   // D columns per work item of the D accumulation
constexpr std::size_t kDLine = 16;    // columns decoded together (one cache line of a row of Y)
constexpr std::size_t kDSpan = 1024;  // output rows per work item of the D accumulation

// s[a] += sum_k d[k] Y[r[k] * ld + a] for a < len, the terms added in the order of k (4 rows per pass over s).
using DAxpy = void (*)(const float*, std::size_t, const std::uint32_t*, const double*, std::size_t, double*,
                       std::size_t);

__attribute__((always_inline)) inline void d_axpy_body(const float* Y, std::size_t ld, const std::uint32_t* r,
                                                       const double* d, std::size_t m, double* s, std::size_t len) {
  std::size_t k = 0;
  for (; k + 4 <= m; k += 4) {
    const float *y0 = Y + r[k] * ld, *y1 = Y + r[k + 1] * ld, *y2 = Y + r[k + 2] * ld, *y3 = Y + r[k + 3] * ld;
    const double d0 = d[k], d1 = d[k + 1], d2 = d[k + 2], d3 = d[k + 3];
    #pragma omp simd
    for (std::size_t a = 0; a < len; ++a) {
      double v = s[a];
      v += d0 * y0[a];
      v += d1 * y1[a];
      v += d2 * y2[a];
      v += d3 * y3[a];
      s[a] = v;
    }
  }
  for (; k < m; ++k) {
    const float* y = Y + r[k] * ld;
    const double dk = d[k];
    #pragma omp simd
    for (std::size_t a = 0; a < len; ++a) s[a] += dk * y[a];
  }
}

void d_axpy_generic(const float* Y, std::size_t ld, const std::uint32_t* r, const double* d, std::size_t m,
                    double* s, std::size_t len) {
  d_axpy_body(Y, ld, r, d, m, s, len);
}

#ifdef GENO_KERNEL_X86
__attribute__((target("avx2,fma")))
void d_axpy_avx2(const float* Y, std::size_t ld, const std::uint32_t* r, const double* d, std::size_t m,
                 double* s, std::size_t len) {
  d_axpy_body(Y, ld, r, d, m, s, len);
}
#endif

// Chosen as select_kernel(): the AVX2 + FMA loop whenever la() picked avx2 or avx512 (an AVX-512 loop
// measured no faster: the loop is bound by reading Y).
DAxpy select_d_axpy() {
#ifdef GENO_KERNEL_X86
  if (std::strcmp(plink4r::la().isa, "generic") != 0) return d_axpy_avx2;
#endif
  return d_axpy_generic;
}

// X'X += H'D + D'H + D'D, as S + S' with S[a, b] = sum_i D_bi (H_ai + D_ai / 2). With cols (0-based, length
// pc), the columns cols of g in that order (out is pc x pc); without, all columns of g. With per-row weights
// wrow, D_bi is multiplied by wrow_i, giving X'WX += H'WD + D'WH + D'WD. Per block of kDRows rows, Y = (H + D / 2)
// / sc of every column (sc = unit / 2) is decoded once into a row-major float buffer (row i contiguous over the
// columns); Y = 2 H / unit + q is an integer below 2^17 for a genotype column and q for an appended one (H = 0),
// so float holds it exactly. Then work items of kDSpan rows x kDCols columns of S, shared out by OpenMP
// (disjoint outputs), add D_bi Y[i, ] to S[, b] for the nonzeros of their columns in the block, in double with
// the span of S and the rows of Y it reads in cache; sc is applied at the end. Every S[a, b] sums its terms in
// ascending row order.
void add_d_xtx(const Geno& g, double* out, int threads, const int* cols = nullptr, std::size_t pc = 0,
               const double* wrow = nullptr) {
  const std::size_t n = g.n, p = cols ? pc : g.p;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  if (!g.dnnz()) return;
  const DAxpy axpy = select_d_axpy();
  const std::size_t RB = std::min(kDRows, n);
  std::vector<double> S(p * p, 0.0), sc(p), val;  // val, row: the nonzeros of the block, column by column
  std::vector<float> Y(RB * p);
  std::vector<std::uint32_t> row;
  for (std::size_t a = 0; a < p; ++a) sc[a] = 0.5 * g.D[col(a)].unit;
  std::vector<std::size_t> e0(p), off(p + 1, 0);  // per column: first nonzero of the block, start in row / val
  const std::size_t ns = (p + kDSpan - 1) / kDSpan, ng = (p + kDCols - 1) / kDCols, nc = (p + kDLine - 1) / kDLine;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> yb(kDLine * RB);
    for (std::size_t i0 = 0; i0 < n; i0 += RB) {
      const std::size_t R = std::min(RB, n - i0), w0 = i0 / 64, w1 = (i0 + R + 63) / 64;
      #pragma omp for schedule(static)
      for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(p); ++aa) {
        const std::size_t a = static_cast<std::size_t>(aa);
        const DCol& c = g.D[col(a)];
        e0[a] = c.mask.empty() ? 0 : d_rank(c, w0);
        off[a + 1] = c.mask.empty() ? 0 : d_rank(c, w1) - e0[a];
      }
      #pragma omp single
      {
        std::partial_sum(off.begin(), off.end(), off.begin());
        if (row.size() < off[p]) { row.resize(off[p]); val.resize(off[p]); }
      }
      const bool any = off[p] > 0;
      #pragma omp barrier
      if (!any) continue;  // every thread has read off[p] before the next block rewrites it
      #pragma omp for schedule(static)
      for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(nc); ++cc) {
        const std::size_t a0 = static_cast<std::size_t>(cc) * kDLine, na = std::min(kDLine, p - a0);
        for (std::size_t t = 0; t < na; ++t) {
          const std::size_t a = a0 + t, j = col(a);
          double* y = yb.data() + t * RB;
          std::fill(y, y + R, 0.0);
          const u64* lo = g.G.lo(j);
          const u64* hi = g.G.hi(j);
          for (std::size_t w = w0; w < w1; ++w) {
            for (u64 bits = lo[w]; bits; bits &= bits - 1) y[64 * (w - w0) + __builtin_ctzll(bits)] = 1.0;
            for (u64 bits = hi[w]; bits; bits &= bits - 1) y[64 * (w - w0) + __builtin_ctzll(bits)] = 2.0;
          }
          const DCol& c = g.D[j];
          const double inv = 1.0 / sc[a];
          for (std::size_t r = 0; r < R; ++r) y[r] *= inv;
          std::size_t e = e0[a], k = off[a];
          for (std::size_t w = w0; w < w1 && off[a + 1] > off[a]; ++w)
            for (u64 bits = c.mask[w]; bits; bits &= bits - 1) {
              const std::size_t r = 64 * (w - w0) + __builtin_ctzll(bits);
              const int q = c.q[e++];
              row[k] = static_cast<std::uint32_t>(r);
              val[k++] = wrow ? q * c.unit * wrow[i0 + r] : q * c.unit;
              y[r] += q;
            }
        }
        for (std::size_t r = 0; r < R; ++r)
          for (std::size_t t = 0; t < na; ++t) Y[r * p + a0 + t] = static_cast<float>(yb[t * RB + r]);
      }
      #pragma omp for schedule(dynamic, 1)
      for (std::ptrdiff_t it = 0; it < static_cast<std::ptrdiff_t>(ns * ng); ++it) {
        const std::size_t s0 = (static_cast<std::size_t>(it) % ns) * kDSpan, len = std::min(kDSpan, p - s0);
        const std::size_t b0 = (static_cast<std::size_t>(it) / ns) * kDCols, b1 = std::min(p, b0 + kDCols);
        for (std::size_t b = b0; b < b1; ++b)
          if (off[b + 1] > off[b])
            axpy(Y.data() + s0, p, row.data() + off[b], val.data() + off[b], off[b + 1] - off[b], S.data() + s0 + b * p,
                 len);
      }  // implicit barrier: Y and the lists are rewritten for the next block after every work item of this one
    }
    #pragma omp for schedule(static)
    for (std::ptrdiff_t bb = 0; bb < static_cast<std::ptrdiff_t>(p); ++bb) {
      const std::size_t b = static_cast<std::size_t>(bb);
      for (std::size_t a = 0; a < p; ++a) out[a + b * p] += sc[a] * S[a + b * p] + sc[b] * S[b + a * p];
    }
  }
}

} // namespace

// [[Rcpp::export]]
SEXP geno_open_bed_cpp(const std::string& path, int n_total, int m_total,
                       const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples,
                       bool mean, int threads) {
  threads = resolve_threads(threads);
  plink4r::BedFile bf(path, static_cast<std::size_t>(n_total), static_cast<std::size_t>(m_total), 1);
  const std::size_t stride = bf.stride();

  Rcpp::XPtr<Geno> g(new Geno, true);
  const std::vector<int> rows = to_rows(samples);
  init_geno(*g, rows.size(), variants.size(), mean);
  static const NibbleTable table({2, -1, 1, 0});  // A1 dosage: 00:2, 01:NA, 10:1, 11:0
  const std::size_t sel_stride = (rows.size() + 3) / 4, W = g->G.words;
  std::vector<unsigned char> buffer(kGroup * stride);
  std::vector<unsigned char> packed(static_cast<std::size_t>(threads) * sel_stride);
  std::vector<u64> missbuf(static_cast<std::size_t>(threads) * W);
  const std::size_t p = g->p;
  for (std::size_t g0 = 0; g0 < p; g0 += kGroup) {
    const std::size_t count = std::min(kGroup, p - g0);
    for (std::size_t k = 0; k < count; ++k) {
      const unsigned char* row = bf.read(static_cast<std::size_t>(variants[g0 + k]), 0);
      if (!row) Rcpp::stop("Could not read BED file: " + path);
      std::memcpy(buffer.data() + k * stride, row, stride);
    }
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const int t = plink4r::thread_id();
      unsigned char* buf = packed.data() + static_cast<std::size_t>(t) * sel_stride;
      gather_packed(buffer.data() + k * stride, rows, buf);
      encode_column(*g, table, buf, g0 + k, missbuf.data() + static_cast<std::size_t>(t) * W);
    }
    Rcpp::checkUserInterrupt();
  }
  finish_stats(*g);
  return g;
}

// [[Rcpp::export]]
SEXP geno_open_pgen_cpp(const std::string& path, int n_total, const Rcpp::IntegerVector& allele_ct,
                        const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples,
                        bool mean, int threads) {
  threads = resolve_threads(threads);
  PgenFile file(path, static_cast<std::size_t>(n_total), std::vector<int>(allele_ct.begin(), allele_ct.end()), threads);
  Rcpp::XPtr<Geno> g(new Geno, true);
  const std::vector<int> rows = to_rows(samples);
  if (file.dosage()) {
    for (int v : variants)
      if (allele_ct[v] > 2)
        Rcpp::stop("Multiallelic variants are not supported in PGEN files with dosages: " + path);
    init_geno(*g, rows.size(), variants.size(), true);
    load_dosage(*g, file, variants, rows, threads, path);
    finish_stats(*g);
    return g;
  }
  init_geno(*g, rows.size(), variants.size(), mean);
  const NibbleTable& table = pgen_table();
  const std::size_t sel_stride = (rows.size() + 3) / 4, W = g->G.words;
  std::vector<unsigned char> packed(static_cast<std::size_t>(threads) * sel_stride);
  std::vector<u64> missbuf(static_cast<std::size_t>(threads) * W);
  const std::size_t p = g->p;
  for (std::size_t g0 = 0; g0 < p; g0 += kGroup) {
    const std::size_t count = std::min(kGroup, p - g0);
    int failed = 0;
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const int t = plink4r::thread_id();
      const unsigned char* column = file.read(static_cast<std::size_t>(variants[g0 + k]), t);
      if (!column) {
        #pragma omp atomic write
        failed = 1;
        continue;
      }
      unsigned char* buf = packed.data() + static_cast<std::size_t>(t) * sel_stride;
      gather_packed(column, rows, buf);
      encode_column(*g, table, buf, g0 + k, missbuf.data() + static_cast<std::size_t>(t) * W);
    }
    if (failed) Rcpp::stop("Could not read PGEN file: " + path);
    Rcpp::checkUserInterrupt();
  }
  finish_stats(*g);
  return g;
}

namespace {

// Sum of the stored values of column j over all cells, in int16 steps of a genotype column's unit (k per allele).
inline std::int64_t step_sum_col(const Geno& g, std::size_t j, int k) {
  std::int64_t tot = k * static_cast<std::int64_t>(g.G.s[j]);
  for (std::int16_t q : g.D[j].q) tot += q;
  return tot;
}

// Sum of the stored values of column j over its observed (not missing) cells. Genotype columns are summed
// exactly in int16 steps, so the result equals the sum of the file's observed calls (hard calls) or 0.01
// dosages; an appended column (no missing cells) gives its column sum.
double obs_sum_col(const Geno& g, std::size_t j) {
  const DCol& c = g.D[j];
  const int k = unit_steps(c.unit);
  if (!k) return static_cast<double>(g.G.s[j]) + g.dsum[j];
  std::int64_t tot = step_sum_col(g, j, k);
  for (std::uint32_t i : g.M[j]) tot -= k * static_cast<std::int64_t>(h_at(g.G, j, i)) + q_at(c, i);
  return static_cast<double>(tot) / k;
}

// Rebuilds D of column j from (row, q) entries in ascending row order and refreshes the column's sums,
// plane statistics, mean and sd. The planes must already hold the final H.
template <class Entries>
void rebuild_col(Geno& g, std::size_t j, double unit, Entries&& each) {
  g.D[j] = DCol();
  g.D[j].unit = unit;
  DBuilder db(g, j);
  each([&](std::size_t i, int q) { if (q) db.add(i, h_at(g.G, j, i), q); });
  db.finish();
  plane_stats(g, j);
  col_stats(g, j);
}

// Exact flip of one stored cell x = h + q / k -> 2 - x (k int16 steps per allele), in the canonical form of a
// direct load (q = +k/2 becomes h + 1, q = -k/2).
inline void flip_cell(unsigned int& h, int& q, int k) {
  h = 2u - h;
  q = -q;
  if (q == k / 2 && h < 2u) { ++h; q = -(k / 2); }
}

}  // namespace

namespace {

// The kernels of the products, writing into caller buffers (column-major); the R functions below and the
// C-callable entry points (bottom of this file) both call them.

// Rows r0..r0+R-1 of stored column j (H + D) as T (float or double) into dst[0..R). grid: the value in steps of
// 0.01 minus 100 (an integer in -100..100), for a column on the 0.01 grid (on_grid).
template <class T>
void decode_rows(const Geno& g, std::size_t j, std::size_t r0, std::size_t R, T* dst, bool grid = false) {
  const T v0 = grid ? T(-kDosSteps) : T(0), v1 = grid ? T(0) : T(1), v2 = grid ? T(kDosSteps) : T(2);
  std::fill(dst, dst + R, v0);
  const u64* lo = g.G.lo(j);
  const u64* hi = g.G.hi(j);
  const std::size_t w0 = r0 / 64, w1 = (r0 + R + 63) / 64;
  for (std::size_t x = w0; x < w1; ++x) {
    for (u64 bits = lo[x]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * x + __builtin_ctzll(bits);
      if (i >= r0 && i < r0 + R) dst[i - r0] = v1;
    }
    for (u64 bits = hi[x]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * x + __builtin_ctzll(bits);
      if (i >= r0 && i < r0 + R) dst[i - r0] = v2;
    }
  }
  const DCol& dc = g.D[j];
  if (!dc.mask.empty()) {
    std::size_t e = d_rank(dc, w0);
    for (std::size_t x = w0; x < w1 && x < dc.mask.size(); ++x)
      for (u64 bits = dc.mask[x]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * x + __builtin_ctzll(bits);
        const double d = grid ? dc.q[e++] : dc.q[e++] * dc.unit;
        if (i >= r0 && i < r0 + R) dst[i - r0] += static_cast<T>(d);
      }
  }
}

// The dense product loop of xtwx_float_into and xtx_dense_into, into o (P x P; P = pc for columns cols, else
// g.p). Per block of RB rows, every column is decoded once into a shared buffer of T (column a at X[a * LD],
// LD = RB + pad: no power-of-two stride, which thrashes the cache); then the output tiles of TS x TS of the
// lower triangle are shared out by OpenMP (disjoint outputs), each adding its block product to o through
// tile(X, LD, r0, R, a0, na, b0, nb) of the thread's own worker from make(). The lower triangle is mirrored.
// grid: the columns are decoded as centred 0.01 steps (decode_rows).
template <class T, class Make>
void dense_tiles(const Geno& g, const int* cols, std::size_t pc, std::size_t RB, std::size_t pad, std::size_t TS,
                 double* o, int threads, Make make, bool grid = false) {
  const std::size_t n = g.n, P = cols ? pc : g.p, LD = RB + pad;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  const std::size_t nt = (P + TS - 1) / TS;
  std::vector<std::pair<std::size_t, std::size_t>> tiles;
  for (std::size_t b = 0; b < nt; ++b) for (std::size_t a = b; a < nt; ++a) tiles.emplace_back(a, b);
  std::vector<T> X(LD * P);
  std::fill(o, o + P * P, 0.0);
  #pragma omp parallel num_threads(threads)
  {
    auto tile = make();
    for (std::size_t r0 = 0; r0 < n; r0 += RB) {
      const std::size_t R = std::min(RB, n - r0);
      #pragma omp for schedule(static)
      for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(P); ++aa)
        decode_rows(g, col(static_cast<std::size_t>(aa)), r0, R, &X[static_cast<std::size_t>(aa) * LD], grid);
      #pragma omp for schedule(dynamic, 1)
      for (std::ptrdiff_t tt = 0; tt < static_cast<std::ptrdiff_t>(tiles.size()); ++tt) {
        const std::size_t a0 = tiles[tt].first * TS, b0 = tiles[tt].second * TS;
        tile(X.data(), LD, r0, R, a0, std::min(TS, P - a0), b0, std::min(TS, P - b0));
      }  // implicit barrier: the buffer is rewritten for the next block after every tile of this one
    }
  }
  for (std::size_t b = 0; b < P; ++b)  // mirror the lower triangle
    for (std::size_t a = b + 1; a < P; ++a) o[b + a * P] = o[a + b * P];
}

constexpr std::size_t kFTile = 256;    // columns per side of an output tile of the float kernel
constexpr std::size_t kFRows = 4096;   // rows per block of the float kernel
constexpr std::size_t kFPad = 16;      // floats after each column of the shared buffer (no power-of-two stride)

// X'WX into o (p x p, or pc x pc for columns cols), W = diag(w), any real w, in float products with double
// sums (dense_tiles). Blocks of kFRows rows, 16 KB per column of the shared buffer (66 MB at p = 4000, 329 MB
// at p = 20000); a tile copies its row group transposed and times w, multiplies it with its column group in
// the shared buffer by the run-time-dispatched serial sgemm and accumulates in double.
void xtwx_float_into(const Geno& g, const double* w, const int* cols, std::size_t pc, double* o, int threads) {
  const std::size_t P = cols ? pc : g.p, RB = std::min(kFRows, g.n);
  const plink4r::LaFns& LA = plink4r::la();
  dense_tiles<float>(g, cols, pc, RB, kFPad, kFTile, o, threads, [&] {
    return [&, At = std::vector<float>(kFTile * RB), Pt = std::vector<float>(kFTile * kFTile)](
               const float* X, std::size_t LD, std::size_t r0, std::size_t R, std::size_t a0, std::size_t na,
               std::size_t b0, std::size_t nb) mutable {
      const float* xa = &X[a0 * LD];
      for (std::size_t r = 0; r < R; ++r) {  // A' : na x R, ld kFTile
        const double wr = w[r0 + r];
        for (std::size_t c = 0; c < na; ++c) At[c + r * kFTile] = static_cast<float>(xa[c * LD + r] * wr);
      }
      LA.sgemm_nn(At.data(), static_cast<int>(na), static_cast<int>(R), static_cast<int>(kFTile), &X[b0 * LD],
                  static_cast<int>(nb), static_cast<int>(LD), Pt.data(), static_cast<int>(kFTile));
      for (std::size_t y = 0; y < nb; ++y)
        for (std::size_t x = 0; x < na; ++x) o[(a0 + x) + (b0 + y) * P] += Pt[x + y * kFTile];
    };
  });
}

constexpr std::size_t kXTile = 256;    // columns per side of an output tile of the double X'X kernel
constexpr std::size_t kXRows = 2048;   // rows per block of the double X'X kernel
constexpr std::size_t kXPad = 8;       // doubles after each column of the shared buffer (no power-of-two stride)

// X'X into o (p x p, or pc x pc for columns cols) in double (dense_tiles): blocks of kXRows rows, 16 KB per
// column of the shared buffer (66 MB at p = 4000, 329 MB at p = 20000); a tile adds X_a'X_b of its row and
// column groups to o by the run-time-dispatched serial dgemm. Columns in 1/16384 steps hold multiples of 2^-14
// in [0, 2], so for n below 2^23 every product and partial sum is exact and the result equals the popcount one.
void xtx_dense_into(const Geno& g, const int* cols, std::size_t pc, double* o, int threads) {
  const std::size_t P = cols ? pc : g.p, RB = std::min(kXRows, g.n);
  const plink4r::LaFns& LA = plink4r::la();
  dense_tiles<double>(g, cols, pc, RB, kXPad, kXTile, o, threads, [&] {
    return [&](const double* X, std::size_t LD, std::size_t, std::size_t R, std::size_t a0, std::size_t na,
               std::size_t b0, std::size_t nb) {
      LA.dgemm_tn(&X[a0 * LD], static_cast<int>(R), static_cast<int>(na), static_cast<int>(LD), &X[b0 * LD],
                  static_cast<int>(nb), static_cast<int>(LD), &o[a0 + b0 * P], static_cast<int>(P), 1);
    };
  });
}

constexpr double kDDenseFrac = 0.01;     // D nonzero fraction above which xtwx_bin_into uses the float kernel
constexpr std::size_t kGRows = 1024;   // rows per block of the grid X'X kernel: partial sums <= 1024 * 100^2 < 2^24

// Whether every column of cols (or all) is on the 0.01 grid: no D, or D in steps of 0.01 (dosage PGEN).
bool on_grid(const Geno& g, const int* cols, std::size_t pc) {
  const std::size_t P = cols ? pc : g.p;
  for (std::size_t a = 0; a < P; ++a) {
    const DCol& c = g.D[cols ? cols[a] : a];
    if (!c.mask.empty() && c.unit != kDosUnit) return false;
  }
  return true;
}

// X'X into o (p x p, or pc x pc for columns cols) for columns on the 0.01 grid (on_grid), exact, in float
// products (dense_tiles, tiles and buffer as xtwx_float_into, blocks of kGRows rows). A column is decoded as
// the integer c = 100 x - 100 in -100..100; within a block every product and partial sum is an integer below
// 2^24, so the serial sgemm is exact in any order; the blocks are summed in double (integers below 2^53) and
// X'X = (C'C + 100 (s_a + s_b) - 1e4 n) / 1e4, s the column sums of 100 x.
void xtx_grid_into(const Geno& g, const int* cols, std::size_t pc, double* o, int threads) {
  const std::size_t P = cols ? pc : g.p, RB = std::min(kGRows, g.n);
  const plink4r::LaFns& LA = plink4r::la();
  dense_tiles<float>(g, cols, pc, RB, kFPad, kFTile, o, threads, [&] {
    return [&, At = std::vector<float>(kFTile * RB), Pt = std::vector<float>(kFTile * kFTile)](
               const float* X, std::size_t LD, std::size_t, std::size_t R, std::size_t a0, std::size_t na,
               std::size_t b0, std::size_t nb) mutable {
      const float* xa = &X[a0 * LD];
      for (std::size_t r = 0; r < R; ++r)  // A' : na x R, ld kFTile
        for (std::size_t c = 0; c < na; ++c) At[c + r * kFTile] = xa[c * LD + r];
      LA.sgemm_nn(At.data(), static_cast<int>(na), static_cast<int>(R), static_cast<int>(kFTile), &X[b0 * LD],
                  static_cast<int>(nb), static_cast<int>(LD), Pt.data(), static_cast<int>(kFTile));
      for (std::size_t y = 0; y < nb; ++y)
        for (std::size_t x = 0; x < na; ++x) o[(a0 + x) + (b0 + y) * P] += Pt[x + y * kFTile];
    };
  }, true);
  std::vector<double> s(P);
  for (std::size_t a = 0; a < P; ++a)
    s[a] = static_cast<double>(step_sum_col(g, cols ? static_cast<std::size_t>(cols[a]) : a, kDosSteps));
  const double k = kDosSteps, nk2 = k * k * static_cast<double>(g.n);
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t b = 0; b < static_cast<std::ptrdiff_t>(P); ++b)
    for (std::size_t a = 0; a < P; ++a) {
      double& v = o[a + static_cast<std::size_t>(b) * P];
      v = (v + k * (s[a] + s[static_cast<std::size_t>(b)]) - nk2) / (k * k);
    }
}

constexpr double kDDenseFracXtx = 0.07;  // D nonzero fraction above which xtx_into uses a dense kernel

// Whether more than frac of the cells of the columns cols (or all) have a nonzero D.
bool d_dense(const Geno& g, const int* cols, std::size_t pc, double frac) {
  const std::size_t P = cols ? pc : g.p;
  double nnz = 0.0;
  for (std::size_t a = 0; a < P; ++a) nnz += static_cast<double>(g.D[cols ? cols[a] : a].q.size());
  return nnz > frac * static_cast<double>(g.n) * static_cast<double>(P);
}

// X'X of the stored matrix into o (p x p), or of columns cols (0-based, length pc) into o (pc x pc): exact
// popcounts of H'H plus add_d_xtx, or, when more than kDDenseFracXtx of the cells of the columns have a nonzero
// D (dosages), a dense kernel: xtx_grid_into for columns on the 0.01 grid, else the double xtx_dense_into.
void xtx_into(const Geno& g, const int* cols, std::size_t pc, double* o, int threads) {
  if (d_dense(g, cols, pc, kDDenseFracXtx)) {
    if (on_grid(g, cols, pc)) xtx_grid_into(g, cols, pc, o, threads);
    else xtx_dense_into(g, cols, pc, o, threads);
    return;
  }
  const TileKernel kernel = select_kernel();
  if (!cols) {
    const std::size_t p = g.p;
    dot_block(g.G, 0, p, g.G, 0, p, true, o, p, threads, kernel);
    add_d_xtx(g, o, threads);
    return;
  }
  Planes sub;
  sub.words = g.G.words;
  sub.bits.resize(2 * pc * sub.words);
  for (std::size_t a = 0; a < pc; ++a) {
    const std::size_t j = static_cast<std::size_t>(cols[a]);
    std::copy(g.G.lo(j), g.G.lo(j) + sub.words, sub.bits.begin() + 2 * a * sub.words);
    std::copy(g.G.hi(j), g.G.hi(j) + sub.words, sub.bits.begin() + (2 * a + 1) * sub.words);
  }
  dot_block(sub, 0, pc, sub, 0, pc, true, o, pc, threads, kernel);
  add_d_xtx(g, o, threads, cols, pc);
}

constexpr std::size_t kBinRows = 64;     // rows per weight bin: one 64-bit word of the gathered planes
constexpr std::size_t kBinGroup = 384;   // bins per streamed group (24576 rows; a multiple of 4)
// X'W~X into o (p x p, or pc x pc for columns cols), W~ = diag(w binned): rows sorted by w (ties by row
// index), each run of kBinRows sorted rows takes the mean w of its rows. With the weight constant within a
// bin, H'W~H = sum_b wbar_b H_b'H_b by exact popcounts: per group of kBinGroup bins the rows of the group
// are gathered for every column (one bin per word), then the 64 x 64 tiles of the upper triangle are shared
// out by OpenMP (disjoint outputs). The D terms come from add_d_xtx with the binned per-row weights. When
// more than kDDenseFrac of the cells of the columns have a nonzero D (dosages), the popcounts do not pay:
// X'WX with w itself comes from xtwx_float_into, and W~ = W. w_bin (n, row order, may be null) receives the
// weights of the result.
void xtwx_bin_into(const Geno& g, const double* w, const int* cols, std::size_t pc, double* o, double* w_bin,
                   int threads) {
  const std::size_t n = g.n, P = cols ? pc : g.p;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  if (d_dense(g, cols, pc, kDDenseFrac)) {
    xtwx_float_into(g, w, cols, pc, o, threads);
    if (w_bin) std::copy(w, w + n, w_bin);
    return;
  }
  std::vector<std::uint32_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0u);
  std::stable_sort(ord.begin(), ord.end(), [w](std::uint32_t a, std::uint32_t b) { return w[a] < w[b]; });
  const std::size_t nbin = (n + kBinRows - 1) / kBinRows, nbin4 = (nbin + 3) / 4 * 4;
  std::vector<double> wbar(nbin4, 0.0), wq(n);  // padding bins (no rows) weigh 0
  for (std::size_t b = 0; b < nbin; ++b) {
    const std::size_t r0 = b * kBinRows, r1 = std::min(n, r0 + kBinRows);
    double s = 0.0;
    for (std::size_t r = r0; r < r1; ++r) s += w[ord[r]];
    wbar[b] = s / static_cast<double>(r1 - r0);
    for (std::size_t r = r0; r < r1; ++r) wq[ord[r]] = wbar[b];
  }
  if (w_bin) std::copy(wq.begin(), wq.end(), w_bin);
  std::fill(o, o + P * P, 0.0);
  const WTileKernel kernel = select_wkernel();
  const std::size_t nt = (P + kTile - 1) / kTile, gw = std::min(kBinGroup, nbin4);
  std::vector<u64> store(2 * P * gw + 3);  // 3 spare words to start the buffer on 32 bytes
  void* start = store.data();
  std::size_t space = store.size() * sizeof(u64);
  u64* buf = static_cast<u64*>(std::align(32, 2 * P * gw * sizeof(u64), start, space));
  for (std::size_t b0 = 0; b0 < nbin4; b0 += kBinGroup) {
    const std::size_t W = std::min(kBinGroup, nbin4 - b0);  // words per plane of this group
    const std::size_t g0 = b0 * kBinRows, g1 = std::min(n, (b0 + W) * kBinRows);
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(P); ++aa) {
      const std::size_t a = static_cast<std::size_t>(aa), j = col(a);
      const u64* lo = g.G.lo(j);
      const u64* hi = g.G.hi(j);
      u64* ol = buf + 2 * a * W;
      u64* oh = ol + W;
      for (std::size_t k = 0; k < W; ++k) {
        u64 l = 0, h = 0;
        const std::size_t r0 = g0 + kBinRows * k, r1 = std::min(g1, r0 + kBinRows);
        for (std::size_t r = r0; r < r1; ++r) {
          const std::uint32_t i = ord[r];
          l |= ((lo[i >> 6] >> (i & 63u)) & 1u) << (r - r0);
          h |= ((hi[i >> 6] >> (i & 63u)) & 1u) << (r - r0);
        }
        ol[k] = l;
        oh[k] = h;
      }
    }
    #pragma omp parallel num_threads(threads)
    {
      std::vector<double> acc(kTile * kTile);
      #pragma omp for schedule(dynamic, 1)
      for (std::ptrdiff_t t = 0; t < static_cast<std::ptrdiff_t>(nt * nt); ++t) {
        const std::size_t ti = static_cast<std::size_t>(t) / nt, tj = static_cast<std::size_t>(t) % nt;
        if (ti > tj) continue;
        const std::size_t i0 = ti * kTile, j0 = tj * kTile, ni = std::min(kTile, P - i0), nj = std::min(kTile, P - j0);
        std::fill(acc.begin(), acc.end(), 0.0);
        kernel(buf, W, i0, ni, j0, nj, wbar.data() + b0, acc.data());
        for (std::size_t x = 0; x < ni; ++x)
          for (std::size_t y = 0; y < nj; ++y) o[(i0 + x) + (j0 + y) * P] += acc[x * kTile + y];
      }
    }
  }
  for (std::size_t b = 0; b < P; ++b)  // mirror the upper triangle
    for (std::size_t a = b + 1; a < P; ++a) o[a + b * P] = o[b + a * P];
  add_d_xtx(g, o, threads, cols, pc, wq.data());
}

// Standardizes a weighted cross product in place: o (pc x pc, columns cols or all when null) becomes
// (X'WX - c s' - s c' + (sum w) c c') / (d d'), s = X'w in double from the planes and D; center c and div d
// are per output column (div may be null: no division).
void xtwx_std(const Geno& g, const double* w, const int* cols, std::size_t pc, const double* center,
              const double* div, double* o, int threads) {
  const std::size_t n = g.n, P = cols ? pc : g.p;
  auto col = [cols](std::size_t a) { return cols ? static_cast<std::size_t>(cols[a]) : a; };
  std::vector<double> sw(P, 0.0);
  double wsum = 0.0;
  for (std::size_t i = 0; i < n; ++i) wsum += w[i];
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t aa = 0; aa < static_cast<std::ptrdiff_t>(P); ++aa) {
    const std::size_t a = static_cast<std::size_t>(aa), j = col(a);
    double acc = 0.0;
    const u64* lo = g.G.lo(j);
    const u64* hi = g.G.hi(j);
    for (std::size_t x = 0; x < g.G.words; ++x) {
      for (u64 bits = lo[x]; bits; bits &= bits - 1) acc += w[64 * x + __builtin_ctzll(bits)];
      for (u64 bits = hi[x]; bits; bits &= bits - 1) acc += 2.0 * w[64 * x + __builtin_ctzll(bits)];
    }
    const DCol& dc = g.D[j];
    std::size_t e = 0;
    for (std::size_t x = 0; x < dc.mask.size(); ++x)
      for (u64 bits = dc.mask[x]; bits; bits &= bits - 1) acc += dc.q[e++] * dc.unit * w[64 * x + __builtin_ctzll(bits)];
    sw[a] = acc;
  }
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t bb = 0; bb < static_cast<std::ptrdiff_t>(P); ++bb) {
    const std::size_t b = static_cast<std::size_t>(bb);
    for (std::size_t a = 0; a < P; ++a) {
      double v = o[a + b * P] - center[a] * sw[b] - sw[a] * center[b] + wsum * center[a] * center[b];
      if (div) v /= div[a] * div[b];
      o[a + b * P] = v;
    }
  }
}

// X'M into o (p x k) for a dense n x k matrix m.
void xtm_into(const Geno& g_, const double* m, std::size_t k, double* o, int threads) {
  const Geno* g = &g_;
  const std::size_t n = g->n, p = g->p, W = g->G.words;
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(p); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    std::vector<double> s1(k, 0.0), s2(k, 0.0), sd(k, 0.0);
    const u64* lo = g->G.lo(j);
    const u64* hi = g->G.hi(j);
    for (std::size_t w = 0; w < W; ++w) {
      for (u64 bits = lo[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) s1[c] += m[i + c * n];
      }
      for (u64 bits = hi[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) s2[c] += m[i + c * n];
      }
    }
    const DCol& dc = g->D[j];
    std::size_t e = 0;
    for (std::size_t w = 0; w < dc.mask.size(); ++w)
      for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) {
        const double d = dc.q[e++] * dc.unit;
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        for (std::size_t c = 0; c < k; ++c) sd[c] += d * m[i + c * n];
      }
    for (std::size_t c = 0; c < k; ++c) o[j + c * p] = s1[c] + 2.0 * s2[c] + sd[c];
  }
}

// X B into o (n x k, overwritten) for a dense p x k matrix b; all-zero rows of b are skipped.
void mx_into(const Geno& g_, const double* b, std::size_t k, double* o, int threads) {
  const Geno* g = &g_;
  const std::size_t n = g->n, p = g->p;
  std::fill(o, o + n * k, 0.0);
  std::vector<char> nonzero(p, 0);
  for (std::size_t j = 0; j < p; ++j)
    for (std::size_t c = 0; c < k; ++c)
      if (b[j + c * p] != 0.0) { nonzero[j] = 1; break; }
  const std::size_t used_words = (n + 63) / 64;
  #pragma omp parallel num_threads(threads)
  {
    const int nt = plink4r::thread_team(), t = plink4r::thread_id();
    const std::size_t w0 = used_words * t / nt, w1 = used_words * (t + 1) / nt;
    std::vector<double> bj(k);
    for (std::size_t j = 0; j < p; ++j) {
      if (!nonzero[j]) continue;
      for (std::size_t c = 0; c < k; ++c) bj[c] = b[j + c * p];
      const u64* lo = g->G.lo(j);
      const u64* hi = g->G.hi(j);
      for (std::size_t w = w0; w < w1; ++w) {
        for (u64 bits = lo[w]; bits; bits &= bits - 1) {
          const std::size_t i = 64 * w + __builtin_ctzll(bits);
          for (std::size_t c = 0; c < k; ++c) o[i + c * n] += bj[c];
        }
        for (u64 bits = hi[w]; bits; bits &= bits - 1) {
          const std::size_t i = 64 * w + __builtin_ctzll(bits);
          for (std::size_t c = 0; c < k; ++c) o[i + c * n] += 2.0 * bj[c];
        }
      }
      const DCol& dc = g->D[j];
      if (!dc.mask.empty() && w0 < w1) {
        std::size_t e = d_rank(dc, w0);
        for (std::size_t w = w0; w < w1; ++w)
          for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) {
            const double d = dc.q[e++] * dc.unit;
            const std::size_t i = 64 * w + __builtin_ctzll(bits);
            for (std::size_t c = 0; c < k; ++c) o[i + c * n] += d * bj[c];
          }
      }
    }
  }
}

// Dense block X[rows, cols] of the stored matrix (0-based) into o (nr x nc).
void dense_into(const Geno& g_, const int* rows, std::size_t nr, const int* cols, std::size_t nc, double* o,
                int threads) {
  const Geno* g = &g_;
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(nc); ++cc) {
    const std::size_t j = static_cast<std::size_t>(cols[cc]);
    const u64* lo = g->G.lo(j);
    const u64* hi = g->G.hi(j);
    double* col = o + static_cast<std::size_t>(cc) * nr;
    for (std::size_t r = 0; r < nr; ++r) {
      const std::size_t i = static_cast<std::size_t>(rows[r]);
      const std::size_t w = i / 64;
      const u64 bit = u64(1) << (i % 64);
      col[r] = (lo[w] & bit) ? 1.0 : (hi[w] & bit) ? 2.0 : 0.0;
    }
    const DCol& dc = g->D[j];
    if (!dc.mask.empty())
      for (std::size_t r = 0; r < nr; ++r) {
        const double d = d_at(dc, static_cast<std::size_t>(rows[r]));
        if (d != 0.0) col[r] += d;
      }
  }
}

// o = (k - s u v') / (da db'), elementwise (r x c); a null da or db is 1. o may be k (in place).
void std_into(const double* k, std::size_t r, std::size_t c, const double* u, const double* v, double s,
              const double* da, const double* db, double* o, int threads) {
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(c); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    const double vj = s * v[j], bj = db ? db[j] : 1.0;
    for (std::size_t i = 0; i < r; ++i) o[i + j * r] = (k[i + j * r] - u[i] * vj) / ((da ? da[i] : 1.0) * bj);
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List geno_info_cpp(SEXP ptr) {
  Rcpp::XPtr<Geno> g(ptr);
  std::vector<double> os(g->p);
  for (std::size_t j = 0; j < g->p; ++j) os[j] = obs_sum_col(*g, j);
  return Rcpp::List::create(Rcpp::Named("n") = static_cast<double>(g->n),
                            Rcpp::Named("p") = static_cast<double>(g->p),
                            Rcpp::Named("mean") = g->mu,
                            Rcpp::Named("sd") = g->sd,
                            Rcpp::Named("n_miss") = g->nmiss,
                            Rcpp::Named("obs_sum") = os,
                            Rcpp::Named("dosage") = g->dosage,
                            Rcpp::Named("dosage_nnz") = static_cast<double>(g->dnnz()),
                            Rcpp::Named("bytes_H") = bytes_h(*g),
                            Rcpp::Named("bytes_D") = bytes_d(*g));
}

// X'X of the stored, unstandardized matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtx_cpp(SEXP ptr, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t p = g->p;
  Rcpp::NumericMatrix out(p, p);
  xtx_into(*g, nullptr, 0, out.begin(), threads);
  return out;
}

// X'M for a dense n x k matrix M.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtm_cpp(SEXP ptr, const Rcpp::NumericMatrix& M, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  const std::size_t k = M.ncol();
  if (static_cast<std::size_t>(M.nrow()) != n) Rcpp::stop("M must have one row per sample.");
  Rcpp::NumericMatrix out(p, k);
  xtm_into(*g, M.begin(), k, out.begin(), threads);
  return out;
}

// X B for a dense p x k matrix B; all-zero rows of B are skipped.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_mx_cpp(SEXP ptr, const Rcpp::NumericMatrix& B, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  const std::size_t k = B.ncol();
  if (static_cast<std::size_t>(B.nrow()) != p) Rcpp::stop("B must have one row per variant.");
  Rcpp::NumericMatrix out(n, k);
  mx_into(*g, B.begin(), k, out.begin(), threads);
  return out;
}

// A new geno = the columns of ptr followed by the k columns of M (n x k, any
// finite values: dense or sparse, 0/1/2 or not). An appended column has H = 0
// and all of x in D, as int16 q = round(x / unit) with its own unit
// max|x| / 32767 (relative error <= 1/65534 of max|x|); only its nonzeros are
// stored. Its mean and sd are kept like those of every column. ptr is unchanged.
// [[Rcpp::export]]
SEXP geno_append_cpp(SEXP ptr, const Rcpp::NumericMatrix& M, int threads) {
  Rcpp::XPtr<Geno> a(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = a->n, p0 = a->p, k = M.ncol();
  if (static_cast<std::size_t>(M.nrow()) != n) Rcpp::stop("M must have one row per sample.");
  if (!k) Rcpp::stop("M has no columns.");
  const double* m = M.begin();
  for (std::size_t e = 0; e < n * k; ++e)
    if (!std::isfinite(m[e])) Rcpp::stop("Appended values must be finite.");
  Rcpp::XPtr<Geno> g(new Geno, true);
  init_geno(*g, n, p0 + k, a->mean);
  g->dosage = a->dosage;
  std::copy(a->G.bits.begin(), a->G.bits.end(), g->G.bits.begin());
  std::copy(a->G.s.begin(), a->G.s.end(), g->G.s.begin());
  std::copy(a->G.v.begin(), a->G.v.end(), g->G.v.begin());
  for (std::size_t j = 0; j < p0; ++j) {
    g->D[j] = a->D[j];
    g->dsum[j] = a->dsum[j];
    g->dsq[j] = a->dsq[j];
    g->dhd[j] = a->dhd[j];
    g->nmiss[j] = a->nmiss[j];
    g->M[j] = a->M[j];
  }
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(k); ++cc) {
    const std::size_t c = static_cast<std::size_t>(cc), j = p0 + c;
    const double* x = m + c * n;
    double amax = 0.0;
    for (std::size_t i = 0; i < n; ++i) amax = std::max(amax, std::fabs(x[i]));
    g->D[j].unit = amax > 0.0 ? amax / 32767.0 : kDUnit;
    const double inv = 1.0 / g->D[j].unit;
    DBuilder db(*g, j);
    for (std::size_t i = 0; i < n; ++i) {
      const int q = static_cast<int>(std::lround(x[i] * inv));
      if (q) db.add(i, 0u, q);
    }
    db.finish();
    plane_stats(*g, j);
  }
  finish_stats(*g);
  return g;
}

// Raw A'B of two geno objects over the same samples, both without D (hard
// calls, median fill): exact popcount on the planes.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_cross_cpp(SEXP ptr_a, SEXP ptr_b, int threads) {
  Rcpp::XPtr<Geno> a(ptr_a), b(ptr_b);
  threads = resolve_threads(threads);
  if (a->n != b->n) Rcpp::stop("A and B must have the same samples.");
  if (a->dnnz() || b->dnnz()) Rcpp::stop("geno_cross_cpp needs hard-call genotypes without D.");
  const TileKernel kernel = select_kernel();
  Rcpp::NumericMatrix out(a->p, b->p);
  dot_block(a->G, 0, a->p, b->G, 0, b->p, false, out.begin(), a->p, threads, kernel);
  return out;
}

// Dense block X[rows, cols] of the stored, unstandardized matrix (0-based).
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_dense_cpp(SEXP ptr, const Rcpp::IntegerVector& rows,
                                   const Rcpp::IntegerVector& cols, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t nr = rows.size(), nc = cols.size();
  Rcpp::NumericMatrix out(nr, nc);
  dense_into(*g, rows.begin(), nr, cols.begin(), nc, out.begin(), threads);
  return out;
}

// Standardized products in C++ (no BLAS from R): out = (K - s u v') / (da db'), elementwise; an empty
// da or db is 1. geno_crossprod / geno_cor give their centered, scaled forms through this.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_std_cpp(const Rcpp::NumericMatrix& K, const Rcpp::NumericVector& u,
                                 const Rcpp::NumericVector& v, double s, const Rcpp::NumericVector& da,
                                 const Rcpp::NumericVector& db, int threads) {
  threads = resolve_threads(threads);
  const std::size_t r = K.nrow(), c = K.ncol();
  if (static_cast<std::size_t>(u.size()) != r || static_cast<std::size_t>(v.size()) != c) Rcpp::stop("u, v do not match K.");
  const bool ha = da.size() > 0, hb = db.size() > 0;
  if ((ha && static_cast<std::size_t>(da.size()) != r) || (hb && static_cast<std::size_t>(db.size()) != c)) Rcpp::stop("da, db do not match K.");
  Rcpp::NumericMatrix out(r, c);
  std_into(K.begin(), r, c, u.begin(), v.begin(), s, ha ? da.begin() : nullptr, hb ? db.begin() : nullptr, out.begin(),
           threads);
  return out;
}

// A'A of a dense n x k matrix (double, columns pairs shared out by OpenMP).
// [[Rcpp::export]]
Rcpp::NumericMatrix dense_xtx_cpp(const Rcpp::NumericMatrix& A, int threads) {
  threads = resolve_threads(threads);
  const std::size_t n = A.nrow(), k = A.ncol();
  Rcpp::NumericMatrix out(k, k);
  const double* a = A.begin();
  double* o = out.begin();
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(k); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    for (std::size_t i = j; i < k; ++i) {
      double acc = 0.0;
      for (std::size_t r = 0; r < n; ++r) acc += a[r + i * n] * a[r + j * n];
      o[i + j * k] = acc; o[j + i * k] = acc;
    }
  }
  return out;
}

// X'WX of the stored matrix, W = diag(w), any real w (xtwx_float_into). With center (and div) given, returns
// the product of the standardized matrix: (X'WX - c s' - s c' + (sum w) c c') / (d d'), s = X'w.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtwx_cpp(SEXP ptr, const Rcpp::NumericVector& w, const Rcpp::NumericVector& center,
                                  const Rcpp::NumericVector& div, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  if (static_cast<std::size_t>(w.size()) != n) Rcpp::stop("w must have one value per sample.");
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(w[i])) Rcpp::stop("w must be finite.");
  if (center.size() > 0 && static_cast<std::size_t>(center.size()) != p) Rcpp::stop("center must have one value per column.");
  Rcpp::NumericMatrix out(p, p);
  xtwx_float_into(*g, w.begin(), nullptr, p, out.begin(), threads);
  if (center.size() > 0)
    xtwx_std(*g, w.begin(), nullptr, p, center.begin(), div.size() > 0 ? div.begin() : nullptr, out.begin(), threads);
  return out;
}

// X'W~X of the stored matrix with w binned in groups of 64 rows, or with w itself when more than 1% of the
// cells have a nonzero D (xtwx_bin_into), standardized as in geno_xtwx_cpp with those weights when center is
// given; w = the weights of the result in row order.
// [[Rcpp::export]]
Rcpp::List geno_xtwx_bin_cpp(SEXP ptr, const Rcpp::NumericVector& w, const Rcpp::NumericVector& center,
                             const Rcpp::NumericVector& div, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  if (static_cast<std::size_t>(w.size()) != n) Rcpp::stop("w must have one value per sample.");
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(w[i])) Rcpp::stop("w must be finite.");
  if (center.size() > 0 && static_cast<std::size_t>(center.size()) != p) Rcpp::stop("center must have one value per column.");
  Rcpp::NumericMatrix out(p, p);
  Rcpp::NumericVector wq(n);
  xtwx_bin_into(*g, w.begin(), nullptr, p, out.begin(), wq.begin(), threads);
  if (center.size() > 0)
    xtwx_std(*g, wq.begin(), nullptr, p, center.begin(), div.size() > 0 ? div.begin() : nullptr, out.begin(), threads);
  return Rcpp::List::create(Rcpp::Named("XtWX") = out, Rcpp::Named("w") = wq);
}

// A new geno = rows x cols of ptr (0-based, any order; repeats allowed), each cell copied as stored (no
// refill); flip[c] (optional) stores 2 - x in output column c, exactly on (H, D): h -> 2 - h, q -> -q
// (q = +8192 -> h + 1, -8192). Missing-cell rows follow the selected rows. ptr is unchanged.
// [[Rcpp::export]]
SEXP geno_take_cpp(SEXP ptr, const Rcpp::IntegerVector& rows, const Rcpp::IntegerVector& cols,
                   const Rcpp::LogicalVector& flip, int threads) {
  Rcpp::XPtr<Geno> a(ptr);
  threads = resolve_threads(threads);
  const std::size_t n0 = a->n, n = rows.size(), p = cols.size();
  if (flip.size() && static_cast<std::size_t>(flip.size()) != p) Rcpp::stop("flip must have one value per column.");
  const std::vector<int> rw(rows.begin(), rows.end()), cl(cols.begin(), cols.end());
  std::vector<char> fl(p, 0);
  for (int r : rw) if (r < 0 || static_cast<std::size_t>(r) >= n0) Rcpp::stop("rows out of range.");
  for (std::size_t c = 0; c < p; ++c) {
    if (cl[c] < 0 || static_cast<std::size_t>(cl[c]) >= a->p) Rcpp::stop("cols out of range.");
    if (flip.size()) {
      if (flip[c] == NA_LOGICAL) Rcpp::stop("flip must not be NA.");
      fl[c] = flip[c] ? 1 : 0;
      if (fl[c] && !unit_steps(a->D[cl[c]].unit)) Rcpp::stop("Appended columns cannot be flipped.");
    }
  }
  Rcpp::XPtr<Geno> g(new Geno, true);
  init_geno(*g, n, p, a->mean);
  g->dosage = a->dosage;
  const std::size_t W = g->G.words;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<std::int16_t> qbuf(n0, 0);
    std::vector<u64> mbuf(a->G.words, 0);
    #pragma omp for schedule(dynamic, 16)
    for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(p); ++cc) {
      const std::size_t c = static_cast<std::size_t>(cc), j = static_cast<std::size_t>(cl[c]);
      const DCol& dc = a->D[j];
      std::size_t e = 0;
      for (std::size_t w = 0; w < dc.mask.size(); ++w)
        for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) qbuf[64 * w + __builtin_ctzll(bits)] = dc.q[e++];
      for (std::uint32_t i : a->M[j]) mbuf[i >> 6] |= u64(1) << (i & 63u);
      u64* lo = g->G.bits.data() + 2 * c * W;
      u64* hi = lo + W;
      g->D[c].unit = dc.unit;
      DBuilder db(*g, c);
      std::vector<std::uint32_t>& Mc = g->M[c];
      for (std::size_t r = 0; r < n; ++r) {
        const std::size_t i = static_cast<std::size_t>(rw[r]);
        unsigned int h = h_at(a->G, j, i);
        int q = qbuf[i];
        if (fl[c]) flip_cell(h, q, unit_steps(dc.unit));
        if (h == 1u) lo[r >> 6] |= u64(1) << (r & 63u);
        else if (h == 2u) hi[r >> 6] |= u64(1) << (r & 63u);
        if (q) db.add(r, h, q);
        if ((mbuf[i >> 6] >> (i & 63u)) & 1u) Mc.push_back(static_cast<std::uint32_t>(r));
      }
      db.finish();
      plane_stats(*g, c);
      g->nmiss[c] = static_cast<double>(Mc.size());
      e = 0;
      for (std::size_t w = 0; w < dc.mask.size(); ++w)
        for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) qbuf[64 * w + __builtin_ctzll(bits)] = 0;
      for (std::uint32_t i : a->M[j]) mbuf[i >> 6] = 0;
    }
  }
  finish_stats(*g);
  return g;
}

namespace {

// A fill value must keep the int16 part in range.
inline bool fill_value_ok(double f) {
  if (!std::isfinite(f)) return false;
  const double hf = std::min(2.0, std::max(0.0, std::floor(f + 0.5)));
  return std::fabs(16384.0 * (f - hf)) <= 32767.0;
}

// Rewrites the missing cells of column j: a cell takes its per-cell value when given (rows ci[k..end),
// ascending, values cg; the first of repeated rows wins), else fc when hasc, else keeps its value.
template <class RowT>
void fill_column(Geno& g, std::size_t j, bool hasc, double fc, const RowT* ci, const double* cg, std::size_t k,
                 std::size_t end) {
  const std::vector<std::uint32_t>& Mj = g.M[j];
  if (Mj.empty() || (!hasc && k == end)) return;
  const std::size_t m = Mj.size();
  std::vector<int> nq(m, 0);
  std::vector<unsigned char> nh(m, 0), chg(m, 0);
  bool any = false;
  for (std::size_t r = 0; r < m; ++r) {
    const std::size_t i = Mj[r];
    while (k < end && static_cast<std::size_t>(ci[k]) < i) ++k;
    double f;
    if (k < end && static_cast<std::size_t>(ci[k]) == i) f = cg[k];
    else if (hasc) f = fc;
    else continue;
    unsigned int h;
    int q;
    quantize_fill(f, g.D[j].unit, h, q);
    nh[r] = static_cast<unsigned char>(h);
    nq[r] = q;
    chg[r] = 1;
    any = true;
  }
  if (!any) return;
  const std::size_t W = g.G.words;
  u64* lo = g.G.bits.data() + 2 * j * W;
  u64* hi = lo + W;
  for (std::size_t r = 0; r < m; ++r) {
    if (!chg[r]) continue;
    const std::size_t i = Mj[r];
    const u64 bit = u64(1) << (i & 63u);
    lo[i >> 6] &= ~bit;
    hi[i >> 6] &= ~bit;
    if (nh[r] == 1) lo[i >> 6] |= bit;
    else if (nh[r] == 2) hi[i >> 6] |= bit;
  }
  const DCol old = std::move(g.D[j]);
  rebuild_col(g, j, old.unit, [&](auto&& emit) {  // old D entries merged with the new missing-cell values, rows ascending
    std::size_t r = 0, e = 0;
    for (std::size_t w = 0; w < old.mask.size(); ++w)
      for (u64 bits = old.mask[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        const int q = old.q[e++];
        for (; r < m && Mj[r] < i; ++r) if (chg[r]) emit(Mj[r], nq[r]);
        if (r < m && Mj[r] == i) { emit(i, chg[r] ? nq[r] : q); ++r; }
        else emit(i, q);
      }
    for (; r < m; ++r) if (chg[r]) emit(Mj[r], nq[r]);
  });
}

}  // namespace

// The one fill routine: writes the missing cells (rows M[j]) of ptr. Column j's missing cells take the
// per-cell value when given (CSR: cells of column j at cell_ptr[j]..cell_ptr[j+1], 0-based row cell_i
// ascending, value cell_g; the first of repeated rows wins), else fill_c[j] when finite, else keep their
// current value. A value f is stored as H = round(f) (halves up, within 0..2) and D = f - H in 1/16384.
// Cells not missing are never touched. inplace = TRUE changes ptr itself (only for an object nobody else
// references yet); otherwise a modified copy is returned and ptr is unchanged.
// [[Rcpp::export]]
SEXP geno_fill_cpp(SEXP ptr, const Rcpp::NumericVector& fill_c, const Rcpp::IntegerVector& cell_ptr,
                   const Rcpp::IntegerVector& cell_i, const Rcpp::NumericVector& cell_g, bool inplace, int threads) {
  Rcpp::XPtr<Geno> src(ptr);
  threads = resolve_threads(threads);
  const std::size_t p = src->p;
  const bool hc = fill_c.size() > 0, hcell = cell_ptr.size() > 0;
  if (hc && static_cast<std::size_t>(fill_c.size()) != p) Rcpp::stop("fill constants must have one value per variant.");
  if (hcell) {
    if (static_cast<std::size_t>(cell_ptr.size()) != p + 1) Rcpp::stop("cell_ptr must have p + 1 entries.");
    if (cell_i.size() != cell_g.size() || cell_ptr[p] != cell_i.size()) Rcpp::stop("cell_i / cell_g do not match cell_ptr.");
  }
  for (std::size_t j = 0; j < static_cast<std::size_t>(fill_c.size()); ++j)
    if (!std::isnan(fill_c[j]) && !fill_value_ok(fill_c[j])) Rcpp::stop("fill values must be finite and within (-2, 4).");
  for (R_xlen_t k = 0; k < cell_g.size(); ++k)
    if (!fill_value_ok(cell_g[k])) Rcpp::stop("cell fill values must be finite and within (-2, 4).");
  Rcpp::XPtr<Geno> out = inplace ? src : Rcpp::XPtr<Geno>(new Geno(*src), true);
  Geno& g = *out;
  const double* fc = hc ? fill_c.begin() : nullptr;
  const int* cp = hcell ? cell_ptr.begin() : nullptr;
  const int* ci = hcell ? cell_i.begin() : nullptr;
  const double* cg = hcell ? cell_g.begin() : nullptr;
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(p); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    const bool hasc = fc && std::isfinite(fc[j]);
    fill_column(g, j, hasc, hasc ? fc[j] : 0.0, ci, cg, cp ? static_cast<std::size_t>(cp[j]) : 0,
                cp ? static_cast<std::size_t>(cp[j + 1]) : 0);
  }
  return out;
}

// Missing cells of every column: CSR ptr (p + 1), 0-based row i (ascending within a column) and the value
// g the cell holds now.
// [[Rcpp::export]]
Rcpp::List geno_miss_cpp(SEXP ptr) {
  Rcpp::XPtr<Geno> g(ptr);
  const std::size_t p = g->p;
  Rcpp::IntegerVector cp(p + 1);
  std::size_t tot = 0;
  for (std::size_t j = 0; j < p; ++j) { cp[j] = static_cast<int>(tot); tot += g->M[j].size(); }
  cp[p] = static_cast<int>(tot);
  Rcpp::IntegerVector ri(tot);
  Rcpp::NumericVector rg(tot);
  std::size_t e = 0;
  for (std::size_t j = 0; j < p; ++j)
    for (std::uint32_t i : g->M[j]) {
      ri[e] = static_cast<int>(i);
      rg[e] = static_cast<double>(h_at(g->G, j, i)) + d_at(g->D[j], i);
      ++e;
    }
  return Rcpp::List::create(Rcpp::Named("ptr") = cp, Rcpp::Named("i") = ri, Rcpp::Named("g") = rg);
}

namespace plink4r {

void geno_write_cells(SEXP ptr, const std::vector<int>& cols, const std::vector<std::size_t>& off,
                      const std::vector<std::uint32_t>& rows, const std::vector<double>& vals, int threads) {
  Rcpp::XPtr<Geno> gp(ptr);
  Geno& g = *gp;
  threads = resolve_threads(threads);
  const std::size_t m = cols.size();
  if (off.size() != m + 1 || off[m] != rows.size() || rows.size() != vals.size())
    Rcpp::stop("geno_write_cells: cols / off / rows / vals do not match.");
  for (int j : cols) if (j < 0 || static_cast<std::size_t>(j) >= g.p) Rcpp::stop("geno_write_cells: column out of range.");
  for (double v : vals) if (!fill_value_ok(v)) Rcpp::stop("cell fill values must be finite and within (-2, 4).");
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(m); ++cc) {
    const std::size_t c = static_cast<std::size_t>(cc);
    fill_column(g, static_cast<std::size_t>(cols[c]), false, 0.0, rows.data(), vals.data(), off[c], off[c + 1]);
  }
}

void geno_genovec(SEXP ptr, const std::vector<int>& cols, std::vector<unsigned char>& gs,
                  std::vector<double>& sum_obs, int& n_out, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  if (g->dosage) Rcpp::stop("LD imputation needs hard calls; this geno object holds dosages.");
  const std::size_t n = g->n, sb = (n + 3) / 4, m = cols.size();
  for (int j : cols) {
    if (j < 0 || static_cast<std::size_t>(j) >= g->p) Rcpp::stop("cols out of range.");
    const DCol& d = g->D[j];
    if (d.unit != kDUnit) Rcpp::stop("LD imputation does not take appended columns.");
    const std::vector<std::uint32_t>& Mj = g->M[j];
    std::size_t r = 0;
    for (std::size_t w = 0; w < d.mask.size(); ++w)
      for (u64 bits = d.mask[w]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * w + __builtin_ctzll(bits);
        while (r < Mj.size() && Mj[r] < i) ++r;
        if (r == Mj.size() || Mj[r] != i) Rcpp::stop("LD imputation needs hard calls at every observed cell.");
      }
  }
  gs.assign(m * sb, 0);
  sum_obs.assign(m, 0.0);
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(m); ++cc) {
    const std::size_t c = static_cast<std::size_t>(cc), j = static_cast<std::size_t>(cols[c]);
    unsigned char* dst = gs.data() + c * sb;
    // 64 rows per step: code = lo | hi << 1 (h_at: lo wins), i.e. the bits of lo and hi interleaved
    const u64* lo = g->G.lo(j);
    const u64* hi = g->G.hi(j);
    auto spread = [](u64 x) {  // bit t of the low 32 bits -> bit 2t
      x &= 0xFFFFFFFFull;
      x = (x | (x << 16)) & 0x0000FFFF0000FFFFull;
      x = (x | (x << 8)) & 0x00FF00FF00FF00FFull;
      x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0Full;
      x = (x | (x << 2)) & 0x3333333333333333ull;
      x = (x | (x << 1)) & 0x5555555555555555ull;
      return x;
    };
    const std::size_t nw = (n + 63) / 64;
    for (std::size_t w = 0; w < nw; ++w) {
      u64 l = lo[w], h = hi[w] & ~l;
      const std::size_t rows = std::min<std::size_t>(64, n - 64 * w);
      if (rows < 64) { const u64 keep = (u64(1) << rows) - 1; l &= keep; h &= keep; }
      const u64 out[2] = {spread(l) | (spread(h) << 1), spread(l >> 32) | (spread(h >> 32) << 1)};
      const std::size_t b0 = 16 * w, nbyt = std::min<std::size_t>(16, sb - b0);
      unsigned char tmp[16];
      for (int k = 0; k < 2; ++k)
        for (int b = 0; b < 8; ++b) tmp[8 * k + b] = static_cast<unsigned char>(out[k] >> (8 * b));
      std::memcpy(dst + b0, tmp, nbyt);
    }
    std::int64_t s = g->G.s[j];
    for (std::uint32_t i : g->M[j]) {
      s -= h_at(g->G, j, i);
      dst[i >> 2] |= static_cast<unsigned char>(3u << (2 * (i & 3u)));
    }
    sum_obs[c] = static_cast<double>(s);
  }
  n_out = static_cast<int>(n);
}

}  // namespace plink4r

// C-callable entry points for other packages' C / C++ code (inst/include/plink4r_api.h; LinkingTo: plink4r).
// Each takes the external pointer of a geno object (X$ptr) and raw column-major buffers, writes its result
// into the caller's buffer (no result is allocated here) and returns 0, or a PLINK4R_API_E* code without
// touching R: 1 not a live geno pointer, 2 bad argument (size / index), 3 failure inside (e.g. memory).
// scale is X$scale; with scale the products are those of the column-standardized matrix, from the
// store's mean and sd (= X$center, X$sd; a column with sd 0 is only centered), as in R/geno.R.
namespace {

const Geno* api_geno(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP) return nullptr;
  return static_cast<const Geno*>(R_ExternalPtrAddr(ptr));
}

bool api_cols_ok(const Geno& g, const int* cols, int nc) {
  if (nc < 0 || (nc > 0 && !cols)) return false;
  for (int a = 0; a < nc; ++a)
    if (cols[a] < 0 || static_cast<std::size_t>(cols[a]) >= g.p) return false;
  return true;
}

// .geno_div(): sd, or 1 for a constant column.
std::vector<double> api_div(const Geno& g) {
  std::vector<double> d(g.p);
  for (std::size_t j = 0; j < g.p; ++j) d[j] = g.sd[j] > 0 ? g.sd[j] : 1.0;
  return d;
}

// R's colSums (long double accumulation in row order) of an r x c column-major matrix.
void api_colsums(const double* x, std::size_t r, std::size_t c, double* out) {
  for (std::size_t j = 0; j < c; ++j) {
    long double s = 0.0;
    for (std::size_t i = 0; i < r; ++i) s += x[i + j * r];
    out[j] = static_cast<double>(s);
  }
}

template <class F>
int api_guard(F&& f) {
  try {
    return f();
  } catch (...) {
    return 3;
  }
}

// n, p, storage (1 = dosage PGEN, 0 = hard calls) and per column mean, sd (n - 1) and missing count of the
// file; any pointer may be NULL (call with NULL arrays first to size them, p values each).
int api_geno_info(SEXP ptr, int* n, int* p, int* dosage, double* mean, double* sd, double* n_miss) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (n) *n = static_cast<int>(g->n);
  if (p) *p = static_cast<int>(g->p);
  if (dosage) *dosage = g->dosage ? 1 : 0;
  if (mean) std::copy(g->mu.begin(), g->mu.end(), mean);
  if (sd) std::copy(g->sd.begin(), g->sd.end(), sd);
  if (n_miss) std::copy(g->nmiss.begin(), g->nmiss.end(), n_miss);
  return 0;
}

// out (p x k) = X'M for M (n x k); geno_crossprod(X, M).
int api_geno_xtm(SEXP ptr, int scale, const double* M, int k, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (k < 0 || (k > 0 && (!M || !out))) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t kk = static_cast<std::size_t>(k);
    xtm_into(*g, M, kk, out, threads);
    if (scale) {
      std::vector<double> cs(kk), d = api_div(*g);
      api_colsums(M, g->n, kk, cs.data());
      std_into(out, g->p, kk, g->mu.data(), cs.data(), 1.0, d.data(), nullptr, out, threads);
    }
    return 0;
  });
}

// out (n x k) = X B for B (p x k); geno_multiply(X, B).
int api_geno_mx(SEXP ptr, int scale, const double* B, int k, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (k < 0 || (k > 0 && (!B || !out))) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t n = g->n, p = g->p, kk = static_cast<std::size_t>(k);
    if (!scale) {
      mx_into(*g, B, kk, out, threads);
      return 0;
    }
    const std::vector<double> d = api_div(*g);
    std::vector<double> Bs(p * kk), Bc(p * kk), cs(kk);
    for (std::size_t c = 0; c < kk; ++c)
      for (std::size_t j = 0; j < p; ++j) {
        Bs[j + c * p] = B[j + c * p] / d[j];
        Bc[j + c * p] = Bs[j + c * p] * g->mu[j];
      }
    api_colsums(Bc.data(), p, kk, cs.data());
    mx_into(*g, Bs.data(), kk, out, threads);
    for (std::size_t c = 0; c < kk; ++c)
      for (std::size_t i = 0; i < n; ++i) out[i + c * n] = out[i + c * n] - cs[c];
    return 0;
  });
}

// out (q x q) = X'X of columns cols (0-based, q = nc), or of all p columns when cols is NULL, standardized in
// place in out; geno_crossprod(X)[cols + 1, cols + 1].
int api_geno_xtx(SEXP ptr, int scale, const int* cols, int nc, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (!out || (cols && !api_cols_ok(*g, cols, nc))) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t q = cols ? static_cast<std::size_t>(nc) : g->p;
    if (!q) return 0;
    xtx_into(*g, cols, q, out, threads);
    if (scale) {
      const std::vector<double> d = api_div(*g);
      std::vector<double> u(q), dq(q);
      for (std::size_t a = 0; a < q; ++a) {
        const std::size_t j = cols ? static_cast<std::size_t>(cols[a]) : a;
        u[a] = g->mu[j];
        dq[a] = d[j];
      }
      std_into(out, q, q, u.data(), u.data(), static_cast<double>(g->n), dq.data(), dq.data(), out, threads);
    }
    return 0;
  });
}

// out (q x q) = X'W~X of columns cols (0-based, q = nc), or of all p columns when cols is NULL, with w (n)
// binned in groups of 64 rows, or w itself when D is dense (xtwx_bin_into), standardized in place with those
// weights when scale; w_bin (n, may be NULL) receives them. geno_wcrossprod_bin(X, w)[cols + 1, cols + 1].
int api_geno_xtwx_bin(SEXP ptr, int scale, const double* w, const int* cols, int nc, double* out, double* w_bin,
                      int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (!out || !w || (cols && !api_cols_ok(*g, cols, nc))) return 2;
  for (std::size_t i = 0; i < g->n; ++i) if (!std::isfinite(w[i])) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t q = cols ? static_cast<std::size_t>(nc) : g->p;
    std::vector<double> wq(g->n);
    xtwx_bin_into(*g, w, cols, q, out, wq.data(), threads);
    if (w_bin) std::copy(wq.begin(), wq.end(), w_bin);
    if (scale) {
      const std::vector<double> d = api_div(*g);
      std::vector<double> u(q), dq(q);
      for (std::size_t a = 0; a < q; ++a) {
        const std::size_t j = cols ? static_cast<std::size_t>(cols[a]) : a;
        u[a] = g->mu[j];
        dq[a] = d[j];
      }
      xtwx_std(*g, wq.data(), cols, q, u.data(), dq.data(), out, threads);
    }
    return 0;
  });
}

// out (n x nc) = X[, cols + 1] (cols 0-based), standardized when scale.
int api_geno_cols(SEXP ptr, int scale, const int* cols, int nc, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (!api_cols_ok(*g, cols, nc) || (nc > 0 && !out)) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t n = g->n, q = static_cast<std::size_t>(nc);
    std::vector<int> rows(n);
    for (std::size_t i = 0; i < n; ++i) rows[i] = static_cast<int>(i);
    dense_into(*g, rows.data(), n, cols, q, out, threads);
    if (scale)
      for (std::size_t a = 0; a < q; ++a) {
        const std::size_t j = static_cast<std::size_t>(cols[a]);
        const double mu = g->mu[j], d = g->sd[j] > 0 ? g->sd[j] : 1.0;
        double* o = out + a * n;
        for (std::size_t i = 0; i < n; ++i) o[i] = o[i] - mu;
        for (std::size_t i = 0; i < n; ++i) o[i] = o[i] / d;
      }
    return 0;
  });
}

}  // namespace

// Called from R_init_plink4r (RcppExports.cpp).
// [[Rcpp::init]]
void plink4r_register_ccallables(DllInfo* dll) {
  (void) dll;
  R_RegisterCCallable("plink4r", "plink4r_geno_info", reinterpret_cast<DL_FUNC>(&api_geno_info));
  R_RegisterCCallable("plink4r", "plink4r_geno_xtm", reinterpret_cast<DL_FUNC>(&api_geno_xtm));
  R_RegisterCCallable("plink4r", "plink4r_geno_mx", reinterpret_cast<DL_FUNC>(&api_geno_mx));
  R_RegisterCCallable("plink4r", "plink4r_geno_xtx", reinterpret_cast<DL_FUNC>(&api_geno_xtx));
  R_RegisterCCallable("plink4r", "plink4r_geno_cols", reinterpret_cast<DL_FUNC>(&api_geno_cols));
  R_RegisterCCallable("plink4r", "plink4r_geno_xtwx_bin", reinterpret_cast<DL_FUNC>(&api_geno_xtwx_bin));
}

// Read-only: the threads and memory budget a call with threads = 0 / mem_gb = 0 would use.
// [[Rcpp::export]]
Rcpp::List resources_cpp(int threads, double mem_gb) {
  const plink4r::Resources r = plink4r::detect_resources(threads, mem_gb);
  return Rcpp::List::create(Rcpp::Named("threads") = r.threads, Rcpp::Named("mem_gb") = r.mem_gb,
                            Rcpp::Named("threads_from") = r.threads_from, Rcpp::Named("mem_from") = r.mem_from);
}
