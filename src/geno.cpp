// geno: a PLINK BED / PGEN genotype matrix held in memory as X = H + D with
// the products SuSiE-style fine-mapping needs: X'X by exact popcount plus the
// D terms, X'M and X B by plane scans, and dense sub-blocks on demand.
// H is a 0/1/2 matrix in 2-bit planes (n p / 4 bytes). D is a sparse
// deviation: per variant a 1-bit nonzero mask and int16 values in units of
// that variant's `unit` (1/16384 for genotypes), allocated only for variants
// with a nonzero D. It holds
//   - dosage PGEN: dosage - H, H the rounded dosage (exact);
//   - every input: the fill of a missing entry, f - H with H = round(f) and
//     f the variant's lower median (impute = "median", D = 0) or observed
//     mean (impute = "mean", and dosage PGEN; quantized to 1/16384);
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
#include <fstream>
#include <string>
#include <utility>
#include "geno_kernel.h"
#include "geno_access.h"
#include "pgen_reader.h"
#include "resources.h"
#include "la_dispatch.h"

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
constexpr double kDUnit = 1.0 / 16384.0;    // value of one int16 step of D

// D of one variant: nonzero mask over the samples (H words), int16 values
// q = D / unit in ascending row order, and the nonzero count before each block
// of kDBlock mask words. All empty when the variant's D is zero. unit is
// 1/16384 for genotype columns; an appended column (H = 0, all of x in D)
// has its own unit max|x| / 32767.
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
// deviation q = round(16384 (f - H)).
inline void quantize_fill(double f, unsigned int& h, int& q) {
  const double hf = std::min(2.0, std::max(0.0, std::floor(f + 0.5)));
  h = static_cast<unsigned int>(hf);
  q = static_cast<int>(std::lround(16384.0 * (f - hf)));
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
  quantize_fill(f, h, q);
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

// Dosage PGEN: H = dosage rounded to the nearest integer (halves up) into the
// planes, D = dosage - H exactly as int16 in [-8192, 8191]. Missing dosages
// take the variant's observed mean dosage, through D; their rows go to M[j].
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
        else obs_sum += v;
      }
      if (bad_range) {
        #pragma omp atomic write
        range = 1;
        continue;
      }
      g.nmiss[j] = static_cast<double>(miss);
      unsigned int hfill = 0;
      int qfill = 0;
      const double fmean = miss < n ? static_cast<double>(obs_sum) / 16384.0 / static_cast<double>(n - miss) : 0.0;
      if (miss) quantize_fill(fmean, hfill, qfill);
      g.M[j].reserve(miss);
      unsigned char* buf = packed.data() + static_cast<std::size_t>(t) * sel_stride;
      std::fill(buf, buf + sel_stride, 0);
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
          h = (static_cast<unsigned int>(v) + 8192u) >> 14;
          q = static_cast<int>(v) - static_cast<int>(h << 14);
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

// X'X += H'D + D'H + D'D, as S + S' with S[a, b] = sum_i D_bi (H_ai + D_ai / 2),
// from a row-interleaved buffer of T columns of H + D / 2 per thread.
void add_d_xtx(const Geno& g, double* out, int threads) {
  const std::size_t n = g.n, p = g.p;
  if (!g.dnnz()) return;
  constexpr std::size_t T = 8;
  std::vector<double> S(p * p, 0.0);
  const std::size_t tiles = (p + T - 1) / T;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> buf(n * T);
    #pragma omp for schedule(dynamic, 1)
    for (std::ptrdiff_t tt = 0; tt < static_cast<std::ptrdiff_t>(tiles); ++tt) {
      const std::size_t a0 = static_cast<std::size_t>(tt) * T, na = std::min(T, p - a0);
      std::fill(buf.begin(), buf.end(), 0.0);
      for (std::size_t t = 0; t < na; ++t) {
        const std::size_t a = a0 + t;
        const u64* lo = g.G.lo(a);
        const u64* hi = g.G.hi(a);
        for (std::size_t w = 0; w < g.G.words; ++w) {
          for (u64 bits = lo[w]; bits; bits &= bits - 1)
            buf[(64 * w + __builtin_ctzll(bits)) * T + t] = 1.0;
          for (u64 bits = hi[w]; bits; bits &= bits - 1)
            buf[(64 * w + __builtin_ctzll(bits)) * T + t] = 2.0;
        }
        const DCol& c = g.D[a];
        std::size_t e = 0;
        for (std::size_t w = 0; w < c.mask.size(); ++w)
          for (u64 bits = c.mask[w]; bits; bits &= bits - 1)
            buf[(64 * w + __builtin_ctzll(bits)) * T + t] += 0.5 * (c.q[e++] * c.unit);
      }
      for (std::size_t b = 0; b < p; ++b) {
        const DCol& c = g.D[b];
        std::size_t e = 0;
        double acc[T] = {};
        for (std::size_t w = 0; w < c.mask.size(); ++w)
          for (u64 bits = c.mask[w]; bits; bits &= bits - 1) {
            const double d = c.q[e++] * c.unit;
            const double* r = buf.data() + (64 * w + __builtin_ctzll(bits)) * T;
            for (std::size_t t = 0; t < T; ++t) acc[t] += d * r[t];
          }
        for (std::size_t t = 0; t < na; ++t) S[a0 + t + b * p] = acc[t];
      }
    }
    #pragma omp barrier
    #pragma omp for schedule(static)
    for (std::ptrdiff_t bb = 0; bb < static_cast<std::ptrdiff_t>(p); ++bb) {
      const std::size_t b = static_cast<std::size_t>(bb);
      for (std::size_t a = 0; a < p; ++a) out[a + b * p] += S[a + b * p] + S[b + a * p];
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

// Sum of the stored values of column j over its observed (not missing) cells. Genotype columns (unit
// 1/16384) are summed exactly in int16 steps, so the result equals the sum of the file's observed calls
// (hard calls) or dosages; an appended column (no missing cells) gives its column sum.
double obs_sum_col(const Geno& g, std::size_t j) {
  const DCol& c = g.D[j];
  if (c.unit != kDUnit) return static_cast<double>(g.G.s[j]) + g.dsum[j];
  std::int64_t tot = 16384 * static_cast<std::int64_t>(g.G.s[j]);
  for (std::int16_t q : c.q) tot += q;
  for (std::uint32_t i : g.M[j]) tot -= 16384 * static_cast<std::int64_t>(h_at(g.G, j, i)) + q_at(c, i);
  return static_cast<double>(tot) / 16384.0;
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

// Exact flip of one stored cell x = h + q / 16384 -> 2 - x, in the canonical form of a direct load
// (q = +8192 becomes h + 1, q = -8192).
inline void flip_cell(unsigned int& h, int& q) {
  h = 2u - h;
  q = -q;
  if (q == 8192 && h < 2u) { ++h; q = -8192; }
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
  const TileKernel kernel = select_kernel();
  const std::size_t p = g->p;
  Rcpp::NumericMatrix out(p, p);
  dot_block(g->G, 0, p, g->G, 0, p, true, out.begin(), p, threads, kernel);
  add_d_xtx(*g, out.begin(), threads);
  return out;
}

// X'M for a dense n x k matrix M.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtm_cpp(SEXP ptr, const Rcpp::NumericMatrix& M, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p, W = g->G.words;
  const std::size_t k = M.ncol();
  if (static_cast<std::size_t>(M.nrow()) != n) Rcpp::stop("M must have one row per sample.");
  Rcpp::NumericMatrix out(p, k);
  const double* m = M.begin();
  double* o = out.begin();
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
  const double* b = B.begin();
  double* o = out.begin();
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
  double* o = out.begin();
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
  const double* k = K.begin();
  double* o = out.begin();
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(c); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    const double vj = s * v[j], bj = hb ? db[j] : 1.0;
    for (std::size_t i = 0; i < r; ++i) o[i + j * r] = (k[i + j * r] - u[i] * vj) / ((ha ? da[i] : 1.0) * bj);
  }
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

namespace {

// Rows r0..r0+R-1 of stored column j (H + D) as floats: dst[r * stride], times w[r0 + r] when w is given.
void decode_rows(const Geno& g, std::size_t j, std::size_t r0, std::size_t R, const double* w, float* dst,
                 std::size_t stride) {
  for (std::size_t r = 0; r < R; ++r) dst[r * stride] = 0.0f;
  const u64* lo = g.G.lo(j);
  const u64* hi = g.G.hi(j);
  const std::size_t w0 = r0 / 64, w1 = (r0 + R + 63) / 64;
  for (std::size_t x = w0; x < w1; ++x) {
    for (u64 bits = lo[x]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * x + __builtin_ctzll(bits);
      if (i >= r0 && i < r0 + R) dst[(i - r0) * stride] = 1.0f;
    }
    for (u64 bits = hi[x]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * x + __builtin_ctzll(bits);
      if (i >= r0 && i < r0 + R) dst[(i - r0) * stride] = 2.0f;
    }
  }
  const DCol& dc = g.D[j];
  if (!dc.mask.empty()) {
    std::size_t e = d_rank(dc, w0);
    for (std::size_t x = w0; x < w1 && x < dc.mask.size(); ++x)
      for (u64 bits = dc.mask[x]; bits; bits &= bits - 1) {
        const std::size_t i = 64 * x + __builtin_ctzll(bits);
        const double d = dc.q[e++] * dc.unit;
        if (i >= r0 && i < r0 + R) dst[(i - r0) * stride] += static_cast<float>(d);
      }
  }
  if (w)
    for (std::size_t r = 0; r < R; ++r) dst[r * stride] = static_cast<float>(dst[r * stride] * w[r0 + r]);
}

}  // namespace

// X'WX of the stored matrix, W = diag(w), any real w. Output tiles of 256 x 256 (lower triangle) are shared
// out by OpenMP within each block of 4096 samples; a tile decodes only its two groups of 256 columns for that
// block into float (w folded into one side), multiplies them with the run-time-dispatched serial sgemm, and
// accumulates in double; ~8 MB per thread whatever n and p. With center (and div) given, returns the
// product of the standardized matrix: (X'WX - c s' - s c' + (sum w) c c') / (d d'), s = X'w.
// [[Rcpp::export]]
Rcpp::NumericMatrix geno_xtwx_cpp(SEXP ptr, const Rcpp::NumericVector& w, const Rcpp::NumericVector& center,
                                  const Rcpp::NumericVector& div, int threads) {
  Rcpp::XPtr<Geno> g(ptr);
  threads = resolve_threads(threads);
  const std::size_t n = g->n, p = g->p;
  if (static_cast<std::size_t>(w.size()) != n) Rcpp::stop("w must have one value per sample.");
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(w[i])) Rcpp::stop("w must be finite.");
  const plink4r::LaFns& LA = plink4r::la();
  constexpr std::size_t TB = 256, RB = 4096;
  const std::size_t nt = (p + TB - 1) / TB;
  std::vector<std::pair<std::size_t, std::size_t>> tiles;
  for (std::size_t b = 0; b < nt; ++b) for (std::size_t a = b; a < nt; ++a) tiles.emplace_back(a, b);
  Rcpp::NumericMatrix out(p, p);
  double* o = out.begin();
  const double* wp = w.begin();
  #pragma omp parallel num_threads(threads)
  {
    std::vector<float> At(TB * RB), B(RB * TB), P(TB * TB);
    for (std::size_t r0 = 0; r0 < n; r0 += RB) {
      const std::size_t R = std::min(RB, n - r0);
      #pragma omp for schedule(dynamic, 1)
      for (std::ptrdiff_t tt = 0; tt < static_cast<std::ptrdiff_t>(tiles.size()); ++tt) {
        const std::size_t a0 = tiles[tt].first * TB, b0 = tiles[tt].second * TB;
        const std::size_t na = std::min(TB, p - a0), nb = std::min(TB, p - b0);
        for (std::size_t c = 0; c < na; ++c) decode_rows(*g, a0 + c, r0, R, wp, &At[c], TB);    // A' : na x R, ld TB
        for (std::size_t c = 0; c < nb; ++c) decode_rows(*g, b0 + c, r0, R, nullptr, &B[c * RB], 1);  // B : R x nb, ld RB
        LA.sgemm_nn(At.data(), static_cast<int>(na), static_cast<int>(R), static_cast<int>(TB), B.data(),
                    static_cast<int>(nb), static_cast<int>(RB), P.data(), static_cast<int>(TB));
        for (std::size_t y = 0; y < nb; ++y)
          for (std::size_t x = 0; x < na; ++x) o[(a0 + x) + (b0 + y) * p] += P[x + y * TB];
      }  // implicit barrier: the next sample block starts after every tile of this one
    }
  }
  for (std::size_t b = 0; b < p; ++b)  // mirror the lower triangle
    for (std::size_t a = b + 1; a < p; ++a) o[b + a * p] = o[a + b * p];
  if (center.size() > 0) {
    if (static_cast<std::size_t>(center.size()) != p) Rcpp::stop("center must have one value per column.");
    const bool hd = div.size() > 0;
    std::vector<double> sw(p, 0.0);
    double wsum = 0.0;
    for (std::size_t i = 0; i < n; ++i) wsum += wp[i];
    #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
    for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(p); ++jj) {  // s = X'w in double
      const std::size_t j = static_cast<std::size_t>(jj);
      double acc = 0.0;
      const u64* lo = g->G.lo(j);
      const u64* hi = g->G.hi(j);
      for (std::size_t x = 0; x < g->G.words; ++x) {
        for (u64 bits = lo[x]; bits; bits &= bits - 1) acc += wp[64 * x + __builtin_ctzll(bits)];
        for (u64 bits = hi[x]; bits; bits &= bits - 1) acc += 2.0 * wp[64 * x + __builtin_ctzll(bits)];
      }
      const DCol& dc = g->D[j];
      std::size_t e = 0;
      for (std::size_t x = 0; x < dc.mask.size(); ++x)
        for (u64 bits = dc.mask[x]; bits; bits &= bits - 1) acc += dc.q[e++] * dc.unit * wp[64 * x + __builtin_ctzll(bits)];
      sw[j] = acc;
    }
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t bb = 0; bb < static_cast<std::ptrdiff_t>(p); ++bb) {
      const std::size_t b = static_cast<std::size_t>(bb);
      for (std::size_t a = 0; a < p; ++a) {
        double v = o[a + b * p] - center[a] * sw[b] - sw[a] * center[b] + wsum * center[a] * center[b];
        if (hd) v /= div[a] * div[b];
        o[a + b * p] = v;
      }
    }
  }
  return out;
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
      if (fl[c] && a->D[cl[c]].unit != kDUnit) Rcpp::stop("Appended columns cannot be flipped.");
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
        if (fl[c]) flip_cell(h, q);
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
  auto ok_value = [](double f) {  // the int16 part of every fill must fit
    if (!std::isfinite(f)) return false;
    const double hf = std::min(2.0, std::max(0.0, std::floor(f + 0.5)));
    return std::fabs(16384.0 * (f - hf)) <= 32767.0;
  };
  for (std::size_t j = 0; j < static_cast<std::size_t>(fill_c.size()); ++j)
    if (!std::isnan(fill_c[j]) && !ok_value(fill_c[j])) Rcpp::stop("fill values must be finite and within (-2, 4).");
  for (R_xlen_t k = 0; k < cell_g.size(); ++k)
    if (!ok_value(cell_g[k])) Rcpp::stop("cell fill values must be finite and within (-2, 4).");
  Rcpp::XPtr<Geno> out = inplace ? src : Rcpp::XPtr<Geno>(new Geno(*src), true);
  Geno& g = *out;
  const double* fc = hc ? fill_c.begin() : nullptr;
  const int* cp = hcell ? cell_ptr.begin() : nullptr;
  const int* ci = hcell ? cell_i.begin() : nullptr;
  const double* cg = hcell ? cell_g.begin() : nullptr;
  const std::size_t W = g.G.words;
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(p); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    const std::vector<std::uint32_t>& Mj = g.M[j];
    if (Mj.empty()) continue;
    const bool hasc = fc && std::isfinite(fc[j]);
    int k = cp ? cp[j] : 0;
    const int end = cp ? cp[j + 1] : 0;
    if (!hasc && k == end) continue;
    const std::size_t m = Mj.size();
    std::vector<int> nq(m, 0);
    std::vector<unsigned char> nh(m, 0), chg(m, 0);
    bool any = false;
    for (std::size_t r = 0; r < m; ++r) {
      const int i = static_cast<int>(Mj[r]);
      while (k < end && ci[k] < i) ++k;
      double f;
      if (k < end && ci[k] == i) f = cg[k];
      else if (hasc) f = fc[j];
      else continue;
      unsigned int h;
      int q;
      quantize_fill(f, h, q);
      nh[r] = static_cast<unsigned char>(h);
      nq[r] = q;
      chg[r] = 1;
      any = true;
    }
    if (!any) continue;
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
    for (std::size_t i = 0; i < n; ++i) dst[i >> 2] |= static_cast<unsigned char>(h_at(g->G, j, i) << (2 * (i & 3u)));
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

// Read-only: the threads and memory budget a call with threads = 0 / mem_gb = 0 would use.
// [[Rcpp::export]]
Rcpp::List resources_cpp(int threads, double mem_gb) {
  const plink4r::Resources r = plink4r::detect_resources(threads, mem_gb);
  return Rcpp::List::create(Rcpp::Named("threads") = r.threads, Rcpp::Named("mem_gb") = r.mem_gb,
                            Rcpp::Named("threads_from") = r.threads_from, Rcpp::Named("mem_from") = r.mem_from);
}
