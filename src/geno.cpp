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

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <utility>
#include "geno_kernel.h"
#include "pgen_file.h"
#include "resources.h"
#include "la_dispatch.h"

namespace {

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

// D of column c at sample i.
inline double d_at(const DCol& c, std::size_t i) {
  if (c.mask.empty()) return 0.0;
  const std::size_t w = i / 64, b = i % 64;
  if (!((c.mask[w] >> b) & 1u)) return 0.0;
  const std::size_t e = d_rank(c, w) + popcount64(c.mask[w] & ((u64(1) << b) - 1));
  return c.q[e] * c.unit;
}

inline int omp_thread() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline int omp_threads() {
#ifdef _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
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

// Fills the missing entries (mask miss) of variant j with fill(i): round(f)
// into the planes, the rest into D. Today fill is the variant's lower median
// or observed mean; a per-entry fill takes the same path.
template <class Fill>
void fill_missing(Geno& g, std::size_t j, const u64* miss, Fill fill, DBuilder& db) {
  const std::size_t W = g.G.words;
  u64* lo = g.G.bits.data() + 2 * j * W;
  u64* hi = lo + W;
  for (std::size_t w = 0; w < W; ++w)
    for (u64 bits = miss[w]; bits; bits &= bits - 1) {
      const std::size_t i = 64 * w + __builtin_ctzll(bits);
      unsigned int h;
      int q;
      quantize_fill(fill(i), h, q);
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

// Externally supplied fills (plink4r::genotype_fill): per-variant constant c[j]
// (NaN = the built-in median / mean) and per-cell values in CSR form, cells of
// variant j at ptr[j]..ptr[j+1] with selected-row index ci (ascending) and value cg.
// A missing cell takes its per-cell value when given, else the constant.
// af[j] (finite) anchors the fill instead (scheme 3): the missing calls of
// variant j all take c = (2 n af - sum of observed) / (number missing),
// clipped to [0, 2], so the filled column has mean 2 af.
struct FillSpec {
  std::vector<double> c, af;
  std::vector<int> ptr, ci;
  std::vector<double> cg;
  bool has_c(std::size_t j) const { return !c.empty() && std::isfinite(c[j]); }
  bool has_cells() const { return !ptr.empty(); }
  bool has_af(std::size_t j) const { return !af.empty() && std::isfinite(af[j]); }
};

inline double anchor_fill(double af, double n, double obs_sum, double miss) {
  return std::min(2.0, std::max(0.0, (2.0 * n * af - obs_sum) / miss));
}

// Fill of missing row i of variant j, cursor k walking the variant's cells (rows ascending).
inline double cell_fill(const FillSpec* fs, std::size_t j, std::size_t i, int& k, double f) {
  if (!fs || !fs->has_cells()) return f;
  const int end = fs->ptr[j + 1];
  while (k < end && fs->ci[k] < static_cast<int>(i)) ++k;
  return (k < end && fs->ci[k] == static_cast<int>(i)) ? fs->cg[k] : f;
}

// Hardcall column: planes of the observed calls, then the missing calls filled
// with the variant's lower median or observed mean, or the supplied fills.
void encode_column(Geno& g, const NibbleTable& table, const unsigned char* buf,
                   std::size_t j, u64* miss, const FillSpec* fs = nullptr) {
  const std::size_t W = g.G.words, n = g.n;
  u64* lo = g.G.bits.data() + 2 * j * W;
  std::int64_t s, v;
  encode_packed(table, buf, n, (n + 3) / 4, W, lo, lo + W, s, v, miss);
  std::int64_t cm = 0;
  for (std::size_t w = 0; w < W; ++w) cm += popcount64(miss[w]);
  const std::int64_t obs = static_cast<std::int64_t>(n) - cm;
  g.nmiss[j] = static_cast<double>(cm);
  DBuilder db(g, j);
  if (cm && (obs || (fs && (fs->has_c(j) || fs->has_cells() || fs->has_af(j))))) {
    std::int64_t c1 = 0, c2 = 0;
    for (std::size_t w = 0; w < W; ++w) {
      c1 += popcount64(lo[w]);
      c2 += popcount64(lo[W + w]);
    }
    const std::int64_t c0 = obs - c1 - c2;
    double f = !obs ? 0.0 : g.mean ? static_cast<double>(c1 + 2 * c2) / static_cast<double>(obs)
      : (2 * c0 >= obs ? 0.0 : (2 * (c0 + c1) >= obs ? 1.0 : 2.0));
    if (fs && fs->has_c(j)) f = fs->c[j];
    else if (fs && fs->has_af(j))
      f = anchor_fill(fs->af[j], static_cast<double>(n), static_cast<double>(c1 + 2 * c2), static_cast<double>(cm));
    int k = fs && fs->has_cells() ? fs->ptr[j] : 0;
    fill_missing(g, j, miss, [fs, j, f, &k](std::size_t i) { return cell_fill(fs, j, i, k, f); }, db);
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
}

void finish_stats(Geno& g) {
  const double n = static_cast<double>(g.n);
  g.mu.assign(g.p, 0.0);
  g.sd.assign(g.p, 0.0);
  for (std::size_t j = 0; j < g.p; ++j) {
    const double s = static_cast<double>(g.G.s[j]);
    const double v = static_cast<double>(g.G.v[j]);
    const double mean = (s + g.dsum[j]) / n;
    const double sumsq = (v + s * s) / n + 2.0 * g.dhd[j] + g.dsq[j];
    g.mu[j] = mean;
    const double var = (sumsq - n * mean * mean) / (n - 1);
    g.sd[j] = var > 0 ? std::sqrt(var) : 0.0;
  }
}

std::vector<int> to_rows(const Rcpp::IntegerVector& x) {
  return std::vector<int>(x.begin(), x.end());
}

// Dosage PGEN: H = dosage rounded to the nearest integer (halves up) into the
// planes, D = dosage - H exactly as int16 in [-8192, 8191]. Missing dosages
// take the variant's observed mean dosage, through D.
void load_dosage(Geno& g, PgenFile& file, const Rcpp::IntegerVector& variants,
                 const std::vector<int>& rows, int threads, const std::string& path,
                 const FillSpec* fs = nullptr) {
  const std::size_t n = g.n, p = g.p, sel_stride = (n + 3) / 4, W = g.G.words;
  g.dosage = true;
  const NibbleTable& table = PgenFile::table();
  const std::uint16_t kMiss = PgenFile::kDosageMissing;
  std::vector<unsigned char> packed(static_cast<std::size_t>(threads) * sel_stride);
  for (std::size_t g0 = 0; g0 < p; g0 += kGroup) {
    const std::size_t count = std::min(kGroup, p - g0);
    int failed = 0, range = 0;
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const int t = omp_thread();
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
      double fmean = miss < n ? static_cast<double>(obs_sum) / 16384.0 / static_cast<double>(n - miss) : 0.0;
      if (fs && fs->has_c(j)) fmean = fs->c[j];
      else if (miss && fs && fs->has_af(j))
        fmean = anchor_fill(fs->af[j], static_cast<double>(n), static_cast<double>(obs_sum) / 16384.0, static_cast<double>(miss));
      if (miss) quantize_fill(fmean, hfill, qfill);
      const bool percell = miss && fs && fs->has_cells() && fs->ptr[j + 1] > fs->ptr[j];
      int kc = percell ? fs->ptr[j] : 0;
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
          if (percell) quantize_fill(cell_fill(fs, j, i, kc, fmean), h, q);
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

FillSpec make_fill(const Rcpp::NumericVector& fill_c, const Rcpp::IntegerVector& cell_ptr,
                   const Rcpp::IntegerVector& cell_i, const Rcpp::NumericVector& cell_g,
                   const Rcpp::NumericVector& fill_af, std::size_t p) {
  FillSpec fs;
  if (fill_af.size()) {
    if (static_cast<std::size_t>(fill_af.size()) != p) Rcpp::stop("af must have one value per variant.");
    fs.af.assign(fill_af.begin(), fill_af.end());
  }
  if (fill_c.size()) {
    if (static_cast<std::size_t>(fill_c.size()) != p) Rcpp::stop("fill constants must have one value per variant.");
    fs.c.assign(fill_c.begin(), fill_c.end());
  }
  if (cell_ptr.size()) {
    if (static_cast<std::size_t>(cell_ptr.size()) != p + 1) Rcpp::stop("cell_ptr must have p + 1 entries.");
    fs.ptr.assign(cell_ptr.begin(), cell_ptr.end());
    fs.ci.assign(cell_i.begin(), cell_i.end());
    fs.cg.assign(cell_g.begin(), cell_g.end());
  }
  return fs;
}

} // namespace

// [[Rcpp::export]]
SEXP geno_open_bed_cpp(const std::string& path, int n_total, int m_total,
                       const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples,
                       bool mean, int threads, Rcpp::NumericVector fill_c, Rcpp::IntegerVector cell_ptr,
                       Rcpp::IntegerVector cell_i, Rcpp::NumericVector cell_g, Rcpp::NumericVector fill_af) {
  threads = resolve_threads(threads);
  const FillSpec fs = make_fill(fill_c, cell_ptr, cell_i, cell_g, fill_af, variants.size());
  const std::size_t stride = (static_cast<std::size_t>(n_total) + 3) / 4;
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input) Rcpp::stop("Cannot open BED file: " + path);
  const std::streamoff actual = input.tellg();
  if (actual < 0 || static_cast<std::uint64_t>(actual) != 3 + static_cast<std::uint64_t>(m_total) * stride)
    Rcpp::stop("BED size does not match its BIM and FAM files: " + path);
  input.seekg(0);
  unsigned char header[3];
  input.read(reinterpret_cast<char*>(header), 3);
  if (header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    Rcpp::stop("BED must be a PLINK SNP-major BED file: " + path);

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
      input.seekg(static_cast<std::streamoff>(3 + static_cast<std::uint64_t>(variants[g0 + k]) * stride));
      input.read(reinterpret_cast<char*>(buffer.data() + k * stride), static_cast<std::streamsize>(stride));
      if (!input) Rcpp::stop("Could not read BED file: " + path);
    }
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const int t = omp_thread();
      unsigned char* buf = packed.data() + static_cast<std::size_t>(t) * sel_stride;
      gather_packed(buffer.data() + k * stride, rows, buf);
      encode_column(*g, table, buf, g0 + k, missbuf.data() + static_cast<std::size_t>(t) * W, &fs);
    }
    Rcpp::checkUserInterrupt();
  }
  finish_stats(*g);
  return g;
}

// [[Rcpp::export]]
SEXP geno_open_pgen_cpp(const std::string& path, int n_total, const Rcpp::IntegerVector& allele_ct,
                        const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples,
                        bool mean, int threads, Rcpp::NumericVector fill_c, Rcpp::IntegerVector cell_ptr,
                        Rcpp::IntegerVector cell_i, Rcpp::NumericVector cell_g, Rcpp::NumericVector fill_af) {
  threads = resolve_threads(threads);
  const FillSpec fs = make_fill(fill_c, cell_ptr, cell_i, cell_g, fill_af, variants.size());
  PgenFile file(path, static_cast<std::size_t>(n_total), allele_ct, threads);
  Rcpp::XPtr<Geno> g(new Geno, true);
  const std::vector<int> rows = to_rows(samples);
  if (file.dosage()) {
    for (int v : variants)
      if (allele_ct[v] > 2)
        Rcpp::stop("Multiallelic variants are not supported in PGEN files with dosages: " + path);
    init_geno(*g, rows.size(), variants.size(), true);
    load_dosage(*g, file, variants, rows, threads, path, &fs);
    finish_stats(*g);
    return g;
  }
  init_geno(*g, rows.size(), variants.size(), mean);
  const NibbleTable& table = PgenFile::table();
  const std::size_t sel_stride = (rows.size() + 3) / 4, W = g->G.words;
  std::vector<unsigned char> packed(static_cast<std::size_t>(threads) * sel_stride);
  std::vector<u64> missbuf(static_cast<std::size_t>(threads) * W);
  const std::size_t p = g->p;
  for (std::size_t g0 = 0; g0 < p; g0 += kGroup) {
    const std::size_t count = std::min(kGroup, p - g0);
    int failed = 0;
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const int t = omp_thread();
      const unsigned char* column = file.read(static_cast<std::size_t>(variants[g0 + k]), t);
      if (!column) {
        #pragma omp atomic write
        failed = 1;
        continue;
      }
      unsigned char* buf = packed.data() + static_cast<std::size_t>(t) * sel_stride;
      gather_packed(column, rows, buf);
      encode_column(*g, table, buf, g0 + k, missbuf.data() + static_cast<std::size_t>(t) * W, &fs);
    }
    if (failed) Rcpp::stop("Could not read PGEN file: " + path);
    Rcpp::checkUserInterrupt();
  }
  finish_stats(*g);
  return g;
}

// [[Rcpp::export]]
Rcpp::List geno_info_cpp(SEXP ptr) {
  Rcpp::XPtr<Geno> g(ptr);
  return Rcpp::List::create(Rcpp::Named("n") = static_cast<double>(g->n),
                            Rcpp::Named("p") = static_cast<double>(g->p),
                            Rcpp::Named("mean") = g->mu,
                            Rcpp::Named("sd") = g->sd,
                            Rcpp::Named("n_miss") = g->nmiss,
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
    const int nt = omp_threads(), t = omp_thread();
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

// Read-only: the threads and memory budget a call with threads = 0 / mem_gb = 0 would use.
// [[Rcpp::export]]
Rcpp::List resources_cpp(int threads, double mem_gb) {
  const plink4r::Resources r = plink4r::detect_resources(threads, mem_gb);
  return Rcpp::List::create(Rcpp::Named("threads") = r.threads, Rcpp::Named("mem_gb") = r.mem_gb,
                            Rcpp::Named("threads_from") = r.threads_from, Rcpp::Named("mem_from") = r.mem_from);
}
