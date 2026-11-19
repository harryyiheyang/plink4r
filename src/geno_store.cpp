// The geno store (see geno_store.h): opening a BED / PGEN, appending columns, row / column subsets with exact
// flips, the fill writers and the views of the missing cells.
#include "geno_store.h"
#include "geno_access.h"

using namespace plink4r::gstore;

namespace {

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

}  // namespace

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

// [[Rcpp::export]]
Rcpp::List geno_info_cpp(SEXP ptr) {
  Rcpp::XPtr<Geno> g(ptr);
  std::vector<double> os(g->p);
  Rcpp::NumericVector af(g->p, NA_REAL);  // observed frequency of the counted allele (NA: nothing observed)
  for (std::size_t j = 0; j < g->p; ++j) {
    os[j] = obs_sum_col(*g, j);
    const double obs = static_cast<double>(g->n) - g->nmiss[j];
    if (obs > 0) af[j] = os[j] / (2.0 * obs);
  }
  return Rcpp::List::create(Rcpp::Named("n") = static_cast<double>(g->n),
                            Rcpp::Named("p") = static_cast<double>(g->p),
                            Rcpp::Named("mean") = g->mu,
                            Rcpp::Named("sd") = g->sd,
                            Rcpp::Named("n_miss") = g->nmiss,
                            Rcpp::Named("obs_sum") = os,
                            Rcpp::Named("af") = af,
                            Rcpp::Named("dosage") = g->dosage,
                            Rcpp::Named("dosage_nnz") = static_cast<double>(g->dnnz()),
                            Rcpp::Named("bytes_H") = bytes_h(*g),
                            Rcpp::Named("bytes_D") = bytes_d(*g));
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

// A new geno merged from the objects ptrs. cols(k, c) / rows(k, r) (0-based, -1 = none) give the column of
// object k that is output column c and the row of object k that is output row r. Cell (r, c) is copied as
// stored (H, D, and whether it is a missing cell) from the first object k with rows(k, r) >= 0 and
// cols(k, c) >= 0. A cell no object holds becomes a missing cell of the result, filled with the observed mean
// of its column (over the cells copied that are not missing cells; 0 when there is none). The objects of a
// column must share the unit of D. A column held by one object in the output row order (rows(k, r) = r, all
// its rows) is copied word by word. Returns list(ptr, uncovered = cells no object holds, per column).
// [[Rcpp::export]]
Rcpp::List geno_merge_cpp(Rcpp::List ptrs, const Rcpp::IntegerMatrix& cols, const Rcpp::IntegerMatrix& rows,
                          int threads) {
  threads = resolve_threads(threads);
  const std::size_t K = ptrs.size(), p = cols.ncol(), n = rows.ncol();
  if (!K || static_cast<std::size_t>(cols.nrow()) != K || static_cast<std::size_t>(rows.nrow()) != K)
    Rcpp::stop("cols and rows must have one row per object.");
  const int* cm = cols.begin();
  const int* rm = rows.begin();
  auto col = [cm, K](std::size_t k, std::size_t c) { return cm[k + c * K]; };
  auto row = [rm, K](std::size_t k, std::size_t r) { return rm[k + r * K]; };
  std::vector<const Geno*> a(K);
  for (std::size_t k = 0; k < K; ++k) a[k] = Rcpp::XPtr<Geno>(static_cast<SEXP>(ptrs[k])).get();
  std::vector<char> ident(K, 1);  // object k holds exactly the output rows, in their order
  bool dosage = false;
  for (std::size_t k = 0; k < K; ++k) {
    dosage = dosage || a[k]->dosage;
    if (a[k]->n != n) ident[k] = 0;
    for (std::size_t r = 0; r < n; ++r) {
      const int i = row(k, r);
      if (i >= 0 && static_cast<std::size_t>(i) >= a[k]->n) Rcpp::stop("rows out of range.");
      if (i != static_cast<int>(r)) ident[k] = 0;
    }
  }
  for (std::size_t c = 0; c < p; ++c) {
    double unit = 0;
    bool first = true;
    for (std::size_t k = 0; k < K; ++k) {
      const int j = col(k, c);
      if (j < 0) continue;
      if (static_cast<std::size_t>(j) >= a[k]->p) Rcpp::stop("cols out of range.");
      const double u = a[k]->D[j].unit;
      if (!first && u != unit) Rcpp::stop("The objects hold column %d in different units (hard calls / dosages).", static_cast<int>(c) + 1);
      unit = u;
      first = false;
    }
    if (first) Rcpp::stop("A column comes from no object.");
  }
  Rcpp::XPtr<Geno> g(new Geno, true);
  init_geno(*g, n, p, a[0]->mean);
  g->dosage = dosage;
  const std::size_t W = g->G.words;
  std::vector<std::size_t> base(K + 1, 0);  // start of object k in the per-thread row buffers
  for (std::size_t k = 0; k < K; ++k) base[k + 1] = base[k] + a[k]->n;
  Rcpp::NumericVector unc(p);
  double* uncp = unc.begin();
  #pragma omp parallel num_threads(threads)
  {
    std::vector<std::int16_t> qbuf, qo;
    std::vector<unsigned char> mbuf, ho, st;  // st: 0 copied, 1 copied missing cell, 2 held by no object
    std::vector<std::size_t> ks;
    #pragma omp for schedule(dynamic, 16)
    for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(p); ++cc) {
      const std::size_t c = static_cast<std::size_t>(cc);
      u64* lo = g->G.bits.data() + 2 * c * W;
      u64* hi = lo + W;
      std::vector<std::uint32_t>& Mc = g->M[c];
      ks.clear();
      for (std::size_t k = 0; k < K; ++k) if (col(k, c) >= 0) ks.push_back(k);
      g->D[c].unit = a[ks[0]]->D[col(ks[0], c)].unit;
      DBuilder db(*g, c);
      if (ident[ks[0]]) {  // the first object holds every row in order: its column as it is
        const Geno& s = *a[ks[0]];
        const std::size_t j = static_cast<std::size_t>(col(ks[0], c));
        const DCol& dc = s.D[j];
        std::copy(s.G.lo(j), s.G.lo(j) + (n + 63) / 64, lo);
        std::copy(s.G.hi(j), s.G.hi(j) + (n + 63) / 64, hi);
        g->D[c].q.reserve(dc.q.size());
        std::size_t e = 0;
        for (std::size_t w = 0; w < dc.mask.size(); ++w)
          for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) {
            const std::size_t i = 64 * w + __builtin_ctzll(bits);
            db.add(i, h_at(s.G, j, i), dc.q[e++]);
          }
        Mc = s.M[j];
      } else {
        if (qbuf.empty()) { qbuf.assign(base[K], 0); mbuf.assign(base[K], 0); qo.resize(n); ho.resize(n); st.resize(n); }
        for (std::size_t k : ks) {
          const DCol& dc = a[k]->D[col(k, c)];
          std::size_t e = 0;
          for (std::size_t w = 0; w < dc.mask.size(); ++w)
            for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) qbuf[base[k] + 64 * w + __builtin_ctzll(bits)] = dc.q[e++];
          for (std::uint32_t i : a[k]->M[col(k, c)]) mbuf[base[k] + i] = 1;
        }
        const double unit = g->D[c].unit;
        const int us = unit_steps(unit);  // observed sum exactly in int16 steps (an appended column: in double)
        std::int64_t ss = 0;
        double so = 0.0;
        std::size_t no = 0, nu = 0;
        for (std::size_t r = 0; r < n; ++r) {
          st[r] = 2;
          for (std::size_t k : ks) {
            const int i = row(k, r);
            if (i < 0) continue;
            ho[r] = static_cast<unsigned char>(h_at(a[k]->G, static_cast<std::size_t>(col(k, c)), static_cast<std::size_t>(i)));
            qo[r] = qbuf[base[k] + i];
            st[r] = mbuf[base[k] + i];
            if (!st[r]) {
              if (us) ss += us * static_cast<std::int64_t>(ho[r]) + qo[r];
              else so += ho[r] + qo[r] * unit;
              ++no;
            }
            break;
          }
          if (st[r] == 2) ++nu;
        }
        unsigned int hf = 0;
        int qf = 0;
        if (us) so = static_cast<double>(ss) / us;
        if (nu) quantize_fill(no ? so / static_cast<double>(no) : 0.0, unit, hf, qf);
        uncp[c] = static_cast<double>(nu);
        for (std::size_t r = 0; r < n; ++r) {
          const unsigned int h = st[r] == 2 ? hf : ho[r];
          const int q = st[r] == 2 ? qf : qo[r];
          if (h == 1u) lo[r >> 6] |= u64(1) << (r & 63u);
          else if (h == 2u) hi[r >> 6] |= u64(1) << (r & 63u);
          if (q) db.add(r, h, q);
          if (st[r]) Mc.push_back(static_cast<std::uint32_t>(r));
        }
        for (std::size_t k : ks) {
          const DCol& dc = a[k]->D[col(k, c)];
          for (std::size_t w = 0; w < dc.mask.size(); ++w)
            for (u64 bits = dc.mask[w]; bits; bits &= bits - 1) qbuf[base[k] + 64 * w + __builtin_ctzll(bits)] = 0;
          for (std::uint32_t i : a[k]->M[col(k, c)]) mbuf[base[k] + i] = 0;
        }
      }
      db.finish();
      plane_stats(*g, c);
      g->nmiss[c] = static_cast<double>(Mc.size());
    }
  }
  finish_stats(*g);
  return Rcpp::List::create(Rcpp::Named("ptr") = g, Rcpp::Named("uncovered") = unc);
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

// Per-column fills: every missing cell (rows M[j]) of column j takes fill_c[j] when finite, else keeps its
// current value. A value f is stored as H = round(f) (halves up, within 0..2) and D = f - H in 1/16384.
// Cells not missing are never touched. inplace = TRUE changes ptr itself (only for an object nobody else
// references yet); otherwise a modified copy is returned and ptr is unchanged.
// [[Rcpp::export]]
SEXP geno_fill_cpp(SEXP ptr, const Rcpp::NumericVector& fill_c, bool inplace, int threads) {
  Rcpp::XPtr<Geno> src(ptr);
  threads = resolve_threads(threads);
  const std::size_t p = src->p;
  if (static_cast<std::size_t>(fill_c.size()) != p) Rcpp::stop("fill values must have one value per variant.");
  for (std::size_t j = 0; j < p; ++j)
    if (!std::isnan(fill_c[j]) && !fill_value_ok(fill_c[j])) Rcpp::stop("fill values must be finite and within (-2, 4).");
  Rcpp::XPtr<Geno> out = inplace ? src : Rcpp::XPtr<Geno>(new Geno(*src), true);
  Geno& g = *out;
  const double* fc = fill_c.begin();
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 16)
  for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(p); ++jj) {
    const std::size_t j = static_cast<std::size_t>(jj);
    if (std::isfinite(fc[j])) fill_column(g, j, true, fc[j], static_cast<const int*>(nullptr), nullptr, 0, 0);
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

// Read-only: the threads and memory budget a call with threads = 0 / mem_gb = 0 would use.
// [[Rcpp::export]]
Rcpp::List resources_cpp(int threads, double mem_gb) {
  const plink4r::Resources r = plink4r::detect_resources(threads, mem_gb);
  return Rcpp::List::create(Rcpp::Named("threads") = r.threads, Rcpp::Named("mem_gb") = r.mem_gb,
                            Rcpp::Named("threads_from") = r.threads_from, Rcpp::Named("mem_from") = r.mem_from);
}
