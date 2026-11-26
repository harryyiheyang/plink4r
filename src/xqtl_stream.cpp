// xqtl_stream_cpp: linear association of many continuous traits, each with the variants of its own window (cis-xQTL:
// expression, protein, methylation ...), streamed in blocks and written to a Parquet file as the scan goes; neither
// the genotype matrix nor a table of the pairs is formed in memory.
//
// The linear test of assoc_stream.cpp (S = 1) for q traits that share their samples and covariates. With V a basis of
// span(X), V'V = I (thin SVD of the covariates, singular values below sqrt(eps) x largest dropped), R = Y - V V'Y the
// residual traits (n x q) and g the centred genotypes of a variant (missing call = mean, i.e. 0 after centring):
//   A_j  = g'g - m'm,  m = V'g                 one number per variant, the same for every trait
//   c_jt = g'R_t                               R_t is orthogonal to V
//   beta_jt = c_jt / A_j,   se_jt = sqrt((R_t'R_t - c_jt^2 / A_j) / (n - k - 1) / A_j)
// which are those of lm(y_t ~ covariates + g). The genotype side (g'g, V'G) never sees a trait and the trait side (R)
// never sees a variant; they meet in one product G'R per block.
//
// Tables. The entry takes the tables as R holds them and does the ordering itself (place): the candidate variants are
// put in position order (chromosome code, then position; a file in that order is only checked), a trait's window
// start - window .. end + window becomes the range lo_t .. hi_t - 1 of those variants by two binary searches, traits
// without a variant in their window are left out, the others ordered by (lo, hi, name), and variants in no window
// dropped from the scan. The traits are read column by column from the table of R straight into R = the residual
// traits, rows taken in the order of the analysed samples (fill_residuals), so no matrix of the traits exists besides.
//
// Windows. Trait t is tested with the scan variants lo_t .. hi_t - 1, the traits ordered by lo. A thread's block of
// variants j0 .. j1 - 1 can only meet the traits first .. last - 1, last = #{t: lo_t < j1}, first = the first t whose
// running maximum of hi exceeds j0: a run of consecutive columns of R, used in place. The block computes G'R over
// that run and keeps the pairs with lo_t <= j < hi_t. A block that meets no trait is not read.
//
// Output. The pairs of a trait wait in memory (Pending: float beta and se at off_t + (j - lo_t), so threads write to
// disjoint places) until the scan has passed its last variant; the trait is then written as one row group of the
// Parquet file (parquet_write.h) and its pairs leave memory. Traits are written in their order, so what waits is the
// pairs of the traits whose windows the scan is inside. A row is one pair: Outcome (the trait), CHR, BP, REF, ALT
// (the counted allele, which BETA is per copy of), BETA, SE (float), N (non-missing calls of the variant). A variant
// without statistics (monomorphic or all missing; or in the span of the covariates) has no row. The file's key /
// value metadata hold what a P value needs besides BETA / SE (df = n - k - 1) and the samples and covariate rank.
//
// G is expanded to doubles kScanRowChunk rows at a time from the 2-bit codes (hard calls) or 0.01 dosages, so a
// thread holds b variants as n/4 (or 2n) bytes each plus kScanRowChunk x b doubles, and (k + active traits) x b
// products. A variant whose A lost more than 1 - 1/kKappaMax of g'g to the subtraction (strong collinearity with the
// covariates) is recomputed by explicit projection, and is in the span of the covariates when what the projection
// leaves of g is below kRankTol of its length, the tolerance that sets the rank of the covariates. The block size b
// follows from mem_gb (shared V, R; then per thread). Products: la().dgemm_tn (SIMD dispatch of la_kernel.cpp).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <limits>
#include "la_dispatch.h"
#include "parquet_write.h"
#include "scan_block.h"
#include "resources.h"

namespace {

using plink4r::kScanRowChunk;
constexpr double kKappaMax = 1e4;        // explicit projection when g'g / A exceeds this
constexpr double kRankTol = 1.4901161193847656e-08;  // sqrt(DBL_EPSILON), as mgcv: rank of X, and g in span(X)

inline double dot(const double* a, const double* b, int len) {
  double s = 0.0;
  for (int i = 0; i < len; ++i) s += a[i] * b[i];
  return s;
}

// What all threads share: the covariate basis V (n x k), the residual traits R (n x q) and rr[t] = R_t'R_t.
struct Pre {
  int n = 0, k = 0, q = 0;
  std::vector<double> V, R, rr;
};

// V: an orthonormal basis of the covariates (thin SVD; singular values below kRankTol x the largest dropped).
void build_basis(Pre& pre, const Rcpp::NumericMatrix& X) {
  const int n = X.nrow(), kX = X.ncol();
  const std::size_t N = n;
  pre.n = n;
  Eigen::MatrixXd A(n, kX);
  std::copy(X.begin(), X.end(), A.data());
  Eigen::BDCSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU);
  const Eigen::VectorXd& sig = svd.singularValues();
  int k = 0;
  while (k < sig.size() && sig[k] > 0 && sig[k] >= kRankTol * sig[0]) ++k;  // min(n, kX) singular values
  if (k == 0) Rcpp::stop("The covariate matrix is zero.");
  pre.k = k;
  pre.V.assign(svd.matrixU().data(), svd.matrixU().data() + N * k);
}

// A column of the table of traits as R holds it: doubles or integers.
struct Column {
  const double* d = nullptr;
  const int* i = nullptr;
};

// R and rr: column t of R is trait cols[t] over the analysed samples (row[i] = its row in the table for sample i),
// then made orthogonal to V. missing[t] = 1 when the trait has a missing or non-finite value among them.
void fill_residuals(Pre& pre, const std::vector<Column>& cols, const int* row, int T, std::vector<char>& missing) {
  const int n = pre.n, k = pre.k, q = static_cast<int>(cols.size());
  const std::size_t N = n;
  pre.q = q;
  pre.R.assign(N * q, 0.0);
  pre.rr.assign(q, 0.0);
  missing.assign(q, 0);
  const double* V = pre.V.data();
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
  {
    std::vector<double> coef(k);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int t = 0; t < q; ++t) {
      double* r = &pre.R[t * N];
      bool ok = true;
      if (cols[t].d) {
        for (int i = 0; i < n; ++i) { r[i] = cols[t].d[row[i]]; ok = ok && std::isfinite(r[i]); }
      } else {
        for (int i = 0; i < n; ++i) {
          const int v = cols[t].i[row[i]];
          ok = ok && v != NA_INTEGER;
          r[i] = v;
        }
      }
      if (!ok) { missing[t] = 1; continue; }
      for (int pass = 0; pass < 2; ++pass) {  // the second pass removes what rounding left of V in r
        for (int c = 0; c < k; ++c) coef[c] = dot(V + c * N, r, n);
        for (int c = 0; c < k; ++c) {
          const double m = coef[c];
          const double* v = V + c * N;
          for (std::size_t i = 0; i < N; ++i) r[i] -= v[i] * m;
        }
      }
      pre.rr[t] = dot(r, r, n);
    }
  }
}

// The variant ranges of the traits (ordered by lo) and the place of every pair among the pairs of all traits.
struct Windows {
  int q = 0;
  std::vector<int> lo, hi, top;    // top[t] = max(hi[0..t])
  std::vector<std::size_t> off;    // off[t] = pairs of the traits before t; off[q] = all pairs
  Windows(const std::vector<int>& lo_, const std::vector<int>& hi_)
      : q(static_cast<int>(lo_.size())), lo(lo_), hi(hi_), top(lo_.size()), off(lo_.size() + 1, 0) {
    int running = 0;
    for (int t = 0; t < q; ++t) {
      running = std::max(running, hi[t]);
      top[t] = running;
      off[t + 1] = off[t] + static_cast<std::size_t>(hi[t] - lo[t]);
    }
  }
  // Traits first .. last - 1: every trait with a variant among j0 .. j1 - 1 is one of them.
  void active(long long j0, long long j1, int& first, int& last) const {
    last = static_cast<int>(std::lower_bound(lo.begin(), lo.end(), j1) - lo.begin());
    first = static_cast<int>(std::upper_bound(top.begin(), top.end(), j0) - top.begin());
    if (first > last) first = last;
  }
  // Most traits a block of b variants meets (blocks start at multiples of b).
  int most_active(int b, int p) const {
    int most = 0, first, last;
    for (long long j0 = 0; j0 < p; j0 += b) {
      active(j0, std::min<long long>(j0 + b, p), first, last);
      most = std::max(most, last - first);
    }
    return most;
  }
};

// The scan as the tables give it: the variants of at least one window in position order, and the traits with a
// variant in their window, ordered by (lo, hi, name), with their ranges of those variants.
struct Scan {
  std::vector<int> row;                          // file row of scan variant j
  std::vector<std::int32_t> chr, bp, ref, alt;   // its chromosome code and position; alleles as codes into alleles
  std::vector<std::string> alleles;
  std::vector<int> trait, lo, hi;                // trait: the place of trait t in the tables given
};

// A place on the genome as one number: chromosome codes and positions are below 2^31.
inline std::int64_t genome_key(std::int64_t chr, std::int64_t bp) { return chr * 4294967296LL + bp; }

// vt: every variant line of the file (variant_text.h); those with two alleles, a chromosome code and a position are
// the candidates. name, o_chr, o_start, o_end: the traits; a trait takes the variants of its chromosome at
// start - window to end + window, both ends included.
Scan place(const plink4r::VariantTable& vt, const std::vector<std::string>& name, const Rcpp::IntegerVector& o_chr,
           const Rcpp::NumericVector& o_start, const Rcpp::NumericVector& o_end, double window) {
  // the candidates in position order: (place, file row); equal places keep the order given
  using Cand = std::pair<std::int64_t, int>;
  std::vector<Cand> cand;
  const int m = static_cast<int>(vt.size());
  cand.reserve(m);
  for (int r = 0; r < m; ++r) {
    if (!vt.biallelic(r) || vt.chr[r] < 0 || vt.bp[r] < 0) continue;
    cand.emplace_back(genome_key(vt.chr[r], vt.bp[r]), r);
  }
  const auto by_place = [](const Cand& a, const Cand& b) { return a.first < b.first; };
  if (!std::is_sorted(cand.begin(), cand.end(), by_place)) std::stable_sort(cand.begin(), cand.end(), by_place);
  const int pc = static_cast<int>(cand.size());
  std::vector<std::int64_t> key(pc);
  for (int j = 0; j < pc; ++j) key[j] = cand[j].first;

  // the range of every trait among the candidates
  const double kPosMax = 2147483647.0;
  const int q_in = static_cast<int>(name.size());
  std::vector<int> lo(q_in, 0), hi(q_in, 0), kept;
  for (int t = 0; t < q_in; ++t) {
    if (o_chr[t] == NA_INTEGER || o_chr[t] < 0 || !std::isfinite(o_start[t]) || !std::isfinite(o_end[t])) continue;
    const double a = std::max(0.0, std::ceil(o_start[t] - window));
    const double b = std::min(kPosMax, std::floor(o_end[t] + window));
    if (!(a <= b)) continue;
    const std::int64_t ka = genome_key(o_chr[t], static_cast<std::int64_t>(a));
    const std::int64_t kb = genome_key(o_chr[t], static_cast<std::int64_t>(b));
    lo[t] = static_cast<int>(std::lower_bound(key.begin(), key.end(), ka) - key.begin());
    hi[t] = static_cast<int>(std::upper_bound(key.begin(), key.end(), kb) - key.begin());
    if (hi[t] > lo[t]) kept.push_back(t);
  }
  std::sort(kept.begin(), kept.end(), [&](int a, int b) {
    if (lo[a] != lo[b]) return lo[a] < lo[b];
    if (hi[a] != hi[b]) return hi[a] < hi[b];
    return name[a] < name[b];
  });

  // the candidates of at least one window, renumbered
  std::vector<int> before(pc + 1, 0);
  for (int t : kept) { ++before[lo[t]]; --before[hi[t]]; }
  Scan sc;
  int inside = 0, count = 0;
  for (int j = 0; j < pc; ++j) {
    inside += before[j];
    before[j] = count;
    if (inside > 0) { sc.row.push_back(cand[j].second); ++count; }
  }
  before[pc] = count;
  sc.trait = kept;
  for (int t : kept) { sc.lo.push_back(before[lo[t]]); sc.hi.push_back(before[hi[t]]); }

  // what the output says of a variant; the alleles as codes into their distinct strings (counted alleles first)
  sc.chr.resize(count); sc.bp.resize(count); sc.ref.resize(count); sc.alt.resize(count);
  std::vector<std::int32_t> code(vt.alleles.size(), -1);
  const auto coded = [&](std::int32_t a) {
    if (code[a] < 0) {
      code[a] = static_cast<std::int32_t>(sc.alleles.size());
      sc.alleles.push_back(vt.alleles[a]);
    }
    return code[a];
  };
  for (int j = 0; j < count; ++j) {
    sc.chr[j] = vt.chr[sc.row[j]]; sc.bp[j] = vt.bp[sc.row[j]];
    sc.alt[j] = coded(vt.alt[sc.row[j]]);
  }
  for (int j = 0; j < count; ++j) sc.ref[j] = coded(vt.ref[sc.row[j]]);
  return sc;
}

// Variants of one thread: their compact genotypes (scan_block.h) and the block products.
struct Worker : plink4r::ScanBlock {
  std::vector<double> G, M, C, gg, z, zt;
  void init(const Pre& pre, int capacity, bool dos, int most_active) {
    ScanBlock::init(capacity, pre.n, dos);
    const std::size_t R = std::min(kScanRowChunk, n);
    G.assign(R * cap, 0.0);
    M.assign(static_cast<std::size_t>(pre.k) * cap, 0.0);
    C.assign(static_cast<std::size_t>(most_active) * cap, 0.0);
    gg.assign(cap, 0.0);
  }
};

// g'g (gg), V'G (M, k x nb) and G'R over traits first .. last - 1 (C, (last - first) x nb) of the nb slots.
void block_products(const Pre& pre, Worker& wk, int nb, int first, int last) {
  const int n = pre.n, k = pre.k, qa = last - first;
  const std::size_t N = n;
  const plink4r::LaFns& LA = plink4r::la();
  if (wk.C.size() < static_cast<std::size_t>(qa) * nb) wk.C.resize(static_cast<std::size_t>(qa) * nb);
  std::fill(wk.gg.begin(), wk.gg.begin() + nb, 0.0);
  for (int row0 = 0; row0 < n; row0 += kScanRowChunk) {
    const int rows = std::min(kScanRowChunk, n - row0);
    for (int s = 0; s < nb; ++s) {
      double* g = &wk.G[static_cast<std::size_t>(s) * rows];
      wk.centred(s, row0, rows, g);
      wk.gg[s] += dot(g, g, rows);
    }
    const int beta = row0 > 0;
    LA.dgemm_tn(&pre.V[row0], rows, k, n, wk.G.data(), nb, rows, wk.M.data(), k, beta);
    LA.dgemm_tn(&pre.R[first * N + row0], rows, qa, n, wk.G.data(), nb, rows, wk.C.data(), qa, beta);
  }
}

// A of slot s by explicit projection; wk.zt is left holding g - V V'g.
double explicit_projection(const Pre& pre, Worker& wk, int s) {
  const int n = pre.n, k = pre.k;
  const std::size_t N = n;
  if (wk.z.empty()) { wk.z.resize(N); wk.zt.resize(N); }
  wk.centred(s, 0, n, wk.z.data());
  wk.zt = wk.z;
  for (int c = 0; c < k; ++c) {
    const double* v = &pre.V[c * N];
    const double m = dot(v, wk.z.data(), n);
    for (std::size_t i = 0; i < N; ++i) wk.zt[i] -= v[i] * m;
  }
  return dot(wk.zt.data(), wk.zt.data(), n);
}

// The pairs of the traits not yet written: beta and se of traits first .., trait by trait.
struct Pending {
  int first = 0;
  std::vector<float> beta, se;
};

// The variants of the scan as the output names them, and the columns of one row group.
struct Rows {
  const std::int32_t* chr; const std::int32_t* bp; const std::int32_t* ref; const std::int32_t* alt;
  const double* n_obs;
  const std::vector<std::string>* alleles;
  std::vector<std::string> outcome;
  std::vector<std::int32_t> c_chr, c_bp, c_ref, c_alt, c_n;
  std::vector<float> c_beta, c_se;
  // Writes trait t (its pairs start at b / s) as one row group: the pairs with statistics. Returns their number.
  std::size_t write(plink4r::ParquetWriter& W, const Windows& win, int t, const float* b, const float* s) {
    c_chr.clear(); c_bp.clear(); c_ref.clear(); c_alt.clear(); c_n.clear(); c_beta.clear(); c_se.clear();
    for (int j = win.lo[t]; j < win.hi[t]; ++j) {
      const int i = j - win.lo[t];
      if (std::isnan(b[i])) continue;
      c_chr.push_back(chr[j]); c_bp.push_back(bp[j]); c_ref.push_back(ref[j]); c_alt.push_back(alt[j]);
      c_n.push_back(static_cast<std::int32_t>(n_obs[j])); c_beta.push_back(b[i]); c_se.push_back(s[i]);
    }
    if (c_bp.empty()) return 0;
    W.begin_group(c_bp.size());
    W.string_const(outcome[t]);
    W.int32_dict(c_chr.data());
    W.int32_delta(c_bp.data());
    W.string_dict(c_ref.data(), *alleles);
    W.string_dict(c_alt.data(), *alleles);
    W.float_plain(c_beta.data());
    W.float_plain(c_se.data());
    W.int32_dict(c_n.data());
    W.end_group();
    return c_bp.size();
  }
};

}  // namespace

// path: .pgen or .bed; bed = true for a BED (A1 = .bim column 5 counted), else PGEN (ALT counted). variant_text: its
// .bim / .pvar, read here (variant_text.h): chromosome, position and the two alleles of every line; the lines with
// two alleles may be scanned. sample_idx: 0-based file rows of the analysed samples. X: n x kX covariates with the
// intercept, rows as sample_idx. table: the table of traits as R holds it (a list of columns); row: for every
// analysed sample its 0-based row in table; col: for every trait its 0-based column in table (double or integer, no
// missing value among the analysed samples). outcome, o_chr, o_start, o_end: name, chromosome code and the two ends
// of every trait, in any order; window: base pairs added on both sides. out: the Parquet file to write. created_by:
// "<application> version <x.y.z>" for the file. block: variants per thread (0 = from mem_gb). Returns the rows
// written, the pairs asked for, which traits were tested (those with a variant in their window), the variants
// scanned, those without statistics (monomorphic: no variation or all missing; singular: in the span of the
// covariates), chromosome code and position of the lines with more than two alleles (not scanned) and, per block
// of the scan, the number of traits it met (active).
// [[Rcpp::export]]
Rcpp::List xqtl_stream_cpp(const std::string& path, bool bed, int n_raw, const std::string& variant_text,
                           Rcpp::IntegerVector sample_idx, Rcpp::NumericMatrix X, Rcpp::List table,
                           Rcpp::IntegerVector row, Rcpp::IntegerVector col, Rcpp::CharacterVector outcome,
                           Rcpp::IntegerVector o_chr, Rcpp::NumericVector o_start, Rcpp::NumericVector o_end,
                           double window, const std::string& out, const std::string& created_by, int threads,
                           double mem_gb, int block) {
  const int T = plink4r::resolve_threads(threads);
  mem_gb = plink4r::resolve_mem_gb(mem_gb, T);
  const int ns = sample_idx.size(), q_in = col.size();
  if (X.nrow() != ns || row.size() != ns) Rcpp::stop("X / row do not match the samples");
  if (outcome.size() != q_in || o_chr.size() != q_in || o_start.size() != q_in || o_end.size() != q_in)
    Rcpp::stop("outcome / chr / start / end do not match the traits");
  if (!(window >= 0)) Rcpp::stop("window must not be negative");

  // ---- the scan: variants in position order, traits with their ranges
  const std::vector<std::string> name = Rcpp::as<std::vector<std::string>>(outcome);
  Scan sc;
  std::vector<int> allele_ct, skip_chr, skip_bp;
  std::size_t m = 0;
  {
    plink4r::VariantTable vt = plink4r::read_variant_table(variant_text, bed);
    m = vt.size();
    bool any = false;
    for (std::size_t r = 0; r < m; ++r) {
      if (vt.biallelic(r)) { any = true; continue; }
      skip_chr.push_back(vt.chr[r] < 0 ? NA_INTEGER : vt.chr[r]);
      skip_bp.push_back(vt.bp[r] < 0 ? NA_INTEGER : vt.bp[r]);
    }
    if (!any) Rcpp::stop("No biallelic variant to test.");
    sc = place(vt, name, o_chr, o_start, o_end, window);
    allele_ct = std::move(vt.allele_ct);
  }
  const int p = static_cast<int>(sc.row.size()), q = static_cast<int>(sc.trait.size());
  if (q == 0) Rcpp::stop("No outcome has a variant in its window.");
  const Windows win(sc.lo, sc.hi);
  plink4r::ScanFile F(path, bed, n_raw, m, allele_ct, sc.row, sample_idx, T);
  const bool dos = F.dosage;

  // ---- what the threads share: the covariate basis, and the traits of the scan read from the table
  Pre pre;
  build_basis(pre, X);
  const int k = pre.k;
  const double df = static_cast<double>(ns) - k - 1;
  if (df < 1) Rcpp::stop("No residual degree of freedom: %d samples, covariate rank %d.", ns, k);
  Rows rows;
  {
    std::vector<Column> cols(q);
    R_xlen_t len = -1;
    for (int t = 0; t < q; ++t) {
      const int c = col[sc.trait[t]];
      rows.outcome.push_back(name[sc.trait[t]]);
      if (c < 0 || c >= table.size()) Rcpp::stop("col is not within the columns of the table");
      const SEXP x = VECTOR_ELT(table, c);
      if (TYPEOF(x) == REALSXP) cols[t].d = REAL(x);
      else if (TYPEOF(x) == INTSXP) cols[t].i = INTEGER(x);
      else Rcpp::stop("Not numeric in outcome: " + rows.outcome.back());
      if (len < 0) len = Rf_xlength(x);
      if (Rf_xlength(x) != len) Rcpp::stop("The columns of the table of outcomes differ in length");
    }
    for (int i = 0; i < ns; ++i)
      if (row[i] < 0 || row[i] >= len) Rcpp::stop("row is not within the rows of the table");
    std::vector<char> missing;
    fill_residuals(pre, cols, row.begin(), T, missing);
    const int kNamed = 5;  // the outcomes with a missing value are named, up to kNamed of them
    std::string names;
    int bad = 0;
    for (int t = 0; t < q; ++t)
      if (missing[t] && bad++ < kNamed) names += (names.empty() ? "" : ", ") + rows.outcome[t];
    if (bad > kNamed) names += ", ... (" + std::to_string(bad) + " in all)";
    if (bad) Rcpp::stop("Missing or not finite among the analysed samples, in outcome: " + names);
  }

  // block size from the memory budget
  const std::size_t N = ns;
  const int b_most = block > 0 ? block : plink4r::kScanBlockMax;
  // shared: V and R, per variant its five output numbers and two counts, X and its copy in the SVD
  const double shared = 8.0 * N * (k + q) + 36.0 * p + 8.0 * N * X.ncol() * 2;
  const double per_var = plink4r::ScanBlock::var_bytes(ns, dos) + 8.0 * std::min<int>(kScanRowChunk, ns) +
                         8.0 * (k + win.most_active(b_most, p) + 1);
  const double per_thread = 8.0 * N * 2 + F.thread_bytes();
  const int b = plink4r::scan_block_size(block, mem_gb, shared, per_thread, per_var, p, ns, T);
  const int most = win.most_active(b, p);

  // the output file
  rows.chr = sc.chr.data(); rows.bp = sc.bp.data(); rows.ref = sc.ref.data(); rows.alt = sc.alt.data();
  rows.alleles = &sc.alleles;
  using PW = plink4r::ParquetWriter;
  PW writer(out,
            {{"Outcome", PW::kString}, {"CHR", PW::kInt32}, {"BP", PW::kInt32}, {"REF", PW::kString},
             {"ALT", PW::kString}, {"BETA", PW::kFloat}, {"SE", PW::kFloat}, {"N", PW::kInt32}},
            {{"df", std::to_string(static_cast<long long>(df))}, {"samples", std::to_string(ns)},
             {"covariate_rank", std::to_string(k)}},
            created_by);

  std::vector<double> nobs(p, 0.0), asum(p, 0.0);
  rows.n_obs = nobs.data();
  Rcpp::IntegerVector active((p + b - 1) / b);
  std::vector<Worker> wk(T);
  for (int t = 0; t < T; ++t) wk[t].init(pre, b, dos, most);
  std::vector<int> bad(T, -1);
  Pending pend;
  const float kNone = std::numeric_limits<float>::quiet_NaN();
  int nexp = 0, n_mono = 0, n_sing = 0;
  double products = 0.0;
  std::size_t written = 0, waiting_most = 0;

  for (long long c0 = 0; c0 < p; c0 += static_cast<long long>(b) * T) {
    const long long c1 = std::min<long long>(p, c0 + static_cast<long long>(b) * T);
    // room for the pairs of every trait this round of blocks can meet
    const int reach = static_cast<int>(std::lower_bound(win.lo.begin(), win.lo.end(), c1) - win.lo.begin());
    const std::size_t base = win.off[pend.first];
    if (reach > pend.first && pend.beta.size() < win.off[reach] - base) {
      pend.beta.resize(win.off[reach] - base, kNone);
      pend.se.resize(win.off[reach] - base, kNone);
    }
    waiting_most = std::max(waiting_most, pend.beta.size());
    float* bt = pend.beta.data();
    float* st = pend.se.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(T) reduction(+ : nexp, n_mono, n_sing, products)
#endif
    for (int t = 0; t < T; ++t) {
      const int th = plink4r::thread_id();
      Worker& W = wk[th];
      const long long j0 = c0 + static_cast<long long>(t) * b;
      if (j0 >= p) continue;
      const int nb = static_cast<int>(std::min<long long>(b, p - j0));
      int first, last;
      win.active(j0, j0 + nb, first, last);
      active[j0 / b] = last - first;
      if (first == last) continue;
      const long long failed = W.read(F, j0, nb, th, &nobs[j0], &asum[j0]);
      if (failed >= 0) bad[th] = static_cast<int>(failed);
      block_products(pre, W, nb, first, last);
      const int qa = last - first;
      products += static_cast<double>(qa) * nb;
      for (int s = 0; s < nb; ++s) {
        const int j = static_cast<int>(j0 + s);
        if (W.mono[s]) { ++n_mono; continue; }
        const double* m = &W.M[static_cast<std::size_t>(s) * k];
        const double* c = &W.C[static_cast<std::size_t>(s) * qa];
        double A = W.gg[s] - dot(m, m, k);
        const bool ill = !(A > W.gg[s] / kKappaMax);
        if (ill) { A = explicit_projection(pre, W, s); ++nexp; }
        if (!(A > kRankTol * kRankTol * W.gg[s])) { ++n_sing; continue; }  // g is in the span of the covariates
        for (int u = first; u < last; ++u) {
          if (j < win.lo[u] || j >= win.hi[u]) continue;
          const double cu = ill ? dot(W.zt.data(), &pre.R[u * N], ns) : c[u - first];
          const std::size_t o = win.off[u] - base + static_cast<std::size_t>(j - win.lo[u]);
          bt[o] = static_cast<float>(cu / A);
          st[o] = static_cast<float>(std::sqrt(std::max(0.0, pre.rr[u] - cu * cu / A) / df / A));
        }
      }
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) Rcpp::stop("Genotype read failed at variant %d of %s", bad[t] + 1, path.c_str());
    // the traits whose last variant the scan has passed are written and leave memory
    const int done = std::min(q, static_cast<int>(std::upper_bound(win.top.begin(), win.top.end(),
                                                                   static_cast<int>(c1)) - win.top.begin()));
    if (done > pend.first) {
      for (int t = pend.first; t < done; ++t)
        written += rows.write(writer, win, t, bt + (win.off[t] - base), st + (win.off[t] - base));
      const std::size_t gone = std::min(pend.beta.size(), win.off[done] - base);
      pend.beta.erase(pend.beta.begin(), pend.beta.begin() + gone);
      pend.se.erase(pend.se.begin(), pend.se.begin() + gone);
      pend.first = done;
    }
    Rcpp::checkUserInterrupt();
  }
  writer.close();

  Rcpp::List res;
  res["rows"] = static_cast<double>(written); res["pairs"] = static_cast<double>(win.off[q]);
  res["bytes"] = static_cast<double>(writer.bytes()); res["row_groups"] = static_cast<double>(writer.groups());
  res["monomorphic"] = n_mono; res["singular"] = n_sing;
  res["waiting_most"] = static_cast<double>(waiting_most);
  res["df"] = df; res["rank"] = k; res["explicit"] = nexp; res["block"] = b; res["active"] = active;
  res["products"] = products; res["threads"] = T; res["mem_gb"] = mem_gb;
  res["dosage"] = dos; res["simd"] = std::string(plink4r::la().isa);
  Rcpp::LogicalVector tested(q_in, false);
  for (int t : sc.trait) tested[t] = true;
  res["tested"] = tested; res["outcomes"] = q; res["variants"] = p;
  res["skipped_chr"] = Rcpp::wrap(skip_chr); res["skipped_bp"] = Rcpp::wrap(skip_bp);
  return res;
}
