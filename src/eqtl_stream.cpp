// eqtl_stream_cpp: linear association of many continuous traits, each with the variants of its own window (cis-eQTL),
// streamed in blocks; neither the genotype matrix nor a variants x traits table of products is formed.
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
// Windows. Trait t is tested with the scan variants lo_t .. hi_t - 1 (a window around an anchor is such a range once
// the variants are in position order; this file knows ranges only). The traits come ordered by lo. A thread's block of
// variants j0 .. j1 - 1 can only meet the traits first .. last - 1, last = #{t: lo_t < j1}, first = the first t whose
// running maximum of hi exceeds j0: a run of consecutive columns of R, used in place. The block computes G'R over
// that run and writes the pairs with lo_t <= j < hi_t. The pairs of trait t are stored together, at off_t + (j - lo_t)
// with off_t = sum_{s < t} (hi_s - lo_s), so every pair has its place before the scan starts and threads write to
// disjoint places. A block that meets no trait is not read.
//
// G is expanded to doubles kScanRowChunk rows at a time from the 2-bit codes (hard calls) or 0.01 dosages, so a
// thread holds b variants as n/4 (or 2n) bytes each plus kScanRowChunk x b doubles, and (k + active traits) x b
// products. A variant whose A lost more than 1 - 1/kKappaMax of g'g to the subtraction (strong collinearity with the
// covariates) is recomputed by explicit projection, and is in the span of the covariates (status 2, no statistics)
// when what the projection leaves of g is below kRankTol of its length, the tolerance that sets the rank of the
// covariates. The block size b follows from mem_gb (shared V, R and the output; then per thread). Products:
// la().dgemm_tn (SIMD dispatch of la_kernel.cpp).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "la_dispatch.h"
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

void build_pre(Pre& pre, const Rcpp::NumericMatrix& X, const Rcpp::NumericMatrix& Y, int T) {
  const int n = X.nrow(), kX = X.ncol(), q = Y.ncol();
  const std::size_t N = n;
  pre.n = n; pre.q = q;
  Eigen::MatrixXd A(n, kX);
  std::copy(X.begin(), X.end(), A.data());
  Eigen::BDCSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU);
  const Eigen::VectorXd& sig = svd.singularValues();
  int k = 0;
  while (k < kX && sig[k] > 0 && sig[k] >= kRankTol * sig[0]) ++k;
  if (k == 0) Rcpp::stop("The covariate matrix is zero.");
  pre.k = k;
  pre.V.assign(svd.matrixU().data(), svd.matrixU().data() + N * k);
  pre.R.assign(Y.begin(), Y.end());
  pre.rr.assign(q, 0.0);
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

// The variant ranges of the traits (ordered by lo) and the place of every pair in the output.
struct Windows {
  int q = 0;
  std::vector<int> lo, hi, top;    // top[t] = max(hi[0..t])
  std::vector<std::size_t> off;    // off[t] = pairs of the traits before t; off[q] = all pairs
  Windows(const Rcpp::IntegerVector& lo_, const Rcpp::IntegerVector& hi_, int p)
      : q(lo_.size()), lo(lo_.begin(), lo_.end()), hi(hi_.begin(), hi_.end()), top(lo_.size()), off(lo_.size() + 1, 0) {
    int running = 0;
    for (int t = 0; t < q; ++t) {
      if (lo[t] < 0 || hi[t] < lo[t] || hi[t] > p) Rcpp::stop("Range of trait %d is not within the variants.", t + 1);
      if (t > 0 && lo[t] < lo[t - 1]) Rcpp::stop("The traits are not ordered by the start of their ranges.");
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

}  // namespace

// path: .pgen or .bed; bed = true for a BED (A1 = .bim column 5 counted), else PGEN (ALT counted).
// var_idx, sample_idx: 0-based file rows; the scan variants are var_idx in the order given. X: n x kX covariates with
// the intercept. Y: n x q traits, no missing value, rows as sample_idx. lo, hi: per trait the 0-based range
// lo .. hi - 1 of scan variants it is tested with; lo not decreasing over the traits. block: variants per thread
// (0 = from mem_gb). Returns beta and se of every pair, trait by trait (trait t: hi_t - lo_t entries, variants
// lo_t .. hi_t - 1), NA for a variant with status 1 (monomorphic / all missing) or 2 (in the span of the covariates);
// per variant n_obs, a1_sum, status (NA when its block met no trait and was not read); per block of the scan the
// number of traits it met (active).
// [[Rcpp::export]]
Rcpp::List eqtl_stream_cpp(const std::string& path, bool bed, int n_raw, Rcpp::IntegerVector allele_ct,
                           Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx, Rcpp::NumericMatrix X,
                           Rcpp::NumericMatrix Y, Rcpp::IntegerVector lo, Rcpp::IntegerVector hi, int threads,
                           double mem_gb, int block) {
  const int T = plink4r::resolve_threads(threads);
  mem_gb = plink4r::resolve_mem_gb(mem_gb, T);
  const int p = var_idx.size(), ns = sample_idx.size(), q = Y.ncol();
  if (X.nrow() != ns || Y.nrow() != ns) Rcpp::stop("X / Y do not match the samples");
  if (lo.size() != q || hi.size() != q) Rcpp::stop("lo / hi do not match the traits");
  for (R_xlen_t i = 0; i < Y.size(); ++i)
    if (!std::isfinite(Y[i])) Rcpp::stop("Y has a missing or non-finite value.");
  const Windows win(lo, hi, p);
  const std::size_t pairs = win.off[q];
  plink4r::ScanFile F(path, bed, n_raw, allele_ct, var_idx, sample_idx, T);
  const bool dos = F.dosage;

  Pre pre;
  build_pre(pre, X, Y, T);
  const int k = pre.k;

  // block size from the memory budget
  const std::size_t N = ns;
  const int b_most = block > 0 ? block : plink4r::kScanBlockMax;
  const double shared = 8.0 * N * (k + q) + 16.0 * pairs + 8.0 * p * 3 + 8.0 * N * X.ncol() * 2;
  const double per_var = plink4r::ScanBlock::var_bytes(ns, dos) + 8.0 * std::min<int>(kScanRowChunk, ns) +
                         8.0 * (k + win.most_active(b_most, p) + 1);
  const double per_thread = 8.0 * N * 2 + F.thread_bytes();
  const int b = plink4r::scan_block_size(block, mem_gb, shared, per_thread, per_var, p, ns, T);
  const int most = win.most_active(b, p);

  // outputs
  Rcpp::NumericVector nobs(p, NA_REAL), asum(p, NA_REAL), beta(pairs, NA_REAL), se(pairs, NA_REAL);
  Rcpp::IntegerVector status(p, NA_INTEGER), active((p + b - 1) / b);
  const double df = static_cast<double>(ns) - k - 1;

  std::vector<Worker> wk(T);
  for (int t = 0; t < T; ++t) wk[t].init(pre, b, dos, most);
  std::vector<int> bad(T, -1);
  double* n_obs = nobs.begin();
  double* a1_sum = asum.begin();
  double* bt = beta.begin();
  double* st = se.begin();
  int nexp = 0;
  double products = 0.0;

  for (long long c0 = 0; c0 < p; c0 += static_cast<long long>(b) * T) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(T) reduction(+ : nexp, products)
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
      const long long failed = W.read(F, j0, nb, th, n_obs, a1_sum);
      if (failed >= 0) bad[th] = static_cast<int>(failed);
      block_products(pre, W, nb, first, last);
      const int qa = last - first;
      products += static_cast<double>(qa) * nb;
      for (int s = 0; s < nb; ++s) {
        const int j = static_cast<int>(j0 + s);
        if (W.mono[s]) { status[j] = 1; continue; }
        const double* m = &W.M[static_cast<std::size_t>(s) * k];
        const double* c = &W.C[static_cast<std::size_t>(s) * qa];
        double A = W.gg[s] - dot(m, m, k);
        const bool ill = !(A > W.gg[s] / kKappaMax);
        if (ill) { A = explicit_projection(pre, W, s); ++nexp; }
        if (!(A > kRankTol * kRankTol * W.gg[s])) { status[j] = 2; continue; }  // g is in the span of the covariates
        status[j] = 0;
        for (int u = first; u < last; ++u) {
          if (j < win.lo[u] || j >= win.hi[u]) continue;
          const double cu = ill ? dot(W.zt.data(), &pre.R[u * N], ns) : c[u - first];
          const std::size_t o = win.off[u] + static_cast<std::size_t>(j - win.lo[u]);
          bt[o] = cu / A;
          st[o] = std::sqrt(std::max(0.0, pre.rr[u] - cu * cu / A) / df / A);
        }
      }
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) Rcpp::stop("Genotype read failed at variant %d of %s", bad[t] + 1, path.c_str());
    Rcpp::checkUserInterrupt();
  }

  Rcpp::List out;
  out["beta"] = beta; out["se"] = se; out["n_obs"] = nobs; out["a1_sum"] = asum; out["status"] = status;
  out["df"] = df; out["rank"] = k; out["explicit"] = nexp; out["block"] = b; out["active"] = active;
  out["products"] = products; out["pairs"] = static_cast<double>(pairs); out["threads"] = T; out["mem_gb"] = mem_gb;
  out["dosage"] = dos; out["simd"] = std::string(plink4r::la().isa);
  return out;
}
