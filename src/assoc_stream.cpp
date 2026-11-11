// assoc_stream_cpp: single-trait association (GWAS, S = 1) and G x E (S > 1) of every selected variant of a
// PGEN or BED, streamed in blocks; the genotype matrix is never formed.
//
// Ported from the StreamBlock scheme of GEM 2 (GPL >= 3; Hong, Chen, Pham, Pan, Salehi Nasab). For variant i let
// Z = [g, g e_1, ..., g e_{S-1}] (n x S), X the covariates (n x k, columns 0..S-1 are the intercept and the e_j),
// W = diag(w) (w = mu (1 - mu) for a logistic null model, 1 for a linear one) and r the null-model residual.
// With a basis V of span(X), V'WV = I (thin SVD of W^{1/2} X, singular values below sqrt(eps) x largest dropped):
//   A = Z'WZ - M'M,  c = Z'r - M'V'r,  M = V'WZ                  (A = Zt'W Zt, c = Zt'r, Zt = Z - V M)
//   B = Z'R2 Z - M2'M - M'M2 + M'(V'R2 V)M,  M2 = V'R2 Z, R2 = diag(r^2)          (robust, linear only)
// and every n-long sum is an entry of L1'G or L2'(G o G), G the centred genotypes of a block (missing call = mean,
// i.e. 0 after centring), L1 = [wV e_j | r2 V e_j | r e_j], L2 = [w e_j e_l | r2 e_j e_l]_{j <= l}. G is expanded
// to doubles kRowChunk rows at a time from the 2-bit codes (hard calls) or 1/16384 dosages, so a thread holds
// b variants as n/4 (or 2n) bytes each plus 2 x kRowChunk x b doubles. A variant whose A_jj lost more than
// 1 - 1/kKappaMax of Z_j'W Z_j to the subtraction (strong collinearity with the covariates) is recomputed by
// explicit projection. The block size b follows from mem_gb (shared L1, L2, V, outputs; then per thread).
// Products: la().dgemm_tn (SIMD dispatch of la_kernel.cpp); small S x S algebra: Eigen.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "la_dispatch.h"
#include "pgen_reader.h"
#include "resources.h"

namespace {

constexpr int kRowChunk = 2048;          // rows of G expanded at a time (multiple of 4)
constexpr int kBlockMax = 4096;          // largest automatic block per thread
constexpr double kKappaMax = 1e4;        // explicit projection when Z_j'WZ_j / A_jj exceeds this
constexpr double kRankTol = 1.4901161193847656e-08;  // sqrt(DBL_EPSILON), as mgcv
constexpr std::uint16_t kDosMiss = 65535;

inline int tid() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline int pair_index(int S, int j, int l) {
  if (j > l) std::swap(j, l);
  return j * S - j * (j - 1) / 2 + (l - j);
}

inline double dot(const double* a, const double* b, int len) {
  double s = 0.0;
  for (int i = 0; i < len; ++i) s += a[i] * b[i];
  return s;
}

// BED byte (4 samples, codes 00 = two A1, 01 missing, 10 one A1, 11 none) -> PGEN-style byte (A1 count, 3 missing)
struct BedTable {
  unsigned char t[256];
  BedTable() {
    static const unsigned char m[4] = {2, 3, 1, 0};
    for (int b = 0; b < 256; ++b) {
      unsigned v = 0;
      for (int k = 0; k < 4; ++k) v |= static_cast<unsigned>(m[(b >> (2 * k)) & 3]) << (2 * k);
      t[b] = static_cast<unsigned char>(v);
    }
  }
};

// Null-model quantities shared by all threads.
struct Pre {
  int n = 0, kX = 0, k = 0, S = 0, P = 0, m1 = 0, m2 = 0;
  bool logistic = false, robust = false;
  std::vector<double> L1, L2, V, E, w, r, Vr, HV;
};

void build_pre(Pre& pre, const Rcpp::NumericMatrix& X, int S, const Rcpp::NumericVector& resid,
               const Rcpp::NumericVector& wt, bool robust) {
  const int n = X.nrow(), kX = X.ncol();
  const std::size_t N = n;
  pre.n = n; pre.kX = kX; pre.S = S; pre.P = S * (S + 1) / 2;
  pre.logistic = wt.size() > 0; pre.robust = robust;
  pre.r.assign(resid.begin(), resid.end());
  if (pre.logistic) pre.w.assign(wt.begin(), wt.end());
  pre.E.resize(N * S);
  for (int j = 0; j < S; ++j)
    for (std::size_t i = 0; i < N; ++i) pre.E[j * N + i] = X(i, j);

  Eigen::MatrixXd A(n, kX);
  for (int c = 0; c < kX; ++c)
    for (int i = 0; i < n; ++i) A(i, c) = X(i, c) * (pre.logistic ? std::sqrt(pre.w[i]) : 1.0);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU);
  const Eigen::VectorXd& sig = svd.singularValues();
  int k = 0;
  while (k < kX && sig[k] > 0 && sig[k] >= kRankTol * sig[0]) ++k;
  if (k == 0) Rcpp::stop("The covariate matrix is zero.");
  pre.k = k;
  const Eigen::MatrixXd& U = svd.matrixU();
  A.resize(0, 0);
  pre.V.resize(N * k);
  std::vector<double> wV(N * k);
  for (int c = 0; c < k; ++c)
    for (std::size_t i = 0; i < N; ++i) {
      const double u = U(i, c), s = pre.logistic ? std::sqrt(pre.w[i]) : 1.0;
      pre.V[c * N + i] = u / s;
      wV[c * N + i] = u * s;
    }

  const int rb = robust ? 2 : 1;
  pre.m1 = k * S * rb + S;
  pre.m2 = pre.P * rb;
  pre.L1.assign(N * pre.m1, 0.0);
  for (int j = 0; j < S; ++j) {
    const double* e = &pre.E[j * N];
    for (int c = 0; c < k; ++c) {
      double* d = &pre.L1[(std::size_t)(j * k + c) * N];
      const double* v = &wV[c * N];
      for (std::size_t i = 0; i < N; ++i) d[i] = v[i] * e[i];
      if (robust) {
        double* d2 = &pre.L1[(std::size_t)(k * S + j * k + c) * N];
        const double* v2 = &pre.V[c * N];
        for (std::size_t i = 0; i < N; ++i) d2[i] = pre.r[i] * pre.r[i] * v2[i] * e[i];
      }
    }
    double* d3 = &pre.L1[(std::size_t)(k * S * rb + j) * N];
    for (std::size_t i = 0; i < N; ++i) d3[i] = pre.r[i] * e[i];
  }
  pre.L2.assign(N * pre.m2, 0.0);
  for (int j = 0; j < S; ++j)
    for (int l = j; l < S; ++l) {
      const int p = pair_index(S, j, l);
      const double* ej = &pre.E[j * N];
      const double* el = &pre.E[l * N];
      double* d = &pre.L2[(std::size_t)p * N];
      for (std::size_t i = 0; i < N; ++i) d[i] = ej[i] * el[i] * (pre.logistic ? pre.w[i] : 1.0);
      if (robust) {
        double* d2 = &pre.L2[(std::size_t)(pre.P + p) * N];
        for (std::size_t i = 0; i < N; ++i) d2[i] = pre.r[i] * pre.r[i] * ej[i] * el[i];
      }
    }
  pre.Vr.assign(k, 0.0);
  for (int c = 0; c < k; ++c) pre.Vr[c] = dot(&pre.V[c * N], pre.r.data(), n);
  if (robust) {
    pre.HV.assign((std::size_t)k * k, 0.0);
    for (int a = 0; a < k; ++a)
      for (int b = a; b < k; ++b) {
        double s = 0.0;
        const double* va = &pre.V[a * N];
        const double* vb = &pre.V[b * N];
        for (std::size_t i = 0; i < N; ++i) s += pre.r[i] * pre.r[i] * va[i] * vb[i];
        pre.HV[a + b * k] = pre.HV[b + a * k] = s;
      }
  }
}

// Variants of one thread: compact genotypes over the analysis samples and the block products.
struct Worker {
  int cap = 0, n = 0;
  bool dosage = false;
  std::size_t sb = 0;
  std::vector<unsigned char> g8;   // cap x sb: 2-bit codes (0/1/2 = count, 3 = missing)
  std::vector<std::uint16_t> g16;  // cap x n: dosage x 16384, kDosMiss = missing
  std::vector<double> mean, G, G2, C1, C2, z, Z, Zt;
  std::vector<int> mono;
  std::vector<unsigned char> raw;  // BED row
  void init(const Pre& pre, int capacity, bool dos) {
    cap = capacity; n = pre.n; dosage = dos;
    sb = (static_cast<std::size_t>(n) + 3) / 4;
    if (dosage) g16.assign(static_cast<std::size_t>(cap) * n, 0); else g8.assign(cap * sb, 0);
    mean.assign(cap, 0.0); mono.assign(cap, 0);
    const std::size_t R = std::min(kRowChunk, n);
    G.assign(R * cap, 0.0); G2.assign(R * cap, 0.0);
    C1.assign(static_cast<std::size_t>(pre.m1) * cap, 0.0);
    C2.assign(static_cast<std::size_t>(pre.m2) * cap, 0.0);
  }
  void centred(int s, int row0, int rows, double* out) const {
    const double mu = mean[s];
    if (dosage) {
      const std::uint16_t* d = &g16[static_cast<std::size_t>(s) * n + row0];
      for (int i = 0; i < rows; ++i) out[i] = d[i] == kDosMiss ? 0.0 : d[i] * (1.0 / 16384.0) - mu;
      return;
    }
    const double lut[4] = {-mu, 1.0 - mu, 2.0 - mu, 0.0};
    const unsigned char* p = &g8[s * sb + (row0 >> 2)];
    int i = 0;
    for (; i + 4 <= rows; i += 4) {
      const unsigned b = p[i >> 2];
      out[i] = lut[b & 3]; out[i + 1] = lut[(b >> 2) & 3]; out[i + 2] = lut[(b >> 4) & 3]; out[i + 3] = lut[b >> 6];
    }
    for (; i < rows; ++i) out[i] = lut[(p[i >> 2] >> ((i & 3) << 1)) & 3];
  }
};

// c (S), A (S x S), B (S x S, robust) of slot s by explicit projection.
void explicit_projection(const Pre& pre, Worker& wk, int s, double* c, double* A, double* B) {
  const int n = pre.n, k = pre.k, S = pre.S;
  const std::size_t N = n;
  if (wk.z.empty()) { wk.z.resize(N); wk.Z.resize(N * S); wk.Zt.resize(N * S); }
  wk.centred(s, 0, n, wk.z.data());
  for (int j = 0; j < S; ++j)
    for (std::size_t i = 0; i < N; ++i) wk.Z[j * N + i] = wk.z[i] * pre.E[j * N + i];
  std::vector<double> M(static_cast<std::size_t>(k) * S);
  for (int j = 0; j < S; ++j)
    for (int c2 = 0; c2 < k; ++c2) M[j * k + c2] = dot(&pre.L1[(std::size_t)(j * k + c2) * N], wk.z.data(), n);
  wk.Zt = wk.Z;
  for (int j = 0; j < S; ++j)
    for (int c2 = 0; c2 < k; ++c2) {
      const double m = M[j * k + c2];
      const double* v = &pre.V[c2 * N];
      double* zt = &wk.Zt[j * N];
      for (std::size_t i = 0; i < N; ++i) zt[i] -= v[i] * m;
    }
  for (int j = 0; j < S; ++j) c[j] = dot(&wk.Zt[j * N], pre.r.data(), n);
  for (int j = 0; j < S; ++j)
    for (int l = j; l < S; ++l) {
      const double* a = &wk.Zt[j * N];
      const double* b = &wk.Zt[l * N];
      double sa = 0.0, sr = 0.0;
      for (std::size_t i = 0; i < N; ++i) {
        const double ab = a[i] * b[i];
        sa += ab * (pre.logistic ? pre.w[i] : 1.0);
        sr += ab * pre.r[i] * pre.r[i];
      }
      A[j + l * S] = A[l + j * S] = sa;
      if (B) B[j + l * S] = B[l + j * S] = sr;
    }
}

// c, A, B of the nb slots of the worker; returns the number of explicit projections.
int block_stats(const Pre& pre, Worker& wk, int nb, double* cc, double* AA, double* BB) {
  const int n = pre.n, k = pre.k, S = pre.S, m1 = pre.m1, m2 = pre.m2, kS = k * S;
  const plink4r::LaFns& LA = plink4r::la();
  for (int row0 = 0; row0 < n; row0 += kRowChunk) {
    const int rows = std::min(kRowChunk, n - row0);
    for (int s = 0; s < nb; ++s) {
      double* g = &wk.G[static_cast<std::size_t>(s) * rows];
      double* g2 = &wk.G2[static_cast<std::size_t>(s) * rows];
      wk.centred(s, row0, rows, g);
      for (int i = 0; i < rows; ++i) g2[i] = g[i] * g[i];
    }
    const int beta = row0 > 0;
    LA.dgemm_tn(&pre.L1[row0], rows, m1, n, wk.G.data(), nb, rows, wk.C1.data(), m1, beta);
    LA.dgemm_tn(&pre.L2[row0], rows, m2, n, wk.G2.data(), nb, rows, wk.C2.data(), m2, beta);
  }
  std::vector<double> HM(static_cast<std::size_t>(k) * S);
  int nexp = 0;
  for (int s = 0; s < nb; ++s) {
    double* c = cc + static_cast<std::size_t>(s) * S;
    double* A = AA + static_cast<std::size_t>(s) * S * S;
    double* B = pre.robust ? BB + static_cast<std::size_t>(s) * S * S : nullptr;
    if (wk.mono[s]) continue;
    const double* c1 = &wk.C1[static_cast<std::size_t>(s) * m1];
    const double* c2 = &wk.C2[static_cast<std::size_t>(s) * m2];
    const double* M = c1;
    const double* M2 = c1 + kS;
    const double* zr = c1 + kS * (pre.robust ? 2 : 1);
    bool ill = false;
    for (int j = 0; j < S; ++j) {
      for (int l = 0; l < S; ++l) A[j + l * S] = c2[pair_index(S, j, l)] - dot(M + j * k, M + l * k, k);
      if (!(A[j + j * S] > c2[pair_index(S, j, j)] / kKappaMax)) ill = true;
    }
    if (ill) { explicit_projection(pre, wk, s, c, A, B); ++nexp; continue; }
    for (int j = 0; j < S; ++j) c[j] = zr[j] - dot(M + j * k, pre.Vr.data(), k);
    if (pre.robust) {
      for (int l = 0; l < S; ++l)
        for (int a = 0; a < k; ++a) HM[l * k + a] = dot(&pre.HV[a * k], M + l * k, k);  // HV symmetric
      for (int j = 0; j < S; ++j)
        for (int l = 0; l < S; ++l)
          B[j + l * S] = c2[pre.P + pair_index(S, j, l)] - dot(M2 + j * k, M + l * k, k) -
                         dot(M + j * k, M2 + l * k, k) + dot(M + j * k, &HM[l * k], k);
    }
  }
  return nexp;
}

}  // namespace

// path: .pgen or .bed; bed = true for a BED (A1 = .bim column 5 counted), else PGEN (ALT counted).
// var_idx, sample_idx: 0-based file rows. X: n x kX covariates, columns 0..S-1 = intercept and the
// interaction variables. resid: null residual y - mu. w: mu (1 - mu) for a logistic null, empty for linear.
// block: variants per thread (0 = from mem_gb). Returns per variant (rows) the statistics below; beta, se in
// columns of Z = [g, g e_1, ...]; status 0 ok, 1 monomorphic / all missing.
// [[Rcpp::export]]
Rcpp::List assoc_stream_cpp(const std::string& path, bool bed, int n_raw, Rcpp::IntegerVector allele_ct,
                            Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx, Rcpp::NumericMatrix X,
                            int S, int n_exp, Rcpp::NumericVector resid, Rcpp::NumericVector w, bool robust,
                            int threads, double mem_gb, int block) {
  const int T = plink4r::resolve_threads(threads);
  mem_gb = plink4r::resolve_mem_gb(mem_gb, T);
  const int p = var_idx.size(), ns = sample_idx.size();
  if (X.nrow() != ns || resid.size() != ns || (w.size() && w.size() != ns)) Rcpp::stop("X / resid / w do not match");
  if (S < 1 || S > X.ncol()) Rcpp::stop("S out of range");
  if (n_exp < 0 || n_exp > S - 1) Rcpp::stop("n_exp out of range");
  if (robust && w.size()) Rcpp::stop("robust standard errors are for the linear model only");
  const std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  const bool ident = plink4r::identity_samples(si.data(), ns, n_raw);

  std::unique_ptr<plink4r::PgenFile> pf;
  if (!bed) pf.reset(new plink4r::PgenFile(path, n_raw, std::vector<int>(allele_ct.begin(), allele_ct.end()), T));
  const bool dos = pf && pf->dosage();
  const std::size_t rawb = (static_cast<std::size_t>(n_raw) + 3) / 4;

  Pre pre;
  build_pre(pre, X, S, resid, w, robust);
  const int k = pre.k, P = pre.P;
  double rr = 0.0;
  for (double v : pre.r) rr += v * v;

  // block size from the memory budget
  const double GiB = 1073741824.0;
  const std::size_t N = ns;
  const double shared = 8.0 * N * (pre.m1 + pre.m2 + k + S + 4) + 8.0 * p * (4 * S + 2 * S * S + 12) +
                        8.0 * N * X.ncol() * 2;
  const double per_var = (dos ? 2.0 * N : (N + 3) / 4.0) + 16.0 * std::min<int>(kRowChunk, ns) +
                         8.0 * (pre.m1 + pre.m2) + 8.0 * (S + 2 * S * S + 4);
  const double per_thread = 8.0 * N * (1 + 2 * S) + (bed ? rawb : (dos ? 3.0 * n_raw : rawb)) + 1048576.0;
  const double budget = mem_gb * GiB * 0.95 - 268435456.0 - shared;
  const int per_thread_max = (p + T - 1) / T;
  int b = block;
  if (b <= 0) {
    const double room = budget / T - per_thread;
    if (room < per_var)
      Rcpp::stop("mem_gb = %.2f is too small: about %.2f GB are needed for %d samples and %d threads.", mem_gb,
                 (shared + T * (per_thread + per_var) + 268435456.0) / 0.95 / GiB, ns, T);
    b = static_cast<int>(std::min<double>(kBlockMax, std::floor(room / per_var)));
  }
  b = std::max(1, std::min(b, std::max(1, per_thread_max)));

  // outputs
  Rcpp::NumericVector nobs(p), asum(p);
  Rcpp::IntegerVector status(p);
  Rcpp::NumericMatrix beta(p, S), se(p, S), se_rob(p, robust ? S : 0);
  Rcpp::NumericVector beta_m(p), se_m(p), se_m_rob(p), chi_int(p), chi_joint(p), chi_int_rob(p), chi_joint_rob(p);
  std::fill(beta.begin(), beta.end(), NA_REAL); std::fill(se.begin(), se.end(), NA_REAL);
  std::fill(se_rob.begin(), se_rob.end(), NA_REAL);
  for (Rcpp::NumericVector* v : {&beta_m, &se_m, &se_m_rob, &chi_int, &chi_joint, &chi_int_rob, &chi_joint_rob})
    std::fill(v->begin(), v->end(), NA_REAL);
  const double df = static_cast<double>(ns) - k - S, df_m = static_cast<double>(ns) - k - 1;

  std::vector<Worker> wk(T);
  std::vector<std::vector<double>> cbuf(T), abuf(T), bbuf(T);
  for (int t = 0; t < T; ++t) {
    wk[t].init(pre, b, dos);
    cbuf[t].resize(static_cast<std::size_t>(b) * S);
    abuf[t].resize(static_cast<std::size_t>(b) * S * S);
    if (robust) bbuf[t].resize(static_cast<std::size_t>(b) * S * S);
  }
  std::vector<std::ifstream> bedf;
  if (bed) {
    bedf.resize(T);
    for (int t = 0; t < T; ++t) {
      bedf[t].open(path, std::ios::binary);
      if (!bedf[t]) Rcpp::stop("Cannot open " + path);
      wk[t].raw.resize(rawb);
    }
  }
  static const BedTable BT;
  std::vector<int> bad(T, -1);
  int nexp = 0;

  for (long long c0 = 0; c0 < p; c0 += static_cast<long long>(b) * T) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(T) reduction(+ : nexp)
#endif
    for (int t = 0; t < T; ++t) {
      const int th = tid();
      Worker& W = wk[th];
      const long long j0 = c0 + static_cast<long long>(t) * b;
      if (j0 >= p) continue;
      const int nb = static_cast<int>(std::min<long long>(b, p - j0));
      // read and compact
      for (int s = 0; s < nb; ++s) {
        const long long j = j0 + s;
        double no = 0, sum = 0, sq = 0;
        if (dos) {
          const std::uint16_t* d = pf->read_dosage(static_cast<std::size_t>(vi[j]), th);
          if (!d) { bad[th] = static_cast<int>(j); W.mono[s] = 1; continue; }
          std::uint16_t* o = &W.g16[static_cast<std::size_t>(s) * ns];
          for (int i = 0; i < ns; ++i) {
            const std::uint16_t v = d[si[i]];
            o[i] = v;
            if (v == kDosMiss) continue;
            const double x = v * (1.0 / 16384.0);
            no += 1; sum += x; sq += x * x;
          }
        } else {
          const unsigned char* gv;
          if (bed) {
            std::ifstream& f = bedf[th];
            f.seekg(3 + static_cast<std::streamoff>(vi[j]) * static_cast<std::streamoff>(rawb));
            f.read(reinterpret_cast<char*>(W.raw.data()), static_cast<std::streamsize>(rawb));
            if (!f) { bad[th] = static_cast<int>(j); W.mono[s] = 1; f.clear(); continue; }
            for (std::size_t y = 0; y < rawb; ++y) W.raw[y] = BT.t[W.raw[y]];
            gv = W.raw.data();
          } else {
            gv = pf->read(static_cast<std::size_t>(vi[j]), th);
            if (!gv) { bad[th] = static_cast<int>(j); W.mono[s] = 1; continue; }
          }
          unsigned char* dst = &W.g8[s * W.sb];
          plink4r::compact_genovec(gv, si.data(), ns, ident, dst);
          long long cnt[4] = {0, 0, 0, 0};
          for (std::size_t y = 0; y < W.sb; ++y) {
            const unsigned v = dst[y];
            ++cnt[v & 3]; ++cnt[(v >> 2) & 3]; ++cnt[(v >> 4) & 3]; ++cnt[v >> 6];
          }
          cnt[0] -= static_cast<long long>(W.sb * 4 - ns);  // padding is code 0
          no = static_cast<double>(cnt[0] + cnt[1] + cnt[2]);
          sum = static_cast<double>(cnt[1] + 2 * cnt[2]);
          sq = static_cast<double>(cnt[1] + 4 * cnt[2]);
        }
        nobs[j] = no; asum[j] = sum;
        const double mu = no > 0 ? sum / no : 0.0, v = no > 0 ? sq / no - mu * mu : 0.0;
        W.mean[s] = mu;
        W.mono[s] = (no < 2 || !(v > 1e-12 * std::max(1.0, sq / std::max(no, 1.0)))) ? 1 : 0;
      }
      double* cc = cbuf[th].data();
      double* AA = abuf[th].data();
      double* BB = robust ? bbuf[th].data() : nullptr;
      nexp += block_stats(pre, W, nb, cc, AA, BB);
      // statistics
      for (int s = 0; s < nb; ++s) {
        const long long j = j0 + s;
        if (W.mono[s]) { status[j] = 1; continue; }
        Eigen::Map<const Eigen::MatrixXd> A(AA + static_cast<std::size_t>(s) * S * S, S, S);
        Eigen::Map<const Eigen::VectorXd> c(cc + static_cast<std::size_t>(s) * S, S);
        const double a00 = A(0, 0), c0v = c[0];
        // marginal G (the model without the interaction terms)
        beta_m[j] = c0v / a00;
        if (pre.logistic) se_m[j] = std::sqrt(1.0 / a00);
        else se_m[j] = std::sqrt(std::max(0.0, rr - c0v * c0v / a00) / df_m / a00);
        if (robust) se_m_rob[j] = std::sqrt(BB[static_cast<std::size_t>(s) * S * S] / (a00 * a00));
        Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
        if (ldlt.info() != Eigen::Success) { status[j] = 2; continue; }
        const Eigen::MatrixXd Ainv = ldlt.solve(Eigen::MatrixXd::Identity(S, S));
        const Eigen::VectorXd bt = Ainv * c;
        const double s2 = pre.logistic ? 1.0 : std::max(0.0, rr - c.dot(bt)) / df;
        const Eigen::MatrixXd Vm = s2 * Ainv;
        for (int q = 0; q < S; ++q) { beta(j, q) = bt[q]; se(j, q) = std::sqrt(Vm(q, q)); }
        Eigen::MatrixXd Vr;
        if (robust) {
          Eigen::Map<const Eigen::MatrixXd> Bm(BB + static_cast<std::size_t>(s) * S * S, S, S);
          Vr = Ainv * Bm * Ainv;
          for (int q = 0; q < S; ++q) se_rob(j, q) = std::sqrt(Vr(q, q));
        }
        if (n_exp > 0) {  // Wald tests of the exposure interactions (terms 1..n_exp) and of G with them (0..n_exp)
          const int ne = n_exp;
          const Eigen::VectorXd bi = bt.segment(1, ne), bj = bt.head(ne + 1);
          chi_int[j] = bi.dot(Vm.block(1, 1, ne, ne).ldlt().solve(bi));
          chi_joint[j] = bj.dot(Vm.topLeftCorner(ne + 1, ne + 1).ldlt().solve(bj));
          if (robust) {
            chi_int_rob[j] = bi.dot(Vr.block(1, 1, ne, ne).ldlt().solve(bi));
            chi_joint_rob[j] = bj.dot(Vr.topLeftCorner(ne + 1, ne + 1).ldlt().solve(bj));
          }
        }
      }
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) Rcpp::stop("Genotype read failed at variant %d of %s", bad[t] + 1, path.c_str());
    Rcpp::checkUserInterrupt();
  }

  Rcpp::List out;
  out["n_obs"] = nobs; out["a1_sum"] = asum; out["status"] = status; out["beta"] = beta; out["se"] = se;
  out["se_robust"] = se_rob; out["beta_marginal"] = beta_m; out["se_marginal"] = se_m;
  out["se_marginal_robust"] = se_m_rob; out["chisq_int"] = chi_int; out["chisq_joint"] = chi_joint;
  out["chisq_int_robust"] = chi_int_rob; out["chisq_joint_robust"] = chi_joint_rob; out["df"] = df;
  out["df_marginal"] = df_m; out["rank"] = k; out["explicit"] = nexp; out["block"] = b; out["threads"] = T;
  out["mem_gb"] = mem_gb; out["dosage"] = dos; out["simd"] = std::string(plink4r::la().isa);
  return out;
}
