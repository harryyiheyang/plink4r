// gwas_stream_cpp, gxe_stream_cpp: single-trait association (GWAS, S = 1) and G x E (S > 1) of every variant
// of a PGEN or BED, streamed in blocks; the genotype matrix is never formed and no table of the results either:
// the rows of a group of blocks are written to the result file as soon as they are computed, each one CHR, BP, A1,
// A2 of the variant as its .bim / .pvar line has them, then its statistics (ResultText, variant_text.h).
//
// Ported from the StreamBlock scheme of GEM 2 (GPL >= 3; Hong, Chen, Pham, Pan, Salehi Nasab). For variant i let
// Z = [g, g e_1, ..., g e_{S-1}] (n x S), X the covariates (n x k, columns 0..S-1 are the intercept and the e_j),
// W = diag(w) (w = mu (1 - mu) for a logistic null model, 1 for a linear one) and r the null-model residual.
// With a basis V of span(X), V'WV = I (thin SVD of W^{1/2} X, singular values below sqrt(eps) x largest dropped):
//   A = Z'WZ - M'M,  c = Z'r - M'V'r,  M = V'WZ                  (A = Zt'W Zt, c = Zt'r, Zt = Z - V M)
//   B = Z'R2 Z - M2'M - M'M2 + M'(V'R2 V)M,  M2 = V'R2 Z, R2 = diag(r^2)          (robust, linear only)
// and every n-long sum is an entry of L1'G or L2'(G o G), G the centred genotypes of a block (missing call = mean,
// i.e. 0 after centring), L1 = [wV e_j | r2 V e_j | r e_j], L2 = [w e_j e_l | r2 e_j e_l]_{j <= l}. G is expanded
// to doubles kScanRowChunk rows at a time from the 2-bit codes (hard calls) or 0.01 dosages, so a thread holds
// b variants as n/4 (or 2n) bytes each plus 2 x kScanRowChunk x b doubles. A variant whose A_jj lost more than
// 1 - 1/kKappaMax of Z_j'W Z_j to the subtraction (strong collinearity with the covariates) is recomputed by
// explicit projection. The block size b follows from mem_gb (shared L1, L2, V; then per thread).
// Products: la().dgemm_tn (SIMD dispatch of la_kernel.cpp); small S x S algebra: Eigen.
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
#include "p_value.h"
#include "scan_block.h"
#include "resources.h"

namespace {

using plink4r::kScanRowChunk;
constexpr double kKappaMax = 1e4;        // explicit projection when Z_j'WZ_j / A_jj exceeds this
constexpr double kRankTol = 1.4901161193847656e-08;  // sqrt(DBL_EPSILON), as mgcv
constexpr double kRowBytes = 512.0;      // text of a variant while its block is written, for the memory budget


inline int pair_index(int S, int j, int l) {
  if (j > l) std::swap(j, l);
  return j * S - j * (j - 1) / 2 + (l - j);
}

inline double dot(const double* a, const double* b, int len) {
  double s = 0.0;
  for (int i = 0; i < len; ++i) s += a[i] * b[i];
  return s;
}

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

// Variants of one thread: their compact genotypes (scan_block.h) and the block products.
struct Worker : plink4r::ScanBlock {
  std::vector<double> G, G2, C1, C2, z, Z, Zt;
  void init(const Pre& pre, int capacity, bool dos) {
    ScanBlock::init(capacity, pre.n, dos);
    const std::size_t R = std::min(kScanRowChunk, n);
    G.assign(R * cap, 0.0); G2.assign(R * cap, 0.0);
    C1.assign(static_cast<std::size_t>(pre.m1) * cap, 0.0);
    C2.assign(static_cast<std::size_t>(pre.m2) * cap, 0.0);
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
  for (int row0 = 0; row0 < n; row0 += kScanRowChunk) {
    const int rows = std::min(kScanRowChunk, n - row0);
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

// What the scan computes of one variant; NaN where a statistic was not computed (no variation among the samples, or
// a singular Z'WZ, in which case the marginal statistics are still there).
struct Stat {
  double n_obs = 0, a1_sum = 0;
  const double* beta = nullptr;     // S: the columns of Z = [g, g e_1, ...]
  const double* se = nullptr;
  const double* se_rob = nullptr;
  double beta_m = NA_REAL, se_m = NA_REAL, se_m_rob = NA_REAL;  // G alone (the model without the interaction terms)
  double chi_int = NA_REAL, chi_joint = NA_REAL, chi_int_rob = NA_REAL, chi_joint_rob = NA_REAL;
};
// What a row needs of the run.
struct Run {
  int S = 1, n_exp = 0;
  bool logistic = false, robust = false;
  double df = 0, df_m = 0;
};

// The scan of every variant of the file, written as it goes (ResultText of variant_text.h): CHR, BP, A1, A2 of each
// variant, then the columns `statistics`: A1Freq, N and what else `row` puts into v for the variant (v:
// statistics.size() - 1 numbers, A1Freq first; N is not one of them).
// file: the genotype file as opened (geno_table.h; the external pointer of geno_file()): a BED counts A1 = .bim
// column 5, a PGEN the ALT allele of each variant. out: the text file to write. sample_idx: 0-based file rows of the analysed samples. X: n x kX
// covariates, columns 0..S-1 = intercept and the interaction variables (the n_exp exposures first, then the
// interaction covariates). resid: null residual y - mu. w: mu (1 - mu) for a logistic null, empty for linear.
// block: variants per thread (0 = from mem_gb).
template <class Row>
Rcpp::List scan(SEXP file, const std::string& out,
                const std::vector<std::string>& statistics, const Rcpp::IntegerVector& sample_idx,
                const Rcpp::NumericMatrix& X, int S, int n_exp, const Rcpp::NumericVector& resid,
                const Rcpp::NumericVector& w, bool robust, int threads, double mem_gb, int block, Row row) {
  const int T = plink4r::resolve_threads(threads);
  mem_gb = plink4r::resolve_mem_gb(mem_gb, T);
  const int ns = sample_idx.size();
  if (X.nrow() != ns || resid.size() != ns || (w.size() && w.size() != ns)) Rcpp::stop("X / resid / w do not match");
  if (S < 1 || S > X.ncol()) Rcpp::stop("S out of range");
  if (n_exp < 0 || n_exp > S - 1) Rcpp::stop("n_exp out of range");
  if (robust && w.size()) Rcpp::stop("robust standard errors are for the linear model only");
  const plink4r::GenoTable& gt = plink4r::geno_table(file);
  plink4r::ScanFile F(gt, sample_idx, T);
  const bool dos = F.dosage;
  const int p = static_cast<int>(F.vi.size());

  Pre pre;
  build_pre(pre, X, S, resid, w, robust);
  const int k = pre.k;
  double rr = 0.0;
  for (double v : pre.r) rr += v * v;
  Run run;
  run.S = S; run.n_exp = n_exp; run.logistic = pre.logistic; run.robust = robust;
  run.df = static_cast<double>(ns) - k - S; run.df_m = static_cast<double>(ns) - k - 1;
  const double df = run.df, df_m = run.df_m;

  // block size from the memory budget
  const std::size_t N = ns;
  const double shared = 8.0 * N * (pre.m1 + pre.m2 + k + S + 4) + 8.0 * N * X.ncol() * 2;
  const double per_var = plink4r::ScanBlock::var_bytes(ns, dos) + 16.0 * std::min<int>(kScanRowChunk, ns) +
                         8.0 * (pre.m1 + pre.m2) + 8.0 * (S + 2 * S * S + 4) + kRowBytes;
  const double per_thread = 8.0 * N * (1 + 2 * S) + F.thread_bytes();
  const int b = plink4r::scan_block_size(block, mem_gb, shared, per_thread, per_var, p, ns, T);

  plink4r::ResultText text(gt.variants, out, statistics);
  const std::size_t nv = statistics.size() - 1;
  std::vector<Worker> wk(T);
  std::vector<std::vector<double>> cbuf(T), abuf(T), bbuf(T), nbuf(T), sbuf(T), qbuf(T), vbuf(T);
  std::vector<std::string> rows(T);
  for (int t = 0; t < T; ++t) {
    wk[t].init(pre, b, dos);
    cbuf[t].resize(static_cast<std::size_t>(b) * S);
    abuf[t].resize(static_cast<std::size_t>(b) * S * S);
    if (robust) bbuf[t].resize(static_cast<std::size_t>(b) * S * S);
    nbuf[t].resize(b); sbuf[t].resize(b);
    qbuf[t].resize(static_cast<std::size_t>(3) * S);
    vbuf[t].resize(nv);
  }
  std::vector<int> bad(T, -1);
  int nexp = 0;
  long long n_mono = 0, n_sing = 0;

  for (long long c0 = 0; c0 < p; c0 += static_cast<long long>(b) * T) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(T) reduction(+ : nexp, n_mono, n_sing)
#endif
    for (int t = 0; t < T; ++t) {
      std::string& o = rows[t];
      o.clear();
      const int th = plink4r::thread_id();
      Worker& W = wk[th];
      const long long j0 = c0 + static_cast<long long>(t) * b;
      if (j0 >= p) continue;
      const int nb = static_cast<int>(std::min<long long>(b, p - j0));
      double* n_obs = nbuf[th].data();
      double* a1_sum = sbuf[th].data();
      const long long failed = W.read(F, j0, nb, th, n_obs, a1_sum);
      if (failed >= 0) bad[th] = static_cast<int>(failed);
      double* cc = cbuf[th].data();
      double* AA = abuf[th].data();
      double* BB = robust ? bbuf[th].data() : nullptr;
      nexp += block_stats(pre, W, nb, cc, AA, BB);
      double* beta = qbuf[th].data();
      double* se = beta + S;
      double* se_rob = se + S;
      double* v = vbuf[th].data();
      for (int s = 0; s < nb; ++s) {
        const long long j = j0 + s;
        Stat st;
        st.n_obs = n_obs[s]; st.a1_sum = a1_sum[s];
        st.beta = beta; st.se = se; st.se_rob = se_rob;
        std::fill(beta, beta + 3 * S, NA_REAL);
        for (;;) {  // the statistics of the variant, as far as they can be computed
          if (W.mono[s]) { ++n_mono; break; }
          Eigen::Map<const Eigen::MatrixXd> A(AA + static_cast<std::size_t>(s) * S * S, S, S);
          Eigen::Map<const Eigen::VectorXd> c(cc + static_cast<std::size_t>(s) * S, S);
          const double a00 = A(0, 0), c0v = c[0];
          // marginal G (the model without the interaction terms)
          st.beta_m = c0v / a00;
          if (pre.logistic) st.se_m = std::sqrt(1.0 / a00);
          else st.se_m = std::sqrt(std::max(0.0, rr - c0v * c0v / a00) / df_m / a00);
          if (robust) st.se_m_rob = std::sqrt(BB[static_cast<std::size_t>(s) * S * S] / (a00 * a00));
          Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
          if (ldlt.info() != Eigen::Success) { ++n_sing; break; }
          const Eigen::MatrixXd Ainv = ldlt.solve(Eigen::MatrixXd::Identity(S, S));
          const Eigen::VectorXd bt = Ainv * c;
          const double s2 = pre.logistic ? 1.0 : std::max(0.0, rr - c.dot(bt)) / df;
          const Eigen::MatrixXd Vm = s2 * Ainv;
          for (int q = 0; q < S; ++q) { beta[q] = bt[q]; se[q] = std::sqrt(Vm(q, q)); }
          Eigen::MatrixXd Vr;
          if (robust) {
            Eigen::Map<const Eigen::MatrixXd> Bm(BB + static_cast<std::size_t>(s) * S * S, S, S);
            Vr = Ainv * Bm * Ainv;
            for (int q = 0; q < S; ++q) se_rob[q] = std::sqrt(Vr(q, q));
          }
          if (n_exp > 0) {  // Wald tests of the exposure interactions (terms 1..n_exp) and of G with them (0..n_exp)
            const int ne = n_exp;
            const Eigen::VectorXd bi = bt.segment(1, ne), bj = bt.head(ne + 1);
            st.chi_int = bi.dot(Vm.block(1, 1, ne, ne).ldlt().solve(bi));
            st.chi_joint = bj.dot(Vm.topLeftCorner(ne + 1, ne + 1).ldlt().solve(bj));
            if (robust) {
              st.chi_int_rob = bi.dot(Vr.block(1, 1, ne, ne).ldlt().solve(bi));
              st.chi_joint_rob = bj.dot(Vr.topLeftCorner(ne + 1, ne + 1).ldlt().solve(bj));
            }
          }
          break;
        }
        row(run, st, v);
        text.row(o, static_cast<std::size_t>(F.vi[j]), static_cast<long long>(st.n_obs), v);
      }
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) Rcpp::stop(F.failure(bad[t]));
    for (int t = 0; t < T; ++t) text.write(rows[t]);
    Rcpp::checkUserInterrupt();
  }
  text.finish();

  Rcpp::List res;
  res["variants"] = p; res["bytes"] = text.bytes();
  res["monomorphic"] = static_cast<double>(n_mono); res["singular"] = static_cast<double>(n_sing);
  res["df"] = df; res["df_marginal"] = df_m; res["rank"] = k; res["explicit"] = nexp; res["block"] = b;
  res["threads"] = T; res["mem_gb"] = mem_gb; res["dosage"] = dos; res["simd"] = std::string(plink4r::la().isa);
  return res;
}

}  // namespace

// The association of one trait with every variant of the file (S = 1): after CHR, BP, A1, A2 of each variant
// A1Freq, N, BETA, SE, T (linear) or Z (logistic), P, and SE_ROBUST, P_ROBUST with robust.
// Arguments as scan().
// [[Rcpp::export]]
Rcpp::List gwas_stream_cpp(SEXP file,
                           const std::string& out, Rcpp::IntegerVector sample_idx, Rcpp::NumericMatrix X,
                           Rcpp::NumericVector resid, Rcpp::NumericVector w, bool robust, int threads, double mem_gb,
                           int block) {
  std::vector<std::string> statistics = {"A1Freq", "N", "BETA", "SE", w.size() ? "Z" : "T", "P"};
  if (robust) { statistics.push_back("SE_ROBUST"); statistics.push_back("P_ROBUST"); }
  return scan(file, out, statistics, sample_idx, X, 1, 0, resid, w, robust, threads, mem_gb,
              block, [](const Run& g, const Stat& s, double* v) {
                const double z = s.beta[0] / s.se[0];
                v[0] = s.a1_sum / (2 * s.n_obs); v[1] = s.beta[0]; v[2] = s.se[0]; v[3] = z;
                v[4] = g.logistic ? plink4r::p_z(z) : plink4r::p_t(z, g.df);
                if (g.robust) { v[5] = s.se_rob[0]; v[6] = plink4r::p_z(s.beta[0] / s.se_rob[0]); }
              });
}

// Gene-environment interaction of one trait with every variant of the file: after CHR, BP, A1, A2 of each variant
// A1Freq, N, the marginal test of G (BETA_MARGINAL, SE_MARGINAL, P_MARGINAL), BETA
// and SE of G and of each interaction term (terms: the names of columns 1..S-1 of X), and the Wald tests of the
// exposure interactions (P_INT) and of G with them (P_JOINT); with robust the robust standard errors and tests
// beside them. Arguments as scan().
// [[Rcpp::export]]
Rcpp::List gxe_stream_cpp(SEXP file,
                          const std::string& out, Rcpp::IntegerVector sample_idx, Rcpp::NumericMatrix X, int S,
                          int n_exp, Rcpp::CharacterVector terms, Rcpp::NumericVector resid, Rcpp::NumericVector w,
                          bool robust, int threads, double mem_gb, int block) {
  if (terms.size() != S - 1) Rcpp::stop("terms do not match S");
  std::vector<std::string> statistics = {"A1Freq", "N", "BETA_MARGINAL", "SE_MARGINAL", "P_MARGINAL"};
  if (robust) { statistics.push_back("SE_MARGINAL_ROBUST"); statistics.push_back("P_MARGINAL_ROBUST"); }
  for (int q = 0; q < S; ++q) {
    const std::string nm = q == 0 ? std::string("G") : "GxE_" + Rcpp::as<std::string>(terms[q - 1]);
    statistics.push_back("BETA_" + nm); statistics.push_back("SE_" + nm);
    if (robust) statistics.push_back("SE_ROBUST_" + nm);
  }
  statistics.push_back("P_INT"); statistics.push_back("P_JOINT");
  if (robust) { statistics.push_back("P_INT_ROBUST"); statistics.push_back("P_JOINT_ROBUST"); }
  return scan(file, out, statistics, sample_idx, X, S, n_exp, resid, w, robust, threads,
              mem_gb, block, [](const Run& g, const Stat& s, double* v) {
                const double zm = s.beta_m / s.se_m;
                *v++ = s.a1_sum / (2 * s.n_obs); *v++ = s.beta_m; *v++ = s.se_m;
                *v++ = g.logistic ? plink4r::p_z(zm) : plink4r::p_t(zm, g.df_m);
                if (g.robust) { *v++ = s.se_m_rob; *v++ = plink4r::p_z(s.beta_m / s.se_m_rob); }
                for (int q = 0; q < g.S; ++q) {
                  *v++ = s.beta[q]; *v++ = s.se[q];
                  if (g.robust) *v++ = s.se_rob[q];
                }
                *v++ = plink4r::p_chisq(s.chi_int, g.n_exp); *v++ = plink4r::p_chisq(s.chi_joint, g.n_exp + 1);
                if (g.robust) {
                  *v++ = plink4r::p_chisq(s.chi_int_rob, g.n_exp); *v++ = plink4r::p_chisq(s.chi_joint_rob, g.n_exp + 1);
                }
              });
}
