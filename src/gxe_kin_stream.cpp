// gxe_kin_stream_cpp: gene-environment interaction of every selected variant of a PGEN or BED against the null
// model of gwas_kin_stream.cpp (the plink4r::KinNull of kin_sparse.h: linear, logistic or Cox, with a sparse
// kinship or with no related sample), streamed in blocks; the genotype matrix is never formed.
//
// The null model is the one of the association test, with the same kinship parameter; the exposures are columns
// of its design. For a variant g (centred, missing call = mean) and the S columns e_0 = 1, e_1, ..., e_{S-1}
// (exposures, then interaction covariates), Z = [g e_0, ..., g e_{S-1}] and, with Omega, P, r, Q of the null
// model (gwas_kin_stream.cpp),
//   c = Z' r (S),   A = Z' P Z = Z' Omega Z - (Q' Z)' (Q' Z) (S x S),   beta = A^-1 c,   Var(beta) = A^-1:
// one step from the null model for all S coefficients (generalised least squares for the linear model). The
// first entries, c_0 and A_00, are the association test of gwas_kin_stream.cpp (the marginal effect of g); the
// interaction test is the Wald test of the exposure coefficients of beta, the joint test that of g with them.
// Z_j' Omega Z_l has two parts: sum_i om_i e_ij e_il g_i^2 and, for the related samples, u_j' u_l with
// u_j = L^-1 (s_r o g_r o e_j) (S sparse triangular solves per variant with the factor of the null model). The
// n-long sums are entries of L1' G and L2' (G o G) per block: L1 = [Q o e_j | r o e_j], L2 = [om o e_j o e_l]
// (j <= l). A variant whose A_jj lost more than 1 - 1 / kKappaMax of Z_j' Omega Z_j to the subtraction is
// recomputed from Z - Vb Q' Z when the null model has Vb; without it (life-table Cox model) A is kept, and the
// variant is singular (status 2) below kSingular of Z_j' Omega Z_j.
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
#include "kin_sparse.h"
#include "la_dispatch.h"
#include "scan_block.h"
#include "resources.h"

namespace {

using plink4r::kScanRowChunk;
constexpr double kKappaMax = 1e4;        // explicit projection when g' Omega g / A exceeds this
constexpr double kSingular = 1e-9;       // without Vb: A below this share of g' Omega g is taken as 0

inline double dot(const double* a, const double* b, int len) {
  double s = 0.0;
  for (int i = 0; i < len; ++i) s += a[i] * b[i];
  return s;
}

inline int pair_index(int S, int j, int l) {
  if (j > l) std::swap(j, l);
  return j * S - j * (j - 1) / 2 + (l - j);
}

// The null model and the interaction columns, shared by all threads.
struct Pre {
  int n = 0, k = 0, S = 0, P = 0, m1 = 0, m2 = 0, nr = 0;
  const plink4r::KinNull* N = nullptr;
  std::vector<double> L1;      // n x (k S + S): Q o e_j (j = 0..S-1), then r o e_j
  std::vector<double> L2;      // n x P: om o e_j o e_l, j <= l
  std::vector<double> E;       // n x S: 1, e_1, ..., e_{S-1}
  Eigen::MatrixXd Er;          // nr x S: s_q e_j of the related samples, in the order of the factor
  const double* Vb = nullptr;  // n x k; nullptr when the null model has no sample-level basis of its columns
  const double* om = nullptr;  // n
};

void build_pre(Pre& pre, const plink4r::KinNull& NM, const Rcpp::NumericMatrix& Ex) {
  const int n = NM.n, k = NM.p, S = Ex.ncol() + 1;
  const std::size_t N = n;
  pre.n = n; pre.k = k; pre.S = S; pre.P = S * (S + 1) / 2; pre.m1 = k * S + S; pre.m2 = pre.P;
  pre.nr = NM.K.nr; pre.N = &NM;
  pre.Vb = NM.Vb.empty() ? nullptr : NM.Vb.data(); pre.om = NM.om.data();
  pre.E.resize(N * S);
  for (std::size_t i = 0; i < N; ++i) pre.E[i] = 1.0;
  for (int j = 1; j < S; ++j)
    for (std::size_t i = 0; i < N; ++i) pre.E[j * N + i] = Ex(i, j - 1);
  const double* Q = NM.L1.data();
  const double* r = Q + static_cast<std::size_t>(k) * N;
  pre.L1.resize(N * pre.m1);
  for (int j = 0; j < S; ++j) {
    const double* e = &pre.E[j * N];
    for (int c = 0; c < k; ++c) {
      double* d = &pre.L1[static_cast<std::size_t>(j * k + c) * N];
      const double* q = Q + c * N;
      for (std::size_t i = 0; i < N; ++i) d[i] = q[i] * e[i];
    }
    double* d = &pre.L1[static_cast<std::size_t>(k * S + j) * N];
    for (std::size_t i = 0; i < N; ++i) d[i] = r[i] * e[i];
  }
  pre.L2.resize(N * pre.m2);
  for (int j = 0; j < S; ++j)
    for (int l = j; l < S; ++l) {
      const double* ej = &pre.E[j * N];
      const double* el = &pre.E[l * N];
      double* d = &pre.L2[static_cast<std::size_t>(pair_index(S, j, l)) * N];
      for (std::size_t i = 0; i < N; ++i) d[i] = pre.om[i] * ej[i] * el[i];
    }
  pre.Er.resize(pre.nr, S);
  for (int q = 0; q < pre.nr; ++q)
    for (int j = 0; j < S; ++j) pre.Er(q, j) = NM.ss[q] * pre.E[j * N + NM.rs[q]];
}

// Variants of one thread: their compact genotypes (scan_block.h) and the block products.
struct Worker : plink4r::ScanBlock {
  std::vector<double> G, G2, C1, C2, z, Zt;
  std::vector<int> sing;
  Eigen::MatrixXd U;               // nr x S: the related samples' scaled columns of Z, then L^-1 of them
  void init(const Pre& pre, int capacity, bool dos) {
    ScanBlock::init(capacity, pre.n, dos);
    sing.assign(cap, 0);
    const std::size_t R = std::min(kScanRowChunk, n);
    G.assign(R * cap, 0.0); G2.assign(R * cap, 0.0);
    C1.assign(static_cast<std::size_t>(pre.m1) * cap, 0.0);
    C2.assign(static_cast<std::size_t>(pre.m2) * cap, 0.0);
    U.resize(pre.nr, pre.S);
  }
  // U' U of U = L^-1 [s_r o g_r o e_j] for slot s (S x S, the related samples' part of Z' Omega Z)
  void related(const Pre& pre, int s, Eigen::MatrixXd& R) {
    const int* rs = pre.N->rs.data();
    for (int q = 0; q < pre.nr; ++q) U.row(q) = value(s, rs[q]) * pre.Er.row(q);
    pre.N->llt.matrixL().solveInPlace(U);
    R.noalias() = U.transpose() * U;
  }
};

// A (S x S) of slot s from Zt = Z - Vb M (M = Q' Z, k x S): Zt' Omega Zt.
void explicit_projection(const Pre& pre, Worker& wk, int s, const double* M, double* A) {
  const int n = pre.n, k = pre.k, S = pre.S;
  const std::size_t N = n;
  if (wk.z.empty()) { wk.z.resize(N); wk.Zt.resize(N * S); }
  wk.centred(s, 0, n, wk.z.data());
  for (int j = 0; j < S; ++j) {
    double* zt = &wk.Zt[j * N];
    const double* e = &pre.E[j * N];
    for (std::size_t i = 0; i < N; ++i) zt[i] = wk.z[i] * e[i];
    for (int c = 0; c < k; ++c) {
      const double m = M[j * k + c];
      const double* v = pre.Vb + c * N;
      for (std::size_t i = 0; i < N; ++i) zt[i] -= v[i] * m;
    }
  }
  for (int j = 0; j < S; ++j)
    for (int l = j; l < S; ++l) {
      const double* a = &wk.Zt[j * N];
      const double* b = &wk.Zt[l * N];
      double sa = 0.0;
      for (std::size_t i = 0; i < N; ++i) sa += pre.om[i] * a[i] * b[i];
      A[j + l * S] = A[l + j * S] = sa;
    }
  if (pre.nr == 0) return;
  for (int q = 0; q < pre.nr; ++q)
    for (int j = 0; j < S; ++j) wk.U(q, j) = pre.N->ss[q] * wk.Zt[j * N + pre.N->rs[q]];
  pre.N->llt.matrixL().solveInPlace(wk.U);
  Eigen::Map<Eigen::MatrixXd> Am(A, S, S);
  Am.noalias() += wk.U.transpose() * wk.U;
}

// c (S) and A (S x S) of the nb slots of the worker; returns the number of explicit projections.
int block_stats(const Pre& pre, Worker& wk, int nb, double* cc, double* AA) {
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
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(S, S);
  std::vector<double> zoz(S);
  int nexp = 0;
  for (int s = 0; s < nb; ++s) {
    wk.sing[s] = 0;
    if (wk.mono[s]) continue;
    double* c = cc + static_cast<std::size_t>(s) * S;
    double* A = AA + static_cast<std::size_t>(s) * S * S;
    const double* M = &wk.C1[static_cast<std::size_t>(s) * m1];
    const double* c2 = &wk.C2[static_cast<std::size_t>(s) * m2];
    if (pre.nr > 0) wk.related(pre, s, R);
    bool ill = false;
    for (int j = 0; j < S; ++j) {
      for (int l = 0; l < S; ++l) A[j + l * S] = c2[pair_index(S, j, l)] + R(j, l) - dot(M + j * k, M + l * k, k);
      zoz[j] = c2[pair_index(S, j, j)] + R(j, j);
      if (!(A[j + j * S] > zoz[j] / kKappaMax)) ill = true;
    }
    for (int j = 0; j < S; ++j) c[j] = M[kS + j];
    if (!ill) continue;
    if (pre.Vb) { explicit_projection(pre, wk, s, M, A); ++nexp; continue; }
    for (int j = 0; j < S; ++j)
      if (!(A[j + j * S] > zoz[j] * kSingular)) wk.sing[s] = 1;
  }
  return nexp;
}

}  // namespace

// path: .pgen or .bed; bed = true for a BED (A1 = .bim column 5 counted), else PGEN (ALT counted).
// var_idx, sample_idx: 0-based file rows. null: external pointer to the plink4r::KinNull of the analysis samples
// (in the order of sample_idx), projected. E: n x (S - 1), the n_exp exposures first, then the interaction
// covariates (columns of the design of the null model, centred or not). block: variants per thread (0 = from
// mem_gb). Returns per variant n_obs, a1_sum, beta and se (columns of Z = [g, g e_1, ...]), beta_marginal,
// se_marginal (the association test of g), chisq_int (exposure terms), chisq_joint (g and the exposure terms),
// status (0 ok, 1 monomorphic / all missing, 2 singular).
// [[Rcpp::export]]
Rcpp::List gxe_kin_stream_cpp(const std::string& path, bool bed, int n_raw, Rcpp::IntegerVector allele_ct,
                             Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                             SEXP null, Rcpp::NumericMatrix E, int n_exp, int threads, double mem_gb, int block) {
  const int T = plink4r::resolve_threads(threads);
  mem_gb = plink4r::resolve_mem_gb(mem_gb, T);
  const Rcpp::XPtr<plink4r::KinNull> np(null);
  const plink4r::KinNull& NM = *np;
  const plink4r::Kin& K = NM.K;
  const int p = var_idx.size(), ns = sample_idx.size(), k = NM.p, S = E.ncol() + 1;
  if (NM.n != ns || NM.L1.size() != static_cast<std::size_t>(ns) * (k + 1) || E.nrow() != ns)
    Rcpp::stop("the null model does not match the samples");
  if (n_exp < 1 || n_exp > S - 1) Rcpp::stop("n_exp out of range");
  plink4r::ScanFile F(path, bed, n_raw, allele_ct.size(), std::vector<int>(allele_ct.begin(), allele_ct.end()),
                      std::vector<int>(var_idx.begin(), var_idx.end()), sample_idx, T);
  const bool dos = F.dosage;

  const std::size_t N = ns;
  Pre pre;
  build_pre(pre, NM, E);

  // block size from the memory budget
  const double shared = 8.0 * N * (pre.m1 + pre.m2 + S + k + 1.0) + 8.0 * p * (2.0 * S + 8.0) +
                        16.0 * (2.0 * K.pairs + K.nr) * 3 + 8.0 * K.nr * S;
  const double per_var = plink4r::ScanBlock::var_bytes(ns, dos) + 16.0 * std::min<int>(kScanRowChunk, ns) +
                         8.0 * (pre.m1 + pre.m2) + 8.0 * (S + S * S + 4);
  const double per_thread = 8.0 * N * (1 + S) + 8.0 * K.nr * S + F.thread_bytes();
  const int b = plink4r::scan_block_size(block, mem_gb, shared, per_thread, per_var, p, ns, T);

  Rcpp::NumericVector nobs(p), asum(p), beta_m(p, NA_REAL), se_m(p, NA_REAL), chi_int(p, NA_REAL), chi_joint(p, NA_REAL);
  Rcpp::NumericMatrix beta(p, S), se(p, S);
  std::fill(beta.begin(), beta.end(), NA_REAL); std::fill(se.begin(), se.end(), NA_REAL);
  Rcpp::IntegerVector status(p);
  std::vector<Worker> wk(T);
  std::vector<std::vector<double>> cbuf(T), abuf(T);
  for (int t = 0; t < T; ++t) {
    wk[t].init(pre, b, dos);
    cbuf[t].resize(static_cast<std::size_t>(b) * S);
    abuf[t].resize(static_cast<std::size_t>(b) * S * S);
  }
  std::vector<int> bad(T, -1);
  double* n_obs = nobs.begin();
  double* a1_sum = asum.begin();
  int nexp = 0;

  for (long long c0 = 0; c0 < p; c0 += static_cast<long long>(b) * T) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(T) reduction(+ : nexp)
#endif
    for (int t = 0; t < T; ++t) {
      const int th = plink4r::thread_id();
      Worker& W = wk[th];
      const long long j0 = c0 + static_cast<long long>(t) * b;
      if (j0 >= p) continue;
      const int nb = static_cast<int>(std::min<long long>(b, p - j0));
      const long long failed = W.read(F, j0, nb, th, n_obs, a1_sum);
      if (failed >= 0) bad[th] = static_cast<int>(failed);
      double* cc = cbuf[th].data();
      double* AA = abuf[th].data();
      nexp += block_stats(pre, W, nb, cc, AA);
      for (int s = 0; s < nb; ++s) {
        const long long j = j0 + s;
        if (W.mono[s]) { status[j] = 1; continue; }
        Eigen::Map<const Eigen::MatrixXd> A(AA + static_cast<std::size_t>(s) * S * S, S, S);
        Eigen::Map<const Eigen::VectorXd> c(cc + static_cast<std::size_t>(s) * S, S);
        if (W.sing[s] || !(A(0, 0) > 0)) { status[j] = 2; continue; }
        beta_m[j] = c[0] / A(0, 0);  // the association test of g (the model without the interaction terms)
        se_m[j] = std::sqrt(1.0 / A(0, 0));
        Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
        if (ldlt.info() != Eigen::Success || !ldlt.isPositive() || !(ldlt.vectorD().minCoeff() > 0)) { status[j] = 2; continue; }
        const Eigen::MatrixXd Ainv = ldlt.solve(Eigen::MatrixXd::Identity(S, S));
        const Eigen::VectorXd bt = Ainv * c;
        for (int q = 0; q < S; ++q) { beta(j, q) = bt[q]; se(j, q) = std::sqrt(Ainv(q, q)); }
        // Wald tests of the exposure interactions (terms 1..n_exp) and of g with them (0..n_exp)
        const Eigen::VectorXd bi = bt.segment(1, n_exp), bj = bt.head(n_exp + 1);
        chi_int[j] = bi.dot(Ainv.block(1, 1, n_exp, n_exp).ldlt().solve(bi));
        chi_joint[j] = bj.dot(Ainv.topLeftCorner(n_exp + 1, n_exp + 1).ldlt().solve(bj));
      }
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) Rcpp::stop("Genotype read failed at variant %d of %s", bad[t] + 1, path.c_str());
    Rcpp::checkUserInterrupt();
  }

  Rcpp::List out;
  out["n_obs"] = nobs; out["a1_sum"] = asum; out["status"] = status; out["beta"] = beta; out["se"] = se;
  out["beta_marginal"] = beta_m; out["se_marginal"] = se_m; out["chisq_int"] = chi_int; out["chisq_joint"] = chi_joint;
  out["explicit"] = nexp; out["block"] = b; out["threads"] = T; out["mem_gb"] = mem_gb; out["dosage"] = dos;
  out["simd"] = std::string(plink4r::la().isa); out["n_related"] = K.nr;
  out["pairs"] = static_cast<double>(K.pairs);
  return out;
}
