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
#include <memory>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "kin_sparse.h"
#include "la_dispatch.h"
#include "pgen_reader.h"
#include "resources.h"

namespace {

constexpr int kRowChunk = 2048;          // rows of G expanded at a time (multiple of 4)
constexpr int kBlockMax = 4096;          // largest automatic block per thread
constexpr double kKappaMax = 1e4;        // explicit projection when g' Omega g / A exceeds this
constexpr double kSingular = 1e-9;       // without Vb: A below this share of g' Omega g is taken as 0
constexpr std::uint16_t kDosMiss = plink4r::PgenFile::kMissing;

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

// Variants of one thread: compact genotypes over the analysis samples and the block products.
struct Worker {
  int cap = 0, n = 0;
  bool dosage = false;
  std::size_t sb = 0;
  std::vector<unsigned char> g8;   // cap x sb: 2-bit codes (0/1/2 = count, 3 = missing)
  std::vector<std::uint16_t> g16;  // cap x n: dosage in steps of 0.01 (dosage_grid), kDosMiss = missing
  std::vector<double> mean, G, G2, C1, C2, z, Zt;
  std::vector<int> mono, sing;
  Eigen::MatrixXd U;               // nr x S: the related samples' scaled columns of Z, then L^-1 of them
  void init(const Pre& pre, int capacity, bool dos) {
    cap = capacity; n = pre.n; dosage = dos;
    sb = (static_cast<std::size_t>(n) + 3) / 4;
    if (dosage) g16.assign(static_cast<std::size_t>(cap) * n, 0); else g8.assign(cap * sb, 0);
    mean.assign(cap, 0.0); mono.assign(cap, 0); sing.assign(cap, 0);
    const std::size_t R = std::min(kRowChunk, n);
    G.assign(R * cap, 0.0); G2.assign(R * cap, 0.0);
    C1.assign(static_cast<std::size_t>(pre.m1) * cap, 0.0);
    C2.assign(static_cast<std::size_t>(pre.m2) * cap, 0.0);
    U.resize(pre.nr, pre.S);
  }
  void centred(int s, int row0, int rows, double* out) const {
    const double mu = mean[s];
    if (dosage) {
      const std::uint16_t* d = &g16[static_cast<std::size_t>(s) * n + row0];
      for (int i = 0; i < rows; ++i) out[i] = d[i] == kDosMiss ? 0.0 : d[i] * 0.01 - mu;
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
  // centred value of sample i of slot s
  double value(int s, int i) const {
    const double mu = mean[s];
    if (dosage) {
      const std::uint16_t v = g16[static_cast<std::size_t>(s) * n + i];
      return v == kDosMiss ? 0.0 : v * 0.01 - mu;
    }
    const unsigned code = (g8[s * sb + (i >> 2)] >> ((i & 3) << 1)) & 3;
    return code == 3 ? 0.0 : static_cast<double>(code) - mu;
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
  const std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  const bool ident = plink4r::identity_samples(si.data(), ns, n_raw);

  std::unique_ptr<plink4r::PgenFile> pf;
  std::unique_ptr<plink4r::BedFile> bf;
  if (bed) bf.reset(new plink4r::BedFile(path, n_raw, allele_ct.size(), T));  // allele_ct: one entry per .bim row
  else pf.reset(new plink4r::PgenFile(path, n_raw, std::vector<int>(allele_ct.begin(), allele_ct.end()), T));
  const bool dos = pf && pf->dosage();
  const std::size_t rawb = (static_cast<std::size_t>(n_raw) + 3) / 4;

  const std::size_t N = ns;
  Pre pre;
  build_pre(pre, NM, E);

  // block size from the memory budget
  const double GiB = 1073741824.0;
  const double shared = 8.0 * N * (pre.m1 + pre.m2 + S + k + 1.0) + 8.0 * p * (2.0 * S + 8.0) +
                        16.0 * (2.0 * K.pairs + K.nr) * 3 + 8.0 * K.nr * S;
  const double per_var = (dos ? 2.0 * N : (N + 3) / 4.0) + 16.0 * std::min<int>(kRowChunk, ns) +
                         8.0 * (pre.m1 + pre.m2) + 8.0 * (S + S * S + 4);
  const double per_thread = 8.0 * N * (1 + S) + 8.0 * K.nr * S + (bed ? rawb : (dos ? 3.0 * n_raw : rawb)) + 1048576.0;
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
      // read and compact
      for (int s = 0; s < nb; ++s) {
        const long long j = j0 + s;
        double no = 0, sum = 0, sq = 0;
        if (dos) {
          const std::uint16_t* dd = pf->read_dosage(static_cast<std::size_t>(vi[j]), th);
          if (!dd) { bad[th] = static_cast<int>(j); W.mono[s] = 1; continue; }
          std::uint16_t* o = &W.g16[static_cast<std::size_t>(s) * ns];
          for (int i = 0; i < ns; ++i) {
            const std::uint16_t v = dd[si[i]];
            if (v == kDosMiss) { o[i] = v; continue; }
            const unsigned int c = plink4r::dosage_grid(v);
            o[i] = static_cast<std::uint16_t>(c);
            const double x = c * 0.01;
            no += 1; sum += x; sq += x * x;
          }
        } else {
          const unsigned char* gv;
          gv = bed ? bf->read_counts(static_cast<std::size_t>(vi[j]), th) : pf->read(static_cast<std::size_t>(vi[j]), th);
          if (!gv) { bad[th] = static_cast<int>(j); W.mono[s] = 1; continue; }
          unsigned char* dst = &W.g8[s * W.sb];
          plink4r::compact_genovec(gv, si.data(), ns, ident, dst);
          long long cnt[4];
          plink4r::count_codes(dst, W.sb, cnt);
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
