// gwas_kin_stream_cpp: the kinship route of gwas(): mixed-model association of every selected variant of a PGEN or
// BED with a sparse kinship
// (fastGWA; Jiang et al. 2019), streamed in blocks; the genotype matrix is never formed.
//
// Null model: the plink4r::KinNull left by reml_tau_gaussian.cpp / reml_tau_binomial.cpp (kin_sparse.h), used as
// it is: Omega = V^-1 = S H^-1 S, P = Omega - Omega X (X' Omega X)^-1 X' Omega, r = P y (y the response, or the
// working response of the logistic model), Q = Omega X La^-T and Vb = X La^-T with X' Omega X = La La'. For a
// variant g (centred, missing call = mean):
//   c = g' r,   A = g' P g = g' Omega g - |Q' g|^2,   beta = c / A,   se = A^-1/2,   chi2 = c^2 / A
// (generalised least squares; for the logistic model on the working response, one step from the null). The exact
// statistic, not the GRAMMAR-GAMMA approximation.
// g' Omega g has two parts. Samples without a relative: sum_i om_i g_i^2. Related samples: |L^-1 (s_r o g_r)|^2,
// one sparse triangular solve per variant with the Cholesky factor of H_rr found by the null model. The n-long
// sums of the first part, Q' g and g' r are entries of L1' G and L2' (G o G) per block of variants, L1 = [Q | r],
// L2 = [om] (0 at related samples), as in assoc_stream.cpp, from which the reading and blocking are taken. A
// variant whose A lost more than 1 - 1 / kKappaMax of g' Omega g to the subtraction is recomputed from
// g - Vb Q' g.
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
constexpr std::uint16_t kDosMiss = plink4r::PgenFile::kMissing;

inline double dot(const double* a, const double* b, int len) {
  double s = 0.0;
  for (int i = 0; i < len; ++i) s += a[i] * b[i];
  return s;
}

// The null model, shared by all threads.
struct Pre {
  int n = 0, k = 0, m1 = 0, nr = 0;
  const plink4r::KinNull* N = nullptr;
  const double* L1 = nullptr;  // n x (k + 1): Q, r
  const double* Vb = nullptr;  // n x k
  const double* om = nullptr;  // n
};

// Variants of one thread: compact genotypes over the analysis samples and the block products.
struct Worker {
  int cap = 0, n = 0;
  bool dosage = false;
  std::size_t sb = 0;
  std::vector<unsigned char> g8;   // cap x sb: 2-bit codes (0/1/2 = count, 3 = missing)
  std::vector<std::uint16_t> g16;  // cap x n: dosage x 16384, kDosMiss = missing
  std::vector<double> mean, G, G2, C1, C2, z;
  std::vector<int> mono;
  Eigen::VectorXd gr;
  void init(const Pre& pre, int capacity, bool dos) {
    cap = capacity; n = pre.n; dosage = dos;
    sb = (static_cast<std::size_t>(n) + 3) / 4;
    if (dosage) g16.assign(static_cast<std::size_t>(cap) * n, 0); else g8.assign(cap * sb, 0);
    mean.assign(cap, 0.0); mono.assign(cap, 0);
    const std::size_t R = std::min(kRowChunk, n);
    G.assign(R * cap, 0.0); G2.assign(R * cap, 0.0);
    C1.assign(static_cast<std::size_t>(pre.m1) * cap, 0.0);
    C2.assign(cap, 0.0);
    gr.resize(pre.nr);
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
  // centred value of sample i of slot s
  double value(int s, int i) const {
    const double mu = mean[s];
    if (dosage) {
      const std::uint16_t v = g16[static_cast<std::size_t>(s) * n + i];
      return v == kDosMiss ? 0.0 : v * (1.0 / 16384.0) - mu;
    }
    const unsigned code = (g8[s * sb + (i >> 2)] >> ((i & 3) << 1)) & 3;
    return code == 3 ? 0.0 : static_cast<double>(code) - mu;
  }
  // the related samples' centred values of slot s, scaled, into gr (in the order of the factor)
  void related(const Pre& pre, int s) {
    const int* rs = pre.N->rs.data();
    const double* ss = pre.N->ss.data();
    for (int q = 0; q < pre.nr; ++q) gr[q] = ss[q] * value(s, rs[q]);
  }
};

// A of slot s from gt = g - Vb M (M = Q' g): gt' Omega gt.
double explicit_projection(const Pre& pre, Worker& wk, int s, const double* M) {
  const int n = pre.n, k = pre.k;
  const std::size_t N = n;
  if (wk.z.empty()) wk.z.resize(N);
  wk.centred(s, 0, n, wk.z.data());
  for (int c = 0; c < k; ++c) {
    const double m = M[c];
    const double* v = pre.Vb + c * N;
    for (std::size_t i = 0; i < N; ++i) wk.z[i] -= v[i] * m;
  }
  double a = 0.0;
  for (std::size_t i = 0; i < N; ++i) a += pre.om[i] * wk.z[i] * wk.z[i];
  if (pre.nr == 0) return a;
  for (int q = 0; q < pre.nr; ++q) wk.gr[q] = pre.N->ss[q] * wk.z[pre.N->rs[q]];
  return a + pre.N->quad_rel(wk.gr);
}

// c and A of the nb slots of the worker; returns the number of explicit projections.
int block_stats(const Pre& pre, Worker& wk, int nb, double* cc, double* AA) {
  const int n = pre.n, k = pre.k, m1 = pre.m1;
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
    LA.dgemm_tn(pre.om + row0, rows, 1, n, wk.G2.data(), nb, rows, wk.C2.data(), 1, beta);
  }
  int nexp = 0;
  for (int s = 0; s < nb; ++s) {
    if (wk.mono[s]) continue;
    const double* M = &wk.C1[static_cast<std::size_t>(s) * m1];
    double gog = wk.C2[s];
    if (pre.nr > 0) {
      wk.related(pre, s);
      gog += pre.N->quad_rel(wk.gr);
    }
    double A = gog - dot(M, M, k);
    if (!(A > gog / kKappaMax)) { A = explicit_projection(pre, wk, s, M); ++nexp; }
    cc[s] = M[k];
    AA[s] = A;
  }
  return nexp;
}

}  // namespace

// path: .pgen or .bed; bed = true for a BED (A1 = .bim column 5 counted), else PGEN (ALT counted).
// var_idx, sample_idx: 0-based file rows. null: external pointer to the plink4r::KinNull of the analysis samples
// (in the order of sample_idx), projected. block: variants per thread (0 = from mem_gb).
// Returns per variant n_obs, a1_sum, beta, se, status (0 ok, 1 monomorphic / all missing, 2 A <= 0).
// [[Rcpp::export]]
Rcpp::List gwas_kin_stream_cpp(const std::string& path, bool bed, int n_raw, Rcpp::IntegerVector allele_ct,
                              Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                              SEXP null, int threads, double mem_gb, int block) {
  const int T = plink4r::resolve_threads(threads);
  mem_gb = plink4r::resolve_mem_gb(mem_gb, T);
  const Rcpp::XPtr<plink4r::KinNull> np(null);
  const plink4r::KinNull& NM = *np;
  const plink4r::Kin& K = NM.K;
  const int p = var_idx.size(), ns = sample_idx.size(), k = NM.p;
  if (NM.n != ns || NM.L1.size() != static_cast<std::size_t>(ns) * (k + 1))
    Rcpp::stop("the null model does not match the samples");
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
  pre.n = ns; pre.k = k; pre.m1 = k + 1; pre.nr = K.nr; pre.N = &NM;
  pre.L1 = NM.L1.data(); pre.Vb = NM.Vb.data(); pre.om = NM.om.data();

  // block size from the memory budget
  const double GiB = 1073741824.0;
  const double shared = 8.0 * N * (pre.m1 + k + 1.0) + 8.0 * p * 6.0 + 16.0 * (2.0 * K.pairs + K.nr) * 3;
  const double per_var = (dos ? 2.0 * N : (N + 3) / 4.0) + 16.0 * std::min<int>(kRowChunk, ns) +
                         8.0 * (pre.m1 + 1) + 8.0 * 4;
  const double per_thread = 8.0 * N + 8.0 * K.nr + (bed ? rawb : (dos ? 3.0 * n_raw : rawb)) + 1048576.0;
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

  Rcpp::NumericVector nobs(p), asum(p), beta(p, NA_REAL), se(p, NA_REAL);
  Rcpp::IntegerVector status(p);
  std::vector<Worker> wk(T);
  std::vector<std::vector<double>> cbuf(T), abuf(T);
  for (int t = 0; t < T; ++t) {
    wk[t].init(pre, b, dos);
    cbuf[t].resize(b);
    abuf[t].resize(b);
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
            o[i] = v;
            if (v == kDosMiss) continue;
            const double x = v * (1.0 / 16384.0);
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
        if (!(AA[s] > 0)) { status[j] = 2; continue; }
        beta[j] = cc[s] / AA[s];
        se[j] = std::sqrt(1.0 / AA[s]);
      }
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) Rcpp::stop("Genotype read failed at variant %d of %s", bad[t] + 1, path.c_str());
    Rcpp::checkUserInterrupt();
  }

  Rcpp::List out;
  out["n_obs"] = nobs; out["a1_sum"] = asum; out["status"] = status; out["beta"] = beta; out["se"] = se;
  out["explicit"] = nexp; out["block"] = b; out["threads"] = T; out["mem_gb"] = mem_gb; out["dosage"] = dos;
  out["simd"] = std::string(plink4r::la().isa); out["n_related"] = K.nr;
  out["pairs"] = static_cast<double>(K.pairs);
  return out;
}
