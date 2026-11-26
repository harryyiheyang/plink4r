// gwas_kin_stream_cpp: the kinship route of gwas(): mixed-model association of every biallelic variant of a PGEN or
// BED with a sparse kinship (fastGWA; Jiang et al. 2019), streamed in blocks; the genotype matrix is never formed
// and no table of the results either: the rows of a group of blocks are written to the result file as soon as they
// are computed, each one CHR, BP, REF, ALT of the variant as its .bim / .pvar line has them, then its statistics
// (ResultText, variant_text.h).
//
// Null model: the plink4r::KinNull left by reml_tau_gaussian.cpp / reml_tau_binomial.cpp / reml_tau_cox.cpp /
// cox_null.cpp (kin_sparse.h), used as
// it is: Omega = V^-1 = S H^-1 S, P = Omega - Omega X (X' Omega X)^-1 X' Omega, r = P y (y the response, or the
// working response of the logistic and Cox models), Q = Omega X La^-T and Vb = X La^-T with X' Omega X = La La'.
// For a variant g (centred, missing call = mean):
//   c = g' r,   A = g' P g = g' Omega g - |Q' g|^2,   beta = c / A,   se = A^-1/2,   chi2 = c^2 / A
// (generalised least squares; for the logistic model on the working response, one step from the null). The exact
// statistic, not the GRAMMAR-GAMMA approximation.
// g' Omega g has two parts. Samples without a relative: sum_i om_i g_i^2. Related samples: |L^-1 (s_r o g_r)|^2,
// one sparse triangular solve per variant with the Cholesky factor of H_rr found by the null model. The n-long
// sums of the first part, Q' g and g' r are entries of L1' G and L2' (G o G) per block of variants, L1 = [Q | r],
// L2 = [om] (0 at related samples), as in assoc_stream.cpp, from which the reading and blocking are taken. A
// variant whose A lost more than 1 - 1 / kKappaMax of g' Omega g to the subtraction is recomputed from
// g - Vb Q' g. The life-table Cox model has no Vb (its columns vary within a sample): there A is kept, and taken
// as 0 (the variant has no statistics) below kSingular of g' Omega g.
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
#include "p_value.h"
#include "scan_block.h"
#include "resources.h"

namespace {

using plink4r::kScanRowChunk;
constexpr double kKappaMax = 1e4;        // explicit projection when g' Omega g / A exceeds this
constexpr double kSingular = 1e-9;       // without Vb: A below this share of g' Omega g is taken as 0
constexpr double kRowBytes = 512.0;      // text of a variant while its block is written, for the memory budget

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
  const double* Vb = nullptr;  // n x k; nullptr when the null model has no sample-level basis of its columns
  const double* om = nullptr;  // n
};

// Variants of one thread: their compact genotypes (scan_block.h) and the block products.
struct Worker : plink4r::ScanBlock {
  std::vector<double> G, G2, C1, C2, z;
  Eigen::VectorXd gr;
  void init(const Pre& pre, int capacity, bool dos) {
    ScanBlock::init(capacity, pre.n, dos);
    const std::size_t R = std::min(kScanRowChunk, n);
    G.assign(R * cap, 0.0); G2.assign(R * cap, 0.0);
    C1.assign(static_cast<std::size_t>(pre.m1) * cap, 0.0);
    C2.assign(cap, 0.0);
    gr.resize(pre.nr);
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
    if (!(A > gog / kKappaMax)) {
      if (pre.Vb) { A = explicit_projection(pre, wk, s, M); ++nexp; }
      else if (!(A > gog * kSingular)) A = 0.0;
    }
    cc[s] = M[k];
    AA[s] = A;
  }
  return nexp;
}

}  // namespace

// path: .pgen or .bed; bed = true for a BED (A1 = .bim column 5 counted), else PGEN (ALT counted). variant_text: its
// .bim / .pvar. out: the text file to write: after CHR, BP, REF, ALT of each line of the .bim / .pvar ALT_FREQ, N,
// BETA, SE, the statistic BETA / SE under the name `statistic` (T or Z) and its normal P; NA for a variant without
// statistics (no variation among the samples, or A <= 0) and for a line with more than two alleles, which is not
// scanned. sample_idx: 0-based file rows of the analysed samples. null: external pointer to the plink4r::KinNull of
// the analysis samples (in the order of sample_idx), projected. block: variants per thread (0 = from mem_gb).
// [[Rcpp::export]]
Rcpp::List gwas_kin_stream_cpp(const std::string& path, bool bed, int n_raw, const std::string& variant_text,
                              const std::string& out, const std::string& statistic, Rcpp::IntegerVector sample_idx,
                              SEXP null, int threads, double mem_gb, int block) {
  const int T = plink4r::resolve_threads(threads);
  mem_gb = plink4r::resolve_mem_gb(mem_gb, T);
  const Rcpp::XPtr<plink4r::KinNull> np(null);
  const plink4r::KinNull& NM = *np;
  const plink4r::Kin& K = NM.K;
  const int ns = sample_idx.size(), k = NM.p;
  if (NM.n != ns || NM.L1.size() != static_cast<std::size_t>(ns) * (k + 1))
    Rcpp::stop("the null model does not match the samples");
  const plink4r::FileVariants fv = plink4r::file_variants(path, bed, n_raw, variant_text);
  plink4r::ScanFile F(path, bed, n_raw, fv.m, fv.allele_ct, fv.rows, sample_idx, T);
  const bool dos = F.dosage;
  const int p = static_cast<int>(F.vi.size());

  const std::size_t N = ns;
  Pre pre;
  pre.n = ns; pre.k = k; pre.m1 = k + 1; pre.nr = K.nr; pre.N = &NM;
  pre.L1 = NM.L1.data(); pre.Vb = NM.Vb.empty() ? nullptr : NM.Vb.data(); pre.om = NM.om.data();

  // block size from the memory budget
  const double shared = 8.0 * N * (pre.m1 + k + 1.0) + 16.0 * (2.0 * K.pairs + K.nr) * 3;
  const double per_var = plink4r::ScanBlock::var_bytes(ns, dos) + 16.0 * std::min<int>(kScanRowChunk, ns) +
                         8.0 * (pre.m1 + 1) + 8.0 * 4 + kRowBytes;
  const double per_thread = 8.0 * N + 8.0 * K.nr + F.thread_bytes();
  const int b = plink4r::scan_block_size(block, mem_gb, shared, per_thread, per_var, p, ns, T);

  plink4r::ResultText text(variant_text, bed, out, {"ALT_FREQ", "N", "BETA", "SE", statistic, "P"});
  std::vector<Worker> wk(T);
  std::vector<std::vector<double>> cbuf(T), abuf(T), nbuf(T), sbuf(T);
  std::vector<std::string> rows(T);
  for (int t = 0; t < T; ++t) {
    wk[t].init(pre, b, dos);
    cbuf[t].resize(b);
    abuf[t].resize(b);
    nbuf[t].resize(b); sbuf[t].resize(b);
  }
  std::vector<int> bad(T, -1);
  int nexp = 0;
  long long n_mono = 0, n_sing = 0;

  for (long long c0 = 0; c0 < p; c0 += static_cast<long long>(b) * T) {
    const long long c1 = std::min<long long>(p, c0 + static_cast<long long>(b) * T);
    text.load(static_cast<long long>(F.vi[c1 - 1]) + 1);
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
      nexp += block_stats(pre, W, nb, cc, AA);
      long long r = j0 == c0 ? text.first() : static_cast<long long>(F.vi[j0 - 1]) + 1;  // the next file row to write
      for (int s = 0; s < nb; ++s) {
        const long long j = j0 + s;
        for (; r < F.vi[j]; ++r) text.row(o, r, 0, nullptr);  // lines with more than two alleles: not scanned
        double v[5] = {a1_sum[s] / (2 * n_obs[s]), NA_REAL, NA_REAL, NA_REAL, NA_REAL};  // ALT_FREQ BETA SE T|Z P
        if (W.mono[s]) ++n_mono;
        else if (!(AA[s] > 0)) ++n_sing;
        else {
          v[1] = cc[s] / AA[s];
          v[2] = std::sqrt(1.0 / AA[s]);
          v[3] = v[1] / v[2];
          v[4] = plink4r::p_z(v[3]);
        }
        text.row(o, r, static_cast<long long>(n_obs[s]), v);
        ++r;
      }
    }
    for (int t = 0; t < T; ++t)
      if (bad[t] >= 0) Rcpp::stop("Genotype read failed at variant %d of %s", bad[t] + 1, path.c_str());
    for (int t = 0; t < T; ++t) text.write(rows[t]);
    Rcpp::checkUserInterrupt();
  }
  text.finish(static_cast<long long>(fv.m));

  Rcpp::List res;
  res["variants"] = static_cast<double>(fv.m); res["tested"] = p; res["bytes"] = text.bytes();
  res["monomorphic"] = static_cast<double>(n_mono); res["singular"] = static_cast<double>(n_sing);
  res["explicit"] = nexp; res["block"] = b; res["threads"] = T; res["mem_gb"] = mem_gb; res["dosage"] = dos;
  res["simd"] = std::string(plink4r::la().isa); res["n_related"] = K.nr;
  res["pairs"] = static_cast<double>(K.pairs);
  return res;
}
