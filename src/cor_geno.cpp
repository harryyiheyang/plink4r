// LD (Pearson correlation) matrix from a PLINK2 .pgen (hardcall or dosage) or a PLINK1 .bed.
// Genotypes are held as uint16 fixed point (ALT/A1 dosage x 16384, exact for pgen dosages; 65535 = missing).
// Missing cells: anchor AF given -> scheme 3 constant c_j = (2 n p_j - sum_obs g)/m_j clipped to [0,2];
// no AF -> observed mean. Then r_jl = sum_i (g_ij - mu_j)(g_il - mu_l) / (n sd_j sd_l), accumulated over
// person chunks as float32 GEMM into a double matrix (OpenMP via Eigen).
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pgen_reader.h"
#include "resources.h"

using Eigen::MatrixXd;
using Eigen::MatrixXf;

namespace {
inline int tidc() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}
inline double secsc(std::chrono::steady_clock::time_point a) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
}

Rcpp::List cor_from_G(const std::vector<std::uint16_t>& G, int ns, int p, const std::vector<double>& af, int T, double t_read) {
  const auto t0 = std::chrono::steady_clock::now();
  const double unit = 1.0 / 16384.0;
  std::vector<double> mu(p), sd(p), fill(p);
  std::vector<int> nmiss(p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
  for (int j = 0; j < p; ++j) {
    const std::uint16_t* g = &G[static_cast<std::size_t>(j) * ns];
    double s = 0, s2 = 0; int m = 0;
    for (int i = 0; i < ns; ++i) {
      if (g[i] == plink4r::PgenFile::kMissing) { ++m; continue; }
      const double v = g[i] * unit; s += v; s2 += v * v;
    }
    double c = (ns - m) > 0 ? s / (ns - m) : 0.0;
    if (m > 0 && std::isfinite(af[j])) c = std::min(2.0, std::max(0.0, (2.0 * ns * af[j] - s) / m));
    const double tot = s + m * c, tot2 = s2 + m * c * c;
    mu[j] = tot / ns;
    sd[j] = std::sqrt(std::max(0.0, tot2 / ns - mu[j] * mu[j]));
    fill[j] = c; nmiss[j] = m;
  }
  MatrixXd XtX = MatrixXd::Zero(p, p);
  const int CH = 16384;
  // X'X by 256-column tiles (lower triangle) shared out by OpenMP; each tile product is a serial Eigen call.
  const int TB = 256, nt = (p + TB - 1) / TB;
  std::vector<std::pair<int, int>> tiles;
  for (int b = 0; b < nt; ++b) for (int a2 = b; a2 < nt; ++a2) tiles.emplace_back(a2, b);
  MatrixXf Z(CH, p);
  for (int a = 0; a < ns; a += CH) {
    const int nb = std::min(CH, ns - a);
    if (nb != Z.rows()) Z.resize(nb, p);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
    for (int j = 0; j < p; ++j) {
      const std::uint16_t* g = &G[static_cast<std::size_t>(j) * ns + a];
      const float f = static_cast<float>(fill[j] - mu[j]), m0 = static_cast<float>(mu[j]);
      for (int i = 0; i < nb; ++i) Z(i, j) = g[i] == plink4r::PgenFile::kMissing ? f : static_cast<float>(g[i] * unit) - m0;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
    for (int t = 0; t < static_cast<int>(tiles.size()); ++t) {
      const int r0 = tiles[t].first * TB, c0 = tiles[t].second * TB;
      const int nr = std::min(TB, p - r0), nc = std::min(TB, p - c0);
      const MatrixXf P = Z.middleCols(r0, nr).transpose() * Z.middleCols(c0, nc);
      XtX.block(r0, c0, nr, nc) += P.cast<double>();
    }
  }
  for (int b = 0; b < nt; ++b)  // mirror the strictly upper tiles
    for (int a2 = b + 1; a2 < nt; ++a2) {
      const int r0 = a2 * TB, c0 = b * TB, nr = std::min(TB, p - r0), nc = std::min(TB, p - c0);
      XtX.block(c0, r0, nc, nr) = XtX.block(r0, c0, nr, nc).transpose();
    }
  Rcpp::NumericMatrix R(p, p);
  for (int l = 0; l < p; ++l)
    for (int j = 0; j < p; ++j) {
      const double d = ns * sd[j] * sd[l];
      R(j, l) = d > 0 ? XtX(j, l) / d : (j == l ? 1.0 : 0.0);
    }
  return Rcpp::List::create(Rcpp::_["R"] = R, Rcpp::_["mean"] = Rcpp::wrap(mu), Rcpp::_["sd"] = Rcpp::wrap(sd),
                            Rcpp::_["n_miss"] = Rcpp::wrap(nmiss), Rcpp::_["fill"] = Rcpp::wrap(fill), Rcpp::_["n"] = ns,
                            Rcpp::_["time"] = Rcpp::NumericVector::create(Rcpp::_["read"] = t_read, Rcpp::_["cor"] = secsc(t0)));
}
}  // namespace

// var_idx / sample_idx 0-based; af = anchor ALT frequency per SNP (NA = mean fill).
// [[Rcpp::export]]
Rcpp::List cor_pgen_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct, Rcpp::IntegerVector var_idx,
                        Rcpp::IntegerVector sample_idx, Rcpp::NumericVector af, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const auto t0 = std::chrono::steady_clock::now();
  const int p = var_idx.size(), ns = sample_idx.size();
  std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  const bool dos = pf.dosage();
  std::vector<std::uint16_t> G(static_cast<std::size_t>(p) * ns);
  std::vector<int> bad(T, -1);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
  for (int j = 0; j < p; ++j) {
    std::uint16_t* dst = &G[static_cast<std::size_t>(j) * ns];
    if (dos) {
      const std::uint16_t* d = pf.read_dosage(static_cast<std::size_t>(vi[j]), tidc());
      if (!d) { bad[tidc()] = j; continue; }
      for (int i = 0; i < ns; ++i) dst[i] = d[si[i]];
    } else {
      const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[j]), tidc());
      if (!gv) { bad[tidc()] = j; continue; }
      static const std::uint16_t code[4] = {0, 16384, 32768, plink4r::PgenFile::kMissing};
      for (int i = 0; i < ns; ++i) { const int r = si[i]; dst[i] = code[(gv[r >> 2] >> (2 * (r & 3))) & 3u]; }
    }
  }
  for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop("PGEN read failed at SNP %d", bad[t] + 1);
  std::vector<double> afv(af.begin(), af.end());
  Rcpp::List out = cor_from_G(G, ns, p, afv, T, secsc(t0));
  out["dosage"] = dos;
  return out;
}

// PLINK1 .bed (SNP-major): counts A1 (bim column 5): 00 -> 2, 10 -> 1, 11 -> 0, 01 -> missing.
// [[Rcpp::export]]
Rcpp::List cor_bed_cpp(const std::string& bed, int n_raw, Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                       Rcpp::NumericVector af, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const auto t0 = std::chrono::steady_clock::now();
  const int p = var_idx.size(), ns = sample_idx.size();
  const std::size_t bpv = (static_cast<std::size_t>(n_raw) + 3) / 4;
  std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  {
    std::ifstream f(bed.c_str(), std::ios::binary);
    unsigned char mg[3];
    if (!f.read(reinterpret_cast<char*>(mg), 3) || mg[0] != 0x6c || mg[1] != 0x1b || mg[2] != 0x01) Rcpp::stop("not a SNP-major PLINK1 .bed: " + bed);
  }
  std::vector<std::uint16_t> G(static_cast<std::size_t>(p) * ns);
  std::vector<int> bad(T, -1);
  static const std::uint16_t code[4] = {32768, plink4r::PgenFile::kMissing, 16384, 0};
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
  {
    std::ifstream f(bed.c_str(), std::ios::binary);
    std::vector<unsigned char> buf(bpv);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int j = 0; j < p; ++j) {
      f.seekg(static_cast<std::streamoff>(3 + static_cast<std::size_t>(vi[j]) * bpv));
      if (!f.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(bpv))) { bad[tidc()] = j; f.clear(); continue; }
      std::uint16_t* dst = &G[static_cast<std::size_t>(j) * ns];
      for (int i = 0; i < ns; ++i) { const int r = si[i]; dst[i] = code[(buf[r >> 2] >> (2 * (r & 3))) & 3u]; }
    }
  }
  for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop(".bed read failed at SNP %d", bad[t] + 1);
  std::vector<double> afv(af.begin(), af.end());
  return cor_from_G(G, ns, p, afv, T, secsc(t0));
}

// Per-variant observed count and ALT-dosage sum over the selected samples (dosage-aware; missing skipped).
// [[Rcpp::export]]
Rcpp::List pgen_stats_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct, Rcpp::IntegerVector var_idx,
                          Rcpp::IntegerVector sample_idx, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const int p = var_idx.size(), ns = sample_idx.size();
  std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  const bool dos = pf.dosage();
  std::vector<double> nobs(p), sum(p);
  std::vector<int> bad(T, -1);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(T)
#endif
  for (int j = 0; j < p; ++j) {
    long long c = 0; double s = 0;
    if (dos) {
      const std::uint16_t* d = pf.read_dosage(static_cast<std::size_t>(vi[j]), tidc());
      if (!d) { bad[tidc()] = j; continue; }
      for (int i = 0; i < ns; ++i) { const std::uint16_t v = d[si[i]]; if (v != plink4r::PgenFile::kMissing) { ++c; s += v; } }
      s /= 16384.0;
    } else {
      const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[j]), tidc());
      if (!gv) { bad[tidc()] = j; continue; }
      for (int i = 0; i < ns; ++i) { const int r = si[i]; const unsigned int x = (gv[r >> 2] >> (2 * (r & 3))) & 3u; if (x != 3u) { ++c; s += x; } }
    }
    nobs[j] = static_cast<double>(c); sum[j] = s;
  }
  for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop("PGEN read failed at SNP %d", bad[t] + 1);
  return Rcpp::List::create(Rcpp::_["n_obs"] = Rcpp::wrap(nobs), Rcpp::_["sum"] = Rcpp::wrap(sum), Rcpp::_["dosage"] = dos);
}
