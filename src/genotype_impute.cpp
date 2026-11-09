// Genotype imputation of hard-call missing cells from an eigen.bin LD block.
// Stage 1 (per person i):
//   z~ = (g - 2p)/sqrt(2p(1-p)) with anchor AF p, 0 at missing;  C = U W U',  W = lam/(lam+delta)
//   z_M = (I - C_MM)^{-1} C_{M,.} z~            ( = R_MO (R_OO + delta I)^{-1} z_O )
//   delta = 0: generalized inverse of (I - C_MM) (eigenvalues <= 1e-6 * max dropped)
// Stage 2 (per SNP j with missing rate >= anchor_min_miss; < 0 = never): shift the missing cells so that
//   sum_obs g + sum_miss g_hat = 2 n p_j, clip [0,2]. Below the cutoff stage 1 is kept: the anchor error is
//   amplified by 2/missing-rate.
//
// Layout for speed (exact, same result as the per-person dot product):
//   with a_lj = C_lj / sd_l and A_j = -sum_l a_lj mu_l (absent rows: a = 0),
//   rhs_ij = C_{.j}' z~_i = A_j + sum_{l in M_i} a_lj mu_l + sum_{l het_i} a_lj + 2 sum_{l homalt_i} a_lj.
//   People are processed in blocks; inside a block the loop runs over SNP columns j (column a_.j / C_.j hot in
//   cache) and over the people missing j, using 1-bit planes (het / homalt) and an 8-bit lookup table, so C is
//   streamed from memory once per person block instead of once per (person, missing SNP). C_MM is gathered
//   from the same hot column.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <vector>
#include <sys/resource.h>
#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"
#include "pgen_reader.h"

using Eigen::MatrixXf;
using Eigen::VectorXf;

namespace {
inline double secs(std::chrono::steady_clock::time_point a) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
}
inline double cpu_secs() {  // process user + system CPU time
  rusage u;
  getrusage(RUSAGE_SELF, &u);
  return u.ru_utime.tv_sec + u.ru_stime.tv_sec + 1e-6 * (u.ru_utime.tv_usec + u.ru_stime.tv_usec);
}
inline int tid() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}
// sum_l col[l] * (het_l + 2 hom_l) over the bit-planes of one person (m padded to a multiple of 16 with zeros)
#if defined(__AVX512F__)
inline float plane_dot(const float* col, const uint8_t* het, const uint8_t* hom, int nbyte) {
  __m512 acc = _mm512_setzero_ps();
  const uint16_t* h16 = reinterpret_cast<const uint16_t*>(het);
  const uint16_t* o16 = reinterpret_cast<const uint16_t*>(hom);
  for (int q = 0; q < nbyte / 2; ++q) {
    const __mmask16 h = h16[q], o = o16[q];
    if ((h | o) == 0) continue;
    const __m512 c = _mm512_loadu_ps(col + 16 * q);
    acc = _mm512_mask_add_ps(acc, h, acc, c);
    acc = _mm512_mask_add_ps(acc, o, acc, _mm512_add_ps(c, c));
  }
  return _mm512_reduce_add_ps(acc);
}

inline float plane_dot_words(const float* col, const uint8_t* het, const uint8_t* hom, const int* wl, int nwl) {
  __m512 acc = _mm512_setzero_ps();
  const uint16_t* h16 = reinterpret_cast<const uint16_t*>(het);
  const uint16_t* o16 = reinterpret_cast<const uint16_t*>(hom);
  for (int t = 0; t < nwl; ++t) {
    const int q = wl[t];
    const __mmask16 h = h16[q], o = o16[q];
    if ((h | o) == 0) continue;
    const __m512 c = _mm512_loadu_ps(col + 16 * q);
    acc = _mm512_mask_add_ps(acc, h, acc, c);
    acc = _mm512_mask_add_ps(acc, o, acc, _mm512_add_ps(c, c));
  }
  return _mm512_reduce_add_ps(acc);
}
#elif defined(__AVX2__)
struct Lut8 {  // byte -> 8 lane masks
  alignas(32) uint32_t v[256][8];
  Lut8() { for (int x = 0; x < 256; ++x) for (int b = 0; b < 8; ++b) v[x][b] = ((x >> b) & 1) ? 0xFFFFFFFFu : 0u; }
};
const Lut8& lut8() { static const Lut8 L; return L; }
inline float plane_dot(const float* col, const uint8_t* het, const uint8_t* hom, int nbyte) {
  const Lut8& L = lut8();
  __m256 acc = _mm256_setzero_ps();
  for (int q = 0; q < nbyte; ++q) {
    const uint8_t h = het[q], o = hom[q];
    if ((h | o) == 0) continue;
    const __m256 c = _mm256_loadu_ps(col + 8 * q);
    acc = _mm256_add_ps(acc, _mm256_and_ps(c, _mm256_load_ps(reinterpret_cast<const float*>(L.v[h]))));
    acc = _mm256_add_ps(acc, _mm256_and_ps(_mm256_add_ps(c, c), _mm256_load_ps(reinterpret_cast<const float*>(L.v[o]))));
  }
  __m128 s = _mm_add_ps(_mm256_castps256_ps128(acc), _mm256_extractf128_ps(acc, 1));
  s = _mm_hadd_ps(s, s); s = _mm_hadd_ps(s, s);
  return _mm_cvtss_f32(s);
}
#else
inline float plane_dot(const float* col, const uint8_t* het, const uint8_t* hom, int nbyte) {
  float acc = 0;
  for (int q = 0; q < nbyte; ++q) {
    const uint8_t h = het[q], o = hom[q];
    if ((h | o) == 0) continue;
    for (int b = 0; b < 8; ++b) acc += col[8 * q + b] * static_cast<float>(((h >> b) & 1) + 2 * ((o >> b) & 1));
  }
  return acc;
}
#endif
#if !defined(__AVX512F__)
inline float plane_dot_words(const float* col, const uint8_t* het, const uint8_t* hom, const int* wl, int nwl) {
  float acc = 0;
  for (int t = 0; t < nwl; ++t) {
    const int q = wl[t];
    acc += plane_dot(col + 16 * q, het + 2 * q, hom + 2 * q, 2);
  }
  return acc;
}
#endif
}  // namespace

// var_idx: 0-based PGEN variant per eigen row (-1 = not in PGEN: missing for everyone, not returned).
// sample_idx: 0-based PSAM rows to impute (anchor counts use these samples only).
// [[Rcpp::export]]
Rcpp::List genotype_impute_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct,
                               Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                               const std::string& eigen_file, Rcpp::NumericVector af, double thresh,
                               double delta, double anchor_min_miss, Rcpp::LogicalVector keep_row, int threads,
                               int nb_people, double mem_gb, double nbr_frac, int nbr_max, Rcpp::IntegerVector rows) {
  const int T = std::max(1, threads);
  const auto t_all = std::chrono::steady_clock::now();
  auto t0 = t_all;

  // ---- LD: whole block C = U W U'; window (rows = subset of eigen rows) C = I - delta (R_WW + delta I)^{-1}
  const plink4r::EigenBlock ld = plink4r::read_eigen_block(eigen_file, thresh);
  const bool window = rows.size() > 0;
  const int m = window ? static_cast<int>(rows.size()) : ld.m, k = ld.k;
  if (window && !(delta > 0)) Rcpp::stop("a window (rows) needs delta > 0");
  if (var_idx.size() != m || af.size() != m || keep_row.size() != m) Rcpp::stop("var_idx/af/keep_row length != eigen m");
  std::vector<int> vi(var_idx.begin(), var_idx.end());
  std::vector<char> kr(m);
  for (int j = 0; j < m; ++j) kr[j] = keep_row[j] ? 1 : 0;
  std::vector<float> mu(m, 0.0f), sd(m, 1.0f), isd(m, 0.0f);
  for (int j = 0; j < m; ++j) {
    if (vi[j] < 0) continue;
    const double p = af[j];
    if (!(p > 0 && p < 1)) Rcpp::stop("af must be in (0,1) for every SNP present in the PGEN (eigen row %d)", j + 1);
    mu[j] = static_cast<float>(2 * p);
    sd[j] = static_cast<float>(std::sqrt(2 * p * (1 - p)));
    isd[j] = 1.0f / sd[j];
  }
  const double t_eig = secs(t0); t0 = std::chrono::steady_clock::now();
  const int mp = (m + 15) / 16 * 16;  // padded rows for the 16-wide plane dot
  VectorXf w(k);
  for (int j = 0; j < k; ++j) w[j] = ld.lambda[j] / (ld.lambda[j] + static_cast<float>(delta));
  Eigen::setNbThreads(T);
  MatrixXf C(m, m);
  if (!window) {
    MatrixXf UW = ld.U * w.asDiagonal();
    C.noalias() = UW * ld.U.transpose();
  } else {
    MatrixXf Uw(m, k);
    for (int r = 0; r < m; ++r) {
      if (rows[r] < 0 || rows[r] >= ld.m) Rcpp::stop("rows out of range");
      Uw.row(r) = ld.U.row(rows[r]);
    }
    MatrixXf Rw = (Uw * ld.lambda.asDiagonal()) * Uw.transpose();
    Rw.diagonal().array() += static_cast<float>(delta);
    Eigen::LLT<MatrixXf> llt(Rw);
    if (llt.info() != Eigen::Success) Rcpp::stop("Cholesky of R_WW + delta I failed");
    C = -static_cast<float>(delta) * llt.solve(MatrixXf::Identity(m, m));
    C.diagonal().array() += 1.0f;
  }
  Eigen::setNbThreads(1);
  MatrixXf Aw = MatrixXf::Zero(mp, m);  // a_lj, padded rows are 0
  std::vector<float> Aj(m);
  // Optional neighbour restriction of the rhs (C_MM stays exact): per SNP j keep the rows l with the largest
  // C_lj^2 until they cover nbr_frac of sum_l C_lj^2 (the rhs error variance from dropping l is ~ C_lj^2),
  // capped at nbr_max; other a_lj are set to 0, and only the 16-SNP words holding a kept row are scanned.
  const bool restrict_nbr = (nbr_frac > 0 && nbr_frac < 1) || nbr_max > 0;
  const int nw = mp / 16;
  std::vector<std::vector<int>> words(restrict_nbr ? m : 0);
  std::vector<int> kept(m, m);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(T)
#endif
  for (int j = 0; j < m; ++j) {
    float s = 0;
    for (int l = 0; l < m; ++l) { const float v = C(l, j) * isd[l]; Aw(l, j) = v; }
    if (restrict_nbr) {
      std::vector<int> idx(m);
      for (int l = 0; l < m; ++l) idx[l] = l;
      std::sort(idx.begin(), idx.end(), [&](int x, int y) { return std::fabs(C(x, j)) > std::fabs(C(y, j)); });
      double tot = 0;
      for (int l = 0; l < m; ++l) tot += static_cast<double>(C(l, j)) * C(l, j);
      const double target = (nbr_frac > 0 && nbr_frac < 1) ? nbr_frac * tot : tot;
      int K = 0;
      double cum = 0;
      while (K < m && (nbr_max <= 0 || K < nbr_max) && cum < target) { cum += static_cast<double>(C(idx[K], j)) * C(idx[K], j); ++K; }
      kept[j] = K;
      std::vector<char> on(m, 0), won(nw, 0);
      for (int q = 0; q < K; ++q) { on[idx[q]] = 1; won[idx[q] >> 4] = 1; }
      for (int l = 0; l < m; ++l) if (!on[l]) Aw(l, j) = 0.0f;
      for (int q = 0; q < nw; ++q) if (won[q]) words[j].push_back(q);
    }
    for (int l = 0; l < m; ++l) s += Aw(l, j) * mu[l];
    Aj[j] = -s;
  }
  double kept_mean = 0, words_frac = 0;
  for (int j = 0; j < m; ++j) { kept_mean += kept[j]; words_frac += restrict_nbr ? static_cast<double>(words[j].size()) / nw : 1.0; }
  kept_mean /= m; words_frac /= m;
  const double t_C = secs(t0); t0 = std::chrono::steady_clock::now();

  // ---- genotypes: SNP-major genovecs of the selected samples (codes 0/1/2 ALT count, 3 missing)
  const int ns = sample_idx.size();
  const int sb = (ns + 3) / 4;
  std::vector<int> sidx(sample_idx.begin(), sample_idx.end());
  std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  std::vector<uint8_t> GS(static_cast<std::size_t>(m) * sb, 0xFF);  // absent rows stay all-missing
  std::vector<double> sum_obs(m, 0.0);
  std::vector<int> bad(T, -1);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(T)
#endif
  for (int j = 0; j < m; ++j) {
    if (vi[j] < 0) continue;
    const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[j]), tid());
    if (!gv) { bad[tid()] = j; continue; }
    uint8_t* dst = &GS[static_cast<std::size_t>(j) * sb];
    std::memset(dst, 0, sb);
    double so = 0;
    for (int s = 0; s < ns; ++s) {
      const int r = sidx[s];
      const unsigned c = (gv[r >> 2] >> (2 * (r & 3))) & 3u;
      dst[s >> 2] |= static_cast<uint8_t>(c << (2 * (s & 3)));
      if (c != 3u) so += c;
    }
    sum_obs[j] = so;
  }
  for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop("PGEN read failed at eigen row %d", bad[t] + 1);
  const double t_geno = secs(t0); t0 = std::chrono::steady_clock::now();

  // ---- bit-planes and missing lists for all people (chunks of 256 people = one cache line per SNP row)
  const int nbyte = mp / 8;
  double c0 = cpu_secs();
  std::vector<uint8_t> hetA(static_cast<std::size_t>(ns) * nbyte, 0), homA(static_cast<std::size_t>(ns) * nbyte, 0);
  std::vector<std::vector<int>> MiA(ns);
  {
    const int nch = (ns + 255) / 256;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
    for (int ch = 0; ch < nch; ++ch) {
      const int i0 = ch * 256, i1 = std::min(ns, i0 + 256);
      for (int lb = 0; lb < nbyte; ++lb) {  // 8 SNPs -> one plane byte per person, branch-free
        const int l0 = 8 * lb, nl = std::min(8, m - l0);
        if (nl <= 0) break;
        const uint8_t* g[8];
        for (int b = 0; b < nl; ++b) g[b] = &GS[static_cast<std::size_t>(l0 + b) * sb];
        for (int s = i0; s < i1; ++s) {
          const int byte = s >> 2, sh = 2 * (s & 3);
          unsigned h = 0, o = 0, x = 0;
          for (int b = 0; b < nl; ++b) {
            const unsigned c = (g[b][byte] >> sh) & 3u;
            h |= static_cast<unsigned>(c == 1u) << b;
            o |= static_cast<unsigned>(c == 2u) << b;
            x |= static_cast<unsigned>(c == 3u) << b;
          }
          hetA[static_cast<std::size_t>(s) * nbyte + lb] = static_cast<uint8_t>(h);
          homA[static_cast<std::size_t>(s) * nbyte + lb] = static_cast<uint8_t>(o);
          while (x) { MiA[s].push_back(l0 + __builtin_ctz(x)); x &= x - 1; }
        }
      }
    }
  }
  std::vector<uint8_t>().swap(GS);
  const double t_planes = secs(t0), u_planes = (cpu_secs() - c0) / (t_planes * T);
  t0 = std::chrono::steady_clock::now();

  // ---- person blocks: fixed size, or as large as the memory budget for S_i, rhs and the SNP lists allows
  std::vector<int> bstart{0};
  {
    const double budget = mem_gb * 1e9;
    double acc = 0;
    for (int s = 0; s < ns; ++s) {
      const double q = static_cast<double>(MiA[s].size());
      const double cost = 4.0 * q * q + 12.0 * q;
      const bool full = nb_people > 0 ? (s - bstart.back() >= nb_people) : (acc + cost > budget && s > bstart.back());
      if (full) { bstart.push_back(s); acc = 0; }
      acc += cost;
    }
    bstart.push_back(ns);
  }
  const int nblk = static_cast<int>(bstart.size()) - 1;

  // ---- stage 1 per block
  std::vector<std::vector<int>> oi(T), oj(T);
  std::vector<std::vector<float>> og(T);
  double t_lists = 0, t_rhs = 0, t_chol = 0, cpu_rhs = 0, cpu_chol = 0;
  long long ns_tot = 0;
  std::vector<int> moff, jcnt(m + 1), jlist_i, jlist_a;
  std::vector<long long> soff;
  std::vector<float> rhs, S;
  for (int bk = 0; bk < nblk; ++bk) {
    const int b0 = bstart[bk], nb = bstart[bk + 1] - b0;
    auto tb = std::chrono::steady_clock::now();
    const uint8_t* het = &hetA[static_cast<std::size_t>(b0) * nbyte];
    const uint8_t* hom = &homA[static_cast<std::size_t>(b0) * nbyte];
    const std::vector<int>* Mi = &MiA[b0];
    // offsets and per-SNP lists of (person, position in M)
    moff.assign(nb + 1, 0);
    soff.assign(nb + 1, 0);
    for (int i = 0; i < nb; ++i) {
      const int q = static_cast<int>(Mi[i].size());
      moff[i + 1] = moff[i] + q;
      soff[i + 1] = soff[i] + static_cast<long long>(q) * q;
    }
    // per-SNP lists in parallel: thread t counts its person range, then fills its own slots (order = person order)
    {
      std::vector<int> tcnt(static_cast<std::size_t>(T) * m, 0);
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
      {
        const int th = tid();
        const int i0 = static_cast<int>(static_cast<long long>(nb) * th / T), i1 = static_cast<int>(static_cast<long long>(nb) * (th + 1) / T);
        int* c = &tcnt[static_cast<std::size_t>(th) * m];
        for (int i = i0; i < i1; ++i) for (int l : Mi[i]) ++c[l];
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
        {
          jcnt[0] = 0;
          for (int j = 0; j < m; ++j) {
            int run = jcnt[j];
            for (int t = 0; t < T; ++t) { const int v = tcnt[static_cast<std::size_t>(t) * m + j]; tcnt[static_cast<std::size_t>(t) * m + j] = run; run += v; }
            jcnt[j + 1] = run;
          }
          jlist_i.resize(moff[nb]); jlist_a.resize(moff[nb]);
        }
        for (int i = i0; i < i1; ++i)
          for (int a = 0; a < static_cast<int>(Mi[i].size()); ++a) {
            const int p = c[Mi[i][a]]++;
            jlist_i[p] = i; jlist_a[p] = a;
          }
      }
    }
    rhs.resize(moff[nb]);                          // fully overwritten below
    S.resize(static_cast<std::size_t>(soff[nb]));  // every column of every S_i is written below
    t_lists += secs(tb); tb = std::chrono::steady_clock::now();
    c0 = cpu_secs();

    // rhs and C_MM, column by column
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(T)
#endif
    for (int j = 0; j < m; ++j) {
      const int e0 = jcnt[j], e1 = jcnt[j + 1];
      if (e0 == e1) continue;
      const float* aj = Aw.col(j).data();
      const float* cj = C.col(j).data();
      for (int e = e0; e < e1; ++e) {
        const int i = jlist_i[e], a = jlist_a[e];
        const std::vector<int>& M = Mi[i];
        const int q = static_cast<int>(M.size());
        const uint8_t* hi = &het[static_cast<std::size_t>(i) * nbyte];
        const uint8_t* oi_ = &hom[static_cast<std::size_t>(i) * nbyte];
        float r = Aj[j] + (restrict_nbr ? plane_dot_words(aj, hi, oi_, words[j].data(), static_cast<int>(words[j].size()))
                                        : plane_dot(aj, hi, oi_, nbyte));
        float* Sc = &S[static_cast<std::size_t>(soff[i]) + static_cast<std::size_t>(a) * q];  // column a of S_i
        for (int bq = 0; bq < q; ++bq) {
          const int l = M[bq];
          r += aj[l] * mu[l];
          Sc[bq] = -cj[l];
        }
        Sc[a] += 1.0f;
        rhs[moff[i] + a] = r;
      }
    }
    t_rhs += secs(tb); tb = std::chrono::steady_clock::now();
    cpu_rhs += cpu_secs() - c0; c0 = cpu_secs();

    // per-person solve
#ifdef _OPENMP
#pragma omp parallel num_threads(T)
#endif
    {
      const int th = tid();
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
      for (int i = 0; i < nb; ++i) {
        const std::vector<int>& M = Mi[i];
        const int q = static_cast<int>(M.size());
        if (q == 0) continue;
        Eigen::Map<MatrixXf> Si(&S[static_cast<std::size_t>(soff[i])], q, q);
        Eigen::Map<VectorXf> ri(&rhs[moff[i]], q);
        VectorXf zm;
        if (delta > 0) {
          Eigen::LLT<Eigen::Ref<MatrixXf>> llt(Si);
          if (llt.info() == Eigen::Success) zm = llt.solve(ri);
          else {
            MatrixXf Sf(q, q);  // LLT overwrote the lower triangle and diagonal; rebuild from the strict upper one
            for (int x = 0; x < q; ++x)
              for (int y = 0; y < q; ++y) Sf(y, x) = y < x ? Si(y, x) : (y > x ? Si(x, y) : 1.0f - C(M[x], M[x]));
            zm = Sf.ldlt().solve(ri);
          }
        } else {
          Eigen::SelfAdjointEigenSolver<MatrixXf> es(Si);
          const VectorXf ev = es.eigenvalues();
          const float cut = 1e-6f * ev.maxCoeff();
          VectorXf y = es.eigenvectors().transpose() * ri;
          for (int x = 0; x < q; ++x) y[x] = ev[x] > cut ? y[x] / ev[x] : 0.0f;
          zm = es.eigenvectors() * y;
        }
        for (int x = 0; x < q; ++x) {
          const int j = M[x];
          if (vi[j] < 0) continue;
          oi[th].push_back(b0 + i); oj[th].push_back(j);
          og[th].push_back(std::min(2.0f, std::max(0.0f, mu[j] + sd[j] * zm[x])));
        }
      }
    }
    t_chol += secs(tb);
    cpu_chol += cpu_secs() - c0;
    for (int i = 0; i < nb; ++i) if (!Mi[i].empty()) ++ns_tot;
  }
  const double t_solve = secs(t0); t0 = std::chrono::steady_clock::now();

  // ---- stage 2: anchor per SNP
  std::vector<double> sum1(m, 0.0), cnt1(m, 0.0);
  for (int t = 0; t < T; ++t)
    for (std::size_t r = 0; r < og[t].size(); ++r) { sum1[oj[t][r]] += og[t][r]; cnt1[oj[t][r]] += 1; }
  std::vector<double> craw(m, NA_REAL), cc(m, NA_REAL), shift(m, 0.0), mean1(m, NA_REAL);
  for (int j = 0; j < m; ++j) {
    if (vi[j] < 0 || cnt1[j] == 0) continue;
    mean1[j] = sum1[j] / cnt1[j];
    craw[j] = (2.0 * ns * af[j] - sum_obs[j]) / cnt1[j];
    cc[j] = std::min(2.0, std::max(0.0, craw[j]));
    if (anchor_min_miss >= 0 && cnt1[j] >= anchor_min_miss * ns) shift[j] = cc[j] - mean1[j];
  }
  std::size_t nout = 0;
  for (int t = 0; t < T; ++t) for (std::size_t r = 0; r < oj[t].size(); ++r) if (kr[oj[t][r]]) ++nout;
  Rcpp::IntegerVector ri(nout), rj(nout);
  Rcpp::NumericVector rg1(nout), rg(nout);
  std::size_t p = 0;
  for (int t = 0; t < T; ++t)
    for (std::size_t r = 0; r < oj[t].size(); ++r) {
      const int j = oj[t][r];
      if (!kr[j]) continue;
      ri[p] = oi[t][r] + 1; rj[p] = j + 1; rg1[p] = og[t][r];
      rg[p] = std::min(2.0, std::max(0.0, og[t][r] + shift[j]));
      ++p;
    }
  const double t_anchor = secs(t0);
  Rcpp::IntegerVector jj(m), mj(m);
  for (int j = 0; j < m; ++j) { jj[j] = j + 1; mj[j] = static_cast<int>(cnt1[j]); }
  return Rcpp::List::create(
    Rcpp::_["fill"] = Rcpp::DataFrame::create(Rcpp::_["i"] = ri, Rcpp::_["j"] = rj, Rcpp::_["g1"] = rg1, Rcpp::_["g"] = rg),
    Rcpp::_["snp"] = Rcpp::DataFrame::create(Rcpp::_["j"] = jj, Rcpp::_["n_miss"] = mj, Rcpp::_["c_raw"] = craw,
                                             Rcpp::_["c"] = cc, Rcpp::_["mean_stage1"] = mean1, Rcpp::_["shift"] = shift),
    Rcpp::_["info"] = Rcpp::List::create(Rcpp::_["m"] = m, Rcpp::_["k"] = k, Rcpp::_["n"] = ns, Rcpp::_["delta"] = delta,
                                         Rcpp::_["persons_solved"] = static_cast<double>(ns_tot), Rcpp::_["threads"] = T,
                                         Rcpp::_["n_blocks"] = nblk, Rcpp::_["people_per_block"] = static_cast<double>(ns) / nblk,
                                         Rcpp::_["cpu_util_planes"] = u_planes,
                                         Rcpp::_["nbr_kept_mean"] = kept_mean, Rcpp::_["nbr_words_frac"] = words_frac,
                                         Rcpp::_["cpu_util_rhs"] = t_rhs > 0 ? cpu_rhs / (t_rhs * T) : NA_REAL,
                                         Rcpp::_["cpu_util_cholesky"] = t_chol > 0 ? cpu_chol / (t_chol * T) : NA_REAL),
    Rcpp::_["time"] = Rcpp::NumericVector::create(Rcpp::_["read_eigen"] = t_eig, Rcpp::_["build_C"] = t_C,
                                                  Rcpp::_["read_geno"] = t_geno, Rcpp::_["planes"] = t_planes, Rcpp::_["lists"] = t_lists,
                                                  Rcpp::_["rhs_CMM"] = t_rhs, Rcpp::_["cholesky"] = t_chol,
                                                  Rcpp::_["solve"] = t_solve, Rcpp::_["anchor"] = t_anchor,
                                                  Rcpp::_["total"] = secs(t_all)));
}
