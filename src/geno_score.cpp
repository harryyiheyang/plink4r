// Polygenic scoring kernels behind geno_score().
//
// score_stream_cpp: plink2 --score style streaming. Every variant carries a 4-entry lookup table per trait
//   (value added to the score for PGEN code 0/1/2 = ALT count, 3 = missing) and a 4-entry allele-count table.
//   R builds the tables (A1 alignment, weights, transforms, the fill of a missing call: plink2 mean, scheme 3
//   anchor, none, or the constant part of an LD-filled cell), so the kernel holds no policy. Variants are read
//   in chunks (parallel over variants), then scored in parallel over sample ranges, so no reduction and no
//   genotype matrix: memory is n x q scores plus one chunk of genovecs.
// score_ld_block_cpp: one eigen.bin LD block through impute_kernel.h (the same stage 1 / stage 2 as
//   genotype_impute), adding B[j, t] * g for every filled cell instead of returning the cells.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include "impute_kernel.h"

namespace {
// Byte decoding: a PGEN byte holds 4 samples (2 bits each, low bits first).
struct ByteTables {
  float x[256][4];        // ALT count, missing -> 0
  unsigned char miss[256];  // 4-bit missing mask
  ByteTables() {
    for (int b = 0; b < 256; ++b) {
      miss[b] = 0;
      for (int k = 0; k < 4; ++k) {
        const int c = (b >> (2 * k)) & 3;
        x[b][k] = c == 3 ? 0.0f : static_cast<float>(c);
        if (c == 3) miss[b] |= static_cast<unsigned char>(1u << k);
      }
    }
  }
};
const ByteTables& byte_tables() { static const ByteTables T; return T; }
}  // namespace

// One pass per chunk of variants: read and compact the genovecs, count each variant's observed calls and ALT
// sum over the selected samples, build its 4-entry tables here (A1 alignment, transform, fill of a missing
// call), then score. So the PGEN is read once, not once for the counts and again for the scores.
// var_idx: 0-based PGEN variants (file order is fastest); sample_idx: 0-based PSAM rows (output row order).
// W: raw weights (variants x q); flip: A1 is REF; pfreq: given A1 frequency (NaN = observed);
// ld_cand: the variant goes through the LD panel when it has missing calls (its missing cells then add only the
// constant part here, the per-cell part comes from score_ld_block_cpp); impute_none: missing calls add nothing;
// transform: 0 none, 1 center, 2 standardize. q = 1: per-variant byte lookup; q >= 2: float32 product.
// Per variant out: n_obs, sum (ALT), status (0 scored, 1 all calls missing, 2 zero variance under
// standardize), fill (0 none, 1 mean, 2 anchor, 3 ld), p1 (A1 frequency used), a (scale).
// [[Rcpp::export]]
Rcpp::List score_stream_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct,
                            Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx, Rcpp::NumericMatrix W,
                            Rcpp::IntegerVector flip, Rcpp::NumericVector pfreq, Rcpp::LogicalVector ld_cand,
                            bool impute_none, int transform, int threads, int chunk) {
  const int T = plink4r::resolve_threads(threads);
  const int p = var_idx.size(), ns = sample_idx.size(), q = W.ncol();
  if (W.nrow() != p || flip.size() != p || pfreq.size() != p || ld_cand.size() != p)
    Rcpp::stop("weights / flip / pfreq / ld_cand do not match the variants");
  const bool gemm = q >= 2;
  const std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  const std::vector<int> act(allele_ct.begin(), allele_ct.end());
  const bool ident = plink4r::identity_samples(si.data(), ns, n_raw);
  const ByteTables& BT = byte_tables();
  const plink4r::LaFns& LA = plink4r::la();
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  if (pf.dosage()) Rcpp::stop("geno_score scores hard calls; this PGEN holds dosages.");
  const std::size_t sb = (static_cast<std::size_t>(ns) + 3) / 4;  // bytes per compacted variant
  const int V = std::max(1, chunk);
  std::vector<unsigned char> buf(static_cast<std::size_t>(V) * sb);
  std::vector<double> S(static_cast<std::size_t>(ns) * q, 0.0);  // row-major: sample, trait
  std::vector<int> ct(ns, 0);
  std::vector<int> bad(T, -1);
  const int sblk = gemm ? 1024 : 4096;  // multiple of 4: a sample block starts on a byte
  const int nsb = (ns + sblk - 1) / sblk;
  std::vector<double> base_sum(q, 0.0);
  std::vector<double> nobs(p), sumv(p), p1v(p, NA_REAL), av(p, 1.0);
  std::vector<int> status(p, 0), fillv(p, 0);
  std::vector<double> tab(static_cast<std::size_t>(V) * 4 * q);  // chunk tables [variant][code][trait]
  std::vector<int> ctab(static_cast<std::size_t>(V) * 4);
  std::vector<double> L1;  // lookup path: per chunk variant, 256 x 4 values (q = 1)
  if (!gemm) L1.resize(static_cast<std::size_t>(V) * 256 * 4);
  for (int c0 = 0; c0 < p; c0 += V) {
    const int nv = std::min(V, p - c0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(T)
#endif
    for (int v = 0; v < nv; ++v) {
      const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[c0 + v]), plink4r::ik::tid());
      if (!gv) { bad[plink4r::ik::tid()] = c0 + v; continue; }
      unsigned char* dst = &buf[static_cast<std::size_t>(v) * sb];
      plink4r::compact_genovec(gv, si.data(), ns, ident, dst);
      long long miss = 0, alt = 0;
      plink4r::count_genovec(dst, sb, &miss, &alt);
      nobs[c0 + v] = static_cast<double>(ns - miss); sumv[c0 + v] = static_cast<double>(alt);
    }
    for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop("PGEN read failed at variant %d", bad[t] + 1);
    // tables of this chunk: value a (d - off) w for A1 dosage d = c0 + s x, x = ALT count of code 0/1/2
    std::fill(tab.begin(), tab.end(), 0.0);
    std::fill(ctab.begin(), ctab.end(), 0);
    for (int v = 0; v < nv; ++v) {
      const int j = c0 + v;
      const double no = nobs[j], so = sumv[j], nm = ns - no;
      if (no == 0) { status[j] = 1; continue; }
      const bool fl = flip[j] != 0;
      const double pg = pfreq[j];
      const bool given = std::isfinite(pg);
      const double palt_obs = so / (2 * no);
      const double p1 = given ? pg : (fl ? 1 - palt_obs : palt_obs);
      p1v[j] = p1;
      if (transform == 2 && !(p1 > 0 && p1 < 1)) { status[j] = 2; continue; }
      const double off = transform ? 2 * p1 : 0, a = transform == 2 ? 1 / std::sqrt(2 * p1 * (1 - p1)) : 1.0;
      av[j] = a;
      const double sg = fl ? -1.0 : 1.0, cc = fl ? 2.0 : 0.0;
      double f_alt = so / no;
      int fill = impute_none ? 0 : 1;
      if (!impute_none && given) {
        fill = 2;
        if (nm > 0) {
          const double pa = fl ? 1 - pg : pg;
          f_alt = std::min(2.0, std::max(0.0, (2.0 * ns * pa - so) / nm));
        }
      }
      if (!impute_none && ld_cand[j] && nm > 0) fill = 3;
      fillv[j] = fill;
      const double dmiss = fill == 0 ? 0.0 : (fill == 3 ? cc - off : cc + sg * f_alt - off);
      double* tv = &tab[static_cast<std::size_t>(v) * 4 * q];
      for (int t = 0; t < q; ++t) {
        const double wa = W(j, t) * a;
        for (int code = 0; code < 3; ++code) tv[code * q + t] = wa * (cc + sg * code - off);
        tv[3 * q + t] = fill == 0 ? 0.0 : wa * dmiss;
      }
      int* cv = &ctab[static_cast<std::size_t>(v) * 4];
      cv[0] = cv[1] = cv[2] = 2; cv[3] = fill == 0 ? 0 : 2;
    }
    int ctsum = 0;
    for (int v = 0; v < nv; ++v) ctsum += ctab[static_cast<std::size_t>(v) * 4];
    if (gemm) {
      Eigen::MatrixXf SlopeF(nv, q);
      Eigen::MatrixXd MissD(nv, q);
      std::vector<int> ct_mis(nv);
      for (int v = 0; v < nv; ++v) {
        const double* tv = &tab[static_cast<std::size_t>(v) * 4 * q];
        for (int t = 0; t < q; ++t) {
          SlopeF(v, t) = static_cast<float>(tv[q + t] - tv[t]);
          MissD(v, t) = tv[3 * q + t] - tv[t];
          base_sum[t] += tv[t];
        }
        ct_mis[v] = ctab[static_cast<std::size_t>(v) * 4 + 3] - ctab[static_cast<std::size_t>(v) * 4];
      }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
      for (int b = 0; b < nsb; ++b) {
        const int s0 = b * sblk, s1 = std::min(ns, s0 + sblk), nr = s1 - s0;
        const int nr4 = (nr + 3) / 4 * 4;
        Eigen::MatrixXf XF(nr4, nv);  // padded to whole bytes; padding rows are dropped below
        std::vector<int> mi, mv;
        for (int v = 0; v < nv; ++v) {
          const unsigned char* g = &buf[static_cast<std::size_t>(v) * sb + (s0 >> 2)];
          float* col = XF.col(v).data();
          for (int y = 0; y < nr4 / 4; ++y) {
            const unsigned byte = g[y];
            std::memcpy(col + 4 * y, BT.x[byte], 4 * sizeof(float));
            unsigned mm = BT.miss[byte];
            while (mm) {
              const int k = __builtin_ctz(mm); mm &= mm - 1;
              const int s = s0 + 4 * y + k;
              if (s < s1) { mi.push_back(s); mv.push_back(v); }
            }
          }
        }
        Eigen::MatrixXf PF(nr, q);
        LA.sgemm_nn(XF.data(), nr, nv, nr4, SlopeF.data(), q, nv, PF.data(), nr);
        for (int s = s0; s < s1; ++s) {
          double* Ss = &S[static_cast<std::size_t>(s) * q];
          for (int t = 0; t < q; ++t) Ss[t] += PF(s - s0, t);
          ct[s] += ctsum;
        }
        for (std::size_t e = 0; e < mi.size(); ++e) {
          double* Ss = &S[static_cast<std::size_t>(mi[e]) * q];
          for (int t = 0; t < q; ++t) Ss[t] += MissD(mv[e], t);
          ct[mi[e]] += ct_mis[mv[e]];
        }
      }
      continue;
    }
    // lookup path (q = 1): per variant a 256 x 4 table byte -> the values of its 4 samples
    for (int v = 0; v < nv; ++v) {
      const double* tv = &tab[static_cast<std::size_t>(v) * 4];
      double* Lv = &L1[static_cast<std::size_t>(v) * 1024];
      for (int byte = 0; byte < 256; ++byte)
        for (int k = 0; k < 4; ++k) Lv[4 * byte + k] = tv[(byte >> (2 * k)) & 3];
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
    for (int b = 0; b < nsb; ++b) {
      const int s0 = b * sblk, s1 = std::min(ns, s0 + sblk);
      const int full = s0 + (s1 - s0) / 4 * 4;  // samples covered by whole bytes
      for (int s = s0; s < s1; ++s) ct[s] += ctsum;  // every call counts; missing ones corrected below
      for (int v = 0; v < nv; ++v) {
        const unsigned char* g = &buf[static_cast<std::size_t>(v) * sb];
        const double* Lv = &L1[static_cast<std::size_t>(v) * 1024];
        const int* cv = &ctab[static_cast<std::size_t>(v) * 4];
        const int dct = cv[3] - cv[0];
        double* Sp = S.data();
        for (int s = s0; s < full; s += 4) {
          const unsigned byte = g[s >> 2];
          const double* e = Lv + 4 * byte;
          Sp[s] += e[0]; Sp[s + 1] += e[1]; Sp[s + 2] += e[2]; Sp[s + 3] += e[3];
          if (dct) for (unsigned mm = BT.miss[byte]; mm; mm &= mm - 1) ct[s + __builtin_ctz(mm)] += dct;
        }
        for (int s = full; s < s1; ++s) {
          const unsigned c = (g[s >> 2] >> (2 * (s & 3))) & 3u;
          Sp[s] += Lv[4 * g[s >> 2] + (s & 3)];
          if (c == 3u) ct[s] += dct;
        }
      }
    }
  }
  if (gemm) for (int s = 0; s < ns; ++s) for (int t = 0; t < q; ++t) S[static_cast<std::size_t>(s) * q + t] += base_sum[t];
  Rcpp::NumericMatrix out(ns, q);
  for (int s = 0; s < ns; ++s) for (int t = 0; t < q; ++t) out(s, t) = S[static_cast<std::size_t>(s) * q + t];
  return Rcpp::List::create(Rcpp::_["score"] = out, Rcpp::_["allele_ct"] = Rcpp::wrap(ct),
                            Rcpp::_["n_obs"] = Rcpp::wrap(nobs), Rcpp::_["sum"] = Rcpp::wrap(sumv),
                            Rcpp::_["status"] = Rcpp::wrap(status), Rcpp::_["fill"] = Rcpp::wrap(fillv),
                            Rcpp::_["p1"] = Rcpp::wrap(p1v), Rcpp::_["a"] = Rcpp::wrap(av));
}

// One LD block over its working rows (panel rows in the PGEN; rows / var_idx / af / sign as in
// genotype_impute_cpp): returns the n x q matrix sum_j B[j, t] * g_ij over the filled cells (g = final ALT
// dosage of a missing call). Rows of B that are all zero are imputed (they condition the others) but add nothing.
// [[Rcpp::export]]
Rcpp::List score_ld_block_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct,
                              Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                              const std::string& eigen_file, Rcpp::NumericVector af, Rcpp::IntegerVector sign, double thresh, double delta,
                              double anchor_min_miss, Rcpp::NumericMatrix B, int threads, int nb_people,
                              double mem_gb, Rcpp::IntegerVector rows) {
  const int T = plink4r::resolve_threads(threads);
  const plink4r::EigenBlock eb = plink4r::read_eigen_block(eigen_file, thresh);
  const std::vector<int> rw(rows.begin(), rows.end());
  const int m = static_cast<int>(rw.size());
  if (var_idx.size() != m || af.size() != m || sign.size() != m || B.nrow() != m)
    Rcpp::stop("rows / var_idx / af / sign / B rows differ in length");
  for (int r : rw) if (r < 0 || r >= eb.m) Rcpp::stop("rows out of range of %s", eigen_file);
  const std::vector<int> sg(sign.begin(), sign.end());
  const int q = B.ncol();
  const std::vector<int> vi(var_idx.begin(), var_idx.end());
  const std::vector<double> afv(af.begin(), af.end());
  std::vector<char> keep(m, 0);
  for (int j = 0; j < m; ++j) for (int t = 0; t < q; ++t) if (B(j, t) != 0) { keep[j] = 1; break; }
  const plink4r::ImputeLD L = plink4r::build_impute_ld(eb, rw, vi, afv, sg, delta, T);
  const std::vector<int> sidx(sample_idx.begin(), sample_idx.end());
  const std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  plink4r::BlockGenotypes G = plink4r::read_block_genotypes(pf, vi, sidx, T);
  const int ns = G.ns;
  const plink4r::PersonPlanes P = plink4r::build_planes(L, G, T);
  const plink4r::Stage1Result R = plink4r::impute_stage1(L, P, delta, nb_people, mem_gb, T);
  const plink4r::Stage2Result A = plink4r::impute_stage2(R, G.sum_obs, afv, ns, anchor_min_miss);
  Rcpp::NumericMatrix S(ns, q);
  double cells = 0;
  plink4r::for_each_cell(R, A, keep, [&](int i, int j, float, double g) {
    for (int t = 0; t < q; ++t) S(i, t) += B(j, t) * g;
    cells += 1;
  });
  return Rcpp::List::create(Rcpp::_["score"] = S, Rcpp::_["cells"] = cells,
                            Rcpp::_["persons_solved"] = static_cast<double>(R.persons_solved));
}

// Instruction set of the SIMD kernels picked at run time (avx512, avx2 or generic; PLINK4R_SIMD caps it).
// [[Rcpp::export]]
std::string simd_isa_cpp() { return plink4r::la().isa; }
