// Genotype imputation of hard-call missing cells from an eigen.bin LD block: R entry point.
// The algorithm (stage 1 conditional mean, stage 2 anchor) lives in impute_kernel.h, shared with geno_score.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>
#include "impute_kernel.h"
#include "geno_access.h"

namespace {

// Working rows of one block / window, checked against the eigen block.
struct ImputeJob {
  plink4r::EigenBlock eb;
  std::vector<int> rw, sg, vi;
  std::vector<double> afv;
  std::vector<char> kr;
  int m = 0;
  double t_read = 0;
};

ImputeJob impute_job(const std::string& eigen_file, double thresh, const Rcpp::IntegerVector& rows,
                     const Rcpp::IntegerVector& var_idx, const Rcpp::NumericVector& af,
                     const Rcpp::IntegerVector& sign, const Rcpp::LogicalVector& keep_row,
                     std::chrono::steady_clock::time_point t_all) {
  ImputeJob J;
  J.eb = plink4r::read_eigen_block(eigen_file, thresh);
  J.t_read = plink4r::ik::secs(t_all);
  J.rw.assign(rows.begin(), rows.end());
  J.m = static_cast<int>(J.rw.size());
  const int m = J.m;
  if (var_idx.size() != m || af.size() != m || sign.size() != m || keep_row.size() != m)
    Rcpp::stop("rows / var_idx / af / sign / keep_row differ in length");
  for (int r : J.rw) if (r < 0 || r >= J.eb.m) Rcpp::stop("rows out of range of %s", eigen_file);
  J.sg.assign(sign.begin(), sign.end());
  J.vi.assign(var_idx.begin(), var_idx.end());
  J.afv.assign(af.begin(), af.end());
  J.kr.assign(m, 0);
  for (int j = 0; j < m; ++j) J.kr[j] = keep_row[j] ? 1 : 0;
  return J;
}

// Stage 1 + 2 on the genotypes of the working rows (shared by the PGEN and the geno-object entries).
struct ImputeOut {
  plink4r::Stage1Result R;
  plink4r::Stage2Result A;
  double t_anchor = 0, t_planes = 0;
  int ns = 0;
};

ImputeOut impute_compute(const ImputeJob& J, const plink4r::ImputeLD& L, plink4r::BlockGenotypes& B, double delta,
                         double anchor_min_miss, int nb_people, double mem_gb, int T) {
  const std::vector<char>& kr = J.kr;
  const int ns = B.ns;
  plink4r::PersonPlanes P = plink4r::build_planes(L, B, T);
  // People whose missing cells all sit in rows that are not returned (flank rows) need no solve: drop their
  // lists. Every person with a missing kept cell keeps the full M_i, so the kept fills are unchanged.
  if (std::find(kr.begin(), kr.end(), 0) != kr.end()) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) num_threads(T)
#endif
    for (int i = 0; i < ns; ++i) {
      std::vector<int>& M = P.Mi[i];
      bool any = false;
      for (int l : M) if (kr[l]) { any = true; break; }
      if (!any) std::vector<int>().swap(M);
    }
  }
  ImputeOut O;
  O.ns = ns;
  O.t_planes = P.t_planes;
  O.R = plink4r::impute_stage1(L, P, delta, nb_people, mem_gb, T);
  const auto t2 = std::chrono::steady_clock::now();
  O.A = plink4r::impute_stage2(O.R, B.sum_obs, J.afv, ns, anchor_min_miss);
  O.t_anchor = plink4r::ik::secs(t2);
  return O;
}

Rcpp::NumericVector na_vec(const std::vector<double>& x) {
  Rcpp::NumericVector out(x.size());
  for (std::size_t j = 0; j < x.size(); ++j) out[j] = std::isnan(x[j]) ? NA_REAL : x[j];
  return out;
}

Rcpp::NumericVector impute_times(const plink4r::ImputeLD& L, const plink4r::BlockGenotypes& B, const ImputeOut& O,
                                 double t_write, std::chrono::steady_clock::time_point t_all) {
  const plink4r::Stage1Result& R = O.R;
  Rcpp::NumericVector t = Rcpp::NumericVector::create(
    Rcpp::_["read_eigen"] = L.t_eig, Rcpp::_["build_C"] = L.t_C, Rcpp::_["read_geno"] = B.t_geno,
    Rcpp::_["planes"] = O.t_planes, Rcpp::_["lists"] = R.t_lists, Rcpp::_["rhs_CMM"] = R.t_rhs,
    Rcpp::_["cholesky"] = R.t_chol, Rcpp::_["solve"] = R.t_solve, Rcpp::_["anchor"] = O.t_anchor);
  if (t_write >= 0) t.push_back(t_write, "write");
  t.push_back(plink4r::ik::secs(t_all), "total");
  return t;
}

// Per-SNP table, info and timings (everything but the cells); t_write < 0 = no write step.
Rcpp::List impute_summary(const ImputeJob& J, const plink4r::ImputeLD& L, const plink4r::BlockGenotypes& B,
                          const ImputeOut& O, double delta, int T, double t_write,
                          std::chrono::steady_clock::time_point t_all) {
  const int m = J.m;
  const plink4r::Stage1Result& R = O.R;
  const plink4r::Stage2Result& A = O.A;
  Rcpp::IntegerVector jj(m), mj(m);
  for (int j = 0; j < m; ++j) { jj[j] = j + 1; mj[j] = static_cast<int>(A.cnt1[j]); }
  const int nblk = R.n_blocks;
  return Rcpp::List::create(
    Rcpp::_["snp"] = Rcpp::DataFrame::create(Rcpp::_["j"] = jj, Rcpp::_["n_miss"] = mj, Rcpp::_["c_raw"] = na_vec(A.craw),
                                             Rcpp::_["c"] = na_vec(A.cc), Rcpp::_["mean_stage1"] = na_vec(A.mean1),
                                             Rcpp::_["shift"] = Rcpp::wrap(A.shift)),
    Rcpp::_["info"] = Rcpp::List::create(Rcpp::_["m"] = m, Rcpp::_["k"] = L.k, Rcpp::_["n"] = O.ns, Rcpp::_["delta"] = delta,
                                         Rcpp::_["persons_solved"] = static_cast<double>(R.persons_solved), Rcpp::_["threads"] = T,
                                         Rcpp::_["n_blocks"] = nblk, Rcpp::_["people_per_block"] = static_cast<double>(O.ns) / nblk),
    Rcpp::_["time"] = impute_times(L, B, O, t_write, t_all));
}

// Cells of the kept rows as an R data frame (i, j = working row, g1, g), plus the summary.
Rcpp::List impute_run(const ImputeJob& J, const plink4r::ImputeLD& L, plink4r::BlockGenotypes& B, double delta,
                      double anchor_min_miss, int nb_people, double mem_gb, int T,
                      std::chrono::steady_clock::time_point t_all) {
  const ImputeOut O = impute_compute(J, L, B, delta, anchor_min_miss, nb_people, mem_gb, T);
  const std::vector<char>& kr = J.kr;
  std::size_t nout = 0;
  for (std::size_t t = 0; t < O.R.oj.size(); ++t) for (std::size_t r = 0; r < O.R.oj[t].size(); ++r) if (kr[O.R.oj[t][r]]) ++nout;
  Rcpp::IntegerVector ri(nout), rj(nout);
  Rcpp::NumericVector rg1(nout), rg(nout);
  std::size_t p = 0;
  plink4r::for_each_cell(O.R, O.A, kr, [&](int i, int j, float g1, double g) {
    ri[p] = i + 1; rj[p] = j + 1; rg1[p] = g1; rg[p] = g; ++p;
  });
  Rcpp::List out = impute_summary(J, L, B, O, delta, T, -1.0, t_all);
  out["fill"] = Rcpp::DataFrame::create(Rcpp::_["i"] = ri, Rcpp::_["j"] = rj, Rcpp::_["g1"] = rg1, Rcpp::_["g"] = rg);
  return out;
}

}  // namespace

// Working rows (the panel rows of the block / window that are in the PGEN, increasing): rows = 0-based eigen
// row, var_idx = 0-based PGEN variant, af = PGEN-ALT frequency, sign = +1 (PGEN ALT = panel alt) or -1.
// sample_idx: 0-based PSAM rows to impute (anchor counts use these samples only).
// [[Rcpp::export]]
Rcpp::List genotype_impute_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct,
                               Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                               const std::string& eigen_file, Rcpp::NumericVector af, Rcpp::IntegerVector sign, double thresh,
                               double delta, double anchor_min_miss, Rcpp::LogicalVector keep_row, int threads,
                               int nb_people, double mem_gb, Rcpp::IntegerVector rows) {
  const int T = plink4r::resolve_threads(threads);
  const auto t_all = std::chrono::steady_clock::now();
  const ImputeJob J = impute_job(eigen_file, thresh, rows, var_idx, af, sign, keep_row, t_all);
  const plink4r::ImputeLD L = plink4r::build_impute_ld(J.eb, J.rw, J.vi, J.afv, J.sg, delta, T, J.t_read);
  const std::vector<int> sidx(sample_idx.begin(), sample_idx.end());
  const std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  plink4r::BlockGenotypes B = plink4r::read_block_genotypes(pf, J.vi, sidx, T);
  return impute_run(J, L, B, delta, anchor_min_miss, nb_people, mem_gb, T, t_all);
}

// Same on a geno object (hard calls): var_idx = 0-based columns of the object, every row of the object is
// imputed; the genotypes come from its planes and its missing-cell rows, in its current counted-allele
// orientation (af and sign refer to that allele).
// [[Rcpp::export]]
Rcpp::List genotype_impute_geno_cpp(SEXP ptr, Rcpp::IntegerVector var_idx, const std::string& eigen_file,
                                    Rcpp::NumericVector af, Rcpp::IntegerVector sign, double thresh, double delta,
                                    double anchor_min_miss, Rcpp::LogicalVector keep_row, int threads, int nb_people,
                                    double mem_gb, Rcpp::IntegerVector rows) {
  const int T = plink4r::resolve_threads(threads);
  const auto t_all = std::chrono::steady_clock::now();
  const ImputeJob J = impute_job(eigen_file, thresh, rows, var_idx, af, sign, keep_row, t_all);
  const plink4r::ImputeLD L = plink4r::build_impute_ld(J.eb, J.rw, J.vi, J.afv, J.sg, delta, T, J.t_read);
  const auto t0 = std::chrono::steady_clock::now();
  plink4r::BlockGenotypes B;
  int ns = 0;
  plink4r::geno_genovec(ptr, J.vi, B.GS, B.sum_obs, ns, T);
  B.m = J.m;
  B.ns = ns;
  B.sb = (ns + 3) / 4;
  B.t_geno = plink4r::ik::secs(t0);
  return impute_run(J, L, B, delta, anchor_min_miss, nb_people, mem_gb, T, t_all);
}

// Observed calls and ALT sum (hard calls, or dosages when the PGEN holds them) of each variant over the selected
// samples, for genotype_fill. var_idx / sample_idx 0-based.
// [[Rcpp::export]]
Rcpp::List pgen_stats_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct, Rcpp::IntegerVector var_idx,
                          Rcpp::IntegerVector sample_idx, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const int p = var_idx.size(), ns = sample_idx.size();
  std::vector<int> vi(var_idx.begin(), var_idx.end()), si(sample_idx.begin(), sample_idx.end());
  std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  const bool dos = pf.dosage();
  const bool ident = plink4r::identity_samples(si.data(), ns, n_raw);
  const std::size_t sb = (static_cast<std::size_t>(ns) + 3) / 4;
  std::vector<unsigned char> cbuf(dos ? 0 : static_cast<std::size_t>(T) * sb);
  std::vector<double> nobs(p), sum(p);
  std::vector<int> bad(T, -1);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(T)
#endif
  for (int j = 0; j < p; ++j) {
    long long c = 0; double s = 0;
    if (dos) {
      const std::uint16_t* d = pf.read_dosage(static_cast<std::size_t>(vi[j]), plink4r::thread_id());
      if (!d) { bad[plink4r::thread_id()] = j; continue; }
      for (int i = 0; i < ns; ++i) { const std::uint16_t v = d[si[i]]; if (v != plink4r::PgenFile::kMissing) { ++c; s += v; } }
      s /= 16384.0;
    } else {
      const unsigned char* gv = pf.read(static_cast<std::size_t>(vi[j]), plink4r::thread_id());
      if (!gv) { bad[plink4r::thread_id()] = j; continue; }
      unsigned char* dst = &cbuf[static_cast<std::size_t>(plink4r::thread_id()) * sb];
      plink4r::compact_genovec(gv, si.data(), ns, ident, dst);
      long long miss = 0, alt = 0;
      plink4r::count_genovec(dst, sb, &miss, &alt);
      c = ns - miss; s = static_cast<double>(alt);
    }
    nobs[j] = static_cast<double>(c); sum[j] = s;
  }
  for (int t = 0; t < T; ++t) if (bad[t] >= 0) Rcpp::stop("PGEN read failed at SNP %d", bad[t] + 1);
  return Rcpp::List::create(Rcpp::_["n_obs"] = Rcpp::wrap(nobs), Rcpp::_["sum"] = Rcpp::wrap(sum), Rcpp::_["dosage"] = dos);
}

// As genotype_impute_geno_cpp, but the final fills of the kept rows are written into the geno object ptr itself
// (in place: give a fresh clone nobody else references) instead of being returned; the genotypes are read from
// the same object (its missing cells, whatever they hold now). Returns the per-SNP table, info and timings.
// [[Rcpp::export]]
Rcpp::List genotype_fill_ld_geno_cpp(SEXP ptr, Rcpp::IntegerVector var_idx, const std::string& eigen_file,
                                     Rcpp::NumericVector af, Rcpp::IntegerVector sign, double thresh, double delta,
                                     double anchor_min_miss, Rcpp::LogicalVector keep_row, int threads, int nb_people,
                                     double mem_gb, Rcpp::IntegerVector rows) {
  const int T = plink4r::resolve_threads(threads);
  const auto t_all = std::chrono::steady_clock::now();
  const ImputeJob J = impute_job(eigen_file, thresh, rows, var_idx, af, sign, keep_row, t_all);
  const plink4r::ImputeLD L = plink4r::build_impute_ld(J.eb, J.rw, J.vi, J.afv, J.sg, delta, T, J.t_read);
  const auto t0 = std::chrono::steady_clock::now();
  plink4r::BlockGenotypes B;
  int ns = 0;
  plink4r::geno_genovec(ptr, J.vi, B.GS, B.sum_obs, ns, T);
  B.m = J.m;
  B.ns = ns;
  B.sb = (ns + 3) / 4;
  B.t_geno = plink4r::ik::secs(t0);
  ImputeOut O = impute_compute(J, L, B, delta, anchor_min_miss, nb_people, mem_gb, T);
  const auto tw = std::chrono::steady_clock::now();
  // CSR over the kept working rows, rows ascending within a column
  const int m = J.m;
  std::vector<int> slot(m, -1), cols;
  for (int j = 0; j < m; ++j) if (J.kr[j] && O.A.cnt1[j] > 0) { slot[j] = static_cast<int>(cols.size()); cols.push_back(J.vi[j]); }
  const std::size_t nc = cols.size();
  std::vector<std::size_t> off(nc + 1, 0);
  for (int j = 0; j < m; ++j) if (slot[j] >= 0) off[slot[j] + 1] = static_cast<std::size_t>(O.A.cnt1[j]);
  for (std::size_t c = 0; c < nc; ++c) off[c + 1] += off[c];
  std::vector<std::uint32_t> ci(off[nc]);
  std::vector<double> cg(off[nc]);
  {
    std::vector<std::size_t> pos(off.begin(), off.end() - 1);
    const int nt = static_cast<int>(O.R.oj.size());
    for (int t = 0; t < nt; ++t) {
      std::vector<int>& oi = O.R.oi[t];
      std::vector<int>& oj = O.R.oj[t];
      std::vector<float>& og = O.R.og[t];
      for (std::size_t r = 0; r < oj.size(); ++r) {
        const int s = slot[oj[r]];
        if (s < 0) continue;
        const std::size_t e = pos[s]++;
        ci[e] = static_cast<std::uint32_t>(oi[r]);
        cg[e] = plink4r::final_fill(og[r], O.A.shift[oj[r]]);
      }
      std::vector<int>().swap(oi); std::vector<int>().swap(oj); std::vector<float>().swap(og);  // release as we go
    }
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(T)
#endif
  for (std::ptrdiff_t cc = 0; cc < static_cast<std::ptrdiff_t>(nc); ++cc) {  // each thread list is in person order; merge-sort by row
    const std::size_t a = off[cc], b = off[cc + 1];
    bool sorted = true;
    for (std::size_t e = a + 1; e < b && sorted; ++e) sorted = ci[e - 1] < ci[e];
    if (sorted) continue;
    std::vector<std::pair<std::uint32_t, double>> v(b - a);
    for (std::size_t e = a; e < b; ++e) v[e - a] = {ci[e], cg[e]};
    std::sort(v.begin(), v.end());
    for (std::size_t e = a; e < b; ++e) { ci[e] = v[e - a].first; cg[e] = v[e - a].second; }
  }
  plink4r::geno_write_cells(ptr, cols, off, ci, cg, T);
  return impute_summary(J, L, B, O, delta, T, plink4r::ik::secs(tw), t_all);
}
