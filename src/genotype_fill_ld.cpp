// The LD fill of genotype_fill(): the missing cells of a geno object (hard calls or dosages) filled from an eigen.bin LD block
// and written into the object. The algorithm (stage 1 conditional mean, stage 2 anchor) lives in impute_kernel.h.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>
#include "impute_kernel.h"
#include "geno_access.h"

// Working rows (the panel rows of the block / window that are columns of the object, increasing): rows = 0-based eigen row,
// var_idx = 0-based column of the object, af = frequency of the object's counted allele, sign = +1 (counted allele = panel
// alt) or -1, keep_row = whether the missing cells of the row are filled (the other rows only condition). Every row of the
// object is used; the genotypes come from its planes, its observed dosage parts (on the dosage_step grid) and its missing-cell
// rows (whatever they hold now), in its current counted-allele orientation. The final fills of the kept rows are written into
// ptr itself (in place: give a fresh copy nobody else references).
// [[Rcpp::export]]
void genotype_fill_ld_cpp(SEXP ptr, Rcpp::IntegerVector var_idx, const std::string& eigen_file, Rcpp::NumericVector af,
                          Rcpp::IntegerVector sign, double thresh, double delta, double anchor_min_miss,
                          Rcpp::LogicalVector keep_row, int threads, int nb_people, double mem_gb, Rcpp::IntegerVector rows,
                          double dosage_step) {
  const int T = plink4r::resolve_threads(threads);
  const plink4r::EigenBlock eb = plink4r::read_eigen_block(eigen_file, thresh);
  const std::vector<int> rw(rows.begin(), rows.end());
  const int m = static_cast<int>(rw.size());
  if (var_idx.size() != m || af.size() != m || sign.size() != m || keep_row.size() != m)
    Rcpp::stop("rows / var_idx / af / sign / keep_row differ in length");
  for (int r : rw) if (r < 0 || r >= eb.m) Rcpp::stop("rows out of range of %s", eigen_file);
  const std::vector<int> sg(sign.begin(), sign.end()), vi(var_idx.begin(), var_idx.end());
  const std::vector<double> afv(af.begin(), af.end());
  std::vector<char> kr(m, 0);
  for (int j = 0; j < m; ++j) kr[j] = keep_row[j] ? 1 : 0;
  const plink4r::ImputeLD L = plink4r::build_impute_ld(eb, rw, vi, afv, sg, delta, T);
  plink4r::BlockGenotypes B;
  int ns = 0;
  B.step = dosage_step;
  plink4r::geno_genovec(ptr, vi, dosage_step, B.GS, B.sum_obs, B.dptr, B.di, B.dd, ns, T);
  B.m = m;
  B.ns = ns;
  B.sb = (ns + 3) / 4;
  plink4r::PersonPlanes P = plink4r::build_planes(L, B, T);
  // People whose missing cells all sit in rows that are not kept (flank rows) need no solve: drop their lists. Every person
  // with a missing kept cell keeps the full M_i, so the kept fills are unchanged.
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
  plink4r::Stage1Result R = plink4r::impute_stage1(L, P, delta, nb_people, mem_gb, T);
  const plink4r::Stage2Result A = plink4r::impute_stage2(R, B.sum_obs, afv, ns, anchor_min_miss);
  // CSR over the kept working rows, rows ascending within a column
  std::vector<int> slot(m, -1), cols;
  for (int j = 0; j < m; ++j) if (kr[j] && A.cnt1[j] > 0) { slot[j] = static_cast<int>(cols.size()); cols.push_back(vi[j]); }
  const std::size_t nc = cols.size();
  std::vector<std::size_t> off(nc + 1, 0);
  for (int j = 0; j < m; ++j) if (slot[j] >= 0) off[slot[j] + 1] = static_cast<std::size_t>(A.cnt1[j]);
  for (std::size_t c = 0; c < nc; ++c) off[c + 1] += off[c];
  std::vector<std::uint32_t> ci(off[nc]);
  std::vector<double> cg(off[nc]);
  {
    std::vector<std::size_t> pos(off.begin(), off.end() - 1);
    const int nt = static_cast<int>(R.oj.size());
    for (int t = 0; t < nt; ++t) {
      std::vector<int>& oi = R.oi[t];
      std::vector<int>& oj = R.oj[t];
      std::vector<float>& og = R.og[t];
      for (std::size_t r = 0; r < oj.size(); ++r) {
        const int s = slot[oj[r]];
        if (s < 0) continue;
        const std::size_t e = pos[s]++;
        ci[e] = static_cast<std::uint32_t>(oi[r]);
        cg[e] = plink4r::final_fill(og[r], A.shift[oj[r]]);
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
}
