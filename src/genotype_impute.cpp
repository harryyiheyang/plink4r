// Genotype imputation of hard-call missing cells from an eigen.bin LD block: R entry point.
// The algorithm (stage 1 conditional mean, stage 2 anchor) lives in impute_kernel.h, shared with geno_score.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
#include <cmath>
#include <vector>
#include "impute_kernel.h"

// Working rows (the panel rows of the block / window that are in the PGEN, increasing): rows = 0-based eigen
// row, var_idx = 0-based PGEN variant, af = PGEN-ALT frequency, sign = +1 (PGEN ALT = panel alt) or -1.
// sample_idx: 0-based PSAM rows to impute (anchor counts use these samples only).
// [[Rcpp::export]]
Rcpp::List genotype_impute_cpp(const std::string& pgen, int n_raw, Rcpp::IntegerVector allele_ct,
                               Rcpp::IntegerVector var_idx, Rcpp::IntegerVector sample_idx,
                               const std::string& eigen_file, Rcpp::NumericVector af, Rcpp::IntegerVector sign, double thresh,
                               double delta, double anchor_min_miss, Rcpp::LogicalVector keep_row, int threads,
                               int nb_people, double mem_gb, double nbr_frac, int nbr_max, Rcpp::IntegerVector rows) {
  namespace ik = plink4r::ik;
  const int T = plink4r::resolve_threads(threads);
  const auto t_all = std::chrono::steady_clock::now();

  const plink4r::EigenBlock eb = plink4r::read_eigen_block(eigen_file, thresh);
  const double t_read = ik::secs(t_all);
  const std::vector<int> rw(rows.begin(), rows.end());
  const int m = static_cast<int>(rw.size());
  if (var_idx.size() != m || af.size() != m || sign.size() != m || keep_row.size() != m)
    Rcpp::stop("rows / var_idx / af / sign / keep_row differ in length");
  for (int r : rw) if (r < 0 || r >= eb.m) Rcpp::stop("rows out of range of %s", eigen_file);
  const std::vector<int> sg(sign.begin(), sign.end());
  const std::vector<int> vi(var_idx.begin(), var_idx.end());
  const std::vector<double> afv(af.begin(), af.end());
  std::vector<char> kr(m);
  for (int j = 0; j < m; ++j) kr[j] = keep_row[j] ? 1 : 0;

  const plink4r::ImputeLD L = plink4r::build_impute_ld(eb, rw, vi, afv, sg, delta, nbr_frac, nbr_max, T, t_read);
  const std::vector<int> sidx(sample_idx.begin(), sample_idx.end());
  const std::vector<int> act(allele_ct.begin(), allele_ct.end());
  plink4r::PgenFile pf(pgen, static_cast<std::size_t>(n_raw), act, T);
  plink4r::BlockGenotypes B = plink4r::read_block_genotypes(pf, vi, sidx, T);
  const int ns = B.ns;
  const plink4r::PersonPlanes P = plink4r::build_planes(L, B, T);
  const plink4r::Stage1Result R = plink4r::impute_stage1(L, P, delta, nb_people, mem_gb, T);
  const auto t2 = std::chrono::steady_clock::now();
  const plink4r::Stage2Result A = plink4r::impute_stage2(R, B.sum_obs, afv, ns, anchor_min_miss);

  std::size_t nout = 0;
  for (std::size_t t = 0; t < R.oj.size(); ++t) for (std::size_t r = 0; r < R.oj[t].size(); ++r) if (kr[R.oj[t][r]]) ++nout;
  Rcpp::IntegerVector ri(nout), rj(nout);
  Rcpp::NumericVector rg1(nout), rg(nout);
  std::size_t p = 0;
  plink4r::for_each_cell(R, A, kr, [&](int i, int j, float g1, double g) {
    ri[p] = i + 1; rj[p] = j + 1; rg1[p] = g1; rg[p] = g; ++p;
  });
  const double t_anchor = ik::secs(t2);
  auto na = [](const std::vector<double>& x) {
    Rcpp::NumericVector out(x.size());
    for (std::size_t j = 0; j < x.size(); ++j) out[j] = std::isnan(x[j]) ? NA_REAL : x[j];
    return out;
  };
  Rcpp::IntegerVector jj(m), mj(m);
  for (int j = 0; j < m; ++j) { jj[j] = j + 1; mj[j] = static_cast<int>(A.cnt1[j]); }
  const int nblk = R.n_blocks;
  return Rcpp::List::create(
    Rcpp::_["fill"] = Rcpp::DataFrame::create(Rcpp::_["i"] = ri, Rcpp::_["j"] = rj, Rcpp::_["g1"] = rg1, Rcpp::_["g"] = rg),
    Rcpp::_["snp"] = Rcpp::DataFrame::create(Rcpp::_["j"] = jj, Rcpp::_["n_miss"] = mj, Rcpp::_["c_raw"] = na(A.craw),
                                             Rcpp::_["c"] = na(A.cc), Rcpp::_["mean_stage1"] = na(A.mean1),
                                             Rcpp::_["shift"] = Rcpp::wrap(A.shift)),
    Rcpp::_["info"] = Rcpp::List::create(Rcpp::_["m"] = m, Rcpp::_["k"] = L.k, Rcpp::_["n"] = ns, Rcpp::_["delta"] = delta,
                                         Rcpp::_["persons_solved"] = static_cast<double>(R.persons_solved), Rcpp::_["threads"] = T,
                                         Rcpp::_["n_blocks"] = nblk, Rcpp::_["people_per_block"] = static_cast<double>(ns) / nblk,
                                         Rcpp::_["cpu_util_planes"] = P.u_planes,
                                         Rcpp::_["nbr_kept_mean"] = L.kept_mean, Rcpp::_["nbr_words_frac"] = L.words_frac,
                                         Rcpp::_["cpu_util_rhs"] = R.t_rhs > 0 ? R.cpu_rhs / (R.t_rhs * T) : NA_REAL,
                                         Rcpp::_["cpu_util_cholesky"] = R.t_chol > 0 ? R.cpu_chol / (R.t_chol * T) : NA_REAL),
    Rcpp::_["time"] = Rcpp::NumericVector::create(Rcpp::_["read_eigen"] = L.t_eig, Rcpp::_["build_C"] = L.t_C,
                                                  Rcpp::_["read_geno"] = B.t_geno, Rcpp::_["planes"] = P.t_planes, Rcpp::_["lists"] = R.t_lists,
                                                  Rcpp::_["rhs_CMM"] = R.t_rhs, Rcpp::_["cholesky"] = R.t_chol,
                                                  Rcpp::_["solve"] = R.t_solve, Rcpp::_["anchor"] = t_anchor,
                                                  Rcpp::_["total"] = ik::secs(t_all)));
}
