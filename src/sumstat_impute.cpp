// Summary-statistics imputation from eigen.bin, ported from SBayesRC.QC2 0.1.4 (Yihe Yang).
// Per block: z_M = R_MO (R_OO + delta I)^{-1} z_O with R = U diag(lam) U', solved in the smallest exact space
// (missing-SNP / Woodbury eigen-space / typed-SNP), float32, OpenMP over blocks.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <chrono>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "eigen_bin.h"

using Eigen::LLT;
using Eigen::MatrixXf;
using Eigen::Ref;
using Eigen::VectorXf;
using plink4r::EigenBlock;

namespace {

std::vector<int> missing_indices(int m, const std::vector<int>& typed) {
  std::vector<int> miss;
  miss.reserve(m - typed.size());
  std::size_t j = 0;
  for (int i = 0; i < m; ++i) {
    if (j < typed.size() && typed[j] == i) ++j;
    else miss.push_back(i);
  }
  return miss;
}

std::string solver_method(int n_typed, int n_missing, int k) {
  if (n_missing <= n_typed && n_missing <= k) return "missing";
  if (k < n_typed) return "eigen";
  return "typed";
}

std::vector<float> impute_missing_space(const EigenBlock& ld, const std::vector<int>& typed, const std::vector<float>& z,
                                        const std::vector<int>& missing, float d) {
  VectorXf zf = VectorXf::Zero(ld.m);
  for (std::size_t i = 0; i < typed.size(); ++i) zf[typed[i]] = z[i];
  const VectorXf dO = ld.U.transpose() * zf;
  const int nm = static_cast<int>(missing.size());
  MatrixXf UM(nm, ld.k);
  for (int i = 0; i < nm; ++i) UM.row(i) = ld.U.row(missing[i]);
  const VectorXf w = ld.lambda.array() / (ld.lambda.array() + d);
  const VectorXf rhs = UM * (w.array() * dO.array()).matrix();
  const VectorXf sw = w.array().sqrt();
  for (int j = 0; j < ld.k; ++j) UM.col(j) *= sw[j];
  MatrixXf S = MatrixXf::Identity(nm, nm);
  S.noalias() -= UM * UM.transpose();
  UM.resize(0, 0);
  LLT<Ref<MatrixXf>, Eigen::Lower> llt(S);
  if (llt.info() != Eigen::Success) throw std::runtime_error("Cholesky failed (missing-SNP system)");
  const VectorXf sol = llt.solve(rhs);
  return std::vector<float>(sol.data(), sol.data() + sol.size());
}

VectorXf eigen_coefficients(MatrixXf& UO, const VectorXf& lam, const VectorXf& z, float d) {
  const int nt = UO.rows(), k = UO.cols();
  if (k < nt) {
    const VectorXf sl = lam.array().sqrt();
    for (int j = 0; j < k; ++j) UO.col(j) *= sl[j];
    MatrixXf S = (UO.transpose() * UO) / d;
    S.diagonal().array() += 1.0f;
    LLT<Ref<MatrixXf>, Eigen::Lower> llt(S);
    if (llt.info() != Eigen::Success) throw std::runtime_error("Cholesky failed (eigen-space system)");
    const VectorXf a = llt.solve((UO.transpose() * z) / d);
    return sl.array() * a.array();
  }
  MatrixXf R = UO * lam.asDiagonal() * UO.transpose();
  R.diagonal().array() += d;
  LLT<Ref<MatrixXf>, Eigen::Lower> llt(R);
  if (llt.info() != Eigen::Success) throw std::runtime_error("Cholesky failed (typed system)");
  const VectorXf s = llt.solve(z);
  return lam.array() * (UO.transpose() * s).array();
}

std::vector<float> impute_block(const EigenBlock& ld, const std::vector<int>& typed, const std::vector<float>& z, float d,
                                std::string& method) {
  const std::vector<int> miss = missing_indices(ld.m, typed);
  method = solver_method(static_cast<int>(typed.size()), static_cast<int>(miss.size()), ld.k);
  if (miss.empty()) return {};
  if (method == "missing") return impute_missing_space(ld, typed, z, miss, d);
  const int nt = static_cast<int>(typed.size());
  MatrixXf UO(nt, ld.k);
  VectorXf zf(nt);
  for (int i = 0; i < nt; ++i) { UO.row(i) = ld.U.row(typed[i]); zf[i] = z[i]; }
  const VectorXf c = eigen_coefficients(UO, ld.lambda, zf, d);
  std::vector<float> out(miss.size());
  for (std::size_t i = 0; i < miss.size(); ++i) out[i] = ld.U.row(miss[i]).dot(c);
  return out;
}

}  // namespace

// typed_index: list of 0-based strictly increasing eigen-row indices per block; z: matching z-scores.
// Returns list of imputed z for the missing rows (all rows not in typed_index, increasing), attr block_log.
// [[Rcpp::export]]
Rcpp::List sumstat_impute_cpp(const std::string& ld_dir, Rcpp::IntegerVector blocks, Rcpp::List typed_index, Rcpp::List z,
                              double thresh, double delta, int threads) {
  const int nb = blocks.size();
  if (typed_index.size() != nb || z.size() != nb) Rcpp::stop("blocks, typed_index and z differ in length");
  std::vector<int> bv(nb);
  std::vector<std::vector<int>> tv(nb);
  std::vector<std::vector<float>> zv(nb);
  for (int i = 0; i < nb; ++i) {
    bv[i] = blocks[i];
    tv[i] = Rcpp::as<std::vector<int>>(typed_index[i]);
    zv[i] = Rcpp::as<std::vector<float>>(z[i]);
  }
  std::vector<std::vector<float>> res(nb);
  std::vector<std::string> err(nb), meth(nb);
  std::vector<int> mm(nb), kk(nb);
  std::vector<double> sec(nb);
  Eigen::setNbThreads(1);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(std::max(1, threads))
#endif
  for (int i = 0; i < nb; ++i) {
    const auto t0 = std::chrono::steady_clock::now();
    try {
      const EigenBlock ld = plink4r::read_eigen_block(plink4r::block_path(ld_dir, bv[i]), thresh);
      mm[i] = ld.m; kk[i] = ld.k;
      for (std::size_t j = 0; j < tv[i].size(); ++j)
        if (tv[i][j] < 0 || tv[i][j] >= ld.m || (j && tv[i][j] <= tv[i][j - 1])) throw std::runtime_error("bad typed_index");
      res[i] = impute_block(ld, tv[i], zv[i], static_cast<float>(delta), meth[i]);
    } catch (const std::exception& e) {
      err[i] = e.what();
    }
    sec[i] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  for (int i = 0; i < nb; ++i)
    if (!err[i].empty()) Rcpp::stop("block %d: %s", bv[i], err[i].c_str());
  Rcpp::List out(nb);
  for (int i = 0; i < nb; ++i) out[i] = Rcpp::wrap(res[i]);
  Rcpp::IntegerVector ntyped(nb);
  for (int i = 0; i < nb; ++i) ntyped[i] = static_cast<int>(tv[i].size());
  out.attr("block_log") = Rcpp::DataFrame::create(Rcpp::_["block"] = bv, Rcpp::_["m"] = mm, Rcpp::_["k"] = kk,
                                                  Rcpp::_["n_typed"] = ntyped, Rcpp::_["method"] = meth,
                                                  Rcpp::_["seconds"] = sec);
  return out;
}
