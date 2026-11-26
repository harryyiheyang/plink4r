// LD-bin list of an eigen LD panel (the file read by ld_bins.h): inside each block, SNPs are grouped greedily by
// r2 > r2_min: the SNP with the most remaining neighbours opens a bin with all its remaining neighbours (ties: the
// first SNP), until every SNP is in a bin (a SNP without neighbours is a bin of its own). r is the correlation of the
// stored (truncated) eigen factor X = U Lambda^{1/2}: r_ij = x_i'x_j / sqrt(x_i'x_i x_j'x_j). The m x m matrix is
// never formed: rows are compared in chunks and only the neighbour lists are kept. OpenMP over blocks, largest first.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <queue>
#include <string>
#include <vector>
#include "eigen_bin.h"
#include "resources.h"

namespace {

std::vector<int> block_bins(const plink4r::EigenBlock& ld, double r2_min, int chunk) {
  const int m = ld.m;
  Eigen::MatrixXf X = ld.U * ld.lambda.cwiseMax(0.f).cwiseSqrt().asDiagonal();
  Eigen::VectorXf d = X.rowwise().squaredNorm();
  std::vector<std::vector<int>> nb(m);
  for (int lo = 0; lo < m; lo += chunk) {
    const int h = std::min(chunk, m - lo);
    const Eigen::MatrixXf G = X.middleRows(lo, h) * X.transpose();  // h x m
    for (int a = 0; a < h; ++a) {
      const int i = lo + a;
      const double di = d[i];
      for (int j = 0; j < m; ++j) {
        if (j == i) continue;
        const double g = G(a, j), den = di * d[j];
        if (den > 0 && g * g > r2_min * den) nb[i].push_back(j);
      }
    }
  }
  std::vector<int> cnt(m), bid(m, -1);
  std::priority_queue<std::pair<int, int>> q;  // (count, -index): most neighbours, then first index
  for (int i = 0; i < m; ++i) { cnt[i] = (int)nb[i].size(); q.push({cnt[i], -i}); }
  int nbins = 0;
  while (!q.empty()) {
    const auto [c, ni] = q.top();
    q.pop();
    const int i = -ni;
    if (bid[i] >= 0 || c != cnt[i]) continue;  // assigned, or a stale count
    std::vector<int> mem{i};
    for (int j : nb[i]) if (bid[j] < 0) mem.push_back(j);
    for (int j : mem) bid[j] = nbins;
    for (int j : mem)
      for (int l : nb[j])
        if (bid[l] < 0) { --cnt[l]; q.push({cnt[l], -l}); }
    ++nbins;
  }
  return bid;
}

}  // namespace

// files: the eigen.bin of each block; m: SNPs (eigen rows) of each block; thresh: eigen truncation (<= 0: every stored
// pair). Returns the bin id of every SNP (0-based inside its block), blocks concatenated.
// [[Rcpp::export]]
Rcpp::IntegerVector ld_bins_cpp(Rcpp::CharacterVector files, Rcpp::IntegerVector m, double thresh, double r2_min,
                                int chunk, int threads) {
  const int T = plink4r::resolve_threads(threads);
  const int nb = files.size();
  if (m.size() != nb) Rcpp::stop("one m per block");
  if (!(r2_min > 0 && r2_min < 1)) Rcpp::stop("r2 must be in (0, 1)");
  if (chunk < 1) Rcpp::stop("chunk must be >= 1");
  const std::vector<std::string> fv = Rcpp::as<std::vector<std::string>>(files);
  std::vector<int64_t> off(nb + 1, 0);
  for (int i = 0; i < nb; ++i) off[i + 1] = off[i] + m[i];
  std::vector<int> out(off[nb], -1), mv(m.begin(), m.end());
  std::vector<std::string> err(nb);
  std::vector<int> order(nb);
  for (int i = 0; i < nb; ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return mv[a] > mv[b]; });
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
#endif
  for (int o = 0; o < nb; ++o) {
    const int i = order[o];
    try {
      const plink4r::EigenBlock ld = plink4r::read_eigen_block(fv[i], thresh);
      if (ld.m != mv[i]) throw std::runtime_error("the block has " + std::to_string(ld.m) + " eigen rows, the panel table " +
                                                  std::to_string(mv[i]));
      const std::vector<int> b = block_bins(ld, r2_min, chunk);
      std::copy(b.begin(), b.end(), out.begin() + off[i]);
    } catch (const std::exception& e) {
      err[i] = fv[i] + ": " + e.what();
    }
  }
  for (const auto& e : err) if (!e.empty()) Rcpp::stop(e);
  return Rcpp::IntegerVector(out.begin(), out.end());
}
