// One eigen.bin block in its int16 form, rows restricted to given SNPs (bottom layer of sbayesrc.cpp).
// plink4r layout: q (m x k, column-major int16), scale (k), lambda (k, double), sumLambda, file cutoff, read as is.
// SBayesRC layout (float32 U): quantized per column exactly as write_eigen_bin() does (scale = max |U_.t| / 32767,
// q = round(U / scale)), so both layouts give the same in-memory form.
// Q is returned row-major over the kept rows (one SNP = ld = k int16 values), c_t = scale_t sqrt(max(lambda_t, 0)),
// so X = U Lambda^{1/2} = Q diag(c).
#ifndef PLINK4R_EIGEN_INT16_H
#define PLINK4R_EIGEN_INT16_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "eigen_bin.h"

namespace plink4r {

struct Block16 {
  int m = 0, ld = 0;           // block SNPs in the file, stored eigen dimension
  int nk = 0;                  // kept rows
  double sum_lambda = 0;       // sum of positive eigenvalues of the block (truncation base)
  float file_cut = 0;          // cutoff the file was written with
  std::vector<int16_t> Q;      // nk x ld, row-major
  std::vector<float> c;        // ld
  std::vector<double> lambda;  // ld
};

inline Block16 read_block16(const std::string& path, const std::vector<int>& keep) {
  std::ifstream in(path.c_str(), std::ios::binary | std::ios::ate);
  if (!in) throw std::runtime_error("cannot open " + path);
  const std::uint64_t size = static_cast<std::uint64_t>(in.tellg());
  in.seekg(0);
  int32_t m = 0, k = 0;
  read_n(in, &m, 1, path);
  read_n(in, &k, 1, path);
  if (m <= 0 || k <= 0 || k > m) throw std::runtime_error("bad dimensions: " + path);
  Block16 B;
  B.m = m; B.ld = k; B.nk = static_cast<int>(keep.size());
  for (std::size_t r = 0; r < keep.size(); ++r)
    if (keep[r] < 0 || keep[r] >= m || (r && keep[r] <= keep[r - 1]))
      throw std::runtime_error("ld_index outside the block or not increasing: " + path);
  B.lambda.resize(k);
  B.c.resize(k);
  B.Q.resize(static_cast<std::size_t>(B.nk) * k);
  std::vector<float> scale(k);
  if (size == 20ull + 12ull * k + 2ull * m * k) {
    read_n(in, &B.sum_lambda, 1, path);
    read_n(in, &B.file_cut, 1, path);
    read_n(in, B.lambda.data(), k, path);
    read_n(in, scale.data(), k, path);
    std::vector<int16_t> col(m);
    for (int t = 0; t < k; ++t) {
      read_n(in, col.data(), m, path);
      for (int r = 0; r < B.nk; ++r) B.Q[static_cast<std::size_t>(r) * k + t] = col[keep[r]];
    }
  } else if (size == 16ull + 4ull * k + 4ull * m * k) {
    float s = 0;
    read_n(in, &s, 1, path);
    B.sum_lambda = s;
    read_n(in, &B.file_cut, 1, path);
    std::vector<float> lf(k);
    read_n(in, lf.data(), k, path);
    for (int t = 0; t < k; ++t) B.lambda[t] = lf[t];
    std::vector<float> col(m);
    for (int t = 0; t < k; ++t) {
      read_n(in, col.data(), m, path);
      float mx = 0;
      for (int i = 0; i < m; ++i) mx = std::max(mx, std::fabs(col[i]));
      const double sc = mx / 32767.;
      scale[t] = static_cast<float>(sc);
      for (int r = 0; r < B.nk; ++r) {
        double q = sc > 0 ? std::nearbyint(col[keep[r]] / sc) : 0.;
        q = std::min(std::max(q, -32767.), 32767.);
        B.Q[static_cast<std::size_t>(r) * k + t] = static_cast<int16_t>(q);
      }
    }
  } else {
    throw std::runtime_error("unknown eigen.bin layout (size mismatch): " + path);
  }
  for (int t = 0; t < k; ++t) B.c[t] = static_cast<float>(scale[t] * std::sqrt(std::max(B.lambda[t], 0.)));
  return B;
}

}  // namespace plink4r
#endif
