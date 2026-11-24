// LD scores from the eigen factors of LD blocks (eigen_bin.h), shared by ld_score.cpp and ld_score_cross.cpp.
// With R1 = U1 diag(lam1) U1' and R2 = U2 diag(lam2) U2' (the stored, possibly truncated factors), the weighted
// LD score of SNP j over a set S of SNPs with weights a_k is
//   l_j = sum_{k in S} a_k R1_jk R2_jk = w1_j' G w2_j,  w_j = diag(lam) u_j,  G = sum_{k in S} a_k u1_k u2_k',
// so a block costs |S| k1 k2 (G) plus (output rows) k1 k2, never the m x m LD matrix. R1 = R2 gives the LD
// score of one panel (a = 1) or a stratified LD score (a = an annotation); two panels give the cross-ancestry
// score. The SNPs of R2 may count the other allele: sign s flips its rows (r2_jk -> s_j s_k r2_jk).
#ifndef PLINK4R_LD_SCORE_KERNEL_H
#define PLINK4R_LD_SCORE_KERNEL_H

#include <RcppEigen.h>
#include <algorithm>
#include <vector>

namespace plink4r {

// sum rows: (r1[k], r2[k], a[k] * s[k]) for k in S; output rows: (o1[j], o2[j], s[j]). out[j] gets l_j.
// chunk: output rows per dense product (memory: chunk x k2 doubles).
inline void weighted_ld_score(const Eigen::MatrixXf& U1, const Eigen::VectorXf& lam1, const Eigen::MatrixXf& U2,
                              const Eigen::VectorXf& lam2, const std::vector<int>& r1, const std::vector<int>& r2,
                              const std::vector<double>& as, const std::vector<int>& o1, const std::vector<int>& o2,
                              const std::vector<double>& so, int chunk, double* out) {
  const int k1 = static_cast<int>(U1.cols()), k2 = static_cast<int>(U2.cols());
  const int ns = static_cast<int>(r1.size()), no = static_cast<int>(o1.size());
  Eigen::MatrixXd A(ns, k1), Bm(ns, k2);
  for (int i = 0; i < ns; ++i) {
    A.row(i) = U1.row(r1[i]).cast<double>() * as[i];
    Bm.row(i) = U2.row(r2[i]).cast<double>();
  }
  const Eigen::MatrixXd G = A.transpose() * Bm;  // k1 x k2
  const Eigen::VectorXd l1 = lam1.cast<double>(), l2 = lam2.cast<double>();
  chunk = std::max(1, chunk);
  for (int j0 = 0; j0 < no; j0 += chunk) {
    const int nr = std::min(chunk, no - j0);
    Eigen::MatrixXd W1(nr, k1), W2(nr, k2);
    for (int j = 0; j < nr; ++j) {
      W1.row(j) = U1.row(o1[j0 + j]).cast<double>().cwiseProduct(l1.transpose());
      W2.row(j) = U2.row(o2[j0 + j]).cast<double>().cwiseProduct(l2.transpose()) * so[j0 + j];
    }
    const Eigen::MatrixXd T = W1 * G;
    for (int j = 0; j < nr; ++j) out[j0 + j] = T.row(j).dot(W2.row(j));
  }
}

}  // namespace plink4r
#endif
