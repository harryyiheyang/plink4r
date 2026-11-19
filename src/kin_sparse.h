// Sparse kinship and the covariance V = diag(d) + tau K of a mixed model over n samples.
//
// Most samples have no relative in the data: their row of K is its diagonal entry alone, so V and its inverse
// Omega are scalars there. The related samples (any off-diagonal entry; a few thousand in a biobank) are held
// together as one sparse symmetric matrix V_rr, factorised by a sparse Cholesky (Eigen SimplicialLLT, fill-
// reducing ordering analysed once, numeric factorisation per tau); no splitting into families.
//   Kin        the kinship: diagonal for everyone, K_rr for the related samples
//   KinFactor  V = diag(d) + tau K factorised: Omega x, x_r' V_rr^-1 x_r, log|V|
//   gls_fit    generalised least squares of y on X under V: beta, Omega (y - X beta) = P y, and the
//              projection factors Q = Omega X L^-T, Vb = X L^-T (X' Omega X = L L'), so that for any g
//              g' P g = g' Omega g - |Q' g|^2.
// Shared by the REML estimators of tau (reml_tau_gaussian.cpp, reml_tau_binomial.cpp) and the association
// stream (gwas_kin_stream.cpp).
#ifndef PLINK4R_KIN_SPARSE_H
#define PLINK4R_KIN_SPARSE_H

#include <RcppEigen.h>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include "la_dispatch.h"

namespace plink4r {

typedef Eigen::SparseMatrix<double> KinSp;

struct Kin {
  int n = 0, nr = 0;
  std::vector<int> rel;     // related samples (0-based, ascending)
  std::vector<int> pos;     // per sample: its row of K_rr, or -1
  std::vector<double> kd;   // diagonal of K, per sample
  KinSp Krr;                // nr x nr, both triangles and the diagonal
  std::size_t pairs = 0;    // off-diagonal pairs
};

// ki, kj (0-based samples, ki != kj, each pair once), kv: the nk off-diagonal entries; kdiag: n diagonal entries.
inline Kin build_kin(int n, const int* ki, const int* kj, const double* kv, std::size_t nk, const double* kdiag) {
  Kin K;
  K.n = n;
  K.kd.assign(kdiag, kdiag + n);
  K.pos.assign(n, -1);
  for (std::size_t e = 0; e < nk; ++e) {
    if (ki[e] < 0 || ki[e] >= n || kj[e] < 0 || kj[e] >= n || ki[e] == kj[e])
      throw std::runtime_error("kinship pair out of range");
    if (kv[e] == 0) continue;
    K.pos[ki[e]] = 0;
    K.pos[kj[e]] = 0;
  }
  for (int i = 0; i < n; ++i)
    if (K.pos[i] == 0) { K.pos[i] = static_cast<int>(K.rel.size()); K.rel.push_back(i); }
  K.nr = static_cast<int>(K.rel.size());
  std::vector<Eigen::Triplet<double>> tr;
  tr.reserve(2 * nk + K.nr);
  for (std::size_t e = 0; e < nk; ++e) {
    if (kv[e] == 0) continue;
    const int a = K.pos[ki[e]], b = K.pos[kj[e]];
    tr.emplace_back(a, b, kv[e]);
    tr.emplace_back(b, a, kv[e]);
    ++K.pairs;
  }
  for (int a = 0; a < K.nr; ++a) tr.emplace_back(a, a, K.kd[K.rel[a]]);
  K.Krr.resize(K.nr, K.nr);
  K.Krr.setFromTriplets(tr.begin(), tr.end());
  K.Krr.makeCompressed();
  return K;
}

struct KinFactor {
  const Kin* K = nullptr;
  std::vector<double> om;   // 1 / (d_i + tau k_ii) for an unrelated sample, 0 for a related one
  Eigen::SimplicialLLT<KinSp> llt;
  bool analysed = false;
  double logdet = 0;

  // V = diag(d) + tau K; false when V is not positive definite.
  bool factor(const Kin& kin, const double* d, double tau) {
    K = &kin;
    const int n = kin.n, nr = kin.nr;
    om.assign(n, 0.0);
    logdet = 0;
    for (int i = 0; i < n; ++i) {
      if (kin.pos[i] >= 0) continue;
      const double v = d[i] + tau * kin.kd[i];
      if (!(v > 0)) return false;
      om[i] = 1.0 / v;
      logdet += std::log(v);
    }
    if (nr == 0) return true;
    KinSp V = kin.Krr * tau;
    for (int a = 0; a < nr; ++a) V.coeffRef(a, a) += d[kin.rel[a]];
    if (!analysed) { llt.analyzePattern(V); analysed = true; }
    llt.factorize(V);
    if (llt.info() != Eigen::Success) return false;
    const KinSp L = llt.matrixL();
    for (int a = 0; a < nr; ++a) {
      const double l = L.coeff(a, a);
      if (!(l > 0)) return false;
      logdet += 2 * std::log(l);
    }
    return true;
  }

  // out = Omega x (n values each).
  void apply(const double* x, double* out) const {
    const int n = K->n, nr = K->nr;
    for (int i = 0; i < n; ++i) out[i] = om[i] * x[i];
    if (nr == 0) return;
    Eigen::VectorXd xr(nr);
    for (int a = 0; a < nr; ++a) xr[a] = x[K->rel[a]];
    const Eigen::VectorXd s = llt.solve(xr);
    for (int a = 0; a < nr; ++a) out[K->rel[a]] = s[a];
  }

  // x_r' V_rr^-1 x_r for the related samples' values xr (nr values, in the order of Kin::rel).
  double quad_rel(const Eigen::VectorXd& xr) const {
    if (K->nr == 0) return 0.0;
    const Eigen::VectorXd z = llt.matrixL().solve(llt.permutationP() * xr);
    return z.squaredNorm();
  }
};

struct GlsFit {
  bool ok = false;
  int p = 0;
  double quad = 0;        // (y - X beta)' Omega (y - X beta) = y' P y
  double logdet_xox = 0;  // log |X' Omega X|
  std::vector<double> beta, r, Q, Vb;  // r = Omega (y - X beta) (n); Q, Vb n x p column-major (when proj)
};

// X: n x p column-major, full column rank; y: n values. proj: also Q and Vb.
inline GlsFit gls_fit(const KinFactor& F, const double* X, int n, int p, const double* y, bool proj) {
  GlsFit out;
  out.p = p;
  const std::size_t N = n;
  std::vector<double> OX(N * p);
  for (int c = 0; c < p; ++c) F.apply(X + c * N, &OX[c * N]);
  Eigen::MatrixXd A(p, p);
  la().dgemm_tn(X, n, p, n, OX.data(), p, n, A.data(), p, 0);
  A = (0.5 * (A + A.transpose())).eval();
  Eigen::LLT<Eigen::MatrixXd> llt(A);
  if (llt.info() != Eigen::Success) return out;
  Eigen::Map<const Eigen::VectorXd> yv(y, n);
  Eigen::Map<const Eigen::MatrixXd> Xm(X, n, p), OXm(OX.data(), n, p);
  const Eigen::VectorXd b = llt.solve(OXm.transpose() * yv);
  const Eigen::MatrixXd L = llt.matrixL();
  for (int c = 0; c < p; ++c) {
    if (!(L(c, c) > 0)) return out;
    out.logdet_xox += 2 * std::log(L(c, c));
  }
  const Eigen::VectorXd res = yv - Xm * b;
  out.beta.assign(b.data(), b.data() + p);
  out.r.resize(N);
  F.apply(res.data(), out.r.data());
  out.quad = res.dot(Eigen::Map<const Eigen::VectorXd>(out.r.data(), n));
  if (proj) {
    out.Q.resize(N * p);
    out.Vb.resize(N * p);
    Eigen::Map<Eigen::MatrixXd> Q(out.Q.data(), n, p), Vb(out.Vb.data(), n, p);
    Q = llt.matrixL().solve(OXm.transpose()).transpose();
    Vb = llt.matrixL().solve(Xm.transpose()).transpose();
  }
  out.ok = true;
  return out;
}

}  // namespace plink4r
#endif
