// The null model with a sparse kinship K over n samples, in the one form all families use:
//   Omega = V^-1 = S H(h)^-1 S,   H(h) = (1 - h) I + h K,   S = diag(s)
// (linear model: s = 1 / sigma; logistic model: s_i = sqrt(w_i) of the working response; Cox model: the same on
// the records of the life table, summed over the records of a sample, see reml_tau_cox.cpp).
//
// Most samples have no relative in the data: their row of K is its diagonal entry kd alone, so H is the scalar
// 1 - h + h kd there. The related samples (any off-diagonal entry; a few thousand in a biobank) are held together
// as one sparse symmetric matrix H_rr, factorised by a sparse Cholesky (Eigen SimplicialLLT, fill-reducing
// ordering analysed once, numeric factorisation per h); no splitting into families.
//
// Nothing of size n is redone per h. With the scaled data Z = [S X | S y] (set_data, once per working data):
//   Gu = Z_u' Z_u (unrelated samples),  Zr = Z_r (related rows),
// and at any h (solve), with c = 1 / (1 - h + h kd) and H_rr = L L':
//   Z' H^-1 Z = c Gu + |L^-1 Zr|^2  ->  X' Omega X, X' Omega y, y' Omega y  ->  beta, y' P y, log|X' Omega X|,
//   log|H| = (n - nr) log(1 - h + h kd) + 2 sum log L_ii.
// project() then writes, once, what the association stream reads (gwas_kin_stream.cpp), with X' Omega X = La La':
//   Vb = X La^-T,  L1 = [Q | r],  Q = Omega Vb,  r = Omega (y - X beta) = P y,  om_i = s_i^2 c (0 if related),
// so that for any g:  g' P g = sum om_i g_i^2 + |L^-1 (s_r o g_r)|^2 - |Q' g|^2.
// One KinNull is made by a null-model estimator (reml_tau_gaussian.cpp, reml_tau_binomial.cpp, reml_tau_cox.cpp,
// cox_null.cpp) and handed to the stream by external pointer: the kinship is built once and its factor at the
// estimate is the one the stream uses. The stream reads om, L1 = [Q | r], rs, ss, the factor (quad_rel) and, when
// the estimator has one, Vb (empty: no sample-level basis of the projected columns; see reml_tau_cox.cpp).
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
  double kd = 1;            // diagonal of K
  std::vector<int> rel;     // related samples (0-based, ascending)
  KinSp Krr;                // nr x nr, both triangles and the diagonal
  std::size_t pairs = 0;    // off-diagonal pairs
};

// ki, kj (0-based samples, ki != kj, each pair once), kv: the nk off-diagonal entries; kd: the diagonal.
inline Kin build_kin(int n, const int* ki, const int* kj, const double* kv, std::size_t nk, double kd) {
  Kin K;
  K.n = n;
  K.kd = kd;
  std::vector<int> pos(n, -1);
  for (std::size_t e = 0; e < nk; ++e) {
    if (ki[e] < 0 || ki[e] >= n || kj[e] < 0 || kj[e] >= n || ki[e] == kj[e])
      throw std::runtime_error("kinship pair out of range");
    if (kv[e] == 0) continue;
    pos[ki[e]] = 0;
    pos[kj[e]] = 0;
  }
  for (int i = 0; i < n; ++i)
    if (pos[i] == 0) { pos[i] = static_cast<int>(K.rel.size()); K.rel.push_back(i); }
  K.nr = static_cast<int>(K.rel.size());
  std::vector<Eigen::Triplet<double>> tr;
  tr.reserve(2 * nk + K.nr);
  for (std::size_t e = 0; e < nk; ++e) {
    if (kv[e] == 0) continue;
    const int a = pos[ki[e]], b = pos[kj[e]];
    tr.emplace_back(a, b, kv[e]);
    tr.emplace_back(b, a, kv[e]);
    ++K.pairs;
  }
  for (int a = 0; a < K.nr; ++a) tr.emplace_back(a, a, kd);
  K.Krr.resize(K.nr, K.nr);
  K.Krr.setFromTriplets(tr.begin(), tr.end());
  K.Krr.makeCompressed();
  return K;
}

struct KinNull {
  Kin K;
  int n = 0, p = 0;

  // ---- H_rr(h) = L L'
  KinSp H;                  // pattern of Krr, values of H_rr(h)
  std::vector<int> dpos;    // position of each diagonal entry among the values of H
  Eigen::SimplicialLLT<KinSp> llt;
  double h = 0, c = 1, logdet = 0;  // c = 1 / (1 - h + h kd); logdet = log|H(h)|

  // ---- scaled data
  Eigen::MatrixXd Gall, Gu, Zr, Wr, G;  // Z'Z, Z_u'Z_u ((p + 1) x (p + 1)); Z_r, L^-1 Z_r (nr x (p + 1))

  // ---- the fit at h
  Eigen::LLT<Eigen::MatrixXd> xox;      // X' S H^-1 S X = La La'
  Eigen::VectorXd beta;
  double quad = 0, logdet_xox = 0;      // (S y)' P_H (S y); log|X' S H^-1 S X|

  // ---- for the association stream
  std::vector<double> om, L1, Vb;       // n; n x (p + 1): Q, r; n x p
  std::vector<int> rs;                  // related samples in the order of the factor
  std::vector<double> ss;               // their scale s

  void init(Kin&& kin) {
    K = std::move(kin);
    n = K.n;
    H = K.Krr;
    dpos.resize(K.nr);
    for (int a = 0; a < K.nr; ++a)
      for (int q = H.outerIndexPtr()[a]; q < H.outerIndexPtr()[a + 1]; ++q)
        if (H.innerIndexPtr()[q] == a) dpos[a] = q;
    if (K.nr > 0) llt.analyzePattern(H);
  }

  // H(h) factorised; false when it is not positive definite.
  bool factor(double hh) {
    h = hh;
    const double v = 1.0 - h + h * K.kd;
    if (!(v > 0)) return false;
    c = 1.0 / v;
    logdet = (n - K.nr) * std::log(v);
    if (K.nr == 0) return true;
    const double* kx = K.Krr.valuePtr();
    double* hx = H.valuePtr();
    const Eigen::Index nz = H.nonZeros();
    for (Eigen::Index q = 0; q < nz; ++q) hx[q] = h * kx[q];
    for (int a = 0; a < K.nr; ++a) hx[dpos[a]] += 1.0 - h;
    llt.factorize(H);
    if (llt.info() != Eigen::Success) return false;
    const KinSp& L = llt.matrixL().nestedExpression();
    for (int a = 0; a < K.nr; ++a) {
      const double l = L.coeff(a, a);
      if (!(l > 0)) return false;
      logdet += 2 * std::log(l);
    }
    return true;
  }

  // The upper end of the search for h: the largest h <= h_max at which H(h) is positive definite. The eigenvalues
  // of H(h) are 1 + h (lambda - 1), lambda those of K, so H is positive definite on [0, 1 / (1 - lambda_min)) and
  // the end is found by bisection on the factorisation (to tol).
  double h_bound(double h_max, double tol) {
    if (factor(h_max)) return h_max;
    double lo = 0.0, hi = h_max;
    while (hi - lo > tol) {
      const double mid = 0.5 * (lo + hi);
      if (factor(mid)) lo = mid; else hi = mid;
    }
    return lo;
  }

  // The scaled data: Xs = S X (n x pp column-major), ys = S y.
  void set_data(const double* Xs, const double* ys, int pp) {
    p = pp;
    const int nr = K.nr, m = p + 1;
    Gall.resize(m, m);
    const LaFns& LA = plink4r::la();
    LA.dgemm_tn(Xs, n, p, n, Xs, p, n, Gall.data(), m, 0);
    LA.dgemm_tn(Xs, n, p, n, ys, 1, n, Gall.data() + static_cast<std::size_t>(p) * m, m, 0);
    double yy = 0;
    for (int i = 0; i < n; ++i) yy += ys[i] * ys[i];
    Gall(p, p) = yy;
    for (int a = 0; a < p; ++a) Gall(p, a) = Gall(a, p);
    Zr.resize(nr, m);
    for (int cc = 0; cc < m; ++cc) {
      const double* col = cc < p ? Xs + static_cast<std::size_t>(cc) * n : ys;
      for (int a = 0; a < nr; ++a) Zr(a, cc) = col[K.rel[a]];
    }
    Gu = Gall;
    if (nr > 0) Gu.noalias() -= Zr.transpose() * Zr;
  }

  // The scaled data given by its cross-products, for a model whose rows are not the samples (reml_tau_cox.cpp):
  // ZtZ = Z'Z over all rows ((pp + 1) x (pp + 1)) and Zrel (nr x (pp + 1)), the related samples' rows of the part of
  // Z on which H acts, so that Z' H^-1 Z = (ZtZ - Zrel' Zrel) + |L^-1 Zrel|^2.
  void set_cross(const Eigen::MatrixXd& ZtZ, const Eigen::MatrixXd& Zrel, int pp) {
    p = pp;
    Gall = ZtZ;
    Zr = Zrel;
    Gu = Gall;
    if (K.nr > 0) Gu.noalias() -= Zr.transpose() * Zr;
  }

  // Columns of X that are not a linear combination of earlier ones (residual sum of squares above tol x its own).
  std::vector<int> independent(double tol) const {
    std::vector<int> keep;
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(p, p);
    Eigen::VectorXd v(p);
    for (int j = 0; j < p; ++j) {
      const int m = static_cast<int>(keep.size());
      double d = Gall(j, j);
      for (int a = 0; a < m; ++a) {
        double s = Gall(keep[a], j);
        for (int b = 0; b < a; ++b) s -= L(a, b) * v[b];
        v[a] = s / L(a, a);
        d -= v[a] * v[a];
      }
      if (!(d > tol * Gall(j, j))) continue;
      for (int a = 0; a < m; ++a) L(m, a) = v[a];
      L(m, m) = std::sqrt(d);
      keep.push_back(j);
    }
    return keep;
  }

  // Generalised least squares of the scaled data at h: beta, quad, logdet_xox (and logdet); false on failure.
  bool solve(double hh) {
    if (!factor(hh)) return false;
    G = c * Gu;
    if (K.nr > 0) {
      Wr = llt.permutationP() * Zr;
      llt.matrixL().solveInPlace(Wr);
      G.noalias() += Wr.transpose() * Wr;
    }
    xox.compute(G.topLeftCorner(p, p));
    if (xox.info() != Eigen::Success) return false;
    logdet_xox = 0;
    for (int a = 0; a < p; ++a) {
      const double l = xox.matrixLLT()(a, a);
      if (!(l > 0)) return false;
      logdet_xox += 2 * std::log(l);
    }
    beta = xox.solve(G.col(p).head(p));
    quad = G(p, p) - G.col(p).head(p).dot(beta);
    return true;
  }

  // After solve(h) at the estimate: om, L1 = [Q | r], Vb for S = a diag(s0) (s0 = nullptr: 1), from the unscaled
  // X (n x p) and y. set_data must have been given diag(s0) X, diag(s0) y.
  void project(const double* X, const double* y, const double* s0, double a) {
    const int nr = K.nr, m = p + 1;
    const std::size_t N = n;
    om.resize(N);
    for (int i = 0; i < n; ++i) om[i] = a * a * c * (s0 ? s0[i] * s0[i] : 1.0);
    rs.resize(nr);
    ss.resize(nr);
    for (int q = 0; q < nr; ++q) {
      const int pos = llt.permutationP().indices()[q];
      rs[pos] = K.rel[q];
      ss[pos] = a * (s0 ? s0[K.rel[q]] : 1.0);
    }
    Vb.resize(N * p);
    L1.resize(N * m);
    Eigen::Map<const Eigen::MatrixXd> Xm(X, n, p);
    Eigen::Map<const Eigen::VectorXd> yv(y, n);
    Eigen::Map<Eigen::MatrixXd> Vm(Vb.data(), n, p), Lm(L1.data(), n, m);
    Eigen::MatrixXd Li = Eigen::MatrixXd::Identity(p, p);
    xox.matrixL().solveInPlace(Li);
    Vm.noalias() = Xm * (Li.transpose() / a);
    Lm.leftCols(p) = Vm;
    Lm.col(p) = yv;
    Lm.col(p).noalias() -= Xm * beta;
    // Omega applied to the columns of L1
    Eigen::MatrixXd T(nr, m);
    for (int cc = 0; cc < m; ++cc) {
      double* col = L1.data() + cc * N;
      for (int q = 0; q < nr; ++q) T(q, cc) = ss[q] * col[rs[q]];
      for (std::size_t i = 0; i < N; ++i) col[i] *= om[i];
    }
    if (nr > 0) {
      llt.matrixL().solveInPlace(T);
      llt.matrixU().solveInPlace(T);
    }
    for (int cc = 0; cc < m; ++cc) {
      double* col = L1.data() + cc * N;
      for (int q = 0; q < nr; ++q) col[rs[q]] = ss[q] * T(q, cc);
    }
    for (int q = 0; q < nr; ++q) om[rs[q]] = 0.0;
  }

  // |L^-1 x|^2 for the related samples' scaled values x (in the order of rs; overwritten).
  double quad_rel(Eigen::VectorXd& x) const {
    llt.matrixL().solveInPlace(x);
    return x.squaredNorm();
  }
};

}  // namespace plink4r
#endif
