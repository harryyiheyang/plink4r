// The null model of a mixed model with a sparse kinship K over n samples, in the one form both families use:
//   V = diag(d) + tau K,   Omega = V^-1
// (linear model: d = s2e for everyone; logistic model on the working scale: d_i = 1 / w_i).
//
// Most samples have no relative in the data: their row of K is its diagonal entry kd alone, so V is the scalar
// d_i + tau kd there. The related samples (any off-diagonal entry; a few thousand in a biobank) are held together
// as one sparse symmetric matrix V_rr = diag(d_r) + tau K_rr, factorised by a sparse Cholesky (Eigen
// SimplicialLLT, fill-reducing ordering analysed once, numeric factorisation per tau); no splitting into families.
//
// With Z = [X | y] (set_data) a fit is generalised least squares from Z' Omega Z = Z_u' diag(om) Z_u +
// |L^-1 Z_r|^2 (om_i = 1 / (d_i + tau kd), V_rr = L L'): beta, y' P y, log|X' Omega X|, log|V|.
//   solve_const(tau, d0)  d the same for everyone: Z_u' Z_u is formed once (cross), so a fit touches the related
//                         samples' rows only;
//   solve(tau, d)         d per sample: the weighted cross-product over all samples is formed at every call.
// project() then writes, once, what the association stream reads (gwas_kin_stream.cpp), with X' Omega X = La La':
//   Vb = X La^-T,  L1 = [Q | r],  Q = Omega Vb,  r = Omega (y - X beta) = P y,  om (0 at related samples),
// so that for any g:  g' P g = sum om_i g_i^2 + |L^-1 g_r|^2 - |Q' g|^2.
// One KinNull is made by a REML estimator (reml_tau_gaussian.cpp, reml_tau_binomial.cpp) and handed to the stream
// by external pointer: the kinship is built once and its factor at the estimate is the one the stream uses.
#ifndef PLINK4R_KIN_SPARSE_H
#define PLINK4R_KIN_SPARSE_H

#include <RcppEigen.h>
#include <algorithm>
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
  int threads = 1;          // of solve()

  // ---- V_rr = diag(d_r) + tau K_rr = L L'
  KinSp V;                  // pattern of Krr
  std::vector<int> dpos;    // position of each diagonal entry among the values of V
  Eigen::SimplicialLLT<KinSp> llt;
  std::vector<int> rs;      // related samples in the order of the factor
  double logdet = 0;        // log|V| (all samples)

  // ---- data Z = [X | y]
  const double* X = nullptr;
  const double* y = nullptr;
  Eigen::MatrixXd Gall, Gu, Zr, Wr, G;  // Z'Z, Z_u'Z_u ((p + 1) x (p + 1)); Z_r, L^-1 Z_r (nr x (p + 1))
  std::vector<double> Zs, sq;           // n x (p + 1), n: rows of Z scaled by sqrt(om) (solve)

  // ---- the fit
  Eigen::LLT<Eigen::MatrixXd> xox;      // X' Omega X = La La'
  Eigen::VectorXd beta;
  double quad = 0, logdet_xox = 0;      // y' P y; log|X' Omega X|
  std::vector<double> om;               // n: 1 / (d_i + tau kd), 0 at related samples

  // ---- for the association stream
  std::vector<double> L1, Vb;           // n x (p + 1): Q, r; n x p

  void init(Kin&& kin) {
    K = std::move(kin);
    n = K.n;
    V = K.Krr;
    dpos.resize(K.nr);
    rs.resize(K.nr);
    for (int a = 0; a < K.nr; ++a)
      for (int q = V.outerIndexPtr()[a]; q < V.outerIndexPtr()[a + 1]; ++q)
        if (V.innerIndexPtr()[q] == a) dpos[a] = q;
    if (K.nr == 0) return;
    llt.analyzePattern(V);
    for (int a = 0; a < K.nr; ++a) rs[llt.permutationP().indices()[a]] = K.rel[a];
  }

  // The data: X (n x pp column-major), y (n); both must stay where they are while the fits are made.
  void set_data(const double* Xp, const double* yp, int pp) {
    X = Xp; y = yp; p = pp;
    const int nr = K.nr, m = p + 1;
    Zr.resize(nr, m);
    for (int cc = 0; cc < m; ++cc) {
      const double* col = cc < p ? X + static_cast<std::size_t>(cc) * n : y;
      for (int a = 0; a < nr; ++a) Zr(a, cc) = col[K.rel[a]];
    }
  }

  // Z'Z and Z_u'Z_u of the data (for independent and solve_const).
  void cross() {
    const int m = p + 1;
    Gall.resize(m, m);
    const LaFns& LA = plink4r::la();
    LA.dgemm_tn(X, n, p, n, X, p, n, Gall.data(), m, 0);
    LA.dgemm_tn(X, n, p, n, y, 1, n, Gall.data() + static_cast<std::size_t>(p) * m, m, 0);
    double yy = 0;
    for (int i = 0; i < n; ++i) yy += y[i] * y[i];
    Gall(p, p) = yy;
    for (int a = 0; a < p; ++a) Gall(p, a) = Gall(a, p);
    Gu = Gall;
    if (K.nr > 0) Gu.noalias() -= Zr.transpose() * Zr;
  }

  // Columns of X that are not a linear combination of earlier ones (residual sum of squares above tol x its own);
  // after cross().
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

  // V = d0 I + tau K: generalised least squares of the data (after cross()); false on failure.
  bool solve_const(double tau, double d0) {
    const double v = d0 + tau * K.kd;
    if (!(v > 0)) return false;
    om.assign(n, 1.0 / v);
    for (int a = 0; a < K.nr; ++a) om[K.rel[a]] = 0.0;
    logdet = (n - K.nr) * std::log(v);
    if (!factor_rr(tau, &d0, 0)) return false;
    G = Gu / v;
    return fit();
  }

  // V = diag(d) + tau K: generalised least squares of the data; false on failure.
  bool solve(double tau, const double* d) {
    const int m = p + 1;
    const std::size_t N = n;
    om.resize(N);
    sq.resize(N);
    Zs.resize(N * m);
    logdet = 0;
    for (int i = 0; i < n; ++i) {
      const double v = d[i] + tau * K.kd;
      if (!(v > 0)) return false;
      om[i] = 1.0 / v;
      logdet += std::log(v);
    }
    for (int a = 0; a < K.nr; ++a) {
      logdet += std::log(om[K.rel[a]]);
      om[K.rel[a]] = 0.0;
    }
    if (!factor_rr(tau, d, 1)) return false;
    // Z_u' diag(om) Z_u: the rows are split between the threads, each scales its rows and forms their product
    const int T = std::max(1, std::min(threads, n / 4096 + 1));
    std::vector<Eigen::MatrixXd> Gt(T, Eigen::MatrixXd(m, m));
    const LaFns& LA = plink4r::la();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(T)
#endif
    for (int t = 0; t < T; ++t) {
      const std::size_t r0 = N * t / T, r1 = N * (t + 1) / T;
      for (std::size_t i = r0; i < r1; ++i) sq[i] = std::sqrt(om[i]);
      for (int cc = 0; cc < m; ++cc) {
        const double* col = cc < p ? X + cc * N : y;
        double* zs = &Zs[cc * N];
        for (std::size_t i = r0; i < r1; ++i) zs[i] = sq[i] * col[i];
      }
      LA.dgemm_tn(&Zs[r0], static_cast<int>(r1 - r0), m, n, &Zs[r0], m, n, Gt[t].data(), m, 0);
    }
    G = Gt[0];
    for (int t = 1; t < T; ++t) G += Gt[t];
    return fit();
  }

  // x <- Omega x (n values), at the last fit.
  void omega(double* x) const {
    const int nr = K.nr;
    Eigen::VectorXd t(nr);
    for (int q = 0; q < nr; ++q) t[q] = x[rs[q]];
    for (int i = 0; i < n; ++i) x[i] *= om[i];
    if (nr == 0) return;
    llt.matrixL().solveInPlace(t);
    llt.matrixU().solveInPlace(t);
    for (int q = 0; q < nr; ++q) x[rs[q]] = t[q];
  }

  // r = Omega (y - X beta) = P y (n values), at the last fit.
  void residual(double* r) const {
    Eigen::Map<const Eigen::MatrixXd> Xm(X, n, p);
    Eigen::Map<Eigen::VectorXd> rv(r, n);
    rv = Eigen::Map<const Eigen::VectorXd>(y, n);
    rv.noalias() -= Xm * beta;
    omega(r);
  }

  // After the fit at the estimate: L1 = [Q | r] and Vb.
  void project() {
    const int m = p + 1;
    const std::size_t N = n;
    Vb.resize(N * p);
    L1.resize(N * m);
    Eigen::Map<const Eigen::MatrixXd> Xm(X, n, p);
    Eigen::Map<Eigen::MatrixXd> Vm(Vb.data(), n, p);
    Eigen::MatrixXd Li = Eigen::MatrixXd::Identity(p, p);
    xox.matrixL().solveInPlace(Li);
    Vm.noalias() = Xm * Li.transpose();
    std::copy(Vb.begin(), Vb.end(), L1.begin());
    for (int cc = 0; cc < p; ++cc) omega(L1.data() + cc * N);
    residual(L1.data() + p * N);
  }

  // |L^-1 x|^2 for the related samples' values x (in the order of rs; overwritten).
  double quad_rel(Eigen::VectorXd& x) const {
    llt.matrixL().solveInPlace(x);
    return x.squaredNorm();
  }

 private:
  // V_rr factorised (d: per sample, stride 1, or one value, stride 0), its log-determinant added to logdet.
  bool factor_rr(double tau, const double* d, int stride) {
    if (K.nr == 0) return true;
    const double* kx = K.Krr.valuePtr();
    double* vx = V.valuePtr();
    const Eigen::Index nz = V.nonZeros();
    for (Eigen::Index q = 0; q < nz; ++q) vx[q] = tau * kx[q];
    for (int a = 0; a < K.nr; ++a) vx[dpos[a]] += d[stride ? K.rel[a] : 0];
    llt.factorize(V);
    if (llt.info() != Eigen::Success) return false;
    const KinSp& L = llt.matrixL().nestedExpression();
    for (int a = 0; a < K.nr; ++a) {
      const double l = L.coeff(a, a);
      if (!(l > 0)) return false;
      logdet += 2 * std::log(l);
    }
    return true;
  }

  // G holds the unrelated samples' part of Z' Omega Z: add the related part, then beta, quad, logdet_xox.
  bool fit() {
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
};

}  // namespace plink4r
#endif
