// Dense linear-algebra kernels compiled once per instruction set and chosen at run time (la_dispatch.h).
//
// Eigen picks its SIMD width at preprocessing time from __AVX2__ / __FMA__ / __AVX512F__, which R's default
// package flags never define, so Eigen products in ordinary translation units run on SSE2. This file is
// compiled three times by src/Makevars (generic flags, -mavx2 -mfma, and AVX-512), each copy with
// PLINK4R_LA_ISA set and with the Eigen namespace renamed to Eigen_<isa>, so that the three copies of every
// Eigen template are distinct symbols and the linker never mixes them (no ODR clash between ISA copies).
// Only Eigen and plain arrays are used here (no std containers, no Rcpp): nothing in this file is an inline
// function shared with the rest of the package.
#ifndef PLINK4R_LA_ISA
#error "PLINK4R_LA_ISA must be defined (generic, avx2 or avx512)"
#endif
#define PLINK4R_CAT2(a, b) a##b
#define PLINK4R_CAT(a, b) PLINK4R_CAT2(a, b)
#define Eigen PLINK4R_CAT(Eigen_, PLINK4R_LA_ISA)
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <cstring>
#include "la_dispatch.h"

namespace {

using MatF = Eigen::MatrixXf;
using VecF = Eigen::VectorXf;
using MapM = Eigen::Map<MatF>;
using MapCM = Eigen::Map<const MatF>;
using MapCV = Eigen::Map<const VecF>;

void set_err(char* err, const char* msg) { if (err) std::snprintf(err, 256, "%s", msg); }

// Column slices of an m-column output, shared out by OpenMP over `threads`; every slice is a serial Eigen
// call (the package is built with EIGEN_DONT_PARALLELIZE, so Eigen never starts threads of its own).
template <class F>
void by_columns(int m, int threads, F f) {
  const int W = 64, ns = (m + W - 1) / W;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
#endif
  for (int s = 0; s < ns; ++s) f(s * W, (s + 1) * W < m ? W : m - s * W);
}

// Copies the strict lower triangle of C (m x m) to the upper, by column slices; each slice writes its own rows.
void mirror_lower(MapM& C, int m, int threads) {
  by_columns(m, threads, [&](int c0, int nc) {
    for (int c = c0; c < c0 + nc; ++c)
      for (int r = c + 1; r < m; ++r) C(c, r) = C(r, c);
  });
}

// C_S = I - delta (R_SS + delta I)^{-1} over the working rows S (rows[0..nrows-1]), R = U diag(lam) U'.
// S = the whole block: U diag(lam / (lam + delta)) U' as before. Otherwise two exact routes, the cheaper one by
// min(|S|, k):
//   k <= |S| (Woodbury): C_S = U_S G U_S', G = (delta Lam^{-1} + U_S'U_S)^{-1} (k x k, double Cholesky
//     G = L^{-T} L^{-1}); V = U_S L^{-T}; C_S = V V' (lower triangle, then mirrored). With every row and an
//     orthonormal U this is U diag(lam / (lam + delta)) U'. delta must be > 0.
//   |S| < k (direct): Cholesky of R_SS + delta I (|S| x |S|), then C = I - delta (R_SS + delta I)^{-1}.
// Products are split into column slices shared out by OpenMP; each slice is a serial Eigen call.
int build_C(const float* U, int ld_m, int k, const float* lambda, double delta, const int* rows, int nrows,
            float* Cout, int threads, char* err) {
  using MatD = Eigen::MatrixXd;
  MapCM Uf(U, ld_m, k);
  MapCV lam(lambda, k);
  const int m = nrows;
  const int T = threads > 0 ? threads : 1;
  if (m <= 0) { set_err(err, "no working rows"); return 1; }
  if (!(delta > 0)) { set_err(err, "delta must be > 0"); return 1; }
  MapM C(Cout, m, m);
  bool whole = m == ld_m;
  for (int r = 0; whole && r < m; ++r) whole = rows[r] == r;
  if (whole) {
    // every row of the block: C = U diag(lam / (lam + delta)) U' (exact for orthonormal U; for an int16 panel
    // the fill moves by ~1e-3, two orders below the imputation RMSE), the cheapest form
    VecF w(k);
    for (int j = 0; j < k; ++j) w[j] = lam[j] / (lam[j] + static_cast<float>(delta));
    const MatF UW = Uf * w.asDiagonal();
    // lower triangle only (half the flops of the full product), then mirrored
    by_columns(m, T, [&](int c0, int nc) {
      C.block(c0, c0, m - c0, nc).noalias() = UW.bottomRows(m - c0) * Uf.middleRows(c0, nc).transpose();
    });
    mirror_lower(C, m, T);
    return 0;
  }
  MatF US(m, k);
  for (int r = 0; r < m; ++r) US.row(r) = Uf.row(rows[r]);
  if (k <= m) {
    // G0 = U_S'U_S (k x k) in double: I - C has entries of order delta / lam_max, so C needs more than float
    // accuracy in G (a float U_S'U_S moved the fill by ~1e-4).
    const MatD USd = US.cast<double>();
    MatD G0(k, k);
    by_columns(k, T, [&](int c0, int nc) { G0.middleCols(c0, nc).noalias() = USd.transpose() * USd.middleCols(c0, nc); });
    MatF V(m, k);
    for (int j = 0; j < k; ++j) { const double l = lam[j]; G0(j, j) += delta / (l > 1e-30 ? l : 1e-30); }
    Eigen::LLT<MatD> llt(G0);
    if (llt.info() != Eigen::Success) { set_err(err, "Cholesky of delta / lam + U_S'U_S failed"); return 1; }
    // V' = L^{-1} U_S', column slices of U_S' (= row slices of U_S)
    by_columns(m, T, [&](int r0, int nr) {
      MatD Y = USd.middleRows(r0, nr).transpose();
      llt.matrixL().solveInPlace(Y);
      V.middleRows(r0, nr) = Y.transpose().cast<float>();
    });
    by_columns(m, T, [&](int c0, int nc) {
      C.block(c0, c0, m - c0, nc).noalias() = V.bottomRows(m - c0) * V.middleRows(c0, nc).transpose();
    });
    mirror_lower(C, m, T);
  } else {
    const MatF UL = US * lam.asDiagonal();
    MatF Rw(m, m);
    by_columns(m, T, [&](int c0, int nc) { Rw.middleCols(c0, nc).noalias() = UL * US.middleRows(c0, nc).transpose(); });
    Rw.diagonal().array() += static_cast<float>(delta);
    Eigen::LLT<MatF> llt(Rw);
    if (llt.info() != Eigen::Success) { set_err(err, "Cholesky of R_SS + delta I failed"); return 1; }
    by_columns(m, T, [&](int c0, int nc) {
      MatF E = MatF::Zero(m, nc);
      for (int c = 0; c < nc; ++c) E(c0 + c, c) = 1.0f;
      llt.solveInPlace(E);
      C.middleCols(c0, nc) = -static_cast<float>(delta) * E;
      for (int c = 0; c < nc; ++c) C(c0 + c, c0 + c) += 1.0f;
    });
  }
  return 0;
}

// Out (M x N, leading dim ldo) = A (M x K, leading dim lda) * B (K x N, leading dim ldb); single thread.
void sgemm_nn(const float* A, int M, int K, int lda, const float* B, int N, int ldb, float* Out, int ldo) {
  Eigen::Map<const MatF, 0, Eigen::OuterStride<>> a(A, M, K, Eigen::OuterStride<>(lda));
  Eigen::Map<const MatF, 0, Eigen::OuterStride<>> b(B, K, N, Eigen::OuterStride<>(ldb));
  Eigen::Map<MatF, 0, Eigen::OuterStride<>> o(Out, M, N, Eigen::OuterStride<>(ldo));
  o.noalias() = a * b;
}

// diag(L^-1 X)'s column squared norms: out[j] = |L^-1 X_j|^2 for a lower Cholesky factor L.
template <class LLT, class M>
void inv_norms(const LLT& llt, const M& X, float* out) {
  MatF Y = X;
  llt.matrixL().solveInPlace(Y);
  for (int j = 0; j < Y.cols(); ++j) out[j] = Y.col(j).squaredNorm();
}

// Summary-statistics imputation of one block (z_M = R_MO (R_OO + delta I)^{-1} z_O, R = U diag(lam) U'),
// solved in the smallest exact space, with the imputation quality
//   r2_j = R_jM (R_OO + delta I)^{-1} R_Oj / R_jj   (fraction of R_jj explained; 0..1).
// typed: increasing 0-based rows (nt); missing: the other rows (nm), increasing. Returns the method
// (1 missing-SNP space, 2 eigen space, 3 typed space), 0 when nothing is missing, -1 on error.
int sumstat_block(const float* U, int m, int k, const float* lambda, const int* typed, int nt, const float* z,
                  const int* miss, int nm, float d, float* zout, float* r2out, char* err) {
  if (nm == 0) return 0;
  MapCM Uf(U, m, k);
  MapCV lam(lambda, k);
  VecF Rjj(nm);
  for (int i = 0; i < nm; ++i) Rjj[i] = (Uf.row(miss[i]).array().square() * lam.transpose().array()).sum();
  const int method = (nm <= nt && nm <= k) ? 1 : (k < nt ? 2 : 3);
  if (method == 1) {
    VecF zf = VecF::Zero(m);
    for (int i = 0; i < nt; ++i) zf[typed[i]] = z[i];
    const VecF dO = Uf.transpose() * zf;
    MatF UM(nm, k);
    for (int i = 0; i < nm; ++i) UM.row(i) = Uf.row(miss[i]);
    const VecF w = lam.array() / (lam.array() + d);
    const VecF rhs = UM * (w.array() * dO.array()).matrix();
    const VecF sw = w.array().sqrt();
    for (int j = 0; j < k; ++j) UM.col(j) *= sw[j];
    MatF S = MatF::Identity(nm, nm);
    S.noalias() -= UM * UM.transpose();
    UM.resize(0, 0);
    Eigen::LLT<Eigen::Ref<MatF>, Eigen::Lower> llt(S);
    if (llt.info() != Eigen::Success) { set_err(err, "Cholesky failed (missing-SNP system)"); return -1; }
    const VecF sol = llt.solve(rhs);
    for (int i = 0; i < nm; ++i) zout[i] = sol[i];
    // S = delta (R + delta I)^{-1}_MM, so R_MO (R_OO + delta I)^{-1} R_OM = R_MM + delta I - delta S^{-1}.
    VecF sinv(nm);
    inv_norms(llt, MatF::Identity(nm, nm), sinv.data());
    for (int i = 0; i < nm; ++i) r2out[i] = Rjj[i] + d - d * sinv[i];
  } else {
    MatF UO(nt, k);
    VecF zo(nt);
    for (int i = 0; i < nt; ++i) { UO.row(i) = Uf.row(typed[i]); zo[i] = z[i]; }
    MatF UMt(k, nm);  // U_M'
    for (int i = 0; i < nm; ++i) UMt.col(i) = Uf.row(miss[i]).transpose();
    VecF c;
    if (method == 2) {
      const VecF sl = lam.array().sqrt();
      for (int j = 0; j < k; ++j) UO.col(j) *= sl[j];
      MatF S = (UO.transpose() * UO) / d;
      S.diagonal().array() += 1.0f;
      Eigen::LLT<Eigen::Ref<MatF>, Eigen::Lower> llt(S);
      if (llt.info() != Eigen::Success) { set_err(err, "Cholesky failed (eigen-space system)"); return -1; }
      const VecF a = llt.solve((UO.transpose() * zo) / d);
      c = sl.array() * a.array();
      // r_j = UO_s v_j with v_j = lam^{1/2} u_j:  quality = |v_j|^2 - |L^{-1} v_j|^2, |v_j|^2 = R_jj.
      const MatF V = sl.asDiagonal() * UMt;
      VecF q(nm);
      inv_norms(llt, V, q.data());
      for (int i = 0; i < nm; ++i) r2out[i] = Rjj[i] - q[i];
    } else {
      MatF R = UO * lam.asDiagonal() * UO.transpose();
      R.diagonal().array() += d;
      Eigen::LLT<Eigen::Ref<MatF>, Eigen::Lower> llt(R);
      if (llt.info() != Eigen::Success) { set_err(err, "Cholesky failed (typed system)"); return -1; }
      const VecF s = llt.solve(zo);
      c = lam.array() * (UO.transpose() * s).array();
      const MatF ROM = UO * (lam.asDiagonal() * UMt);  // R_OM
      VecF q(nm);
      inv_norms(llt, ROM, q.data());
      for (int i = 0; i < nm; ++i) r2out[i] = q[i];
    }
    for (int i = 0; i < nm; ++i) zout[i] = Uf.row(miss[i]).dot(c);
  }
  for (int i = 0; i < nm; ++i) {
    float r = Rjj[i] > 0 ? r2out[i] / Rjj[i] : 0.0f;
    r2out[i] = r < 0 ? 0.0f : (r > 1 ? 1.0f : r);
  }
  return method;
}

}  // namespace

#define PLINK4R_LA_TABLE PLINK4R_CAT(plink4r_la_, PLINK4R_LA_ISA)
#define PLINK4R_STR2(x) #x
#define PLINK4R_STR(x) PLINK4R_STR2(x)
extern const plink4r::LaFns PLINK4R_LA_TABLE;
const plink4r::LaFns PLINK4R_LA_TABLE = {build_C, sgemm_nn, sumstat_block, PLINK4R_STR(PLINK4R_LA_ISA)};
