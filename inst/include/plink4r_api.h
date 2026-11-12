/* C / C++ interface to the genotype kernels of plink4r, for other packages' compiled code.
 *
 * In the dependent package: DESCRIPTION `LinkingTo: plink4r` and `Imports: plink4r`, NAMESPACE
 * `importFrom(plink4r, geno_open)` (or any import, so that plink4r is loaded), then in C / C++
 *   #include <Rcpp.h>            // or <Rinternals.h>; before this header
 *   #include <plink4r_api.h>
 * and pass the geno object's external pointer X$ptr and scale flag X$scale from R. The functions are
 * fetched from plink4r with R_GetCCallable on first use (do that first call on the R main thread); the
 * dependent package is built against this header and must be rebuilt when plink4r changes it.
 *
 * Conventions: all matrices are column-major double; column indices are 0-based; every function writes its
 * result into the caller's buffer (nothing is allocated for the result) and returns PLINK4R_API_OK or an
 * error code (no R error is raised); threads is the OpenMP thread count of plink4r's kernels (0 = from
 * plink4r_resources()). With scale != 0 the result is that of the column-standardized matrix (mean 0, sd 1
 * with n - 1; constant columns only centered), with the same values as the R functions named below.
 */
#ifndef PLINK4R_API_H
#define PLINK4R_API_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define PLINK4R_API_OK 0
#define PLINK4R_API_EPTR 1   /* not a live geno external pointer */
#define PLINK4R_API_EARG 2   /* bad size, NULL buffer or column index out of range */
#define PLINK4R_API_EFAIL 3  /* failure inside plink4r (e.g. out of memory) */

#ifdef __cplusplus
extern "C" {
#endif

/* n samples, p columns, dosage (1 = dosage PGEN, 0 = hard calls), and per column mean, sd (n - 1) and
 * number of missing calls in the file (p values each; = X$center, X$sd, X$n_miss). Any pointer may be
 * NULL: call with NULL arrays first to get n and p. */
static inline int plink4r_geno_info(SEXP ptr, int* n, int* p, int* dosage, double* mean, double* sd,
                                    double* n_miss) {
  typedef int (*fn_t)(SEXP, int*, int*, int*, double*, double*, double*);
  static fn_t fn = NULL;
  if (!fn) fn = (fn_t) R_GetCCallable("plink4r", "plink4r_geno_info");
  return fn(ptr, n, p, dosage, mean, sd, n_miss);
}

/* out (p x k) = X'M for M (n x k); geno_crossprod(X, M). */
static inline int plink4r_geno_xtm(SEXP ptr, int scale, const double* M, int k, double* out, int threads) {
  typedef int (*fn_t)(SEXP, int, const double*, int, double*, int);
  static fn_t fn = NULL;
  if (!fn) fn = (fn_t) R_GetCCallable("plink4r", "plink4r_geno_xtm");
  return fn(ptr, scale, M, k, out, threads);
}

/* out (n x k) = X B for B (p x k); geno_multiply(X, B). */
static inline int plink4r_geno_mx(SEXP ptr, int scale, const double* B, int k, double* out, int threads) {
  typedef int (*fn_t)(SEXP, int, const double*, int, double*, int);
  static fn_t fn = NULL;
  if (!fn) fn = (fn_t) R_GetCCallable("plink4r", "plink4r_geno_mx");
  return fn(ptr, scale, B, k, out, threads);
}

/* out (q x q) = X'X of the q = nc columns cols, or of all p columns when cols is NULL (nc ignored),
 * standardized in place in out; geno_crossprod(X)[cols + 1, cols + 1]. */
static inline int plink4r_geno_xtx(SEXP ptr, int scale, const int* cols, int nc, double* out, int threads) {
  typedef int (*fn_t)(SEXP, int, const int*, int, double*, int);
  static fn_t fn = NULL;
  if (!fn) fn = (fn_t) R_GetCCallable("plink4r", "plink4r_geno_xtx");
  return fn(ptr, scale, cols, nc, out, threads);
}

/* out (n x nc) = the columns cols as double; X[, cols + 1]. */
static inline int plink4r_geno_cols(SEXP ptr, int scale, const int* cols, int nc, double* out, int threads) {
  typedef int (*fn_t)(SEXP, int, const int*, int, double*, int);
  static fn_t fn = NULL;
  if (!fn) fn = (fn_t) R_GetCCallable("plink4r", "plink4r_geno_cols");
  return fn(ptr, scale, cols, nc, out, threads);
}

#ifdef __cplusplus
}
#endif

#endif
