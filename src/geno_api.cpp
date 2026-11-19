// C-callable entry points for other packages' C / C++ code (inst/include/plink4r_api.h; LinkingTo: plink4r).
// Each takes the external pointer of a geno object (X$ptr) and raw column-major buffers, writes its result
// into the caller's buffer (no result is allocated here) and returns 0, or a PLINK4R_API_E* code without
// touching R: 1 not a live geno pointer, 2 bad argument (size / index), 3 failure inside (e.g. memory).
// scale is X$scale; with scale the products are those of the column-standardized matrix, from the
// store's mean and sd (= X$center, X$sd; a column with sd 0 is only centered), as in R/geno.R.
#include "geno_store.h"
#include "geno_products.h"

using namespace plink4r::gstore;

namespace {

const Geno* api_geno(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP) return nullptr;
  return static_cast<const Geno*>(R_ExternalPtrAddr(ptr));
}

bool api_cols_ok(const Geno& g, const int* cols, int nc) {
  if (nc < 0 || (nc > 0 && !cols)) return false;
  for (int a = 0; a < nc; ++a)
    if (cols[a] < 0 || static_cast<std::size_t>(cols[a]) >= g.p) return false;
  return true;
}

// .geno_div(): sd, or 1 for a constant column.
std::vector<double> api_div(const Geno& g) {
  std::vector<double> d(g.p);
  for (std::size_t j = 0; j < g.p; ++j) d[j] = g.sd[j] > 0 ? g.sd[j] : 1.0;
  return d;
}

// R's colSums (long double accumulation in row order) of an r x c column-major matrix.
void api_colsums(const double* x, std::size_t r, std::size_t c, double* out) {
  for (std::size_t j = 0; j < c; ++j) {
    long double s = 0.0;
    for (std::size_t i = 0; i < r; ++i) s += x[i + j * r];
    out[j] = static_cast<double>(s);
  }
}

template <class F>
int api_guard(F&& f) {
  try {
    return f();
  } catch (...) {
    return 3;
  }
}

// n, p, storage (1 = dosage PGEN, 0 = hard calls) and per column mean, sd (n - 1) and missing count of the
// file; any pointer may be NULL (call with NULL arrays first to size them, p values each).
int api_geno_info(SEXP ptr, int* n, int* p, int* dosage, double* mean, double* sd, double* n_miss) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (n) *n = static_cast<int>(g->n);
  if (p) *p = static_cast<int>(g->p);
  if (dosage) *dosage = g->dosage ? 1 : 0;
  if (mean) std::copy(g->mu.begin(), g->mu.end(), mean);
  if (sd) std::copy(g->sd.begin(), g->sd.end(), sd);
  if (n_miss) std::copy(g->nmiss.begin(), g->nmiss.end(), n_miss);
  return 0;
}

// out (p x k) = X'M for M (n x k); geno_crossprod(X, M).
int api_geno_xtm(SEXP ptr, int scale, const double* M, int k, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (k < 0 || (k > 0 && (!M || !out))) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t kk = static_cast<std::size_t>(k);
    xtm_into(*g, M, kk, out, threads);
    if (scale) {
      std::vector<double> cs(kk), d = api_div(*g);
      api_colsums(M, g->n, kk, cs.data());
      std_into(out, g->p, kk, g->mu.data(), cs.data(), 1.0, d.data(), nullptr, out, threads);
    }
    return 0;
  });
}

// out (n x k) = X B for B (p x k); geno_multiply(X, B).
int api_geno_mx(SEXP ptr, int scale, const double* B, int k, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (k < 0 || (k > 0 && (!B || !out))) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t n = g->n, p = g->p, kk = static_cast<std::size_t>(k);
    if (!scale) {
      mx_into(*g, B, kk, out, threads);
      return 0;
    }
    const std::vector<double> d = api_div(*g);
    std::vector<double> Bs(p * kk), Bc(p * kk), cs(kk);
    for (std::size_t c = 0; c < kk; ++c)
      for (std::size_t j = 0; j < p; ++j) {
        Bs[j + c * p] = B[j + c * p] / d[j];
        Bc[j + c * p] = Bs[j + c * p] * g->mu[j];
      }
    api_colsums(Bc.data(), p, kk, cs.data());
    mx_into(*g, Bs.data(), kk, out, threads);
    for (std::size_t c = 0; c < kk; ++c)
      for (std::size_t i = 0; i < n; ++i) out[i + c * n] = out[i + c * n] - cs[c];
    return 0;
  });
}

// out (q x q) = X'X of columns cols (0-based, q = nc), or of all p columns when cols is NULL, standardized in
// place in out; geno_crossprod(X)[cols + 1, cols + 1].
int api_geno_xtx(SEXP ptr, int scale, const int* cols, int nc, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (!out || (cols && !api_cols_ok(*g, cols, nc))) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t q = cols ? static_cast<std::size_t>(nc) : g->p;
    if (!q) return 0;
    xtx_into(*g, cols, q, out, threads);
    if (scale) {
      const std::vector<double> d = api_div(*g);
      std::vector<double> u(q), dq(q);
      for (std::size_t a = 0; a < q; ++a) {
        const std::size_t j = cols ? static_cast<std::size_t>(cols[a]) : a;
        u[a] = g->mu[j];
        dq[a] = d[j];
      }
      std_into(out, q, q, u.data(), u.data(), static_cast<double>(g->n), dq.data(), dq.data(), out, threads);
    }
    return 0;
  });
}

// out (q x q) = X'W~X of columns cols (0-based, q = nc), or of all p columns when cols is NULL, with w (n)
// binned in groups of 64 rows, or w itself when D is dense (xtwx_bin_into), standardized in place with those
// weights when scale; w_bin (n, may be NULL) receives them. geno_wcrossprod_bin(X, w)[cols + 1, cols + 1].
int api_geno_xtwx_bin(SEXP ptr, int scale, const double* w, const int* cols, int nc, double* out, double* w_bin,
                      int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (!out || !w || (cols && !api_cols_ok(*g, cols, nc))) return 2;
  for (std::size_t i = 0; i < g->n; ++i) if (!std::isfinite(w[i])) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t q = cols ? static_cast<std::size_t>(nc) : g->p;
    std::vector<double> wq(g->n);
    xtwx_bin_into(*g, w, cols, q, out, wq.data(), threads);
    if (w_bin) std::copy(wq.begin(), wq.end(), w_bin);
    if (scale) {
      const std::vector<double> d = api_div(*g);
      std::vector<double> u(q), dq(q);
      for (std::size_t a = 0; a < q; ++a) {
        const std::size_t j = cols ? static_cast<std::size_t>(cols[a]) : a;
        u[a] = g->mu[j];
        dq[a] = d[j];
      }
      xtwx_std(*g, wq.data(), cols, q, u.data(), dq.data(), out, threads);
    }
    return 0;
  });
}

// out (n x nc) = X[, cols + 1] (cols 0-based), standardized when scale.
int api_geno_cols(SEXP ptr, int scale, const int* cols, int nc, double* out, int threads) {
  const Geno* g = api_geno(ptr);
  if (!g) return 1;
  if (!api_cols_ok(*g, cols, nc) || (nc > 0 && !out)) return 2;
  return api_guard([&] {
    threads = resolve_threads(threads);
    const std::size_t n = g->n, q = static_cast<std::size_t>(nc);
    std::vector<int> rows(n);
    for (std::size_t i = 0; i < n; ++i) rows[i] = static_cast<int>(i);
    dense_into(*g, rows.data(), n, cols, q, out, threads);
    if (scale)
      for (std::size_t a = 0; a < q; ++a) {
        const std::size_t j = static_cast<std::size_t>(cols[a]);
        const double mu = g->mu[j], d = g->sd[j] > 0 ? g->sd[j] : 1.0;
        double* o = out + a * n;
        for (std::size_t i = 0; i < n; ++i) o[i] = o[i] - mu;
        for (std::size_t i = 0; i < n; ++i) o[i] = o[i] / d;
      }
    return 0;
  });
}

}  // namespace

// Called from R_init_plink4r (RcppExports.cpp).
// [[Rcpp::init]]
void plink4r_register_ccallables(DllInfo* dll) {
  (void) dll;
  R_RegisterCCallable("plink4r", "plink4r_geno_info", reinterpret_cast<DL_FUNC>(&api_geno_info));
  R_RegisterCCallable("plink4r", "plink4r_geno_xtm", reinterpret_cast<DL_FUNC>(&api_geno_xtm));
  R_RegisterCCallable("plink4r", "plink4r_geno_mx", reinterpret_cast<DL_FUNC>(&api_geno_mx));
  R_RegisterCCallable("plink4r", "plink4r_geno_xtx", reinterpret_cast<DL_FUNC>(&api_geno_xtx));
  R_RegisterCCallable("plink4r", "plink4r_geno_cols", reinterpret_cast<DL_FUNC>(&api_geno_cols));
  R_RegisterCCallable("plink4r", "plink4r_geno_xtwx_bin", reinterpret_cast<DL_FUNC>(&api_geno_xtwx_bin));
}
