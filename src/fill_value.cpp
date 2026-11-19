// fill_value_cpp: per-column fill value of the missing cells from the column counts (fill_value.h).
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include "fill_value.h"

// For each column: cnst when finite; else the anchored value clip((2 n af - s_obs) / m, lo, hi) when af is
// finite and m > 0; else the observed mean s_obs / (n - m) (`none` when nothing is observed) when
// mean_default, else NaN (= keep the current fill). af and cnst may be empty (= all NA).
// [[Rcpp::export]]
Rcpp::NumericVector fill_value_cpp(double n, const Rcpp::NumericVector& s_obs, const Rcpp::NumericVector& m,
                                   const Rcpp::NumericVector& af, const Rcpp::NumericVector& cnst,
                                   bool mean_default, double none = 0.0, double lo = 0.0, double hi = 2.0) {
  const R_xlen_t p = s_obs.size();
  if (m.size() != p || (af.size() && af.size() != p) || (cnst.size() && cnst.size() != p))
    Rcpp::stop("fill_value_cpp: lengths differ.");
  Rcpp::NumericVector f(p);
  for (R_xlen_t j = 0; j < p; ++j) {
    double v = mean_default ? plink4r::observed_mean(n, s_obs[j], m[j], none) : std::numeric_limits<double>::quiet_NaN();
    if (af.size() && std::isfinite(af[j]) && m[j] > 0)
      v = plink4r::clip_dosage(plink4r::anchor_raw(n, s_obs[j], m[j], af[j]), lo, hi);
    if (cnst.size() && std::isfinite(cnst[j])) v = cnst[j];
    f[j] = v;
  }
  return f;
}
