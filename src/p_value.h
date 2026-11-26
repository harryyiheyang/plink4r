// The P values that the association scans write beside their statistics: R's own distribution functions (Rmath),
// called where a row is written. They are called from several threads at once, so every argument that would make R
// raise a warning is answered here: a missing statistic or no degree of freedom gives NaN.
#ifndef PLINK4R_P_VALUE_H
#define PLINK4R_P_VALUE_H

#include <Rcpp.h>
#include <cmath>
#include <limits>

namespace plink4r {

// Two-sided P of a t statistic with df degrees of freedom: 2 * pt(-|t|, df).
inline double p_t(double t, double df) {
  if (std::isnan(t) || !(df > 0)) return std::numeric_limits<double>::quiet_NaN();
  return 2.0 * R::pt(-std::fabs(t), df, 1, 0);
}
// Two-sided P of a z statistic: 2 * pnorm(-|z|).
inline double p_z(double z) {
  if (std::isnan(z)) return std::numeric_limits<double>::quiet_NaN();
  return 2.0 * R::pnorm(-std::fabs(z), 0.0, 1.0, 1, 0);
}
// Upper tail of a chi-square statistic with df degrees of freedom: pchisq(x, df, lower.tail = FALSE).
inline double p_chisq(double x, double df) {
  if (std::isnan(x) || !(df > 0)) return std::numeric_limits<double>::quiet_NaN();
  return R::pchisq(x, df, 0, 0);
}

}  // namespace plink4r

#endif
