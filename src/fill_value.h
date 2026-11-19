// The per-column fill value of the missing cells of one variant (the anchor): the one place this arithmetic
// lives. Used by geno refills (fill_value_cpp), the file route of genotype_fill and stage 2 of the LD
// imputation (impute_kernel.h).
#ifndef PLINK4R_FILL_VALUE_H
#define PLINK4R_FILL_VALUE_H

#include <algorithm>

namespace plink4r {

// Raw anchored value: with n rows, observed sum s_obs and m cells to fill, the value every filled cell takes
// so that the column sum is 2 n af.
inline double anchor_raw(double n, double s_obs, double m, double af) { return (2.0 * n * af - s_obs) / m; }

// A dosage clipped to [lo, hi].
inline double clip_dosage(double x, double lo, double hi) { return std::min(hi, std::max(lo, x)); }

// Observed mean over the n - m observed rows (`none` when no row is observed).
inline double observed_mean(double n, double s_obs, double m, double none) { return m < n ? s_obs / (n - m) : none; }

}  // namespace plink4r

#endif
