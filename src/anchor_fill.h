// The anchored fill of a column: the one value f its missing cells take so that the column's allele count is the
// anchor's, sum of the observed cells + n_miss f = 2 n af, clipped to the range of a dosage [0, 2]. Used by
// geno_fill_cpp (genotype_fill(af = )) and stage 2 of the LD fill (impute_kernel.h).
#ifndef PLINK4R_ANCHOR_FILL_H
#define PLINK4R_ANCHOR_FILL_H

#include <algorithm>

namespace plink4r {

// n = rows, af = anchor frequency of the counted allele, sum_obs = sum of the observed cells, n_miss > 0 = missing cells.
inline double anchor_fill(double n, double af, double sum_obs, double n_miss) {
  return std::min(2.0, std::max(0.0, (2.0 * n * af - sum_obs) / n_miss));
}

}  // namespace plink4r

#endif
