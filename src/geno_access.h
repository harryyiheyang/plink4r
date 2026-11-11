// Read access to a geno object (struct Geno stays in geno.cpp) for the other compilation units.
#ifndef PLINK4R_GENO_ACCESS_H
#define PLINK4R_GENO_ACCESS_H

#include <Rcpp.h>
#include <vector>

namespace plink4r {

// Columns cols (0-based) of the geno object ptr as PGEN-style genovecs over all its rows: 2 bits per row,
// first row in the low bits, (n + 3) / 4 bytes per column; code = count of the counted allele (0 / 1 / 2)
// or 3 for a missing cell (whatever it holds now). sum_obs = counted-allele sum over the observed cells.
// Stops for a dosage object, an appended column, or a non-integer value at an observed cell.
void geno_genovec(SEXP ptr, const std::vector<int>& cols, std::vector<unsigned char>& gs,
                  std::vector<double>& sum_obs, int& n, int threads);

}  // namespace plink4r

#endif
