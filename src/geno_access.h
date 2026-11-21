// Access to a geno object (struct Geno stays in geno.cpp) for the other compilation units.
#ifndef PLINK4R_GENO_ACCESS_H
#define PLINK4R_GENO_ACCESS_H

#include <Rcpp.h>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace plink4r {

// Columns cols (0-based) of the geno object ptr as PGEN-style genovecs over all its rows: 2 bits per row,
// first row in the low bits, (n + 3) / 4 bytes per column; code = count of the counted allele (0 / 1 / 2)
// or 3 for a missing cell (whatever it holds now); an observed dosage gives its H = round(g). The rest of an
// observed dosage, on the grid step: d = round((g - H) / step) (int8), the nonzero ones per column c as rows
// di[dptr[c] .. dptr[c + 1]) (ascending) with values dd. sum_obs = counted-allele sum over the observed cells
// (H + step d). Hard calls and dosages alike; stops for an appended column.
void geno_genovec(SEXP ptr, const std::vector<int>& cols, double step, std::vector<unsigned char>& gs,
                  std::vector<double>& sum_obs, std::vector<long long>& dptr, std::vector<int>& di,
                  std::vector<std::int8_t>& dd, int& n, int threads);

// Writes per-cell fills into the geno object ptr in place (ptr must be referenced by nobody else, e.g. a fresh
// clone): for column cols[c] (0-based), rows rows[off[c]..off[c+1]) (0-based, ascending, missing cells only;
// others are ignored) take vals. Columns must be distinct. Same storage rule as geno_fill_cpp.
void geno_write_cells(SEXP ptr, const std::vector<int>& cols, const std::vector<std::size_t>& off,
                      const std::vector<std::uint32_t>& rows, const std::vector<double>& vals, int threads);

}  // namespace plink4r

#endif
