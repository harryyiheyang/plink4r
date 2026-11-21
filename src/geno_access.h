// Access to a geno object (struct Geno stays in geno.cpp) for the other compilation units.
#ifndef PLINK4R_GENO_ACCESS_H
#define PLINK4R_GENO_ACCESS_H

#include <Rcpp.h>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace plink4r {

// For LD imputation on columns cols (0-based) of the geno object ptr: sum_obs = ALT sum over the observed cells
// of each column (H + D, so dosages count exactly), q = missing cells per row over cols, set to 0 for a row with
// no missing cell in a column with keep set (no solve needed). Hard calls and dosages alike; stops for an
// appended column.
void geno_ld_prepare(SEXP ptr, const std::vector<int>& cols, const std::vector<char>& keep,
                     std::vector<double>& sum_obs, std::vector<int>& q, int& n);

// Rows [r0, r1) x cols of the geno object as a person-major uint16 matrix (row stride ld >= cols.size(),
// entries past cols.size() written 0): the stored value H + D in 1/16384 units, 65535 for a missing cell
// (whatever it holds now).
void geno_ld_rows(SEXP ptr, const std::vector<int>& cols, int r0, int r1, std::uint16_t* out, int ld, int threads);

// A geno object of the PGEN variants (0-based) and samples (0-based PSAM rows), hard calls or dosages as the
// file holds them (geno_open_pgen_cpp).
SEXP geno_from_pgen(const std::string& pgen, int n_raw, const std::vector<int>& allele_ct,
                    const std::vector<int>& variants, const std::vector<int>& samples, int threads);

// Writes per-cell fills into the geno object ptr in place (ptr must be referenced by nobody else, e.g. a fresh
// clone): for column cols[c] (0-based), rows rows[off[c]..off[c+1]) (0-based, ascending, missing cells only;
// others are ignored) take vals. Columns must be distinct. Same storage rule as geno_fill_cpp.
void geno_write_cells(SEXP ptr, const std::vector<int>& cols, const std::vector<std::size_t>& off,
                      const std::vector<std::uint32_t>& rows, const std::vector<double>& vals, int threads);

}  // namespace plink4r

#endif
