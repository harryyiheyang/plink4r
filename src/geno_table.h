// A genotype file as opened once: where it is, what kind (BED or PGEN), its variants (variant_text.h: a .pvar line
// with several ALT alleles is one variant per ALT) and its samples. geno_file() of R builds it and every entry that
// reads the file takes it behind an external pointer, so the .bim / .pvar and the .fam / .psam are read once.
#ifndef PLINK4R_GENO_TABLE_H
#define PLINK4R_GENO_TABLE_H

#include <Rcpp.h>
#include <string>
#include "variant_text.h"

namespace plink4r {

struct GenoTable {
  std::string path;   // the .bed / .pgen
  bool bed = false;
  VariantTable variants;
  SampleTable samples;
};

// The table behind the external pointer of a geno_file() result; an error when the pointer is not of this session.
inline const GenoTable& geno_table(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP || !R_ExternalPtrAddr(ptr))
    Rcpp::stop("This geno_file() result is not valid in this R session (it holds a pointer, which is neither saved "
               "nor sent to another process): call geno_file() again.");
  return *Rcpp::XPtr<GenoTable>(ptr);
}

}  // namespace plink4r

#endif
