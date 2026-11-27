// The entry of geno_file(): a genotype file opened once (geno_table.h), and its two tables as R shows them.
#include <Rcpp.h>
#include <string>
#include <vector>
#include "geno_table.h"

// path: the .bed / .pgen; bed: whether a BED; variant_text: its .bim / .pvar (plain text); sample_text: its .fam /
// .psam. Reads the two texts and returns the table behind an external pointer.
// [[Rcpp::export]]
SEXP geno_table_open_cpp(const std::string& path, bool bed, const std::string& variant_text,
                         const std::string& sample_text) {
  Rcpp::XPtr<plink4r::GenoTable> t(new plink4r::GenoTable, true);
  t->path = path;
  t->bed = bed;
  t->variants = plink4r::read_variant_table(variant_text, bed);
  t->samples = plink4r::read_sample_table(sample_text, bed);
  if (t->variants.size() == 0) Rcpp::stop("No variant in " + variant_text);
  if (t->samples.size() == 0) Rcpp::stop("No sample in " + sample_text);
  return t;
}

// The variants of an opened file as columns: CHR (chromosome code, NA without one), BP (NA when not a position),
// and A1, A2 as 1-based codes into `alleles`, the distinct allele texts of the file in upper case.
// [[Rcpp::export]]
Rcpp::List geno_table_variants_cpp(SEXP ptr) {
  const plink4r::VariantTable& v = plink4r::geno_table(ptr).variants;
  const R_xlen_t m = static_cast<R_xlen_t>(v.size());
  Rcpp::IntegerVector chr(m), bp(m), a1(m), a2(m);
  for (R_xlen_t r = 0; r < m; ++r) {
    const std::int32_t c = v.code(r);
    chr[r] = c < 0 ? NA_INTEGER : c;
    bp[r] = v.bp[r] < 0 ? NA_INTEGER : v.bp[r];
    a1[r] = v.a1[r] + 1; a2[r] = v.a2[r] + 1;
  }
  Rcpp::CharacterVector alleles(v.alleles.size());
  for (std::size_t k = 0; k < v.alleles.size(); ++k) {
    std::string u = v.alleles[k];
    for (char& ch : u) ch = plink4r::ascii_upper(ch);
    alleles[k] = u;
  }
  return Rcpp::List::create(Rcpp::Named("CHR") = chr, Rcpp::Named("BP") = bp, Rcpp::Named("A1") = a1,
                            Rcpp::Named("A2") = a2, Rcpp::Named("alleles") = alleles);
}

// The samples of an opened file: fid and iid, in file order.
// [[Rcpp::export]]
Rcpp::List geno_table_samples_cpp(SEXP ptr) {
  const plink4r::SampleTable& s = plink4r::geno_table(ptr).samples;
  return Rcpp::List::create(Rcpp::Named("fid") = Rcpp::wrap(s.fid), Rcpp::Named("iid") = Rcpp::wrap(s.iid));
}
