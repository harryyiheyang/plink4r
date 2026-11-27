// geno_haplotype_cpp: the two haplotypes of every sample at the variants of a PGEN file with phased hard calls, as
// two 0 / 1 matrices (samples x variants): h1 = ALT on the first haplotype, h2 = ALT on the second.
//
// A variant is read with its phase (PgenFile::read_phased) and each sample takes, by its hard call:
//   0 ALT        h1 = 0, h2 = 0
//   2 ALT        h1 = 1, h2 = 1
//   missing      h1 = NA, h2 = NA
//   1 ALT, phase in the file     ALT on the first haplotype (1|0): h1 = 1, h2 = 0; on the second (0|1): h1 = 0, h2 = 1
//   1 ALT, no phase in the file  h1 = NA, h2 = NA, counted in n_unphased of the variant
// so a count over the rows of h1 and h2 is a count of haplotypes and never takes a heterozygote of unknown phase for
// one. n_missing and n_unphased tell the two kinds of NA of a variant apart. A multiallelic variant counts every
// non-REF allele as ALT. Variants are independent, one OpenMP loop over them; every sample of the file is read and
// the requested rows are taken in the order asked for.
#include <Rcpp.h>
#include <cstdint>
#include <vector>
#include "pgen_reader.h"
#include "resources.h"

// [[Rcpp::export]]
Rcpp::List geno_haplotype_cpp(const std::string& path, int n_total, const Rcpp::IntegerVector& allele_ct,
                              const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples, int threads) {
  threads = plink4r::resolve_threads(threads);
  plink4r::PgenFile file(path, static_cast<std::size_t>(n_total),
                         std::vector<int>(allele_ct.begin(), allele_ct.end()), threads);
  const std::ptrdiff_t n = samples.size(), p = variants.size();
  Rcpp::IntegerMatrix h1(n, p), h2(n, p);
  Rcpp::IntegerVector n_unphased(p), n_missing(p);
  int* a1 = h1.begin();
  int* a2 = h2.begin();
  int* unph = n_unphased.begin();
  int* miss = n_missing.begin();
  const int* rows = samples.begin();
  const int* vars = variants.begin();
  int failed = 0;
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t j = 0; j < p; ++j) {
    plink4r::PgenFile::Phased ph;
    if (!file.read_phased(static_cast<std::size_t>(vars[j]), plink4r::thread_id(), ph)) {
      #pragma omp atomic write
      failed = 1;
      continue;
    }
    int* c1 = a1 + j * n;
    int* c2 = a2 + j * n;
    int nu = 0, nm = 0;
    for (std::ptrdiff_t i = 0; i < n; ++i) {
      const std::uint32_t r = static_cast<std::uint32_t>(rows[i]);
      const unsigned code = static_cast<unsigned>(ph.genovec[r >> 5] >> (2 * (r & 31u))) & 3u;
      if (code == 0) { c1[i] = 0; c2[i] = 0; }
      else if (code == 2) { c1[i] = 1; c2[i] = 1; }
      else if (code == 3) { c1[i] = NA_INTEGER; c2[i] = NA_INTEGER; ++nm; }
      else if (ph.present_ct && ((ph.present[r >> 6] >> (r & 63u)) & 1u)) {
        const int first = static_cast<int>((ph.info[r >> 6] >> (r & 63u)) & 1u);
        c1[i] = first; c2[i] = 1 - first;
      } else { c1[i] = NA_INTEGER; c2[i] = NA_INTEGER; ++nu; }
    }
    unph[j] = nu;
    miss[j] = nm;
  }
  if (failed) Rcpp::stop("Could not read the phased genotypes of " + path);
  return Rcpp::List::create(Rcpp::Named("h1") = h1, Rcpp::Named("h2") = h2,
                            Rcpp::Named("n_unphased") = n_unphased, Rcpp::Named("n_missing") = n_missing);
}
