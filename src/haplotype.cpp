// geno_haplotype_cpp: the two haplotypes of every sample at the variants of a PGEN file with phased hard calls, as
// two 0 / 1 matrices (samples x variants): h1 = the variant's counted allele (A1: one ALT allele of its line) on the
// first haplotype, h2 = on the second.
//
// A line is read with its phase (PgenFile::read_haplotypes: the allele on each haplotype of every sample) and a
// variant, which counts the ALT allele a of that line, takes for each sample:
//   both alleles known and placed (a homozygote, or a heterozygote whose phase the file holds)
//                                  h1 = 1 if the first haplotype carries a, else 0; h2 the same for the second
//   missing call                   h1 = NA, h2 = NA, counted in n_missing of the variant
//   heterozygote without phase     neither of its alleles is a: h1 = 0, h2 = 0;
//                                  one of them is a: h1 = NA, h2 = NA, counted in n_unphased of the variant
// so a count over the rows of h1 and h2 is a count of haplotypes and never takes a heterozygote of unknown phase for
// one. n_missing and n_unphased tell the two kinds of NA of a variant apart. Variants are independent, one OpenMP
// loop over them (the variants of one line follow one another and are served by one read); every sample of the
// file is read and the requested rows are taken in the order asked for.
#include <Rcpp.h>
#include <cstdint>
#include <vector>
#include "geno_table.h"
#include "pgen_reader.h"
#include "resources.h"

// file: the genotype file as opened (geno_table.h), a PGEN; variants: 0-based variants of its table, the columns
// wanted; samples: 0-based file rows, the rows wanted.
// [[Rcpp::export]]
Rcpp::List geno_haplotype_cpp(SEXP file, const Rcpp::IntegerVector& variants, const Rcpp::IntegerVector& samples,
                              int threads) {
  using plink4r::PgenFile;
  threads = plink4r::resolve_threads(threads);
  const plink4r::GenoTable& gt = plink4r::geno_table(file);
  if (gt.bed) Rcpp::stop("A BED file holds no phase: " + gt.path);
  PgenFile pf(gt.path, gt.samples.size(), gt.variants.line_alleles, threads);
  const std::ptrdiff_t n = samples.size(), p = variants.size();
  Rcpp::IntegerMatrix h1(n, p), h2(n, p);
  Rcpp::IntegerVector n_unphased(p), n_missing(p);
  int* a1 = h1.begin();
  int* a2 = h2.begin();
  int* unph = n_unphased.begin();
  int* miss = n_missing.begin();
  const int* rows = samples.begin();
  const int* vars = variants.begin();
  int failed = -1;
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (std::ptrdiff_t j = 0; j < p; ++j) {
    const std::int32_t line = gt.variants.line[vars[j]];
    const unsigned char a = static_cast<unsigned char>(gt.variants.alt[vars[j]]);
    PgenFile::Haplotypes hp;
    if (!pf.read_haplotypes(static_cast<std::size_t>(line), plink4r::thread_id(), hp)) {
      #pragma omp atomic write
      failed = line;
      continue;
    }
    int* c1 = a1 + j * n;
    int* c2 = a2 + j * n;
    int nu = 0, nm = 0;
    for (std::ptrdiff_t i = 0; i < n; ++i) {
      const int r = rows[i];
      const unsigned char st = hp.state[r];
      if (st == PgenFile::kPhased) { c1[i] = hp.h1[r] == a; c2[i] = hp.h2[r] == a; }
      else if (st == PgenFile::kNoCall) { c1[i] = NA_INTEGER; c2[i] = NA_INTEGER; ++nm; }
      else if (hp.h1[r] != a && hp.h2[r] != a) { c1[i] = 0; c2[i] = 0; }
      else { c1[i] = NA_INTEGER; c2[i] = NA_INTEGER; ++nu; }
    }
    unph[j] = nu;
    miss[j] = nm;
  }
  if (failed >= 0) Rcpp::stop(pf.failure(static_cast<std::size_t>(failed)));
  return Rcpp::List::create(Rcpp::Named("h1") = h1, Rcpp::Named("h2") = h2,
                            Rcpp::Named("n_unphased") = n_unphased, Rcpp::Named("n_missing") = n_missing);
}
