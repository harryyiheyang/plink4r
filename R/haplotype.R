#' Read the two haplotypes of every sample from a phased PGEN file
#'
#' For the variants and samples asked for, returns two 0 / 1 matrices (samples by variants): `h1` is 1 where the
#' first haplotype of the sample carries ALT, `h2` the same for the second. Stacked, they are the 2 n haplotypes of
#' the samples, so carriers, the four combinations of two alleles and the rank of a set of columns are counted on
#' haplotypes. A homozygote gives the same value in both; a missing call gives `NA` in both; a heterozygote whose
#' phase the file does not hold also gives `NA` in both (its two alleles cannot be placed) and is counted in
#' `n_unphased`. A multiallelic variant counts every non-REF allele as ALT. The phase is the file's own: nothing is
#' phased here. A BED file holds no phase.
#'
#' @param pgenfile PLINK 2 PGEN path or prefix (`.pvar` or `.pvar.zst`, and `.psam` alongside).
#' @param snp_vec Variants, as IDs or 1-based file rows, in the column order wanted; `NULL` takes all. For a
#'   window, rows from [geno_file()], e.g. `gf$variants[chr == 16 & pos >= lo & pos <= hi, row]`.
#' @param sample_vec Samples, as IIDs (or `FID_IID`) or 1-based file rows, in the row order wanted; `NULL` takes
#'   all.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @param file A [geno_file()] result of the same file, used instead of reading the variant and sample tables
#'   (then leave `pgenfile` `NULL`).
#' @return A list: `h1`, `h2` (integer matrices, rows named by IID, columns by variant; 1 = ALT, `NA` = missing
#'   call or heterozygote of unknown phase); `n_unphased` and `n_missing` (per variant, the samples of each kind
#'   of `NA`); `phased` (`TRUE` when no heterozygote asked for lacks its phase, i.e. `all(n_unphased == 0)`);
#'   `variants` = data.table(snp, chr, pos, a1, a2, row) with `a1` = ALT, the allele `h1` and `h2` count, and
#'   `a2` = REF; and `sample`, `fid`.
#' @export
geno_haplotype <- function(pgenfile = NULL, snp_vec = NULL, sample_vec = NULL, threads = 4L, file = NULL) {
  threads <- .threads(threads)
  gf <- .as_geno_file(file, NULL, pgenfile)
  if (gf$type != "pgen") stop("geno_haplotype() reads PGEN files; a BED file holds no phase.", call. = FALSE)
  v <- gf$variants
  s <- gf$samples
  ids <- if ("id" %in% names(v)) v$id else NULL
  vidx <- if (is.null(ids) && is.character(snp_vec)) .geno_vid_index(v, snp_vec) else
    .geno_index(snp_vec, if (is.null(ids)) v$row else ids, NULL, "snp_vec")
  sidx <- .geno_index(sample_vec, s$iid, paste(s$fid, s$iid, sep = "_"), "sample_vec")
  out <- geno_haplotype_cpp(gf$path, nrow(s), v$allele_ct, vidx - 1L, sidx - 1L, threads)
  nm <- if (is.null(ids)) .geno_vid(v[vidx]) else ids[vidx]
  dimnames(out$h1) <- dimnames(out$h2) <- list(s$iid[sidx], nm)
  list(h1 = out$h1, h2 = out$h2, n_unphased = out$n_unphased, n_missing = out$n_missing,
       phased = all(out$n_unphased == 0L),
       variants = data.table::data.table(snp = nm, chr = v$chr[vidx], pos = v$pos[vidx], a1 = v$a1[vidx],
                                         a2 = v$a2[vidx], row = vidx),
       sample = s$iid[sidx], fid = s$fid[sidx])
}
