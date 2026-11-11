#' LD correlation matrix of selected variants of a PLINK2 .pgen or PLINK1 .bed
#'
#' The genotypes are held by [geno_open()] (2-bit planes H plus a sparse int16
#' deviation D) and the correlations come from [geno_cor()]: X'X by popcount
#' on the planes plus the D terms, centred and scaled in C++ from the stored
#' means and sds, so no n x p dense matrix is formed. Values are ALT (PGEN) or
#' bim A1 (BED) dosages; a PGEN with dosages is read as dosages.
#'
#' Missing calls: with `af` (the counted allele's anchor frequency) every
#' missing call of the variant takes c = (2 n af - sum of observed) /
#' (number missing), clipped to `[0, 2]` (scheme 3), so the filled column has
#' mean 2 af; `NA` entries of `af`, or `af = NULL`, use `impute`. The sd is
#' always the sample sd after filling.
#'
#' @param pfile,bfile PLINK2 / PLINK1 prefix or file path.
#' @param snps Variants, as IDs or 1-based file row indices (use indices when
#'   IDs repeat), in output order. `NULL` = all.
#' @param sample IIDs or 1-based indices; `NULL` = all samples.
#' @param af Anchor frequency of the counted allele per variant, or `NULL`.
#' @param impute Fill for variants without an anchor: `"mean"` or `"median"`.
#' @param threads OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return list(R = p x p correlation matrix with dimnames (rows and columns of
#'   variants with sd 0 are `NA`), mean, sd, n_miss (missing calls per variant
#'   before filling), n, a1 (counted allele), a2, dosage, time)
#' @export
cor_pgen <- function(pfile, snps = NULL, sample = NULL, af = NULL, threads = 4L, impute = c("mean", "median")) {
  .cor_file(list(pgenfile = pfile), snps, sample, af, match.arg(impute), threads)
}

#' @rdname cor_pgen
#' @export
cor_bed <- function(bfile, snps = NULL, sample = NULL, af = NULL, threads = 4L, impute = c("mean", "median")) {
  .cor_file(list(bedfile = bfile), snps, sample, af, match.arg(impute), threads)
}

#' @rdname cor_pgen
#' @details `pgen_cor()` and `bed_cor()` are other names of `cor_pgen()` and `cor_bed()` (same arguments and
#'   return value).
#' @export
pgen_cor <- cor_pgen

#' @rdname cor_pgen
#' @export
bed_cor <- cor_bed

.cor_file <- function(file, snps, sample, af, impute, threads, max_cells = Inf) {
  threads <- .threads(threads)
  t0 <- proc.time()[["elapsed"]]
  X <- do.call(geno_open, c(file, list(snp_vec = snps, sample_vec = sample, impute = impute,
                                       scale = FALSE, threads = threads, af = af)))
  t1 <- proc.time()[["elapsed"]]
  if (as.numeric(X$p) * X$p > max_cells) {
    stop("Requested correlation matrix is too large; subset variants first.", call. = FALSE)
  }
  R <- geno_cor(X, NULL, threads)
  info <- geno_info_cpp(X$ptr)
  list(R = R, mean = X$center, sd = X$sd, n_miss = info$n_miss, n = X$n, a1 = X$a1, a2 = X$a2,
       dosage = X$dosage, time = c(read = t1 - t0, cor = proc.time()[["elapsed"]] - t1))
}

#' Eigen-decompose an LD matrix, truncate at `thresh` of the positive eigenvalue sum, write plink4r eigen.bin
#'
#' Layout: int32 m; int32 k; float64 sumPos; float32 cutoff; float64 lam[k]; float32 scale[k]; int16 U[m*k]
#' (column-major, u = q * scale, scale_j = max|U_j| / 32767).
#' @return list(m, k, sumpos) invisibly
#' @export
write_eigen_bin <- function(R, path, thresh = 0.999) {
  e <- eigen(R, symmetric = TRUE)
  lam <- e$values; lp <- pmax(lam, 0); sumpos <- sum(lp)
  k <- min(which(cumsum(lp) / sumpos >= thresh)[1], length(lam))
  U <- e$vectors[, seq_len(k), drop = FALSE]
  sc <- apply(abs(U), 2, max) / 32767
  q <- round(sweep(U, 2, sc, "/"))
  q[q > 32767] <- 32767; q[q < -32767] <- -32767
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.integer(c(nrow(R), k)), con, size = 4L, endian = "little")
  writeBin(sumpos, con, size = 8L, endian = "little")
  writeBin(thresh, con, size = 4L, endian = "little")
  writeBin(lam[seq_len(k)], con, size = 8L, endian = "little")
  writeBin(sc, con, size = 4L, endian = "little")
  writeBin(as.integer(q), con, size = 2L, endian = "little")
  invisible(list(m = nrow(R), k = k, sumpos = sumpos))
}
