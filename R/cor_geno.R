.read_pvar <- function(pfile) {
  f <- paste0(pfile, ".pvar")
  if (file.exists(f)) return(data.table::fread(f, select = c("ID", "ALT"), skip = "#CHROM"))
  if (file.exists(paste0(f, ".zst")))
    return(data.table::fread(cmd = paste("zstdcat", shQuote(paste0(f, ".zst")), "| grep -v '^##'"), select = c("ID", "ALT")))
  stop("no .pvar or .pvar.zst for ", pfile)
}

#' LD correlation matrix from a PLINK2 .pgen (hardcall or dosage)
#'
#' Genotypes are read as ALT dosage (pgen dosages kept at 1/16384 resolution). Missing cells: with `af` the
#' scheme-3 constant c_j = (2 n p_j - sum_obs g) / m_j (clipped to [0, 2]); without `af` the observed mean.
#' @param pfile PLINK2 prefix. @param snps variant IDs (.pvar ID), output order.
#' @param sample IIDs (NULL = all). @param af anchor ALT frequency per SNP (NULL or NA = mean fill).
#' @return list(R = p x p matrix with dimnames, mean, sd, n_miss, fill, n, dosage, time)
#' @export
cor_pgen <- function(pfile, snps, sample = NULL, af = NULL, threads = 4L) {
  pvar <- .read_pvar(pfile)
  psam <- data.table::fread(paste0(pfile, ".psam"))
  iid <- as.character(psam[[intersect(c("IID", "#IID"), names(psam))[1]]])
  vi <- match(snps, pvar$ID)
  if (anyNA(vi)) stop(sum(is.na(vi)), " of `snps` not in ", pfile, ".pvar")
  si <- if (is.null(sample)) seq_along(iid) else match(as.character(sample), iid)
  if (anyNA(si)) stop(sum(is.na(si)), " samples not in .psam")
  act <- as.integer(2L + lengths(regmatches(pvar$ALT, gregexpr(",", pvar$ALT))))
  a <- if (is.null(af)) rep(NA_real_, length(snps)) else as.numeric(af)
  out <- cor_pgen_cpp(normalizePath(paste0(pfile, ".pgen")), length(iid), act, as.integer(vi - 1L), as.integer(si - 1L), a, as.integer(threads))
  dimnames(out$R) <- list(snps, snps)
  out
}

#' LD correlation matrix from a PLINK1 .bed (counts A1; missing as in [cor_pgen()])
#' @export
cor_bed <- function(bfile, snps, sample = NULL, af = NULL, threads = 4L) {
  bim <- data.table::fread(paste0(bfile, ".bim"), header = FALSE, select = 2L)[[1]]
  fam <- as.character(data.table::fread(paste0(bfile, ".fam"), header = FALSE, select = 2L)[[1]])
  vi <- match(snps, bim)
  if (anyNA(vi)) stop(sum(is.na(vi)), " of `snps` not in ", bfile, ".bim")
  si <- if (is.null(sample)) seq_along(fam) else match(as.character(sample), fam)
  if (anyNA(si)) stop(sum(is.na(si)), " samples not in .fam")
  a <- if (is.null(af)) rep(NA_real_, length(snps)) else as.numeric(af)
  out <- cor_bed_cpp(normalizePath(paste0(bfile, ".bed")), length(fam), as.integer(vi - 1L), as.integer(si - 1L), a, as.integer(threads))
  dimnames(out$R) <- list(snps, snps)
  out
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
