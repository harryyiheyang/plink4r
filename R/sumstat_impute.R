#' Impute summary statistics from eigen.bin LD blocks (SBayesRC setting, fast solver)
#'
#' z_M = R_MO (R_OO + delta I)^{-1} z_O per block, R = U diag(lam) U' truncated at `thresh`,
#' solved in the smallest exact space (missing-SNP / eigen / typed); OpenMP over blocks, largest first;
#' the dense algebra uses AVX-512 / AVX2 / generic code chosen at run time.
#' Output effects are standardized: BETA = Z / sqrt(N + Z^2), SE = 1 / sqrt(N + Z^2).
#'
#' `tidy` must have SNPs on every chromosome 1 to 22 (an error otherwise). X is optional: the output covers the
#' panel SNPs on the chromosomes of `tidy`, so without X in `tidy` it has no X rows.
#'
#' Imputation quality: `r2` is 1 for a typed SNP and, for an imputed SNP,
#' r2 = R_jO (R_OO + delta I)^{-1} R_Oj / R_jj in [0, 1) (the share of the SNP's LD variance explained by the
#' typed SNPs; 0 in a block without typed SNPs, where Z = 0). The imputed Z is the conditional mean
#' E[z_j | z_O], whose null variance is r2: it is already shrunk like sqrt(r2) * N(0, 1), so Z needs no
#' further adjustment. Only rows of `tidy` with `r2 == 1` (or without an
#' `r2` column) are used as typed, so re-running on an imputed file re-imputes its imputed rows from the typed ones.
#' @param tidy output of [summarystat_tidy()] with the same panel (data.table or file); its rows are placed in
#'   the panel by Block and Index.
#' @param ref panel SNP table (CHR, BP37 and/or BP38, A1, A2, A1Freq, Block, Index) as a data.frame or parquet
#'   path; NULL (default) reads `snpinfo.parquet` of `ld_dir`.
#' @param ld_dir the LD panel folder: block<N>.eigen.bin (SBayesRC or plink4r layout) and `snpinfo.parquet`.
#' @param thresh eigenvalue variance kept (SBayesRC default 0.995); <= 0 keeps every stored pair.
#' @param delta diagonal added to R (SBayesRC: 0.1); must be > 0.
#' @param output optional output path (.parquet or text).
#' @return data.table over the panel SNPs on the chromosomes of `tidy`, by Block and Index: CHR, BP37, BP38, A1,
#'   A2, A1Freq (of the GWAS where typed, else the panel's), Z, BETA, SE, P, N, r2, Block, Index.
#' @export
sumstat_impute <- function(tidy, ref = NULL, ld_dir, thresh = 0.995, delta = 0.1, threads = 4L, output = NULL) {
  threads <- .threads(threads)
  .check_delta(delta)
  x <- .read_table(tidy)
  miss <- setdiff(c("CHR", "A1", "A2", "A1Freq", "Z", "N", "Block", "Index"), names(x))
  if (length(miss)) stop("tidy lacks ", paste(miss, collapse = ", "), ": give the output of summarystat_tidy()", call. = FALSE)
  if (!"r2" %in% names(x)) x[, r2 := 1]
  n_prev <- sum(is.finite(x$Z) & !(x$r2 >= 1))
  if (n_prev) message(n_prev, " rows with r2 < 1 (imputed before) are re-imputed, not used as typed")
  x <- x[is.finite(Z) & r2 >= 1]
  x[, `:=`(Block = as.integer(Block), Index = as.integer(Index))]
  chr <- unique(.chr_code(x$CHR))
  .need_autosomes(chr, "tidy")
  r <- .with_bp(.ld_panel(ref, ld_dir))
  r <- r[CHR %in% chr]
  data.table::setorder(r, Block, Index)
  i <- r[x, on = .(Block, Index), which = TRUE]
  if (anyNA(i) || anyDuplicated(i) || !isTRUE(all(r$A1[i] == as.character(x$A1) & r$A2[i] == as.character(x$A2))))
    stop("tidy does not belong to this panel (its Block, Index, A1, A2 must name panel SNPs, each once): ",
         "make it with summarystat_tidy() and this panel", call. = FALSE)
  N_med <- stats::median(x$N)
  r[, `:=`(Z = NA_real_, N = NA_real_)]
  data.table::set(r, i, "Z", as.numeric(x$Z))
  data.table::set(r, i, "N", as.numeric(x$N))
  data.table::set(r, i, "A1Freq", as.numeric(x$A1Freq))
  r[, r2 := data.table::fifelse(is.finite(Z), 1, NA_real_)]
  bl <- r[, .(typed = list(Index[is.finite(Z)]), z = list(Z[is.finite(Z)]), rows = list(.I)), by = Block]
  blocks <- bl$Block; typed <- bl$typed; zt <- bl$z
  run <- lengths(typed) > 0
  message(sum(run), " of ", length(blocks), " blocks have typed SNPs; ", sum(lengths(typed)), " typed SNPs")
  t0 <- Sys.time()
  res <- sumstat_impute_cpp(.ld_block_file(ld_dir, blocks[run]), as.integer(blocks[run]), lapply(typed[run], as.integer),
                            zt[run], thresh = thresh, delta = delta, threads = as.integer(threads))
  log <- data.table::as.data.table(attr(res, "block_log"))
  message("imputation: ", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s (",
          log$isa[1], ")")
  for (k in seq_along(res)) {  # res[[k]]: z and r2 of every eigen row of block blocks[run][k]
    br <- bl$rows[[which(run)[k]]]
    br <- br[!is.finite(r$Z[br])]
    li <- r$Index[br] + 1L
    data.table::set(r, br, "Z", res[[k]]$z[li])
    data.table::set(r, br, "r2", pmin(res[[k]]$r2[li], 1 - 1e-6))
  }
  r[is.na(r2), `:=`(Z = 0, r2 = 0)]
  r[is.na(N), N := N_med]
  r[, `:=`(BETA = Z / sqrt(N + Z^2), SE = 1 / sqrt(N + Z^2), P = 2 * stats::pnorm(-abs(Z)))]
  out <- r[, .(CHR, BP37, BP38, A1, A2, A1Freq, Z, BETA, SE, P, N, r2, Block, Index)]
  message(sum(out$r2 >= 1), " typed, ", sum(out$r2 > 0 & out$r2 < 1), " imputed (median r2 ",
          format(stats::median(out$r2[out$r2 > 0 & out$r2 < 1]), digits = 3), "), ",
          sum(out$r2 == 0), " with r2 = 0 (incl. blocks without typed SNPs)")
  nthr <- threads; log[, threads := nthr]
  attr(out, "block_log") <- log
  if (!is.null(output)) {
    if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
  }
  out[]
}
