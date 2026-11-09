#' Impute summary statistics from eigen.bin LD blocks (SBayesRC setting, fast solver)
#'
#' z_M = R_MO (R_OO + delta I)^{-1} z_O per block, R = U diag(lam) U' truncated at `thresh`,
#' solved in the smallest exact space (missing-SNP / eigen / typed); OpenMP over blocks.
#' Output effects are standardized: BETA = Z / sqrt(N + Z^2), SE = 1 / sqrt(N + Z^2).
#' @param tidy output of [summarystat_tidy()] (data.table or file).
#' @param ref reference parquet / data.frame (ID, chr, bp37, bp38, ref, alt, AF_PanUKB, block, ld_index).
#' @param ld_dir directory with block<N>.eigen.bin (SBayesRC or plink4r layout).
#' @param thresh eigenvalue variance kept (SBayesRC default 0.995); <= 0 keeps every stored pair.
#' @param delta diagonal added to R (SBayesRC: 0.1).
#' @param output optional output path (.parquet or text).
#' @return data.table over all `ref` SNPs: ID, chr, bp37, bp38, A1, A2, freq, Z, BETA, SE, P, N, r2
#'   (r2 = 1 typed, 0 imputed, -1 block without typed SNPs).
#' @export
summarystat_impute <- function(tidy, ref, ld_dir, thresh = 0.995, delta = 0.1, threads = 4L, output = NULL) {
  x <- .read_table(tidy)
  r <- .read_table(ref)[, .(ID, chr, bp37, bp38, A1 = alt, A2 = ref, AF_PanUKB, block, ld_index)]
  data.table::setorder(r, block, ld_index)
  N_med <- stats::median(x$N)
  r <- x[, .(ID, freq, Z, N)][r, on = "ID"]
  r[, r2 := data.table::fifelse(is.finite(Z), 1, NA_real_)]
  bl <- r[, .(typed = list(ld_index[is.finite(Z)]), z = list(Z[is.finite(Z)]), rows = list(.I)), by = block]
  blocks <- bl$block; typed <- bl$typed; zt <- bl$z
  run <- lengths(typed) > 0
  message(sum(run), " of ", length(blocks), " blocks have typed SNPs; ", sum(lengths(typed)), " typed SNPs")
  t0 <- Sys.time()
  res <- sumstat_impute_cpp(ld_dir, as.integer(blocks[run]), lapply(typed[run], as.integer), zt[run],
                            thresh = thresh, delta = delta, threads = as.integer(threads))
  log <- data.table::as.data.table(attr(res, "block_log"))
  message("imputation: ", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s")
  for (a in which(run)) {
    i <- match(a, which(run))
    m <- log$m[i]
    miss <- setdiff(0:(m - 1L), typed[[a]])
    zi <- res[[i]]
    br <- bl$rows[[a]]
    w <- match(miss, r$ld_index[br])
    rows <- br[w[!is.na(w)]]
    data.table::set(r, rows, "Z", zi[!is.na(w)])
    data.table::set(r, rows, "r2", 0)
  }
  r[is.na(r2), `:=`(Z = 0, r2 = -1)]
  r[is.na(N), N := N_med]
  r[is.na(freq), freq := AF_PanUKB]
  r[, `:=`(BETA = Z / sqrt(N + Z^2), SE = 1 / sqrt(N + Z^2), P = 2 * stats::pnorm(-abs(Z)))]
  out <- r[, .(ID, chr, bp37, bp38, A1, A2, freq, Z, BETA, SE, P, N, r2)]
  message(sum(out$r2 == 1), " typed, ", sum(out$r2 == 0), " imputed, ", sum(out$r2 == -1), " in blocks without typed SNPs")
  attr(out, "block_log") <- log
  if (!is.null(output)) {
    if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
  }
  out[]
}
