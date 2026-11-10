#' Impute summary statistics from eigen.bin LD blocks (SBayesRC setting, fast solver)
#'
#' z_M = R_MO (R_OO + delta I)^{-1} z_O per block, R = U diag(lam) U' truncated at `thresh`,
#' solved in the smallest exact space (missing-SNP / eigen / typed); OpenMP over blocks, largest first;
#' the dense algebra uses AVX-512 / AVX2 / generic code chosen at run time.
#' Output effects are standardized: BETA = Z / sqrt(N + Z^2), SE = 1 / sqrt(N + Z^2).
#'
#' Imputation quality: `r2` is 1 for a typed SNP and, for an imputed SNP,
#' r2 = R_jO (R_OO + delta I)^{-1} R_Oj / R_jj in [0, 1) (the share of the SNP's LD variance explained by the
#' typed SNPs; 0 in a block without typed SNPs, where Z = 0). The imputed Z is the conditional mean
#' E[z_j | z_O], whose null variance is r2: it is already shrunk like sqrt(r2) * N(0, 1), so Z needs no
#' further adjustment. Only rows of `tidy` with `r2 == 1` (or without an
#' `r2` column) are used as typed, so re-running on an imputed file re-imputes its imputed rows from the typed ones.
#' @param tidy output of [summarystat_tidy()] (data.table or file).
#' @param ref reference parquet / data.frame (ID, chr, bp37, bp38, ref, alt, AF_PanUKB, block, ld_index).
#' @param ld_dir directory with block<N>.eigen.bin (SBayesRC or plink4r layout).
#' @param thresh eigenvalue variance kept (SBayesRC default 0.995); <= 0 keeps every stored pair.
#' @param delta diagonal added to R (SBayesRC: 0.1).
#' @param output optional output path (.parquet or text).
#' @return data.table over all `ref` SNPs: ID, chr, bp37, bp38, A1, A2, freq, Z, BETA, SE, P, N, r2.
#' @export
sumstat_impute <- function(tidy, ref, ld_dir, thresh = 0.995, delta = 0.1, threads = 4L, output = NULL) {
  threads <- .threads(threads)
  x <- .read_table(tidy)
  if (!"r2" %in% names(x)) x[, r2 := 1]
  n_prev <- sum(is.finite(x$Z) & !(x$r2 >= 1))
  if (n_prev) message(n_prev, " rows with r2 < 1 (imputed before) are re-imputed, not used as typed")
  x <- x[is.finite(Z) & r2 >= 1]
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
  message("imputation: ", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s (",
          log$isa[1], ")")
  for (a in which(run)) {
    i <- match(a, which(run))
    m <- log$m[i]
    miss <- setdiff(0:(m - 1L), typed[[a]])
    zi <- res[[i]]$z; qi <- res[[i]]$r2
    br <- bl$rows[[a]]
    w <- match(miss, r$ld_index[br])
    rows <- br[w[!is.na(w)]]
    data.table::set(r, rows, "Z", zi[!is.na(w)])
    data.table::set(r, rows, "r2", pmin(qi[!is.na(w)], 1 - 1e-6))
  }
  r[is.na(r2), `:=`(Z = 0, r2 = 0)]
  r[is.na(N), N := N_med]
  r[is.na(freq), freq := AF_PanUKB]
  r[, `:=`(BETA = Z / sqrt(N + Z^2), SE = 1 / sqrt(N + Z^2), P = 2 * stats::pnorm(-abs(Z)))]
  out <- r[, .(ID, chr, bp37, bp38, A1, A2, freq, Z, BETA, SE, P, N, r2)]
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

#' @rdname sumstat_impute
#' @details `summarystat_impute()` is the former name, kept as an alias.
#' @export
summarystat_impute <- function(...) sumstat_impute(...)
