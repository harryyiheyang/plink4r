#' Predict the weights of SNPs a target cohort lacks onto its own SNPs through the LD of an eigen.bin panel
#'
#' A score `S = X b` built on the panel SNPs (e.g. SBayesRC weights on 7M UKB SNPs) is to be computed in a cohort
#' (e.g. AoU) with its own SNP set. SNPs in both keep their weight as it is. The weights of the SNPs the target lacks
#' (UKB-specific, M) are predicted onto the target SNPs that carry no weight (target-specific, O): per LD block,
#'
#'   w_O = D_O^{-1} (R_OO + delta I)^{-1} R_OM D_M b_M,
#'
#' where `D = diag(sqrt(2 p (1 - p)))` (p = panel A1Freq) and `R = U diag(lam) U'` is the block's eigen LD truncated
#' at `thresh`. With the shared weights fixed, this `w_O` minimises the variance of `X_M b_M - X_O w_O`, the part of
#' `S` still missing. Target-specific SNPs must be in the panel to receive weight. The solve is the transposed
#' summary-statistics imputation of [sumstat_impute()], in the smallest exact space (missing-SNP / eigen / typed);
#' every panel row outside O is in M, with weight 0 where it carries none. The whole computation is one compiled call.
#'
#' Matching is by (CHR, BP, A1, A2), alleles in either order (a swapped pair flips the sign), no strand flipping;
#' the target and the weights are each matched on the panel build (BP37 or BP38) that places more of their rows.
#' Weights outside the panel are matched to the target directly (same build needed) and keep their value there;
#' those found in neither are lost. A variant repeated in the target or the weights is an error.
#'
#' Per block, with `total` = b' R b the block's score variance (standardized scale), `kept` = 1 - lost / total with
#' `lost` = b_M' (R_MM - R_MO (R_OO + delta I)^{-1} R_OM) b_M, and `kept_drop` = 1 - b_M' R_MM b_M / total, the
#' share kept when the missing SNPs are dropped.
#'
#' @param target the target SNPs: data.frame or file with CHR, BP, A1, A2; other columns are carried to the output.
#' @param weights data.frame or file with CHR, BP, A1, A2 and one or more numeric weight columns (per-allele
#'   effect of A1).
#' @param ld_dir the LD panel folder (block<N>.eigen.bin and `snpinfo.parquet`).
#' @param ref panel SNP table in place of `snpinfo.parquet` (data.frame or parquet path); NULL reads the folder's.
#' @param score_cols weight columns; NULL takes every numeric column of `weights` other than CHR and BP.
#' @param thresh eigenvalue variance kept (the plink4r panels are cut at 0.999); <= 0 keeps every stored pair.
#' @param delta diagonal added to R; must be > 0.
#' @param slice rows of U_O converted to double at a time in the eigen-space and missing-space solves.
#' @param threads OpenMP threads over blocks (0 = all detected).
#' @param output optional output path (.parquet or text).
#' @return data.table of the target rows: its columns, then per score column `<col>` (the weight w), `<col>_orig`
#'   (the weight it had, 0 if none) and `<col>_delta` (the predicted weight, 0 on SNPs with a
#'   weight of their own); `in_panel`, `in_weights`. Attributes: `block_log` (per block with weights to move: m,
#'   k, n_obs (receiving rows), n_miss, n_miss_weighted, method, seconds, and per score `kept_<col>` and
#'   `kept_drop_<col>`), `weight_status` (per weight row: target, moved, no_target_in_block, lost).
#' @export
weight_condition <- function(target, weights, ld_dir, ref = NULL, score_cols = NULL, thresh = 0.999, delta = 0.01,
                             slice = 1024L, threads = 4L, output = NULL) {
  threads <- .threads(threads)
  .check_delta(delta)
  tg <- .read_table(target)
  wt <- .read_table(weights)
  for (x in list(list(tg, "target"), list(wt, "weights"))) {
    miss <- setdiff(c("CHR", "BP", "A1", "A2"), names(x[[1]]))
    if (length(miss)) stop(x[[2]], " lacks ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (is.null(score_cols)) score_cols <- setdiff(names(wt)[vapply(wt, is.numeric, TRUE)], c("CHR", "BP"))
  if (!length(score_cols) || !all(score_cols %in% names(wt)))
    stop("score_cols must name numeric columns of weights", call. = FALSE)
  W <- as.matrix(wt[, score_cols, with = FALSE])
  storage.mode(W) <- "double"
  if (anyNA(W)) stop("weights have NA in ", paste(score_cols, collapse = ", "), call. = FALSE)
  r <- .ld_panel(ref, ld_dir)
  bp <- function(b) if (b %in% names(r)) r[[b]] else rep(NA_integer_, nrow(r))
  blocks <- sort(unique(r$Block))
  code <- function(x) { v <- .chr_code(x); v[is.na(v)] <- 0L; v }
  t0 <- Sys.time()
  res <- weight_condition_cpp(r$CHR, bp("BP37"), bp("BP38"), r$A1, r$A2, r$A1Freq, r$Block, r$Index, blocks,
                              .ld_block_file(ld_dir, blocks), code(tg$CHR), as.integer(tg$BP), as.character(tg$A1),
                              as.character(tg$A2), code(wt$CHR), as.integer(wt$BP), as.character(wt$A1),
                              as.character(wt$A2), W, thresh = thresh, delta = delta, slice = as.integer(slice),
                              threads = as.integer(threads))
  out <- data.table::copy(tg)
  for (c in seq_along(score_cols)) {
    s <- score_cols[c]
    data.table::set(out, j = s, value = res$w[, c])
    data.table::set(out, j = paste0(s, "_orig"), value = res$beta[, c])
    data.table::set(out, j = paste0(s, "_delta"), value = res$delta[, c])
  }
  out[, `:=`(in_panel = res$target_panel, in_weights = res$target_weight)]
  log <- data.table::as.data.table(res$block_log)
  for (c in seq_along(score_cols)) {
    data.table::set(log, j = paste0("kept_", score_cols[c]), value = res$kept[, c])
    data.table::set(log, j = paste0("kept_drop_", score_cols[c]), value = res$kept_drop[, c])
  }
  st <- table(factor(res$weight_status, c("target", "moved", "no_target_in_block", "lost")))
  message("weights: ", st[["target"]], " on target SNPs, ", st[["moved"]], " moved through LD, ",
          st[["no_target_in_block"]], " in blocks without target SNPs that lack a weight, ", st[["lost"]],
          " in neither panel nor target; ", nrow(log), " blocks solved in ",
          format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s (", res$isa, ")")
  attr(out, "block_log") <- log
  attr(out, "weight_status") <- res$weight_status
  if (!is.null(output)) {
    if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
  }
  out[]
}
