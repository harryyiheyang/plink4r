#' Move score weights onto a target SNP set through the LD of an eigen.bin panel (conditional weights)
#'
#' A score `S = X b` built on the panel SNPs (e.g. SBayesRC weights on 7M SNPs) is to be computed in a cohort that
#' has only part of them (the target set). Per LD block, with O = the block's panel SNPs in the target and M = all
#' its other panel SNPs, the best linear predictor of `S` from `X_O` is `X_O w_O` with
#'
#'   w_O = b_O + Delta_O,   Delta_O = D_O^{-1} (R_OO + delta I)^{-1} R_OM D_M b_M,
#'
#' where `D = diag(sqrt(2 p (1 - p)))` (p = panel A1Freq) and `R = U diag(lam) U'` is the block's eigen LD truncated
#' at `thresh`. This is the score of the conditional-mean genotypes `E[X_M | X_O]` (the LD fill of
#' [genotype_fill()] when a whole SNP is missing for everyone), folded into the weights. A SNP in both the weights and
#' the target keeps its weight `b_O` and receives `Delta_O` on top; a target SNP without a weight can receive
#' `Delta_O` alone. The solve is the transposed summary-statistics imputation of [sumstat_impute()], in the smallest
#' exact space (missing-SNP / eigen / typed); the whole computation is one compiled call.
#'
#' Matching is by (CHR, BP, A1, A2), alleles in either order (a swapped pair flips the sign), no strand flipping;
#' the target and the weights are each matched on the panel build (BP37 or BP38) that places more of their rows.
#' Weights outside the panel are matched to the target directly (same build needed) and keep their value there;
#' those found in neither are lost. A variant repeated in the target or the weights is an error.
#'
#' Per block, `lost` = b_M' (R_MM - R_MO (R_OO + delta I)^{-1} R_OM) b_M and `total` = b' R b (standardized scale)
#' give the share of the block's score variance kept, `1 - lost / total`.
#'
#' @param target the target SNPs: data.frame or file with CHR, BP, A1, A2; other columns are carried to the output.
#' @param weights data.frame or file with CHR, BP, A1, A2 and one or more numeric weight columns (per-allele
#'   effect of A1).
#' @param ld_dir the LD panel folder (block<N>.eigen.bin and `snpinfo.parquet`).
#' @param ref panel SNP table in place of `snpinfo.parquet` (data.frame or parquet path); NULL reads the folder's.
#' @param score_cols weight columns; NULL takes every numeric column of `weights` other than CHR and BP.
#' @param thresh eigenvalue variance kept (SBayesRC default 0.995); <= 0 keeps every stored pair.
#' @param delta diagonal added to R; must be > 0.
#' @param slice rows of U_O converted to double at a time in the eigen-space and missing-space solves.
#' @param threads OpenMP threads over blocks (0 = all detected).
#' @param output optional output path (.parquet or text).
#' @return data.table of the target rows: its columns, then per score column `<col>` (the weight w), `<col>_orig`
#'   (the weight it had, 0 if none) and `<col>_delta`; `in_panel`, `in_weights`. Attributes: `block_log` (per block
#'   with weights to move: m, k, n_obs, n_miss, n_miss_weighted, method, seconds, and per score `kept_<col>` =
#'   1 - lost / total), `weight_status` (per weight row: target, moved, no_target_in_block, lost).
#' @export
weight_condition <- function(target, weights, ld_dir, ref = NULL, score_cols = NULL, thresh = 0.995, delta = 0.01,
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
  for (c in seq_along(score_cols))
    data.table::set(log, j = paste0("kept_", score_cols[c]), value = res$kept[, c])
  st <- table(factor(res$weight_status, c("target", "moved", "no_target_in_block", "lost")))
  message("weights: ", st[["target"]], " on target SNPs, ", st[["moved"]], " moved through LD, ",
          st[["no_target_in_block"]], " in blocks without target SNPs, ", st[["lost"]],
          " in neither panel nor target; ", nrow(log), " blocks solved in ",
          format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s (", res$isa, ")")
  attr(out, "block_log") <- log
  attr(out, "weight_status") <- res$weight_status
  if (!is.null(output)) {
    if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
  }
  out[]
}
