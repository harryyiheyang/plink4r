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
#' share kept when the missing SNPs are dropped. Both are variance ratios (not r2); `kept_drop` can be negative when
#' the dropped SNPs and the kept ones are negatively correlated in their contributions.
#'
#' Duplicates: a variant (CHR, BP, unordered allele pair) repeated within the target or the weights is left out with
#' every copy (`dup = "drop"`, status "duplicate") or stops the call (`dup = "error"`). No target or no weight SNP in
#' the panel stops the call (wrong build or chromosome coding). Off-panel weights reach the target only when the two
#' tables are on the same build.
#'
#' The output is keyed by the target's own columns, among them `CHR`, `BP`, `A1`, `A2`: it is the `weights` table
#' of [geno_score()] as it is.
#'
#' @param target the target SNPs: data.frame or file with CHR, BP, A1, A2; other columns are carried to the output.
#' @param weights data.frame or file with CHR, BP, A1, A2 and the weight columns (per-allele effect of A1).
#' @param ld_dir the LD panel folder (block<N>.eigen.bin and `snpinfo.parquet`).
#' @param score_cols the weight columns of `weights` (finite numbers; e.g. "BETA").
#' @param ref panel SNP table in place of `snpinfo.parquet` (data.frame or parquet path); NULL reads the folder's.
#' @param thresh eigenvalue variance kept; <= 0 (default) keeps every stored pair, whatever the panel's cut.
#' @param delta diagonal added to R; must be > 0.
#' @param dup "drop" or "error", for variants repeated within the target or the weights.
#' @param slice rows of U_O converted to double at a time in the eigen-space and missing-space solves.
#' @param threads OpenMP threads over blocks (0 = all detected).
#' @param output optional output path (.parquet or text).
#' @return data.table of the target rows: its columns, then per score column `<col>` (the weight w), `<col>_orig`
#'   (the weight it had, 0 if none) and `<col>_delta` (the predicted weight, 0 on SNPs with a weight of their own);
#'   `in_panel`, `in_weights`. Attributes: `block_log` (per block with weights to move: m, k, n_obs (receiving rows),
#'   n_miss, n_miss_weighted, method, seconds, and per score `kept_<col>` and `kept_drop_<col>`), `weight_status`
#'   (per weight row: target, moved, no_target_in_block, no_freq, duplicate, lost), `match` (panel build used for
#'   the target and the weights, rows matched, duplicates left out).
#' @export
weight_condition <- function(target, weights, ld_dir, score_cols, ref = NULL, thresh = 0, delta = 0.01,
                             dup = c("drop", "error"), slice = 1024L, threads = 4L, output = NULL) {
  dup <- match.arg(dup)
  threads <- .threads(threads)
  .check_delta(delta)
  tg <- .read_table(target)
  wt <- .read_table(weights)
  for (x in list(list(tg, "target"), list(wt, "weights"))) {
    miss <- setdiff(c("CHR", "BP", "A1", "A2"), names(x[[1]]))
    if (length(miss)) stop(x[[2]], " lacks ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!length(score_cols) || !all(score_cols %in% names(wt)) ||
      !all(vapply(wt[, score_cols, with = FALSE], is.numeric, TRUE)))
    stop("score_cols must name numeric columns of weights", call. = FALSE)
  W <- as.matrix(wt[, score_cols, with = FALSE])
  storage.mode(W) <- "double"
  bad <- colSums(!is.finite(W))
  if (any(bad > 0)) stop("weights have non-finite values in ", paste(score_cols[bad > 0], collapse = ", "), call. = FALSE)
  clash <- intersect(names(tg), c(score_cols, paste0(score_cols, "_orig"), paste0(score_cols, "_delta"), "in_panel",
                                  "in_weights"))
  if (length(clash)) stop("target already has the output columns ", paste(clash, collapse = ", "), call. = FALSE)
  r <- .ld_panel(ref, ld_dir)
  bp <- function(b) if (b %in% names(r)) r[[b]] else rep(NA_integer_, nrow(r))
  blocks <- sort(unique(r$Block))
  code <- function(x) { v <- .chr_code(x); v[is.na(v)] <- 0L; v }
  t0 <- Sys.time()
  res <- weight_condition_cpp(r$CHR, bp("BP37"), bp("BP38"), r$A1, r$A2, r$A1Freq, r$Block, r$Index, blocks,
                              .ld_block_file(ld_dir, blocks), code(tg$CHR), as.integer(tg$BP), as.character(tg$A1),
                              as.character(tg$A2), code(wt$CHR), as.integer(wt$BP), as.character(wt$A1),
                              as.character(wt$A2), W, thresh = thresh, delta = delta, slice = as.integer(slice),
                              dup = as.integer(dup == "drop"), threads = as.integer(threads))
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
  st <- table(factor(res$weight_status, c("target", "moved", "no_target_in_block", "no_freq", "duplicate", "lost")))
  b <- c("BP37", "BP38")
  match <- list(build_target = b[res$build_target + 1L], build_weights = b[res$build_weight + 1L],
                n_target_panel = res$n_target_panel, n_weights_panel = res$n_weight_panel,
                target_duplicate = sum(res$target_dup), weights_duplicate = st[["duplicate"]])
  message("target: ", match$n_target_panel, " of ", nrow(tg), " in the panel (", match$build_target, "); weights: ",
          match$n_weights_panel, " of ", nrow(wt), " in the panel (", match$build_weights, ")",
          if (match$target_duplicate + match$weights_duplicate > 0)
            paste0("; left out as repeated: ", match$target_duplicate, " target, ", match$weights_duplicate, " weight rows"))
  message("weights: ", st[["target"]], " on target SNPs, ", st[["moved"]], " moved through LD, ",
          st[["no_target_in_block"]], " in blocks without target SNPs that lack a weight, ", st[["no_freq"]],
          " without a panel frequency, ", st[["lost"]], " in neither panel nor target; ", nrow(log), " blocks solved in ",
          format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s (", res$isa, ")")
  attr(out, "block_log") <- log
  attr(out, "weight_status") <- res$weight_status
  attr(out, "match") <- match
  if (!is.null(output)) {
    if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
  }
  out[]
}
