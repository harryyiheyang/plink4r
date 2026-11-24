# Inputs shared by the LD score regression modules (ldsc_univ.R, ldsc_bicov.R, ldsc_trans.R, ldsc_strat.R): the
# summary statistics (output of summarystat_tidy() or sumstat_impute()) and the LD scores (output of ld_score())
# of one panel are matched by Block and Index, the panel's own key, so their alleles are the panel's on both sides.

# gwas as a data.table with Z, N, Block, Index (and A1, A2 when present); `what` names it in messages.
.ldsc_gwas <- function(gwas, what = "gwas") {
  g <- .read_table(gwas)
  miss <- setdiff(c("Z", "N", "Block", "Index"), names(g))
  if (length(miss)) stop(what, " lacks ", paste(miss, collapse = ", "), ": give the output of summarystat_tidy()",
                         call. = FALSE)
  g <- g[is.finite(Z) & is.finite(N) & N > 0]
  g[, `:=`(Block = as.integer(Block), Index = as.integer(Index))]
  if (anyDuplicated(g, by = c("Block", "Index"))) stop(what, " names a panel SNP (Block, Index) more than once", call. = FALSE)
  g
}

# The rows of ld (an ld_score() table, or any table keyed by Block and Index) at the rows of g, NA where absent;
# A1 / A2 must agree where both tables carry them.
.ldsc_at <- function(ld, g, what) {
  i <- ld[g, on = .(Block, Index), which = TRUE]
  if (all(c("A1", "A2") %in% names(ld)) && all(c("A1", "A2") %in% names(g))) {
    ok <- !is.na(i)
    if (!isTRUE(all(ld$A1[i[ok]] == g$A1[ok] & ld$A2[i[ok]] == g$A2[ok])))
      stop(what, " and the summary statistics come from different panels (A1 / A2 differ at the same Block, Index)",
           call. = FALSE)
  }
  i
}

# Bootstrap group counts as ldscR draws them: one row per replicate, B = the number of contiguous groups of M SNPs.
.ldsc_counts <- function(M, nblock, n_boot) {
  if (n_boot <= 0) return(matrix(0, 0L, 0L))
  B <- ldsc_groups_cpp(as.integer(M), as.integer(nblock))
  cnt <- matrix(0, n_boot, B, dimnames = list(NULL, paste0("block", seq_len(B))))
  for (r in seq_len(n_boot)) cnt[r, ] <- tabulate(sample.int(B, B, replace = TRUE), nbins = B)
  cnt
}

# The point estimates as a one-row data.frame; with a bootstrap the replicates are kept in $bootstrap.
.ldsc_result <- function(res, cols, counts) {
  out <- data.frame(res[cols], check.names = FALSE)
  if (nrow(counts)) {
    est <- cols[!grepl("\\.se$", cols) & cols != "M"]
    b <- stats::setNames(lapply(seq_along(est), function(k) res$boot[, k]), est)
    b$block.counts <- counts
    attr(out, "bootstrap") <- b
    class(out) <- c("ldsc_bootstrap_result", class(out))
  }
  out
}

#' @export
`$.ldsc_bootstrap_result` <- function(x, name) {
  if (identical(name, "bootstrap")) return(attr(x, "bootstrap", exact = TRUE))
  NextMethod("$")
}
