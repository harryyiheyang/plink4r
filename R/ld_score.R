#' LD scores of an eigen LD panel
#'
#' For every panel SNP j, \eqn{L2_j = \sum_k r_{jk}^2} over the SNPs k of its LD block, and the
#' stratified scores \eqn{L2\_c_j = \sum_k r_{jk}^2 A_{kc}}, r the correlations of the block's stored eigen factor
#' (`block<N>.eigen.bin`, truncated at `thresh`). The sum runs over the panel's LD block, not over a 1 cM window,
#' and \eqn{r^2} is not corrected for the reference sample size. All numerical work is in `src/ld_score.cpp`
#' (a block costs \eqn{m k^2} per annotation, k the eigen pairs kept, never the m x m LD matrix).
#'
#' With `snps`, k runs over the given SNPs only and only their rows are returned, so the LD scores and the
#' number of SNPs M that the LD score regression takes from its input refer to the same SNPs. Rows with an `r2`
#' column below 1 (SNPs filled by [sumstat_impute()]) are not given SNPs.
#' @param ld_dir the LD panel folder.
#' @param ref panel SNP table in place of `snpinfo.parquet` of `ld_dir` (data.frame or parquet path).
#' @param annot annotations (required): CHR, the position column of `build` (BP37 or BP38) and one numeric
#'   column per annotation; matched to the panel by chromosome and position, panel SNPs without a row get 0. The
#'   all-ones LD score L2 is always computed as well (it enters the regression weights).
#' @param snps NULL (default: every panel SNP) or the summary statistics the LD scores are for (output of
#'   [summarystat_tidy()] or [sumstat_impute()], data.table or file; Block and Index, and `r2` when present).
#' @param build 37 or 38: the position column `annot` is matched on.
#' @param thresh eigenvalue variance kept; <= 0 (default) keeps every stored pair.
#' @param chunk output rows per dense product (memory per thread: chunk x k doubles, plus the block's factor).
#' @param threads OpenMP threads over blocks; 0 takes them from [plink4r_resources()].
#' @return data.table by Block and Index: CHR, BP37, BP38, A1, A2, Block, Index, L2, and one `L2_<name>` and one
#'   `A_<name>` (the annotation value) per annotation column.
#' @export
ld_score <- function(ld_dir, annot, ref = NULL, snps = NULL, build = 38, thresh = 0, chunk = 1024L, threads = 0L) {
  threads <- .threads(threads)
  r <- .with_bp(.ld_panel(ref, ld_dir))
  data.table::setorder(r, Block, Index)
  m <- r[, .(m = .N, ok = identical(Index, seq.int(0L, .N - 1L))), by = Block]
  if (!all(m$ok)) stop("the panel table must hold every eigen row of each block (Index 0..m-1)", call. = FALSE)
  a <- .read_table(annot)
  bp <- paste0("BP", build)
  if (!all(c("CHR", bp) %in% names(a))) stop("annot needs CHR and ", bp, call. = FALSE)
  an <- setdiff(names(a), c("CHR", "BP37", "BP38", "SNP", "A1", "A2"))
  if (!length(an)) stop("annot has no annotation column", call. = FALSE)
  a <- a[, c("CHR", bp, an), with = FALSE]
  data.table::setnames(a, bp, "BP")
  a[, `:=`(CHR = .chr_code(CHR), BP = as.integer(BP))]
  a <- unique(a, by = c("CHR", "BP"))
  i <- a[data.table::data.table(CHR = r$CHR, BP = r[[bp]]), on = .(CHR, BP), which = TRUE]
  message(sum(!is.na(i)), " of ", nrow(r), " panel SNPs have an annotation row (the others get 0)")
  Aa <- as.matrix(a[i, an, with = FALSE])
  storage.mode(Aa) <- "double"
  Aa[is.na(Aa)] <- 0
  A <- cbind(1, Aa)
  keep <- rep(TRUE, nrow(r))
  if (!is.null(snps)) {
    s <- .read_table(snps)
    if (!all(c("Block", "Index") %in% names(s))) stop("snps needs Block and Index", call. = FALSE)
    if ("r2" %in% names(s)) s <- s[r2 >= 1]
    j <- r[data.table::data.table(Block = as.integer(s$Block), Index = as.integer(s$Index)), on = .(Block, Index),
           which = TRUE]
    if (anyNA(j)) stop(sum(is.na(j)), " rows of snps are not panel SNPs (Block, Index)", call. = FALSE)
    keep <- seq_len(nrow(r)) %in% j
    A[!keep, ] <- 0
    message(sum(keep), " of ", nrow(r), " panel SNPs given: LD scores over them only")
  }
  t0 <- Sys.time()
  L <- ld_score_cpp(.ld_block_file(ld_dir, m$Block), c(0L, cumsum(m$m)), A, thresh, as.integer(chunk),
                    as.integer(threads))
  message("LD scores of ", nrow(r), " SNPs in ", nrow(m), " blocks: ",
          format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s")
  out <- r[, .(CHR, BP37, BP38, A1, A2, Block, Index)]
  out[, L2 := L[, 1L]]
  for (k in seq_along(an)) {
    data.table::set(out, j = paste0("L2_", an[k]), value = L[, k + 1L])
    data.table::set(out, j = paste0("A_", an[k]), value = A[, k + 1L])
  }
  out[keep]
}
