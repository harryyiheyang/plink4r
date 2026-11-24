#' Cross-ancestry LD scores of two eigen LD panels
#'
#' The SNPs of panel 1 and panel 2 are paired by chromosome, the `build` position and the allele pair (either
#' order; `sign` is -1 where panel 2 counts the other allele). For a paired SNP j,
#' \eqn{L2\_cross_j = \sum_k r^{(1)}_{jk} r^{(2)}_{jk}} over the paired SNPs k in the same block of panel 1 and the
#' same block of panel 2, \eqn{r^{(2)}} turned to the alleles of panel 1. All numerical work is in
#' `src/ld_score_cross.cpp`.
#' @param ld_dir1,ld_dir2 the two panel folders.
#' @param ref1,ref2 their SNP tables in place of `snpinfo.parquet`.
#' @param build 37 or 38: the position both panels are paired on.
#' @param thresh1,thresh2 eigenvalue variance kept in each panel; <= 0 keeps every stored pair.
#' @inheritParams ld_score
#' @return data.table by Block1 and Index1: CHR, BP37, BP38, A1, A2 (of panel 1), Block1, Index1, Block2, Index2,
#'   sign, L2_cross.
#' @export
ld_score_cross <- function(ld_dir1, ld_dir2, ref1 = NULL, ref2 = NULL, build = 38, thresh1 = 0, thresh2 = 0,
                           chunk = 1024L, threads = 0L) {
  threads <- .threads(threads)
  bp <- paste0("BP", build)
  r1 <- .with_bp(.ld_panel(ref1, ld_dir1))
  r2 <- .with_bp(.ld_panel(ref2, ld_dir2))
  k2 <- r2[!is.na(get(bp)), .(CHR, BP = get(bp), A1, A2, Block2 = Block, Index2 = Index)]
  k1 <- r1[!is.na(get(bp)), .(CHR, BP37, BP38, BP = get(bp), A1, A2, Block1 = Block, Index1 = Index)]
  same <- k2[k1, on = .(CHR, BP, A1, A2), nomatch = NULL][, sign := 1]
  swap <- k2[k1, on = .(CHR, BP, A1 = A2, A2 = A1), nomatch = NULL][, sign := -1]
  x <- rbind(same, swap, use.names = TRUE)
  x <- x[!duplicated(x[, .(Block1, Index1)]) & !duplicated(x[, .(Block2, Index2)])]
  if (!nrow(x)) stop("no SNP of panel 1 found in panel 2 (by chromosome, ", bp, " and alleles)", call. = FALSE)
  data.table::setorder(x, Block1, Block2, Index1)
  grp <- x[, .N, by = .(Block1, Block2)]
  message(nrow(x), " SNPs paired (", sum(x$sign < 0), " with A1 / A2 swapped) in ", nrow(grp), " block pairs")
  t0 <- Sys.time()
  L <- ld_score_cross_cpp(.ld_block_file(ld_dir1, grp$Block1), .ld_block_file(ld_dir2, grp$Block2),
                          c(0L, cumsum(grp$N)), x$Index1, x$Index2, as.numeric(x$sign), thresh1, thresh2,
                          as.integer(chunk), as.integer(threads))
  message("cross LD scores: ", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)), " s")
  out <- x[, .(CHR, BP37, BP38, A1, A2, Block1, Index1, Block2, Index2, sign)]
  out[, L2_cross := L]
  data.table::setorder(out, Block1, Index1)
  out[]
}
