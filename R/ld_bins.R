#' LD-bin list of an eigen LD panel
#'
#' Groups the SNPs of every LD block into bins of high LD, for the LD-bin suppression of [viprs_annot()] and
#' [sbayesrc()]: inside a block, the SNP with the most remaining neighbours (r2 > `r2`) opens a bin with all its
#' remaining neighbours, until every SNP is in a bin (a SNP without neighbours is a bin of its own). r is the
#' correlation of the block's stored eigen factor (truncated at `thresh`); the m x m LD matrix is never formed. All
#' numerical work is in `src/ld_bins.cpp`.
#' @param ld_dir the LD panel folder.
#' @param r2 r2 threshold of a bin.
#' @param ref panel SNP table in place of `snpinfo.parquet` of `ld_dir` (data.frame or parquet path); it must hold
#'   every eigen row of each of its blocks.
#' @param thresh eigenvalue variance kept; <= 0 (default) keeps every stored pair.
#' @param output the bin file: one int16 per SNP of the panel table in the order Block, Index, its bin id inside its
#'   block (the file read by [viprs_annot()] and [sbayesrc()]); default `ld_dir`/ldbin<r2>.bin. NULL writes nothing.
#' @param chunk rows compared per dense product (memory per thread: chunk x m floats, plus the block's factor).
#' @param threads OpenMP threads over blocks; 0 takes them from [plink4r_resources()].
#' @return data.table by Block and Index: Block, Index, BIN (invisibly when written).
#' @export
ld_bins <- function(ld_dir, r2 = 0.98, ref = NULL, thresh = 0, output = file.path(ld_dir, paste0("ldbin", r2, ".bin")),
                    chunk = 1024L, threads = 0L) {
  threads <- .threads(threads)
  r <- .ld_panel(ref, ld_dir)
  data.table::setorder(r, Block, Index)
  m <- r[, .(m = .N, ok = identical(Index, seq.int(0L, .N - 1L))), by = Block]
  if (!all(m$ok)) stop("the panel table must hold every eigen row of each block (Index 0..m-1)", call. = FALSE)
  t0 <- Sys.time()
  b <- ld_bins_cpp(.ld_block_file(ld_dir, m$Block), as.integer(m$m), thresh, r2, as.integer(chunk), as.integer(threads))
  if (max(b) > 32767L) stop("more than 32768 bins in a block: the int16 bin file cannot hold them", call. = FALSE)
  out <- r[, .(Block, Index)][, BIN := b]
  sz <- out[, .N, by = .(Block, BIN)]
  message(sprintf("LD bins (r2 > %g) of %d SNPs in %d blocks: %d bins, %d with more than one SNP (%d SNPs), %.1f s",
                  r2, nrow(out), nrow(m), nrow(sz), sum(sz$N > 1), sum(sz$N[sz$N > 1]),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (is.null(output)) return(out[])
  con <- file(output, "wb")
  on.exit(close(con))
  writeBin(as.integer(b), con, size = 2L, endian = "little")
  message("written: ", output)
  invisible(out[])
}
