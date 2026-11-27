# The LD panel in one place: what a panel folder holds, its SNP table, the file of a block, and the alignment of
# the panel to a genotype file.
#
# A panel is a folder with
#   snpinfo.parquet     one row per panel SNP: CHR (1..22, 23 = X), BP37, BP38 (NA where a build has no position),
#                       A1 (the allele the LD counts), A2, A1Freq (frequency of A1), Block (block number, from 1),
#                       Index (0-based eigen row inside the block)
#   ldminfo.bed         one header line, then row N = block N: chr, start, end (a block holds the SNPs with
#                       start <= position < end), then any further columns
#   block<N>.eigen.bin  the eigen factor of the LD of block N (src/eigen_bin.h)
# A user's file need not share the panel's allele order, IDs or build: .ld_align() matches every panel SNP to
# the file by chromosome code (R/chr.R), position and the allele pair, and records the orientation as a sign
# (+1: the file counts the panel's A1, -1: it counts A2). Panel SNPs absent from the file are not carried into an
# imputation at all.

.panel_cols <- c("CHR", "BP37", "BP38", "A1", "A2", "A1Freq", "Block", "Index")
.panel_snpinfo <- "snpinfo.parquet"
.panel_ldminfo <- "ldminfo.bed"

# The eigen.bin of blocks of a panel folder (the one place that names these files; "~" expanded for the C++ readers).
.ld_block_file <- function(ld_dir, block) file.path(path.expand(ld_dir), paste0("block", block, ".eigen.bin"))

# Panel SNP table, optionally only some blocks. ref: a data.frame or parquet path in place of the folder's table;
# NULL (or a folder) reads snpinfo.parquet of ld_dir. With `blocks` only those blocks are scanned (the table is
# sorted by block, so whole row groups are skipped).
.ld_panel <- function(ref = NULL, ld_dir = NULL, blocks = NULL) {
  if (is.character(ref) && length(ref) == 1L && dir.exists(ref)) { ld_dir <- ref; ref <- NULL }
  if (is.null(ref)) {
    if (is.null(ld_dir)) stop("give ld_dir (or ref)", call. = FALSE)
    ref <- file.path(ld_dir, .panel_snpinfo)
    if (!file.exists(ref)) stop("no ", .panel_snpinfo, " in ", ld_dir, call. = FALSE)
  }
  if (is.data.frame(ref)) {
    r <- data.table::as.data.table(ref)
    r <- r[, intersect(.panel_cols, names(r)), with = FALSE]
  } else {
    ds <- arrow::open_dataset(ref)
    cols <- intersect(.panel_cols, names(ds))
    r <- if (is.null(blocks)) data.table::as.data.table(arrow::read_parquet(ref, col_select = cols)) else {
      inb <- arrow::Expression$create("is_in", arrow::Expression$field_ref("Block"),
                                      options = list(value_set = arrow::Array$create(as.integer(blocks)), skip_nulls = TRUE))
      data.table::as.data.table(as.data.frame(ds$NewScan()$Filter(inb)$Project(cols)$Finish()$ToTable()))
    }
  }
  need <- c("CHR", "A1", "A2", "A1Freq", "Block", "Index")
  if (length(setdiff(need, names(r)))) stop("the panel table lacks ", paste(setdiff(need, names(r)), collapse = ", "), call. = FALSE)
  builds <- intersect(c("BP37", "BP38"), names(r))
  if (!length(builds)) stop("the panel table has neither BP37 nor BP38", call. = FALSE)
  code <- .chr_code(r$CHR)
  if (anyNA(code)) stop("the panel table has a CHR without a chromosome code: ",
                        paste(utils::head(unique(r$CHR[is.na(code)]), 5L), collapse = ", "), call. = FALSE)
  r[, `:=`(CHR = code, A1 = .by_unique(as.character(A1), toupper), A2 = .by_unique(as.character(A2), toupper),
           A1Freq = as.numeric(A1Freq), Block = as.integer(Block), Index = as.integer(Index))]
  for (b in builds) data.table::set(r, j = b, value = as.integer(r[[b]]))
  if (!is.null(blocks)) r <- r[which(r$Block %in% blocks)]
  if (anyNA(r$Block) || anyNA(r$Index) || anyDuplicated(r, by = c("Block", "Index")))
    stop("Block and Index of the panel table must be complete and name each SNP once", call. = FALSE)
  r
}

# Adds vi (0-based variant of the file) and sign (+1 / -1) to the panel rows found in the file; the others are
# dropped. variants: CHR (chromosome code), BP, A1 (the counted allele), A2, row, as geno_file()$variants. A panel
# SNP is matched by chromosome, position and allele pair, in either orientation; a variant that shares its CHR, BP
# and A2 with another (one ALT allele of several, .several()) only in its own orientation. The position column is
# the panel's only build column (BP37 or BP38: the file must be on that build); with both, every panel row on the
# file's chromosomes is matched on each build and the build that places more SNPs is taken. A row without a
# position on that build is matched to nothing. A variant of the file is used once. No panel SNP in the file is an
# error.
.ld_align <- function(panel, variants) {
  builds <- intersect(c("BP37", "BP38"), names(panel))
  pv <- data.table::data.table(CHR = variants$CHR, POS = as.integer(variants$BP), REF = variants$A2,
                               ALT = variants$A1, row = variants$row, several = .several(variants))
  pv <- pv[!is.na(CHR) & !is.na(POS)]
  q <- data.table::data.table(i = seq_len(nrow(panel)), CHR = panel$CHR, A1 = panel$A1, A2 = panel$A2)
  q <- q[CHR %in% unique(pv$CHR)]
  hits <- function(b) {
    qb <- data.table::copy(q)[, POS := panel[[b]][i]][!is.na(POS)]
    same <- pv[qb, on = .(CHR, POS, REF = A2, ALT = A1), nomatch = NULL, .(i = i.i, row = x.row, s = 1L)]
    swap <- pv[!several][qb, on = .(CHR, POS, REF = A1, ALT = A2), nomatch = NULL, .(i = i.i, row = x.row, s = -1L)]
    h <- rbind(same, swap)
    h[!duplicated(i) & !duplicated(row)]
  }
  best <- hits(builds[1L])
  if (length(builds) == 2L) {
    other <- hits(builds[2L])
    if (nrow(other) > nrow(best)) best <- other
  }
  if (!nrow(best)) stop("no panel SNP found in the genotype file (by chromosome, ", paste(builds, collapse = " or "),
                        " and alleles): check the build", call. = FALSE)
  data.table::setorder(best, i)
  p <- panel[best$i]
  p[, `:=`(vi = best$row - 1L, sign = best$s)]
  p
}

# Working rows of one block: the aligned panel rows of block `blk` (eigen rows lo..hi when given), by eigen row.
.ld_rows <- function(aligned, blk, lo = NULL, hi = NULL) {
  a <- aligned[which(aligned$Block == blk)]
  if (!is.null(lo)) a <- a[which(a$Index >= lo & a$Index <= hi)]
  data.table::setorder(a, Index)
  a
}

# The dosage grid of the LD input: the step count round(D / step), |D| <= 0.5, is int8.
.check_dosage_step <- function(step) {
  if (!is.numeric(step) || length(step) != 1L || !is.finite(step) || step < 0.5 / 127)
    stop("dosage_step must be a single number >= 0.5 / 127", call. = FALSE)
  invisible(step)
}

# The ridge delta of every LD imputation must be > 0 (delta = 0 has no stable C = I - delta (R + delta I)^-1).
.check_delta <- function(delta) {
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0)
    stop("delta must be a single number > 0", call. = FALSE)
  invisible(delta)
}

# BP37 / BP38 present in every output table; a build the panel lacks is NA.
.with_bp <- function(r) {
  for (b in c("BP37", "BP38")) if (!b %in% names(r)) data.table::set(r, j = b, value = NA_integer_)
  r
}
