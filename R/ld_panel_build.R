# Build a blockwise eigen LD panel (the folder of R/ld_panel.R: snpinfo.parquet, ldminfo.bed, block<N>.eigen.bin)
# from PLINK genotypes. Offline tool: the columns of a block are opened as a geno object (geno_open(), 2-bit
# store), its missing calls are filled and the correlations come from geno_cor() (threads);
# eigendecomposition by base R eigen() (R's BLAS / LAPACK, as the user decided for panel building),
# blocks one after another.

#' Build an eigen LD panel from PLINK BED / PGEN genotypes
#'
#' The output folder is a panel as read by the package (`ld_dir`):
#' `snpinfo.parquet` (CHR, BP37, BP38, A1, A2, A1Freq, Block, Index),
#' `ldminfo.bed` (a header line, then row N = block N: chr (chr1..chr22,
#' chrX), start, end, NumSNP) and `block<N>.eigen.bin`.
#'
#' 1. Rows of `ref` (CHR, BP37 and/or BP38, A1, A2; optional A1Freq =
#'    frequency of A1) on chromosomes 1 to 22 and X are matched to the `.pvar` /
#'    `.bim` of `geno` by chromosome, the `build` position and alleles: a file
#'    whose counted allele (ALT for PGEN, bim column 5 for BED) is `A1` keeps
#'    its sign, one whose counted allele is `A2` is flipped (not one ALT
#'    allele of several at a position, which is matched by its own `A1` only:
#'    see [geno_file()]). The LD of the panel counts `A1`.
#' 2. A variant is placed in the block of `blocks` with start <= position < end. Blocks
#'    without variants are dropped; a block with fewer than `min_snp` variants
#'    is merged into its neighbour on the same chromosome with fewer variants,
#'    repeatedly, before anything is computed. Blocks are numbered 1..B.
#' 3. The plan (one row per block) goes to `out_dir/build/plan.tsv`;
#'    `dry_run = TRUE` stops here.
#' 4. Block by block: the block's variants are opened with [geno_open()], their
#'    missing calls filled as by [genotype_fill()] (with A1Freq the anchored
#'    value, so the filled column has mean 2 A1Freq; otherwise the observed
#'    mean; sd is the sample sd after filling), LD by [geno_cor()], signs
#'    flipped to count A1, variants with sd 0 in the sample dropped, then
#'    [write_eigen_bin()]. A finished block is skipped on a rerun unless
#'    `overwrite`.
#' 5. The SNP table: the matched rows of `ref` with `Block` and `Index`
#'    (0-based eigen row in the block), sorted by Block and Index; CHR is the
#'    integer chromosome code, a build `ref` lacks is `NA`, and A1Freq is kept
#'    from `ref` or, when absent, is the sample frequency of A1.
#'
#' @param geno PLINK prefix(es) (or `.pgen` / `.bed` paths); each file's
#'   chromosomes are read from its own `.pvar` (uncompressed) or `.bim`.
#' @param ref Parquet path or data.frame of candidate variants (see 1.).
#' @param out_dir Output folder (created if needed).
#' @param blocks Block file on the same build as the genotypes: a header line
#'   (names not checked), then whitespace-separated columns chr, start, end
#'   first on every line; further columns are not read, so the `ldminfo.bed`
#'   of a panel is such a file. `NULL` uses the bundled file for `build`
#'   ([ld_blocks_default()]).
#' @param build 38 or 37: the build of the genotypes (picks BP38 or BP37).
#' @param sample IIDs to use (`NULL` = all samples of each file).
#' @param thresh Eigenvalue truncation (proportion of the positive sum).
#' @param min_snp Blocks with fewer matched variants are merged (see 2.).
#' @param threads OpenMP threads for the correlations; 0 takes them from
#'   [plink4r_resources()]. The eigendecomposition uses R's BLAS threads.
#' @param dry_run Only write and return the plan.
#' @param overwrite Recompute blocks whose output already exists.
#' @return Invisibly, list(plan, log, panel) (only `plan` when `dry_run`).
#' @export
ld_panel_build <- function(geno, ref, out_dir, blocks = NULL, build = 38, sample = NULL,
                           thresh = 0.999, min_snp = 100L, threads = 4L, dry_run = FALSE, overwrite = FALSE) {
  threads <- .threads(threads)
  build <- as.integer(build)
  if (!build %in% c(37L, 38L)) stop("build must be 37 or 38.", call. = FALSE)
  if (is.null(blocks)) blocks <- ld_blocks_default(build)
  wd <- file.path(out_dir, "build")
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)

  r <- .ldb_ref(ref, build)
  gv <- .ldb_geno_variants(geno)
  m <- .ldb_match(r$key, gv)
  bl <- .ldb_read_blocks(blocks)
  m <- .ldb_assign(m, bl)
  plan <- .ldb_plan(m, bl, as.integer(min_snp))
  m <- plan$variants
  plan <- plan$blocks
  n_s <- vapply(gv$files, function(f) if (is.null(sample)) nrow(f$file$samples) else length(sample), numeric(1))
  plan[, est_gb := round((n_s[file] * NumSNP / 4 + 5 * 8 * as.numeric(NumSNP)^2) / 1e9, 2)]
  plan[, prefix := vapply(gv$files, `[[`, "", "prefix")[file]]
  data.table::fwrite(plan, file.path(wd, "plan.tsv"), sep = "\t")
  message(sprintf("%d of %d variants matched (%d flipped), %d blocks; largest block %d variants (about %.1f GB)",
                  nrow(m), nrow(r$key), sum(m$sgn < 0), nrow(plan), max(plan$NumSNP), max(plan$est_gb)))
  if (dry_run) return(invisible(list(plan = plan)))

  logf <- file.path(wd, "log.tsv")
  has_af <- "AF" %in% names(m)
  for (b in plan$Block) {
    eig <- .ld_block_file(out_dir, b)
    rf <- file.path(wd, paste0("block", b, ".rows.tsv"))
    if (!overwrite && .ldb_done(eig, rf)) next
    v <- m[Block == b]
    fi <- plan$file[plan$Block == b]
    t0 <- proc.time()[["elapsed"]]
    af <- if (has_af) ifelse(v$sgn > 0, v$AF, 1 - v$AF) else NULL
    X <- geno_open(file = gv$files[[fi]]$file, snp_vec = v$row, sample_vec = sample, scale = FALSE, threads = threads)
    # missing calls: the anchored value where af is finite, else the observed mean (in place: X is local)
    X <- .geno_refill(X, af = af, inplace = TRUE, impute = if (has_af) "af" else "mean", threads = threads)
    t_read <- proc.time()[["elapsed"]]
    keep <- X$sd > 0
    if (!any(keep)) stop("Block ", b, ": no variant varies in the sample.", call. = FALSE)
    R <- geno_cor(X, threads = threads)
    t_cor <- proc.time()[["elapsed"]]
    R <- R[keep, keep, drop = FALSE]
    s <- v$sgn[keep]
    if (any(s < 0)) {
      R[s < 0, ] <- -R[s < 0, ]
      R[, s < 0] <- -R[, s < 0]
    }
    diag(R) <- 1
    dimnames(R) <- NULL
    t1 <- proc.time()[["elapsed"]]
    e <- write_eigen_bin(R, paste0(eig, ".tmp"), thresh)
    rm(R)
    file.rename(paste0(eig, ".tmp"), eig)
    t2 <- proc.time()[["elapsed"]]
    p_alt <- ifelse(s > 0, X$center[keep] / 2, 1 - X$center[keep] / 2)
    data.table::fwrite(data.table::data.table(rid = v$rid[keep], af_sample = p_alt), paste0(rf, ".tmp"), sep = "\t")
    file.rename(paste0(rf, ".tmp"), rf)
    lg <- data.table::data.table(Block = b, NumSNP = nrow(v), kept = sum(keep), k = e$k,
                                 n = X$n, read_s = t_read - t0, cor_s = t_cor - t_read,
                                 eigen_s = t2 - t1, total_s = t2 - t0)
    data.table::fwrite(lg, logf, sep = "\t", append = file.exists(logf))
    message(sprintf("block %d/%d: m = %d (%d dropped), k = %d, %.1f s", b, nrow(plan), sum(keep),
                    sum(!keep), e$k, t2 - t0))
    rm(X)
    gc(FALSE)
  }

  rows <- data.table::rbindlist(lapply(plan$Block, function(b) {
    x <- data.table::fread(file.path(wd, paste0("block", b, ".rows.tsv")))
    x[, `:=`(Block = b, Index = seq_len(.N) - 1L)]
  }))
  panel <- .with_bp(r$table[rows$rid])
  panel[, `:=`(Block = rows$Block, Index = rows$Index)]
  if (!"A1Freq" %in% names(panel)) panel[, A1Freq := rows$af_sample]
  data.table::setcolorder(panel, .panel_cols)
  arrow::write_parquet(panel, file.path(out_dir, .panel_snpinfo))
  data.table::fwrite(plan[, .(chr = paste0("chr", ifelse(Chrom == 23L, "X", Chrom)), start = StartBP, end = EndBP, NumSNP)],
                     file.path(out_dir, .panel_ldminfo), sep = "\t")
  invisible(list(plan = plan, log = if (file.exists(logf)) data.table::fread(logf) else NULL, panel = panel))
}

#' Eigen-decompose an LD matrix, truncate at `thresh` of the positive eigenvalue sum, write plink4r eigen.bin
#'
#' Layout: int32 m; int32 k; float64 sumPos; float32 cutoff; float64 lam[k]; float32 scale[k]; int16 U[m*k]
#' (column-major, u = q * scale, scale_j = max|U_j| / 32767).
#' @return list(m, k, sumpos) invisibly
#' @export
write_eigen_bin <- function(R, path, thresh = 0.999) {
  e <- eigen(R, symmetric = TRUE)
  lam <- e$values; lp <- pmax(lam, 0); sumpos <- sum(lp)
  k <- min(which(cumsum(lp) / sumpos >= thresh)[1], length(lam))
  U <- e$vectors[, seq_len(k), drop = FALSE]
  sc <- apply(abs(U), 2, max) / 32767
  q <- round(sweep(U, 2, sc, "/"))
  q[q > 32767] <- 32767; q[q < -32767] <- -32767
  con <- file(path, "wb"); on.exit(close(con))
  writeBin(as.integer(c(nrow(R), k)), con, size = 4L, endian = "little")
  writeBin(sumpos, con, size = 8L, endian = "little")
  writeBin(thresh, con, size = 4L, endian = "little")
  writeBin(lam[seq_len(k)], con, size = 8L, endian = "little")
  writeBin(sc, con, size = 4L, endian = "little")
  writeBin(as.integer(q), con, size = 2L, endian = "little")
  invisible(list(m = nrow(R), k = k, sumpos = sumpos))
}

#' Bundled LD block files
#'
#' GRCh38: the blocks of the Pan-UKB EUR panel, its `ldminfo.bed` as is
#' (1,341 blocks on chr1-22, contiguous from 0, from pyrho with gnomAD / deCODE
#' gap filling; 30 blocks on chrX, with gaps between them: X variants in a gap
#' are in no block). Columns chr, start, end, source, n_gnomAD_snv.
#' GRCh37: LDetect EUR blocks (Berisa and Pickrell 2016; 1,702 blocks,
#' chr1-22), columns chr, start, end.
#' @param build 38 or 37.
#' @return Path to the block file.
#' @export
ld_blocks_default <- function(build = 38) {
  f <- switch(as.character(build),
              "38" = "ldminfo_GRCh38_PanUKB_EUR.bed",
              "37" = "blocks_GRCh37_LDetect_EUR.txt",
              stop("build must be 37 or 38.", call. = FALSE))
  system.file("extdata", f, package = "plink4r", mustWork = TRUE)
}

# ref table (parquet or data.frame): its panel columns (chromosome code, alleles in upper case), and the key CHR /
# BP / A1 / A2 (/ AF) of its rows on the panel chromosomes (1 to 22, X) with a position on the build, with row id rid.
.ldb_ref <- function(x, build) {
  d <- if (is.data.frame(x)) data.table::as.data.table(x) else data.table::as.data.table(arrow::read_parquet(x))
  bp <- paste0("BP", build)
  miss <- setdiff(c("CHR", bp, "A1", "A2"), names(d))
  if (length(miss)) stop("ref lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  d <- d[, intersect(setdiff(.panel_cols, c("Block", "Index")), names(d)), with = FALSE]
  d[, `:=`(CHR = .chr_code(CHR), A1 = .by_unique(as.character(A1), toupper), A2 = .by_unique(as.character(A2), toupper))]
  for (b in intersect(c("BP37", "BP38"), names(d))) data.table::set(d, j = b, value = as.integer(d[[b]]))
  key <- data.table::data.table(rid = seq_len(nrow(d)), CHR = d$CHR, BP = d[[bp]], A1 = d$A1, A2 = d$A2)
  if ("A1Freq" %in% names(d)) key[, AF := as.numeric(d$A1Freq)]
  list(table = d, key = key[CHR %in% .panel_chr & !is.na(BP)])
}

# Variants of every genotype file: file, 1-based row, chromosome code, position, counted / other allele of its
# variants on the panel chromosomes (1 to 22, X; a .pvar's PAR1 / PAR2 / XY rows are X), and several (.several()).
# files[[f]] = list(prefix as given, file = its geno_file()), from which every block of the file is opened.
.ldb_geno_variants <- function(geno) {
  files <- list()
  vars <- list()
  for (f in seq_along(geno)) {
    pre <- sub("\\.(pgen|bed)$", "", geno[f])
    gf <- if (file.exists(paste0(pre, ".pgen"))) geno_file(pgenfile = pre) else
      if (file.exists(paste0(pre, ".bed"))) geno_file(bedfile = pre) else
        stop("No .pgen or .bed for ", geno[f], call. = FALSE)
    files[[f]] <- list(prefix = pre, file = gf)
    v <- gf$variants
    v <- data.table::data.table(file = f, row = v$row, CHR = v$CHR, BP = v$BP, cnt = v$A1, oth = v$A2, several = .several(v))
    vars[[f]] <- v[CHR %in% .panel_chr]
  }
  vars <- data.table::rbindlist(vars)
  fc <- unique(vars[, .(file, CHR)])
  if (anyDuplicated(fc$CHR)) {
    stop("Chromosome ", fc$CHR[anyDuplicated(fc$CHR)], " appears in more than one genotype file.", call. = FALSE)
  }
  list(files = files, variants = vars)
}

.ldb_match <- function(ref, gv) {
  v <- gv$variants
  same <- v[ref, on = .(CHR, BP, cnt = A1, oth = A2), nomatch = NULL]
  same[, sgn := 1L]
  v <- v[!v$several]  # one ALT allele of several is matched in its own orientation only
  data.table::setnames(same, c("cnt", "oth"), c("A1", "A2"))
  swap <- v[ref, on = .(CHR, BP, cnt = A2, oth = A1), nomatch = NULL]
  swap[, sgn := -1L]
  data.table::setnames(swap, c("cnt", "oth"), c("A2", "A1"))
  m <- data.table::rbindlist(list(same, swap), use.names = TRUE)[, several := NULL]
  data.table::setorder(m, CHR, BP, A2, A1, file, row)
  dup <- duplicated(m, by = c("CHR", "BP", "A2", "A1"))
  if (any(dup)) {
    warning(sum(dup), " ref variants match more than one genotype row; the first is used.", call. = FALSE)
    m <- m[!dup]
  }
  if (!nrow(m)) stop("No ref variant matches the genotype files.", call. = FALSE)
  m
}

# Block file: first line is a header (names not checked), then on every line the whitespace-separated columns chr
# start end first (further columns are not read: the ldminfo.bed of a panel is such a file). The blocks of
# the panel chromosomes (1 to 22, X), with the chromosome code.
.ldb_read_blocks <- function(path) {
  x <- trimws(readLines(path))
  x <- x[nzchar(x)]
  f <- strsplit(x, "[[:space:]]+")
  nf <- lengths(f)
  if (length(f) < 2L || any(nf < 3L)) {
    bad <- which(nf < 3L)
    stop("Block file must be a header line plus lines that begin with the columns chr start end",
         if (length(bad)) paste0(" (line ", bad[1L], " has ", nf[bad[1L]], ")"), ": ", path, call. = FALSE)
  }
  body <- vapply(f[-1L], function(a) a[1:3], character(3))
  bl <- data.table::data.table(src = seq_len(ncol(body)), CHR = .chr_code(body[1L, ]),
                               start = suppressWarnings(as.numeric(body[2L, ])),
                               end = suppressWarnings(as.numeric(body[3L, ])))
  if (anyNA(bl$start) || anyNA(bl$end)) stop("start and end must be numbers: ", path, call. = FALSE)
  bl <- bl[CHR %in% .panel_chr]
  data.table::setorder(bl, CHR, start)
  bl
}

# Each variant to the block with start <= BP < end on its chromosome; others dropped.
.ldb_assign <- function(m, bl) {
  m[, src := NA_integer_]
  for (ch in intersect(unique(m$CHR), unique(bl$CHR))) {
    b <- bl[CHR == ch]
    ii <- which(m$CHR == ch)
    k <- findInterval(m$BP[ii], b$start)
    ok <- k > 0L
    ok[ok] <- m$BP[ii][ok] < b$end[k[ok]]
    data.table::set(m, ii[ok], "src", b$src[k[ok]])
  }
  out <- m[!is.na(src)]
  if (nrow(out) < nrow(m)) message(nrow(m) - nrow(out), " matched variants fall outside every block and are dropped.")
  out
}

# Drop empty blocks, merge blocks below min_snp into the smaller neighbour, number 1..B.
.ldb_plan <- function(m, bl, min_snp) {
  cnt <- m[, .N, by = src]
  b <- bl[cnt, on = "src", nomatch = NULL]
  data.table::setorder(b, CHR, start)
  b[, members := as.character(src)]
  res <- list()
  for (ch in unique(b$CHR)) {
    x <- b[CHR == ch]
    while (nrow(x) > 1L && any(x$N < min_snp)) {
      i <- which.min(x$N)
      nb <- c(i - 1L, i + 1L)
      nb <- nb[nb >= 1L & nb <= nrow(x)]
      j <- nb[which.min(x$N[nb])]
      lo <- min(i, j); hi <- max(i, j)
      x$end[lo] <- x$end[hi]
      x$N[lo] <- x$N[lo] + x$N[hi]
      x$members[lo] <- paste(x$members[lo], x$members[hi], sep = ",")
      m[src == x$src[hi], src := x$src[lo]]
      x <- x[-hi]
    }
    res[[length(res) + 1L]] <- x
  }
  b <- data.table::rbindlist(res)
  b <- b[order(CHR, start)]
  b[, Block := seq_len(.N)]
  m <- m[b[, .(src, Block)], on = "src"]
  m <- m[order(Block, BP, A2, A1)]
  fl <- m[, .(file = file[1L], nfile = data.table::uniqueN(file), n_flip = sum(sgn < 0)), by = Block]
  plan <- b[fl, on = "Block"][, .(Block, Chrom = CHR, StartBP = start, EndBP = end, NumSNP = N,
                                  n_flip, src_blocks = members, file)]
  list(blocks = plan, variants = m)
}

.ldb_done <- function(eig, vf) {
  if (!file.exists(eig) || !file.exists(vf)) return(FALSE)
  con <- file(eig, "rb")
  on.exit(close(con))
  m <- readBin(con, integer(), n = 1L, size = 4L, endian = "little")
  length(m) == 1L && m == nrow(data.table::fread(vf, select = 1L))
}
