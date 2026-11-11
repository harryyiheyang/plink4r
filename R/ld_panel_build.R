# Build a blockwise eigen LD panel (one SNP-table parquet + block<N>.eigen.bin) from
# PLINK genotypes. Offline tool: correlations by cor_pgen() / cor_bed() (2-bit
# store, threads), eigendecomposition by base R eigen() (R's BLAS / LAPACK, as
# the user decided for panel building), blocks one after another.

#' Build an eigen LD panel from PLINK BED / PGEN genotypes
#'
#' The output folder is a panel as read by the package (`ld_dir`): one parquet
#' SNP table (ID, chr, bp37 and/or bp38, ref, alt, AF_PanUKB, block, ld_index)
#' and `block<N>.eigen.bin`.
#'
#' 1. Rows of `ref` (chr, bp37 and/or bp38, ref, alt; optional ID and
#'    AF_PanUKB = alt frequency) are matched to the `.pvar` / `.bim` of `geno`
#'    by chromosome, the `build` position and alleles: a file whose counted
#'    allele (ALT for PGEN, bim column 5 for BED) is `alt` keeps its sign, one
#'    whose counted allele is `ref` is flipped; multiallelic `.pvar` rows are
#'    skipped. The panel is in the `ref` / `alt` orientation of `ref`.
#' 2. Variants are placed in the blocks of `blocks` by `[start, stop)`. Blocks
#'    without variants are dropped; a block with fewer than `min_snp` variants
#'    is merged into its neighbour on the same chromosome with fewer variants,
#'    repeatedly, before anything is computed. Blocks are numbered 1..B.
#' 3. The plan (one row per block) goes to `out_dir/build/plan.tsv`;
#'    `dry_run = TRUE` stops here.
#' 4. Block by block: LD by [cor_pgen()] / [cor_bed()] (with AF_PanUKB the
#'    missing calls take the anchored constant, so the filled column has mean
#'    2 AF; otherwise the observed mean; sd is the sample sd after filling),
#'    signs flipped to the alt orientation, variants with sd 0 in the sample
#'    dropped, then [write_eigen_bin()]. A finished block is skipped on a rerun
#'    unless `overwrite`.
#' 5. The SNP table: the matched rows of `ref` with `block` and `ld_index`
#'    (0-based eigen row in the block), sorted by block and ld_index; ID is
#'    kept from `ref` or made as `chr:pos:ref:alt`; AF_PanUKB is kept from
#'    `ref` or, when absent, the sample alt frequency.
#'
#' @param geno PLINK prefix(es) (or `.pgen` / `.bed` paths); each file's
#'   chromosomes are read from its own `.pvar` (uncompressed) or `.bim`.
#' @param ref Parquet path or data.frame of candidate variants (see 1.).
#' @param out_dir Output folder (created if needed).
#' @param blocks Block file on the same build as the genotypes: a header line
#'   (names not checked), then exactly three whitespace-separated columns chr,
#'   start, stop; anything else is an error. `NULL` uses the bundled file for
#'   `build` ([ld_blocks_default()]).
#' @param build 38 or 37: the build of the genotypes (picks bp38 or bp37).
#' @param sample IIDs to use (`NULL` = all samples of each file).
#' @param thresh Eigenvalue truncation (proportion of the positive sum).
#' @param min_snp Blocks with fewer matched variants are merged (see 2.).
#' @param threads OpenMP threads for the correlations; 0 takes them from
#'   [plink4r_resources()]. The eigendecomposition uses R's BLAS threads.
#' @param name File name of the SNP table.
#' @param dry_run Only write and return the plan.
#' @param overwrite Recompute blocks whose output already exists.
#' @return Invisibly, list(plan, log, panel) (only `plan` when `dry_run`).
#' @export
ld_panel_build <- function(geno, ref, out_dir, blocks = NULL, build = 38, sample = NULL,
                           thresh = 0.999, min_snp = 100L, threads = 4L,
                           name = "ld_panel.parquet", dry_run = FALSE, overwrite = FALSE) {
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
  n_s <- vapply(gv$files, function(f) if (is.null(sample)) f$n else length(sample), numeric(1))
  plan[, est_gb := round((n_s[file] * NumSNP / 4 + 5 * 8 * as.numeric(NumSNP)^2) / 1e9, 2)]
  plan[, prefix := vapply(gv$files, `[[`, "", "prefix")[file]]
  data.table::fwrite(plan, file.path(wd, "plan.tsv"), sep = "\t")
  message(sprintf("%d of %d variants matched (%d flipped), %d blocks; largest block %d variants (about %.1f GB)",
                  nrow(m), nrow(r$key), sum(m$sgn < 0), nrow(plan), max(plan$NumSNP), max(plan$est_gb)))
  if (dry_run) return(invisible(list(plan = plan)))

  logf <- file.path(wd, "log.tsv")
  has_af <- "AF" %in% names(m)
  for (b in plan$Block) {
    eig <- file.path(out_dir, paste0("block", b, ".eigen.bin"))
    rf <- file.path(wd, paste0("block", b, ".rows.tsv"))
    if (!overwrite && .ldb_done(eig, rf)) next
    v <- m[Block == b]
    fi <- plan$file[plan$Block == b]
    t0 <- proc.time()[["elapsed"]]
    af <- if (has_af) ifelse(v$sgn > 0, v$AF, 1 - v$AF) else NULL
    cor_fun <- if (gv$files[[fi]]$type == "pgen") cor_pgen else cor_bed
    res <- cor_fun(gv$files[[fi]]$prefix, snps = v$row, sample = sample, af = af,
                   threads = threads, impute = "mean")
    keep <- res$sd > 0
    if (!any(keep)) stop("Block ", b, ": no variant varies in the sample.", call. = FALSE)
    R <- res$R[keep, keep, drop = FALSE]
    res$R <- NULL
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
    p_alt <- ifelse(s > 0, res$mean[keep] / 2, 1 - res$mean[keep] / 2)
    data.table::fwrite(data.table::data.table(rid = v$rid[keep], af_sample = p_alt), paste0(rf, ".tmp"), sep = "\t")
    file.rename(paste0(rf, ".tmp"), rf)
    lg <- data.table::data.table(Block = b, NumSNP = nrow(v), kept = sum(keep), k = e$k,
                                 n = res$n, read_s = res$time[["read"]], cor_s = res$time[["cor"]],
                                 eigen_s = t2 - t1, total_s = t2 - t0)
    data.table::fwrite(lg, logf, sep = "\t", append = file.exists(logf))
    message(sprintf("block %d/%d: m = %d (%d dropped), k = %d, %.1f s", b, nrow(plan), sum(keep),
                    sum(!keep), e$k, t2 - t0))
    rm(res)
    gc(FALSE)
  }

  rows <- data.table::rbindlist(lapply(plan$Block, function(b) {
    x <- data.table::fread(file.path(wd, paste0("block", b, ".rows.tsv")))
    x[, `:=`(block = b, ld_index = seq_len(.N) - 1L)]
  }))
  panel <- r$table[rows$rid]
  panel[, `:=`(block = rows$block, ld_index = rows$ld_index)]
  if (!"AF_PanUKB" %in% names(panel)) panel[, AF_PanUKB := rows$af_sample]
  if (!"ID" %in% names(panel)) panel[, ID := paste(sub("^chr", "", chr, ignore.case = TRUE),
                                                   get(paste0("bp", build)), ref, alt, sep = ":")]
  data.table::setcolorder(panel, intersect(.ldb_cols, names(panel)))
  for (f in list.files(out_dir, pattern = "\\.parquet$", full.names = TRUE)) if (basename(f) != name) file.remove(f)
  arrow::write_parquet(panel, file.path(out_dir, name))
  invisible(list(plan = plan, log = if (file.exists(logf)) data.table::fread(logf) else NULL, panel = panel))
}

#' Bundled LD block files
#'
#' GRCh38: pyrho EUR blocks (1,336 blocks, chr1-22). GRCh37: LDetect EUR
#' blocks (Berisa and Pickrell 2016; 1,702 blocks, chr1-22). Both are a header
#' line plus chr, start, stop.
#' @param build 38 or 37.
#' @return Path to the block file.
#' @export
ld_blocks_default <- function(build = 38) {
  f <- switch(as.character(build),
              "38" = "blocks_GRCh38_pyrho_EUR.txt",
              "37" = "blocks_GRCh37_LDetect_EUR.txt",
              stop("build must be 37 or 38.", call. = FALSE))
  system.file("extdata", f, package = "plink4r", mustWork = TRUE)
}

.ldb_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
.ldb_auto <- as.character(1:22)

.ldb_cols <- c("ID", "chr", "bp37", "bp38", "ref", "alt", "AF_PanUKB", "block", "ld_index")

# ref table (parquet or data.frame): its panel columns, and a key CHR / BP / REF / ALT (/ AF) with row id rid.
.ldb_ref <- function(x, build) {
  d <- if (is.data.frame(x)) data.table::as.data.table(x) else data.table::as.data.table(arrow::read_parquet(x))
  bp <- paste0("bp", build)
  miss <- setdiff(c("chr", bp, "ref", "alt"), names(d))
  if (length(miss)) stop("ref lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  d <- d[, intersect(setdiff(.ldb_cols, c("block", "ld_index")), names(d)), with = FALSE]
  key <- data.table::data.table(rid = seq_len(nrow(d)), CHR = .ldb_chr(d$chr), BP = as.integer(d[[bp]]),
                                REF = toupper(d$ref), ALT = toupper(d$alt))
  if ("AF_PanUKB" %in% names(d)) key[, AF := as.numeric(d$AF_PanUKB)]
  list(table = d, key = key[CHR %in% .ldb_auto & !is.na(BP)])
}

# Variants of every genotype file: file, 1-based row, chromosome, position, counted / other allele.
.ldb_geno_variants <- function(geno) {
  files <- list()
  vars <- list()
  for (f in seq_along(geno)) {
    pre <- sub("\\.(pgen|bed)$", "", geno[f])
    if (file.exists(paste0(pre, ".pgen"))) {
      pv <- .read_pvar5(pre)
      files[[f]] <- list(prefix = pre, type = "pgen", n = length(.psam_ids(pre)))
      v <- data.table::data.table(file = f, row = pv$row, CHR = .ldb_chr(pv$CHROM),
                                  BP = pv$POS, cnt = pv$ALT, oth = pv$REF)
      v <- v[pv$allele_ct == 2L]
    } else if (file.exists(paste0(pre, ".bed"))) {
      bim <- data.table::fread(paste0(pre, ".bim"), header = FALSE, colClasses = "character")
      fam <- data.table::fread(paste0(pre, ".fam"), header = FALSE, select = 1L)
      files[[f]] <- list(prefix = pre, type = "bed", n = nrow(fam))
      v <- data.table::data.table(file = f, row = seq_len(nrow(bim)), CHR = .ldb_chr(bim[[1L]]),
                                  BP = as.integer(bim[[4L]]), cnt = toupper(bim[[5L]]), oth = toupper(bim[[6L]]))
    } else {
      stop("No .pgen or .bed for ", geno[f], call. = FALSE)
    }
    vars[[f]] <- v[CHR %in% .ldb_auto]
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
  same <- v[ref, on = .(CHR, BP, cnt = ALT, oth = REF), nomatch = NULL]
  same[, sgn := 1L]
  data.table::setnames(same, c("cnt", "oth"), c("ALT", "REF"))
  swap <- v[ref, on = .(CHR, BP, cnt = REF, oth = ALT), nomatch = NULL]
  swap[, sgn := -1L]
  data.table::setnames(swap, c("cnt", "oth"), c("REF", "ALT"))
  m <- data.table::rbindlist(list(same, swap), use.names = TRUE)
  data.table::setorder(m, CHR, BP, REF, ALT, file, row)
  dup <- duplicated(m, by = c("CHR", "BP", "REF", "ALT"))
  if (any(dup)) {
    warning(sum(dup), " ref variants match more than one genotype row; the first is used.", call. = FALSE)
    m <- m[!dup]
  }
  if (!nrow(m)) stop("No ref variant matches the genotype files.", call. = FALSE)
  m
}

# Block file: first line is a header (names not checked), then exactly three
# whitespace-separated columns chr start stop; anything else is an error.
.ldb_read_blocks <- function(path) {
  x <- trimws(readLines(path))
  x <- x[nzchar(x)]
  f <- strsplit(x, "[[:space:]]+")
  nf <- lengths(f)
  if (length(f) < 2L || any(nf != 3L)) {
    bad <- which(nf != 3L)
    stop("Block file must be a header line plus three whitespace-separated columns chr start stop",
         if (length(bad)) paste0(" (line ", bad[1L], " has ", nf[bad[1L]], ")"), ": ", path, call. = FALSE)
  }
  body <- do.call(rbind, f[-1L])
  bl <- data.table::data.table(src = seq_len(nrow(body)), CHR = .ldb_chr(body[, 1L]),
                               start = suppressWarnings(as.numeric(body[, 2L])),
                               end = suppressWarnings(as.numeric(body[, 3L])))
  if (anyNA(bl$start) || anyNA(bl$end)) stop("start and stop must be numbers: ", path, call. = FALSE)
  bl <- bl[CHR %in% .ldb_auto]
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
    res[[ch]] <- x
  }
  b <- data.table::rbindlist(res)
  b <- b[order(as.integer(CHR), start)]
  b[, Block := seq_len(.N)]
  m <- m[b[, .(src, Block)], on = "src"]
  m <- m[order(Block, BP, REF, ALT)]
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
