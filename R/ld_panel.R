# The LD panel in one place: reading its SNP table and aligning it to a PGEN.
#
# The panel shipped with plink4r is the Pan-UKB EUR LD (block<N>.eigen.bin, GRCh38 eigen rows) with its SNP
# table (a parquet: ID chr:pos:ref:alt, chr, bp37, bp38, ref, alt, AF_PanUKB = frequency of alt, block,
# ld_index = 0-based eigen row inside the block). A user's PGEN need not share the panel's REF/ALT orientation,
# IDs or build: .ld_align() matches every panel SNP to the PGEN by ID and alleles, else by chromosome, position
# and the allele pair, and records the orientation as a sign (+1: PGEN ALT is the panel alt, -1: swapped).
# Panel SNPs absent from the PGEN are not carried into the imputation at all.

.panel_cols <- c("ID", "chr", "bp37", "bp38", "ref", "alt", "AF_PanUKB", "block", "ld_index")

# Panel SNP table, optionally only some blocks / IDs. ref: parquet path or data.frame; NULL (or a directory)
# reads the one .parquet in ld_dir (where the released panel keeps it). With `blocks` only those blocks are
# scanned; otherwise the whole table (panel columns only) is read.
.ld_panel <- function(ref = NULL, ld_dir = NULL, blocks = NULL, ids = NULL) {
  if (is.character(ref) && length(ref) == 1L && dir.exists(ref)) { ld_dir <- ref; ref <- NULL }
  if (is.null(ref)) {
    if (is.null(ld_dir)) stop("give ld_dir (or ref)", call. = FALSE)
    f <- list.files(ld_dir, pattern = "\\.parquet$", full.names = TRUE)
    if (length(f) != 1L) stop("expected one .parquet SNP table in ", ld_dir, ", found ", length(f), call. = FALSE)
    ref <- f
  }
  if (is.data.frame(ref)) {
    r <- data.table::as.data.table(ref)
  } else {
    ds <- arrow::open_dataset(ref)
    cols <- intersect(.panel_cols, names(ds))
    r <- if (is.null(blocks)) data.table::as.data.table(arrow::read_parquet(ref, col_select = cols)) else {
      # only the requested blocks (the released parquet is sorted by block, so whole row groups are skipped)
      inb <- arrow::Expression$create("is_in", arrow::Expression$field_ref("block"),
                                      options = list(value_set = arrow::Array$create(as.integer(blocks)), skip_nulls = TRUE))
      data.table::as.data.table(as.data.frame(ds$NewScan()$Filter(inb)$Project(cols)$Finish()$ToTable()))
    }
  }
  need <- c("ID", "block", "ld_index", "AF_PanUKB")
  if (length(setdiff(need, names(r)))) stop("the panel table lacks ", paste(setdiff(need, names(r)), collapse = ", "), call. = FALSE)
  r <- r[, intersect(.panel_cols, names(r)), with = FALSE]
  if (!is.null(blocks)) r <- r[which(r$block %in% blocks)]
  if (!is.null(ids)) r <- r[which(r$ID %in% ids)]
  r
}

# CHROM POS ID REF ALT of a .pvar (or .pvar.zst), with row (1-based) and allele_ct.
.read_pvar5 <- function(pfile) {
  pf <- paste0(pfile, ".pvar")
  pvar <- if (file.exists(pf)) data.table::fread(pf, skip = "#CHROM", select = 1:5, colClasses = "character") else
    if (file.exists(paste0(pf, ".zst")))
      data.table::fread(cmd = paste("zstdcat", shQuote(paste0(pf, ".zst")), "| grep -v '^##'"), select = 1:5,
                        colClasses = "character") else stop("no .pvar or .pvar.zst for ", pfile, call. = FALSE)
  data.table::setnames(pvar, c("CHROM", "POS", "ID", "REF", "ALT"))
  pvar[, `:=`(row = .I, CHROM = sub("^chr", "", CHROM, ignore.case = TRUE), POS = as.integer(POS),
              REF = toupper(REF), ALT = toupper(ALT))]
  pvar[, allele_ct := 2L + lengths(regmatches(ALT, gregexpr(",", ALT, fixed = TRUE)))]
  pvar
}

# Adds vi (0-based PGEN variant) and sign (+1 / -1) to the panel rows found in the PGEN; others are dropped.
# 1. by ID: kept when the panel has no ref/alt columns (sign +1), else the alleles decide the sign;
# 2. the rest by chr + position + allele pair. The position column is the panel's only build column (bp37 or
#    bp38: the PGEN must be on that build); with both, the build that places more SNPs.
# A PGEN variant is used once; multiallelic variants are not aligned. No panel SNP in the PGEN is an error.
.ld_align <- function(panel, pvar) {
  p <- data.table::copy(panel)
  p[, `:=`(vi = NA_integer_, sign = NA_integer_)]
  bi <- pvar$allele_ct == 2L
  has_al <- all(c("ref", "alt") %in% names(p))
  if (has_al) p[, `:=`(ref_u = toupper(ref), alt_u = toupper(alt))]
  k <- which(bi)[match(p$ID, pvar$ID[bi])]
  if (has_al) {
    sg <- ifelse(is.na(k), NA_integer_,
                 ifelse(pvar$ALT[k] == p$alt_u & pvar$REF[k] == p$ref_u, 1L,
                        ifelse(pvar$ALT[k] == p$ref_u & pvar$REF[k] == p$alt_u, -1L, NA_integer_)))
    k[is.na(sg)] <- NA_integer_
  } else sg <- ifelse(is.na(k), NA_integer_, 1L)
  p[, `:=`(vi = k - 1L, sign = sg)]
  builds <- intersect(c("bp37", "bp38"), names(p))
  rest <- which(is.na(p$vi))
  if (length(rest) && has_al && "chr" %in% names(p) && length(builds)) {
    used <- rep(FALSE, nrow(pvar)); used[p$vi[!is.na(p$vi)] + 1L] <- TRUE
    pv <- pvar[bi & !used, .(CHROM, POS, REF, ALT, row)]
    q <- p[rest, .(i = rest, CHROM = sub("^chr", "", as.character(chr), ignore.case = TRUE), ref_u, alt_u)]
    best <- NULL
    for (b in builds) {
      qb <- data.table::copy(q)[, POS := as.integer(p[[b]][i])]
      same <- pv[qb, on = .(CHROM, POS, REF = ref_u, ALT = alt_u), nomatch = NULL, .(i = i.i, row = x.row, s = 1L)]
      swap <- pv[qb, on = .(CHROM, POS, REF = alt_u, ALT = ref_u), nomatch = NULL, .(i = i.i, row = x.row, s = -1L)]
      h <- rbind(same, swap)
      h <- h[!duplicated(i) & !duplicated(row)]
      if (is.null(best) || nrow(h) > nrow(best)) best <- h
    }
    if (nrow(best)) p[best$i, `:=`(vi = best$row - 1L, sign = best$s)]
  }
  p <- p[!is.na(vi)]
  if (!nrow(p)) stop("no panel SNP found in the PGEN (by ID, or by chr / ", paste(builds, collapse = " or "),
                     " / alleles): check the build and the variant IDs", call. = FALSE)
  if (has_al) p[, c("ref_u", "alt_u") := NULL]
  p
}

# Working rows of one block: the panel rows of `block` (eigen rows lo..hi when given) that are in the PGEN.
# Returns rows (0-based eigen row), vi, sign, af (PGEN-ALT frequency) and the PGEN ID, by eigen row.
.ld_rows <- function(aligned, blk, lo = NULL, hi = NULL) {
  a <- aligned[which(aligned$block == blk)]
  if (!is.null(lo)) a <- a[which(a$ld_index >= lo & a$ld_index <= hi)]
  data.table::setorder(a, ld_index)
  a[!duplicated(ld_index)]
}

# The ridge delta of every LD imputation must be > 0 (delta = 0 has no stable C = I - delta (R + delta I)^-1).
.check_delta <- function(delta) {
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0)
    stop("delta must be a single number > 0", call. = FALSE)
  invisible(delta)
}

# bp37 / bp38 present in every output table; a build the panel lacks is NA.
.with_bp <- function(r) {
  for (b in c("bp37", "bp38")) if (!b %in% names(r)) data.table::set(r, j = b, value = NA_integer_)
  r
}
