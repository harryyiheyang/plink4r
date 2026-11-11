# The LD panel in one place: reading its SNP table and aligning it to a PGEN.
#
# The panel shipped with plink4r is the Pan-UKB EUR LD (block<N>.eigen.bin, GRCh38 eigen rows) with its SNP
# table (a parquet: ID chr:pos:ref:alt, chr, bp37, bp38, ref, alt, AF_PanUKB = frequency of alt, block,
# ld_index = 0-based eigen row inside the block). A user's PGEN need not share the panel's REF/ALT orientation,
# IDs or build: .ld_align() matches every panel SNP to the PGEN by ID and alleles, else by chromosome, position
# and the allele pair, and records the orientation as a sign (+1: PGEN ALT is the panel alt, -1: swapped).
# Panel SNPs absent from the PGEN are not carried into the imputation at all.

.panel_cols <- c("ID", "chr", "bp37", "bp38", "ref", "alt", "AF_PanUKB", "block", "ld_index")

# Panel SNP table, optionally only some blocks. ref: parquet path or data.frame; NULL (or a directory)
# reads the one .parquet in ld_dir (where the released panel keeps it). With `blocks` only those blocks are
# scanned; otherwise the whole table (panel columns only) is read.
.ld_panel <- function(ref = NULL, ld_dir = NULL, blocks = NULL) {
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
  # required: chr, at least one of bp37 / bp38, ref, alt, AF_PanUKB, block, ld_index; ID is optional
  need <- c("chr", "ref", "alt", "AF_PanUKB", "block", "ld_index")
  if (length(setdiff(need, names(r)))) stop("the panel table lacks ", paste(setdiff(need, names(r)), collapse = ", "), call. = FALSE)
  if (!any(c("bp37", "bp38") %in% names(r))) stop("the panel table has neither bp37 nor bp38", call. = FALSE)
  r <- r[, intersect(.panel_cols, names(r)), with = FALSE]
  if (is.character(r$chr)) r[, chr := sub("^chr", "", chr, ignore.case = TRUE)]
  r[, `:=`(ref = toupper(ref), alt = toupper(alt))]
  if (!"ID" %in% names(r)) r[, ID := paste(chr, if ("bp38" %in% names(r)) bp38 else bp37, ref, alt, sep = ":")]
  if (!is.null(blocks)) r <- r[which(r$block %in% blocks)]
  r
}

# Adds vi (0-based PGEN variant) and sign (+1 / -1) to the panel rows found in the PGEN; others are dropped.
# 1. by ID, kept only when the alleles agree (the sign comes from them);
# 2. the rest by chr + position + allele pair. The position column is the panel's only build column (bp37 or
#    bp38: the PGEN must be on that build); with both, the build that places more SNPs.
# A PGEN variant is used once; multiallelic variants are not aligned. No panel SNP in the PGEN is an error.
.ld_align <- function(panel, pvar) {
  p <- data.table::copy(panel)
  bi <- pvar$allele_ct == 2L
  k <- which(bi)[match(p$ID, pvar$ID[bi])]
  sg <- ifelse(is.na(k), NA_integer_,
               ifelse(pvar$ALT[k] == p$alt & pvar$REF[k] == p$ref, 1L,
                      ifelse(pvar$ALT[k] == p$ref & pvar$REF[k] == p$alt, -1L, NA_integer_)))
  k[is.na(sg)] <- NA_integer_
  p[, `:=`(vi = k - 1L, sign = sg)]
  builds <- intersect(c("bp37", "bp38"), names(p))
  rest <- which(is.na(p$vi))
  if (length(rest)) {
    used <- rep(FALSE, nrow(pvar)); used[p$vi[!is.na(p$vi)] + 1L] <- TRUE
    pv <- pvar[bi & !used, .(CHROM, POS, REF, ALT, row)]
    q <- p[rest, .(i = rest, CHROM = as.character(chr), ref, alt)]
    best <- NULL
    for (b in builds) {
      qb <- data.table::copy(q)[, POS := as.integer(p[[b]][i])]
      same <- pv[qb, on = .(CHROM, POS, REF = ref, ALT = alt), nomatch = NULL, .(i = i.i, row = x.row, s = 1L)]
      swap <- pv[qb, on = .(CHROM, POS, REF = alt, ALT = ref), nomatch = NULL, .(i = i.i, row = x.row, s = -1L)]
      h <- rbind(same, swap)
      h <- h[!duplicated(i) & !duplicated(row)]
      if (is.null(best) || nrow(h) > nrow(best)) best <- h
    }
    if (nrow(best)) p[best$i, `:=`(vi = best$row - 1L, sign = best$s)]
  }
  p <- p[!is.na(vi)]
  if (!nrow(p)) stop("no panel SNP found in the PGEN (by ID, or by chr / ", paste(builds, collapse = " or "),
                     " / alleles): check the build and the variant IDs", call. = FALSE)
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
