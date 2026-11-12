#' Read the variant and sample tables of a BED or PGEN file once
#'
#' Reads `.bim` / `.fam` (BED) or `.pvar` / `.psam` (PGEN) once, without any genotypes. Give the result to
#' [geno_open()] as `file =` to open many windows or sample sets of the same file without re-reading the
#' tables (a 10-million-row `.bim` takes seconds; an open from the tables takes milliseconds).
#'
#' @param bedfile PLINK 1 BED path or prefix (`.bim` and `.fam` alongside).
#' @param pgenfile PLINK 2 PGEN path or prefix (`.pvar` or `.pvar.zst`, and `.psam` alongside).
#' @param id Whether to read the variant ID column. `FALSE` saves time and memory on large tables (IDs are
#'   often duplicated, e.g. "."); objects opened from it name their columns `chr:pos:a1:a2` and `snp_vec`
#'   then selects by those names or by row.
#' @return A list of class `geno_file`: `type` ("bed" / "pgen"), `path` (the .bed / .pgen), `prefix`,
#'   `variants` = data.table(row, chr, pos, id, a1, a2, allele_ct) with `a1` the counted allele (bim
#'   column 5 / ALT) and `a2` the other (bim column 6 / REF), `chr` without a leading "chr"; and
#'   `samples` = data.table(row, fid, iid). Select windows by row, e.g.
#'   `gf$variants[chr == "16" & pos >= lo & pos <= hi, row]`, and pass them as `snp_vec`.
#' @export
geno_file <- function(bedfile = NULL, pgenfile = NULL, id = TRUE) {
  if (is.null(bedfile) == is.null(pgenfile)) stop("Give exactly one of bedfile and pgenfile.", call. = FALSE)
  if (!is.null(bedfile)) {
    bed <- .geno_bed_path(bedfile)
    base <- sub("\\.bed$", "", bed)
    sel <- if (id) c(1L, 2L, 4L, 5L, 6L) else c(1L, 4L, 5L, 6L)
    v <- data.table::fread(paste0(base, ".bim"), header = FALSE, select = sel,
                           colClasses = list(character = intersect(c(1L, 2L, 5L, 6L), sel), integer = 4L))
    data.table::setnames(v, if (id) c("chr", "id", "pos", "a1", "a2") else c("chr", "pos", "a1", "a2"))
    v[, `:=`(row = .I, chr = sub("^chr", "", chr, ignore.case = TRUE), allele_ct = 2L)]
    fam <- data.table::fread(paste0(base, ".fam"), header = FALSE, select = 1:2, colClasses = "character")
    s <- data.table::data.table(row = seq_len(nrow(fam)), fid = fam[[1L]], iid = fam[[2L]])
    type <- "bed"; path <- bed
  } else {
    files <- .geno_pgen_files(pgenfile)
    base <- files$prefix
    pv <- .read_pvar5(base, id = id)
    v <- data.table::setnames(pv, c("CHROM", "POS", "ID", "ALT", "REF"), c("chr", "pos", "id", "a1", "a2"), skip_absent = TRUE)
    ps <- .psam_ids(base, fid = TRUE)
    s <- data.table::data.table(row = seq_along(ps$iid), fid = ps$fid, iid = ps$iid)
    type <- "pgen"; path <- files$pgen
  }
  data.table::setcolorder(v, intersect(c("row", "chr", "pos", "id", "a1", "a2", "allele_ct"), names(v)))
  structure(list(type = type, path = path, prefix = base, variants = v, samples = s), class = "geno_file")
}

# The geno_file behind a geno_open call: `file` as given, else read from bedfile / pgenfile.
.as_geno_file <- function(file = NULL, bedfile = NULL, pgenfile = NULL) {
  if (is.null(file)) return(geno_file(bedfile = bedfile, pgenfile = pgenfile))
  if (!inherits(file, "geno_file")) stop("file must be a geno_file() result.", call. = FALSE)
  if (!is.null(bedfile) || !is.null(pgenfile)) stop("Give file or bedfile / pgenfile, not both.", call. = FALSE)
  file
}

# Column names of a geno_file's variants: the file ID, or chr:pos:a1:a2 without IDs.
.geno_vid <- function(v) if ("id" %in% names(v)) v$id else paste(v$chr, v$pos, v$a1, v$a2, sep = ":")

# Rows of chr:pos:a1:a2 names x in a variants table without IDs, as match(x, .geno_vid(v)) but naming only the rows
# at the positions in x (when no chr contains ":", the second field of a name is its row's pos, so every row with
# that name is among them); otherwise, or if a name is not found, through all names.
.geno_vid_index <- function(v, x) {
  pos <- suppressWarnings(as.integer(sub("^[^:]*:([^:]*):.*$", "\\1", x)))
  if (!anyNA(pos) && !any(grepl(":", unique(v$chr), fixed = TRUE))) {
    cand <- which(v$pos %in% pos)
    idx <- cand[match(x, .geno_vid(v[cand]))]
    if (!anyNA(idx)) return(idx)
  }
  .geno_index(x, .geno_vid(v), NULL, "snp_vec")
}

# Column of X (or NA) for each variant given by chr, pos and its two alleles a / b (either orientation):
# chromosome without "chr", alleles upper-case; a column whose counted allele is a matches first, then one
# whose counted allele is b. IDs and column names are not used.
.geno_match_vars <- function(X, chr, pos, a, b) {
  key <- function(c, p, x, y) paste(sub("^chr", "", as.character(c), ignore.case = TRUE), as.integer(p),
                                    toupper(x), toupper(y), sep = ":")
  kx <- key(X$chr, X$pos, X$a1, X$a2)
  kx[is.na(X$chr) | is.na(X$pos)] <- NA_character_  # appended columns match nothing
  j <- match(key(chr, pos, a, b), kx, incomparables = NA)
  miss <- is.na(j)
  j[miss] <- match(key(chr, pos, b, a)[miss], kx, incomparables = NA)
  j
}
