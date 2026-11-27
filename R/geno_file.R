#' Open a BED or PGEN file: its variant and sample tables, read once
#'
#' Reads `.bim` / `.fam` (BED) or `.pvar` / `.psam` (PGEN) once, in compiled code, without any genotypes. Every
#' function of the package that reads a genotype file takes the result as `file =` and then reads neither table
#' again; a path given in its place is opened by the same reader, once in that call.
#'
#' A `.pvar` line that names several ALT alleles is one variant per ALT allele, in the order the line names them:
#' each counts that one allele (a sample that carries another ALT of the line counts 0 of it), so every variant of
#' the table has one counted allele `A1` and one other allele `A2`, and every function works on it as on a line
#' with one ALT. A file without such lines has one variant per line. No variant ID is read: a variant is `CHR`,
#' `BP`, `A1`, `A2`.
#'
#' @param bedfile PLINK 1 BED path or prefix (`.bim` and `.fam` alongside).
#' @param pgenfile PLINK 2 PGEN path or prefix (`.pvar` or `.pvar.zst`, and `.psam` alongside).
#' @return A list of class `geno_file`: `type` ("bed" / "pgen"), `path` (the .bed / .pgen), `prefix`, `ptr` (the
#'   tables as the compiled code holds them: an external pointer, valid in this R session only, so the object is
#'   neither saved with `saveRDS` nor sent to parallel workers), `variants` = data.table(row, CHR, BP, A1, A2) with
#'   `A1` the counted allele (bim column 5 / one ALT allele) and `A2` the other (bim column 6 / REF), both in upper
#'   case, `CHR` the integer chromosome code: 1..22, 23 = X (a file's X, 23, PAR1, PAR2, XY and 25 are all 23),
#'   24 = Y, 26 = MT, `NA` for any other label, and `row` the place of the variant in this table; and `samples` =
#'   data.table(row, fid, iid). Select windows by row, e.g. `gf$variants[CHR == 16 & BP >= lo & BP <= hi, row]`,
#'   and pass them as `snp_vec`.
#' @export
geno_file <- function(bedfile = NULL, pgenfile = NULL) {
  if (is.null(bedfile) == is.null(pgenfile)) stop("Give exactly one of bedfile and pgenfile.", call. = FALSE)
  tmp <- NULL
  if (!is.null(bedfile)) {
    path <- .geno_bed_path(bedfile)
    base <- sub("\\.bed$", "", path)
    text <- paste0(base, c(".bim", ".fam"))
    type <- "bed"
  } else {
    files <- .geno_pgen_files(pgenfile)
    base <- files$prefix; path <- files$pgen
    pvar <- paste0(base, ".pvar")
    if (!file.exists(pvar)) {
      if (!file.exists(paste0(pvar, ".zst"))) stop("no .pvar or .pvar.zst for ", base, call. = FALSE)
      pvar <- tmp <- .zst_expand(paste0(pvar, ".zst"))
      on.exit(unlink(tmp), add = TRUE)
    }
    text <- c(pvar, paste0(base, ".psam"))
    type <- "pgen"
  }
  for (f in text) if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
  ptr <- geno_table_open_cpp(path, type == "bed", normalizePath(text[1L], winslash = "/"),
                             normalizePath(text[2L], winslash = "/"))
  v <- geno_table_variants_cpp(ptr)
  s <- geno_table_samples_cpp(ptr)
  structure(list(type = type, path = path, prefix = base, ptr = ptr,
                 variants = data.table::data.table(row = seq_along(v$BP), CHR = v$CHR, BP = v$BP, A1 = v$alleles[v$A1],
                                                   A2 = v$alleles[v$A2]),
                 samples = data.table::data.table(row = seq_along(s$iid), fid = s$fid, iid = s$iid)),
            class = "geno_file")
}

# The geno_file behind a call: `file` as given, else opened from bedfile / pgenfile.
.as_geno_file <- function(file = NULL, bedfile = NULL, pgenfile = NULL) {
  if (is.null(file)) return(geno_file(bedfile = bedfile, pgenfile = pgenfile))
  if (!inherits(file, "geno_file")) stop("file must be a geno_file() result.", call. = FALSE)
  if (!is.null(bedfile) || !is.null(pgenfile)) stop("Give file or a path, not both.", call. = FALSE)
  file
}

# The name of a variant: CHR:BP:A1:A2.
.geno_vid <- function(v) paste(v$CHR, v$BP, v$A1, v$A2, sep = ":")

# Rows of the names x (CHR:BP:A1:A2) in a variants table, as match(x, .geno_vid(v)) but naming only the rows at the
# positions in x (the second field of a name is its row's BP, so every row with that name is among them); if a name
# is not found, through all names.
.geno_vid_index <- function(v, x) {
  pos <- suppressWarnings(as.integer(sub("^[^:]*:([^:]*):.*$", "\\1", x)))
  if (!anyNA(pos)) {
    cand <- which(v$BP %in% pos)
    idx <- cand[match(x, .geno_vid(v[cand]))]
    if (!anyNA(idx)) return(idx)
  }
  .geno_index(x, .geno_vid(v), NULL, "snp_vec")
}

# The variants (rows of the table) that snp_vec names: NULL = all, names CHR:BP:A1:A2, or rows.
.geno_variants <- function(gf, snp_vec) {
  v <- gf$variants
  if (is.character(snp_vec)) .geno_vid_index(v, snp_vec) else .geno_index(snp_vec, v$row, NULL, "snp_vec")
}

# The samples (file rows) that sample_vec names: NULL = all, IIDs (or FID_IID), or rows.
.geno_samples <- function(gf, sample_vec) {
  s <- gf$samples
  .geno_index(sample_vec, s$iid, paste(s$fid, s$iid, sep = "_"), "sample_vec")
}

# Per row of a variants table (CHR, BP, A1, A2): whether another row has its CHR, BP and A2. Such rows are the ALT
# alleles of one position (a .pvar line with several ALT alleles, or the same written on several lines): 2 - x of
# one of them counts "not A1", which is not the count of A2, so it is never matched in the swapped orientation.
.several <- function(v) {
  k <- data.table::data.table(CHR = v$CHR, BP = v$BP, A2 = v$A2)
  duplicated(k) | duplicated(k, fromLast = TRUE)
}
