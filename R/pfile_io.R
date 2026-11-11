# Reading the .pvar and .psam of a PLINK2 fileset, shared by genotype_impute, genotype_fill and geno_score.

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

# IIDs of a .psam, in file order (the #IID or IID column, read as character); with fid, list(iid, fid), fid "0"
# when the .psam has no FID column.
.psam_ids <- function(pfile, fid = FALSE) {
  psam <- data.table::fread(paste0(pfile, ".psam"), colClasses = "character")
  col <- intersect(c("#IID", "IID"), names(psam))
  if (!length(col)) stop(pfile, ".psam has no IID column", call. = FALSE)
  iid <- psam[[col[1L]]]
  if (!fid) return(iid)
  fc <- intersect(c("#FID", "FID"), names(psam))
  list(iid = iid, fid = if (length(fc)) psam[[fc[1L]]] else rep("0", length(iid)))
}

# 1-based .psam rows of `sample` (all rows when NULL), in the order of `sample`.
.sample_rows <- function(sample, iid) {
  si <- if (is.null(sample)) seq_along(iid) else match(as.character(sample), iid)
  if (anyNA(si)) stop(sum(is.na(si)), " samples not in .psam", call. = FALSE)
  si
}
