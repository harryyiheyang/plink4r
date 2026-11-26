# Reading the .pvar and .psam of a PLINK2 fileset, for geno_file().

# CHROM (the chromosome code, R/chr.R) POS ID REF ALT of a .pvar (or .pvar.zst), with row (1-based) and allele_ct;
# id = FALSE skips the ID column.
.read_pvar5 <- function(pfile, id = TRUE) {
  pf <- paste0(pfile, ".pvar")
  sel <- if (id) 1:5 else c(1L, 2L, 4L, 5L)
  cc <- list(character = setdiff(sel, 2L), integer = 2L)
  pvar <- if (file.exists(pf)) data.table::fread(pf, skip = "#CHROM", select = sel, colClasses = cc) else
    if (file.exists(paste0(pf, ".zst"))) .read_pvar_zst(paste0(pf, ".zst"), sel, cc) else
      stop("no .pvar or .pvar.zst for ", pfile, call. = FALSE)
  data.table::setnames(pvar, c("CHROM", "POS", "ID", "REF", "ALT")[sel])
  # per distinct value: a 22M-row .pvar has few distinct CHROM / REF / ALT strings
  pvar[, `:=`(row = .I, CHROM = .chr_code(CHROM), POS = as.integer(POS), REF = .by_unique(REF, toupper),
              ALT = .by_unique(ALT, toupper))]
  ac <- rep(2L, nrow(pvar))
  w <- grep(",", pvar$ALT, fixed = TRUE)
  if (length(w)) ac[w] <- 2L + lengths(regmatches(pvar$ALT[w], gregexpr(",", pvar$ALT[w], fixed = TRUE)))
  pvar[, allele_ct := ac]
  pvar
}

# f(x) for an elementwise string function f, computed once per distinct value.
.by_unique <- function(x, f) {
  u <- unique(x)
  fu <- f(u)
  if (identical(fu, u)) x else fu[match(x, u)]
}

# A .pvar.zst read through R's zstd connection (R >= 4.5) into a temporary file, else through zstdcat.
.read_pvar_zst <- function(path, sel, cc) {
  if (!exists("zstdfile", baseenv())) {
    return(data.table::fread(cmd = paste("zstdcat", shQuote(path)), skip = "#CHROM", select = sel, colClasses = cc))
  }
  tmp <- tempfile(fileext = ".pvar")
  on.exit(unlink(tmp))
  con <- get("zstdfile", baseenv())(path, "rb")
  out <- file(tmp, "wb")
  repeat {
    b <- readBin(con, "raw", 64L * 1024L^2)
    if (!length(b)) break
    writeBin(b, out)
  }
  close(out); close(con)
  data.table::fread(tmp, skip = "#CHROM", select = sel, colClasses = cc)
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
