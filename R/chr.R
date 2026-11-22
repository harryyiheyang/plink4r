# The chromosome code, in one place.
#
# Every variant table of the package holds its chromosome as this integer, and every match between two tables is
# made on it: 1..22 for the autosomes; 23 for X including its pseudo-autosomal parts (which files write as X, 23,
# PAR1, PAR2, XY or 25, with their positions on X); 24 for Y; 26 for the mitochondrion (MT, M). A leading "chr"
# (any case) is ignored. Any other label has no code (NA) and is matched to nothing.

# The code of chromosome labels or numbers.
.chr_code <- function(chr) {
  u <- unique(chr)
  x <- toupper(sub("^chr", "", as.character(u), ignore.case = TRUE))
  code <- rep(NA_integer_, length(x))
  num <- grepl("^[0-9]{1,2}$", x)
  code[num] <- as.integer(x[num])
  named <- c(X = 23L, PAR1 = 23L, PAR2 = 23L, XY = 23L, Y = 24L, MT = 26L, M = 26L)
  code[!num] <- unname(named[x[!num]])
  code[code %in% 25L] <- 23L
  code[!code %in% c(1:24, 26L)] <- NA_integer_
  code[match(chr, u)]
}

# Summary statistics must cover chromosomes 1 to 22; X (23) may be absent. chr: the codes of the SNPs of `what`.
.need_autosomes <- function(chr, what) {
  miss <- setdiff(1:22, chr)
  if (length(miss))
    stop(what, " has no SNP on chromosome ", paste(miss, collapse = ", "),
         ": chromosomes 1 to 22 are required (X is optional).", call. = FALSE)
  invisible(NULL)
}
