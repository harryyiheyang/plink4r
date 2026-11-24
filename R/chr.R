# The chromosome code, in one place.
#
# Every variant table of the package holds its chromosome as this integer, and every match between two tables is
# made on it: 1..22 for the autosomes; 23 for X (X or 23; the LD panel and summary statistics write X this way);
# 24 for Y; 25 for the pseudo-autosomal parts of X (XY, PAR1, PAR2 or 25), kept apart from 23 because a PGEN codes
# males differently there (diploid) than on X (haploid); 26 for the mitochondrion (MT, M). A leading "chr" (any
# case) is ignored. Any other label has no code (NA) and is matched to nothing.

# The code of chromosome labels or numbers.
.chr_code <- function(chr) {
  u <- unique(chr)
  x <- toupper(sub("^chr", "", as.character(u), ignore.case = TRUE))
  code <- rep(NA_integer_, length(x))
  num <- grepl("^[0-9]{1,2}$", x)
  code[num] <- as.integer(x[num])
  named <- c(X = 23L, Y = 24L, XY = 25L, PAR1 = 25L, PAR2 = 25L, MT = 26L, M = 26L)
  code[!num] <- unname(named[x[!num]])
  code[!code %in% 1:26] <- NA_integer_
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
