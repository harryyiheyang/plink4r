# The sample ID as text, in one place.
#
# Sample IDs are matched as text with those of the genotype file. as.character() writes a round whole number held as
# a double in exponent form (100000 as "1e+05"), which is no ID of any file: a table whose reader holds whole numbers
# as doubles would lose those samples without a word.

# The IDs of a column as text: a double written out in full (its 15 significant digits), anything else by
# as.character().
.id_text <- function(x) {
  if (!is.double(x) || is.object(x)) return(as.character(x))
  out <- formatC(x, digits = 15L, width = 1L, format = "fg")
  out[is.na(x)] <- NA_character_
  out
}
