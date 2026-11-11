# Editing geno objects in memory: row / column subsets, allele harmonisation (exact flips), missing cells.
# Every result is a new object (new external pointer); the input is never changed.

#' Rows and columns of a geno object as a new geno object
#'
#' The selected cells are copied as stored (missing cells keep their current fill; nothing is refilled) and
#' the rows of the missing cells follow the selection. `center`, `sd`, `n_miss` and `af` are those of the
#' subset. Use [genotype_fill()] afterwards to refill the missing cells from the subset's own samples.
#'
#' @param X A `geno` object.
#' @param sample_vec Rows to keep, as IIDs, `FID_IID` or 1-based row numbers of `X`, in the order wanted.
#'   `NULL` keeps all.
#' @param snp_vec Columns to keep, as column names (`X$snp`) or 1-based column numbers, in the order
#'   wanted. `NULL` keeps all.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A `geno` object.
#' @export
geno_subset <- function(X, sample_vec = NULL, snp_vec = NULL, threads = X$threads) {
  if (!inherits(X, "geno")) stop("X must be a geno object.", call. = FALSE)
  rows <- .geno_index(sample_vec, X$sample, if (!is.null(X$fid)) paste(X$fid, X$sample, sep = "_"), "sample_vec")
  cols <- .geno_index(snp_vec, X$snp, NULL, "snp_vec")
  .geno_take(X, rows, cols, NULL, threads)
}

# Rows x cols (1-based) of X, columns with flip TRUE stored as 2 - x (a1 / a2 swapped, `flipped` toggled).
.geno_take <- function(X, rows, cols, flip = NULL, threads = X$threads) {
  ptr <- geno_take_cpp(X$ptr, as.integer(rows) - 1L, as.integer(cols) - 1L,
                       if (is.null(flip)) logical(0) else as.logical(flip), .threads(threads))
  cl <- .geno_cols(X, cols)
  if (!is.null(flip) && any(flip)) {
    a1 <- cl$a1
    cl$a1[flip] <- cl$a2[flip]; cl$a2[flip] <- a1[flip]
    cl$flipped <- xor(cl$flipped, flip)
  }
  .geno_new(ptr, cl, list(sample = X$sample[rows], fid = X$fid[rows]), impute = X$impute, scale = X$scale,
            threads = X$threads, file = X$file, appended = sort(which(cols %in% X$appended)))
}

#' Put geno objects on one counted allele per variant
#'
#' Variants are matched across objects by chromosome, position and the unordered allele pair (IDs are not
#' used). Without `ref`, the counted allele of each variant is that of the first object holding it, and
#' every other object is flipped there exactly (stored x -> 2 - x on the 2-bit and int16 parts, so the
#' values, missing cells and fills of a flipped column are 2 minus the old ones). With `ref`, only `ref`
#' decides: columns matching a `ref` row are put on its A1, and columns absent from `ref` (or whose allele
#' pair does not match) are left as they are in every object (remove them with [geno_subset()]; the
#' `ref_row` column of the `variants` table is `NA` for them). A variant whose unordered pair occurs twice
#' at one position in any object or in `ref` (e.g. TC/T and T/TC) is matched by its ordered alleles and
#' never flipped. Nothing is dropped and no sample or fill is changed; strand is not considered
#' (palindromic variants are only marked).
#'
#' @param X A named list of `geno` objects (or one object), in order of priority.
#' @param ref Optional SNP table giving the counted allele A1 of its variants. Columns as in the file route
#'   of [genotype_fill()]: `CHR` (or chr / CHROM), `BP37` or `BP38` (or pos / POS / BP; the build of the
#'   objects' positions), `A1`, `A2`, optional `A1Freq` (frequency of A1); "chr" prefixes are removed and
#'   alleles upper-cased.
#' @param build 37 or 38: which of `BP37` / `BP38` to use when `ref` has both.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A list of `geno` objects like `X`, with column names `chr:pos:a1:a2` (a1 = counted allele)
#'   and `a1`, `a2`, `flipped` updated; objects without a flip keep their store. Attribute `variants`: one
#'   row per column name (key, chr, pos, a1, a2, anchor = name of the object (or "ref") fixing a1,
#'   ref_row = row of `ref` (NA if absent), A1Freq = `ref`'s A1Freq (NA if absent), present_<name>,
#'   af_<name> = observed frequency of a1 in that object, palindromic). Map an object's columns to it with
#'   `match(G$snp, variants$key)`.
#' @export
geno_harmonise <- function(X, ref = NULL, build = NULL, threads = NULL) {
  if (inherits(X, "geno")) X <- list(X)
  if (!is.list(X) || !length(X) || !all(vapply(X, inherits, NA, "geno"))) {
    stop("X must be a geno object or a list of geno objects.", call. = FALSE)
  }
  nm <- names(X)
  if (is.null(nm)) nm <- rep("", length(X))
  nm[!nzchar(nm)] <- paste0("G", which(!nzchar(nm)))
  if (anyDuplicated(nm)) stop("The names of X must be unique.", call. = FALSE)
  metas <- lapply(seq_along(X), function(k) {
    G <- X[[k]]
    if (length(G$appended)) stop(nm[k], " has appended columns; harmonise before geno_append().", call. = FALSE)
    if (is.null(G$chr) || anyNA(G$chr) || anyNA(G$pos)) stop(nm[k], " has no chr / pos per column.", call. = FALSE)
    data.table::data.table(chr = sub("^chr", "", as.character(G$chr), ignore.case = TRUE), pos = as.integer(G$pos),
                           A1 = toupper(G$a1), A2 = toupper(G$a2))
  })
  r <- if (!is.null(ref)) .harm_ref(ref, build)
  all_m <- c(if (!is.null(r)) list(r), metas)
  for (m in all_m) m[, cid := paste(chr, pos, pmin(A1, A2), pmax(A1, A2), sep = ":")]
  amb <- unique(unlist(lapply(all_m, function(m) m$cid[duplicated(m$cid)])))
  for (m in all_m) m[, mk := ifelse(cid %in% amb, paste(chr, pos, A1, A2, sep = ":"), cid)]
  # target allele per column: from ref (if given; columns not in ref keep their own), else from the first object
  tgt <- if (!is.null(r)) {
    r[, ref_row := .I]
    lapply(seq_along(X), function(k) {
      m <- metas[[k]]; i <- match(m$mk, r$mk); hit <- !is.na(i)
      data.table::data.table(A1 = ifelse(hit, r$A1[i], m$A1), A2 = ifelse(hit, r$A2[i], m$A2),
                             anchor = ifelse(hit, "ref", nm[k]), ref_row = i, A1Freq = r$A1Freq[i])
    })
  } else {
    all <- data.table::rbindlist(lapply(seq_along(metas), function(k) metas[[k]][, list(mk, A1, A2, k = k)]))
    U <- all[!duplicated(mk)]
    lapply(seq_along(X), function(k) {
      i <- match(metas[[k]]$mk, U$mk)
      data.table::data.table(A1 = U$A1[i], A2 = U$A2[i], anchor = nm[U$k[i]], ref_row = NA_integer_, A1Freq = NA_real_)
    })
  }
  out <- lapply(seq_along(X), function(k) {
    G <- X[[k]]; m <- metas[[k]]; t <- tgt[[k]]
    flip <- m$A1 != t$A1
    H <- if (any(flip)) .geno_take(G, seq_len(G$n), seq_len(G$p), flip, if (is.null(threads)) G$threads else threads) else {
      H <- G; H$cache <- new.env(parent = emptyenv())
      if (!is.null(G$cache$xtx)) H$cache$xtx <- G$cache$xtx
      H
    }
    H$snp <- paste(m$chr, m$pos, t$A1, t$A2, sep = ":"); H$a1 <- t$A1; H$a2 <- t$A2
    H
  })
  names(out) <- nm
  U <- data.table::rbindlist(lapply(seq_along(X), function(k) {
    t <- tgt[[k]]; m <- metas[[k]]
    data.table::data.table(snp_key = out[[k]]$snp, chr = m$chr, pos = m$pos, a1 = t$A1, a2 = t$A2, anchor = t$anchor,
                           ref_row = t$ref_row, A1Freq = t$A1Freq)
  }))
  data.table::setnames(U, "snp_key", "key")  # (data.table()'s own `key` argument would set a sort key)
  U <- U[!duplicated(key)]
  for (k in seq_along(out)) {
    H <- out[[k]]
    i <- match(U$key, H$snp)
    data.table::set(U, j = paste0("present_", nm[k]), value = !is.na(i))
    data.table::set(U, j = paste0("af_", nm[k]), value = H$af[i])
  }
  U[, palindromic := paste(a1, a2) %in% c("A T", "T A", "C G", "G C")]
  data.table::setorder(U, chr, pos, a1, a2)
  attr(out, "variants") <- U[]
  out
}

# A SNP table as given to genotype_fill(pfile, snp) (CHR, BP37 / BP38, A1, A2, optional A1Freq; or chr, pos,
# a1, a2), cleaned the same way: "chr" prefix removed, alleles upper-cased, positions integer.
.harm_ref <- function(ref, build = NULL) {
  r <- data.table::as.data.table(ref)
  pick <- function(cands, what) {
    k <- intersect(cands, names(r))
    if (!length(k)) stop("ref has no ", what, " column.", call. = FALSE)
    r[[k[1L]]]
  }
  pc <- intersect(c("pos", "POS", "BP", "bp"), names(r))
  bc <- intersect(c("BP37", "BP38"), names(r))
  pcol <- if (!is.null(build)) paste0("BP", build) else if (length(pc)) pc[1L] else if (length(bc) == 1L) bc else
    stop("ref needs exactly one of BP37 / BP38 (or pos), or set build", call. = FALSE)
  if (!pcol %in% names(r)) stop("ref has no column ", pcol, call. = FALSE)
  data.table::data.table(chr = sub("^chr", "", as.character(pick(c("chr", "CHR", "CHROM", "#CHROM"), "chr")), ignore.case = TRUE),
                         pos = as.integer(r[[pcol]]),
                         A1 = toupper(pick(c("A1", "a1"), "A1")), A2 = toupper(pick(c("A2", "a2"), "A2")),
                         A1Freq = if ("A1Freq" %in% names(r)) as.numeric(r$A1Freq) else NA_real_)
}

#' Missing cells of a geno object
#'
#' The cells that were missing in the file (and are carried through [geno_subset()] and
#' [geno_harmonise()]), with the value each holds now (its fill, in the object's current counted-allele
#' orientation). The inverse of `genotype_fill(X, cells = )`.
#'
#' @param X A `geno` object.
#' @return data.table(snp, IID, g).
#' @export
geno_missing <- function(X) {
  if (!inherits(X, "geno")) stop("X must be a geno object.", call. = FALSE)
  mc <- geno_miss_cpp(X$ptr)
  j <- rep.int(seq_len(X$p), diff(mc$ptr))
  data.table::data.table(snp = X$snp[j], IID = X$sample[mc$i + 1L], g = mc$g)
}

# The object's variants as a .pvar-like table for the LD panel code: ALT = counted allele, REF = the other.
.geno_vartab <- function(X) {
  if (is.null(X$chr) || is.null(X$pos)) stop("the geno object has no chr / pos per column.", call. = FALSE)
  alt <- toupper(X$a1)
  act <- 2L + lengths(regmatches(alt, gregexpr(",", alt, fixed = TRUE)))
  act[is.na(alt) | seq_len(X$p) %in% X$appended] <- 0L  # appended columns: matched to nothing
  data.table::data.table(CHROM = sub("^chr", "", as.character(X$chr), ignore.case = TRUE), POS = as.integer(X$pos),
                         ID = X$snp, REF = toupper(X$a2), ALT = alt, row = seq_len(X$p), allele_ct = act)
}
