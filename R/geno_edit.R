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
#' used). The counted allele of each variant is that of the first object holding it (`ref` first when
#' given), and every other object is flipped there exactly (stored x -> 2 - x on the 2-bit and int16 parts,
#' so the values, missing cells and fills of a flipped column are 2 minus the old ones). A variant whose
#' unordered pair occurs twice at one position in any object (e.g. TC/T and T/TC) is matched by its ordered
#' alleles and never flipped. Nothing is dropped and no sample or fill is changed; strand is not
#' considered (palindromic variants are only marked).
#'
#' @param X A named list of `geno` objects (or one object), in order of priority.
#' @param ref Optional variant table with columns chr, pos, a1, a2 (also CHR / CHROM, POS / BP, A1, A2)
#'   giving the counted allele `a1` of its variants, ahead of every object. This is how an object is put
#'   on the A1 of a SNP table.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A list of `geno` objects like `X`, with column names `chr:pos:a1:a2` (a1 = counted allele)
#'   and `a1`, `a2`, `flipped` updated; objects without a flip keep their store. Attribute `variants`: the
#'   union table (key, chr, pos, a1, a2, anchor = name of the object (or "ref") fixing a1, present_<name>,
#'   af_<name> = observed frequency of a1 in that object, palindromic).
#' @export
geno_harmonise <- function(X, ref = NULL, threads = NULL) {
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
  src <- nm
  if (!is.null(ref)) {
    metas <- c(list(.harm_ref(ref)), metas)
    src <- c("ref", nm)
  }
  for (m in metas) m[, cid := paste(chr, pos, pmin(A1, A2), pmax(A1, A2), sep = ":")]
  amb <- unique(unlist(lapply(metas, function(m) m$cid[duplicated(m$cid)])))
  for (m in metas) m[, mk := ifelse(cid %in% amb, paste(chr, pos, A1, A2, sep = ":"), cid)]
  all <- data.table::rbindlist(lapply(seq_along(metas), function(k) metas[[k]][, list(mk, chr, pos, A1, A2, k = k)]))
  U <- all[!duplicated(mk)]
  off <- if (is.null(ref)) 0L else 1L
  U <- U[mk %in% unlist(lapply(metas[off + seq_along(X)], `[[`, "mk"))]
  U[, `:=`(key = paste(chr, pos, A1, A2, sep = ":"), anchor = src[k])]
  out <- lapply(seq_along(X), function(k) {
    G <- X[[k]]; m <- metas[[off + k]]
    j <- match(m$mk, U$mk)
    flip <- m$A1 != U$A1[j]
    H <- if (any(flip)) .geno_take(G, seq_len(G$n), seq_len(G$p), flip, if (is.null(threads)) G$threads else threads) else {
      H <- G; H$cache <- new.env(parent = emptyenv())
      if (!is.null(G$cache$xtx)) H$cache$xtx <- G$cache$xtx
      H
    }
    H$snp <- U$key[j]; H$a1 <- U$A1[j]; H$a2 <- U$A2[j]
    H
  })
  names(out) <- nm
  for (k in seq_along(out)) {
    H <- out[[k]]
    i <- match(U$key, H$snp)
    data.table::set(U, j = paste0("present_", nm[k]), value = !is.na(i))
    data.table::set(U, j = paste0("af_", nm[k]), value = H$af[i])
  }
  U[, palindromic := paste(A1, A2) %in% c("A T", "T A", "C G", "G C")]
  U[, c("mk", "k") := NULL]
  data.table::setnames(U, c("A1", "A2"), c("a1", "a2"))
  data.table::setcolorder(U, c("key", "chr", "pos", "a1", "a2", "anchor"))
  data.table::setorder(U, chr, pos, a1, a2)
  attr(out, "variants") <- U[]
  out
}

.harm_ref <- function(ref) {
  r <- data.table::as.data.table(ref)
  pick <- function(cands, what) {
    k <- intersect(cands, names(r))
    if (!length(k)) stop("ref has no ", what, " column.", call. = FALSE)
    r[[k[1L]]]
  }
  data.table::data.table(chr = sub("^chr", "", as.character(pick(c("chr", "CHR", "CHROM", "#CHROM"), "chr")), ignore.case = TRUE),
                         pos = as.integer(pick(c("pos", "POS", "BP", "bp"), "pos")),
                         A1 = toupper(pick(c("a1", "A1"), "a1")), A2 = toupper(pick(c("a2", "A2"), "a2")))
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
