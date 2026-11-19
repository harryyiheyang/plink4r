# Binding geno objects: columns side by side (geno_cbind) and rows one below the other (geno_rbind).
# Both copy every cell as stored through one C++ routine (geno_bind_cpp) and carry the per-column and
# per-row fields (the object's stand-ins for the .bim / .pvar and .fam / .psam) along.

# Objects given as geno_xbind(A, B) or geno_xbind(list(A, B)).
.bind_objects <- function(...) {
  X <- list(...)
  if (length(X) == 1L && !inherits(X[[1L]], "geno") && is.list(X[[1L]])) X <- X[[1L]]
  if (!length(X) || !all(vapply(X, inherits, NA, "geno"))) stop("give geno objects, or a list of them.", call. = FALSE)
  X
}

# Rank of a chromosome label in PLINK order: 1, 2, ... by number, then X, Y, XY, MT, then the other labels.
.chr_rank <- function(chr) {
  x <- toupper(sub("^chr", "", as.character(chr), ignore.case = TRUE))
  num <- suppressWarnings(as.integer(x))
  sex <- match(x, c("X", "Y", "XY", "PAR1", "PAR2", "MT", "M"))
  r <- ifelse(!is.na(num), num, 1e6 + sex)
  other <- is.na(r) & !is.na(x)
  r[other] <- 2e6 + match(x[other], sort(unique(x[other])))
  r
}

#' Bind geno objects by columns
#'
#' The objects must hold the same samples in the same order (IID and FID). The columns of the result are
#' those of all objects ordered by chromosome (PLINK order), position, counted allele `a1`, then `a2` (the
#' order of [geno_harmonise()]), whatever the order of the objects; appended columns ([geno_append()], no
#' position) come last in the order given. A variant (chr, pos, a1, a2) held by two objects is an
#' error. Every cell is copied as stored: the fills of the missing cells are kept.
#'
#' @param ... `geno` objects, or one list of them (e.g. the result of [genotype_fill()]).
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A `geno` object; `scale`, `impute` and `file` are those of the first object (`impute = "mixed"`,
#'   `file = NULL` when the objects differ).
#' @export
geno_cbind <- function(..., threads = 4L) {
  X <- .bind_objects(...)
  K <- length(X)
  for (k in seq_len(K)[-1L])
    if (!identical(X[[k]]$sample, X[[1L]]$sample) || !identical(X[[k]]$fid, X[[1L]]$fid))
      stop("object ", k, " does not hold the samples of object 1 in the same order (see geno_subset()).", call. = FALSE)
  src <- rep(seq_len(K), vapply(X, function(x) x$p, 1L, USE.NAMES = FALSE))
  col <- unlist(lapply(X, function(x) seq_len(x$p)), use.names = FALSE)
  cat_f <- function(f, na) unlist(lapply(X, function(x) if (is.null(x[[f]])) rep(na, x$p) else x[[f]]), use.names = FALSE)
  cols <- list(snp = cat_f("snp", NA_character_), id = cat_f("id", NA_character_), chr = cat_f("chr", NA_character_),
               pos = cat_f("pos", NA_integer_), a1 = cat_f("a1", NA_character_), a2 = cat_f("a2", NA_character_),
               row = cat_f("row", NA_integer_), flipped = cat_f("flipped", NA))
  app <- unlist(lapply(X, function(x) seq_len(x$p) %in% x$appended), use.names = FALSE)
  key <- paste(cols$chr, cols$pos, cols$a1, cols$a2)[!app]
  if (anyDuplicated(key)) stop("a variant (chr, pos, a1, a2) is held by two objects: ", key[anyDuplicated(key)], call. = FALSE)
  o <- order(app, .chr_rank(cols$chr), cols$pos, cols$a1, cols$a2, method = "radix")
  M <- matrix(-1L, K, length(o))
  M[cbind(src[o], seq_along(o))] <- col[o] - 1L
  ptr <- geno_bind_cpp(lapply(X, `[[`, "ptr"), M, integer(0), integer(0), .threads(threads))
  same <- function(f) all(vapply(X, function(x) identical(x[[f]], X[[1L]][[f]]), NA))
  .geno_new(ptr, lapply(cols, `[`, o), list(sample = X[[1L]]$sample, fid = X[[1L]]$fid),
            impute = if (same("impute")) X[[1L]]$impute else "mixed", scale = X[[1L]]$scale, threads = X[[1L]]$threads,
            file = if (same("file")) X[[1L]]$file, appended = which(app[o]))
}

#' Bind geno objects by rows
#'
#' The objects must hold the same columns in the same order and orientation (chr, pos, a1, a2; bring them
#' there with [geno_harmonise()] and [geno_subset()]), all hard calls or all dosages.
#' * `sample = NULL`: the rows of the objects one below the other, in the order given; of a sample (FID and
#'   IID) held more than once only the first row is kept.
#' * `sample` given: the rows of the samples of `sample` that the objects hold, in the order of `sample`
#'   (first row of a sample held more than once); samples not in `sample` are dropped.
#' Every cell is copied as stored: a missing cell keeps the fill it had in its own object ([geno_fill()]
#' refills them from the bound samples).
#'
#' @param ... `geno` objects, or one list of them.
#' @param sample `NULL`, or the samples to keep and their order: a data.frame whose first two columns are FID
#'   and IID, or the path of such a file (.fam, .psam, or any table with FID and IID first).
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A `geno` object with the column fields of the first object; `impute` and `file` as in [geno_cbind()].
#' @export
geno_rbind <- function(..., sample = NULL, threads = 4L) {
  X <- .bind_objects(...)
  K <- length(X)
  A <- X[[1L]]
  if (length(A$appended)) stop("objects with appended columns cannot be bound by rows.", call. = FALSE)
  for (k in seq_len(K)[-1L]) {
    B <- X[[k]]
    if (B$p != A$p || !identical(as.character(B$chr), as.character(A$chr)) || !identical(B$pos, A$pos) ||
        !identical(B$a1, A$a1) || !identical(B$a2, A$a2) || length(B$appended))
      stop("object ", k, " does not hold the columns of object 1 in the same order and orientation ",
           "(see geno_harmonise(), geno_subset()).", call. = FALSE)
  }
  n <- vapply(X, function(x) x$n, 1L, USE.NAMES = FALSE)
  iid <- as.character(unlist(lapply(X, `[[`, "sample"), use.names = FALSE))
  fid <- as.character(unlist(lapply(X, function(x) if (is.null(x$fid)) x$sample else x$fid), use.names = FALSE))
  obj <- rep(seq_len(K), n); row <- sequence(n)
  key <- paste(fid, iid, sep = "\r")
  if (is.null(sample)) {
    keep <- which(!duplicated(key))
  } else {
    if (is.character(sample) && length(sample) == 1L) sample <- data.table::fread(sample, colClasses = "character", data.table = FALSE)
    if (!is.data.frame(sample) || ncol(sample) < 2L) stop("sample must be a data.frame (FID, IID) or the path of such a file.", call. = FALSE)
    keep <- match(unique(paste(sample[[1L]], sample[[2L]], sep = "\r")), key)
    keep <- keep[!is.na(keep)]
    if (length(keep) < 2L) stop("fewer than two samples of `sample` are in the objects.", call. = FALSE)
  }
  whole <- length(keep) == length(key) && identical(keep, seq_along(key))
  M <- matrix(seq_len(A$p) - 1L, K, A$p, byrow = TRUE)
  ptr <- geno_bind_cpp(lapply(X, `[[`, "ptr"), M, if (whole) integer(0) else obj[keep] - 1L,
                       if (whole) integer(0) else row[keep] - 1L, .threads(threads))
  same <- function(f) all(vapply(X, function(x) identical(x[[f]], A[[f]]), NA))
  .geno_new(ptr, .geno_cols(A), list(sample = iid[keep], fid = fid[keep]), impute = if (same("impute")) A$impute else "mixed",
            scale = A$scale, threads = A$threads, file = if (same("file")) A$file)
}
