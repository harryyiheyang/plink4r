# Merging geno objects: one tool for binding columns, binding rows, and objects that differ in both.
# The unit is the cell (sample x variant), as in plink2 --pmerge: the result holds the samples and the
# variants of all objects, and a cell no object holds is a missing cell. The per-column and per-row fields
# (the object's stand-ins for the .bim / .pvar and .fam / .psam) follow.

#' Merge geno objects
#'
#' The result holds the variants of all objects and, by default, the samples of all objects:
#' * variants: the union over the objects of (chr, pos, a1, a2), ordered by chromosome code (1..22, 23 = X,
#'   24 = Y, 25 = XY (pseudo-autosomal), 26 = MT; variants without a code last), position, `a1`, then `a2`, whatever the order of the
#'   objects. The same variant in the opposite orientation is another variant: bring the objects to one
#'   orientation first with [geno_harmonise()].
#' * samples: `sample = NULL`, the union over the objects of (FID, IID) in the order of first appearance;
#'   `sample` given, the samples of `sample` that some object holds, in the order of `sample`.
#' * cells: a cell (sample, variant) is copied as stored (a missing cell keeps its fill) from the first
#'   object, in the order given, that holds both the sample and the variant. A cell no object holds becomes
#'   a missing cell filled with the observed mean of its column; [genotype_fill()] refills the missing cells.
#'
#' So objects with the same samples are bound by columns, objects with the same variants by rows, and a
#' sample or variant held twice is taken from the first object. All objects hold hard calls, or all dosages.
#'
#' @param ... `geno` objects, or one list of them (e.g. the result of [geno_harmonise()]).
#' @param sample `NULL`, or the samples to keep and their order: a data.frame whose first two columns are FID
#'   and IID, or the path of such a file (.fam, .psam, or any table with FID and IID first).
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A `geno` object. `scale` is that of the first object, `impute` and `file` those of the objects
#'   when they agree (else `"mixed"` and `NULL`), the column fields of a variant those of the first object
#'   holding it, and `uncovered` the number of cells per column that no object held.
#' @export
geno_merge <- function(..., sample = NULL, threads = 4L) {
  X <- list(...)
  if (length(X) == 1L && !inherits(X[[1L]], "geno") && is.list(X[[1L]])) X <- X[[1L]]
  if (!length(X) || !all(vapply(X, inherits, NA, "geno"))) stop("give geno objects, or a list of them.", call. = FALSE)
  K <- length(X)
  skey <- lapply(X, function(x) paste(if (is.null(x$fid)) x$sample else x$fid, x$sample, sep = "\r"))
  vkey <- lapply(X, function(x) paste(x$chr, x$pos, x$a1, x$a2, sep = "\r"))
  # samples of the result
  s_obj <- rep(seq_len(K), lengths(skey)); s_all <- unlist(skey, use.names = FALSE)
  if (is.null(sample)) {
    s_first <- which(!duplicated(s_all))
  } else {
    if (is.character(sample) && length(sample) == 1L) sample <- data.table::fread(sample, colClasses = "character", data.table = FALSE)
    if (!is.data.frame(sample) || ncol(sample) < 2L) stop("sample must be a data.frame (FID, IID) or the path of such a file.", call. = FALSE)
    s_first <- match(unique(paste(sample[[1L]], sample[[2L]], sep = "\r")), s_all)
    s_first <- s_first[!is.na(s_first)]
  }
  if (length(s_first) < 2L) stop("fewer than two samples to merge.", call. = FALSE)
  # variants of the result, in chr / pos / a1 / a2 order
  v_obj <- rep(seq_len(K), lengths(vkey)); v_all <- unlist(vkey, use.names = FALSE)
  v_col <- sequence(lengths(vkey))
  v_first <- which(!duplicated(v_all))
  f <- function(name, na) unlist(lapply(X, function(x) if (is.null(x[[name]])) rep(na, x$p) else x[[name]]), use.names = FALSE)[v_first]
  cols <- list(snp = f("snp", NA_character_), id = f("id", NA_character_), chr = f("chr", NA_integer_), pos = f("pos", NA_integer_),
               a1 = f("a1", NA_character_), a2 = f("a2", NA_character_), row = f("row", NA_integer_), flipped = f("flipped", NA))
  o <- order(cols$chr, cols$pos, cols$a1, cols$a2, method = "radix")
  v_first <- v_first[o]
  C <- matrix(-1L, K, length(v_first)); R <- matrix(-1L, K, length(s_first))
  for (k in seq_len(K)) {
    C[k, ] <- match(v_all[v_first], vkey[[k]], nomatch = 0L) - 1L
    R[k, ] <- match(s_all[s_first], skey[[k]], nomatch = 0L) - 1L
  }
  res <- geno_merge_cpp(lapply(X, `[[`, "ptr"), C, R, .threads(threads))
  same <- function(name) all(vapply(X, function(x) identical(x[[name]], X[[1L]][[name]]), NA))
  first_s <- function(name, alt) unlist(lapply(X, function(x) if (is.null(x[[name]])) x[[alt]] else x[[name]]), use.names = FALSE)[s_first]
  out <- .geno_new(res$ptr, lapply(cols, `[`, o), list(sample = first_s("sample", "sample"), fid = first_s("fid", "sample")),
                   impute = if (same("impute")) X[[1L]]$impute else "mixed", scale = X[[1L]]$scale, threads = X[[1L]]$threads,
                   file = if (same("file")) X[[1L]]$file)
  out$uncovered <- res$uncovered
  out
}
