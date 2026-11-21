# geno: a PLINK BED or PGEN genotype matrix kept in memory as 2-bit planes
# (n p / 4 bytes) instead of an n x p double matrix. The object behaves like a numeric matrix
# for dim(), dimnames() and [ , and gives the products fine-mapping needs:
# geno_crossprod() (X'X; X'M), geno_multiply() (X B),
# geno_wcrossprod() and geno_cor() /
# cor_bed() / cor_pgen() (correlations). Values are A1 (BED) or ALT (PGEN)
# allele dosages, as in BEDMatrix. Every column is held as X = H + D: H a
# 0/1/2 matrix in 2-bit planes, D a sparse int16 part in the column's own
# unit (1/16384 for hard-call inputs, 0.01 for a PGEN with dosages) holding dosage - H for a PGEN with dosages,
# f - round(f) for a missing entry filled with f (the variant's lower median,
# impute = "median", or observed mean, "mean"); with scale = TRUE every product is that of the
# column-standardized matrix (mean 0, sd 1 with n - 1; constant columns are
# only centered).

#' Open a BED or PGEN file as a compact genotype matrix
#'
#' Keeps the genotypes of a PLINK BED or PGEN file in 2-bit form, so R
#' never holds an n by p double matrix in R. Values are A1 (BED) or ALT (PGEN)
#' allele counts as in BEDMatrix. Multiallelic PGEN variants are counted as REF
#' versus all non-REF alleles. To choose the counted direction, set the target
#' allele as REF upstream (`plink2 --ref-allele force target.txt 1 2
#' --make-pgen`); the counted allele is then everything else.
#'
#' @param bedfile PLINK 1 BED path or prefix (`.bim` and `.fam` alongside).
#' @param pgenfile PLINK 2 PGEN path or prefix (uncompressed `.pvar` and
#'   `.psam` alongside). Give either `bedfile` or `pgenfile`. A PGEN whose
#'   header flags dosages (e.g. `--make-pgen` from imputed data) is detected
#'   automatically and read as ALT dosages taken to the nearest 0.01 (the
#'   precision of imputed dosages; the file holds 1/16384 steps): values,
#'   `center`, `sd` and all products are those of that dosage matrix, stored as
#'   the rounded hard call in 2-bit planes plus the dosage - hard call
#'   correction (a 1-bit nonzero mask and int16 values in 0.01 units). Such
#'   files must be biallelic; a missing dosage takes the variant's observed
#'   mean dosage (to 0.01) whatever `impute` is.
#' @param snp_vec Variants to keep, as IDs or 1-based file indices, in the
#'   column order wanted. `NULL` keeps all in file order. Selection by ID takes
#'   the first match, so use indices when IDs are duplicated (e.g. ".").
#' @param sample_vec Samples to keep, as IIDs (or `FID_IID`) or 1-based file
#'   indices, in the row order wanted. `NULL` keeps all in file order. The rows
#'   of `y`, `Z` and `status` passed to `SuSiE4I()` must be in this order (or in
#'   `.fam`/`.psam` order when `NULL`). IIDs repeated across families need the
#'   `FID_IID` form.
#' @param impute `"median"` (default, the variant's lower median) or `"mean"`
#'   (the variant's observed mean, stored to 1/32768: its rounded value in the
#'   2-bit planes, the rest as an int16 deviation at the missing entries).
#'   Hardcall input only.
#' @param scale Whether products and dense blocks are of the standardized
#'   matrix.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @param af Optional anchor frequency of the counted allele (bim column 5 for
#'   BED, ALT for PGEN), one value per selected variant; `NA` entries keep
#'   `impute`. With a finite `af`, every missing call of that variant takes
#'   c = (2 n af - sum of observed) / (number missing), clipped to `[0, 2]`
#'   (stored as round(c) in the planes and the rest in D), so the filled
#'   column has mean 2 af and `sd` is the sample sd after filling: the same as
#'   [genotype_fill()] with `af` on the opened object. `NULL` (default): no
#'   anchor.
#' @param file A [geno_file()] result of the same file, used instead of reading the variant and sample
#'   tables (then leave `bedfile` and `pgenfile` `NULL`).
#' @return An object of class `geno` with fields `n`, `p`, `snp`,
#'   `a1`, `a2` (aligned with `snp`; `a1` is the counted allele: bim column 5
#'   for BED, ALT for PGEN, and `a2` the other allele: bim column 6 or REF),
#'   `sample` (IIDs in row order), `impute` (how the missing cells are filled now: "median", "mean", "af",
#'   "ld" or "custom"), `scale`, `center`, `sd`, `dosage`
#'   (whether the PGEN holds dosages) and `dosage_nnz` (nonzero entries of the
#'   int16 deviation D); per column also `id` (file ID; `chr:pos:a1:a2` when opened from
#'   `geno_file(id = FALSE)`, built once at open for the opened rows only), `chr`, `pos`, `row` (file row), `flipped` (stored
#'   value = 2 - file value), `n_miss` (missing calls) and `af` (observed frequency of `a1`); per row `fid`;
#'   and `file` = list(type, path). The rows of the missing cells are kept (see [geno_missing()]). It holds an
#'   external pointer, so it is valid only in the current R session (it cannot
#'   be saved with `saveRDS` or sent to parallel workers).
#' @export
geno_open <- function(bedfile = NULL, pgenfile = NULL, snp_vec = NULL,
                      sample_vec = NULL, impute = c("median", "mean"),
                      scale = TRUE, threads = 4L, af = NULL, file = NULL) {
  impute <- match.arg(impute)
  threads <- .threads(threads)
  gf <- .as_geno_file(file, bedfile, pgenfile)
  v <- gf$variants
  s <- gf$samples
  # without IDs, chr:pos:a1:a2 names only of the opened rows (and of the candidate rows of character snp_vec)
  ids <- if ("id" %in% names(v)) v$id else NULL
  vidx <- if (is.null(ids) && is.character(snp_vec)) .geno_vid_index(v, snp_vec) else
    .geno_index(snp_vec, if (is.null(ids)) v$row else ids, NULL, "snp_vec")
  sidx <- .geno_index(sample_vec, s$iid, paste(s$fid, s$iid, sep = "_"), "sample_vec")
  ptr <- if (gf$type == "bed") {
    geno_open_bed_cpp(gf$path, nrow(s), nrow(v), vidx - 1L, sidx - 1L, impute == "mean", threads)
  } else {
    geno_open_pgen_cpp(gf$path, nrow(s), v$allele_ct, vidx - 1L, sidx - 1L, impute == "mean", threads)
  }
  # names built once here, for the opened rows only, and kept in the object (snp and id)
  nm <- if (is.null(ids)) .geno_vid(v[vidx]) else ids[vidx]
  X <- .geno_new(ptr, list(snp = nm, id = nm,
                           chr = v$chr[vidx], pos = v$pos[vidx], a1 = v$a1[vidx], a2 = v$a2[vidx], row = vidx,
                           flipped = rep(FALSE, length(vidx))),
                 list(sample = s$iid[sidx], fid = s$fid[sidx]), impute = impute, scale = scale, threads = threads,
                 file = list(type = gf$type, path = gf$path))
  if (!is.null(af))  # on the new object, which nothing else references yet: in place
    X <- .geno_refill(X, af = af, default = "keep", inplace = TRUE, impute = "af", threads = threads)
  X
}

# The one constructor of a geno object: per-column fields `cols` (snp, id, chr, pos, a1, a2, row, flipped),
# per-row fields `rows` (sample, fid); center, sd, n_miss and af come from the store.
.geno_new <- function(ptr, cols, rows, impute, scale, threads, file = NULL) {
  info <- geno_info_cpp(ptr)
  n <- as.integer(info$n); p <- as.integer(info$p)
  for (k in names(cols)) if (length(cols[[k]]) != p) stop("internal: column field ", k, " has the wrong length", call. = FALSE)
  for (k in names(rows)) if (length(rows[[k]]) != n) stop("internal: row field ", k, " has the wrong length", call. = FALSE)
  structure(list(ptr = ptr, n = n, p = p, snp = cols$snp, a1 = cols$a1, a2 = cols$a2, sample = rows$sample,
                 impute = impute, scale = scale, center = info$mean, sd = info$sd, threads = threads,
                 dosage = isTRUE(info$dosage), dosage_nnz = info$dosage_nnz,
                 id = cols$id, chr = cols$chr, pos = cols$pos, row = cols$row, flipped = cols$flipped,
                 n_miss = info$n_miss, af = info$af, fid = rows$fid, file = file),
            class = "geno")
}

.geno_cols <- function(X, j = seq_len(X$p)) {
  list(snp = X$snp[j], id = X$id[j], chr = X$chr[j], pos = X$pos[j], a1 = X$a1[j], a2 = X$a2[j],
       row = X$row[j], flipped = X$flipped[j])
}

.geno_index <- function(x, ids, alt_ids, arg) {
  if (is.null(x)) return(seq_along(ids))
  if (is.character(x)) {
    idx <- match(x, ids)
    if (!is.null(alt_ids) && anyNA(idx)) idx[is.na(idx)] <- match(x[is.na(idx)], alt_ids)
    if (anyNA(idx)) {
      stop(sum(is.na(idx)), " entries of ", arg, " are not in the file, e.g. ",
           x[which(is.na(idx))[1L]], call. = FALSE)
    }
    return(idx)
  }
  idx <- as.integer(x)
  if (anyNA(idx) || any(idx < 1L | idx > length(ids))) {
    stop(arg, " indices must lie between 1 and ", length(ids), ".", call. = FALSE)
  }
  idx
}

#' @export
dim.geno <- function(x) c(x$n, x$p)

#' @export
dimnames.geno <- function(x) list(NULL, x$snp)

#' @export
print.geno <- function(x, ...) {
  if (isTRUE(x$dosage)) {
    cat(sprintf("geno matrix: %d samples x %d variants, dosage (%.1f%% of entries off the hard call), scale = %s (%.2f GB in memory)\n",
                x$n, x$p, 100 * x$dosage_nnz / (x$n * x$p), x$scale,
                sum(unlist(geno_info_cpp(x$ptr)[c("bytes_H", "bytes_D")])) / 1e9))
    if (!is.null(x$n_miss)) cat(sprintf("missing cells: %.0f in %d variants\n", sum(x$n_miss), sum(x$n_miss > 0)))
    return(invisible(x))
  }
  cat(sprintf("geno matrix: %d samples x %d variants, impute = %s, scale = %s (%.2f GB in memory)\n",
              x$n, x$p, x$impute, x$scale,
              sum(unlist(geno_info_cpp(x$ptr)[c("bytes_H", "bytes_D")])) / 1e9))
  if (!is.null(x$n_miss)) cat(sprintf("missing cells: %.0f in %d variants\n", sum(x$n_miss), sum(x$n_miss > 0)))
  invisible(x)
}

#' @export
`[.geno` <- function(x, i, j, drop = TRUE) {
  rows <- if (missing(i)) seq_len(x$n) else seq_len(x$n)[i]
  cols <- if (missing(j)) seq_len(x$p) else if (is.character(j)) match(j, x$snp) else seq_len(x$p)[j]
  if (anyNA(rows) || anyNA(cols)) stop("subscript out of bounds", call. = FALSE)
  G <- geno_dense_cpp(x$ptr, rows - 1L, cols - 1L, x$scale, x$threads)
  colnames(G) <- x$snp[cols]
  if (drop) drop(G) else G
}

#' @export
as.matrix.geno <- function(x, ...) x[, , drop = FALSE]

#' Read a BED or PGEN file into a dense numeric matrix
#'
#' Same reader and arguments as [geno_open()], returned as an ordinary n x p
#' double matrix (rows = samples, named by IID; columns = variants). Values are
#' A1 (BED) or ALT (PGEN) allele counts, or dosages for a dosage PGEN (to
#' 0.01); missing hardcalls are filled by `impute`. With `scale = TRUE`
#' every column is standardized (mean 0, sd 1 with n - 1). For a large region
#' prefer [geno_open()], which never forms this matrix.
#'
#' @inheritParams geno_open
#' @return A numeric matrix with attributes `a1` and `a2` (counted and other
#'   allele per column).
#' @export
geno_read <- function(bedfile = NULL, pgenfile = NULL, snp_vec = NULL, sample_vec = NULL,
                      impute = c("median", "mean"), scale = FALSE, threads = 4L) {
  X <- geno_open(bedfile = bedfile, pgenfile = pgenfile, snp_vec = snp_vec, sample_vec = sample_vec,
                 impute = impute, scale = scale, threads = threads)
  G <- X[, , drop = FALSE]
  rownames(G) <- X$sample
  attr(G, "a1") <- X$a1
  attr(G, "a2") <- X$a2
  G
}

#' Cross products of a geno matrix
#'
#' `geno_crossprod(X)` is `crossprod(X)` by exact popcount, or by products of
#' decoded row blocks when more than 7% of the cells have an int16 part
#' (dosages; exact for the 0.01 dosages of a PGEN); `geno_crossprod(X, M)` is `crossprod(X, M)` for a dense `M`
#' with one row per sample. With `X$scale = TRUE` both are those of the
#' column-standardized matrix. Every call computes its product; the object keeps none, so keep the result when
#' it is needed again.
#'
#' @param X A `geno` object from [geno_open()].
#' @param M Optional dense matrix (or vector) with one row per sample.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A p by p matrix when `M` is `NULL`, otherwise a p by `ncol(M)` matrix.
#' @export
geno_crossprod <- function(X, M = NULL, threads = X$threads) {
  if (is.null(M)) {
    out <- geno_xtx_cpp(X$ptr, X$scale, threads)
    dimnames(out) <- list(X$snp, X$snp)
    return(out)
  }
  M <- as.matrix(M) + 0
  out <- geno_xtm_cpp(X$ptr, M, X$scale, threads)
  rownames(out) <- X$snp
  out
}

#' Multiply a geno matrix by a dense matrix or vector
#'
#' `geno_multiply(X, B)` is `X %*% B` for `B` with one row per variant; a
#' vector `B` gives a vector. With `X$scale = TRUE` it is the product of the
#' column-standardized matrix.
#'
#' @param X A `geno` object from [geno_open()].
#' @param B Dense matrix with one row per variant, or a numeric vector.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return An n by `ncol(B)` matrix, or a numeric vector when `B` is a vector.
#' @export
geno_multiply <- function(X, B, threads = X$threads) {
  vec <- is.null(dim(B))
  B <- as.matrix(B) + 0
  out <- geno_mx_cpp(X$ptr, B, X$scale, threads)
  if (vec) as.numeric(out) else out
}

# crossprod(X, X * w) and crossprod(X, M) (of the standardized matrix when X$scale), all in C++:
# X'WX by 256 x 256 tiles over blocks of 4096 samples (float products, double sums; geno_xtwx_cpp), X'M by
# geno_crossprod().
geno_wcrossprod <- function(X, w, M = NULL, threads = X$threads) {
  w <- as.numeric(w)
  if (length(w) != X$n) stop("w must have one value per sample.", call. = FALSE)
  XtWX <- geno_xtwx_cpp(X$ptr, w, X$scale, threads)
  XtM <- if (is.null(M)) matrix(0, X$p, 0L) else geno_crossprod(X, M, threads = threads)
  dimnames(XtWX) <- list(X$snp, X$snp)
  rownames(XtM) <- X$snp
  list(XtWX = XtWX, XtM = XtM)
}

# crossprod(X, X * wq) (of the standardized matrix when X$scale, centered with wq), wq = w binned: rows sorted
# by w, each run of 64 sorted rows shares the mean of its w, so X'diag(wq)X sums exact popcount products per
# bin times the bin weight (geno_xtwx_bin_cpp). When more than 1% of the cells have a nonzero D (dosages,
# mean-filled missing calls) the popcounts do not pay: wq = w and the product is
# geno_wcrossprod's. Returns the matrix and wq, the weights it corresponds to.
geno_wcrossprod_bin <- function(X, w, threads = X$threads) {
  w <- as.numeric(w)
  if (length(w) != X$n) stop("w must have one value per sample.", call. = FALSE)
  out <- geno_xtwx_bin_cpp(X$ptr, w, X$scale, threads)
  dimnames(out$XtWX) <- list(X$snp, X$snp)
  list(XtWX = out$XtWX, w = out$w)
}

.geno_bed_path <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("bedfile must be one BED file path.", call. = FALSE)
  }
  if (!grepl("\\.bed$", x)) x <- paste0(x, ".bed")
  if (!file.exists(x)) stop("BED file not found: ", x, call. = FALSE)
  normalizePath(x, winslash = "/")
}

.geno_pgen_files <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("pgenfile must be a single PGEN file path or prefix.", call. = FALSE)
  }
  prefix <- sub("\\.pgen$", "", x, ignore.case = TRUE)
  pgen <- paste0(prefix, ".pgen")
  for (path in c(pgen, paste0(prefix, ".psam"))) {
    if (!file.exists(path)) stop("Required PGEN companion file is missing: ",
                                 path, call. = FALSE)
  }
  list(pgen = normalizePath(pgen, mustWork = TRUE), prefix = prefix)
}

#' Correlations of the columns of geno matrices
#'
#' `geno_cor(A)` is `cor(A)`, `geno_cor(A, B)` is `cor(A, B)` for two geno
#' objects over the same samples in the same order, from the raw cross
#' products (exact popcount when neither has an int16 part). Constant columns
#' give `NA`.
#'
#' @param A,B `geno` objects.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A `ncol(A)` by `ncol(B)` matrix (B = A when `NULL`).
#' @export
geno_cor <- function(A, B = NULL, threads = A$threads) {
  threads <- .threads(threads)
  if (!is.null(B) && (!inherits(B, "geno") || B$n != A$n || !identical(A$sample, B$sample)))
    stop("A and B must be geno objects with the same samples in the same order.", call. = FALSE)
  R <- geno_cor_cpp(A$ptr, if (is.null(B)) NULL else B$ptr, threads)
  dimnames(R) <- list(A$snp, if (is.null(B)) A$snp else B$snp)
  R
}

