# geno: a PLINK BED or PGEN genotype matrix kept in memory as 2-bit planes
# (n p / 4 bytes) instead of an n x p double matrix (moved here from SuSiE4I,
# which now uses these functions). The object behaves like a numeric matrix
# for dim(), dimnames() and [ , and gives the products fine-mapping needs:
# geno_crossprod() (X'X, cached per object; X'M), geno_multiply() (X B),
# geno_wcrossprod(), geno_append() (extra columns) and geno_cor() /
# bed_cor() / pgen_cor() (correlations). Values are A1 (BED) or ALT (PGEN)
# allele dosages, as in BEDMatrix. Every column is held as X = H + D: H a
# 0/1/2 matrix in 2-bit planes, D a sparse int16 part in the column's own
# unit (1/16384 for genotypes) holding dosage - H for a PGEN with dosages,
# f - round(f) for a missing entry filled with f (the variant's lower median,
# impute = "median", or observed mean, "mean"), and all of x (H = 0) for an
# appended column; with scale = TRUE every product is that of the
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
#'   automatically and read as ALT dosages: values, `center`, `sd` and all
#'   products are those of the dosage matrix, stored as the rounded hard call
#'   in 2-bit planes plus the exact dosage - hard call correction (a 1-bit
#'   nonzero mask and int16 values in 1/16384 units). Such files must be
#'   biallelic; a missing dosage takes the variant's observed mean dosage
#'   (to 1/32768) whatever `impute` is.
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
#' @param fill Externally supplied fills (used by [genotype_fill()]).
#' @param af Optional anchor frequency of the counted allele (bim column 5 for
#'   BED, ALT for PGEN), one value per selected variant; `NA` entries keep
#'   `impute`. With a finite `af`, every missing call of that variant takes
#'   c = (2 n af - sum of observed) / (number missing), clipped to `[0, 2]`
#'   (stored as round(c) in the planes and the rest in D), so the filled
#'   column has mean 2 af and `sd` is the sample sd after filling. `NULL`
#'   (default) leaves the behaviour unchanged.
#' @return An object of class `geno` with fields `n`, `p`, `snp`,
#'   `a1`, `a2` (aligned with `snp`; `a1` is the counted allele: bim column 5
#'   for BED, ALT for PGEN, and `a2` the other allele: bim column 6 or REF),
#'   `sample` (IIDs in row order), `impute`, `scale`, `center`, `sd`, `dosage`
#'   (whether the PGEN holds dosages) and `dosage_nnz` (nonzero entries of the
#'   int16 deviation D). It holds an
#'   external pointer, so it is valid only in the current R session (it cannot
#'   be saved with `saveRDS` or sent to parallel workers).
#' @export
geno_open <- function(bedfile = NULL, pgenfile = NULL, snp_vec = NULL,
                      sample_vec = NULL, impute = c("median", "mean"),
                      scale = TRUE, threads = 4L, fill = NULL, af = NULL) {
  impute <- match.arg(impute)
  threads <- .threads(threads)
  if (is.null(bedfile) == is.null(pgenfile)) {
    stop("Give exactly one of bedfile and pgenfile.", call. = FALSE)
  }
  if (!is.null(bedfile)) {
    bed <- .geno_bed_path(bedfile)
    base <- sub("\\.bed$", "", bed)
    bim <- data.table::fread(paste0(base, ".bim"), header = FALSE, colClasses = "character")
    fam <- data.table::fread(paste0(base, ".fam"), header = FALSE, colClasses = "character")
    snp <- bim[[2L]]
    a1 <- bim[[5L]]
    a2 <- bim[[6L]]
    fid <- fam[[1L]]
    iid <- fam[[2L]]
  } else {
    files <- .geno_pgen_files(pgenfile)
    pv <- data.table::fread(files$pvar, skip = "#CHROM", colClasses = "character")
    snp <- pv[["ID"]]
    a1 <- pv[["ALT"]]
    a2 <- pv[["REF"]]
    allele_ct <- 1L + nchar(pv$ALT) - nchar(gsub(",", "", pv$ALT, fixed = TRUE))
    ps <- data.table::fread(files$psam, colClasses = "character")
    iid <- ps[[if ("#IID" %in% names(ps)) "#IID" else "IID"]]
    fid_col <- intersect(c("#FID", "FID"), names(ps))
    fid <- if (length(fid_col)) ps[[fid_col[1L]]] else rep("0", length(iid))
  }
  vidx <- .geno_index(snp_vec, snp, NULL, "snp_vec")
  sidx <- .geno_index(sample_vec, iid, paste(fid, iid, sep = "_"), "sample_vec")
  fs <- .geno_fill_spec(fill, length(vidx), iid[sidx])
  if (is.null(af)) {
    af <- numeric(0)
  } else {
    af <- as.numeric(af)
    if (length(af) != length(vidx)) stop("af needs one value per selected variant.", call. = FALSE)
    if (any(af < 0 | af > 1, na.rm = TRUE)) stop("af must lie in [0, 1].", call. = FALSE)
  }
  ptr <- if (!is.null(bedfile)) {
    geno_open_bed_cpp(bed, length(iid), length(snp), vidx - 1L, sidx - 1L,
                      impute == "mean", threads, fs$c, fs$ptr, fs$i, fs$g, af)
  } else {
    geno_open_pgen_cpp(files$pgen, length(iid), allele_ct,
                       vidx - 1L, sidx - 1L, impute == "mean", threads, fs$c, fs$ptr, fs$i, fs$g, af)
  }
  info <- geno_info_cpp(ptr)
  structure(list(ptr = ptr, n = length(sidx), p = length(vidx),
                 snp = snp[vidx],
                 a1 = a1[vidx], a2 = a2[vidx], sample = iid[sidx],
                 impute = impute, scale = scale,
                 center = info$mean, sd = info$sd, threads = threads,
                 dosage = isTRUE(info$dosage), dosage_nnz = info$dosage_nnz,
                 appended = integer(0), cache = new.env(parent = emptyenv())),
            class = "geno")
}

#' Append columns to a geno matrix
#'
#' Appends the columns of a dense n x k matrix `M` (one row per sample, in the
#' row order of `X`) after the columns of `X`, in the same store: the 2-bit
#' part of an appended column is 0 and all of its values sit in the int16 part,
#' with the column's own unit max|x| / 32767 (relative error at most
#' 1/65534 of max|x|), so a dense or a sparse, 0/1/2 or continuous column is
#' held the same way. Means and sds are kept per column. `X` is not changed.
#' The raw X'X of `X` is computed once per object and reused, so only X'M and
#' M'M are added.
#'
#' @param X A `geno` object.
#' @param M Numeric matrix (or vector) with one row per sample.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A `geno` object with `p + ncol(M)` columns; `appended` gives the
#'   positions of all appended columns and `a1`, `a2` are `NA` there.
#' @export
geno_append <- function(X, M, threads = X$threads) {
  if (!inherits(X, "geno")) stop("X must be a geno object.", call. = FALSE)
  M <- as.matrix(M) + 0
  if (nrow(M) != X$n) stop("M must have one row per sample.", call. = FALSE)
  if (any(!is.finite(M))) stop("Appended columns must be finite.", call. = FALSE)
  if (is.null(colnames(M))) colnames(M) <- paste0("A", X$p + seq_len(ncol(M)))
  if (anyDuplicated(c(X$snp, colnames(M)))) stop("Appended column names must differ from the column names of X.", call. = FALSE)
  threads <- .threads(threads)
  ptr <- geno_append_cpp(X$ptr, M, threads)
  info <- geno_info_cpp(ptr)
  p0 <- X$p; k <- ncol(M); p <- p0 + k
  out <- X
  out$ptr <- ptr; out$p <- p
  out$snp <- c(X$snp, colnames(M))
  out$a1 <- c(X$a1, rep(NA_character_, k)); out$a2 <- c(X$a2, rep(NA_character_, k))
  out$center <- info$mean; out$sd <- info$sd
  out$dosage_nnz <- info$dosage_nnz
  out$appended <- c(X$appended, p0 + seq_len(k))
  out$cache <- new.env(parent = emptyenv())
  A <- geno_dense_cpp(ptr, seq_len(X$n) - 1L, p0 + seq_len(k) - 1L, threads)  # stored values of the appended columns
  K01 <- geno_xtm_cpp(X$ptr, A, threads)
  K <- matrix(0, p, p)
  K[seq_len(p0), seq_len(p0)] <- .geno_raw_xtx(X, threads)
  K[seq_len(p0), p0 + seq_len(k)] <- K01
  K[p0 + seq_len(k), seq_len(p0)] <- t(K01)
  K[p0 + seq_len(k), p0 + seq_len(k)] <- dense_xtx_cpp(A, threads)
  out$cache$xtx <- K
  out
}

# Raw (unstandardized) X'X of a geno, computed once per object.
.geno_raw_xtx <- function(X, threads = X$threads) {
  if (!is.null(X$cache) && !is.null(X$cache$xtx)) return(X$cache$xtx)
  K <- geno_xtx_cpp(X$ptr, as.integer(threads))
  if (!is.null(X$cache)) X$cache$xtx <- K
  K
}

# fill = list(const = one value per selected variant (NA = built-in fill),
#             cells = data.frame(j = column, IID, g)) -> CSR over the selected rows.
.geno_fill_spec <- function(fill, p, iid) {
  out <- list(c = numeric(0), ptr = integer(0), i = integer(0), g = numeric(0))
  if (is.null(fill)) return(out)
  if (!is.null(fill$const)) {
    if (length(fill$const) != p) stop("fill$const needs one value per selected variant.", call. = FALSE)
    out$c <- as.numeric(fill$const)
  }
  cl <- fill$cells
  if (!is.null(cl) && nrow(cl)) {
    i <- match(as.character(cl$IID), iid)
    j <- as.integer(cl$j)
    ok <- !is.na(i) & !is.na(j)
    o <- order(j[ok], i[ok])
    j <- j[ok][o]; i <- i[ok][o]; g <- as.numeric(cl$g)[ok][o]
    out$ptr <- c(0L, cumsum(tabulate(j, nbins = p)))
    out$i <- i - 1L
    out$g <- g
  }
  out
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

# Scale divisor: constant columns are centered but not divided (as scale()).
.geno_div <- function(X) ifelse(X$sd > 0, X$sd, 1)

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
    return(invisible(x))
  }
  cat(sprintf("geno matrix: %d samples x %d variants, impute = %s, scale = %s (%.2f GB in memory)\n",
              x$n, x$p, x$impute, x$scale,
              sum(unlist(geno_info_cpp(x$ptr)[c("bytes_H", "bytes_D")])) / 1e9))
  invisible(x)
}

#' @export
`[.geno` <- function(x, i, j, drop = TRUE) {
  rows <- if (missing(i)) seq_len(x$n) else seq_len(x$n)[i]
  cols <- if (missing(j)) seq_len(x$p) else if (is.character(j)) match(j, x$snp) else seq_len(x$p)[j]
  if (anyNA(rows) || anyNA(cols)) stop("subscript out of bounds", call. = FALSE)
  G <- geno_dense_cpp(x$ptr, rows - 1L, cols - 1L, x$threads)
  if (x$scale) {
    G <- sweep(sweep(G, 2L, x$center[cols], "-"), 2L, .geno_div(x)[cols], "/")
  }
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
#' 1/16384); missing hardcalls are filled by `impute`. With `scale = TRUE`
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
#' `geno_crossprod(X)` is `crossprod(X)` by exact popcount; `geno_crossprod(X, M)`
#' is `crossprod(X, M)` for a dense `M` with one row per sample. With
#' `X$scale = TRUE` both are those of the column-standardized matrix.
#'
#' @param X A `geno` object from [geno_open()].
#' @param M Optional dense matrix (or vector) with one row per sample.
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A p by p matrix when `M` is `NULL`, otherwise a p by `ncol(M)` matrix.
#' @export
geno_crossprod <- function(X, M = NULL, threads = X$threads) {
  if (is.null(M)) {
    out <- .geno_raw_xtx(X, threads)
    if (X$scale) {
      d <- .geno_div(X)
      out <- geno_std_cpp(out, X$center, X$center, X$n, d, d, .threads(threads))
    }
    dimnames(out) <- list(X$snp, X$snp)
    return(out)
  }
  M <- as.matrix(M) + 0
  out <- geno_xtm_cpp(X$ptr, M, threads)
  if (X$scale) out <- geno_std_cpp(out, X$center, colSums(M), 1, .geno_div(X), numeric(0), .threads(threads))
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
  if (X$scale) {
    B <- B / .geno_div(X)
    out <- geno_mx_cpp(X$ptr, B, threads)
    out <- sweep(out, 2L, colSums(B * X$center), "-")
  } else {
    out <- geno_mx_cpp(X$ptr, B, threads)
  }
  if (vec) as.numeric(out) else out
}

# crossprod(X, X * w) and crossprod(X, M) (of the standardized matrix when X$scale), all in C++:
# X'WX by 256 x 256 tiles over blocks of 4096 samples (float products, double sums; geno_xtwx_cpp), X'M by
# geno_crossprod(). block_size is no longer used (kept so existing calls still work).
geno_wcrossprod <- function(X, w, M = NULL, block_size = 10000L, threads = X$threads) {
  threads <- .threads(threads)
  w <- as.numeric(w)
  if (length(w) != X$n) stop("w must have one value per sample.", call. = FALSE)
  XtWX <- if (X$scale) geno_xtwx_cpp(X$ptr, w, X$center, .geno_div(X), threads) else
    geno_xtwx_cpp(X$ptr, w, numeric(0), numeric(0), threads)
  XtM <- if (is.null(M)) matrix(0, X$p, 0L) else geno_crossprod(X, M, threads = threads)
  dimnames(XtWX) <- list(X$snp, X$snp)
  rownames(XtM) <- X$snp
  list(XtWX = XtWX, XtM = XtM)
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
  pvar <- paste0(prefix, ".pvar")
  psam <- paste0(prefix, ".psam")
  for (path in c(pgen, pvar, psam)) {
    if (!file.exists(path)) stop("Required PGEN companion file is missing: ",
                                 path, call. = FALSE)
  }
  list(pgen = normalizePath(pgen, mustWork = TRUE),
       pvar = normalizePath(pvar, mustWork = TRUE),
       psam = normalizePath(psam, mustWork = TRUE))
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
  n <- A$n
  if (is.null(B)) {
    K <- .geno_raw_xtx(A, threads)
    B <- A
  } else {
    if (!inherits(B, "geno") || B$n != n || !identical(A$sample, B$sample)) {
      stop("A and B must be geno objects with the same samples in the same order.", call. = FALSE)
    }
    if (A$dosage_nnz == 0 && B$dosage_nnz == 0) {
      K <- geno_cross_cpp(A$ptr, B$ptr, threads)
    } else {
      K <- matrix(0, A$p, B$p)
      rows <- seq_len(n) - 1L
      for (s in seq.int(1L, B$p, by = 64L)) {
        cols <- s:min(B$p, s + 63L)
        K[, cols] <- geno_xtm_cpp(A$ptr, geno_dense_cpp(B$ptr, rows, cols - 1L, threads), threads)
      }
    }
  }
  R <- geno_std_cpp(K, A$center, B$center, n, (n - 1) * A$sd, B$sd, threads)
  R[A$sd == 0, ] <- NA_real_
  R[, B$sd == 0] <- NA_real_
  R <- pmin(pmax(R, -1), 1)
  dimnames(R) <- list(A$snp, B$snp)
  R
}

#' Correlations of all variants of PLINK files
#'
#' `bed_cor(A)` / `pgen_cor(A)`: correlations of all variants of a BED / PGEN
#' file; with `B`, of the variants of A against those of B (same samples in
#' the same order). Missing calls take the variant's lower median; a PGEN with
#' dosages is read as dosages. Moved here from CppMatrix (same arguments); the
#' genotypes are held by [geno_open()] and the products are those of
#' [geno_cor()]. These are thin wrappers: `bed_cor(A)` is
#' `cor_bed(A, impute = "median")$R` (and likewise for `pgen_cor`); use
#' [cor_pgen()] / [cor_bed()] to select variants or samples or to anchor the fill.
#'
#' @param A,B BED paths or prefixes (`bed_cor`), PGEN paths or prefixes
#'   (`pgen_cor`).
#' @param threads Number of OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return A correlation matrix with variant IDs as dimnames.
#' @export
bed_cor <- function(A, B = NULL, threads = 4L) {
  .geno_file_cor(A, B, threads, "bedfile")
}

#' @rdname bed_cor
#' @export
pgen_cor <- function(A, B = NULL, threads = 4L) {
  .geno_file_cor(A, B, threads, "pgenfile")
}

.geno_file_cor <- function(A, B, threads, arg) {
  threads <- .threads(threads)
  if (is.null(B)) {
    file <- stats::setNames(list(A), arg)
    return(.cor_file(file, NULL, NULL, NULL, "median", threads, max_cells = 1e8)$R)
  }
  open1 <- function(f) do.call(geno_open, stats::setNames(list(f, "median", FALSE, as.integer(threads)),
                                                          c(arg, "impute", "scale", "threads")))
  XA <- open1(A)
  XB <- open1(B)
  if (as.numeric(XA$p) * XB$p > 1e8) {
    stop("Requested correlation matrix is too large; subset variants first.", call. = FALSE)
  }
  geno_cor(XA, XB, threads)
}
