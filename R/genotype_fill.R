#' Fill the missing genotype cells of a geno object
#'
#' Returns a new object whose missing cells (the cells missing in the file, see [geno_missing()]) are
#' rewritten; the observed cells and `x` itself are unchanged. For a column j with m_j missing cells,
#' n rows and observed sum S_j (counted allele = `x$a1`):
#' - default: f_j = S_j / (n - m_j), the observed mean of `x`'s own rows;
#' - `af` (finite entries): f_j = clip((2 n af_j - S_j) / m_j, 0, 2), so the column mean is 2 af_j;
#' - `const` (finite entries): f_j = const_j;
#' - `ld`: per-cell conditional means from the eigen.bin panel ([genotype_impute()] on `x`, anchored to `af`
#'   when given, else to the panel AF); columns absent from the panel keep f_j;
#' - `cells`: the given value of that cell.
#' Priority cells > ld > const > af > mean.
#'
#' `snps` restricts the call to a subset of columns: only their missing cells are rewritten (every other column
#' keeps its current fill, and its `af`, `const` and `cells` entries are ignored), and the LD fill conditions on
#' their observed calls alone (the `cond` of [genotype_impute()]), so columns outside `snps`, e.g. wrong calls
#' that are being replaced, carry no information. Nothing is copied for the subset: the kernel reads the
#' selected columns of `x` in place. A value f is stored as round(f) in the 2-bit planes and the
#' rest in the int16 part (1/16384; 0.01 for a dosage PGEN). All values refer to the counted allele `x$a1` (after [geno_harmonise()],
#' the harmonised one). Hard calls and dosages are taken alike: observed dosages condition the LD fill and
#' missing dosages are filled.
#'
#' @param x A `geno` object.
#' @param af Anchor frequency of the counted allele, one per column (`NA` = not anchored).
#' @param ld NULL, or list(dir = folder of block<N>.eigen.bin and the panel SNP .parquet, ref = optional
#'   parquet path or data.frame in place of that table).
#' @param const Fill constant per column (`NA` = not given).
#' @param snps NULL (every column), or the columns (names, or numbers) to fill and, for `ld`, to condition on.
#' @param cells data.frame(snp = column name (or j = column number), IID, g) of per-cell fills; cells that
#'   are not missing are ignored.
#' @param para named list overriding [genotype_fill_para()].
#' @return A `geno` object (new external pointer), `impute` = "mean", "af", "custom" or "ld".
#' @export
genotype_fill <- function(x, af = NULL, ld = NULL, const = NULL, cells = NULL, snps = NULL, para = list()) {
  if (!inherits(x, "geno")) stop("x must be a geno object.", call. = FALSE)
  P <- genotype_fill_para()
  bad <- setdiff(names(para), names(P))
  if (length(bad)) stop("unknown para: ", paste(bad, collapse = ", "), call. = FALSE)
  P[names(para)] <- para
  P$threads <- .threads(P$threads)
  p <- x$p
  sel <- rep(TRUE, p)
  if (!is.null(snps)) {
    js <- if (is.numeric(snps)) as.integer(snps) else match(as.character(snps), x$snp)
    bad <- is.na(js) | js < 1L | js > p
    if (any(bad)) stop(sum(bad), " of `snps` are not columns of x.", call. = FALSE)
    sel <- seq_len(p) %in% js
  }
  if (!is.null(const)) {
    const <- as.numeric(const)
    if (length(const) != p) stop("const needs one value per column.", call. = FALSE)
  }
  cl <- NULL
  if (!is.null(cells) && nrow(cells)) {
    j <- if (!is.null(cells$j)) as.integer(cells$j) else match(as.character(cells$snp), x$snp)
    cl <- data.frame(j = j, IID = as.character(cells$IID), g = as.numeric(cells$g))
    cl <- cl[!is.na(cl$j) & sel[pmax(cl$j, 1L)], , drop = FALSE]
    if (!nrow(cl)) cl <- NULL
  }
  lab <- if (!is.null(cl)) "custom" else if (!is.null(ld)) "ld" else if (any(is.finite(const))) "custom" else
    if (any(is.finite(af))) "af" else "mean"
  ld_write <- NULL
  if (!is.null(ld) && any(x$n_miss > 0)) {
    P$mem_gb <- .mem_gb(P$mem_gb, P$threads)
    .check_delta(P$delta)
    pan <- .ld_panel(ld$ref, ld$dir)
    al <- .ld_align(pan, .geno_vartab(x))
    use <- sort(intersect(al$vi + 1L, which(sel & x$n_miss > 0 & !(seq_len(p) %in% x$appended))))
    if (length(use)) {
      a <- if (!is.null(af)) { ua <- which(is.finite(af)); stats::setNames(af[ua], x$snp[ua]) }
      cond <- if (!all(sel)) x$snp[sel] else P$cond
      # the kernel writes the LD fills straight into the new object (no per-cell table in R)
      ld_write <- function(ptr) .genotype_impute(x, pan, ld$dir, NULL, x$snp[use], P$flank, NULL, P$delta,
                                                 P$anchor_min_miss, P$thresh, P$threads, NULL, P$mem_gb, a, cond,
                                                 write_ptr = ptr, dosage_step = P$dosage_step, aligned = al)
    }
  }
  .geno_refill(x, af = af, const = const, cells = .geno_cells_csr(cl, x),
               default = "mean", inplace = FALSE, impute = lab, threads = P$threads, cols = if (!all(sel)) sel,
               ld_write = ld_write)
}

# Refill the missing cells of X through geno_fill_cpp, which computes the per-column value f: default ("mean":
# observed mean, 0 when every cell is missing; "keep": the current fill), then af (anchored), then const; cells
# (CSR from .geno_cells_csr) win over all. inplace only for an object nothing else references. ld_write(ptr), when given,
# writes the LD fills into the new store after the per-column values and before the cells (so cells > LD > f).
.geno_refill <- function(X, af = NULL, const = NULL, cells = NULL, default = c("mean", "keep"), inplace = FALSE,
                         impute = X$impute, threads = X$threads, cols = NULL, ld_write = NULL) {
  default <- match.arg(default)
  p <- X$p
  if (!is.null(af)) {
    af <- as.numeric(af)
    if (length(af) != p) stop("af needs one value per selected variant.", call. = FALSE)
    if (any(af < 0 | af > 1, na.rm = TRUE)) stop("af must lie in [0, 1].", call. = FALSE)
  }
  if (length(const)) {
    const <- as.numeric(const)
    if (length(const) != p) stop("const needs one value per selected variant.", call. = FALSE)
  }
  um <- default == "mean"
  af <- if (is.null(af)) numeric(0) else af
  const <- if (is.null(const)) numeric(0) else const
  cols <- if (is.null(cols)) logical(0) else as.logical(cols)
  if (is.null(cells) || !length(cells$ptr)) cells <- list(ptr = integer(0), i = integer(0), g = numeric(0))
  th <- .threads(threads)
  if (is.null(ld_write)) {
    ptr <- geno_fill_cpp(X$ptr, um, af, const, cols, cells$ptr, cells$i, cells$g, inplace, th)
  } else {
    ptr <- geno_fill_cpp(X$ptr, um, af, const, cols, integer(0), integer(0), numeric(0), inplace, th)
    ld_write(ptr)
    if (length(cells$ptr))
      ptr <- geno_fill_cpp(ptr, FALSE, numeric(0), numeric(0), logical(0), cells$ptr, cells$i, cells$g, TRUE, th)
  }
  .geno_new(ptr, .geno_cols(X), list(sample = X$sample, fid = X$fid), impute = impute, scale = X$scale,
            threads = X$threads, file = X$file, appended = X$appended)
}

# Per-cell fills cl = data.frame(j = column, IID, g) as CSR over the rows of X: cells matched by IID, unmatched
# cells dropped, rows ascending within a column (geno_fill_cpp takes the first of repeated cells).
.geno_cells_csr <- function(cl, X) {
  if (is.null(cl) || !nrow(cl)) return(NULL)
  i <- match(as.character(cl$IID), X$sample)
  j <- as.integer(cl$j)
  ok <- !is.na(i) & !is.na(j) & j >= 1L & j <= X$p
  o <- order(j[ok], i[ok])
  list(ptr = c(0L, cumsum(tabulate(j[ok], nbins = X$p))), i = i[ok][o] - 1L, g = as.numeric(cl$g)[ok][o])
}

#' Defaults of [genotype_fill()]
#'
#' delta 0.01 ridge of stage 1; anchor_min_miss 0.01 (stage 2 only at or above this missing rate); flank 2500
#' panel rows each side; thresh 0 = all stored eigenpairs; threads; mem_gb memory budget for the people batches
#' of stage 1; cond NULL = condition on every panel SNP of the window, else the column names of the conditioning
#' set (see `cond` of [genotype_impute()]; `snps` of [genotype_fill()] sets it); dosage_step 0.01 grid of the
#' dosage part of an observed cell in the LD input (see [genotype_impute()]). The LD fill uses all of them, the
#' other fills `threads` alone.
#' @export
genotype_fill_para <- function() {
  list(delta = 0.01, anchor_min_miss = 0.01, flank = 2500L, thresh = 0, threads = 4L, mem_gb = 2, cond = NULL,
       dosage_step = 0.01)
}
