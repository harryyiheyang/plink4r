#' Fill the missing genotype cells of a geno object
#'
#' Returns a new object whose missing cells (the cells missing in the file, see [geno_missing()]) are
#' rewritten; the observed cells and `x` itself are unchanged. For a column j with m_j missing cells,
#' n rows and observed sum S_j (counted allele = `x$A1`):
#' - default: f_j = S_j / (n - m_j), the observed mean of `x`'s own rows;
#' - `af` (finite entries): f_j = clip((2 n af_j - S_j) / m_j, 0, 2), so the column mean is 2 af_j;
#' - `const` (finite entries): f_j = const_j;
#' - `ld`: per-cell conditional means from the eigen.bin panel (the LD fill below), anchored to `af` when given,
#'   else to the panel AF; columns absent from the panel keep f_j;
#' - `cells`: the given value of that cell.
#' Priority cells > ld > const > af > mean.
#'
#' `snps` restricts the call to a subset of columns: only their missing cells are rewritten (every other column
#' keeps its current fill, and its `af`, `const` and `cells` entries are ignored), and the LD fill conditions on
#' their observed calls alone, so columns outside `snps`, e.g. wrong calls
#' that are being replaced, carry no information. Nothing is copied for the subset: the kernel reads the
#' selected columns of `x` in place. A value f is stored as round(f) in the 2-bit planes and the
#' rest in the int16 part (1/16384; 0.01 for a dosage PGEN). All values refer to the counted allele `x$A1` (after [geno_harmonise()],
#' the harmonised one). Hard calls and dosages are taken alike: observed dosages condition the LD fill and
#' missing dosages are filled.
#'
#' The LD fill works on the panel SNPs that are columns of `x` (panel SNPs absent from `x` are not carried). The
#' panel need not share the object's allele orientation or build: SNPs are aligned by chromosome, position
#' and allele pair, and a swapped pair enters with a negative sign.
#' Stage 1 (OpenMP over people): z_M = R_MO (R_OO + delta I)^{-1} z_O per person on the working set S,
#' standardized with the anchor AF (`af`, else the panel's A1Freq in the object's orientation), from
#' C_S = I - delta (R_SS + delta I)^{-1}, built by Woodbury in the k eigenpairs or directly in |S|, whichever is
#' smaller; delta must be > 0. S is the panel rows of the columns to fill plus `flank` eigen rows on each side
#' inside the same LD block (blocks are independent in the panel, so a region spanning blocks is filled block by
#' block), restricted to the conditioning set when one is given. Observed dosages enter z_O as they are (no
#' hard-call threshold), on the `dosage_step` grid.
#' Stage 2: for columns with missing rate >= anchor_min_miss, the filled cells are shifted so that observed +
#' filled allele count = 2 n p, then clipped to `[0, 2]`; below the cutoff stage 1 is kept. anchor_min_miss < 0
#' disables it: `para = list(anchor_min_miss = -1)` and then [geno_missing()] give the stage 1 values.
#'
#' @param x A `geno` object.
#' @param af Anchor frequency of the counted allele, one per column (`NA` = not anchored).
#' @param ld NULL, or list(dir = the LD panel folder (`snpinfo.parquet` and block<N>.eigen.bin), ref = optional
#'   parquet path or data.frame in place of `snpinfo.parquet`).
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
    .check_dosage_step(P$dosage_step)
    al <- .ld_align(.ld_panel(ld$ref, ld$dir), .geno_vartab(x))
    use <- sort(intersect(al$vi + 1L, which(sel & x$n_miss > 0)))
    if (length(use)) {
      cond <- if (!all(sel)) which(sel) else if (!is.null(P$cond)) which(x$snp %in% P$cond)
      # the kernel writes the LD fills straight into the new object (no per-cell table in R)
      ld_write <- function(ptr) .genotype_fill_ld(ptr, al, ld$dir, use, af, cond, P)
    }
  }
  .geno_refill(x, af = af, const = const, cells = .geno_cells_csr(cl, x), inplace = FALSE, impute = lab,
               threads = P$threads, cols = if (!all(sel)) sel, ld_write = ld_write)
}

# Refill the missing cells of X through geno_fill_cpp, which computes the per-column value f: the observed mean (0
# when every cell is missing), then af (anchored), then const; cells (CSR from .geno_cells_csr) win over all.
# inplace only for an object nothing else references. ld_write(ptr), when given, writes the LD fills into the new
# store after the per-column values and before the cells (so cells > LD > f).
.geno_refill <- function(X, af = NULL, const = NULL, cells = NULL, inplace = FALSE, impute = X$impute,
                         threads = X$threads, cols = NULL, ld_write = NULL) {
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
  af <- if (is.null(af)) numeric(0) else af
  const <- if (is.null(const)) numeric(0) else const
  cols <- if (is.null(cols)) logical(0) else as.logical(cols)
  if (is.null(cells) || !length(cells$ptr)) cells <- list(ptr = integer(0), i = integer(0), g = numeric(0))
  th <- .threads(threads)
  if (is.null(ld_write)) {
    ptr <- geno_fill_cpp(X$ptr, TRUE, af, const, cols, cells$ptr, cells$i, cells$g, inplace, th)
  } else {
    ptr <- geno_fill_cpp(X$ptr, TRUE, af, const, cols, integer(0), integer(0), numeric(0), inplace, th)
    ld_write(ptr)
    if (length(cells$ptr))
      ptr <- geno_fill_cpp(ptr, FALSE, numeric(0), numeric(0), logical(0), cells$ptr, cells$i, cells$g, TRUE, th)
  }
  .geno_new(ptr, .geno_cols(X), list(sample = X$sample, fid = X$fid), impute = impute, scale = X$scale,
            threads = X$threads, file = X$file)
}

# The LD fill of genotype_fill(): writes the final fills of the missing cells of the columns `use` into the store ptr
# (a fresh copy of the object's store; genotype_fill_ld_cpp reads the genotypes from it and writes into it), one LD
# block after another. al = the panel aligned to the columns of the object (.ld_align: vi = 0-based column, sign = +1
# when the counted allele is the panel's A1), af = NULL or the anchor frequency of the counted allele per column (finite
# entries replace the panel's A1Freq), cond = NULL or the columns of the conditioning set, P = the genotype_fill_para()
# values (threads and mem_gb resolved). The working rows of a block are the panel rows of its columns to fill plus
# P$flank eigen rows on each side, with an AF in (0, 1) and in cond.
.genotype_fill_ld <- function(ptr, al, ld_dir, use, af, cond, P) {
  if (!is.null(af)) {  # anchors are frequencies of the counted allele; the panel column is in panel orientation
    pa <- as.numeric(af)[al$vi + 1L]
    ok <- is.finite(pa)
    al$A1Freq[ok] <- ifelse(al$sign[ok] > 0, pa[ok], 1 - pa[ok])
  }
  rs <- al[(al$vi + 1L) %in% use]
  for (blk in sort(unique(rs$Block))) {
    li <- rs$Index[rs$Block == blk]
    keep <- rs$vi[rs$Block == blk]
    w <- .ld_rows(al, blk, min(li) - P$flank, max(li) + P$flank)
    w <- w[A1Freq > 0 & A1Freq < 1]
    if (!is.null(cond)) {
      inc <- (w$vi + 1L) %in% cond
      if (any(lost <- !inc & w$vi %in% keep))
        warning(sum(lost), " column(s) to fill in block ", blk, " are not in `cond` and get no LD fill", call. = FALSE)
      w <- w[inc]
    }
    if (!nrow(w)) stop("no SNP of block ", blk, " found in the geno object", call. = FALSE)
    genotype_fill_ld_cpp(ptr, as.integer(w$vi), .ld_block_file(ld_dir, blk),
                         ifelse(w$sign > 0, w$A1Freq, 1 - w$A1Freq), as.integer(w$sign), P$thresh, P$delta,
                         P$anchor_min_miss, w$vi %in% keep, as.integer(P$threads), 0L, P$mem_gb,
                         as.integer(w$Index), P$dosage_step)
  }
  invisible(NULL)
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
#' panel rows each side; thresh 0 = all stored eigenpairs (else the eigenvalue variance kept); threads; mem_gb
#' memory budget for the people batches of stage 1 (0 takes it from [plink4r_resources()]); cond NULL = condition
#' on every panel SNP of the window, else the column names of the conditioning set: a person's missing cells are
#' then predicted from their observed calls at these columns alone, and a column to fill that is not among them
#' gets no LD fill (warning; `snps` of [genotype_fill()] sets it); dosage_step 0.01 grid of the dosage part of an
#' observed cell in the LD input: an observed dosage g enters stage 1 as
#' round(g) + dosage_step * round((g - round(g)) / dosage_step) (hard calls are exact; only this input is
#' quantized, observed values are never changed; must be >= 0.5 / 127, the step count is int8). The LD fill uses
#' all of them, the other fills `threads` alone.
#' @export
genotype_fill_para <- function() {
  list(delta = 0.01, anchor_min_miss = 0.01, flank = 2500L, thresh = 0, threads = 4L, mem_gb = 2, cond = NULL,
       dosage_step = 0.01)
}
