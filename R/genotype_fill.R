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
#' Priority cells > ld > const > af > mean. A value f is stored as round(f) in the 2-bit planes and the
#' rest in the int16 part (1/16384). All values refer to the counted allele `x$a1` (after [geno_harmonise()],
#' the harmonised one).
#'
#' The file route of earlier versions is kept: when `x` is a character PLINK2 prefix the call is
#' `genotype_fill(pfile, snp, ld, sample, para)` and behaves as before (match `snp` to the .pvar, QC, per-SNP
#' fill constants and LD cells, for `geno_open(fill = )`). It stays the route for LD fills of a PGEN with
#' dosages (the object route takes hard calls only).
#'
#' @param x A `geno` object; or (old route) a PLINK2 prefix.
#' @param ... Old route only: `snp`, `ld`, `sample` by position or name (see `.genotype_fill_file`).
#' @param af Anchor frequency of the counted allele, one per column (`NA` = not anchored).
#' @param ld NULL, or list(dir = folder of block<N>.eigen.bin and the panel SNP .parquet, ref = optional
#'   parquet path or data.frame in place of that table).
#' @param const Fill constant per column (`NA` = not given).
#' @param cells data.frame(snp = column name (or j = column number), IID, g) of per-cell fills; cells that
#'   are not missing are ignored.
#' @param para named list overriding [genotype_fill_para()] (used: delta, anchor_min_miss, flank, thresh,
#'   threads, mem_gb; the old route uses all).
#' @return A `geno` object (new external pointer), `impute` = "mean", "af", "custom" or "ld". Old route: as
#'   before, list(snp, cells, dropped, method, para, time).
#' @export
genotype_fill <- function(x, ..., af = NULL, ld = NULL, const = NULL, cells = NULL, para = list()) {
  if (missing(x)) {  # old route called as genotype_fill(pfile = , ...)
    a <- list(...)
    if (is.null(a$pfile)) stop("x is missing.", call. = FALSE)
    x <- a$pfile
    a$pfile <- NULL
    return(do.call(genotype_fill, c(list(x), a, if (!missing(ld)) list(ld = ld), if (!missing(para)) list(para = para))))
  }
  if (is.character(x)) {
    a <- list(...)
    an <- names(a)
    if (is.null(an)) an <- rep("", length(a))
    bad <- setdiff(an[nzchar(an)], c("snp", "sample"))
    if (length(bad)) stop("unknown argument: ", paste(bad, collapse = ", "), call. = FALSE)
    if (!is.null(af) || !is.null(const) || !is.null(cells)) stop("af / const / cells need a geno object.", call. = FALSE)
    args <- list(snp = a$snp, ld = if (!missing(ld)) ld, sample = a$sample, para = if (!missing(para)) para)
    given <- c(snp = "snp" %in% an, ld = !missing(ld), sample = "sample" %in% an, para = !missing(para))
    for (v in a[!nzchar(an)]) {  # positional: snp, ld, sample, para in that order, skipping named ones
      slot <- names(given)[!given][1L]
      if (is.na(slot)) stop("too many arguments for genotype_fill(pfile, snp, ld, sample, para).", call. = FALSE)
      if (!is.null(v)) args[[slot]] <- v
      given[slot] <- TRUE
    }
    if (is.null(args$snp)) stop("genotype_fill(pfile, ...) needs `snp`.", call. = FALSE)
    return(.genotype_fill_file(x, args$snp, args$ld, args$sample, if (is.null(args$para)) list() else args$para))
  }
  if (!inherits(x, "geno")) stop("x must be a geno object (or a PLINK2 prefix).", call. = FALSE)
  if (length(list(...))) stop("extra arguments are only for genotype_fill(pfile, snp, ld, sample, para).", call. = FALSE)
  P <- genotype_fill_para()
  bad <- setdiff(names(para), names(P))
  if (length(bad)) stop("unknown para: ", paste(bad, collapse = ", "), call. = FALSE)
  P[names(para)] <- para
  P$threads <- .threads(P$threads)
  p <- x$p
  if (!is.null(const)) {
    const <- as.numeric(const)
    if (length(const) != p) stop("const needs one value per column.", call. = FALSE)
  }
  cl <- NULL
  if (!is.null(cells) && nrow(cells)) {
    j <- if (!is.null(cells$j)) as.integer(cells$j) else match(as.character(cells$snp), x$snp)
    cl <- data.frame(j = j, IID = as.character(cells$IID), g = as.numeric(cells$g))
  }
  lab <- if (!is.null(cl)) "custom" else if (!is.null(ld)) "ld" else if (any(is.finite(const))) "custom" else
    if (any(is.finite(af))) "af" else "mean"
  if (!is.null(ld) && any(x$n_miss > 0)) {
    P$mem_gb <- .mem_gb(P$mem_gb, P$threads)
    .check_delta(P$delta)
    if (x$dosage) stop("LD fill on an object needs hard calls; for a PGEN with dosages use the file route ",
                       "genotype_fill(pfile, snp, ld, sample) and geno_open(pgenfile =, fill = ).", call. = FALSE)
    pan <- .ld_panel(ld$ref, ld$dir)
    al <- .ld_align(pan, .geno_vartab(x))
    use <- sort(intersect(al$vi + 1L, which(x$n_miss > 0 & !(seq_len(p) %in% x$appended))))
    if (length(use)) {
      a <- if (!is.null(af)) { ua <- which(is.finite(af)); stats::setNames(af[ua], x$snp[ua]) }
      gi <- genotype_impute(x, pan, ld$dir, snps = x$snp[use], flank = P$flank, delta = P$delta,
                            anchor_min_miss = P$anchor_min_miss, thresh = P$thresh, threads = P$threads,
                            mem_gb = P$mem_gb, af = a)
      lc <- data.frame(j = gi$fill$j, IID = gi$fill$IID, g = gi$fill$g, i = gi$fill$i)
      if (!is.null(cl)) cl$i <- NA_integer_
      cl <- if (is.null(cl)) lc else rbind(cl, lc)  # user cells first: they win
    }
  }
  .geno_refill(x, af = af, const = const, cells = .geno_fill_spec(list(cells = cl), x$snp, x$sample),
               default = "mean", inplace = FALSE, impute = lab, threads = P$threads)
}

# Refill the missing cells of X through geno_fill_cpp. Per column f = default ("mean": observed mean, 0 when
# every cell is missing; "keep": the current fill), then af (anchored), then const; cells (CSR from
# .geno_fill_spec) win over all. inplace only for an object nothing else references.
.geno_refill <- function(X, af = NULL, const = NULL, cells = NULL, default = c("mean", "keep"), inplace = FALSE,
                         impute = X$impute, threads = X$threads) {
  default <- match.arg(default)
  p <- X$p; n <- X$n
  if (!is.null(af)) {
    af <- as.numeric(af)
    if (length(af) != p) stop("af needs one value per selected variant.", call. = FALSE)
    if (any(af < 0 | af > 1, na.rm = TRUE)) stop("af must lie in [0, 1].", call. = FALSE)
  }
  info <- geno_info_cpp(X$ptr)
  m <- info$n_miss; S <- info$obs_sum
  f <- if (default == "mean") ifelse(m < n, S / (n - m), 0) else rep(NaN, p)
  if (length(af)) {
    ok <- is.finite(af) & m > 0
    f[ok] <- pmin(2, pmax(0, (2 * n * af[ok] - S[ok]) / m[ok]))
  }
  if (length(const)) {
    ok <- is.finite(const)
    f[ok] <- const[ok]
  }
  f[m == 0] <- NaN
  if (is.null(cells) || !length(cells$ptr)) cells <- list(ptr = integer(0), i = integer(0), g = numeric(0))
  ptr <- geno_fill_cpp(X$ptr, f, cells$ptr, cells$i, cells$g, inplace, .threads(threads))
  .geno_new(ptr, .geno_cols(X), list(sample = X$sample, fid = X$fid), impute = impute, scale = X$scale,
            threads = X$threads, file = X$file, appended = X$appended)
}

# fill = list(const = one value per selected variant (NA = built-in fill),
#             cells = data.frame(j = column, IID, g)), or a genotype_fill() file-route result -> CSR over the
# selected rows (cells matched by IID; unmatched cells dropped; the first of repeated cells wins).
.geno_fill_spec <- function(fill, ids, iid) {
  p <- length(ids)
  out <- list(c = numeric(0), ptr = integer(0), i = integer(0), g = numeric(0))
  if (is.null(fill)) return(out)
  if (!is.null(fill$snp)) {  # a genotype_fill() result: by variant ID
    cl <- fill$cells
    fill <- list(const = fill$snp$fill_const[match(ids, fill$snp$ID)],
                 cells = if (!is.null(cl) && nrow(cl)) data.frame(j = match(cl$ID, ids), IID = cl$IID, g = cl$g))
  }
  if (!is.null(fill$const)) {
    if (length(fill$const) != p) stop("fill$const needs one value per selected variant.", call. = FALSE)
    out$c <- as.numeric(fill$const)
  }
  cl <- fill$cells
  if (!is.null(cl) && nrow(cl)) {
    i <- match(as.character(cl$IID), iid)
    if (!is.null(cl$i)) i[!is.na(cl$i)] <- as.integer(cl$i[!is.na(cl$i)])  # row numbers when known
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

# The file route of genotype_fill() (unchanged from plink4r 0.1.3): locate SNPs in a PGEN by CHR / BP / A1 /
# A2 and compute the fills of their missing genotypes, for geno_open(fill = ) (SuSiE4I calls it for
# X = list(pfile, snp, ld, para)):
# 1. SNPs of `snp` are matched to the .pvar by chromosome, position and the allele pair
#    ({A1, A2} = {ALT, REF}, either orientation; values are ALT dosages).
# 2. Per-SNP observed count and ALT-dosage sum over `sample` (dosage-aware) -> missing rate and observed AF;
#    SNPs with missing rate > `miss_max` or observed MAF < `maf` are dropped.
# 3. Fill of a missing cell: default observed mean dosage; `A1Freq` in `snp` (scheme 3): c_j = (2 n p_j -
#    sum_obs g) / m_j clipped to [0, 2], p_j the ALT frequency; `ld` (scheme 4): per-cell conditional mean from
#    the eigen.bin panel (genotype_impute()), anchored to `A1Freq` when given, else to the panel AF; SNPs
#    absent from the panel keep the per-SNP constant.
# pfile PLINK2 prefix; snp data.frame with CHR, BP37 or BP38, A1, A2, optional A1Freq; ld NULL or list(dir,
# ref); sample IIDs (NULL = all); para overrides genotype_fill_para().
# Returns list(snp = one row per kept SNP: ID (pvar), idx (1-based pvar row), CHR, POS, REF, ALT, A1, A2,
# a1_is_alt, n_miss, miss_rate, af_alt_obs, fill_const, in_ld; cells = data.table(ID, IID, g) per-cell fills
# (scheme 4 only, ALT dosage); dropped = SNPs not used with the reason; method; para; time).
.genotype_fill_file <- function(pfile, snp, ld = NULL, sample = NULL, para = list()) {
  t0 <- proc.time()[["elapsed"]]
  P <- genotype_fill_para()
  bad <- setdiff(names(para), names(P))
  if (length(bad)) stop("unknown para: ", paste(bad, collapse = ", "))
  P[names(para)] <- para
  P$threads <- .threads(P$threads); P$mem_gb <- .mem_gb(P$mem_gb, P$threads)
  if (!is.null(ld)) .check_delta(P$delta)
  snp <- data.table::copy(data.table::as.data.table(snp))
  bpcol <- if (!is.null(P$build)) paste0("BP", P$build) else intersect(c("BP37", "BP38"), names(snp))
  if (length(bpcol) != 1L) stop("`snp` needs exactly one of BP37 / BP38, or set para$build")
  miss_col <- setdiff(c("CHR", bpcol, "A1", "A2"), names(snp))
  if (length(miss_col)) stop("`snp` lacks ", paste(miss_col, collapse = ", "))
  chr_std <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

  pvar <- .read_pvar5(pfile)
  act <- pvar$allele_ct

  q <- data.table::data.table(k = seq_len(nrow(snp)), CHROM = chr_std(snp$CHR), POS = as.integer(snp[[bpcol]]),
                              A1 = toupper(snp$A1), A2 = toupper(snp$A2),
                              A1Freq = if ("A1Freq" %in% names(snp)) as.numeric(snp$A1Freq) else NA_real_)
  m <- merge(q, pvar[, .(CHROM, POS, ID, REF, ALT, row)], by = c("CHROM", "POS"), allow.cartesian = TRUE)
  m <- m[(A1 == ALT & A2 == REF) | (A1 == REF & A2 == ALT)]
  m[, a1_is_alt := A1 == ALT]
  m <- m[!duplicated(k)]
  dropped <- data.table::data.table(k = setdiff(q$k, m$k), reason = "not in pvar (CHR/BP/alleles)")
  data.table::setorder(m, k)

  iid <- .psam_ids(pfile)
  si <- .sample_rows(sample, iid)
  n <- length(si)
  pg <- normalizePath(paste0(pfile, ".pgen"))
  st <- pgen_stats_cpp(pg, length(iid), act, as.integer(m$row - 1L), as.integer(si - 1L), as.integer(P$threads))
  m[, `:=`(n_miss = n - st$n_obs, miss_rate = 1 - st$n_obs / n, af_alt_obs = st$sum / (2 * pmax(st$n_obs, 1)), sum_obs = st$sum)]
  r_miss <- m$miss_rate > P$miss_max
  r_maf <- !r_miss & pmin(m$af_alt_obs, 1 - m$af_alt_obs) < P$maf
  dropped <- rbind(dropped, data.table::data.table(k = m$k[r_miss], reason = paste0("miss_rate > ", P$miss_max)),
                   data.table::data.table(k = m$k[r_maf], reason = paste0("MAF < ", P$maf)))
  m <- m[!(r_miss | r_maf)]
  if (!nrow(m)) stop("no SNP left after matching and QC")

  m[, p_alt := data.table::fifelse(a1_is_alt, A1Freq, 1 - A1Freq)]
  m[, fill_const := data.table::fifelse(is.finite(p_alt) & n_miss > 0,
                                        pmin(2, pmax(0, (2 * n * p_alt - sum_obs) / pmax(n_miss, 1))),
                                        sum_obs / pmax(n - n_miss, 1))]
  method <- if (all(is.finite(m$p_alt))) "scheme3" else if (any(is.finite(m$p_alt))) "scheme3_partial" else "mean"

  cells <- data.table::data.table(ID = character(0), IID = character(0), g = numeric(0))
  m[, in_ld := FALSE]
  if (!is.null(ld)) {
    pan <- .ld_panel(ld$ref, ld$dir)
    al <- .ld_align(pan, pvar)
    m[, in_ld := (row - 1L) %in% al$vi]
    if (any(m$in_ld)) {
      ua <- m[in_ld == TRUE & is.finite(p_alt)]
      gi <- genotype_impute(pfile, pan, ld$dir, snps = m$ID[m$in_ld], flank = P$flank, sample = iid[si],
                            delta = P$delta, anchor_min_miss = P$anchor_min_miss, thresh = P$thresh,
                            threads = P$threads, mem_gb = P$mem_gb, af = stats::setNames(ua$p_alt, ua$ID))
      cells <- gi$fill[, .(ID = snp, IID, g)]
      method <- paste0("scheme4+", method)
    }
  }
  dropped <- merge(dropped, snp[, `:=`(k = .I)][, c("k", "CHR", bpcol, "A1", "A2"), with = FALSE], by = "k")
  list(snp = m[, .(ID, idx = row, CHR = CHROM, POS, REF, ALT, A1, A2, a1_is_alt, n_miss, miss_rate, af_alt_obs,
                   fill_const, in_ld)],
       cells = cells, dropped = dropped[, !"k"], method = method, para = P,
       time = proc.time()[["elapsed"]] - t0)
}

#' Defaults of [genotype_fill()]
#'
#' build NULL = the one BP column present; maf 0.01 (observed); miss_max 1 = keep all (drop SNPs whose
#' missing rate exceeds it); delta 0.01 ridge of stage 1; anchor_min_miss 0.01 (stage 2 only at or above
#' this missing rate); flank 2500 panel rows each side; thresh 0 = all stored eigenpairs; threads; mem_gb
#' memory budget for the people batches of stage 1.
#' @export
genotype_fill_para <- function() {
  list(build = NULL, maf = 0.01, miss_max = 1, delta = 0.01, anchor_min_miss = 0.01, flank = 2500L,
       thresh = 0, threads = 4L, mem_gb = 2)
}
