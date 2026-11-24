#' Association of many continuous traits, each with the variants around its own anchor (cis-eQTL)
#'
#' Tests every trait of `pheno` (e.g. the expression of a gene) with the variants within `window` base pairs of its
#' anchor (e.g. the transcription start site), adjusting for the covariates of `formula`: per pair the linear test of
#' [gwas()], `lm(trait ~ covariates + g)`. The traits share their samples and covariates, so the genotype side of
#' the test is computed once per variant and the trait side once per trait; the genotype file is streamed once, in
#' position order, and a block of variants meets only the traits whose window it falls in (see `info`). Variants in
#' no window are not read.
#'
#' Continuous traits only. A missing genotype call is replaced by the variant's mean among the analysed samples.
#'
#' @param pheno Numeric matrix, one row per sample and one column per trait, with the sample IDs as row names and
#'   the trait IDs as column names. No missing value among the analysed samples.
#' @param anchor data.frame with one row per trait: `ID` (a column name of `pheno`), `CHR` and `POS` (the anchor).
#'   Traits of `pheno` without a row here, or on a chromosome without a variant in their window, are not tested.
#' @param formula The covariates, one-sided as for [stats::model.matrix()]: e.g. `~ age + sex + PC1`; `~ 1` without
#'   covariates. Every variable must be a column of `data`.
#' @param data data.frame with one row per sample: `id_col` and the columns named in `formula`; `NULL` with `~ 1`.
#'   Samples with a missing covariate, or not in `pheno` or the genotype file, are left out.
#' @param window Half-width of the window in base pairs: a trait is tested with the variants at `POS - window` to
#'   `POS + window` of its chromosome.
#' @inheritParams gwas
#' @param block Variants per thread and block. A block computes the products of its variants with every trait it
#'   meets, so a block much longer than a window computes products that no pair uses (`products` against `pairs` in
#'   `info`); 0 = the largest that fits `mem_gb` (at most 4096), as [gwas()].
#' @return A list of class `eqtl`, which holds every pair once and in compact form: `beta` and `se`, one entry per
#'   trait-variant pair, trait by trait; `traits`, a data.table with one row per trait tested in the order of their
#'   windows along the genome: `PHENO`, `CHR`, `POS` (the anchor), `FIRST` (the row of `variants` of its first
#'   variant) and `N_VAR` (its variants, which are consecutive rows of `variants`; its pairs are the next `N_VAR`
#'   entries of `beta` and `se` after those of the traits before it); `variants`, a data.table with one row per
#'   variant read, in position order: `CHR`, `POS`, `ID`, `A1` (counted allele), `A2`, `N`, `A1_FREQ`; `df`, the
#'   `n - k - 1` degrees of freedom of the t statistics; and `info`: samples, traits tested, variants read, pairs,
#'   covariate rank, block size, threads, memory budget, SIMD set, `active` (per block of variants the number of
#'   traits it met), `products` (trait-variant products computed; at least `pairs`), monomorphic and singular
#'   variants, traits not tested, skipped multiallelic variants. Variants with no variation among the analysed
#'   samples or in the span of the covariates have `NA` statistics. [eqtl_pairs()] gives the pairs as a table.
#' @export
eqtl <- function(pheno, anchor, formula = ~ 1, data = NULL, pfile = NULL, bfile = NULL, window = 1e6,
                 id_col = "IID", snp_vec = NULL, threads = 0L, mem_gb = 0, block = 256L) {
  if (!is.matrix(pheno) || !is.numeric(pheno) || is.null(rownames(pheno)) || is.null(colnames(pheno)))
    stop("pheno must be a numeric matrix with sample IDs as row names and trait IDs as column names.", call. = FALSE)
  if (anyDuplicated(rownames(pheno)) || anyDuplicated(colnames(pheno)))
    stop("Sample IDs or trait IDs repeated in pheno.", call. = FALSE)
  anchor <- as.data.frame(anchor)
  miss <- setdiff(c("ID", "CHR", "POS"), names(anchor))
  if (length(miss)) stop("anchor lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyDuplicated(anchor$ID)) stop("Trait IDs repeated in anchor.", call. = FALSE)
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 0)
    stop("window must be a non-negative number of base pairs.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 2L)
    stop("formula must be one-sided: ~ covariates.", call. = FALSE)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  gv <- .geno_variants(pfile, bfile, snp_vec)
  var <- gv$var; iid <- gv$iid; vrow <- gv$vrow

  # ---- samples and covariates: the rules of gwas() (.pheno_samples), the trait matrix standing for the phenotype
  d <- if (is.null(data)) data.frame(ids = rownames(pheno)) else as.data.frame(data)
  if (is.null(data)) names(d) <- id_col
  if (!id_col %in% names(d)) stop("data lacks column: ", id_col, call. = FALSE)
  d <- d[as.character(d[[id_col]]) %in% rownames(pheno), , drop = FALSE]
  d[[".plink4r_y"]] <- rep(0, nrow(d))
  sm <- .pheno_samples(d, stats::update(formula, .plink4r_y ~ .), "gaussian", id_col, iid)
  d <- sm$d; srow <- sm$srow
  X <- stats::model.matrix(formula, d)
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL

  # ---- windows as ranges of the scan variants, which are taken in position order
  vrow <- vrow[!is.na(var$CHR[vrow]) & !is.na(var$POS[vrow])]
  vrow <- vrow[order(var$CHR[vrow], var$POS[vrow])]
  span <- 1e10  # above any position: chromosome code x span + position orders the genome
  key <- var$CHR[vrow] * span + var$POS[vrow]
  a <- anchor[anchor$ID %in% colnames(pheno), c("ID", "CHR", "POS"), drop = FALSE]
  a$CHR <- .chr_code(a$CHR)
  a <- a[!is.na(a$CHR) & is.finite(a$POS), , drop = FALSE]
  lo <- findInterval(a$CHR * span + pmax(a$POS - window, 0), key, left.open = TRUE)
  hi <- findInterval(a$CHR * span + pmin(a$POS + window, span - 1), key)
  keep <- hi > lo
  not_tested <- setdiff(colnames(pheno), a$ID[keep])
  a <- a[keep, , drop = FALSE]; lo <- lo[keep]; hi <- hi[keep]
  if (!nrow(a)) stop("No trait has a variant in its window.", call. = FALSE)
  o <- order(lo, hi)
  a <- a[o, , drop = FALSE]; lo <- lo[o]; hi <- hi[o]
  # variants in no window are left out, and the ranges renumbered
  p <- length(vrow)
  covered <- cumsum(tabulate(lo + 1L, p + 1L) - tabulate(hi + 1L, p + 1L))[seq_len(p)] > 0L
  before <- c(0L, cumsum(covered))
  lo <- before[lo + 1L]; hi <- before[hi + 1L]
  vrow <- vrow[covered]

  Y <- pheno[match(as.character(d[[id_col]]), rownames(pheno)), a$ID, drop = FALSE]
  storage.mode(Y) <- "double"
  if (anyNA(Y)) stop("pheno has a missing value among the analysed samples.", call. = FALSE)

  res <- eqtl_stream_cpp(gv$path, is.null(pfile), length(iid), as.integer(var$allele_ct), as.integer(vrow - 1L),
                         as.integer(srow - 1L), X, Y, as.integer(lo), as.integer(hi), threads, mem_gb,
                         as.integer(block))
  variants <- var[vrow, list(CHR, POS, ID, A1, A2)]
  variants[, `:=`(N = as.integer(res$n_obs), A1_FREQ = res$a1_sum / (2 * res$n_obs))]
  traits <- data.table::data.table(PHENO = a$ID, CHR = a$CHR, POS = a$POS, FIRST = lo + 1L, N_VAR = hi - lo)
  info <- list(samples = nrow(d), traits = nrow(a), variants = length(vrow), pairs = res$pairs,
               covariate_rank = res$rank, block = res$block, threads = res$threads, mem_gb = res$mem_gb,
               simd = res$simd, dosage = res$dosage, active = res$active, products = res$products,
               explicit_projection = res$explicit, monomorphic = sum(res$status == 1L),
               singular = sum(res$status == 2L), traits_not_tested = not_tested, skipped_multiallelic = gv$skipped)
  structure(list(beta = res$beta, se = res$se, traits = traits, variants = variants, df = res$df, info = info),
            class = "eqtl")
}

#' The pairs of an eqtl() result as a table
#'
#' @param x An [eqtl()] result.
#' @param pheno Trait IDs (`NULL` = every trait tested), given in the order wanted.
#' @param p_max Keep the pairs with `P <= p_max` (1 = every pair, those without statistics included).
#' @return data.table with one row per pair, trait by trait and variants in position order: `PHENO`, `CHR`, `POS`,
#'   `ID`, `A1`, `A2`, `N`, `A1_FREQ`, `BETA`, `SE`, `T`, `P` (t distribution, `x$df` degrees of freedom).
#' @export
eqtl_pairs <- function(x, pheno = NULL, p_max = 1) {
  if (!inherits(x, "eqtl")) stop("x must be an eqtl() result.", call. = FALSE)
  if (!is.numeric(p_max) || length(p_max) != 1L || is.na(p_max) || p_max < 0 || p_max > 1)
    stop("p_max must be a number between 0 and 1.", call. = FALSE)
  tr <- x$traits
  i <- if (is.null(pheno)) seq_len(nrow(tr)) else match(as.character(pheno), tr$PHENO)
  if (anyNA(i)) stop("Not a trait of x: ", paste(pheno[is.na(i)], collapse = ", "), call. = FALSE)
  len <- tr$N_VAR[i]
  k <- sequence(len, from = (cumsum(as.numeric(tr$N_VAR)) - tr$N_VAR)[i] + 1)  # entries of beta / se
  vi <- sequence(len, from = tr$FIRST[i])                                      # rows of variants
  ti <- rep(i, len)
  beta <- x$beta[k]; se <- x$se[k]
  stat <- beta / se
  pv <- 2 * stats::pt(-abs(stat), x$df)
  if (p_max < 1) {
    keep <- which(pv <= p_max)
    vi <- vi[keep]; ti <- ti[keep]; beta <- beta[keep]; se <- se[keep]; stat <- stat[keep]; pv <- pv[keep]
  }
  out <- x$variants[vi]
  out[, `:=`(PHENO = tr$PHENO[ti], BETA = beta, SE = se, T = stat, P = pv)]
  data.table::setcolorder(out, "PHENO")
  out[]
}

#' @export
print.eqtl <- function(x, ...) {
  cat(sprintf("eqtl: %d traits x the variants of their windows = %.0f pairs (%d samples, %d variants read, df %g)
",
              nrow(x$traits), x$info$pairs, x$info$samples, nrow(x$variants), x$df))
  cat("beta, se: one entry per pair; traits, variants: their tables; eqtl_pairs(x, pheno, p_max): the pairs as a table
")
  invisible(x)
}
