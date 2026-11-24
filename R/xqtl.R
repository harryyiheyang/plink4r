#' Association of many continuous traits, each with the variants around its own anchor (cis-xQTL), written to disk
#'
#' Tests every trait of `pheno` (the expression of a gene, a protein, a methylation site ...) with the variants
#' within `window` base pairs of its anchor (e.g. the transcription start site), adjusting for the covariates of
#' `formula`: per pair the linear test of [gwas()], `lm(trait ~ covariates + g)`. The traits share their samples and
#' covariates, so the genotype side of the test is computed once per variant and the trait side once per trait; the
#' genotype file is streamed once, in position order, and a block of variants meets only the traits whose window it
#' falls in. Variants in no window are not read. Nothing is returned in memory: each trait is written to the Parquet
#' file `out` when the scan has passed its window, as one row group.
#'
#' Continuous traits only, and no missing value in `pheno`: run the tissues (sample sets) one at a time. A missing
#' genotype call is replaced by the variant's mean among the analysed samples.
#'
#' The file has one row per trait-variant pair, traits in the order of their windows along the genome and variants in
#' position order: `Outcome` (the trait ID), `CHR` (integer code: 23 = X), `BP`, `REF`, `ALT` (the counted allele,
#' which `BETA` is per copy of: ALT of a PGEN, .bim column 5 of a BED), `BETA`, `SE` (single precision), `N`
#' (non-missing calls of the variant among the analysed samples). A variant with no variation among the analysed
#' samples, or in the span of the covariates, has no row. The P value of a pair is
#' `2 * pt(-abs(BETA / SE), df)` with `df = n - k - 1` (`n` samples, `k` the rank of the covariates with the
#' intercept); `df`, `samples` and `covariate_rank` are in the key-value metadata of the file and in the value
#' returned here. Read one trait with `arrow::open_dataset(out) |> dplyr::filter(Outcome == "...")`: the row groups
#' carry the minimum and maximum of `Outcome`, `CHR` and `BP`, so only that trait is read.
#'
#' @param pheno Numeric matrix, one row per sample and one column per trait, with the sample IDs as row names and
#'   the trait IDs as column names. No missing value among the analysed samples.
#' @param anchor data.frame with one row per trait: `ID` (a column name of `pheno`), `CHR` and `POS` (the anchor).
#'   Traits of `pheno` without a row here, or on a chromosome without a variant in their window, are not tested.
#' @param out Path of the Parquet file to write (replaced when it exists).
#' @param formula The covariates, one-sided as for [stats::model.matrix()]: e.g. `~ age + sex + PC1`; `~ 1` without
#'   covariates. Every variable must be a column of `data`.
#' @param data data.frame with one row per sample: `id_col` and the columns named in `formula`; `NULL` with `~ 1`.
#'   Samples with a missing covariate, or not in `pheno` or the genotype file, are left out.
#' @param window Half-width of the window in base pairs: a trait is tested with the variants at `POS - window` to
#'   `POS + window` of its chromosome.
#' @inheritParams gwas
#' @param block Variants per thread and block. A block computes the products of its variants with every trait it
#'   meets, so a block much longer than a window computes products that no pair uses (`products` against `pairs` in
#'   the value); 0 = the largest that fits `mem_gb` (at most 4096).
#' @return Invisibly, a list: `out` (the file), `rows` (pairs written), `bytes` (size of the file), `df`, and the
#'   run: samples, traits tested, variants read, pairs asked for, covariate rank, block size, threads, memory budget,
#'   SIMD set, `active` (per block of variants the number of traits it met), `products` (trait-variant products
#'   computed; at least `pairs`), `waiting_most` (the most pairs held in memory at once), monomorphic and singular
#'   variants, traits not tested, skipped multiallelic variants.
#' @export
xqtl <- function(pheno, anchor, out, formula = ~ 1, data = NULL, pfile = NULL, bfile = NULL, window = 1e6,
                 id_col = "IID", snp_vec = NULL, threads = 0L, mem_gb = 0, block = 256L) {
  if (!is.matrix(pheno) || !is.numeric(pheno) || is.null(rownames(pheno)) || is.null(colnames(pheno)))
    stop("pheno must be a numeric matrix with sample IDs as row names and trait IDs as column names.", call. = FALSE)
  if (anyDuplicated(rownames(pheno)) || anyDuplicated(colnames(pheno)))
    stop("Sample IDs or trait IDs repeated in pheno.", call. = FALSE)
  anchor <- as.data.frame(anchor)
  miss <- setdiff(c("ID", "CHR", "POS"), names(anchor))
  if (length(miss)) stop("anchor lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyDuplicated(anchor$ID)) stop("Trait IDs repeated in anchor.", call. = FALSE)
  if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out))
    stop("out must be the path of the Parquet file to write.", call. = FALSE)
  if (!dir.exists(dirname(out))) stop("No such directory: ", dirname(out), call. = FALSE)
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

  # the alleles as codes into their distinct strings: REF = the other allele, ALT = the counted one
  alleles <- unique(c(var$A1[vrow], var$A2[vrow]))
  res <- xqtl_stream_cpp(gv$path, is.null(pfile), length(iid), as.integer(var$allele_ct), as.integer(vrow - 1L),
                         as.integer(srow - 1L), X, Y, as.integer(lo), as.integer(hi), path.expand(out),
                         as.character(a$ID), as.integer(var$CHR[vrow]), as.integer(var$POS[vrow]),
                         match(var$A2[vrow], alleles) - 1L, match(var$A1[vrow], alleles) - 1L, alleles,
                         paste("plink4r version", utils::packageVersion("plink4r")), threads, mem_gb,
                         as.integer(block))
  invisible(list(out = out, rows = res$rows, bytes = res$bytes, df = res$df, samples = nrow(d), traits = nrow(a),
                 variants = length(vrow), pairs = res$pairs, covariate_rank = res$rank, block = res$block,
                 threads = res$threads, mem_gb = res$mem_gb, simd = res$simd, dosage = res$dosage,
                 active = res$active, products = res$products, waiting_most = res$waiting_most,
                 explicit_projection = res$explicit, monomorphic = res$monomorphic, singular = res$singular,
                 traits_not_tested = not_tested, skipped_multiallelic = gv$skipped))
}
