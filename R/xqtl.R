#' Association of many continuous outcomes, each with the variants around its own position (cis-xQTL), written to disk
#'
#' Three tables go in. `outcome` holds the outcomes of the samples (the expression of the genes, the proteins, the
#' methylation sites ...), one column each; `position` says where each outcome lies on the genome; `covariate`
#' holds the covariates of the samples, the same for every outcome. An outcome is tested with the variants from
#' `window` base pairs before its `Start` to `window` base pairs after its `End`: per pair the linear test of
#' [gwas()], `lm(outcome ~ covariates + g)`. The outcomes share their samples and covariates, so the genotype side
#' of the test is computed once per variant and the outcome side once per outcome; the genotype file is streamed
#' once, in position order, and a block of variants meets only the outcomes whose window it falls in. Variants in
#' no window are not read. Nothing is returned in memory: each outcome is written to the Parquet file `out` when the
#' scan has passed its window, as one row group.
#'
#' Continuous outcomes only, and no missing value in `outcome`: run the tissues (sample sets) one at a time. A
#' missing genotype call is replaced by the variant's mean among the analysed samples. The analysed samples are
#' those in `outcome`, in `covariate` with no missing covariate, and in the genotype file.
#'
#' The file has one row per outcome-variant pair, outcomes in the order of their windows along the genome and
#' variants in position order: `Outcome` (its name), `CHR` (integer code: 23 = X), `BP`, `REF`, `ALT` (the counted
#' allele, which `BETA` is per copy of: ALT of a PGEN, .bim column 5 of a BED), `BETA`, `SE` (single precision), `N`
#' (non-missing calls of the variant among the analysed samples). A variant with no variation among the analysed
#' samples, or in the span of the covariates, has no row. The P value of a pair is
#' `2 * pt(-abs(BETA / SE), df)` with `df = n - k - 1` (`n` samples, `k` the rank of the covariates with the
#' intercept); `df`, `samples` and `covariate_rank` are in the key-value metadata of the file and in the value
#' returned here. Read one outcome with `arrow::open_dataset(out) |> dplyr::filter(Outcome == "...")`: the row
#' groups carry the minimum and maximum of `Outcome`, `CHR` and `BP`, so only that outcome is read.
#'
#' @param outcome data.frame (data.table) with one row per sample: the sample IDs (.psam / .fam IID) in its first
#'   column, then one numeric column per outcome, named by the outcome. No missing value among the analysed samples.
#' @param position data.frame with one row per outcome: `Outcome` (a column name of `outcome`), `CHR`, `Start` and
#'   `End`, positions on `CHR` in the coordinates of the genotype file with `Start <= End`: the two ends of a gene
#'   or region, or `Start = End` for a single anchor (e.g. the transcription start site). Rows naming no column of
#'   `outcome` are ignored. Outcomes without a row here, without a position, or with no variant in their window are
#'   not tested (`outcomes_not_tested` in the value).
#' @param out Path of the Parquet file to write (replaced when it exists).
#' @param covariate `NULL` (no covariate: the intercept alone), or a data.frame with one row per sample: the sample
#'   IDs in its first column, then the covariates, every column of which is used for every outcome. Numeric columns
#'   enter as they are, factors and characters expanded by base R with the current contrasts; the intercept is
#'   always added. Samples with a missing covariate are left out.
#' @param window Base pairs added on both sides of an outcome's position: it is tested with the variants at
#'   `Start - window` to `End + window` of its chromosome, both ends included. With 0, `Start` to `End` is the
#'   tested region itself.
#' @inheritParams gwas
#' @param block Variants per thread and block. A block computes the products of its variants with every outcome it
#'   meets, so a block much longer than a window computes products that no pair uses (`products` against `pairs` in
#'   the value); 0 = the largest that fits `mem_gb` (at most 4096).
#' @return Invisibly, a list: `out` (the file), `rows` (pairs written), `bytes` (size of the file), `df`, and the
#'   run: samples, outcomes tested, variants read, pairs asked for, covariate rank, block size, threads, memory
#'   budget, SIMD set, `active` (per block of variants the number of outcomes it met), `products` (outcome-variant
#'   products computed; at least `pairs`), `waiting_most` (the most pairs held in memory at once), monomorphic and
#'   singular variants, outcomes not tested, skipped multiallelic variants.
#' @export
xqtl <- function(outcome, position, out, covariate = NULL, pfile = NULL, bfile = NULL, window = 1e6, snp_vec = NULL,
                 threads = 0L, mem_gb = 0, block = 256L) {
  some <- function(x) paste0(paste(utils::head(x, 5L), collapse = ", "), if (length(x) > 5L) ", ...")
  if (!is.data.frame(outcome) || ncol(outcome) < 2L)
    stop("outcome must be a data.frame: the sample IDs in its first column, then one column per outcome.",
         call. = FALSE)
  oname <- names(outcome)[-1L]
  if (anyNA(oname) || anyDuplicated(oname)) stop("Column names of outcome missing or repeated.", call. = FALSE)
  oid <- as.character(.subset2(outcome, 1L))
  if (anyNA(oid) || anyDuplicated(oid))
    stop("Sample IDs (the first column of outcome) missing or repeated.", call. = FALSE)
  num <- vapply(seq_along(oname) + 1L, function(j) is.numeric(.subset2(outcome, j)), NA)
  if (!all(num)) stop("Not numeric in outcome: ", some(oname[!num]), call. = FALSE)
  if (!is.data.frame(position)) stop("position must be a data.frame: Outcome, CHR, Start, End.", call. = FALSE)
  miss <- setdiff(c("Outcome", "CHR", "Start", "End"), names(position))
  if (length(miss)) stop("position lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (!is.null(covariate) && (!is.data.frame(covariate) || ncol(covariate) < 2L))
    stop("covariate must be a data.frame: the sample IDs in its first column, then one column per covariate ",
         "(NULL without covariates).", call. = FALSE)
  if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out))
    stop("out must be the path of the Parquet file to write.", call. = FALSE)
  if (!dir.exists(dirname(out))) stop("No such directory: ", dirname(out), call. = FALSE)
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 0)
    stop("window must be a non-negative number of base pairs.", call. = FALSE)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  gv <- .geno_variants(pfile, bfile, snp_vec)
  var <- gv$var; iid <- gv$iid; vrow <- gv$vrow

  # ---- samples and covariates: the rules of gwas() (.pheno_samples), every column of covariate a term of the model
  if (is.null(covariate)) {
    d <- data.frame(IID = oid, stringsAsFactors = FALSE); term <- "1"
  } else {
    d <- as.data.frame(covariate)
    if (anyNA(names(d)) || anyDuplicated(names(d)) || any(grepl("`", names(d), fixed = TRUE)))
      stop("Column names of covariate missing, repeated or with a backtick.", call. = FALSE)
    term <- sprintf("`%s`", names(d)[-1L])
    cid <- as.character(d[[1L]])
    if (anyDuplicated(cid[!is.na(cid)])) stop("Sample IDs (the first column of covariate) repeated.", call. = FALSE)
    d <- d[cid %in% oid, , drop = FALSE]
    if (!nrow(d)) stop("No sample of covariate is in outcome.", call. = FALSE)
  }
  id_col <- names(d)[1L]
  d[[".plink4r_y"]] <- rep(0, nrow(d))
  sm <- .pheno_samples(d, stats::reformulate(term, response = ".plink4r_y"), "gaussian", id_col, iid)
  d <- sm$d; srow <- sm$srow
  X <- stats::model.matrix(stats::reformulate(term), d)
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL

  # ---- windows as ranges of the scan variants, which are taken in position order
  vrow <- vrow[!is.na(var$CHR[vrow]) & !is.na(var$POS[vrow])]
  vrow <- vrow[order(var$CHR[vrow], var$POS[vrow])]
  span <- 1e10  # above any position: chromosome code x span + position orders the genome
  key <- var$CHR[vrow] * span + var$POS[vrow]
  a <- as.data.frame(position)[, c("Outcome", "CHR", "Start", "End")]
  a$Outcome <- as.character(a$Outcome)
  a <- a[a$Outcome %in% oname, , drop = FALSE]
  if (anyDuplicated(a$Outcome))
    stop("Repeated in position: ", some(unique(a$Outcome[duplicated(a$Outcome)])), call. = FALSE)
  if (!is.numeric(a$Start) || !is.numeric(a$End)) stop("Start and End of position must be numeric.", call. = FALSE)
  a$CHR <- .chr_code(a$CHR)
  a <- a[!is.na(a$CHR) & is.finite(a$Start) & is.finite(a$End), , drop = FALSE]
  if (any(a$Start > a$End)) stop("Start after End in position: ", some(a$Outcome[a$Start > a$End]), call. = FALSE)
  lo <- findInterval(a$CHR * span + pmax(a$Start - window, 0), key, left.open = TRUE)
  hi <- findInterval(a$CHR * span + pmin(a$End + window, span - 1), key)
  keep <- hi > lo
  not_tested <- setdiff(oname, a$Outcome[keep])
  a <- a[keep, , drop = FALSE]; lo <- lo[keep]; hi <- hi[keep]
  if (!nrow(a)) stop("No outcome has a variant in its window.", call. = FALSE)
  o <- order(lo, hi)
  a <- a[o, , drop = FALSE]; lo <- lo[o]; hi <- hi[o]
  # variants in no window are left out, and the ranges renumbered
  p <- length(vrow)
  covered <- cumsum(tabulate(lo + 1L, p + 1L) - tabulate(hi + 1L, p + 1L))[seq_len(p)] > 0L
  before <- c(0L, cumsum(covered))
  lo <- before[lo + 1L]; hi <- before[hi + 1L]
  vrow <- vrow[covered]

  # the tested outcomes of the analysed samples, column by column: the table itself is not copied
  orow <- match(as.character(d[[id_col]]), oid); ocol <- match(a$Outcome, oname) + 1L
  Y <- vapply(ocol, function(j) as.numeric(.subset2(outcome, j)[orow]), numeric(length(orow)))
  dim(Y) <- c(length(orow), length(ocol))
  bad <- which(colSums(!is.finite(Y)) > 0)
  if (length(bad))
    stop("Missing or not finite among the analysed samples, in outcome: ", some(a$Outcome[bad]), call. = FALSE)

  # the alleles as codes into their distinct strings: REF = the other allele, ALT = the counted one
  alleles <- unique(c(var$A1[vrow], var$A2[vrow]))
  res <- xqtl_stream_cpp(gv$path, is.null(pfile), length(iid), as.integer(var$allele_ct), as.integer(vrow - 1L),
                         as.integer(srow - 1L), X, Y, as.integer(lo), as.integer(hi), path.expand(out),
                         a$Outcome, as.integer(var$CHR[vrow]), as.integer(var$POS[vrow]),
                         match(var$A2[vrow], alleles) - 1L, match(var$A1[vrow], alleles) - 1L, alleles,
                         paste("plink4r version", utils::packageVersion("plink4r")), threads, mem_gb,
                         as.integer(block))
  invisible(list(out = out, rows = res$rows, bytes = res$bytes, df = res$df, samples = nrow(d), outcomes = nrow(a),
                 variants = length(vrow), pairs = res$pairs, covariate_rank = res$rank, block = res$block,
                 threads = res$threads, mem_gb = res$mem_gb, simd = res$simd, dosage = res$dosage,
                 active = res$active, products = res$products, waiting_most = res$waiting_most,
                 explicit_projection = res$explicit, monomorphic = res$monomorphic, singular = res$singular,
                 outcomes_not_tested = not_tested, skipped_multiallelic = gv$skipped))
}
