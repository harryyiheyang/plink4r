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
#'   column, whatever its name, then one numeric column per outcome, named by the outcome. Every column is an
#'   outcome to test: give the columns wanted. No missing value among the analysed samples.
#' @param position data.frame with the position of every outcome: `Outcome` (the column name in `outcome`), `CHR`,
#'   `Start` and `End`, positions on `CHR` in the coordinates of the genotype file with `Start <= End`: the two ends
#'   of a gene or region, or `Start = End` for a single anchor (e.g. the transcription start site). A column of
#'   `outcome` without a row here is an error; rows naming no column of `outcome` are ignored (the table may be that
#'   of the whole genome). An outcome whose row has no chromosome that the package reads (1 to 22, X, Y, MT), no
#'   `Start` or `End`, or no variant in its window (e.g. on a chromosome that the genotype file does not hold) is
#'   not tested (`outcomes_not_tested` in the value).
#' @param covariate `NULL` (no covariate: the intercept alone), or a data.frame with one row per sample: the sample
#'   IDs in its first column, whatever its name, then the covariates, every column of which is used for every
#'   outcome. Numeric columns enter as they are, factors and characters expanded by base R with the current
#'   contrasts; the intercept is always added. Samples with a missing covariate are left out.
#' @param out Path of the Parquet file to write (replaced when it exists).
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
#'   singular variants, outcomes not tested, and `skipped_multiallelic`: the `CHR` and `BP` of the variants left out
#'   for having more than two alleles.
#' @export
xqtl <- function(outcome, position, covariate = NULL, out, pfile = NULL, bfile = NULL, window = 1e6, snp_vec = NULL,
                 threads = 0L, mem_gb = 0, block = 256L) {
  some <- function(x) paste0(paste(utils::head(x, 5L), collapse = ", "),
                             if (length(x) > 5L) sprintf(", ... (%d in all)", length(x)))
  if (!is.data.frame(outcome) || ncol(outcome) < 2L)
    stop("outcome must be a data.frame: the sample IDs in its first column, then one column per outcome.",
         call. = FALSE)
  oname <- names(outcome)[-1L]
  if (anyNA(oname) || anyDuplicated(oname)) stop("Column names of outcome missing or repeated.", call. = FALSE)
  oid <- .id_text(.subset2(outcome, 1L))
  if (anyNA(oid) || anyDuplicated(oid))
    stop("Sample IDs (the first column of outcome) missing or repeated.", call. = FALSE)
  num <- vapply(seq_along(oname) + 1L, function(j) is.numeric(.subset2(outcome, j)), NA)
  if (!all(num)) stop("Not numeric in outcome: ", some(oname[!num]), call. = FALSE)
  if (!is.null(covariate) && (!is.data.frame(covariate) || ncol(covariate) < 2L))
    stop("covariate must be a data.frame: the sample IDs in its first column, then one column per covariate ",
         "(NULL without covariates).", call. = FALSE)
  if (!is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out))
    stop("out must be the path of the Parquet file to write.", call. = FALSE)
  if (!dir.exists(dirname(out))) stop("No such directory: ", dirname(out), call. = FALSE)
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 0)
    stop("window must be a non-negative number of base pairs.", call. = FALSE)

  # ---- position: its rows for the columns of outcome, every one of which must have its row
  if (!is.data.frame(position)) stop("position must be a data.frame: Outcome, CHR, Start, End.", call. = FALSE)
  miss <- setdiff(c("Outcome", "CHR", "Start", "End"), names(position))
  if (length(miss)) stop("position lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  a <- as.data.frame(position)[, c("Outcome", "CHR", "Start", "End")]
  a$Outcome <- as.character(a$Outcome)
  a <- a[a$Outcome %in% oname, , drop = FALSE]
  if (anyDuplicated(a$Outcome))
    stop("Repeated in position: ", some(unique(a$Outcome[duplicated(a$Outcome)])), call. = FALSE)
  lack <- setdiff(oname, a$Outcome)
  if (length(lack)) stop("No row of position for outcome column(s): ", some(lack), call. = FALSE)
  if (!is.numeric(a$Start) || !is.numeric(a$End)) stop("Start and End of position must be numeric.", call. = FALSE)
  a$CHR <- .chr_code(a$CHR)
  a <- a[!is.na(a$CHR) & is.finite(a$Start) & is.finite(a$End), , drop = FALSE]
  if (!nrow(a))
    stop("No outcome has a chromosome (1 to 22, X, Y, MT; with or without chr), a Start and an End in position.",
         call. = FALSE)
  if (any(a$Start > a$End)) stop("Start after End in position: ", some(a$Outcome[a$Start > a$End]), call. = FALSE)

  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  gv <- .geno_variants(pfile, bfile, snp_vec)
  var <- gv$var; iid <- gv$iid; vrow <- gv$vrow

  # ---- samples and covariates: the rules of gwas() (.pheno_samples), every column of covariate a term of the model.
  # .pheno_samples takes its samples with a phenotype: here a column of zeros under a name that covariate does not have
  if (is.null(covariate)) {
    d <- data.frame(IID = oid, stringsAsFactors = FALSE)
  } else {
    d <- as.data.frame(covariate)
    if (anyNA(names(d)) || !all(nzchar(names(d))) || anyDuplicated(names(d)))
      stop("Column names of covariate missing or repeated.", call. = FALSE)
    d[[1L]] <- .id_text(d[[1L]])
    d <- d[d[[1L]] %in% oid, , drop = FALSE]
    if (!nrow(d)) stop("No sample of covariate is in outcome.", call. = FALSE)
  }
  id_col <- names(d)[1L]
  rhs <- Reduce(function(x, y) call("+", x, y), lapply(names(d)[-1L], as.name), 1)
  y0 <- make.unique(c(names(d), "y"))[ncol(d) + 1L]
  d[[y0]] <- rep(0, nrow(d))
  sm <- .pheno_samples(d, eval(call("~", as.name(y0), rhs)), "gaussian", id_col, iid)
  d <- sm$d; srow <- sm$srow
  X <- stats::model.matrix(eval(call("~", rhs)), d)
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL

  # ---- the tables go to C++ as they are: it puts the variants in position order, places the windows among them and
  # reads the outcomes of the analysed samples from the columns of outcome
  res <- xqtl_stream_cpp(gv$path, is.null(pfile), length(iid), as.integer(var$allele_ct), as.integer(vrow - 1L),
                         as.integer(var$CHR), as.integer(var$POS), var$A2, var$A1, as.integer(srow - 1L), X,
                         outcome, match(d[[id_col]], oid) - 1L, match(a$Outcome, oname), a$Outcome, a$CHR,
                         as.numeric(a$Start), as.numeric(a$End), window, path.expand(out),
                         paste("plink4r version", utils::packageVersion("plink4r")), threads, mem_gb,
                         as.integer(block))
  invisible(list(out = out, rows = res$rows, bytes = res$bytes, df = res$df, samples = nrow(d),
                 outcomes = res$outcomes, variants = res$variants, pairs = res$pairs, covariate_rank = res$rank, block = res$block,
                 threads = res$threads, mem_gb = res$mem_gb, simd = res$simd, dosage = res$dosage,
                 active = res$active, products = res$products, waiting_most = res$waiting_most,
                 explicit_projection = res$explicit, monomorphic = res$monomorphic, singular = res$singular,
                 outcomes_not_tested = setdiff(oname, a$Outcome[res$tested]),
                 skipped_multiallelic = data.table::setnames(gv$skipped, "POS", "BP")))
}
