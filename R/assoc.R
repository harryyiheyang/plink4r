#' Genome-wide association of one trait (linear, logistic or Cox)
#'
#' Tests every biallelic variant of a PGEN or BED for association with one phenotype, adjusting for covariates:
#' `formula + g` with a linear (`family = "gaussian"`) or logistic (`"binomial"`) model. The genotype matrix is
#' never formed: the null model `formula` is fitted once with [stats::glm()] (which calls glm.fit), and its numeric design (`x = TRUE`: factors and characters expanded by base R with the current contrasts) is the covariate matrix used below; then variants are streamed
#' from the file in blocks whose size follows from `mem_gb` and reduced in C++ to a few numbers each (the
#' StreamBlock scheme of GEM 2: two products with fixed n x m matrices per block), so memory is about
#' `8 n (m1 + m2 + k + S)` bytes shared (the basis `V` and the fixed matrices `L1`, `L2`: `m1 = k S (1 + robust) + S`,
#' `m2 = S (S + 1) / 2 (1 + robust)` columns, `k` covariate rank, `S` = 1 + interaction terms; e.g. 0.7 GB for
#' n = 500,000, k = 20, S = 4, robust) plus, per thread, `n / 4` bytes per variant of the block (2-bit calls;
#' `2 n` for a PGEN with dosages) and two 2048-row double chunks.
#'
#' Linear: `BETA`, `SE`, `T` and `P` are those of `lm(y ~ covar + g)` (t distribution, `n - k - 1` degrees of
#' freedom, `k` the rank of the covariates with the intercept), with a missing call replaced by the variant's
#' mean among the analysed samples. `robust = TRUE` adds `SE_ROBUST` and `P_ROBUST` (normal), the sandwich
#' `A^-1 B A^-1` with `B` from the null-model residuals, as GEM. Logistic: one Newton step from the null model,
#' `BETA = U / I`, `SE = I^-1/2` (`U`, `I` the score and information of `g` given the covariates), `Z` and its
#' normal `P`; `Z^2` is the score test. No robust SE for the logistic model.
#'
#' Common variants: every test here is a large-sample test (normal or chi-square reference), meant for variants
#' with a minor allele frequency above 0.01. For the logistic and Cox models the reference also needs minor
#' alleles among the cases (events), about `2 x cases x MAF >= 30`, which MAF > 0.01 gives from 1,500 cases on.
#' Variants are not filtered here: every biallelic variant of the file is tested, in file order. To test a part of
#' them, write that part to a file of its own first (e.g. with PLINK).
#'
#' The result is written to the text file `out` as the scan goes; R holds nothing of it and reads no line of the
#' `.bim` / `.pvar`. That file is the skeleton of the result: every line of it is written as it is, in file order,
#' followed by the statistics of its variant (tab-separated), so the columns of the variant file (chromosome, ID,
#' position, alleles, and whatever else a `.pvar` holds) are in the result without having been read. The first line
#' names the columns: the `#CHROM` line of the `.pvar`, or `CHR ID CM POS A1 A2` for a `.bim` (which has none),
#' then the statistics: `N` (non-missing calls among the analysed samples), `A1_FREQ` (frequency of the counted
#' allele: ALT of a PGEN, .bim column 5 of a BED), `BETA`, `SE`, `T` (linear) or `Z` (logistic, Cox), `P`, and
#' `SE_ROBUST`, `P_ROBUST` when `robust`. A variant with no variation among the analysed samples or a singular
#' `Z'WZ` has `NA` after `A1_FREQ`; a line with more than two alleles (PGEN) is not tested and has `NA` throughout.
#' The file has one row per line of the `.bim` / `.pvar`, in its order, so any other table of the same variants can
#' be pasted beside it. Numbers are written with the digits that read back as the same double.
#'
#' @param formula The null model, `phenotype ~ covariates`, as for [stats::glm()]: e.g. `y ~ age + sex + PC1`,
#'   `y ~ 1` without covariates; transformations, factors and interactions among covariates as in base R. Logistic:
#'   the phenotype is coded 0/1 (or 1/2, PLINK case/control). Every variable must be a column of `data`; `.` and
#'   `offset()` are not accepted.
#' @param data data.frame with one row per sample: `id_col` and the columns named in `formula`. Samples with a
#'   missing value in any of them, or not in the genotype file, are left out. A covariate that is not finite after
#'   its transformation (e.g. `log(0)`) in one of the analysed samples is an error.
#' @param out Path of the text file to write (replaced when it exists); read it with `data.table::fread(out)`.
#' @param pfile PLINK2 prefix (.pgen/.pvar/.psam; .pvar may be .pvar.zst); or give `bfile`.
#' @param bfile PLINK1 prefix (.bed/.bim/.fam).
#' @param family `"gaussian"` (linear), `"binomial"` (logistic) or `"cox"` (proportional hazards). Cox: the
#'   phenotype of `formula` is `cbind(time, status)` (or `survival::Surv(time,
#'   status)`), `status` 1 for an event and 0 for censoring. The null model is fitted once and every variant is
#'   then tested one step from it (`BETA` a log hazard ratio; `Z` and its normal `P`). Without `kinship` (or with
#'   no related pair among the analysed samples) the null model is the Cox partial-likelihood fit with the Breslow
#'   baseline cumulative hazard `Lambda0`, and a variant is tested in the Poisson model of the event indicator
#'   with `log Lambda0(time)` as an offset. With `kinship` the null model is that of the life table of the data
#'   (binary records of the samples in `para$n_interval` intervals, complementary log-log link, interval
#'   intercepts as the baseline, the kinship correlation between the records of relatives; see [reml_tau()]), and
#'   a variant is tested in that model of the records.
#' @param id_col Column of `data` with the sample IDs (.psam / .fam IID).
#' @param robust Linear only: add sandwich standard errors (not with `kinship`).
#' @param kinship `NULL` (default: the samples are taken as unrelated), or the related pairs as a three-column
#'   data.frame: ID, ID, kinship coefficient (see [reml_tau()]). With `kinship` the null model accounts for the
#'   similarity of relatives through a sparse kinship matrix `K` and is fitted once ([reml_tau()]). Linear: the
#'   mixed model `formula + b`, `b ~ N(0, tau K)`, by REML (as fastGWA, Jiang et al. 2019). Logistic and Cox: the
#'   marginal model in which the Pearson residuals of relatives have correlation `h K_ij` (Cox: on the records of
#'   the life table). With `V` the covariance of the (working) response under that null model and
#'   `P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1` every variant gets the generalised least
#'   squares `BETA = g' P y / g' P g`, `SE = (g' P g)^-1/2` (exact, no GRAMMAR-GAMMA approximation), `T` or `Z` =
#'   `BETA / SE` and its normal `P`. Samples without a listed relative enter as in the ordinary test; with no
#'   related pair among the analysed samples the logistic result is that of `kinship = NULL`, and the linear one
#'   differs only by using the null-model variance and the normal distribution. The value then also holds `tau`,
#'   `s2e`, `h2`, `lrt`, `p_lrt`, `n_related`, `pairs`.
#' @param para With `kinship`: named list overriding [reml_tau_para()].
#' @param threads OpenMP threads; 0 takes them from [plink4r_resources()].
#' @param mem_gb Memory budget in GB that sets the block size; 0 takes it from [plink4r_resources()].
#' @param block Variants per thread and block; 0 = the largest that fits `mem_gb` (at most 4096). A thread expands
#'   its block of genotypes to doubles 2048 samples at a time, and a few hundred variants keep that within the
#'   cache: 256 took 0.4 to 0.6 of the time of 4096 for 1,000 and 20,000 samples (linear, logistic, kinship, Cox;
#'   0.5 to 0.65 for [gxe()]).
#' @return Invisibly, a list: `out` (the file), `variants` (its rows: the lines of the `.bim` / `.pvar`), `tested`
#'   (the biallelic ones), `bytes` (size of the file), and the run: samples, covariate rank, whether the null model
#'   converged, `df` (of the t distribution of the linear model), block size, threads, memory budget, SIMD set,
#'   whether dosages were read, variants recomputed by explicit projection, `monomorphic` and `singular` (variants
#'   without statistics) and `multiallelic` (lines not tested).
#' @export
gwas <- function(formula, data, out, pfile = NULL, bfile = NULL, family = c("gaussian", "binomial", "cox"),
                 id_col = "IID", robust = FALSE, threads = 0L, mem_gb = 0, block = 256L,
                 kinship = NULL, para = list()) {
  family <- match.arg(family)
  out <- .result_path(out)
  if (family == "cox" && is.null(kinship))  # the Cox model has one route: its null model with no related pair
    kinship <- data.frame(ID1 = character(0), ID2 = character(0), kinship = numeric(0))
  if (!is.null(kinship)) {
    if (robust) stop("Robust standard errors are not available with kinship or for the Cox model.", call. = FALSE)
    return(.gwas_kinship(formula, data, out, pfile, bfile, kinship, family, id_col, threads, mem_gb, block, para))
  }
  if (length(para)) stop("para is used with kinship only.", call. = FALSE)
  gf <- .geno_files(pfile, bfile)
  if (!is.null(gf$tmp)) on.exit(unlink(gf$tmp), add = TRUE)
  nl <- .assoc_null(formula, data, gf$iid, exposure = NULL, int_covar = NULL, family, id_col, robust, center = FALSE,
                    threads, mem_gb)
  res <- gwas_stream_cpp(gf$path, gf$bed, length(gf$iid), gf$text, out, nl$srow, nl$X, nl$resid, nl$w, robust,
                         nl$threads, nl$mem_gb, as.integer(block))
  invisible(.assoc_value(out, nl, res))
}

#' Gene-environment interaction of one trait (linear, logistic or Cox)
#'
#' For every biallelic variant fits `formula + g + g:E + g:int_covar` as GEM does, from one null-model fit
#' (`formula`, which contains the exposures and interaction covariates) and the same streaming as [gwas()]. The
#' null model is the one of [gwas()] called with the same `formula`, `family`, `kinship` and `para`: without
#' `kinship` the linear or logistic fit of [stats::glm()] (design taken from the fit); with `kinship`, and for the
#' Cox model, the null model described there, with the same kinship parameter for the association and the
#' interaction tests. Several exposures are tested together: `exposure` may
#' name any number of terms of `formula`, factors and character columns become indicator columns (one interaction
#' term each), and `int_covar` names terms whose interaction with `g` is adjusted for but not tested. `center = TRUE` centres the exposures and interaction
#' covariates (numeric ones and indicators alike), so `BETA_G` is the effect of `g` at their mean; the
#' interaction estimates do not change.
#'
#' Written per variant to `out`, after its line of the `.bim` / `.pvar` and `N`, `A1_FREQ` as [gwas()] writes them: the marginal effect from the model without interactions (`BETA_MARGINAL`,
#' `SE_MARGINAL`, `P_MARGINAL`, as [gwas()]); the main effect `BETA_G`, `SE_G` and every interaction
#' `BETA_GxE_<term>`, `SE_GxE_<term>` (terms of `exposure` and `int_covar`, after expansion); `P_INT`, the Wald test
#' of all exposure interactions (`df` = number of exposure terms), and `P_JOINT`, of `g` and those interactions
#' together (`df` + 1). Linear: model-based variances use the full-model residual variance (`lm`), `robust` (default
#' `TRUE` for linear, as GEM) adds `SE_ROBUST_*`, `P_INT_ROBUST`, `P_JOINT_ROBUST` and `SE_MARGINAL_ROBUST`.
#' Logistic: one Newton step from the null model (score-based), no robust SE. Chi-square tests throughout.
#'
#' With `kinship`, and for the Cox model: with `V` the covariance of the (working) response under the null model,
#' `P` as in [gwas()] and `Z = [g, g e_1, ..., g e_{S-1}]` (the exposure terms, then the interaction covariates), a
#' variant gets `c = Z' P y` and `A = Z' P Z`. The marginal effect is `c_1 / A_11` with `SE = A_11^-1/2`, the test
#' of [gwas()]; `A^-1 c` with covariance `A^-1` gives `BETA_G`, the interactions and the Wald tests `P_INT` and
#' `P_JOINT` (generalised least squares for the linear model, one step from the null model for the logistic and
#' Cox models; normal and chi-square references, model-based variances, no robust SE). The exposures of a Cox
#' model are fixed in time.
#'
#' Common variants, as [gwas()]: use variants with a minor allele frequency above 0.01.
#'
#' @inheritParams gwas
#' @param exposure Names of the exposures (one or more): terms of `formula` as written there, e.g. `"smoke"` for
#'   `y ~ age + smoke`, `"log(bmi)"` for `y ~ age + log(bmi)`.
#' @param int_covar Names of interaction covariates (adjusted for, not tested), terms of `formula`, or `NULL`.
#' @param center Centre the exposures and interaction covariates.
#' @param robust Linear model without `kinship` only: add sandwich standard errors (default `TRUE` there).
#' @return Invisibly, a list as that of [gwas()], with `df_marginal` (of the t distribution of `P_MARGINAL`, linear
#'   model) and the columns of the design used as `exposure_terms` and `int_covar_terms`.
#' @export
gxe <- function(formula, data, exposure, out, int_covar = NULL, pfile = NULL, bfile = NULL,
                family = c("gaussian", "binomial", "cox"), id_col = "IID", center = TRUE, robust = NULL,
                threads = 0L, mem_gb = 0, block = 256L, kinship = NULL, para = list()) {
  family <- match.arg(family)
  if (missing(exposure) || !length(exposure)) stop("Give at least one exposure.", call. = FALSE)
  out <- .result_path(out)
  if (family == "cox" && is.null(kinship))  # the Cox model has one route: its null model with no related pair
    kinship <- data.frame(ID1 = character(0), ID2 = character(0), kinship = numeric(0))
  if (!is.null(kinship)) {
    if (isTRUE(robust)) stop("Robust standard errors are not available with kinship or for the Cox model.", call. = FALSE)
    return(.gxe_kinship(formula, data, exposure, out, int_covar, pfile, bfile, kinship, family, id_col, center,
                        threads, mem_gb, block, para))
  }
  if (length(para)) stop("para is used with kinship only.", call. = FALSE)
  if (is.null(robust)) robust <- family == "gaussian"
  gf <- .geno_files(pfile, bfile)
  if (!is.null(gf$tmp)) on.exit(unlink(gf$tmp), add = TRUE)
  nl <- .assoc_null(formula, data, gf$iid, exposure, int_covar, family, id_col, robust, center, threads, mem_gb)
  res <- gxe_stream_cpp(gf$path, gf$bed, length(gf$iid), gf$text, out, nl$srow, nl$X, nl$S, nl$n_exp, nl$terms,
                        nl$resid, nl$w, robust, nl$threads, nl$mem_gb, as.integer(block))
  invisible(c(.assoc_value(out, nl, res),
              list(df_marginal = res$df_marginal, exposure_terms = nl$terms[seq_len(nl$n_exp)],
                   int_covar_terms = nl$terms[-seq_len(nl$n_exp)])))
}

# The path of the result file of gwas() / gxe(), checked.
.result_path <- function(out) {
  if (missing(out) || !is.character(out) || length(out) != 1L || is.na(out) || !nzchar(out))
    stop("out must be the path of the text file to write.", call. = FALSE)
  if (!dir.exists(dirname(out))) stop("No such directory: ", dirname(out), call. = FALSE)
  path.expand(out)
}

# The value of gwas() / gxe() without kinship: the file and the run (nl: .assoc_null, res: the C++ stream).
.assoc_value <- function(out, nl, res) {
  list(out = out, variants = res$variants, tested = res$tested, bytes = res$bytes, samples = nl$n,
       covariate_rank = res$rank, null_converged = nl$converged, df = res$df, block = res$block,
       threads = res$threads, mem_gb = res$mem_gb, simd = res$simd, dosage = res$dosage,
       explicit_projection = res$explicit, monomorphic = res$monomorphic, singular = res$singular,
       multiallelic = res$variants - res$tested)
}

# The files of an association run: the genotype file (path; bed = whether a BED), its variant text, which the compiled
# code reads (text: the .bim, or the .pvar; a .pvar.zst is expanded into a temporary file, tmp, which the caller
# removes), and the sample IDs in file order (iid). R reads no line of the variant text.
.geno_files <- function(pfile, bfile) {
  if (is.null(pfile) == is.null(bfile)) stop("Give exactly one of pfile and bfile.", call. = FALSE)
  if (!is.null(bfile)) {
    bed <- .geno_bed_path(bfile)
    base <- sub("\\.bed$", "", bed)
    for (f in paste0(base, c(".bim", ".fam"))) if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
    fam <- data.table::fread(paste0(base, ".fam"), header = FALSE, select = 2L, colClasses = "character")
    return(list(path = bed, bed = TRUE, text = paste0(base, ".bim"), tmp = NULL, iid = fam[[1L]]))
  }
  files <- .geno_pgen_files(pfile)
  pvar <- paste0(files$prefix, ".pvar")
  tmp <- NULL
  if (!file.exists(pvar)) {
    if (!file.exists(paste0(pvar, ".zst"))) stop("no .pvar or .pvar.zst for ", files$prefix, call. = FALSE)
    pvar <- tmp <- .zst_expand(paste0(pvar, ".zst"))
  }
  list(path = files$pgen, bed = FALSE, text = normalizePath(pvar, winslash = "/"), tmp = tmp,
       iid = .psam_ids(files$prefix))
}

# The analysis samples of the null model `formula` (phenotype ~ covariates): rows of d with an ID and no missing
# value in the variables of the formula (and, with iid, present in the genotype file, in file order; srow = their
# file rows), a finite phenotype after any transformation, and that phenotype (0/1 for a logistic model; the
# follow-up time, with status the 0/1 event indicator, for a Cox model). Every route takes its samples here, so
# what no route handles is refused here: offset() in the formula, and a covariate that is not finite after its
# transformation in an analysed sample.
.pheno_samples <- function(d, formula, family, id_col, iid = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("formula must be two-sided: phenotype ~ covariates.", call. = FALSE)
  vars <- all.vars(formula)
  if ("." %in% vars) stop("Write the covariates out in the formula ('.' is not accepted).", call. = FALSE)
  tm <- stats::terms(formula)
  if (attr(tm, "intercept") != 1L) stop("The formula must keep its intercept.", call. = FALSE)
  if (!is.null(attr(tm, "offset"))) stop("offset() is not accepted in the formula.", call. = FALSE)
  miss <- setdiff(c(id_col, vars), names(d))
  if (length(miss)) stop("data lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  ids <- as.character(d[[id_col]])
  if (anyDuplicated(ids[!is.na(ids)])) stop("IDs repeated in data.", call. = FALSE)
  ok <- !is.na(ids) & stats::complete.cases(d[, vars, drop = FALSE])
  srow <- NULL
  if (!is.null(iid)) {
    srow <- match(ids, iid)
    ok <- ok & !is.na(srow)
  }
  y <- eval(formula[[2L]], d, environment(formula))
  status <- NULL
  if (family == "cox") {  # time and event indicator: cbind(time, status), or survival::Surv(time, status)
    if (!is.matrix(y) || ncol(y) != 2L || nrow(y) != nrow(d) || !(is.numeric(y) || is.logical(y)))
      stop("A Cox model takes cbind(time, status) (or Surv(time, status)) as its phenotype.", call. = FALSE)
    status <- as.numeric(y[, 2L]); y <- as.numeric(y[, 1L])
    ok <- ok & is.finite(status)
  }
  if ((!is.numeric(y) && !is.logical(y)) || length(y) != nrow(d)) stop("The phenotype must be numeric.", call. = FALSE)
  y <- as.numeric(y)
  ok <- ok & is.finite(y)
  if (!any(ok)) stop("No sample of data has complete data", if (!is.null(iid)) " and is in the genotype file", ".",
                     call. = FALSE)
  d <- d[ok, , drop = FALSE]; y <- y[ok]; status <- status[ok]
  if (!is.null(iid)) {
    srow <- srow[ok]
    o <- order(srow); d <- d[o, , drop = FALSE]; srow <- srow[o]; y <- y[o]; status <- status[o]
  }
  if (family == "binomial") {
    u <- sort(unique(y))
    if (all(u %in% c(1, 2))) y <- y - 1 else if (!all(u %in% c(0, 1)))
      stop("A logistic phenotype must be coded 0/1 or 1/2.", call. = FALSE)
    if (length(unique(y)) < 2L) stop("The phenotype has one class only.", call. = FALSE)
  }
  if (family == "cox") {
    if (!all(status %in% c(0, 1))) stop("The event indicator of a Cox model must be coded 0/1.", call. = FALSE)
    if (!any(status == 1)) stop("No event among the analysed samples.", call. = FALSE)
    if (any(y < 0)) stop("Follow-up times must not be negative.", call. = FALSE)
  }
  for (cl in all.vars(formula[[3L]]))
    if ((is.factor(d[[cl]]) || is.character(d[[cl]])) && length(unique(d[[cl]])) < 2L)
      stop("Column ", cl, " has a single level among the analysed samples.", call. = FALSE)
  mf <- tryCatch(stats::model.frame(stats::delete.response(tm), d, na.action = stats::na.pass),
                 error = function(e) stop("A term of the formula cannot be evaluated on the analysed samples: ",
                                          conditionMessage(e), call. = FALSE))
  bad <- names(mf)[!vapply(mf, function(v) !is.numeric(v) || all(is.finite(v)), NA)]
  if (length(bad))
    stop("Not finite in an analysed sample (after its transformation): ", paste(bad, collapse = ", "), ".",
         call. = FALSE)
  list(d = d, srow = srow, y = y, status = status)
}

# Shared part of gwas() and gxe() without kinship: the analysis samples (iid: the sample IDs of the genotype file),
# the design and the null model, as the C++ streams take them: X (columns: intercept, exposures, interaction
# covariates, the rest), S (1 + interaction columns), n_exp, the names of the interaction columns (terms), the
# null residual, the weights of a logistic null, and the 0-based file rows of the samples (srow).
.assoc_null <- function(formula, data, iid, exposure, int_covar, family, id_col, robust, center, threads, mem_gb) {
  if (robust && family == "binomial") stop("Robust standard errors are for the linear model only.", call. = FALSE)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  sm <- .pheno_samples(as.data.frame(data), formula, family, id_col, iid)
  d <- sm$d; srow <- sm$srow; y <- sm$y
  n <- nrow(d)

  # ---- null model with base glm() on the user's formula: factors and characters are expanded by R's own
  # model.matrix (contrasts as options("contrasts")), and the numeric design is taken from the fit (x = TRUE).
  # Columns, reordered: intercept, exposures, interaction covariates (the S interaction columns of Z), the rest.
  d[[".plink4r_y"]] <- y  # the phenotype as analysed (transformed, 0/1)
  fml <- stats::update(formula, .plink4r_y ~ .)
  fam <- if (family == "gaussian") stats::gaussian() else stats::binomial()
  fit <- stats::glm(fml, family = fam, data = d, x = TRUE)
  if (!fit$converged) warning("The null model did not converge.", call. = FALSE)
  X <- fit$x
  asg <- attr(X, "assign")
  lab <- gsub("`", "", attr(stats::terms(fit), "term.labels"))
  nt <- setdiff(gsub("`", "", c(exposure, int_covar)), lab)
  if (length(nt))
    stop("Not a term of the formula: ", paste(nt, collapse = ", "), ". Its terms: ", paste(lab, collapse = ", "), ".",
         call. = FALSE)
  cols_of <- function(v) which(asg %in% match(gsub("`", "", v), lab))
  ie <- cols_of(exposure); ii <- cols_of(int_covar)
  X <- X[, c(1L, ie, ii, setdiff(seq_len(ncol(X)), c(1L, ie, ii))), drop = FALSE]
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL
  S <- 1L + length(ie) + length(ii)
  if (S > 1L) {
    E <- X[, 2:S, drop = FALSE]
    if (any(apply(E, 2L, stats::var) == 0))
      stop("An exposure or interaction covariate is constant among the analysed samples.", call. = FALSE)
    # centring changes the parametrisation only (the null fit and its residuals are the same)
    if (center) X[, 2:S] <- sweep(E, 2L, colMeans(E))
  }
  terms <- colnames(X)[seq_len(S)[-1L]]
  mu <- fit$fitted.values
  resid <- y - mu
  w <- if (family == "binomial") mu * (1 - mu) else numeric(0)

  list(X = X, S = S, n_exp = length(ie), terms = terms, resid = resid, w = w, srow = as.integer(srow - 1L), n = n,
       converged = fit$converged, threads = threads, mem_gb = mem_gb)
}
