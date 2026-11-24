#' Genome-wide association of one trait (linear, logistic or Cox)
#'
#' Tests every selected variant of a PGEN or BED for association with one phenotype, adjusting for covariates:
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
#' Variants are not filtered here: give the variants to test in `snp_vec`.
#'
#' @param formula The null model, `phenotype ~ covariates`, as for [stats::glm()]: e.g. `y ~ age + sex + PC1`,
#'   `y ~ 1` without covariates; transformations, factors and interactions among covariates as in base R. Logistic:
#'   the phenotype is coded 0/1 (or 1/2, PLINK case/control). Every variable must be a column of `data`; `.` and
#'   `offset()` are not accepted.
#' @param data data.frame with one row per sample: `id_col` and the columns named in `formula`. Samples with a
#'   missing value in any of them, or not in the genotype file, are left out. A covariate that is not finite after
#'   its transformation (e.g. `log(0)`) in one of the analysed samples is an error.
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
#' @param snp_vec Variant IDs to test (`NULL` = all biallelic variants), tested in file order.
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
#'   differs only by using the null-model variance and the normal distribution. The attribute `info` then also
#'   holds `tau`, `s2e`, `h2`, `lrt`, `p_lrt`, `n_related`, `pairs`.
#' @param para With `kinship`: named list overriding [reml_tau_para()].
#' @param threads OpenMP threads; 0 takes them from [plink4r_resources()].
#' @param mem_gb Memory budget in GB that sets the block size; 0 takes it from [plink4r_resources()].
#' @param block Variants per thread and block; 0 = the largest that fits `mem_gb` (at most 4096). A thread expands
#'   its block of genotypes to doubles 2048 samples at a time, and a few hundred variants keep that within the
#'   cache: 256 took 0.4 to 0.6 of the time of 4096 for 1,000 and 20,000 samples (linear, logistic, kinship, Cox;
#'   0.5 to 0.65 for [gxe()]).
#' @return data.table with one row per variant: `CHR`, `POS`, `ID`, `A1` (counted allele: ALT of a PGEN, .bim
#'   column 5 of a BED), `A2`, `N` (non-missing calls), `A1_FREQ`, `BETA`, `SE`, `T` (linear) or `Z`
#'   (logistic), `P`, and `SE_ROBUST`, `P_ROBUST` when `robust`. Variants with no variation among the analysed
#'   samples (`monomorphic` in `info`) or a singular `Z'WZ` (`singular`) have `NA` statistics. Attribute `info`:
#'   samples, covariate rank, block size, threads, memory budget, SIMD set, variants recomputed by explicit
#'   projection, monomorphic and singular counts, skipped multiallelic variants.
#' @export
gwas <- function(formula, data, pfile = NULL, bfile = NULL, family = c("gaussian", "binomial", "cox"),
                 id_col = "IID", snp_vec = NULL, robust = FALSE, threads = 0L, mem_gb = 0, block = 256L,
                 kinship = NULL, para = list()) {
  family <- match.arg(family)
  if (family == "cox" && is.null(kinship))  # the Cox model has one route: its null model with no related pair
    kinship <- data.frame(ID1 = character(0), ID2 = character(0), kinship = numeric(0))
  if (!is.null(kinship)) {
    if (robust) stop("Robust standard errors are not available with kinship or for the Cox model.", call. = FALSE)
    return(.gwas_kinship(formula, data, pfile, bfile, kinship, family, id_col, snp_vec, threads, mem_gb, block,
                         para))
  }
  if (length(para)) stop("para is used with kinship only.", call. = FALSE)
  r <- .assoc_run(formula, data, pfile, bfile, exposure = NULL, int_covar = NULL, family, id_col, snp_vec,
                  robust, center = FALSE, threads, mem_gb, block)
  res <- r$res; out <- r$out
  stat <- if (family == "gaussian") "T" else "Z"
  out[, `:=`(BETA = res$beta[, 1L], SE = res$se[, 1L])]
  out[, (stat) := BETA / SE]
  out[, P := if (family == "gaussian") 2 * stats::pt(-abs(get(stat)), res$df) else 2 * stats::pnorm(-abs(get(stat)))]
  if (robust) {
    out[, SE_ROBUST := res$se_robust[, 1L]]
    out[, P_ROBUST := 2 * stats::pnorm(-abs(BETA / SE_ROBUST))]
  }
  data.table::setattr(out, "info", r$info)
  out[]
}

#' Gene-environment interaction of one trait (linear, logistic or Cox)
#'
#' For every selected variant fits `formula + g + g:E + g:int_covar` as GEM does, from one null-model fit
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
#' Reported per variant: the marginal effect from the model without interactions (`BETA_MARGINAL`,
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
#' @return data.table, one row per variant, as described above; attribute `info` as [gwas()], with the columns
#'   of the design used as `exposure_terms` and `int_covar_terms`.
#' @export
gxe <- function(formula, data, exposure, int_covar = NULL, pfile = NULL, bfile = NULL,
                family = c("gaussian", "binomial", "cox"), id_col = "IID", center = TRUE, robust = NULL,
                snp_vec = NULL, threads = 0L, mem_gb = 0, block = 256L, kinship = NULL, para = list()) {
  family <- match.arg(family)
  if (missing(exposure) || !length(exposure)) stop("Give at least one exposure.", call. = FALSE)
  if (family == "cox" && is.null(kinship))  # the Cox model has one route: its null model with no related pair
    kinship <- data.frame(ID1 = character(0), ID2 = character(0), kinship = numeric(0))
  if (!is.null(kinship)) {
    if (isTRUE(robust)) stop("Robust standard errors are not available with kinship or for the Cox model.", call. = FALSE)
    return(.gxe_kinship(formula, data, exposure, int_covar, pfile, bfile, kinship, family, id_col, center, snp_vec,
                        threads, mem_gb, block, para))
  }
  if (length(para)) stop("para is used with kinship only.", call. = FALSE)
  if (is.null(robust)) robust <- family == "gaussian"
  r <- .assoc_run(formula, data, pfile, bfile, exposure, int_covar, family, id_col, snp_vec, robust, center,
                  threads, mem_gb, block)
  res <- r$res; out <- r$out; terms <- r$terms; ne <- r$n_exp
  lin <- family == "gaussian"
  pz <- function(b, s) if (lin) 2 * stats::pt(-abs(b / s), res$df_marginal) else 2 * stats::pnorm(-abs(b / s))
  out[, `:=`(BETA_MARGINAL = res$beta_marginal, SE_MARGINAL = res$se_marginal)]
  out[, P_MARGINAL := pz(BETA_MARGINAL, SE_MARGINAL)]
  if (robust) {
    out[, SE_MARGINAL_ROBUST := res$se_marginal_robust]
    out[, P_MARGINAL_ROBUST := 2 * stats::pnorm(-abs(BETA_MARGINAL / SE_MARGINAL_ROBUST))]
  }
  nm <- c("G", paste0("GxE_", terms))
  for (q in seq_along(nm)) {
    out[, (paste0("BETA_", nm[q])) := res$beta[, q]]
    out[, (paste0("SE_", nm[q])) := res$se[, q]]
    if (robust) out[, (paste0("SE_ROBUST_", nm[q])) := res$se_robust[, q]]
  }
  out[, `:=`(P_INT = stats::pchisq(res$chisq_int, ne, lower.tail = FALSE),
             P_JOINT = stats::pchisq(res$chisq_joint, ne + 1L, lower.tail = FALSE))]
  if (robust)
    out[, `:=`(P_INT_ROBUST = stats::pchisq(res$chisq_int_robust, ne, lower.tail = FALSE),
               P_JOINT_ROBUST = stats::pchisq(res$chisq_joint_robust, ne + 1L, lower.tail = FALSE))]
  data.table::setattr(out, "info", c(r$info, list(exposure_terms = terms[seq_len(ne)],
                                                  int_covar_terms = terms[-seq_len(ne)])))
  out[]
}

# The genotype file of an association run, from geno_file(): its variants (var), sample IDs (iid), path, the file
# rows of the biallelic variants to test (vrow) and the IDs of the multiallelic ones left out (skipped).
.geno_variants <- function(pfile, bfile, snp_vec) {
  if (is.null(pfile) == is.null(bfile)) stop("Give exactly one of pfile and bfile.", call. = FALSE)
  gf <- geno_file(bedfile = bfile, pgenfile = pfile)
  v <- gf$variants
  var <- data.table::data.table(CHR = v$chr, POS = v$pos, ID = v$id, A1 = v$a1, A2 = v$a2, allele_ct = v$allele_ct)
  iid <- gf$samples$iid
  path <- gf$path
  vrow <- if (is.null(snp_vec)) seq_len(nrow(var)) else {
    i <- match(as.character(snp_vec), var$ID)
    if (anyNA(i)) warning(sum(is.na(i)), " of snp_vec not in the genotype file; dropped.", call. = FALSE)
    sort(unique(i[!is.na(i)]))
  }
  skipped <- var$ID[vrow[var$allele_ct[vrow] > 2L]]
  vrow <- vrow[var$allele_ct[vrow] <= 2L]
  if (!length(vrow)) stop("No biallelic variant to test.", call. = FALSE)
  list(var = var, iid = iid, path = path, vrow = vrow, skipped = skipped)
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

# Shared part of gwas() and gxe(): samples, design, null model, variants, and the C++ stream.
.assoc_run <- function(formula, data, pfile, bfile, exposure, int_covar, family, id_col, snp_vec, robust,
                       center, threads, mem_gb, block) {
  if (robust && family == "binomial") stop("Robust standard errors are for the linear model only.", call. = FALSE)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  gv <- .geno_variants(pfile, bfile, snp_vec)
  var <- gv$var; iid <- gv$iid; path <- gv$path; vrow <- gv$vrow; skipped <- gv$skipped
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

  res <- assoc_stream_cpp(path, is.null(pfile), length(iid), as.integer(var$allele_ct), as.integer(vrow - 1L),
                          as.integer(srow - 1L), X, S, length(ie), resid, w, robust, threads, mem_gb,
                          as.integer(block))
  out <- var[vrow, list(CHR, POS, ID, A1, A2)]
  out[, `:=`(N = as.integer(res$n_obs), A1_FREQ = res$a1_sum / (2 * res$n_obs))]
  info <- list(samples = n, covariate_rank = res$rank, null_converged = fit$converged, block = res$block,
               threads = res$threads, mem_gb = res$mem_gb, simd = res$simd, dosage = res$dosage,
               explicit_projection = res$explicit, monomorphic = sum(res$status == 1L),
               singular = sum(res$status == 2L),
               skipped_multiallelic = skipped)
  list(res = res, out = out, info = info, terms = terms, n_exp = length(ie))
}
