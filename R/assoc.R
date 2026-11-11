#' Genome-wide association of one trait (linear or logistic)
#'
#' Tests every selected variant of a PGEN or BED for association with one phenotype, adjusting for covariates:
#' `y ~ covar + g` with a linear (`family = "gaussian"`) or logistic (`"binomial"`) model. The genotype matrix is
#' never formed: the null model `y ~ covar` is fitted once with [stats::glm()] (which calls glm.fit), and its numeric design (`x = TRUE`: factors and characters expanded by base R with the current contrasts) is the covariate matrix used below; then variants are streamed
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
#' @param pfile PLINK2 prefix (.pgen/.pvar/.psam; .pvar may be .pvar.zst); or give `bfile`.
#' @param bfile PLINK1 prefix (.bed/.bim/.fam).
#' @param data data.frame with one row per sample: `id_col`, the phenotype and covariate columns.
#' @param pheno Name of the phenotype column. Logistic: coded 0/1 (or 1/2, PLINK case/control).
#' @param covar Names of covariate columns (numeric, or factors / characters, expanded by base glm()). An
#'   intercept is always added.
#' @param family `"gaussian"` (linear) or `"binomial"` (logistic).
#' @param id_col Column of `data` with the sample IDs (.psam / .fam IID).
#' @param snp_vec Variant IDs to test (`NULL` = all biallelic variants), tested in file order.
#' @param robust Linear only: add sandwich standard errors.
#' @param threads OpenMP threads; 0 takes them from [plink4r_resources()].
#' @param mem_gb Memory budget in GB that sets the block size; 0 takes it from [plink4r_resources()].
#' @param block Variants per thread and block; 0 = the largest that fits `mem_gb` (at most 4096).
#' @return data.table with one row per variant: `CHR`, `POS`, `ID`, `A1` (counted allele: ALT of a PGEN, .bim
#'   column 5 of a BED), `A2`, `N` (non-missing calls), `A1_FREQ`, `BETA`, `SE`, `T` (linear) or `Z`
#'   (logistic), `P`, and `SE_ROBUST`, `P_ROBUST` when `robust`. Variants with no variation among the analysed
#'   samples (`monomorphic` in `info`) or a singular `Z'WZ` (`singular`) have `NA` statistics. Attribute `info`:
#'   samples, covariate rank, block size, threads, memory budget, SIMD set, variants recomputed by explicit
#'   projection, monomorphic and singular counts, skipped multiallelic variants.
#' @export
gwas <- function(pfile = NULL, bfile = NULL, data, pheno, covar = NULL, family = c("gaussian", "binomial"),
                 id_col = "IID", snp_vec = NULL, robust = FALSE, threads = 0L, mem_gb = 0, block = 0L) {
  family <- match.arg(family)
  r <- .assoc_run(pfile, bfile, data, pheno, covar, exposure = NULL, int_covar = NULL, family, id_col, snp_vec,
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

#' Gene-environment interaction of one trait (GEM, linear or logistic)
#'
#' For every selected variant fits `y ~ covar + E + int_covar + g + g:E + g:int_covar` (linear or logistic) as GEM
#' does, from one null-model fit (`y ~ covar + E + int_covar`, [stats::glm()], design taken from the fit) and the same streaming as
#' [gwas()]. Several exposures are tested together: `exposure` may name any number of columns, factors and
#' character columns become indicator columns (one interaction term each), and `int_covar` names covariates whose
#' interaction with `g` is adjusted for but not tested. `center = TRUE` centres the exposures and interaction
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
#' @inheritParams gwas
#' @param exposure Names of the exposure columns of `data` (one or more).
#' @param int_covar Names of interaction covariates (adjusted for, not tested), or `NULL`.
#' @param center Centre the exposures and interaction covariates.
#' @param robust Linear only: add sandwich standard errors (default `TRUE` for linear).
#' @return data.table, one row per variant, as described above; attribute `info` as [gwas()].
#' @export
gxe <- function(pfile = NULL, bfile = NULL, data, pheno, exposure, int_covar = NULL, covar = NULL,
                family = c("gaussian", "binomial"), id_col = "IID", center = TRUE, robust = NULL,
                snp_vec = NULL, threads = 0L, mem_gb = 0, block = 0L) {
  family <- match.arg(family)
  if (is.null(robust)) robust <- family == "gaussian"
  if (missing(exposure) || !length(exposure)) stop("Give at least one exposure.", call. = FALSE)
  r <- .assoc_run(pfile, bfile, data, pheno, covar, exposure, int_covar, family, id_col, snp_vec, robust, center,
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

# Shared part of gwas() and gxe(): samples, design, null model, variants, and the C++ stream.
.assoc_run <- function(pfile, bfile, data, pheno, covar, exposure, int_covar, family, id_col, snp_vec, robust,
                       center, threads, mem_gb, block) {
  if (is.null(pfile) == is.null(bfile)) stop("Give exactly one of pfile and bfile.", call. = FALSE)
  if (robust && family == "binomial") stop("Robust standard errors are for the linear model only.", call. = FALSE)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  d <- as.data.frame(data)
  need <- c(id_col, pheno, covar, exposure, int_covar)
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("data lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (length(pheno) != 1L) stop("Give one phenotype column.", call. = FALSE)

  # ---- genotype file
  if (!is.null(pfile)) {
    pvar <- .read_pvar5(pfile)
    var <- data.table::data.table(CHR = pvar$CHROM, POS = pvar$POS, ID = pvar$ID, A1 = pvar$ALT, A2 = pvar$REF,
                                  allele_ct = pvar$allele_ct)
    iid <- .psam_ids(pfile)
    path <- normalizePath(paste0(pfile, ".pgen"))
  } else {
    bim <- data.table::fread(paste0(bfile, ".bim"), header = FALSE, colClasses = "character")
    var <- data.table::data.table(CHR = sub("^chr", "", bim[[1L]], ignore.case = TRUE), POS = as.integer(bim[[4L]]),
                                  ID = bim[[2L]], A1 = toupper(bim[[5L]]), A2 = toupper(bim[[6L]]), allele_ct = 2L)
    iid <- data.table::fread(paste0(bfile, ".fam"), header = FALSE, colClasses = "character")[[2L]]
    path <- normalizePath(paste0(bfile, ".bed"))
  }
  vrow <- if (is.null(snp_vec)) seq_len(nrow(var)) else {
    i <- match(as.character(snp_vec), var$ID)
    if (anyNA(i)) warning(sum(is.na(i)), " of snp_vec not in the genotype file; dropped.", call. = FALSE)
    sort(unique(i[!is.na(i)]))
  }
  skipped <- var$ID[vrow[var$allele_ct[vrow] > 2L]]
  vrow <- vrow[var$allele_ct[vrow] <= 2L]
  if (!length(vrow)) stop("No biallelic variant to test.", call. = FALSE)

  # ---- samples: in the genotype file, complete phenotype / covariates / exposures; file order
  ids <- as.character(d[[id_col]])
  if (anyDuplicated(ids[!is.na(ids)])) stop("IDs repeated in data.", call. = FALSE)
  srow <- match(ids, iid)
  ok <- !is.na(srow) & stats::complete.cases(d[, setdiff(need, id_col), drop = FALSE])
  if (!any(ok)) stop("No sample of data is in the genotype file with complete data.", call. = FALSE)
  d <- d[ok, , drop = FALSE]; srow <- srow[ok]
  o <- order(srow); d <- d[o, , drop = FALSE]; srow <- srow[o]
  n <- nrow(d)

  y <- d[[pheno]]
  if (!is.numeric(y) && !is.logical(y)) stop("The phenotype must be numeric.", call. = FALSE)
  y <- as.numeric(y)
  if (family == "binomial") {
    u <- sort(unique(y))
    if (all(u %in% c(1, 2))) y <- y - 1 else if (!all(u %in% c(0, 1)))
      stop("A logistic phenotype must be coded 0/1 or 1/2.", call. = FALSE)
    if (length(unique(y)) < 2L) stop("The phenotype has one class only.", call. = FALSE)
  }

  # ---- null model with base glm(): factors and characters are expanded by R's own model.matrix (contrasts
  # as options("contrasts")), and the numeric design is taken from the fit (x = TRUE). Columns: intercept,
  # exposures, interaction covariates (the S interaction columns of Z), then the other covariates.
  for (cl in c(exposure, int_covar, covar))
    if ((is.factor(d[[cl]]) || is.character(d[[cl]])) && length(unique(d[[cl]])) < 2L)
      stop("Column ", cl, " has a single level among the analysed samples.", call. = FALSE)
  d[[pheno]] <- y
  rhs <- c(exposure, int_covar, covar)
  fml <- stats::reformulate(if (length(rhs)) sprintf("`%s`", rhs) else "1", response = as.name(pheno))
  fam <- if (family == "gaussian") stats::gaussian() else stats::binomial()
  fit <- stats::glm(fml, family = fam, data = d, x = TRUE)
  if (!fit$converged) warning("The null model did not converge.", call. = FALSE)
  X <- fit$x
  asg <- attr(X, "assign")
  lab <- gsub("`", "", attr(stats::terms(fit), "term.labels"))
  cols_of <- function(v) which(asg %in% match(v, lab))
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
