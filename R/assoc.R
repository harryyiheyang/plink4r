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
#' mean among the analysed samples. `robust = TRUE` adds `SE_ROBUST` and `P_ROBUST`, the sandwich
#' `A^-1 B A^-1` with `B` from the null-model residuals, as GEM. With `mac_adjust = TRUE` (default) `P_ROBUST` is
#' the small-count version: `T^2 = (BETA / SE_ROBUST)^2` is mapped to `H = (M - 1) T^2 / (M - T^2)` and referred
#' to `F(1, M - 1)`, where `M = N (1 - (1 - p)^2)` is the expected number of carriers under Hardy-Weinberg
#' (`p` the minor allele frequency; hard calls and dosages alike), and `P_ROBUST` is the larger of that P and the
#' normal one. This removes the inflated median of the normal P at low minor-allele counts without making any P
#' smaller than before (the tail keeps the normal P, which is conservative there); when the minor allele count
#' `2 N p` is below 5, or with `mac_adjust = FALSE`, `P_ROBUST` is the normal P. Logistic: one Newton step from the null model,
#' `BETA = U / I`, `SE = I^-1/2` (`U`, `I` the score and information of `g` given the covariates), `Z = BETA / SE`
#' (`Z^2` is the score test) and `P`. For variants whose normal P is below `spa_cut`, `P` is the saddlepoint P of
#' the score (fastSPA: the cumulant generating function is exact over the carriers of the minor allele, dosage at
#' least 0.05, and normal over the other samples), which keeps rare variants and unbalanced case-control ratios
#' calibrated; elsewhere it is the normal P. No robust SE for the logistic model. `BETA` and `SE` are not penalised
#' (use them for meta-analysis); with `shrink = TRUE` (default) the logistic result also has `BETA_POST` and
#' `SE_POST`, the posterior mean and sd of the per-allele log-OR under a mixture prior (null, three normal scales
#' that grow as the minor allele frequency falls, and a wide normal for large effects; see `prior`). For variants
#' with at most `carrier_max` expected carriers the likelihood is first refined over the carriers alone (penalised
#' Newton from the one-step ridge, to convergence), so large effects of rare variants are not overestimated by the
#' one-step score. The posterior is shrunk towards 0 by design and must not be meta-analysed.
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
#' @param mac_adjust Refer one-degree-of-freedom robust tests to `F(1, M - 1)` with `M` the expected carrier
#'   count (see Details); `FALSE` keeps the normal (chi-square) P.
#' @param shrink Logistic only: add `BETA_POST` and `SE_POST`.
#' @param prior Logistic only: list with any of `sigma0sq` (scale: the per-allele prior variance of component k
#'   is `sigma0sq * grid[k] * (2 p (1 - p))^alpha`), `weights` (of the `grid` components), `alpha` (default -0.38),
#'   `grid` (default `c(0, 1, 10, 100)`), `wide_sd` (default 1.5, log-OR, not scaled by MAF) and `wide_weight`
#'   (default 0.02). `sigma0sq` and `weights` not given are learned from the common variants (MAF >= 1%) of the
#'   same run: `sigma0sq` from the strongest variant of each LD block with P < 5e-8 (normal likelihood of the z
#'   values truncated at the threshold, which removes the winner's curse; the middle grid component is put at the
#'   hits' variance), the weights by EM over one random common variant per block. With fewer than 10 hits or 50
#'   blocks the corresponding defaults are kept (sd 0.05 at MAF 0.25; weights 0.7, 0.15, 0.1, 0.05); the prior used
#'   is in `info$prior`. The carrier step uses the ridge penalty `1 / wide_sd^2`.
#' @param robust_min_mac With `mac_adjust`, `P_ROBUST` is the model-based `P` for variants whose minor allele count
#'   `2 N p` is below this (default 5): the sandwich variance carries no information there (real data: lambda 2.2
#'   at MAC 1-5 against 1.08 for the model P). Those rows are recognisable from `N` and `A1_FREQ`.
#' @param min_mac `P` (and `P_ROBUST`) are `NA` for variants whose minor allele count `2 N p` is below this (e.g.
#'   dosage-only variants with MAC < 1, where neither the normal nor the saddlepoint reference is reliable); 0 =
#'   never.
#' @param blocks Logistic with `shrink`: LD blocks for prior learning, 38 or 37 ([ld_blocks_default()]) or a block
#'   file (header, then chr start stop).
#' @param carrier_max Logistic only: refine the likelihood over the carriers for variants with at most this many
#'   expected carriers (0 = never).
#' @param spa_cut Logistic only: saddlepoint P for the variants whose normal P is below this; 0 = normal P
#'   throughout.
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
                 id_col = "IID", snp_vec = NULL, robust = FALSE, threads = 0L, mem_gb = 0, block = 0L,
                 mac_adjust = TRUE, spa_cut = 0.05, shrink = TRUE, prior = NULL, carrier_max = 2000,
                 blocks = 38, min_mac = 1, robust_min_mac = 5) {
  family <- match.arg(family)
  shrink <- family == "binomial" && shrink
  pr <- if (shrink) .rv_prior(prior) else NULL
  r <- .assoc_run(pfile, bfile, data, pheno, covar, exposure = NULL, int_covar = NULL, family, id_col, snp_vec,
                  robust, center = FALSE, threads, mem_gb, block, spa_cut,
                  carrier_max = if (shrink) carrier_max else 0, carrier_lambda = if (shrink) 1 / pr$wide_sd^2 else 1)
  res <- r$res; out <- r$out
  stat <- if (family == "gaussian") "T" else "Z"
  out[, `:=`(BETA = res$beta[, 1L], SE = res$se[, 1L])]
  out[, (stat) := BETA / SE]
  out[, P := if (family == "gaussian") 2 * stats::pt(-abs(get(stat)), res$df) else
    data.table::fifelse(is.na(res$p_spa), 2 * stats::pnorm(-abs(get(stat))), res$p_spa)]
  if (shrink) {
    pr <- .rv_learn(out, pr, names(prior), blocks)
    .rv_shrink(out, res, pr)
    r$info$prior <- pr
  }
  if (robust) {
    out[, SE_ROBUST := res$se_robust[, 1L]]
    out[, P_ROBUST := .p_robust_1df((BETA / SE_ROBUST)^2, N, A1_FREQ, mac_adjust, P, robust_min_mac)]
  }
  if (min_mac > 0) {
    lowp <- which(!(2 * out$N * pmin(out$A1_FREQ, 1 - out$A1_FREQ) >= min_mac))
    if (length(lowp)) for (cl in grep("^P($|_)", names(out), value = TRUE)) data.table::set(out, lowp, cl, NA_real_)
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
#' `BETA_GxE_<term>`, `SE_GxE_<term>` (terms of `exposure` and `int_covar`, after expansion); `P_INT`, the Wald test of
#' all exposure interactions (`df` = number of exposure terms; with one term and a logistic model, the saddlepoint
#' P of its efficient score below `spa_cut`), and `P_JOINT`, of `g` and those interactions
#' together (`df` + 1). Linear: model-based variances use the full-model residual variance (`lm`), `robust` (default
#' `TRUE` for linear, as GEM) adds `SE_ROBUST_*`, `P_INT_ROBUST`, `P_JOINT_ROBUST` and `SE_MARGINAL_ROBUST`.
#' With two or more exposure terms, `P_INT_CAUCHY` (and `P_INT_ROBUST_CAUCHY`) is an additional test: the Cauchy
#' combination of the single-term P values (half of the weight, shared equally; t, the small-count robust
#' reference, or the logistic saddlepoint P) and a variance-component test (the other half: the sum of the squared
#' single-term z values, weighted chi-square null by saddlepoint). It gains power when the effect is spread over
#' terms but loses it when the exposure terms are strongly correlated (e.g. BMI and BMI^2), where the Wald `P_INT`
#' is the better test; both are reported.
#' Logistic: one Newton step from the null model (score-based), no robust SE; `P_MARGINAL` is the saddlepoint P of
#' [gwas()] below `spa_cut`. With `shrink = TRUE` (default) the logistic result also has `BETA_POST_G`, `SE_POST_G`,
#' `BETA_POST_GxE_<term>` and `SE_POST_GxE_<term>`: the joint posterior of the S coefficients, each pair of a
#' main-effect component (the mixture of [gwas()], learned from the marginal results of the same run unless given)
#' and an interaction component (null; variance `int_ratio * sigma0sq * (2p(1-p))^alpha` with weight `int_weight`,
#' default 10 and 0.015; the wide N(0, wide_sd^2) with weight `int_wide_weight`, default 0.005) giving a prior
#' `N(0, D C D)`, where `D` holds the component sds and `C = I / (tr(I) / S)` is the shape of the variant's
#' information matrix (its eigenvalues raised to `gamma`, default 1). The shape keeps the prior invariant to how
#' the exposures are coded and shrinks directions the carriers cannot separate; its size still follows the MAF.
#' Variants with at most `carrier_max` expected carriers get the S-dimensional carrier step first. Not for
#' meta-analysis. Chi-square tests throughout, except
#' that with `mac_adjust = TRUE` the one-degree-of-freedom robust tests (`P_MARGINAL_ROBUST`, and `P_INT_ROBUST`
#' when there is one exposure term) use the small-count `F(1, M - 1)` reference of [gwas()]; tests with two or more
#' degrees of freedom keep the chi-square.
#'
#' @inheritParams gwas
#' @param exposure Names of the exposure columns of `data` (one or more).
#' @param int_covar Names of interaction covariates (adjusted for, not tested), or `NULL`.
#' @param center Centre the exposures and interaction covariates.
#' @param robust Linear only: add sandwich standard errors (default `TRUE` for linear).
#' @param min_mac Every P column is `NA` for variants whose minor allele count `2 N p` is below this (e.g.
#'   dosage-only variants with MAC < 1, where no reference distribution is reliable); 0 = never.
#' @param min_mac_int Interaction tests are `NA` below this minor allele count: for the linear model all of them
#'   (default 10), for the logistic model the Wald chi-square ones that have no saddlepoint version (`P_INT` with two
#'   or more exposure terms, `P_JOINT`; default 100, as they were seen at 1e-12 to 1e-52 under the null up to MAC
#'   1000 in real data). **With a logistic model and two or more exposure terms, use `P_INT_CAUCHY`**, which is
#'   built from saddlepoint P values and stays calibrated at low counts.
#'   Note that the logistic interaction estimates are one Newton step from the null: at moderate counts a Wald
#'   chi-square can still be extreme (real data: MAC 161, `P_INT` 2.8e-12 with `BETA_GxE` -41 and a posterior near
#'   0). For logistic interactions trust `P_INT_CAUCHY` and the posterior columns over the Wald `P_INT` /
#'   `P_JOINT`.
#' @param min_mac_int_spa Logistic only: floor for the saddlepoint-based interaction tests (`P_INT` with one term,
#'   `P_INT_CAUCHY`), default 10 (real data: P down to 6e-14 at MAC 1-5).
#' @param robust_min_mac As [gwas()], for `P_MARGINAL_ROBUST` and the one-df `P_INT_ROBUST`.
#' @param eig_tol Interaction and joint tests are `NA` when the information matrix of G and its interactions,
#'   scaled to unit diagonal, has its smallest eigenvalue below this (near singular: the carriers barely vary in
#'   the exposures); `info$int_singular` counts them.
#' @param shrink,prior,carrier_max,blocks Logistic only, as [gwas()]; `prior` also takes `int_ratio`, `int_weight`,
#'   `int_wide_weight` and `gamma` (see Details).
#' @return data.table, one row per variant, as described above; attribute `info` as [gwas()].
#' @export
gxe <- function(pfile = NULL, bfile = NULL, data, pheno, exposure, int_covar = NULL, covar = NULL,
                family = c("gaussian", "binomial"), id_col = "IID", center = TRUE, robust = NULL,
                snp_vec = NULL, threads = 0L, mem_gb = 0, block = 0L, mac_adjust = TRUE, spa_cut = 0.05,
                shrink = TRUE, prior = NULL, carrier_max = 2000, blocks = 38, min_mac = 1, min_mac_int = NULL,
                robust_min_mac = 5, eig_tol = 1e-6, min_mac_int_spa = 10) {
  family <- match.arg(family)
  if (is.null(min_mac_int)) min_mac_int <- if (family == "binomial") 100 else 10
  shrink <- family == "binomial" && shrink
  pr <- if (shrink) .rv_prior(prior) else NULL
  if (is.null(robust)) robust <- family == "gaussian"
  if (missing(exposure) || !length(exposure)) stop("Give at least one exposure.", call. = FALSE)
  r <- .assoc_run(pfile, bfile, data, pheno, covar, exposure, int_covar, family, id_col, snp_vec, robust, center,
                  threads, mem_gb, block, spa_cut, carrier_max = if (shrink) carrier_max else 0,
                  carrier_lambda = if (shrink) 1 / pr$wide_sd^2 else 1)
  res <- r$res; out <- r$out; terms <- r$terms; ne <- r$n_exp
  lin <- family == "gaussian"
  pz <- function(b, s) if (lin) 2 * stats::pt(-abs(b / s), res$df_marginal) else 2 * stats::pnorm(-abs(b / s))
  out[, `:=`(BETA_MARGINAL = res$beta_marginal, SE_MARGINAL = res$se_marginal)]
  out[, P_MARGINAL := pz(BETA_MARGINAL, SE_MARGINAL)]
  if (!lin) out[!is.na(res$p_spa), P_MARGINAL := res$p_spa[!is.na(res$p_spa)]]
  if (robust) {
    out[, SE_MARGINAL_ROBUST := res$se_marginal_robust]
    out[, P_MARGINAL_ROBUST := .p_robust_1df((BETA_MARGINAL / SE_MARGINAL_ROBUST)^2, N, A1_FREQ, mac_adjust,
                                             P_MARGINAL, robust_min_mac)]
  }
  nm <- c("G", paste0("GxE_", terms))
  for (q in seq_along(nm)) {
    out[, (paste0("BETA_", nm[q])) := res$beta[, q]]
    out[, (paste0("SE_", nm[q])) := res$se[, q]]
    if (robust) out[, (paste0("SE_ROBUST_", nm[q])) := res$se_robust[, q]]
  }
  if (shrink) {
    pr <- .rv_learn(data.table::data.table(CHR = out$CHR, POS = out$POS, A1_FREQ = out$A1_FREQ,
                                           BETA = out$BETA_MARGINAL, SE = out$SE_MARGINAL), pr, names(prior), blocks)
    .rv_shrink_gxe(out, res, pr, nm, r$info$threads)
    r$info$prior <- pr
  }
  mac <- 2 * out$N * pmin(out$A1_FREQ, 1 - out$A1_FREQ)
  p_int <- stats::pchisq(res$chisq_int, ne, lower.tail = FALSE)
  p_joint <- stats::pchisq(res$chisq_joint, ne + 1L, lower.tail = FALSE)
  ki <- 1L + seq_len(ne)
  # single-term P of the exposure interactions: t (linear), saddlepoint below spa_cut else normal (logistic)
  z1 <- res$beta[, ki, drop = FALSE] / res$se[, ki, drop = FALSE]
  p1 <- if (lin) 2 * stats::pt(-abs(z1), res$df) else 2 * stats::pnorm(-abs(z1))
  if (!lin) { k <- !is.na(res$p_spa_int); p1[k] <- res$p_spa_int[k] }
  if (!lin && ne == 1L) p_int <- p1[, 1L]
  out[, `:=`(P_INT = p_int, P_JOINT = p_joint)]
  if (robust) {
    p1r <- vapply(seq_len(ne), function(q) .p_robust_1df((res$beta[, 1L + q] / res$se_robust[, 1L + q])^2, out$N,
                                                          out$A1_FREQ, mac_adjust, p1[, q], robust_min_mac),
                  numeric(nrow(out)))
    p1r <- matrix(p1r, ncol = ne)
    out[, `:=`(P_INT_ROBUST = if (ne == 1L) p1r[, 1L] else stats::pchisq(res$chisq_int_robust, ne, lower.tail = FALSE),
               P_JOINT_ROBUST = stats::pchisq(res$chisq_joint_robust, ne + 1L, lower.tail = FALSE))]
  }
  if (ne >= 2L) {
    # z values carried by the single-term P (so the VC part sees the same calibration as the single terms)
    cz <- function(pp, V) { zz <- stats::qnorm(pmax(pp, 1e-300) / 2, lower.tail = FALSE) * sign(z1)
      zz * sqrt(pmax(V[, (ki - 2L) * ne + (ki - 1L), drop = FALSE], 0)) }
    out[, P_INT_CAUCHY := int_cauchy_cpp(cz(p1, res$v_int), res$v_int, p1, r$info$threads)]
    if (robust)
      out[, P_INT_ROBUST_CAUCHY := int_cauchy_cpp(cz(p1r, res$v_int_robust), res$v_int_robust, p1r, r$info$threads)]
  }
  # Interaction tests not reported (NA): (a) below min_mac_int, every chi-square (Wald) interaction / joint test,
  # and for the linear model also the one-df and Cauchy interaction tests (logistic one-df and Cauchy tests use the
  # saddlepoint and are kept); (b) the information matrix of (G, interactions) near singular: smallest eigenvalue of
  # its unit-diagonal scaling below eig_tol, every interaction / joint test.
  ic <- intersect(c("P_INT", "P_JOINT", "P_INT_ROBUST", "P_JOINT_ROBUST", "P_INT_CAUCHY", "P_INT_ROBUST_CAUCHY"),
                  names(out))
  chi <- c("P_JOINT", "P_JOINT_ROBUST", if (ne >= 2L) c("P_INT", "P_INT_ROBUST"))
  lowc <- which(!(mac >= min_mac_int))
  if (length(lowc)) for (cl in if (lin) ic else intersect(chi, ic)) data.table::set(out, lowc, cl, NA_real_)
  if (!lin) {  # the saddlepoint-based logistic interaction tests have their own, lower floor
    lows <- which(!(mac >= min_mac_int_spa))
    if (length(lows)) for (cl in setdiff(ic, chi)) data.table::set(out, lows, cl, NA_real_)
  }
  sing <- which(!(res$eig_min >= eig_tol))
  if (length(sing)) for (cl in ic) data.table::set(out, sing, cl, NA_real_)
  r$info$int_low_mac <- length(lowc); r$info$int_singular <- length(sing); r$info$eig_min <- res$eig_min
  if (min_mac > 0) {
    pc <- grep("^P_", names(out), value = TRUE)
    lowp <- which(!(mac >= min_mac))
    if (length(lowp)) for (cl in pc) data.table::set(out, lowp, cl, NA_real_)
  }
  data.table::setattr(out, "info", c(r$info, list(exposure_terms = terms[seq_len(ne)],
                                                  int_covar_terms = terms[-seq_len(ne)])))
  out[]
}

# P of a one-df robust Wald statistic W (= T^2). Small-count reference (rare_variant_calibration.md, sec. 10):
# M = N (1 - (1 - p)^2), the HWE-expected carrier count from the minor allele frequency p; W is mapped to
# H = (M - 1) W / (M - W) and referred to F(1, M - 1). The mapping is only used to make a P larger, never smaller:
# P = max(P_F, P_normal). Its job is the inflated median of the normal P at low counts (W < 1); in the tail it
# would shrink P, which real data showed to be unsafe (dosages and non-HWE counts give W near or above M, where
# H -> infinity: a MAC-5 dosage variant went from 0.02 to 1e-19), so the tail keeps the normal P, which the
# simulations found conservative there. W >= M therefore gives the normal P, and a non-integer M is fine (F with
# real df). The floor "M >= 5" is checked on the minor allele count 2 N p (M = MAC - MAC^2 / 4N is a hair below 5
# at MAC 5); below it chi-square(1). With p_model given, variants with a minor allele count below min_mac_model take
# the model-based P instead (the robust variance has no information there). adjust = FALSE: chi-square(1) as before.
.p_robust_1df <- function(W, N, freq, adjust = TRUE, p_model = NULL, min_mac_model = 5) {
  p <- stats::pchisq(W, 1, lower.tail = FALSE)
  if (!adjust) return(p)
  maf <- rep_len(pmin(freq, 1 - freq), length(W)); N <- rep_len(N, length(W))
  if (!is.null(p_model)) {  # very low counts: the sandwich is unreliable (real data: lambda 2.2 at MAC 1-5)
    lo <- which(2 * N * maf < min_mac_model)
    p[lo] <- rep_len(p_model, length(W))[lo]
  }
  M <- N * (1 - (1 - maf)^2)
  k <- which(is.finite(W) & is.finite(M) & 2 * N * maf >= max(5 - 1e-8, min_mac_model) & W < M)
  if (length(k)) {
    Mk <- M[k]; Wk <- W[k]
    p[k] <- pmax(p[k], stats::pf((Mk - 1) * Wk / (Mk - Wk), 1, Mk - 1, lower.tail = FALSE))
  }
  p
}

# Shared part of gwas() and gxe(): samples, design, null model, variants, and the C++ stream.
.assoc_run <- function(pfile, bfile, data, pheno, covar, exposure, int_covar, family, id_col, snp_vec, robust,
                       center, threads, mem_gb, block, spa_cut = 0, carrier_max = 0,
                       carrier_lambda = 1) {
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
                          as.integer(block), if (family == "binomial") as.numeric(spa_cut) else 0,
                          as.numeric(carrier_max), as.numeric(carrier_lambda))
  out <- var[vrow, list(CHR, POS, ID, A1, A2)]
  out[, `:=`(N = as.integer(res$n_obs), A1_FREQ = res$a1_sum / (2 * res$n_obs))]
  info <- list(samples = n, covariate_rank = res$rank, null_converged = fit$converged, block = res$block,
               threads = res$threads, mem_gb = res$mem_gb, simd = res$simd, dosage = res$dosage,
               explicit_projection = res$explicit, spa = sum(!is.na(res$p_spa)), monomorphic = sum(res$status == 1L),
               singular = sum(res$status == 2L),
               skipped_multiallelic = skipped)
  list(res = res, out = out, info = info, terms = terms, n_exp = length(ie))
}
