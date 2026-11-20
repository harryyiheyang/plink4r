#' REML estimate of the kinship variance of a (generalised) linear mixed model
#'
#' Fits the null model `formula + b`, `b ~ N(0, tau K)` with a sparse kinship matrix `K` and returns `tau`.
#' Linear (`family = "gaussian"`): `y = X beta + b + e`, `e ~ N(0, s2e I)`, so `V = s2e I + tau K`; one Brent
#' search over `h = tau / (tau + s2e)` of the profiled REML likelihood. Logistic (`"binomial"`):
#' `logit(mu) = X beta + b`, fitted on the IRLS working scale by penalised quasi-likelihood (as GMMAT and
#' fastGWA-GLMM): the working response has covariance `V = W^-1 + tau K`, `W = diag(mu (1 - mu))`; `tau` by one
#' Brent search of the working REML likelihood, alternated with the update of `beta` and `b` until both settle.
#' `tau` is then a variance on the logit scale.
#'
#' Samples without a relative in `kinship` enter with the diagonal entry of `K` alone (a scalar variance); the
#' related samples are held together as one sparse matrix and factorised by sparse Cholesky at every value of the
#' variance parameter. Linear: the cross-products over all samples are formed once, and a likelihood evaluation
#' touches the related samples' rows only. Logistic: the weights differ between samples, so every evaluation also
#' forms the weighted cross-products over all samples (`n p^2`). With `tau = 0` the fit is ordinary least squares /
#' logistic regression.
#'
#' @param formula The fixed part of the model, `phenotype ~ covariates`, as in [gwas()] (factors and characters
#'   expanded by `model.matrix`). Logistic: the phenotype is coded 0/1 (or 1/2). An externally estimated polygenic
#'   score is given here as a covariate.
#' @param data data.frame with one row per sample: `id_col` and the columns named in `formula`.
#' @param kinship The related pairs: a data.frame with exactly three columns, in this order: the two sample IDs
#'   of a pair and their kinship coefficient as biobanks release it (0.25 for full siblings or parent and child,
#'   0.125 for second-degree relatives). Column names are not used. The matrix `K` is built here: `K_ij =
#'   para$kin_scale` (2) times the coefficient of the pair, `K_ii = para$diag` (1), 0 for every pair not listed.
#'   Pairs with a sample outside the analysis are dropped.
#' @param family `"gaussian"` or `"binomial"`.
#' @param id_col Column of `data` with the sample IDs.
#' @param para named list overriding [reml_tau_para()].
#' @param threads OpenMP threads (logistic model); 0 takes them from [plink4r_resources()].
#' @return list(tau, s2e and h2 = tau / (tau + s2e) (linear model; `NA` for the logistic model), lrt =
#'   likelihood-ratio statistic of tau = 0 (on the working scale for the logistic model), p_lrt (half chi-square
#'   1), beta, n, n_related, pairs, converged, iterations, evals, time).
#' @export
reml_tau <- function(formula, data, kinship, family = c("gaussian", "binomial"), id_col = "IID", para = list(),
                     threads = 0L) {
  family <- match.arg(family)
  t0 <- proc.time()[["elapsed"]]
  P <- .reml_para(para)
  D <- .mm_design(formula, data, family, id_col, NULL)
  nl <- .reml_null(D, kinship, family, P, proj = FALSE, .threads(threads))
  c(nl$est, list(time = proc.time()[["elapsed"]] - t0))
}

#' Defaults of [reml_tau()] and of [gwas()] with `kinship`
#'
#' h_max 0.999 upper end of the search for h (linear); tau_max 10 upper end of the search for tau (logistic);
#' tol 1e-6 Brent tolerance and (relative) change of tau at convergence; eta_tol 1e-6 largest change of the linear
#' predictor at convergence (logistic); max_iter 50 outer iterations (logistic); max_eval 100 likelihood
#' evaluations per Brent search; mu_eps 1e-10 bound of the fitted probabilities; kin_scale 2: K_ij = kin_scale x
#' the pairwise kinship coefficient (0.25 for full siblings gives 0.5); diag 1 diagonal of K; rank_tol 1e-10 a
#' covariate column is aliased (dropped with a warning) when its residual sum of squares on the earlier columns
#' is below rank_tol x its own sum of squares.
#' @export
reml_tau_para <- function() {
  list(h_max = 0.999, tau_max = 10, tol = 1e-6, eta_tol = 1e-6, max_iter = 50L, max_eval = 100L, mu_eps = 1e-10,
       kin_scale = 2, diag = 1, rank_tol = 1e-10)
}

# The kinship route of gwas(): the null model formula + b, b ~ N(0, tau K) is fitted once (.reml_null), then
# variants are streamed (gwas_kin_stream_cpp) against that same null-model object. With V the covariance of the
# (working) response under the null model and P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1, a variant g gets
# BETA = g' P y / g' P g and SE = (g' P g)^-1/2 (generalised least squares; for the logistic model on the working
# response, one step from the null model). g' P g is exact (no GRAMMAR-GAMMA approximation). Same columns as the
# route without kinship.
.gwas_kinship <- function(formula, data, pfile, bfile, kinship, family, id_col, snp_vec, threads, mem_gb, block,
                          para) {
  t0 <- proc.time()[["elapsed"]]
  P <- .reml_para(para)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  gv <- .geno_variants(pfile, bfile, snp_vec)
  D <- .mm_design(formula, data, family, id_col, gv$iid)
  nl <- .reml_null(D, kinship, family, P, proj = TRUE, threads)
  t1 <- proc.time()[["elapsed"]]
  res <- gwas_kin_stream_cpp(gv$path, is.null(pfile), length(gv$iid), as.integer(gv$var$allele_ct),
                             as.integer(gv$vrow - 1L), as.integer(D$srow - 1L), nl$null, threads, mem_gb,
                             as.integer(block))
  out <- gv$var[gv$vrow, list(CHR, POS, ID, A1, A2)]
  out[, `:=`(N = as.integer(res$n_obs), A1_FREQ = res$a1_sum / (2 * res$n_obs), BETA = res$beta, SE = res$se)]
  stat <- if (family == "gaussian") "T" else "Z"
  out[, (stat) := BETA / SE]
  out[, P := 2 * stats::pnorm(-abs(get(stat)))]
  t2 <- proc.time()[["elapsed"]]
  data.table::setattr(out, "info", c(nl$est[setdiff(names(nl$est), "beta")],
                                     list(samples = nl$est$n, family = family, covariate_rank = length(nl$est$beta),
                                          block = res$block, threads = res$threads, mem_gb = res$mem_gb,
                                          simd = res$simd, dosage = res$dosage, explicit_projection = res$explicit,
                                          monomorphic = sum(res$status == 1L), singular = sum(res$status == 2L),
                                          skipped_multiallelic = gv$skipped,
                                          time = c(null = t1 - t0, stream = t2 - t1, total = t2 - t0))))
  out[]
}

.reml_para <- function(para) {
  P <- reml_tau_para()
  bad <- setdiff(names(para), names(P))
  if (length(bad)) stop("unknown para: ", paste(bad, collapse = ", "), call. = FALSE)
  P[names(para)] <- para
  P
}

# Samples and design of a mixed null model: the samples of .pheno_samples, their phenotype and the numeric design
# X (intercept first; factors expanded by model.matrix). Aliased columns are found and left out by the fit.
.mm_design <- function(formula, data, family, id_col, iid = NULL) {
  s <- .pheno_samples(as.data.frame(data), formula, family, id_col, iid)
  X <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data = s$d)
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL
  list(ids = as.character(s$d[[id_col]]), y = s$y, X = X, srow = s$srow)
}

# K over the analysis samples `ids` from the related pairs (a three-column data.frame: ID, ID, kinship
# coefficient): off-diagonal pairs (i, j 0-based, each pair once) with K_ij = scale x coefficient. Pairs with a
# sample outside `ids`, self pairs and zero coefficients are dropped; a repeated pair keeps its first row.
.kin_prepare <- function(kinship, ids, scale) {
  if (!is.data.frame(kinship) || ncol(kinship) != 3L)
    stop("kinship must be a data.frame with three columns: ID, ID, kinship coefficient.", call. = FALSE)
  v <- as.numeric(kinship[[3L]])
  if (anyNA(v)) stop("kinship has missing coefficients.", call. = FALSE)
  i <- match(as.character(kinship[[1L]]), ids); j <- match(as.character(kinship[[2L]]), ids)
  ok <- !is.na(i) & !is.na(j) & i != j & v != 0
  i <- i[ok]; j <- j[ok]; v <- v[ok]
  lo <- pmin(i, j); hi <- pmax(i, j)
  keep <- !duplicated(data.frame(lo, hi))
  list(i = as.integer(lo[keep] - 1L), j = as.integer(hi[keep] - 1L), v = as.numeric(scale) * v[keep])
}

# Null model of reml_tau() and of gwas() with kinship: the estimates (est) and, with proj, the null-model object of
# the association stream (null: external pointer, made and read in compiled code only).
.reml_null <- function(D, kinship, family, P, proj, threads) {
  kin <- .kin_prepare(kinship, D$ids, P$kin_scale)
  if (family == "gaussian") {
    f <- reml_tau_gaussian_cpp(D$y, D$X, kin$i, kin$j, kin$v, P$diag, P$h_max, P$tol, as.integer(P$max_eval),
                               P$rank_tol, proj)
    est <- list(tau = f$tau, s2e = f$s2e, h2 = f$h, converged = TRUE, iterations = 1L)
  } else {
    f <- reml_tau_binomial_cpp(D$y, D$X, kin$i, kin$j, kin$v, P$diag, P$tau_max, P$tol, P$eta_tol,
                               as.integer(P$max_iter), as.integer(P$max_eval), P$mu_eps, P$rank_tol, proj,
                               as.integer(threads))
    if (!f$converged) warning("The logistic null model with kinship did not converge.", call. = FALSE)
    est <- list(tau = f$tau, s2e = NA_real_, h2 = NA_real_, converged = f$converged, iterations = f$iterations)
  }
  if (length(f$keep) < ncol(D$X))
    warning("Aliased covariate column(s) dropped: ", paste(colnames(D$X)[-f$keep], collapse = ", "), call. = FALSE)
  lrt <- max(0, f$m2l_0 - f$m2l)
  est <- c(est, list(lrt = lrt, p_lrt = if (lrt > 0) 0.5 * stats::pchisq(lrt, 1, lower.tail = FALSE) else 1,
                     beta = stats::setNames(f$beta, colnames(D$X)[f$keep]), n = length(D$y),
                     n_related = f$n_related, pairs = f$pairs, evals = f$evals))
  list(est = est, null = f$null)
}
