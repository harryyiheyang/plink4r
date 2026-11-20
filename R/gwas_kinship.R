#' REML estimate of the kinship variance of a (generalised) linear mixed model
#'
#' Fits the null model `formula + b`, `b ~ N(0, tau K)` with a sparse kinship matrix `K` and returns `tau`.
#' Linear (`family = "gaussian"`): `y = X beta + b + e`, `e ~ N(0, s2e I)`; one Brent search over
#' `h = tau / (tau + s2e)` of the profiled REML likelihood. Logistic (`"binomial"`): marginal
#' form, `logit(mu) = X beta` and `Var(y) = W^1/2 ((1 - h) I + h K) W^1/2`, `W = diag(mu (1 - mu))`, i.e. the
#' Pearson residuals `e = (y - mu) / sqrt(mu (1 - mu))` have covariance `H = (1 - h) I + h K`. Nested: for a given
#' `h`, `beta_hat(h)` solves `X' W^1/2 H^-1 e = 0` by Fisher scoring to convergence; `h` minimises
#' `f(h) = log|H| + sum log w + log|X' W^1/2 H^-1 W^1/2 X| + e' H^-1 e` at `beta_hat(h)` (the Gaussian restricted
#' quasi-likelihood of `y` with these two moments) by a Brent search started from a few steps of the alternating
#' scheme. `tau = h` and `s2e = 1 - h` are then on the Pearson-residual scale. Cox (`"cox"`; phenotype
#' `cbind(time, status)`): the same marginal form on the life table of the data. Follow-up is cut into
#' `para$n_interval` intervals at quantiles of the event times and a sample contributes one binary record to every
#' interval at whose start it is at risk (1 when its event falls in the interval), with
#' `P(event) = 1 - exp(-exp(alpha_j + X beta))` (complementary log-log link; the interval intercepts `alpha_j` are
#' the baseline). The records of one sample are uncorrelated, and the Pearson residuals of the records of two
#' relatives have correlation `h K_ij` whatever the two intervals. `h` and `(alpha, beta)` are found in turn until
#' the criterion settles: at the current working weights `h` minimises the working restricted likelihood by a
#' Brent search, then `(alpha, beta)` are updated by generalised least squares. With no related pair among the
#' analysed samples the fit is the Cox partial-likelihood fit (Breslow ties).
#'
#' Samples without a relative in `kinship` enter with the diagonal entry of `K` alone (a scalar variance, which
#' does not depend on `h` when that entry is 1, so they carry no information on `h`); the related samples are
#' held together as one sparse matrix and factorised by sparse Cholesky at every value of `h`, so the cost per
#' likelihood evaluation is one small sparse factorisation and solves with the related samples' rows only (the
#' cross-products over all samples are formed once per scoring step). With no related pair, or at `h = 0`, the
#' fit is ordinary least squares / logistic regression.
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
#' @param family `"gaussian"`, `"binomial"` or `"cox"`.
#' @param id_col Column of `data` with the sample IDs.
#' @param para named list overriding [reml_tau_para()].
#' @return list(tau, s2e, h2 = tau / (tau + s2e), lrt = likelihood-ratio statistic of tau = 0 (of the
#'   restricted quasi-likelihood for the logistic model), p_lrt (half chi-square 1), beta, n, n_related, pairs, converged, iterations,
#'   evals, h_bound = upper end of the search actually used, time; Cox: also baseline = data.frame(time,
#'   cumhaz), the cumulative baseline hazard at the end of every interval of the life table, whose log increments
#'   are the first entries of beta (with no related pair: at the event times, and beta's intercept is the log of
#'   the factor by which the Poisson fit rescales that baseline)).
#' @export
reml_tau <- function(formula, data, kinship, family = c("gaussian", "binomial", "cox"), id_col = "IID",
                     para = list()) {
  family <- match.arg(family)
  t0 <- proc.time()[["elapsed"]]
  P <- .reml_para(para)
  D <- .mm_design(formula, data, family, id_col, NULL)
  nl <- .reml_null(D, kinship, family, P, proj = FALSE)
  c(nl$est, list(time = proc.time()[["elapsed"]] - t0))
}

#' Defaults of [reml_tau()] and of [gwas()] with `kinship`
#'
#' h_max 0.999 upper end of the search for h (the search stops earlier, at `h_bound`, when `(1 - h) I + h K` is
#' not positive definite up to h_max, e.g. with a pair whose K_ij exceeds the diagonal);
#' tol 1e-6 Brent tolerance on h; max_iter 50 scoring steps per value of h (logistic); max_eval 100 evaluations
#' per Brent search;
#' mu_eps 1e-10 bound of the fitted probabilities; start_iter 2 steps of the alternating scheme that give the
#' start value of h (logistic); h_window 0.1 half-width of the search around it (the whole interval is searched
#' when the minimum falls on the edge of the window); f_tol 1e-8 convergence of the logistic fit is that of its criterion f: the
#' scoring stops when a step changes f by less than f_tol (1 + |f|), the search over h when its three best
#' points differ by less than f_tol (1 + |f(0)|), and h is exactly 0 unless it improves f(0) by more than that;
#' kin_scale 2: K_ij = kin_scale x
#' the pairwise kinship coefficient (0.25 for full siblings gives 0.5); diag 1 diagonal of K; rank_tol 1e-10 a
#' covariate column is aliased (dropped with a warning) when its residual sum of squares on the earlier columns
#' is below rank_tol x its own sum of squares; n_interval 10 intervals of the life table of the Cox model with
#' kinship (cut at quantiles of the event times; fewer when there are fewer distinct event times). For the Cox
#' model h is a correlation of records, of the order of 1 / n_interval, and max_iter counts the iterations.
#' @export
reml_tau_para <- function() {
  list(h_max = 0.999, tol = 1e-6, max_iter = 50L, max_eval = 100L, mu_eps = 1e-10,
       start_iter = 2L, h_window = 0.1, f_tol = 1e-8,
       kin_scale = 2, diag = 1, rank_tol = 1e-10, n_interval = 10L)
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
  nl <- .reml_null(D, kinship, family, P, proj = TRUE)
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
  data.table::setattr(out, "info", c(nl$est[setdiff(names(nl$est), c("beta", "baseline"))],
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
  list(ids = as.character(s$d[[id_col]]), y = s$y, status = s$status, X = X, srow = s$srow)
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
.reml_null <- function(D, kinship, family, P, proj) {
  kin <- .kin_prepare(kinship, D$ids, P$kin_scale)
  if (family == "gaussian") {
    f <- reml_tau_gaussian_cpp(D$y, D$X, kin$i, kin$j, kin$v, P$diag, P$h_max, P$tol, as.integer(P$max_eval),
                               P$rank_tol, proj)
    est <- list(tau = f$tau, s2e = f$s2e, h2 = f$h, converged = TRUE, iterations = 1L)
  } else if (family == "cox" && !length(kin$i)) {  # no related pair: the Cox partial-likelihood fit
    f <- cox_null_cpp(D$y, D$status, D$X, as.integer(P$max_iter), P$mu_eps, P$f_tol, P$rank_tol, proj)
    if (!f$converged) warning("The Cox null model did not converge.", call. = FALSE)
    f <- c(f, list(tau = 0, m2l_0 = f$m2l, n_related = 0L, pairs = 0, evals = 0L, h_bound = NA_real_))
    est <- list(tau = 0, s2e = 1, h2 = 0, converged = f$converged, iterations = f$iterations,
                baseline = data.frame(time = f$base_time, cumhaz = f$base_cumhaz))
  } else if (family == "cox") {
    if (P$diag != 1) stop("The Cox model with kinship takes diag = 1.", call. = FALSE)
    f <- reml_tau_cox_cpp(D$y, D$status, D$X, kin$i, kin$j, kin$v, as.integer(P$n_interval), P$h_max, P$tol,
                          as.integer(P$max_iter), as.integer(P$max_eval), P$mu_eps, P$f_tol, P$rank_tol, proj)
    if (!f$converged) warning("The Cox null model with kinship did not converge.", call. = FALSE)
    est <- list(tau = f$tau, s2e = 1 - f$tau, h2 = f$tau, converged = f$converged, iterations = f$iterations,
                baseline = data.frame(time = f$base_time, cumhaz = f$base_cumhaz))
  } else {
    f <- reml_tau_binomial_cpp(D$y, D$X, kin$i, kin$j, kin$v, P$diag, P$h_max, P$tol,
                               as.integer(P$max_iter), as.integer(P$max_eval), P$mu_eps, as.integer(P$start_iter),
                               P$h_window, P$f_tol, P$rank_tol, proj)
    if (!f$converged) warning("The logistic null model with kinship did not converge.", call. = FALSE)
    est <- list(tau = f$tau, s2e = 1 - f$tau, h2 = f$tau, converged = f$converged, iterations = f$iterations)
  }
  if (length(f$keep) < ncol(D$X))
    warning("Aliased covariate column(s) dropped: ", paste(colnames(D$X)[-f$keep], collapse = ", "), call. = FALSE)
  lrt <- max(0, f$m2l_0 - f$m2l)
  bn <- colnames(D$X)[f$keep]
  if (!is.null(f$n_interval))  # life table: the interval intercepts take the place of the intercept
    bn <- c(paste0("(Interval ", seq_len(f$n_interval), ")"), bn[-1L])
  est <- c(est, list(lrt = lrt, p_lrt = if (lrt > 0) 0.5 * stats::pchisq(lrt, 1, lower.tail = FALSE) else 1,
                     beta = stats::setNames(f$beta, bn), n = length(D$y),
                     n_related = f$n_related, pairs = f$pairs, evals = f$evals, h_bound = f$h_bound))
  list(est = est, null = f$null)
}
