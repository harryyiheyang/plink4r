#' REML estimate of the kinship variance of a (generalised) linear mixed model
#'
#' Fits the null model `y ~ covar + b`, `b ~ N(0, tau K)` with a sparse kinship matrix `K` and returns `tau`.
#' Linear (`family = "gaussian"`): `y = X beta + b + e`, `e ~ N(0, s2e I)`; one Brent search over
#' `h = tau / (tau + s2e)` of the profiled REML likelihood. Logistic (`"binomial"`): on the IRLS working scale,
#' in marginal form: `logit(mu) = X beta` and the Pearson residuals `(y - mu) / sqrt(mu (1 - mu))` have covariance
#' `(1 - h) I + h K`, so the working response has covariance `W^-1/2 ((1 - h) I + h K) W^-1/2`,
#' `W = diag(mu (1 - mu))`; `h` by one Brent search of the working REML likelihood, alternated with the
#' generalised-least-squares update of `beta` until both settle. `tau = h` and `s2e = 1 - h` are then on the
#' Pearson-residual scale.
#'
#' Samples without a relative in `kinship` enter with the diagonal entry of `K` alone (a scalar variance, which
#' does not depend on `h` when that entry is 1, so they carry no information on `h`); the related samples are
#' held together as one sparse matrix and factorised by sparse Cholesky at every value of `h`, so the cost per
#' likelihood evaluation is one small sparse factorisation plus `n p^2`. With no related pair, or at `h = 0`, the
#' fit is ordinary least squares / logistic regression.
#'
#' @param data data.frame with one row per sample: `id_col`, the phenotype and covariate columns.
#' @param pheno Name of the phenotype column. Logistic: coded 0/1 (or 1/2).
#' @param covar Names of covariate columns (numeric, or factors / characters, expanded by `model.matrix`); an
#'   intercept is always added. An externally estimated polygenic score is given here as a covariate.
#' @param kinship The related pairs: a data.frame with exactly three columns, in this order: the two sample IDs
#'   of a pair and their kinship coefficient as biobanks release it (0.25 for full siblings or parent and child,
#'   0.125 for second-degree relatives). Column names are not used. The matrix `K` is built here: `K_ij =
#'   para$kin_scale` (2) times the coefficient of the pair, `K_ii = para$diag` (1), 0 for every pair not listed.
#'   Pairs with a sample outside the analysis are dropped.
#' @param family `"gaussian"` or `"binomial"`.
#' @param id_col Column of `data` with the sample IDs.
#' @param para named list overriding [reml_tau_para()].
#' @return list(tau, s2e, h2 = tau / (tau + s2e), lrt = likelihood-ratio statistic of tau = 0 (on the working
#'   scale for the logistic model), p_lrt (half chi-square 1), beta, n, n_related, pairs, converged, iterations,
#'   evals, time).
#' @export
reml_tau <- function(data, pheno, covar = NULL, kinship, family = c("gaussian", "binomial"), id_col = "IID",
                     para = list()) {
  family <- match.arg(family)
  t0 <- proc.time()[["elapsed"]]
  P <- .reml_para(para)
  D <- .mm_design(data, pheno, covar, family, id_col, NULL)
  nl <- .reml_null(D, kinship, family, P)
  c(nl$est, list(time = proc.time()[["elapsed"]] - t0))
}

#' Defaults of [reml_tau()] and of [gwas()] with `kinship`
#'
#' h_max 0.999 upper end of the search for h;
#' tol 1e-6 Brent tolerance and change of h at convergence; eta_tol 1e-6 largest change of the linear
#' predictor at convergence (logistic); max_iter 50 outer iterations (logistic); max_eval 100 likelihood
#' evaluations per Brent search; mu_eps 1e-10 bound of the fitted probabilities; kin_scale 2: K_ij = kin_scale x
#' the pairwise kinship coefficient (0.25 for full siblings gives 0.5); diag 1 diagonal of K.
#' @export
reml_tau_para <- function() {
  list(h_max = 0.999, tol = 1e-6, eta_tol = 1e-6, max_iter = 50L, max_eval = 100L, mu_eps = 1e-10,
       kin_scale = 2, diag = 1)
}

# The kinship route of gwas(): the null model y ~ covar + b, b ~ N(0, tau K) is fitted once (.reml_null), then
# variants are streamed (gwas_kin_stream_cpp). With V the covariance of the (working) response under the null
# model and P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1, a variant g gets BETA = g' P y / g' P g and SE = (g' P g)^-1/2
# (generalised least squares; for the logistic model on the working response, one step from the null model).
# g' P g is exact (no GRAMMAR-GAMMA approximation). Same columns as the route without kinship.
.gwas_kinship <- function(pfile, bfile, data, pheno, covar, kinship, family, id_col, snp_vec, threads, mem_gb, block,
                          para) {
  t0 <- proc.time()[["elapsed"]]
  if (is.null(pfile) == is.null(bfile)) stop("Give exactly one of pfile and bfile.", call. = FALSE)
  P <- .reml_para(para)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)

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

  D <- .mm_design(data, pheno, covar, family, id_col, iid)
  nl <- .reml_null(D, kinship, family, P)
  t1 <- proc.time()[["elapsed"]]
  f <- nl$fit; kin <- nl$kin
  res <- gwas_kin_stream_cpp(path, is.null(pfile), length(iid), as.integer(var$allele_ct), as.integer(vrow - 1L),
                            as.integer(D$srow - 1L), f$Q, f$Vb, f$r, f$d, f$tau, kin$i, kin$j, kin$v, kin$diag,
                            threads, mem_gb, as.integer(block))
  out <- var[vrow, list(CHR, POS, ID, A1, A2)]
  out[, `:=`(N = as.integer(res$n_obs), A1_FREQ = res$a1_sum / (2 * res$n_obs), BETA = res$beta, SE = res$se)]
  stat <- if (family == "gaussian") "T" else "Z"
  out[, (stat) := BETA / SE]
  out[, P := 2 * stats::pnorm(-abs(get(stat)))]
  t2 <- proc.time()[["elapsed"]]
  data.table::setattr(out, "info", c(nl$est[setdiff(names(nl$est), "beta")],
                                     list(samples = nl$est$n, family = family, covariate_rank = ncol(D$X), block = res$block,
                                          threads = res$threads, mem_gb = res$mem_gb, simd = res$simd,
                                          dosage = res$dosage, explicit_projection = res$explicit,
                                          monomorphic = sum(res$status == 1L), singular = sum(res$status == 2L),
                                          skipped_multiallelic = skipped,
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

# Samples and design of a mixed null model: rows of `data` with complete phenotype and covariates (and, with
# `iid`, present in the genotype file, in file order; srow = their file rows), the phenotype (0/1 for a logistic
# model) and the numeric design X (intercept first; factors expanded by model.matrix; aliased columns dropped).
.mm_design <- function(data, pheno, covar, family, id_col, iid = NULL) {
  d <- as.data.frame(data)
  need <- c(id_col, pheno, covar)
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("data lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (length(pheno) != 1L) stop("Give one phenotype column.", call. = FALSE)
  ids <- as.character(d[[id_col]])
  if (anyDuplicated(ids[!is.na(ids)])) stop("IDs repeated in data.", call. = FALSE)
  ok <- !is.na(ids) & stats::complete.cases(d[, setdiff(need, id_col), drop = FALSE])
  srow <- NULL
  if (!is.null(iid)) {
    srow <- match(ids, iid)
    ok <- ok & !is.na(srow)
  }
  if (!any(ok)) stop("No sample of data has complete data", if (!is.null(iid)) " and is in the genotype file", ".",
                     call. = FALSE)
  d <- d[ok, , drop = FALSE]
  if (!is.null(iid)) {
    srow <- srow[ok]
    o <- order(srow); d <- d[o, , drop = FALSE]; srow <- srow[o]
  }
  y <- d[[pheno]]
  if (!is.numeric(y) && !is.logical(y)) stop("The phenotype must be numeric.", call. = FALSE)
  y <- as.numeric(y)
  if (family == "binomial") {
    u <- sort(unique(y))
    if (all(u %in% c(1, 2))) y <- y - 1 else if (!all(u %in% c(0, 1)))
      stop("A logistic phenotype must be coded 0/1 or 1/2.", call. = FALSE)
    if (length(unique(y)) < 2L) stop("The phenotype has one class only.", call. = FALSE)
  }
  fml <- stats::reformulate(if (length(covar)) sprintf("`%s`", covar) else "1")
  X <- stats::model.matrix(fml, data = d)
  attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL
  q <- qr(X)
  if (q$rank < ncol(X)) {
    drop <- colnames(X)[q$pivot[-seq_len(q$rank)]]
    warning("Aliased covariate column(s) dropped: ", paste(drop, collapse = ", "), call. = FALSE)
    X <- X[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  }
  list(ids = as.character(d[[id_col]]), y = y, X = X, srow = srow)
}

# K over the analysis samples `ids` from the related pairs (a three-column data.frame: ID, ID, kinship
# coefficient): off-diagonal pairs (i, j 0-based, each pair once) with K_ij = scale x coefficient, and the diagonal
# (diag for everyone). Pairs with a sample outside `ids`, self pairs and zero coefficients are dropped; a repeated
# pair keeps its first row.
.kin_prepare <- function(kinship, ids, scale, diag) {
  if (!is.data.frame(kinship) || ncol(kinship) != 3L)
    stop("kinship must be a data.frame with three columns: ID, ID, kinship coefficient.", call. = FALSE)
  v <- as.numeric(kinship[[3L]])
  if (anyNA(v)) stop("kinship has missing coefficients.", call. = FALSE)
  i <- match(as.character(kinship[[1L]]), ids); j <- match(as.character(kinship[[2L]]), ids)
  ok <- !is.na(i) & !is.na(j) & i != j & v != 0
  i <- i[ok]; j <- j[ok]; v <- v[ok]
  lo <- pmin(i, j); hi <- pmax(i, j)
  keep <- !duplicated(data.frame(lo, hi))
  list(i = as.integer(lo[keep] - 1L), j = as.integer(hi[keep] - 1L), v = as.numeric(scale) * v[keep],
       diag = rep(as.numeric(diag), length(ids)))
}

# Null model of reml_tau() and of gwas() with kinship: the estimates (est) and the quantities of the association stream (fit:
# tau, d, r, Q, Vb) with the kinship of its V = diag(d) + tau K (kin; scaled by W^-1/2 for the logistic model).
.reml_null <- function(D, kinship, family, P) {
  kin <- .kin_prepare(kinship, D$ids, P$kin_scale, P$diag)
  n <- length(D$y)
  if (family == "gaussian") {
    f <- reml_tau_gaussian_cpp(D$y, D$X, kin$i, kin$j, kin$v, kin$diag, P$h_max, P$tol, as.integer(P$max_eval))
    lrt <- max(0, f$m2l_0 - f$m2l)
    est <- list(tau = f$tau, s2e = f$s2e, h2 = f$h, lrt = lrt, converged = TRUE, iterations = 1L)
  } else {
    f <- reml_tau_binomial_cpp(D$y, D$X, kin$i, kin$j, kin$v, kin$diag, P$h_max, P$tol, P$eta_tol,
                               as.integer(P$max_iter), as.integer(P$max_eval), P$mu_eps)
    if (!f$converged) warning("The logistic null model with kinship did not converge.", call. = FALSE)
    lrt <- max(0, f$m2l_0 - f$m2l)
    est <- list(tau = f$tau, s2e = 1 - f$tau, h2 = f$tau, lrt = lrt, converged = f$converged,
                iterations = f$iterations)
    kin$v <- f$ks; kin$diag <- f$ksdiag  # the stream's V = diag(d) + tau W^-1/2 K W^-1/2
  }
  est <- c(est, list(p_lrt = if (lrt > 0) 0.5 * stats::pchisq(lrt, 1, lower.tail = FALSE) else 1,
                     beta = stats::setNames(f$beta, colnames(D$X)), n = n, n_related = f$n_related,
                     pairs = f$pairs, evals = f$evals))
  list(est = est, fit = f, kin = kin)
}
