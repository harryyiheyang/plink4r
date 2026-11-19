#' Polygenic score weights by SBayesRC on the eigen LD panel
#'
#' SBayesRC's model (Zheng et al. 2024) fitted as fastSBayesRC (Yihe Yang): short Gibbs chains for the effects and
#' stochastic-approximation EM for the annotation coefficients, on the int16 eigen panel. Per LD block the data are
#' w = Lambda^{-1} X' bhat with X = U Lambda^{1/2}; beta_j ~ sum_k pi_jk N(0, gam_k sigma^2) with gam_0 = 0, and
#' annotations enter only the prior weights through SBayesRC's nested structure
#' P(delta_j >= k | delta_j >= k - 1) = logistic(f_j' theta_k), f_j = (1, A_j), with the adaptive spike-and-slab lasso
#' on the slopes. Blocks are sampled in parallel.
#'
#' Steps (fast_chain.py):
#' 1. SBayesRC pseudo summaries: bhat_t = bhat + sqrt(1 / n_t - 1 / N) X z (z ~ N(0, I)), n_t = round(train_frac N),
#'    bhat_v = (N bhat - n_t bhat_t) / (N - n_t).
#' 2. Warm-start chain over `cutoffs` (eigen variance kept per block): `chain_iter` iterations per level on bhat_t,
#'    each level started from the previous one; pseudo-validation r = beta' bhat_v / sqrt(beta' R beta).
#' 3. SBayesRC rule: base = the largest cutoff; another level is chosen when its |r / r_base| > `rel_r_min` (r > 0), or
#'    when r_base < 0; choosing the smallest cutoff is an error (widen the grid).
#' 4. Main fit on bhat at the chosen cutoff, warm from that level: `main_iter` iterations, posterior mean over the
#'    iterations after `main_burn`.
#'
#' @param sumstat output of [summarystat_impute()] (or [summarystat_tidy()]): ID, A1, A2, Z, N (r2 optional; rows with
#'   r2 <= 0 are dropped). bhat = Z / sqrt(N + Z^2); N = median N. The original fits every panel SNP, so the input
#'   should cover the panel (the output of [summarystat_impute()]); otherwise only the panel rows with data are used.
#' @param ld_dir directory with block<N>.eigen.bin (plink4r int16 layout, or SBayesRC layout quantized to int16 on
#'   reading) and the panel SNP .parquet.
#' @param ref panel SNP table in place of the .parquet in `ld_dir`.
#' @param annot NULL, or a data.frame with ID and one numeric column per annotation (continuous columns as given);
#'   panel SNPs absent from it get 0.
#' @param control [sbayesrc_control()].
#' @param threads OpenMP threads (0: [plink4r_resources()]).
#' @param output optional output path (.parquet or text).
#' @return data.table: ID, chr, bp37, bp38, A1, A2, freq, BETA (per A1 allele), BETA_STD; attributes `chain` (one row
#'   per cutoff), `cutoff` (chosen), `theta` (levels x (intercept, annotations)), `sigma2`, `nnz` (mean number of
#'   nonzero effects after burn-in) and `trace` of the main fit.
#' @export
sbayesrc <- function(sumstat, ld_dir, ref = NULL, annot = NULL, control = sbayesrc_control(), threads = 4L,
                     output = NULL) {
  threads <- .threads(threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  ctl <- do.call(sbayesrc_control, control)
  A <- matrix(0, nrow(d$x), 0)
  if (!is.null(annot)) {
    an <- data.table::as.data.table(annot)
    if (!"ID" %in% names(an)) stop("annot needs an ID column", call. = FALSE)
    cols <- setdiff(names(an), "ID")
    A <- as.matrix(an[d$x[, .(ID)], on = "ID", ..cols])
    A[is.na(A)] <- 0
    storage.mode(A) <- "double"
    message(length(cols), " annotation columns; ", sum(!d$x$ID %in% an$ID), " SNPs without annotation (set to 0)")
  }
  t0 <- Sys.time()
  res <- sbayesrc_cpp(ld_dir, d$blocks, d$keep, d$x$Z, d$x$N, d$x$freq, A, control = ctl, replay = ctl$replay,
                      threads = threads)
  message(sprintf("sbayesrc: %d SNPs, cutoff %g, sigma2 %.3g, %.0f nonzero effects, %.1f s", nrow(d$x), res$cutoff,
                  res$sigma2, res$nnz, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out <- d$x[, .(ID, chr, bp37, bp38, A1, A2, freq)]
  out[, `:=`(BETA = res$BETA, BETA_STD = res$BETA_STD)]
  th <- res$theta
  dimnames(th) <- list(paste0("level", seq_len(nrow(th))), c("intercept", colnames(A)))
  attr(out, "chain") <- data.table::as.data.table(res$chain)
  attr(out, "cutoff") <- res$cutoff
  attr(out, "theta") <- th
  attr(out, "sigma2") <- res$sigma2
  attr(out, "nnz") <- res$nnz
  attr(out, "trace") <- data.table::as.data.table(res$trace)
  .prs_write(out, output)
  out[]
}

#' Control settings of [sbayesrc()]
#'
#' Defaults are those of fastSBayesRC (sbrc_mc.fit, fast_chain.py).
#' @param cutoffs eigen variance cutoffs of the chain, run in this order (the largest is the SBayesRC base).
#' @param chain_iter,chain_burn iterations per chain level and burn-in of its posterior mean.
#' @param main_iter,main_burn iterations and burn-in of the main fit.
#' @param rel_r_min SBayesRC rule threshold on |r / r_base|.
#' @param train_frac,pseudo_seed pseudo training proportion and its random seed.
#' @param seed seed of the Gibbs sampler (every fit restarts it, as the original).
#' @param gam,pi0 mixture variances (relative to sigma^2; 0 = point mass) and their starting weights.
#' @param ve residual variance of the eigen likelihood (precision n / ve).
#' @param nu prior degrees of freedom of sigma^2 (scaled inverse chi^2 centred at the current value).
#' @param hess_every refresh the Hessian of the annotation M-step every this many iterations.
#' @param slope_start annotation slopes enter after this many iterations (intercepts only before, unless warm-started).
#' @param warm iterations after `slope_start` with the ridge prior on the slopes before the SSL starts.
#' @param ridge_prec ridge precision of the slopes in the warm phase.
#' @param ssl_ratio,ssl_l1,ssl_beta,ssl_gamma,ssl_th_range adaptive SSL: lam0 = ssl_ratio lam1, start lam1, Beta mixture
#'   prior of the slab weight (rows a, b), Gamma(shape, rate) prior of lam1, clip of the slab weight.
#' @param newton_steps Newton steps of the intercept-only M-step.
#' @param pen_intercept ridge precision of the intercepts.
#' @param cd_maxit,cd_tol,irls_tol coordinate descent iterations / tolerance, Newton step tolerance.
#' @param eta_clip,w_floor,h_floor clip of the linear predictor, minimum soft-count weight, minimum Newton weight.
#' @param het_floor floor of 2p(1-p) in BETA = BETA_STD / sqrt(2p(1-p)).
#' @param lam_eps eigenvalues at or below are dropped from w.
#' @param vg_floor floor of the starting genetic variance.
#' @param replay NULL; a list(normal, seed, chisq) of recorded random draws replays them (validation only).
#' @return a list.
#' @export
sbayesrc_control <- function(cutoffs = c(0.9, 0.95, 0.99, 0.995, 0.999), chain_iter = 10L, chain_burn = NULL,
                             main_iter = 80L, main_burn = 30L, rel_r_min = 1.25, train_frac = 0.9, pseudo_seed = 22,
                             seed = 1, gam = c(0, 1e-3, 1e-2, 1e-1, 1), pi0 = NULL, ve = 1, nu = 4, hess_every = 5L,
                             slope_start = 5L, warm = 10L, ridge_prec = 1, ssl_ratio = 10, ssl_l1 = 0.5,
                             ssl_beta = rbind(c(1.5, 8.5), c(7.5, 42.5), c(30, 170)), ssl_gamma = c(2, 1),
                             ssl_th_range = c(1e-3, 0.999), newton_steps = 3L, pen_intercept = 1e-8, cd_maxit = 200L,
                             cd_tol = 1e-7, irls_tol = 1e-7, eta_clip = 30, w_floor = 1e-12, h_floor = 1e-12,
                             het_floor = 1e-6, lam_eps = 1e-8, vg_floor = 1e-6, replay = NULL) {
  K <- length(gam)
  if (is.null(pi0)) pi0 <- c(0.99, rep(0.01 / (K - 1), K - 1))
  if (is.null(chain_burn)) chain_burn <- chain_iter %/% 2L
  list(cutoffs = as.numeric(cutoffs), chain_iter = as.integer(chain_iter), chain_burn = as.integer(chain_burn),
       main_iter = as.integer(main_iter), main_burn = as.integer(main_burn), rel_r_min = as.numeric(rel_r_min),
       train_frac = as.numeric(train_frac), pseudo_seed = as.numeric(pseudo_seed), seed = as.numeric(seed),
       gam = as.numeric(gam), pi0 = as.numeric(pi0), ve = as.numeric(ve), nu = as.numeric(nu),
       hess_every = as.integer(hess_every), slope_start = as.integer(slope_start), warm = as.integer(warm),
       ridge_prec = as.numeric(ridge_prec), ssl_ratio = as.numeric(ssl_ratio), ssl_l1 = as.numeric(ssl_l1),
       ssl_beta = matrix(as.numeric(ssl_beta), ncol = 2), ssl_gamma = as.numeric(ssl_gamma),
       ssl_th_range = as.numeric(ssl_th_range), newton_steps = as.integer(newton_steps),
       pen_intercept = as.numeric(pen_intercept), cd_maxit = as.integer(cd_maxit), cd_tol = as.numeric(cd_tol),
       irls_tol = as.numeric(irls_tol), eta_clip = as.numeric(eta_clip), w_floor = as.numeric(w_floor),
       h_floor = as.numeric(h_floor), het_floor = as.numeric(het_floor), lam_eps = as.numeric(lam_eps),
       vg_floor = as.numeric(vg_floor), replay = replay)
}
