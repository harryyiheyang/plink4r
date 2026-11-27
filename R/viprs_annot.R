#' Annotation-partitioned polygenic score weights by VIPRS (PartitionPRS)
#'
#' [viprs()] with an additive annotation model of the per-SNP effect variance, and the partition of each effect into
#' annotation-specific components (tissue-CRE-specific PRS weights; Research Design D.1.2). Model:
#' beta_j = sum_k beta_j^(k), beta_j^(k) independent with variance A_jk sigma2_k (A_jk >= 0), so
#' h2_j = Var(beta_j) = sigma2_0 + sum_k A_jk sigma2_k, with an intercept sigma2_0 on every SNP (A_j0 = 1, the Base
#' of S-LDSC). Prior on beta_j:
#' VIPRS spike and slab with pi fixed on the grid; the slab is a mixture of small / medium / large effects with fixed
#' shares (`mix_prop`, default 6 : 3 : 1) and variances `mix_var` (relative, default 0.01 : 0.1 : 1) times
#' h2_j / pi, scaled so that Var(beta_j) = h2_j. M-step: sigma2_k maximize the ELBO, whose sigma2-dependent part is
#' -1/2 sum_j (gamma_j log h2_j + pi zeta_j / h2_j) (gamma_j the slab probability, zeta_j = gamma_j E[beta_j^2 / kappa_c]
#' the slab second moment on the scale of h2_j; SNPs not suppressed): stationarity
#' sum_j A_jk (pi zeta_j / h2_j^2 - gamma_j / h2_j) = 0, by the monotone MM step
#' sigma2_k <- sigma2_k sqrt(sum_j A_jk pi zeta_j / h2_j^2 / sum_j A_jk gamma_j / h2_j) (10 per iteration); the intercept is
#' kept above 1e-20. sigma2_k = 0 is a fixed point, so the start (see `start_*` of [viprs_annot_control()]) gives every
#' column, the intercept included, a positive value.
#' Output: the MMSE allocation of the posterior mean, BETA_k = w_jk BETA (k = intercept and each annotation); the
#' columns sum to BETA for every SNP and the annotation-specific scores sum to the total PRS for every individual.
#' Hyperparameters by VI on one PUMAS split (the same pseudo training / validation data for every fit), each fit
#' scored by r = eta' valid / sqrt(eta' R eta). Stage 1, the eigen cutoff (as [sbayesrc()]): a chain over `cutoffs`
#' (default 0.995, 0.999) from the lowest up, each level warm-started from the one below, at the fixed pi `cut_pi`
#' (default 0.01); the largest cutoff unless a smaller one has r / r_largest > `rel_r_min` (SBayesRC rule). Stage 2, pi
#' at the chosen cutoff, warm-started from its stage-1 fit: the pi grid of [viprs()] ascending on
#' all data, each point warm-started from the last (sigma2_k carried over); the converged point of the largest ELBO is
#' refitted to full convergence (`pi_by = "pseudo_r"`: the grid on the PUMAS training half, scored like stage 1, the
#' best point refitted on all data). pi and sigma2_k are not chosen jointly by one criterion (weakly identified).
#' `sigma2_init` sets the sigma2_k of the first fit of each stage, e.g. the stratified LDSC estimates of [ldsc_strat()].
#' Final fit at the chosen pi, started from the VI refit: a short Gibbs sampler (`gibbs_burn` + `gibbs_draws` sweeps)
#' with Rao-Blackwell means; sigma2_k by the same M-step on stochastic-approximation averages of gamma_j and zeta_j,
#' sigma by stochastic approximation. BETA is the mean of the conditional posterior means m_j over the draws and
#' BETA_k the mean of w_jk m_j (w_jk of each sweep), so the columns still sum to BETA.
#'
#' LD-bin suppression (as [sbayesrc()]): inside an LD bin of the list `bin_file` (bins of r2 > `bin_r2`, default 0.98,
#' made by [ld_bins()]) with at least one annotated member (nonzero in any of the columns `skip_annot`), the members
#' without annotation are fixed at beta = 0 (never updated) and are left out of the
#' M-step.
#'
#' @inheritParams viprs
#' @param annot annotations keyed by Block and Index (or CHR and BP37 / BP38), one numeric column per annotation
#'   (binary 0/1 or continuous score >= 0);
#'   SNPs not in `annot` get 0.
#' @param sigma2_init NULL (all sigma2_k equal, sum_j h2_j = the naive LDSC h2), the output of [ldsc_strat()] with the
#'   same annotations (tau_c; its Base row, named Residual, is the intercept), or a numeric vector named `intercept`
#'   and the annotation names; negative values set to 0, then the start rule of `start_scale`, `start_zero_share`,
#'   `start_total_min` (a zero column, the intercept included, gets a positive start).
#' @param bin_file LD-bin list (default `ld_dir`/ldbin<bin_r2>.bin, from [ld_bins()]; must exist).
#' @param skip_annot annotation columns that make a SNP annotated (any nonzero value) for the LD-bin suppression:
#'   NULL (default) every column, character(0) no suppression.
#' @param control [viprs_annot_control()].
#' @return data.table as [viprs()] (BETA, BETA_STD, PIP), plus H2_SNP (h2_j), SUPPRESSED and one column
#'   BETA_<annotation> per annotation and BETA_intercept (BETA_STD scale: BETA_STD_<annotation>, ...); per SNP the
#'   BETA_<...> columns sum to BETA. Attributes `sigma2` (sigma2_k), `h2_annot` (model h2_k = (sum_j A_jk) sigma2_k
#'   over the SNPs not suppressed; k = intercept: M sigma2_0) and `h2_share`; `h2_shapley` (heritability of each
#'   column as its Shapley value in beta' R beta: the mean over the draws of beta^(k)' R beta, cross terms split evenly
#'   between the two columns, so the values add up to the total) with `h2_shapley_sd` (sd over the draws);
#'   `var_shapley` (the same for the output PRS, b_k' R b: they add up to b' R b, the PRS variance on the panel),
#'   `grid` (one row per fit: stage, cutoff, pi, its pseudo-validation r (NaN for ELBO fits on all data) and ELBO),
#'   `grid_sigma2` (sigma2_k at each grid point), `cutoff`, `best_pseudo_r` (stage 1, the chosen cutoff) and
#'   `best_pi_crit` (stage 2: ELBO or pseudo r of the chosen pi), `n_bins`, `n_suppressed`, and those of [viprs()].
#' @export
viprs_annot <- function(sumstat, ld_dir, annot, ref = NULL, bin_file = NULL, skip_annot = NULL, sigma2_init = NULL,
                        fp16 = TRUE, control = viprs_annot_control(), threads = 4L, output = NULL) {
  threads <- .threads(threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  ctl <- do.call(viprs_annot_control, control)
  A <- .prs_annot(annot, d$x)
  cols <- colnames(A)
  if ("intercept" %in% cols) stop("the annotation name 'intercept' is reserved", call. = FALSE)
  if (any(A < 0)) stop("annotations must be >= 0", call. = FALSE)
  big <- cols[apply(A, 2, max) > 1]
  if (length(big)) warning("annotation values above 1 in ", paste(big, collapse = ", "),
                           ": 0-100 scores are read as probabilities only after dividing by 100", call. = FALSE)
  s2 <- .sigma2_init(sigma2_init, cols)
  if (is.null(skip_annot)) skip_annot <- cols
  miss <- setdiff(skip_annot, cols)
  if (length(miss)) stop("skip_annot names no annotation column: ", paste(miss, collapse = ", "), call. = FALSE)
  if (is.null(bin_file)) bin_file <- file.path(ld_dir, paste0("ldbin", ctl$bin_r2, ".bin"))
  if (!file.exists(bin_file)) stop("LD-bin list not found: ", bin_file, call. = FALSE)
  t0 <- Sys.time()
  res <- viprs_annot_cpp(.ld_block_file(ld_dir, d$blocks), d$blocks, d$keep, d$x$Z, d$x$N, d$x$A1Freq, A, bin_file, d$start, d$n_panel,
                         sigma2_init = s2, skip_cols = match(skip_annot, cols) - 1L,
                         fp16 = fp16, control = ctl, threads = threads)
  n_sup <- sum(res$extra$suppressed)
  hyp <- res$hyp
  message(sprintf("viprs_annot: %d SNPs (%d suppressed in %.0f LD bins), sum k %d; pi %.4g, h2 %.4f, sigma %.4f, %d iterations (%s), %.1f s",
                  as.integer(hyp$n_snps), n_sup, res$extra$n_bins, as.integer(hyp$sum_k), hyp$pi, hyp$h2, hyp$sigma,
                  as.integer(hyp$iters), hyp$message, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  nm <- c("intercept", cols)
  out <- d$x[, .(CHR, BP37, BP38, A1, A2, A1Freq, Block, Index)]
  out[, `:=`(BETA = res$BETA, BETA_STD = res$BETA_STD, PIP = res$PIP, H2_SNP = res$H2_SNP,
             SUPPRESSED = res$extra$suppressed)]
  sc <- ifelse(res$BETA_STD != 0, res$BETA / res$BETA_STD, 0)  # 1 / sqrt(2 p (1 - p)) per SNP
  ord <- c(cols, "intercept")
  for (k in ord) {
    i <- match(k, nm)
    out[, (paste0("BETA_", k)) := res$PART_STD[, i] * sc]
  }
  for (k in ord) out[, (paste0("BETA_STD_", k)) := res$PART_STD[, match(k, nm)]]
  attr(out, "hyp") <- hyp
  attr(out, "sigma2") <- stats::setNames(res$sigma2, nm)
  attr(out, "h2_annot") <- stats::setNames(res$h2_annot, nm)
  attr(out, "h2_share") <- stats::setNames(res$h2_annot / sum(res$h2_annot), nm)
  attr(out, "h2_shapley") <- stats::setNames(res$h2_shapley, nm)
  attr(out, "h2_shapley_sd") <- stats::setNames(res$h2_shapley_sd, nm)
  attr(out, "var_shapley") <- stats::setNames(res$var_shapley, nm)
  attr(out, "n_bins") <- res$extra$n_bins
  attr(out, "n_suppressed") <- n_sup
  attr(out, "grid") <- data.table::as.data.table(res$extra$grid)
  gs <- res$extra$grid_sigma2
  colnames(gs) <- nm
  attr(out, "grid_sigma2") <- gs
  attr(out, "best_pseudo_r") <- res$extra$best_pseudo_r
  attr(out, "best_pi_crit") <- res$extra$best_pi_crit
  attr(out, "ssl") <- res$extra$ssl
  attr(out, "mash") <- res$extra$mash
  attr(out, "cutoff") <- res$extra$cutoff
  attr(out, "pumas_xi") <- res$extra$pumas_xi
  .prs_write(out, output)
  out[]
}

#' Control settings of [viprs_annot()]
#'
#' [viprs_control()] settings plus those of the annotation variance model.
#' @param ... passed to [viprs_control()].
#' @param bin_r2 r2 threshold of the default LD-bin list `ld_dir`/ldbin<bin_r2>.bin.
#' @param mix_prop prior shares of the slab classes (small, medium, large), fixed; normalized to sum 1.
#' @param mix_var relative variances of the slab classes; one class (`mix_prop = 1, mix_var = 1`) is the Gaussian
#'   slab of [viprs()].
#' @param cutoffs eigen cutoffs (eigenvalue variance kept per block) compared at pi = `cut_pi`; one value fixes the
#'   cutoff; a value at or above the panel's stored cutoff keeps every stored pair.
#' @param rel_r_min a smaller cutoff is taken when its r / r_largest exceeds this.
#' @param cut_pi the fixed pi of the cutoff choice (stage 1).
#' @param pi_by criterion of the pi grid at the chosen cutoff (stage 2): "elbo" (all data) or "pseudo_r".
#' @param start_scale,start_zero_share,start_total_min with `sigma2_init` (e.g. ldsc_strat()), h the LDSC h2: negative
#'   values set to 0; every sigma2_k times `start_scale`; a column at 0 (the intercept included) set so that its h2
#'   (sum_j A_jk sigma2_k) is `start_zero_share` h; then all scaled up so that the start h2 is at least
#'   `start_total_min` h (sigma2_k = 0 is a fixed point of the M-step, so no column may start at 0).
#' @param gibbs_burn,gibbs_draws sweeps of the final Gibbs fit before and after burn-in; `gibbs_draws = 0` keeps the
#'   VI refit as the final fit.
#' @param fix_sigma2 TRUE keeps sigma2_k at the start (no M-step; a diagnostic with the true values as `sigma2_init`,
#'   `start_scale = 1`, `start_total_min = 0`).
#' @param select annotation selection on sigma2_k (not the intercept): "ssl" (default), the spike-and-slab lasso on
#'   s_k = sigma2_k / b (b = LDSC h2 / SNPs not suppressed): prior theta Exp(lambda1) + (1 - theta) Exp(lambda0),
#'   lambda0 = `ssl_ratio` lambda1; M-step s_k <- s_k sqrt(N_k / (D_k + 2 lambda*_k / b)), lambda*_k =
#'   p*_k lambda1 + (1 - p*_k) lambda0, p*_k the posterior slab probability of s_k; exact zeros by the KKT condition
#'   (a column at 0 stays there while b (N_k - D_k) <= 2 lambda*_k(0)); lambda1 and theta by EB after each M-step,
#'   lambda1 <- (sum p* + g1 - 1) / (sum p* s + g2) with Gamma(`ssl_gamma`) and theta <- (sum p* + a - 1) / (K + a + b - 2)
#'   with Beta(`ssl_beta`). The result has attribute `ssl` (lambda1, theta, scale b).
#'   "mash": atomic prior sigma2_k in `mash_grid` x sigma2_0 (default 0, 1, 10, 100 times the intercept) with fixed
#'   weights `mash_w` (default 0.5, 0.3, 0.1, 0.1); sigma2_0 by the ELBO (MM step of the column
#'   a_j = 1 + sum_k A_jk grid[m_k]); each
#'   column's atom m_k by the largest F(sigma2_k = grid[m] b) + log w_m in VI (MAP; atom 0 is an exact zero), drawn from
#'   p_km proportional to exp(F + log w_m) in the Gibbs fit. Attribute `mash`: scale (sigma2_0), the final atoms (0-based), pzero
#'   (the Gibbs mean of p(m_k = 0)), draws.
#' @return a list.
#' @export
viprs_annot_control <- function(..., bin_r2 = 0.98, mix_prop = c(0.6, 0.3, 0.1), mix_var = c(0.01, 0.1, 1),
                                cutoffs = c(0.995, 0.999), rel_r_min = 1.25, cut_pi = 0.01, pi_by = c("elbo", "pseudo_r"),
                                start_scale = 1.3, start_zero_share = 0.05, start_total_min = 1.3,
                                gibbs_burn = 100L, gibbs_draws = 300L, fix_sigma2 = FALSE,
                                select = c("ssl", "mash"), ssl_ratio = 10, ssl_gamma = c(2, 1), ssl_beta = c(1, 1),
                                mash_grid = c(0, 1, 10, 100), mash_w = c(0.5, 0.3, 0.1, 0.1)) {
  c(viprs_control(...), list(bin_r2 = as.numeric(bin_r2), mix_prop = as.numeric(mix_prop),
                             mix_var = as.numeric(mix_var), cutoffs = as.numeric(cutoffs),
                             rel_r_min = as.numeric(rel_r_min), cut_pi = as.numeric(cut_pi), pi_by = match.arg(pi_by),
                             start_scale = as.numeric(start_scale),
                             start_zero_share = as.numeric(start_zero_share), start_total_min = as.numeric(start_total_min),
                             gibbs_burn = as.integer(gibbs_burn),
                             gibbs_draws = as.integer(gibbs_draws),
                             fix_sigma2 = isTRUE(fix_sigma2), select = match.arg(select), ssl_ratio = as.numeric(ssl_ratio),
                             ssl_gamma = as.numeric(ssl_gamma), ssl_beta = as.numeric(ssl_beta),
                             mash_grid = as.numeric(mash_grid), mash_w = as.numeric(mash_w)))
}

# sigma2_k of the first grid point (intercept first, then cols), or numeric(0) for the LDSC start
.sigma2_init <- function(x, cols) {
  if (is.null(x)) return(numeric(0))
  if (is.list(x) && !is.null(x$Parameter_summary)) {
    p <- x$Parameter_summary
    tau <- stats::setNames(p$Estimate, p$Parameter)
    miss <- setdiff(c("Residual", cols), names(tau))
    if (length(miss)) stop("ldsc_strat() result lacks ", paste(miss, collapse = ", "), call. = FALSE)
    return(c(tau[["Residual"]], unname(tau[cols])))
  }
  if (!is.numeric(x) || is.null(names(x))) stop("sigma2_init: NULL, a ldsc_strat() result or a named numeric vector",
                                                call. = FALSE)
  miss <- setdiff(c("intercept", cols), names(x))
  if (length(miss)) stop("sigma2_init lacks ", paste(miss, collapse = ", "), call. = FALSE)
  unname(x[c("intercept", cols)])
}
