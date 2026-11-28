#' Annotation-partitioned polygenic score weights by VIPRS (PartitionPRS)
#'
#' [viprs()] with an additive annotation model of the per-SNP effect variance, and the partition of each effect into
#' annotation-specific components (tissue-CRE-specific PRS weights; Research Design D.1.2). Model:
#' beta_j = sum_k beta_j^(k), beta_j^(k) independent with variance A_jk sigma2_k (A_jk >= 0), so
#' h2_j = Var(beta_j) = sigma2_0 + sum_k A_jk sigma2_k, with an intercept sigma2_0 on every SNP (A_j0 = 1, the Base
#' of S-LDSC). Prior on beta_j:
#' VIPRS spike and slab with one pi, chosen on the grid; the slab N(0, h2_j / pi) (default `mix_prop = 1`,
#' `mix_var = 1`; a fixed mixture of classes is allowed), so E[beta_j^2] = h2_j. M-step (Research Design D.1.2):
#' sigma2_k <- S_k / sum_j A_jk, S_k = sum_j [w_jk^2 E[beta_j^2] + w_jk (1 - w_jk) h2_j] (E[beta_j^2] the posterior
#' second moment, spike included), whose fixed point sum_j A_jk (E[beta_j^2] - h2_j) / h2_j^2 = 0 is the stationary point
#' of Fm = -1/2 sum_j (log h2_j + E[beta_j^2] / h2_j) over the SNPs not suppressed. MASH prior on the annotation
#' variances: sigma2_k = a_k sigma2_0, a_k in `mash_grid` (0, 1, 10, 100) with fixed weights `mash_w` (0.5, 0.3, 0.1,
#' 0.1) common to every annotation; VI averages over a_k with q(a_k = m) proportional to w_m exp(Fm(a_k = m)) (the other
#' columns at their mean) and uses E_q[a_k]; sigma2_0 by monotone MM steps on Fm; the Gibbs fit keeps sigma2_0 and draws
#' a_k from q.
#' Output: the MMSE allocation of the posterior mean, BETA_k = w_jk BETA (k = intercept and each annotation); the
#' columns sum to BETA for every SNP and the annotation-specific scores sum to the total PRS for every individual.
#' Hyperparameters by VI. Stage 1, the eigen cutoff (as [sbayesrc()]) on one PUMAS split (the same pseudo training /
#' validation data for every fit), each fit scored by r = eta' valid / sqrt(eta' R eta): a chain over `cutoffs`
#' (default 0.995, 0.999) from the lowest up, each level warm-started from the one below, at the fixed pi `cut_pi`
#' (default 0.01); the largest cutoff unless a smaller one has r / r_largest > `rel_r_min` (SBayesRC rule). Stage 2, pi
#' at the chosen cutoff, warm-started from its stage-1 fit: the pi grid of [viprs()] ascending on all data, each point
#' warm-started from the last; the converged point of the largest ELBO + E_q[log w(a)] - E_q[log q(a)] is refitted to
#' full convergence. The output sigma2_k of the Gibbs fit is its mean over the draws.
#' `sigma2_init` sets the sigma2_k of the first fit (stage 1, lowest cutoff; later fits are warm-started), e.g. the stratified LDSC estimates of [ldsc_strat()].
#' Final fit at the chosen pi, started from the VI refit: a short Gibbs sampler (`gibbs_burn` + `gibbs_draws` sweeps)
#' with Rao-Blackwell means; the atoms a_k drawn from q on stochastic-approximation averages of E[beta_j^2], sigma by
#' stochastic approximation. BETA is the mean of the conditional posterior means m_j over the draws and
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
#'   and the annotation names; negative values set to 0, then the start rule of `start_scale` and `start_zero_share`
#'   (a zero column, the intercept included, gets a positive start).
#' @param bin_file LD-bin list (default `ld_dir`/ldbin<bin_r2>.bin, from [ld_bins()]; must exist).
#' @param skip_annot annotation columns that make a SNP annotated (any nonzero value) for the LD-bin suppression:
#'   NULL (default) every column, character(0) no suppression.
#' @param control [viprs_annot_control()].
#' @return data.table as [viprs()] (BETA, BETA_STD, PIP), plus H2_SNP (h2_j), SUPPRESSED and one column
#'   BETA_<annotation> per annotation and BETA_intercept (BETA_STD scale: BETA_STD_<annotation>, ...); per SNP the
#'   BETA_<...> columns sum to BETA. Attributes `sigma2` (sigma2_k), `h2_annot` (model h2_k = (sum_j A_jk) sigma2_k
#'   over the SNPs not suppressed; k = intercept: M sigma2_0) and `h2_share`; `h2_shapley` (heritability of each
#'   column as its Shapley value in beta' R beta: the mean over the draws of beta^(k)' R beta, cross terms split evenly
#'   between the two columns, so the values add up to the total; VI: b_k' R b), `mash` (sigma2_0, q: q(a_k) of the VI
#'   refit, annotations x atoms; mean_atom: E_q[a_k]; pzero_gibbs: the mean over the draws of q(a_k = 0); draws),
#'   `grid` (one row per fit: stage, cutoff, pi, its pseudo-validation r (NaN for ELBO fits on all data) and ELBO),
#'   `grid_sigma2` (sigma2_k at each grid point), `cutoff`, `best_pseudo_r` (stage 1, the chosen cutoff) and
#'   `best_pi_crit` (stage 2: ELBO plus the MASH part at the chosen pi), `n_bins`, `n_suppressed`, and those of [viprs()].
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
  attr(out, "n_bins") <- res$extra$n_bins
  attr(out, "n_suppressed") <- n_sup
  attr(out, "grid") <- data.table::as.data.table(res$extra$grid)
  gs <- res$extra$grid_sigma2
  colnames(gs) <- nm
  attr(out, "grid_sigma2") <- gs
  attr(out, "best_pseudo_r") <- res$extra$best_pseudo_r
  attr(out, "best_pi_crit") <- res$extra$best_pi_crit
  mash <- res$extra$mash
  dimnames(mash$q) <- list(nm, paste0("a", ctl$mash_grid))
  if (length(mash$mean_atom) == length(nm)) names(mash$mean_atom) <- nm
  if (length(mash$pzero_gibbs) == length(nm)) names(mash$pzero_gibbs) <- nm
  attr(out, "mash") <- mash
  attr(out, "cutoff") <- res$extra$cutoff
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
#' @param start_scale,start_zero_share with `sigma2_init` (e.g. ldsc_strat()), h the LDSC h2: negative
#'   values set to 0; every sigma2_k times `start_scale`; a column at 0 (the intercept included) set so that its h2
#'   (sum_j A_jk sigma2_k) is `start_zero_share` h.
#' @param gibbs_burn,gibbs_draws sweeps of the final Gibbs fit before and after burn-in; `gibbs_draws = 0` keeps the
#'   VI refit as the final fit.
#' @param fix_sigma2 TRUE keeps sigma2_k at the start (no M-step; a diagnostic with the true values as `sigma2_init`,
#'   `start_scale = 1`).
#' @param mash_grid,mash_w MASH prior on the annotation variances: sigma2_k = a_k sigma2_0, a_k in `mash_grid`
#'   (0 first, then increasing) with fixed weights `mash_w`, common to every annotation; averaged in VI (q(a_k)),
#'   drawn in the Gibbs fit.
#' @return a list.
#' @export
viprs_annot_control <- function(..., bin_r2 = 0.98, mix_prop = 1, mix_var = 1,
                                cutoffs = c(0.995, 0.999), rel_r_min = 1.25, cut_pi = 0.01,
                                start_scale = 1.3, start_zero_share = 0.05,
                                gibbs_burn = 100L, gibbs_draws = 300L, fix_sigma2 = FALSE,
                                mash_grid = c(0, 1, 10, 100), mash_w = c(0.5, 0.3, 0.1, 0.1)) {
  mash_grid <- as.numeric(mash_grid); mash_w <- as.numeric(mash_w)
  if (length(mash_grid) < 2 || length(mash_grid) != length(mash_w) || mash_grid[1] != 0 || is.unsorted(mash_grid, strictly = TRUE) ||
      any(mash_w <= 0)) stop("mash_grid: 0 first, then strictly increasing positive atoms; mash_w: one positive weight per atom")
  c(viprs_control(...), list(bin_r2 = as.numeric(bin_r2), mix_prop = as.numeric(mix_prop),
                             mix_var = as.numeric(mix_var), cutoffs = as.numeric(cutoffs),
                             rel_r_min = as.numeric(rel_r_min), cut_pi = as.numeric(cut_pi),
                             start_scale = as.numeric(start_scale),
                             start_zero_share = as.numeric(start_zero_share),
                             gibbs_burn = as.integer(gibbs_burn),
                             gibbs_draws = as.integer(gibbs_draws),
                             fix_sigma2 = isTRUE(fix_sigma2),
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
