#' Annotation-partitioned polygenic score weights by VIPRS (PartitionPRS)
#'
#' [viprs()] with an additive annotation model of the per-SNP effect variance, and the partition of each effect into
#' annotation-specific components (tissue-CRE-specific PRS weights; Research Design D.1.2). Model:
#' beta_j = sum_k beta_j^(k), beta_j^(k) independent with variance A_jk sigma2_k (A_jk >= 0), so
#' h2_j = Var(beta_j) = sigma2_0 + sum_k A_jk sigma2_k, with an intercept sigma2_0 on every SNP (A_j0 = 1, the Base
#' of S-LDSC). Prior on beta_j:
#' VIPRS spike and slab with pi fixed on the grid; the slab is a mixture of small / medium / large effects with fixed
#' shares (`mix_prop`, default 6 : 3 : 1) and variances `mix_var` (relative, default 0.01 : 0.1 : 1) times
#' h2_j / pi, scaled so that Var(beta_j) = h2_j. M-step (D.1.2): with
#' w_jk = A_jk sigma2_k / h2_j and E[(beta_j^(k))^2] = w_jk^2 zeta_j + w_jk (1 - w_jk) h2_j (zeta_j the posterior
#' second moment of beta_j on the scale of h2_j),
#' sigma2_k <- (S_k + 2 b0) / (sum_j A_jk + 2 (a0 + 1)), S_k = sum_j E[(beta_j^(k))^2] (Inverse-Gamma(a0, b0) prior).
#' No floor: an annotation at sigma2_k = 0 (e.g. a negative LDSC start cut at 0) stays 0 (b0 = 0) and its SNPs fall
#' back on the intercept.
#' Output: the MMSE allocation of the posterior mean, BETA_k = w_jk BETA (k = intercept and each annotation); the
#' columns sum to BETA for every SNP and the annotation-specific scores sum to the total PRS for every individual.
#' Search, PUMAS pseudo-validation and refit as [viprs()] (pi grid ascending, each point warm-started from the last);
#' `sigma2_init` sets the sigma2_k of the first grid point, e.g. the stratified LDSC estimates of [ldsc_strat()].
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
#'   and the annotation names. Annotation values are cut at 0; the intercept must be > 0 (else the LDSC start).
#' @param bin_file LD-bin list (default `ld_dir`/ldbin<bin_r2>.bin, from [ld_bins()]; must exist).
#' @param skip_annot annotation columns that make a SNP annotated (any nonzero value) for the LD-bin suppression:
#'   NULL (default) every column, character(0) no suppression.
#' @param control [viprs_annot_control()].
#' @return data.table as [viprs()] (BETA, BETA_STD, PIP), plus H2_SNP (h2_j), SUPPRESSED and one column
#'   BETA_<annotation> per annotation and BETA_intercept (BETA_STD scale: BETA_STD_<annotation>, ...); per SNP the
#'   BETA_<...> columns sum to BETA. Attributes `sigma2` (sigma2_k), `h2_annot` (h2_k = (sum_j A_jk) sigma2_k over
#'   the SNPs not suppressed; k = intercept: M sigma2_0) and `h2_share`, `h2_post` (heritability of each partitioned PRS,
#'   E[beta^(k)' R beta^(k)] = b_k' R b_k + sum_j R_jj Var(beta_j^(k)) with the posterior variance; the PRS add up,
#'   these do not: the columns have cross terms),
#'   `grid_sigma2` (sigma2_k at each grid point), `n_bins`,
#'   `n_suppressed`, and those of [viprs()].
#' @export
viprs_annot <- function(sumstat, ld_dir, annot, ref = NULL, bin_file = NULL, skip_annot = NULL, sigma2_init = NULL,
                        thresh = 0, fp16 = TRUE, control = viprs_annot_control(), threads = 4L, output = NULL) {
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
                         thresh = thresh, fp16 = fp16, control = ctl, threads = threads)
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
  attr(out, "h2_post") <- stats::setNames(res$h2_post, nm)
  attr(out, "n_bins") <- res$extra$n_bins
  attr(out, "n_suppressed") <- n_sup
  attr(out, "grid") <- data.table::as.data.table(res$extra$grid)
  gs <- res$extra$grid_sigma2
  colnames(gs) <- nm
  attr(out, "grid_sigma2") <- gs
  attr(out, "best_pseudo_r2") <- res$extra$best_pseudo_r2
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
#' @param a0,b0 Inverse-Gamma(a0, b0) prior of each sigma2_k: sigma2_k <- (S_k + 2 b0) / (sum_j A_jk + 2 (a0 + 1)).
#'   a0 = -1, b0 = 0 is the flat prior (with no annotation column this is [viprs()]).
#' @return a list.
#' @export
viprs_annot_control <- function(..., bin_r2 = 0.98, mix_prop = c(0.6, 0.3, 0.1), mix_var = c(0.01, 0.1, 1),
                                a0 = 0, b0 = 0) {
  c(viprs_control(...), list(bin_r2 = as.numeric(bin_r2), mix_prop = as.numeric(mix_prop),
                             mix_var = as.numeric(mix_var),
                             a0 = as.numeric(a0), b0 = as.numeric(b0)))
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
