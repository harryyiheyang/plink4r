#' Annotation-partitioned polygenic score weights by VIPRS (PartitionPRS)
#'
#' [viprs()] with an additive annotation model of the per-SNP effect variance, and the partition of each effect into
#' annotation-specific components (tissue-CRE-specific PRS weights; Research Design D.1.2). Model:
#' beta_j = sum_k beta_j^(k), beta_j^(k) independent with variance A_jk sigma2_k (A_jk >= 0), so
#' h2_j = Var(beta_j) = sum_k A_jk sigma2_k; a residual annotation A_R = 1 for every SNP is added. Prior on beta_j:
#' VIPRS spike and slab with pi fixed on the grid; the slab is a mixture of small / medium / large effects with fixed
#' shares (`mix_prop`, default 6 : 3 : 1) and variances `mix_var` (relative, default 0.01 : 0.1 : 1) times
#' h2_j / pi, scaled so that Var(beta_j) = h2_j. M-step (D.1.2): with
#' w_jk = A_jk sigma2_k / h2_j and E[(beta_j^(k))^2] = w_jk^2 zeta_j + w_jk (1 - w_jk) h2_j (zeta_j the posterior
#' second moment of beta_j on the scale of h2_j),
#' sigma2_k <- (S_k + 2 b0) / (sum_j A_jk + 2 (a0 + 1)), S_k = sum_j E[(beta_j^(k))^2] (Inverse-Gamma(a0, b0) prior).
#' The fitted posterior mean is allocated as beta_j^(k) = w_jk beta_j (MMSE allocation), so the annotation columns
#' sum to BETA for every SNP and the annotation-specific scores sum to the total PRS for every individual:
#' `geno_score(pfile, out, score_cols = grep("^BETA_", names(out), value = TRUE))`.
#' Search, PUMAS pseudo-validation and refit as [viprs()] (pi grid ascending, each point warm-started from the last);
#' `sigma2_init` sets the sigma2_k of the first grid point, e.g. the stratified LDSC estimates of [ldsc_strat()].
#'
#' LD-bin suppression (as [sbayesrc()]): inside an LD bin of the list `bin_file` (bins of r2 > `bin_r2`, default 0.98, see
#' [sbayesrc()] for the format) with at least one annotated member (any nonzero binary 0/1 annotation column), the
#' members without annotation are fixed at beta = 0 (never updated, so they never change the residual) and are left
#' out of the M-step. Without binary annotation columns nothing is suppressed.
#'
#' @inheritParams viprs
#' @param annot annotations keyed by Block and Index (or CHR and BP37 / BP38), one numeric column per annotation
#'   (binary 0/1 or continuous score >= 0);
#'   SNPs not in `annot` get 0.
#' @param resid_all residual annotation A_R = 1 for every SNP (TRUE) or only for SNPs without annotation (FALSE).
#' @param sigma2_init NULL (all sigma2_k equal, sum_j h2_j = the naive LDSC h2), the output of [ldsc_strat()] with the
#'   same annotations (tau_c; with `resid_all = FALSE` the annotation values are tau_0 + tau_k, exact for
#'   non-overlapping 0/1 annotations), or a numeric vector named `residual` and the annotation names. Values are cut
#'   at 0 and then raised to at least `init_floor` times the mean h2_j (a sigma2_k at 0 would stay 0).
#' @param bin_file LD-bin list (default `ld_dir`/ldbin<bin_r2>.bin; must exist).
#' @param control [viprs_annot_control()].
#' @return data.table as [viprs()] (BETA, BETA_STD, PIP), plus H2_SNP (h2_j), SUPPRESSED and one column
#'   BETA_<annotation> per annotation and BETA_residual (BETA_STD scale: BETA_STD_<annotation>, ...); per SNP the
#'   BETA_<...> columns sum to BETA. Attributes `sigma2` (sigma2_k), `h2_annot` (h2_k = (sum_j A_jk) sigma2_k over the
#'   SNPs not suppressed) and `h2_share` (h2_k / sum_k h2_k), `grid_sigma2` (sigma2_k at each grid point), `n_bins`,
#'   `n_suppressed`, and those of [viprs()].
#' @export
viprs_annot <- function(sumstat, ld_dir, annot, ref = NULL, bin_file = NULL, resid_all = TRUE, sigma2_init = NULL,
                        thresh = 0, fp16 = TRUE, control = viprs_annot_control(), threads = 4L, output = NULL) {
  threads <- .threads(threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  ctl <- do.call(viprs_annot_control, control)
  A <- .prs_annot(annot, d$x)
  cols <- colnames(A)
  if ("residual" %in% cols) stop("the annotation name 'residual' is reserved", call. = FALSE)
  if (any(A < 0)) stop("annotations must be >= 0", call. = FALSE)
  s2 <- .sigma2_init(sigma2_init, cols, resid_all)
  if (is.null(bin_file)) bin_file <- file.path(ld_dir, paste0("ldbin", ctl$bin_r2, ".bin"))
  if (!file.exists(bin_file)) stop("LD-bin list not found: ", bin_file, call. = FALSE)
  t0 <- Sys.time()
  res <- viprs_annot_cpp(.ld_block_file(ld_dir, d$blocks), d$blocks, d$keep, d$x$Z, d$x$N, d$x$A1Freq, A, bin_file, d$start, d$n_panel,
                         resid_all = resid_all, sigma2_init = s2, thresh = thresh, fp16 = fp16, control = ctl, threads = threads)
  n_sup <- sum(res$extra$suppressed)
  hyp <- res$hyp
  message(sprintf("viprs_annot: %d SNPs (%d suppressed in %.0f LD bins), sum k %d; pi %.4g, h2 %.4f, sigma %.4f, %d iterations (%s), %.1f s",
                  as.integer(hyp$n_snps), n_sup, res$extra$n_bins, as.integer(hyp$sum_k), hyp$pi, hyp$h2, hyp$sigma,
                  as.integer(hyp$iters), hyp$message, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  nm <- c("residual", cols)
  out <- d$x[, .(CHR, BP37, BP38, A1, A2, A1Freq, Block, Index)]
  out[, `:=`(BETA = res$BETA, BETA_STD = res$BETA_STD, PIP = res$PIP, H2_SNP = res$H2_SNP,
             SUPPRESSED = res$extra$suppressed)]
  sc <- ifelse(res$BETA_STD != 0, res$BETA / res$BETA_STD, 0)  # 1 / sqrt(2 p (1 - p)) per SNP
  ord <- c(cols, "residual")
  for (k in ord) {
    i <- match(k, nm)
    out[, (paste0("BETA_", k)) := res$PART_STD[, i] * sc]
  }
  for (k in ord) out[, (paste0("BETA_STD_", k)) := res$PART_STD[, match(k, nm)]]
  attr(out, "hyp") <- hyp
  attr(out, "sigma2") <- stats::setNames(res$sigma2, nm)
  attr(out, "h2_annot") <- stats::setNames(res$h2_annot, nm)
  attr(out, "h2_share") <- stats::setNames(res$h2_annot / sum(res$h2_annot), nm)
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
#' @param init_floor smallest start value of a sigma2_k given by `sigma2_init`, as a fraction of the mean h2_j.
#' @param a0,b0 Inverse-Gamma(a0, b0) prior of each sigma2_k: sigma2_k <- (S_k + 2 b0) / (sum_j A_jk + 2 (a0 + 1)).
#'   a0 = -1, b0 = 0 is the flat prior (with no annotation column this is [viprs()]).
#' @return a list.
#' @export
viprs_annot_control <- function(..., bin_r2 = 0.98, mix_prop = c(0.6, 0.3, 0.1), mix_var = c(0.01, 0.1, 1),
                                init_floor = 0.01, a0 = 0, b0 = 0) {
  c(viprs_control(...), list(bin_r2 = as.numeric(bin_r2), mix_prop = as.numeric(mix_prop),
                             mix_var = as.numeric(mix_var), init_floor = as.numeric(init_floor), a0 = as.numeric(a0), b0 = as.numeric(b0)))
}

# sigma2_k of the first grid point (residual first, then cols), or numeric(0) for the LDSC start
.sigma2_init <- function(x, cols, resid_all) {
  if (is.null(x)) return(numeric(0))
  if (is.list(x) && !is.null(x$Parameter_summary)) {
    p <- x$Parameter_summary
    tau <- stats::setNames(p$Estimate, p$Parameter)
    miss <- setdiff(c("Residual", cols), names(tau))
    if (length(miss)) stop("ldsc_strat() result lacks ", paste(miss, collapse = ", "), call. = FALSE)
    t0 <- tau[["Residual"]]
    return(c(t0, unname(tau[cols]) + if (resid_all) 0 else t0))
  }
  if (!is.numeric(x) || is.null(names(x))) stop("sigma2_init: NULL, a ldsc_strat() result or a named numeric vector",
                                                call. = FALSE)
  miss <- setdiff(c("residual", cols), names(x))
  if (length(miss)) stop("sigma2_init lacks ", paste(miss, collapse = ", "), call. = FALSE)
  unname(x[c("residual", cols)])
}
