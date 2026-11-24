#' Polygenic score weights by VIPRS with an annotation model of the per-SNP variance
#'
#' [viprs()] with the slab variance of each SNP given by its annotations (LDSC-type variance model, log scale):
#' beta_j | gamma_j = 1 ~ N(0, v_j), v_j = exp(alpha_0 + sum_k A_jk alpha_k). With pi fixed at each grid point, the
#' coefficients alpha are estimated by EM together with sigma: the M-step maximizes
#' sum_j [-pi/2 log v_j - zeta_j / (2 v_j)] (zeta_j = E[gamma_j b_j^2], weight pi per SNP as VIPRS 0.1.5), which is
#' concave in alpha (Newton steps with step halving). With no annotation column this is [viprs()].
#' Search, PUMAS pseudo-validation and refit as [viprs()]; a move of pi along the grid shifts alpha_0 by
#' log(prev pi / pi) (keeps the implied h2).
#'
#' @inheritParams viprs
#' @param annot data.frame with ID and one numeric column per annotation (binary or continuous, used as given); SNPs
#'   not in `annot` get 0.
#' @param control [viprs_annot_control()].
#' @return as [viprs()], plus column VAR (fitted slab variance v_j) and attributes `coef` (alpha: intercept and one
#'   per annotation), `enrichment` (per annotation: annotation-weighted mean of v_j over the mean of all SNPs) and
#'   `grid_coef` (alpha at each grid point).
#' @export
viprs_annot <- function(sumstat, ld_dir, annot, ref = NULL, thresh = 0, fp16 = TRUE, control = viprs_annot_control(),
                        threads = 4L, output = NULL) {
  threads <- .threads(threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  ctl <- do.call(viprs_annot_control, control)
  an <- data.table::as.data.table(annot)
  if (!"ID" %in% names(an)) stop("annot needs an ID column", call. = FALSE)
  cols <- setdiff(names(an), "ID")
  A <- if (length(cols)) as.matrix(an[d$x[, .(ID)], on = "ID", ..cols]) else matrix(0, nrow(d$x), 0)
  A[is.na(A)] <- 0
  storage.mode(A) <- "double"
  message(length(cols), " annotation columns; ", sum(!d$x$ID %in% an$ID), " SNPs without annotation (set to 0)")
  t0 <- Sys.time()
  res <- viprs_annot_cpp(ld_dir, d$blocks, d$keep, d$x$Z, d$x$N, d$x$freq, A, thresh = thresh, fp16 = fp16,
                         control = ctl, threads = threads)
  hyp <- res$hyp
  message(sprintf("viprs_annot: %d SNPs, sum k %d; pi %.4g, h2 %.4f, sigma %.4f, %d iterations (%s), %.1f s",
                  as.integer(hyp$n_snps), as.integer(hyp$sum_k), hyp$pi, hyp$h2, hyp$sigma, as.integer(hyp$iters),
                  hyp$message, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out <- d$x[, .(ID, chr, bp37, bp38, A1, A2, freq)]
  out[, `:=`(BETA = res$BETA, BETA_STD = res$BETA_STD, PIP = res$PIP, VAR = res$VAR)]
  attr(out, "hyp") <- hyp
  attr(out, "coef") <- stats::setNames(res$coef, c("intercept", cols))
  attr(out, "enrichment") <- stats::setNames(res$enrichment, cols)
  attr(out, "grid") <- data.table::as.data.table(res$extra$grid)
  gc <- res$extra$grid_coef
  colnames(gc) <- c("intercept", cols)
  attr(out, "grid_coef") <- gc
  attr(out, "best_pseudo_r2") <- res$extra$best_pseudo_r2
  attr(out, "pumas_xi") <- res$extra$pumas_xi
  .prs_write(out, output)
  out[]
}

#' Control settings of [viprs_annot()]
#'
#' [viprs_control()] settings plus those of the M-step of the variance model.
#' @param ... passed to [viprs_control()].
#' @param mstep_max_iter,mstep_tol Newton steps per M-step at most, and stop when the largest coefficient change is
#'   below `mstep_tol`.
#' @param mstep_max_halving step halvings per Newton step at most (a step that does not raise the M-step objective
#'   after these is not taken).
#' @return a list.
#' @export
viprs_annot_control <- function(..., mstep_max_iter = 20L, mstep_tol = 1e-8, mstep_max_halving = 40L) {
  c(viprs_control(...), list(mstep_max_iter = as.integer(mstep_max_iter), mstep_tol = as.numeric(mstep_tol),
                             mstep_max_halving = as.integer(mstep_max_halving)))
}
