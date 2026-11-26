#' Annotation mixture partition of polygenic score weights
#'
#' Each SNP's effect comes from one class: null, or (group g, size c) for a group g the SNP belongs to (annotations
#' with A_jg > 0, plus a residual group) and a size c in (large, medium, small) with variance gamma_c sigma2 shared
#' across groups. P(z_j = (g, c)) = A_jg w_gc / (1 + s_j), P(null) = 1 / (1 + s_j), s_j = sum_gc A_jg w_gc. Fit by
#' mean-field coordinate ascent on the eigen panel (b_hat ~ N(R beta, R / n), R never formed), sigma2 by EM and w by
#' an MM update with a spike-and-slab lasso type penalty on each annotation's total weight (residual group not
#' penalised). The weight of annotation g is the posterior mean restricted to g's classes; per SNP the BETA_<group>
#' columns sum to BETA, so [geno_score()] on them gives annotation-specific scores that sum to the total score.
#'
#' @inheritParams viprs
#' @param annot annotations keyed by Block and Index (or CHR and BP37 / BP38), one numeric column per annotation (>= 0);
#'   SNPs not in `annot` get 0.
#' @param resid_all residual group on every SNP (TRUE) or only on SNPs without annotation.
#' @param control [mixpart_control()].
#' @return data.table: CHR, BP37, BP38, A1, A2, A1Freq, Block, Index, BETA, BETA_STD, BETA_<annotation>, BETA_residual and the
#'   BETA_STD_<...> columns; attributes `h2_share` (posterior), `h2_share_prior`, `w` (groups x size classes),
#'   `sigma2`, `iters`, `converged`.
#' @export
mixpart <- function(sumstat, ld_dir, annot, ref = NULL, resid_all = TRUE, thresh = 0, fp16 = TRUE,
                    control = mixpart_control(), threads = 4L, output = NULL) {
  threads <- .threads(threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  ctl <- do.call(mixpart_control, control)
  A <- .prs_annot(annot, d$x)
  cols <- colnames(A)
  if ("residual" %in% cols) stop("the annotation name 'residual' is reserved", call. = FALSE)
  t0 <- Sys.time()
  res <- mixpart_cpp(.ld_block_file(ld_dir, d$blocks), d$blocks, d$keep, d$x$Z, d$x$N, A, resid_all, thresh, fp16, ctl, threads)
  message(sprintf("mixpart: %d SNPs, %d annotations, %d iterations (%s), %.1f s", nrow(d$x), length(cols),
                  as.integer(res$iters), if (res$converged) "converged" else "max_iter reached",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  nm <- c(cols, "residual")
  p <- pmin(pmax(d$x$A1Freq, 1e-6), 1 - 1e-6)
  sc <- 1 / sqrt(2 * p * (1 - p))
  out <- d$x[, .(CHR, BP37, BP38, A1, A2, A1Freq, Block, Index)]
  out[, `:=`(BETA = res$BETA_STD * sc, BETA_STD = res$BETA_STD)]
  for (i in seq_along(nm)) out[, (paste0("BETA_", nm[i])) := res$PART_STD[, i] * sc]
  for (i in seq_along(nm)) out[, (paste0("BETA_STD_", nm[i])) := res$PART_STD[, i]]
  attr(out, "h2_share") <- stats::setNames(res$h2_post / sum(res$h2_post), nm)
  attr(out, "h2_share_prior") <- stats::setNames(res$h2_prior / sum(res$h2_prior), nm)
  w <- res$w
  dimnames(w) <- list(nm, c("large", "medium", "small"))
  attr(out, "w") <- w
  attr(out, "sigma2") <- res$sigma2
  attr(out, "iters") <- res$iters
  attr(out, "converged") <- res$converged
  .prs_write(out, output)
  out[]
}

#' Control settings of [mixpart()]
#'
#' @param gamma variances of the three size classes relative to sigma2.
#' @param w_init starting weight of each group (split evenly over the size classes).
#' @param h2_init heritability implied by the starting sigma2.
#' @param lam0_per_snp,lam1,theta penalty on each annotation's total weight: spike rate lam0_per_snp * M, slab rate
#'   lam1, prior slab probability theta.
#' @param penalty use the penalty (FALSE: plain MM update).
#' @param warm iterations before the penalty starts.
#' @param mm_iter MM updates of w per iteration.
#' @param max_iter,tol stop when the largest change of a posterior mean is below tol * max |posterior mean|.
#' @return a list.
#' @export
mixpart_control <- function(gamma = c(1, 0.1, 0.01), w_init = 0.02, h2_init = 0.3, lam0_per_snp = 0.5, lam1 = 1,
                            theta = 0.5, penalty = TRUE, warm = 15L, mm_iter = 20L, max_iter = 400L, tol = 2e-5) {
  list(gamma = as.numeric(gamma), w_init = as.numeric(w_init), h2_init = as.numeric(h2_init),
       lam0_per_snp = as.numeric(lam0_per_snp), lam1 = as.numeric(lam1), theta = as.numeric(theta),
       penalty = isTRUE(penalty), warm = as.integer(warm), mm_iter = as.integer(mm_iter),
       max_iter = as.integer(max_iter), tol = as.numeric(tol))
}
