#' Annotation-partitioned polygenic score weights by VIPRS with a small slab mixture (variance route)
#'
#' [viprs_annot()] with the slab a mixture of `C` normals: spike (1 - pi) and slab pi, slab class c with fixed
#' proportion p_c and variance s_c h2_j / pi, s_c = g_c / sum_c p_c g_c (so E[beta_j^2] = h2_j). The annotation
#' variance model h2_j = sum_k A_jk sigma2_k, its M-step and the MMSE allocation of the posterior mean into
#' annotation-specific components are those of [viprs_annot()]. pi is searched on the VIPRS grid (ascending, warm
#' start, PUMAS pseudo-validation, best converged point refitted), the sigma2_k carried over between grid points.
#' The sigma2_k start at `sigma2_init` (e.g. stratified LDSC tau, truncated at 0) or, without it, all equal with
#' sum_j h2_j = h2 of the naive LDSC.
#'
#' @inheritParams viprs_annot
#' @param control [viprs_annot_mix_control()].
#' @return data.table as [viprs_annot()]; attribute `mix` (proportions p_c and relative variances s_c).
#' @export
viprs_annot_mix <- function(sumstat, ld_dir, annot, ref = NULL, bin_file = NULL, resid_all = TRUE, thresh = 0,
                            fp16 = TRUE, control = viprs_annot_mix_control(), threads = 4L, output = NULL) {
  threads <- .threads(threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  ctl <- do.call(viprs_annot_mix_control, control)
  an <- data.table::as.data.table(annot)
  if (!"ID" %in% names(an)) stop("annot needs an ID column", call. = FALSE)
  cols <- setdiff(names(an), "ID")
  if ("residual" %in% cols) stop("the annotation name 'residual' is reserved", call. = FALSE)
  A <- if (length(cols)) as.matrix(an[d$x[, .(ID)], on = "ID", ..cols]) else matrix(0, nrow(d$x), 0)
  A[is.na(A)] <- 0
  storage.mode(A) <- "double"
  if (any(A < 0)) stop("annotations must be >= 0", call. = FALSE)
  message(length(cols), " annotation columns; ", sum(!d$x$ID %in% an$ID), " SNPs without annotation (set to 0)")
  nm <- c("residual", cols)
  if (!is.null(ctl$sigma2_init)) {
    s0 <- ctl$sigma2_init
    if (!is.null(names(s0))) {
      if (!setequal(names(s0), nm)) stop("sigma2_init names must be 'residual' and the annotation columns", call. = FALSE)
      s0 <- s0[nm]
    } else if (length(s0) != length(nm)) stop("sigma2_init needs ", length(nm), " values (residual first)", call. = FALSE)
    ctl$sigma2_init <- pmax(as.numeric(s0), 0)
  }
  if (is.null(bin_file)) bin_file <- file.path(ld_dir, paste0("ldbin", ctl$bin_r2, ".bin"))
  if (!file.exists(bin_file)) stop("LD-bin list not found: ", bin_file, call. = FALSE)
  t0 <- Sys.time()
  res <- viprs_annot_mix_cpp(ld_dir, d$blocks, d$keep, d$x$Z, d$x$N, d$x$freq, A, bin_file, d$start, d$n_panel,
                             resid_all = resid_all, thresh = thresh, fp16 = fp16, control = ctl, threads = threads)
  n_sup <- sum(res$extra$suppressed)
  hyp <- res$hyp
  message(sprintf("viprs_annot_mix: %d SNPs (%d suppressed in %.0f LD bins), sum k %d; pi %.4g, h2 %.4f, sigma %.4f, %d iterations (%s), %.1f s",
                  as.integer(hyp$n_snps), n_sup, res$extra$n_bins, as.integer(hyp$sum_k), hyp$pi, hyp$h2, hyp$sigma,
                  as.integer(hyp$iters), hyp$message, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out <- d$x[, .(ID, chr, bp37, bp38, A1, A2, freq)]
  out[, `:=`(BETA = res$BETA, BETA_STD = res$BETA_STD, PIP = res$PIP, H2_SNP = res$H2_SNP,
             SUPPRESSED = res$extra$suppressed)]
  sc <- ifelse(res$BETA_STD != 0, res$BETA / res$BETA_STD, 0)  # 1 / sqrt(2 p (1 - p)) per SNP
  ord <- c(cols, "residual")
  for (k in ord) out[, (paste0("BETA_", k)) := res$PART_STD[, match(k, nm)] * sc]
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
  p <- ctl$mix_p / sum(ctl$mix_p)
  attr(out, "mix") <- data.table::data.table(class = names(ctl$mix_p), p = p,
                                             s = ctl$mix_g / sum(p * ctl$mix_g))
  .prs_write(out, output)
  out[]
}

#' Control settings of [viprs_annot_mix()]
#'
#' [viprs_annot_control()] settings plus those of the slab mixture and of the start.
#' @param ... passed to [viprs_annot_control()].
#' @param mix_p proportions of the slab classes (normalized to sum 1).
#' @param mix_g relative variances of the slab classes (rescaled so that sum p_c g_c = 1).
#' @param sigma2_init start of sigma2_k: named (`residual` and the annotation columns) or ordered (residual first)
#'   numeric vector, values < 0 set to 0; `NULL`: all equal with sum_j h2_j = naive LDSC h2.
#' @param h2_floor lower bound of the per-SNP variance h2_j.
#' @return a list.
#' @export
viprs_annot_mix_control <- function(..., mix_p = c(small = 0.5, medium = 0.3, large = 0.2),
                                    mix_g = c(small = 0.01, medium = 0.1, large = 1), sigma2_init = NULL,
                                    h2_floor = 1e-12) {
  c(viprs_annot_control(...), list(mix_p = mix_p, mix_g = as.numeric(mix_g),
                                   sigma2_init = if (is.null(sigma2_init)) NULL else sigma2_init,
                                   h2_floor = as.numeric(h2_floor)))
}
