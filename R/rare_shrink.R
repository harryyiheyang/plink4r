# Effect-size shrinkage of logistic association results (rare_variant_calibration.md, sec. 11-14).
#
# Per-allele prior of a variant with minor allele frequency p: a mixture of N(0, sigma_k^2(p)) with
# sigma_k^2(p) = sigma0sq * grid_k * (2 p (1 - p))^alpha, grid = (0, 1, 10, 100), weights `weights`, plus a wide
# component N(0, wide_sd^2) (not scaled by MAF) with fixed weight `wide_weight` (the grid weights are scaled by
# 1 - wide_weight). The likelihood enters as a Gaussian summary bhat ~ N(beta, s^2).

.rv_prior <- function(prior = NULL) {
  d <- list(sigma0sq = 0.05^2 / (2 * 0.25 * 0.75)^-0.38,  # sd 0.05 at MAF 0.25 (fallback, nothing learned)
            weights = c(0.7, 0.15, 0.1, 0.05), alpha = -0.38, grid = c(0, 1, 10, 100),
            wide_sd = 1.5, wide_weight = 0.02, learned = FALSE)
  if (is.null(prior)) return(d)
  if (!is.list(prior)) stop("prior must be a list (see gwas()).", call. = FALSE)
  bad <- setdiff(names(prior), names(d))
  if (length(bad)) stop("Unknown prior field(s): ", paste(bad, collapse = ", "), call. = FALSE)
  d[names(prior)] <- prior
  if (length(d$weights) != length(d$grid) || any(d$weights < 0) || sum(d$weights) <= 0)
    stop("prior$weights must be non-negative, one per grid component.", call. = FALSE)
  d$weights <- d$weights / sum(d$weights)
  if (!(d$wide_weight >= 0 && d$wide_weight < 1)) stop("prior$wide_weight must be in [0, 1).", call. = FALSE)
  d
}

# Prior sds (n x K, the wide component last) and weights (K) at minor allele frequencies maf.
.rv_prior_sd <- function(prior, maf) {
  sc <- (2 * maf * (1 - maf))^prior$alpha
  sd <- sqrt(outer(sc * prior$sigma0sq, prior$grid))
  w <- prior$weights * (1 - prior$wide_weight)
  if (prior$wide_weight > 0) { sd <- cbind(sd, prior$wide_sd); w <- c(w, prior$wide_weight) }
  list(sd = sd, w = w)
}

# Mixture posterior mean and sd from bhat ~ N(beta, s2), prior sds sd (n x K), weights w (K):
# lbf_k = log N(bhat; 0, s2 + sd_k^2) - log N(bhat; 0, s2), w~_k propto w_k exp(lbf_k), component posterior
# N(bhat sd_k^2 / (sd_k^2 + s2), sd_k^2 s2 / (sd_k^2 + s2)); mean = sum w~_k m_k, var by the law of total variance.
.post_mix <- function(bhat, s2, sd, w) {
  n <- length(bhat); K <- ncol(sd)
  v <- sd^2
  tot <- v + s2
  lw <- matrix(log(w), n, K, byrow = TRUE) - 0.5 * log(tot / s2) - 0.5 * bhat^2 * (1 / tot - 1 / s2)
  lw <- lw - apply(lw, 1L, max)
  pw <- exp(lw); pw <- pw / rowSums(pw)
  shr <- v / tot
  m <- bhat * shr
  vk <- shr * s2
  mean <- rowSums(pw * m)
  var <- pmax(rowSums(pw * (vk + m^2)) - mean^2, 0)
  bad <- !is.finite(bhat) | !is.finite(s2) | s2 <= 0
  mean[bad] <- NA_real_; var[bad] <- NA_real_
  list(mean = mean, sd = sqrt(var))
}

# BETA_POST / SE_POST of a logistic gwas() result: the likelihood is the carrier step's (beta_carrier,
# 1 / info_carrier) where it was computed, else the one-step (BETA, SE^2).
.rv_shrink <- function(out, res, prior) {
  bh <- out$BETA; s2 <- out$SE^2
  k <- which(!is.na(res$beta_carrier))
  bh[k] <- res$beta_carrier[k]; s2[k] <- 1 / res$info_carrier[k]
  maf <- pmin(out$A1_FREQ, 1 - out$A1_FREQ)
  ps <- .rv_prior_sd(prior, maf)
  pm <- .post_mix(bh, s2, ps$sd, ps$w)
  out[, `:=`(BETA_POST = pm$mean, SE_POST = pm$sd)]
  invisible(out)
}
