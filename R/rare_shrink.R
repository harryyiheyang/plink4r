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

# Block index of each variant (NA outside the blocks / non-autosomal). blocks: 37 / 38 (ld_blocks_default()) or
# a block file (header plus chr start stop).
.rv_block_of <- function(chr, pos, blocks) {
  f <- if (is.numeric(blocks)) ld_blocks_default(blocks) else blocks
  bl <- .ldb_read_blocks(f)
  data.table::setorder(bl, CHR, start)
  chr <- .ldb_chr(chr)
  b <- rep(NA_integer_, length(chr))
  for (cc in unique(bl$CHR)) {
    k <- which(chr == cc); if (!length(k)) next
    bc <- bl[CHR == cc]
    i <- findInterval(pos[k], bc$start)
    ok <- i > 0L
    ok[ok] <- pos[k][ok] <= bc$end[i[ok]]
    b[k[ok]] <- bc$src[i[ok]]
  }
  b
}

# Learn sigma0sq and the grid weights of the prior from the common variants (MAF >= common_maf) of a gwas()
# result (rare_variant_calibration.md sec. 13-14). Fields already in `fixed` (user-given) are kept.
# - sigma0sq: in each LD block the variant with the largest |z| (z = BETA / SE) and P < 5e-8 (|z| > zc); with
#   c_i = (2 p_i q_i)^alpha / SE_i^2, z_i ~ N(0, 1 + tau2 c_i) truncated to |z_i| > zc, tau2 by one-dimensional
#   maximum likelihood (the Gaussian log-likelihood minus log P(|z_i| > zc) per hit). The middle grid component
#   is put at the hits' variance: sigma0sq = tau2 / grid[3].
# - weights: one random common variant per block, EM on the fixed grid (wide component left out).
.rv_learn <- function(out, prior, fixed, blocks, common_maf = 0.01, zc = stats::qnorm(5e-8 / 2, lower.tail = FALSE),
                      min_hits = 10L, min_blocks = 50L) {
  info <- list(hits = 0L, random = 0L)
  maf <- pmin(out$A1_FREQ, 1 - out$A1_FREQ)
  z <- out$BETA / out$SE
  ok <- which(is.finite(z) & maf >= common_maf & out$SE > 0)
  if (!length(ok)) return(c(prior, list(learn = info)))
  blk <- .rv_block_of(out$CHR[ok], out$POS[ok], blocks)
  ok <- ok[!is.na(blk)]; blk <- blk[!is.na(blk)]
  if (!("sigma0sq" %in% fixed)) {
    o <- order(blk, -abs(z[ok]))
    top <- ok[o][!duplicated(blk[o])]
    top <- top[abs(z[top]) > zc]
    info$hits <- length(top)
    if (length(top) >= min_hits) {
      zz <- z[top]; ci <- (2 * maf[top] * (1 - maf[top]))^prior$alpha / out$SE[top]^2
      nll <- function(lt) { v <- 1 + exp(lt) * ci
        -sum(-0.5 * log(v) - 0.5 * zz^2 / v - (log(2) + stats::pnorm(zc / sqrt(v), lower.tail = FALSE, log.p = TRUE))) }
      s0 <- max(mean(zz^2) - 1, 1e-6) / stats::median(ci)
      lt <- stats::optimize(nll, log(s0) + c(-15, 8))$minimum
      prior$sigma0sq <- exp(lt) / prior$grid[min(3L, length(prior$grid))]
      info$tau2 <- exp(lt)
    }
  }
  if (!("weights" %in% fixed)) {
    # one variant per block, chosen at random without touching the caller's RNG stream
    if (exists(".Random.seed", globalenv())) { old <- get(".Random.seed", globalenv())
      on.exit(assign(".Random.seed", old, globalenv()))
    } else on.exit(rm(".Random.seed", envir = globalenv()))
    set.seed(20261116L)
    sh <- sample.int(length(ok))
    pick <- ok[sh][!duplicated(blk[sh])]
    info$random <- length(pick)
    if (length(pick) >= min_blocks) {
      sd <- sqrt(outer((2 * maf[pick] * (1 - maf[pick]))^prior$alpha * prior$sigma0sq, prior$grid))
      s2 <- out$SE[pick]^2; bh <- out$BETA[pick]
      L <- -0.5 * log(sd^2 + s2) - 0.5 * bh^2 / (sd^2 + s2)
      w <- rep(1 / ncol(L), ncol(L))
      for (it in seq_len(1000L)) {
        lw <- L + matrix(log(pmax(w, 1e-300)), nrow(L), ncol(L), byrow = TRUE)
        lw <- exp(lw - apply(lw, 1L, max)); lw <- lw / rowSums(lw)
        wn <- colMeans(lw)
        if (max(abs(wn - w)) < 1e-8) { w <- wn; break }
        w <- wn
      }
      prior$weights <- w
    }
  }
  prior$learned <- info$hits >= min_hits || info$random >= min_blocks
  c(prior, list(learn = info))
}
