#' Annotation-partitioned polygenic scores (PartitionPRS): variational selection, short Gibbs with Rao-Blackwellization
#'
#' Research Design D.1.2. Model of [viprs_annot()]: h2_j = sigma2_0 + sum_k A_jk sigma2_k (sigma2_0 the intercept on
#' every SNP; its component is the residual PRS), spike and slab with global pi, the slab a mixture with fixed shares
#' (`mix_prop`, small / medium / large = 6 : 3 : 1) and relative variances (`mix_var`, 1 : 10 : 100).
#' Step 1, selection: [viprs_annot()] (variational, pi grid ascending with warm start, PUMAS pseudo-validation, the
#' moment M-step of D.1.2) chooses pi and gives the start of sigma2_k, sigma and beta.
#' Step 2, fit: a short Gibbs sampler at that pi, warm-started there (`iter` sweeps, the first `burn` discarded);
#' sigma2_k and sigma by SAEM (the moment M-step of D.1.2 on the averaged Rao-Blackwellized statistics). The weights
#' are the Rao-Blackwellized posterior means b_j; the partitioned weights are the MMSE allocation
#' b_j^(k) = w_jk b_j, w_jk = A_jk sigma2_k / h2_j at the final sigma2, so per SNP the BETA_<...> columns sum to BETA
#' and the partitioned scores of [geno_score()] sum to the total score.
#' Heritability: `h2_annot` = (sum_j A_jk) sigma2_k (the model heritability of D.1.2, additive) and `h2_shapley` =
#' b^(k)' R b (the Shapley value of the PRS variance v(S) = b_S' R b_S; the columns sum to b' R b).
#' LD-bin skip as [viprs_annot()] (r2 > `bin_r2`, unannotated members of a bin with an annotated member fixed at 0).
#'
#' @inheritParams viprs_annot
#' @param control [partprs_control()] (the Gibbs step).
#' @param vi_control [viprs_annot_control()] (the selection step).
#' @return data.table as [viprs_annot()] (BETA, BETA_STD, PIP, H2_SNP, SUPPRESSED, BETA_<annotation>,
#'   BETA_residual and their BETA_STD_ columns); attributes `sigma2`, `h2_annot`, `h2_shapley`, `hyp`, `trace`
#'   (sigma2_k, sigma and the beta' R beta of each draw) and `vi` (the attributes of the selection fit).
#' @export
partprs <- function(sumstat, ld_dir, annot, ref = NULL, bin_file = NULL, skip_annot = NULL, sigma2_init = NULL,
                    thresh = 0, fp16 = TRUE, control = partprs_control(), vi_control = viprs_annot_control(exact_em = FALSE),
                    threads = 4L, output = NULL) {
  threads <- .threads(threads)
  ctl <- do.call(partprs_control, control)
  vctl <- do.call(viprs_annot_control, vi_control)
  vctl$mix_prop <- ctl$mix_prop; vctl$mix_var <- ctl$mix_var
  vi <- viprs_annot(sumstat, ld_dir, annot, ref = ref, bin_file = bin_file, skip_annot = skip_annot,
                    sigma2_init = sigma2_init, thresh = thresh, fp16 = fp16, control = vctl, threads = threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  A <- .prs_annot(annot, d$x)
  cols <- colnames(A)
  hv <- attr(vi, "hyp")
  t0 <- Sys.time()
  res <- partprs_cpp(.ld_block_file(ld_dir, d$blocks), d$blocks, d$keep, d$x$Z, d$x$N, d$x$A1Freq, A,
                     suppressed = vi$SUPPRESSED, pi = hv$pi, sigma2 = unname(attr(vi, "sigma2")), sigma = hv$sigma,
                     beta = vi$BETA_STD, h2_ref = hv$h2_ldsc_init, thresh = thresh, fp16 = fp16, control = ctl,
                     threads = threads)
  hyp <- res$hyp
  message(sprintf("partprs: %d SNPs, pi %.4g (viprs_annot), Gibbs %d sweeps (%d burn-in), h2 %.4f, sigma %.4f, %.1f s",
                  as.integer(hyp$n_snps), hyp$pi, as.integer(hyp$iter), as.integer(hyp$burn), hyp$h2, hyp$sigma,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  nm <- c("residual", cols)
  out <- d$x[, .(CHR, BP37, BP38, A1, A2, A1Freq, Block, Index)]
  out[, `:=`(BETA = res$BETA, BETA_STD = res$BETA_STD, PIP = res$PIP, H2_SNP = res$H2_SNP, SUPPRESSED = vi$SUPPRESSED)]
  sc <- ifelse(res$BETA_STD != 0, res$BETA / res$BETA_STD, 0)
  ord <- c(cols, "residual")
  for (k in ord) out[, (paste0("BETA_", k)) := res$PART_STD[, match(k, nm)] * sc]
  for (k in ord) out[, (paste0("BETA_STD_", k)) := res$PART_STD[, match(k, nm)]]
  tr <- res$trace
  colnames(tr) <- c(nm, "sigma", "bRb")
  attr(out, "hyp") <- hyp
  attr(out, "sigma2") <- stats::setNames(res$sigma2, nm)
  attr(out, "h2_annot") <- stats::setNames(res$h2_annot, nm)
  attr(out, "h2_shapley") <- stats::setNames(res$h2_shapley, nm)
  attr(out, "trace") <- data.table::as.data.table(tr)
  attr(out, "vi") <- attributes(vi)[c("hyp", "sigma2", "h2_annot", "grid", "grid_sigma2", "best_pseudo_r2", "n_bins",
                                      "n_suppressed")]
  .prs_write(out, output)
  out[]
}

#' Control settings of [partprs()] (the Gibbs step)
#'
#' @param iter,burn Gibbs sweeps and the first `burn` discarded (SAEM step 1 during burn-in).
#' @param mix_prop,mix_var slab classes small / medium / large: prior shares and relative variances (also used by the
#'   selection step).
#' @param a0,b0 Inverse-Gamma(a0, b0) prior of each sigma2_k in the M-step of D.1.2.
#' @param var_floor lower bound of every sigma2_k and h2_j, as a fraction of the mean h2_j of the naive LDSC h2.
#' @param sigma_fixed residual variance used when 1 - 2 beta'b + beta'R beta of a draw is not positive.
#' @param af_clip allele frequencies clipped to [af_clip, 1 - af_clip] for BETA = BETA_STD / sqrt(2 p (1 - p)).
#' @param seed seed of the sampler.
#' @return a list.
#' @export
partprs_control <- function(iter = 400L, burn = 100L, mix_prop = c(0.6, 0.3, 0.1), mix_var = c(0.01, 0.1, 1),
                            a0 = 0, b0 = 0, var_floor = 1e-3, sigma_fixed = 0.95, af_clip = 1e-6, seed = 7209) {
  list(iter = as.integer(iter), burn = as.integer(burn), mix_prop = as.numeric(mix_prop), mix_var = as.numeric(mix_var),
       a0 = as.numeric(a0), b0 = as.numeric(b0), var_floor = as.numeric(var_floor), sigma_fixed = as.numeric(sigma_fixed),
       af_clip = as.numeric(af_clip), seed = as.numeric(seed))
}
