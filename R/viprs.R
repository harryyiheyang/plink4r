#' Polygenic score weights by VIPRS on the eigen LD panel
#'
#' Variational spike-and-slab regression (VIPRS 0.1.5 updates) on the SBayesRC eigen likelihood, with one set of
#' hyperparameters (pi, tau, sigma) for the genome; ported from viprs_eigen (Yihe Yang). Per LD block the data are
#' w = Lambda^{-1/2} U' b and the LD is R = V V', V = U Lambda^{1/2} (rows of the SNPs with data only), so R is never
#' formed; each coordinate update costs O(k) (k = eigen dimension of the block). Blocks are swept in parallel.
#'
#' Hyperparameter search: grid search with warm start (viprs_eigen `"path"`). The pi grid is run ascending (strong to
#' weak shrinkage), each point warm-started from the previous one, tau and sigma by EM, scored by PUMAS pseudo-validation
#' (training / validation split of the summary statistics); the search stops after `stop_after` consecutive drops of the
#' pseudo R2. A point that does not converge is skipped in the choice; the best converged point is refitted on all data.
#'
#' The pi grid is VIPRS's official one (`pi_steps` values, log10 spaced from max(10 / M, 1e-5) to min(1e4 / M, 0.2),
#' extended upward by `grid_extend` points), or `pi_grid` when given (the cgrid setting: `pi_grid = c * r^(kmin:kmax)`).
#'
#' @param sumstat output of [summarystat_tidy()] or [summarystat_impute()] (data.table or file): columns ID, A1, A2,
#'   Z, N (optionally r2: rows with r2 <= 0 carry no data and are dropped). Rows are matched to the panel by ID; A1 / A2
#'   swapped relative to the panel alt / ref flips Z; other rows are dropped.
#' @param ld_dir directory with block<N>.eigen.bin (SBayesRC or plink4r layout) and the panel SNP .parquet.
#' @param ref panel SNP table in place of the .parquet in `ld_dir` (data.frame or parquet path): ID, chr, bp37 and / or
#'   bp38, ref, alt, AF_PanUKB, block, ld_index.
#' @param thresh eigenvalue variance kept per block; <= 0 (default) keeps every stored eigenpair (as viprs_eigen).
#' @param fp16 store V as IEEE float16 (half the memory; the data terms are computed from the rounded rows).
#' @param control [viprs_control()].
#' @param threads OpenMP threads (0: [plink4r_resources()]).
#' @param output optional output path (.parquet or text).
#' @return data.table over the fitted SNPs: ID, chr, bp37, bp38, A1, A2, freq (A1 frequency of the panel), BETA (per
#'   A1 allele, BETA_STD / sqrt(2 freq (1 - freq))), BETA_STD, PIP; attributes `hyp` (hyperparameters, iterations,
#'   counts), `grid` (one row per grid point), `best_pseudo_r2` and `pumas_xi`. ID, A1 and BETA are
#'   the columns [geno_score()] reads.
#' @export
viprs <- function(sumstat, ld_dir, ref = NULL, thresh = 0, fp16 = TRUE, control = viprs_control(), threads = 4L,
                  output = NULL) {
  threads <- .threads(threads)
  d <- .prs_inputs(sumstat, ref, ld_dir)
  ctl <- do.call(viprs_control, control)
  t0 <- Sys.time()
  res <- viprs_cpp(ld_dir, d$blocks, d$keep, d$x$Z, d$x$N, d$x$freq, thresh = thresh, fp16 = fp16, control = ctl,
                   threads = threads)
  hyp <- res$hyp
  message(sprintf("viprs: %d SNPs, sum k %d; pi %.4g, h2 %.4f, sigma %.4f, %d iterations (%s), %.1f s",
                  as.integer(hyp$n_snps), as.integer(hyp$sum_k), hyp$pi, hyp$h2, hyp$sigma, as.integer(hyp$iters),
                  hyp$message, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out <- d$x[, .(ID, chr, bp37, bp38, A1, A2, freq)]
  out[, `:=`(BETA = res$BETA, BETA_STD = res$BETA_STD, PIP = res$PIP)]
  attr(out, "hyp") <- hyp
  attr(out, "grid") <- data.table::as.data.table(res$extra$grid)
  attr(out, "best_pseudo_r2") <- res$extra$best_pseudo_r2
  attr(out, "pumas_xi") <- res$extra$pumas_xi
  .prs_write(out, output)
  out[]
}

#' Control settings of [viprs()]
#'
#' Defaults are those of viprs_eigen (engine.Fitter, viprs_eigen_fast.py).
#' @param max_iter iterations of the final refit.
#' @param grid_max_iter iterations per grid point.
#' @param stop_after stop a pi path after this many consecutive drops of the pseudo R2 (0 = never).
#' @param pi_steps points of the official pi grid.
#' @param grid_extend extra grid points above the official upper end, same log spacing (0 = official grid).
#' @param pi_grid user pi grid (replaces the official grid and `grid_extend`), e.g. `c * sqrt(3)^(-6:3)`.
#' @param prop_train,seed PUMAS training proportion and random seed.
#' @param use_active,full_every,act_tol active set: partial sweeps skip SNPs that moved by <= act_tol at the last
#'   full sweep; a full sweep every `full_every` iterations.
#' @param refresh_every recompute s = V' eta exactly every this many iterations.
#' @param x_tol,f_tol,min_iter,patience convergence (VIPRS 0.1.5): |dELBO| <= f_tol, max |d eta| < x_tol, or the
#'   LD-weighted effects stable for more than `patience` iterations; only after `min_iter` iterations.
#' @param sigma_fixed residual variance used when the sigma update is <= 0.
#' @param h2_min,h2_max clip of the LDSC start h2 = M (mean chi2 - 1) / (N mean l).
#' @param pi_lo_count,pi_lo_min,pi_hi_count,pi_hi_max,pi_cap official grid ends max(pi_lo_count / M, pi_lo_min),
#'   min(pi_hi_count / M, pi_hi_max); extended points above pi_cap are dropped.
#' @return a list.
#' @export
viprs_control <- function(max_iter = 1000L, grid_max_iter = 200L, stop_after = 2L, pi_steps = 10L, grid_extend = 2L,
                          pi_grid = NULL, prop_train = 0.8, seed = 7209, use_active = TRUE, full_every = 5L,
                          act_tol = 1e-7, refresh_every = 20L, x_tol = 1e-6, f_tol = 1e-6, min_iter = 3L,
                          patience = 5L, sigma_fixed = 0.95, h2_min = 0.01, h2_max = 0.99, pi_lo_count = 10,
                          pi_lo_min = 1e-5, pi_hi_count = 1e4, pi_hi_max = 0.2, pi_cap = 0.5) {
  list(max_iter = as.integer(max_iter), grid_max_iter = as.integer(grid_max_iter), stop_after = as.integer(stop_after),
       pi_steps = as.integer(pi_steps), grid_extend = as.integer(grid_extend),
       pi_grid = if (is.null(pi_grid)) NULL else as.numeric(pi_grid),
       prop_train = as.numeric(prop_train), seed = as.numeric(seed), use_active = isTRUE(use_active),
       full_every = as.integer(full_every), act_tol = as.numeric(act_tol), refresh_every = as.integer(refresh_every),
       x_tol = as.numeric(x_tol), f_tol = as.numeric(f_tol), min_iter = as.integer(min_iter),
       patience = as.integer(patience), sigma_fixed = as.numeric(sigma_fixed), h2_min = as.numeric(h2_min),
       h2_max = as.numeric(h2_max), pi_lo_count = as.numeric(pi_lo_count), pi_lo_min = as.numeric(pi_lo_min),
       pi_hi_count = as.numeric(pi_hi_count), pi_hi_max = as.numeric(pi_hi_max), pi_cap = as.numeric(pi_cap))
}
