#' Stratified LD score regression on the squared-effect scale
#'
#' Ported from ldscR::ldsc.strat. With \eqn{\hat\beta_j^2 = Z_j^2 / N_j},
#' \deqn{E[\hat\beta_j^2] = \sum_c \tau_c \ell_{jc} + a / N_j,}
#' \eqn{\ell_{jc}} the LD score of annotation c (the all-ones column first) from [ld_score()] with `annot`; the
#' LDSC intercept is the coefficient \eqn{a} of \eqn{1/N_j} (1 without confounding). Z is truncated to
#' `[-z_max, z_max]`. Weights \eqn{1 / \max(1, \ell_{j1}) / (1/N_j + h^2 \ell_{j1} / M)^2}, first at `h2_init`, then
#' `n_reweight` times at the current \eqn{h^2 = \sum_c \tau_c \sum_j A_{jc}}, then fixed. With `ridge`, the
#' \eqn{\tau_c} get the penalty \eqn{\lambda u \sum_c v_c \tau_c^2} (\eqn{v_c} the variance of annotation c over
#' the SNPs, \eqn{u} the median of \eqn{(X^\top W X)_{cc} / v_c}); `lambda = NULL` chooses \eqn{\lambda} from
#' `lambda_grid` by leave-one-block-out cross-validation. Standard errors: leave-one-block-out jackknife over the
#' panel's LD blocks, the bound on \eqn{a} applied in every refit (under ridge, for the chosen \eqn{\lambda}).
#' All numerical work is in `src/ldsc_strat.cpp`.
#' @param gwas summary statistics (as in [ldsc_univ()]).
#' @param ldscore output of [ld_score()] with `annot`: L2, L2_<name> and the annotation values A_<name>.
#' @param intercept TRUE: \eqn{a} estimated within `[intercept_lower, intercept_upper]`; FALSE: fixed at 1.
#' @param intercept_lower,intercept_upper bounds for \eqn{a}.
#' @param z_max Z-scores are truncated to `[-z_max, z_max]` (not removed).
#' @param ridge ridge penalty on the \eqn{\tau_c}.
#' @param lambda ridge parameter (on the scale of `lambda_grid`); `NULL` chooses it by cross-validation.
#' @param lambda_grid candidate values for the cross-validation.
#' @inheritParams ldsc_univ
#' @return list: `H2_summary` (total and per-annotation heritability \eqn{\tau_c \sum_j A_{jc}}),
#'   `Parameter_summary` (`Intercept` = \eqn{a}, then \eqn{\tau_c}), `SE_method`, `n_truncated`, `lambda`, `cv`
#'   (cross-validation error minus its minimum, ridge with `lambda = NULL` only) and `annot_var`.
#' @export
ldsc_strat <- function(gwas, ldscore, intercept = TRUE, intercept_lower = 0.9, intercept_upper = 2, z_max = 50,
                       ridge = FALSE, lambda = NULL, lambda_grid = c(0, 10^seq(-6, 2, length.out = 49)),
                       h2_init = 0.1, n_reweight = 3L, tol = 1e-12) {
  g <- .ldsc_gwas(gwas)
  ld <- .read_table(ldscore)
  an <- sub("^A_", "", grep("^A_", names(ld), value = TRUE))
  lcol <- c("L2", paste0("L2_", an))
  if (length(setdiff(lcol, names(ld)))) stop("ldscore lacks ", paste(setdiff(lcol, names(ld)), collapse = ", "),
                                             ": give the output of ld_score() with annot", call. = FALSE)
  i <- .ldsc_at(ld, g, "ldscore")
  keep <- !is.na(i)
  g <- g[keep]; i <- i[keep]
  o <- order(g$Block, g$Index)
  g <- g[o]; i <- i[o]
  L <- as.matrix(ld[i, lcol, with = FALSE])
  A <- cbind(1, as.matrix(ld[i, paste0("A_", an), with = FALSE]))
  keep <- rowSums(!is.finite(L)) == 0
  blk <- match(g$Block[keep], unique(g$Block[keep])) - 1L
  res <- ldsc_strat_cpp(g$Z[keep], as.numeric(g$N[keep]), L[keep, , drop = FALSE], A[keep, , drop = FALSE], blk,
                        isTRUE(intercept), intercept_lower, intercept_upper, z_max, h2_init, as.integer(n_reweight),
                        isTRUE(ridge), if (is.null(lambda)) NaN else lambda, as.numeric(lambda_grid), tol)
  cats <- c("Residual", an)
  list(H2_summary = data.frame(Category = c("Total", cats), Estimate = c(res$h2, res$tau * res$annot_cov),
                               SE = c(res$h2.se, res$tau.se * res$annot_cov)),
       Parameter_summary = data.frame(Parameter = c("Intercept", cats), Estimate = c(res$intercept, res$tau),
                                      SE = c(res$intercept.se, res$tau.se)),
       SE_method = paste0("LOBO jackknife, ", res$B, " blocks"), n_truncated = res$n_truncated,
       lambda = res$lambda, cv = if (isTRUE(ridge) && is.null(lambda)) data.frame(lambda = lambda_grid, cv = res$cv),
       annot_var = res$annot_var)
}
