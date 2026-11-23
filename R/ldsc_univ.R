#' Univariate LD score regression (SNP heritability)
#'
#' Ported from ldscR::ldsc.univ. The summary statistics and the LD scores are matched by the panel key (Block,
#' Index); with \eqn{M} SNPs, \eqn{\ell_j = L2_j N_j / M} and \eqn{Z_j^2} is regressed on \eqn{(1, \ell_j)} by
#' weighted least squares, weights \eqn{1 / \max(1, L2_j) / (1 + \ell_j h^2)^2}, first at `h2_init` and then
#' `n_reweight` times at the current estimate. Standard errors: block bootstrap over `nblock` contiguous groups
#' of SNPs (in panel order), `n_boot` replicates. All numerical work is in `src/ldsc_univ.cpp`.
#' @param gwas output of [summarystat_tidy()] or [sumstat_impute()] (data.table or file): Z, N, Block, Index.
#' @param ldscore output of [ld_score()] on the same panel (its column `L2`).
#' @param n_boot bootstrap replicates; 0 gives point estimates only.
#' @param nblock contiguous groups for the bootstrap.
#' @param h2_init heritability in the first weights.
#' @param n_reweight reweighting rounds after the first fit.
#' @param tol rank tolerance of the weighted least-squares solve.
#' @return data.frame: intercept, intercept.se, h2, h2.se, M. With `n_boot > 0`, `$bootstrap` holds the
#'   replicate estimates (`intercept`, `h2`) and `block.counts` (replicates x groups).
#' @export
ldsc_univ <- function(gwas, ldscore, n_boot = 0L, nblock = 500L, h2_init = 0.1, n_reweight = 3L, tol = 1e-12) {
  g <- .ldsc_gwas(gwas)
  ld <- .read_table(ldscore)
  g[, L2 := ld$L2[.ldsc_at(ld, g, "ldscore")]]
  g <- g[is.finite(L2)]
  data.table::setorder(g, Block, Index)
  counts <- .ldsc_counts(nrow(g), nblock, n_boot)
  res <- ldsc_univ_cpp(g$Z, as.numeric(g$N), g$L2, as.integer(nblock), counts, h2_init, as.integer(n_reweight), tol)
  .ldsc_result(res, c("intercept", "intercept.se", "h2", "h2.se", "M"), counts)
}
