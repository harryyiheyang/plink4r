#' Bivariate LD score regression (genetic covariance of two traits)
#'
#' Ported from ldscR::ldsc.bicov. Both traits' summary statistics come from the same panel and are matched with
#' its LD scores by (Block, Index). With \eqn{\ell_j = L2_j \sqrt{N_{1j} N_{2j}} / M}, \eqn{z_{1j} z_{2j}} is
#' regressed on \eqn{(1, \ell_j)}: weights \eqn{1 / \max(1, L2_j) / v_j},
#' \eqn{v_j = (1 + L2_j N_{1j} h_{21} / M)(1 + L2_j N_{2j} h_{22} / M)}, then `n_reweight` times
#' \eqn{1 / \max(1, L2_j) / (v_j + 2 \mu_j^2)}, \eqn{\mu_j} the current fitted value. Standard errors: block
#' bootstrap as [ldsc_univ()]. All numerical work is in `src/ldsc_bicov.cpp`.
#' @param gwas1,gwas2 summary statistics of the two traits (as in [ldsc_univ()]).
#' @param ldscore output of [ld_score()] on the panel.
#' @param h21,h22 heritabilities of the two traits (in the weights), e.g. from [ldsc_univ()].
#' @inheritParams ldsc_univ
#' @return data.frame: ecov (intercept), ecov.se, gcov (genetic covariance), gcov.se, M; `$bootstrap` as in
#'   [ldsc_univ()].
#' @export
ldsc_bicov <- function(gwas1, gwas2, ldscore, h21, h22, n_boot = 0L, nblock = 500L, n_reweight = 2L, tol = 1e-12) {
  g1 <- .ldsc_gwas(gwas1, "gwas1")
  g2 <- .ldsc_gwas(gwas2, "gwas2")
  ld <- .read_table(ldscore)
  i2 <- .ldsc_at(g2, g1, "gwas2")
  g <- g1[, .(Block, Index, z1 = Z, N1 = as.numeric(N), z2 = g2$Z[i2], N2 = as.numeric(g2$N[i2]))]
  g[, L2 := ld$L2[.ldsc_at(ld, g1, "ldscore")]]
  g <- g[is.finite(z2) & is.finite(L2)]
  data.table::setorder(g, Block, Index)
  counts <- .ldsc_counts(nrow(g), nblock, n_boot)
  res <- ldsc_bicov_cpp(g$z1, g$z2, g$N1, g$N2, g$L2, h21, h22, as.integer(nblock), counts, as.integer(n_reweight), tol)
  .ldsc_result(res, c("ecov", "ecov.se", "gcov", "gcov.se", "M"), counts)
}
