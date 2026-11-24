#' Trans-ancestry LD score regression (genetic covariance across two ancestries)
#'
#' Ported from ldscR::ldsc.trans. Trait 1 is summarised on panel 1 (ancestry 1), trait 2 on panel 2; the SNPs of
#' both panels are paired by [ld_score_cross()], whose `sign` turns trait 2's z-scores to the alleles of panel 1.
#' With the cross LD score \eqn{L_{12}} and \eqn{\ell_j = L_{12,j} \sqrt{N_{1j} N_{2j}} / M}, \eqn{z_{1j} z_{2j}}
#' is regressed on \eqn{(1, \ell_j)}: weights \eqn{1 / v_j},
#' \eqn{v_j = (1 + L_{1j} N_{1j} h_{21} / M)(1 + L_{2j} N_{2j} h_{22} / M)}, then `n_reweight` times
#' \eqn{1 / (v_j + 2 \mu_j^2)}. Standard errors: block bootstrap as [ldsc_univ()]. All numerical work is in
#' `src/ldsc_trans.cpp`.
#' @param gwas1 summary statistics of trait 1 on panel 1; `gwas2` of trait 2 on panel 2.
#' @param ldscore1,ldscore2 [ld_score()] of panel 1 and of panel 2.
#' @param ldscore_cross [ld_score_cross()] of panel 1 and panel 2.
#' @param h21,h22 heritabilities of the two traits (in the weights).
#' @inheritParams ldsc_univ
#' @return data.frame: ecov, ecov.se, gcov, gcov.se, M; `$bootstrap` as in [ldsc_univ()].
#' @export
ldsc_trans <- function(gwas1, gwas2, ldscore1, ldscore2, ldscore_cross, h21, h22, n_boot = 0L, nblock = 500L,
                       n_reweight = 2L, tol = 1e-12) {
  g1 <- .ldsc_gwas(gwas1, "gwas1")
  g2 <- .ldsc_gwas(gwas2, "gwas2")
  l1 <- .read_table(ldscore1)
  l2 <- .read_table(ldscore2)
  x <- .read_table(ldscore_cross)
  miss <- setdiff(c("Block1", "Index1", "Block2", "Index2", "sign", "L2_cross"), names(x))
  if (length(miss)) stop("ldscore_cross lacks ", paste(miss, collapse = ", "), ": give the output of ld_score_cross()",
                         call. = FALSE)
  k1 <- x[, .(Block = as.integer(Block1), Index = as.integer(Index1))]
  k2 <- x[, .(Block = as.integer(Block2), Index = as.integer(Index2))]
  i1 <- g1[k1, on = .(Block, Index), which = TRUE]
  i2 <- g2[k2, on = .(Block, Index), which = TRUE]
  g <- data.table::data.table(Block = k1$Block, Index = k1$Index, z1 = g1$Z[i1], N1 = as.numeric(g1$N[i1]),
                              z2 = x$sign * g2$Z[i2], N2 = as.numeric(g2$N[i2]),
                              L1 = l1$L2[l1[k1, on = .(Block, Index), which = TRUE]],
                              L2 = l2$L2[l2[k2, on = .(Block, Index), which = TRUE]], L12 = x$L2_cross)
  g <- g[is.finite(z1) & is.finite(z2) & is.finite(L1) & is.finite(L2) & is.finite(L12)]
  data.table::setorder(g, Block, Index)
  counts <- .ldsc_counts(nrow(g), nblock, n_boot)
  res <- ldsc_trans_cpp(g$z1, g$z2, g$N1, g$N2, g$L1, g$L2, g$L12, h21, h22, as.integer(nblock), counts,
                        as.integer(n_reweight), tol)
  .ldsc_result(res, c("ecov", "ecov.se", "gcov", "gcov.se", "M"), counts)
}
