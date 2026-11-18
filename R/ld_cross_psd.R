#' Replace the cross block of an LD matrix and keep the matrix positive semi-definite
#'
#' For an LD matrix whose variants split into a kept part 1 and a replaced part 2: R11 stays as it is, R22 is the
#' new block (e.g. `cor()` of the corrected genotypes) and `C12` an estimate of the new cross block (e.g. the
#' correlation of the observed or filled part-1 calls with the corrected part-2 genotypes). An estimated C12 does
#' not in general make `[R11, C12; t(C12), R22]` PSD. With R11 = W W' (W = U1 diag(sqrt(lambda)), the part-1 rows
#' of the eigen.bin factor, or from `R11`) and B = R22^(1/2), every cross block C = W G B with ||G||_2 <= 1 gives
#' `[R11, C; C', R22] = diag(W, B) [I, G; G', I] diag(W, B)' >= 0`, and every PSD completion has this form. The
#' function takes G = W^+ C12 B^+, clips its singular values at 1 - eps and returns C = W G B. R11 and R22 are not
#' changed; the cost is one SVD of W (p1 x k), one eigen decomposition of R22 and one SVD of G (k x p2).
#' @param C12 p1 x p2 estimate of the cross block.
#' @param R22 p2 x p2 new block (symmetric PSD).
#' @param W p1 x k factor with R11 = W W' (preferred: no p1 x p1 work); or give `R11`.
#' @param R11 p1 x p1 kept block, used when `W` is NULL.
#' @param eps singular values of G are clipped at 1 - eps (eps > 0 keeps the result positive definite on the
#'   joint range).
#' @param tol relative cut for the singular values of W and the eigenvalues of R22 treated as zero.
#' @param full TRUE also returns the assembled (p1 + p2) x (p1 + p2) matrix.
#' @return list(C12 = the PSD cross block, sv_max = largest singular value of G before clipping,
#'   n_clipped = number of clipped singular values, change = Frobenius norm of C12 - input, R = full matrix
#'   when `full`).
#' @export
ld_cross_psd <- function(C12, R22, W = NULL, R11 = NULL, eps = 1e-6, tol = 1e-8, full = FALSE) {
  C12 <- as.matrix(C12); R22 <- as.matrix(R22)
  if (is.null(W)) {
    if (is.null(R11)) stop("give W or R11.", call. = FALSE)
    e1 <- eigen(as.matrix(R11), symmetric = TRUE)
    k1 <- e1$values > tol * max(e1$values)
    W <- e1$vectors[, k1, drop = FALSE] %*% diag(sqrt(e1$values[k1]), sum(k1))
  }
  W <- as.matrix(W)
  p1 <- nrow(W); p2 <- nrow(R22)
  if (!identical(dim(C12), c(p1, p2))) stop("C12 must be ", p1, " x ", p2, ".", call. = FALSE)
  if (ncol(R22) != p2) stop("R22 must be square.", call. = FALSE)
  if (!(eps > 0 && eps < 1)) stop("eps must lie in (0, 1).", call. = FALSE)
  sw <- svd(W)
  kw <- sw$d > tol * max(sw$d)
  e2 <- eigen((R22 + t(R22)) / 2, symmetric = TRUE)
  k2 <- e2$values > tol * max(e2$values)
  V2 <- e2$vectors[, k2, drop = FALSE]; l2 <- sqrt(e2$values[k2])
  # coordinates: W = A D Bw', R22^(1/2) = V2 L V2'; work in the reduced bases
  A <- sw$u[, kw, drop = FALSE]; d <- sw$d[kw]
  H <- (crossprod(A, C12) %*% V2) / outer(d, l2)  # = Bw' G V2 in the reduced bases
  sh <- svd(H)
  s_new <- pmin(sh$d, 1 - eps)
  Ht <- sh$u %*% (s_new * t(sh$v))
  Cn <- A %*% ((d * Ht) * rep(l2, each = length(d))) %*% t(V2)
  dimnames(Cn) <- dimnames(C12)
  out <- list(C12 = Cn, sv_max = if (length(sh$d)) max(sh$d) else 0, n_clipped = sum(sh$d > 1 - eps),
              change = sqrt(sum((Cn - C12)^2)))
  if (full) {
    R11v <- tcrossprod(W)
    out$R <- rbind(cbind(R11v, Cn), cbind(t(Cn), R22))
  }
  out
}
