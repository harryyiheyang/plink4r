# Fast checks of the rare-variant calibration helpers (no genotype file needed).
library(plink4r)
ns <- asNamespace("plink4r")
f <- ns$.p_robust_1df

# 1. small-count robust P: never below the normal P, any W (also W >= M) and non-integer M
N <- 10000
for (m in c(5, 5.1, 7.3, 19.5, 190.2)) {
  W <- c(seq(0, 300, 0.37), m - 1e-6, m, m + 1e-6)
  p <- f(W, N, m / (2 * N))
  stopifnot(all(p >= pchisq(W, 1, lower.tail = FALSE) - 1e-15), all(p <= 1), !anyNA(p))
  stopifnot(all(p[W >= m * (1 - 1e-9)] == pchisq(W[W >= m * (1 - 1e-9)], 1, lower.tail = FALSE)))
}
# the median correction is kept: at W = 0.455 (chi-square median) the P is above 0.5 for small M
stopifnot(f(0.455, N, 5 / (2 * N)) > 0.5)
# below MAC 5 and with adjust = FALSE: the normal P
stopifnot(f(3, N, 4 / (2 * N)) == pchisq(3, 1, lower.tail = FALSE), f(3, N, 0.01, FALSE) == pchisq(3, 1, lower.tail = FALSE))

# 2. variance-component saddlepoint: with V = I it is chi-square(q); Cauchy of identical P returns that P
set.seed(1); q <- 3; b <- matrix(rnorm(15, 0, 2), 5)
V <- matrix(rep(as.vector(diag(q)), each = 5), 5)
pv <- pchisq(rowSums(b^2), q, lower.tail = FALSE)
got <- ns$int_cauchy_cpp(b, V, matrix(pv, 5, q), 1L)
stopifnot(all(abs(log(got) - log(pv)) < 0.02))

# 3. mixture posterior: the multivariate C++ code with S = 1 equals the scalar closed form
pr <- ns$.rv_prior(); I <- rexp(50, 1 / 3); bh <- rnorm(50, 0, 1.5)
ps <- ns$.rv_prior_sd(pr, rep(0.01, 50)); base <- pr$sigma0sq * (2 * 0.01 * 0.99)^pr$alpha
r1 <- ns$.post_mix(bh, 1 / I, ps$sd, ps$w)
r2 <- ns$post_mix_mv_cpp(matrix(bh * I), matrix(I), rep(base, 50), c(pr$grid, 0), c(pr$grid * 0, pr$wide_sd^2),
                         c(pr$weights * (1 - pr$wide_weight), pr$wide_weight), 0, 0, 1, 1, 1L)
stopifnot(max(abs(r1$mean - r2$mean)) < 1e-10, max(abs(r1$sd - r2$sd)) < 1e-10)
cat("rare-variant checks passed\n")
