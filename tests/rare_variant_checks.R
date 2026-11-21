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
# 4. below robust_min_mac the model-based P is used
pm <- c(0.3, 0.02, 0.5)
stopifnot(identical(f(c(4, 9, 1), N, c(2, 4, 50) / (2 * N), TRUE, pm, 5)[1:2], pm[1:2]),
          f(1, N, 50 / (2 * N), TRUE, pm[3], 5) != pm[3])
# 5. default reporting floors on a tiny BED: GWAS P only for MAC > 20, interaction P only for MAF >= min_maf_int
td <- tempfile("p4r"); dir.create(td); pre <- file.path(td, "t")
n <- 2000; macs <- c(10, 20, 21, 40)
G <- sapply(macs, function(m) { g <- integer(n); g[seq_len(m) * 7L] <- 1L; g })
code <- matrix(c(3L, 2L, 0L)[G + 1L], n)                        # BED: 11 = 0 copies of A1, 10 = 1, 00 = 2
cm <- array(code, c(4L, n / 4L, length(macs)))
con <- file(paste0(pre, ".bed"), "wb"); writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
writeBin(as.raw(as.vector(cm[1, , ] + 4L * cm[2, , ] + 16L * cm[3, , ] + 64L * cm[4, , ])), con); close(con)
writeLines(sprintf("1\tv%d\t0\t%d\tA\tG", seq_along(macs), seq_along(macs)), paste0(pre, ".bim"))
writeLines(sprintf("i%d\ti%d\t0\t0\t0\t-9", 1:n, 1:n), paste0(pre, ".fam"))
set.seed(2); d <- data.frame(IID = paste0("i", 1:n), c1 = rnorm(n), e = rnorm(n), y = rnorm(n))
a <- gwas(y ~ c1, data = d, bfile = pre, robust = TRUE, threads = 1)
stopifnot(identical(is.na(a$P), macs <= 20), identical(is.na(a$P_ROBUST), macs <= 20), !anyNA(a$BETA))
b <- gxe(y ~ e, data = d, exposure = "e", bfile = pre, threads = 1, min_maf_int = 30 / (2 * n))
stopifnot(identical(is.na(b$P_INT_ROBUST), macs < 30), identical(is.na(b$P_MARGINAL_ROBUST), macs <= 20))
b <- gxe(y ~ e, data = d, exposure = "e", bfile = pre, threads = 1)   # default MAF 0.001 = MAC 4 here
stopifnot(!anyNA(b$P_INT_ROBUST))
cat("rare-variant checks passed\n")
