#' Polygenic scores from a BED or PGEN file (plink2 `--score` in R)
#'
#' Computes `S = D W` for every sample, with `D` the dosage of each weight's
#' A1 allele and `W` one column of weights per trait, as `plink2 --score` with
#' `--score-col-nums` does. The genotypes are read once into the compact
#' [geno_open()] store (missing calls filled in C++), and all traits are scored
#' in one [geno_multiply()] pass.
#'
#' Allele matching follows plink2: a weight's `A1` must equal the file's
#' counted allele (bim column 5 / pvar ALT) or the other allele (bim column 6 /
#' REF); when it is the other allele the dosage is `2 - x`. There is no strand
#' flipping. Variants whose ID is absent, whose `A1` matches neither allele,
#' that are multiallelic, or whose ID is repeated in the file (`dup_ids =
#' "skip"`) are dropped and listed in `skipped`.
#'
#' Missing genotypes (`impute`):
#' * `"mean"` (default, plink2's default): a missing call is `2 p`, `p` the A1
#'   frequency from `freq` (`--read-freq`) or, without `freq`, observed in the
#'   selected samples (stored to 1/32768).
#' * `"median"`: the variant's lower median of the observed calls.
#' * `"none"` (`no-mean-imputation`): a missing call contributes nothing and
#'   lowers the sample's `allele_ct` by 2. Reads the file twice.
#'
#' `transform = "center"` scores `d - 2p` and `"standardize"` scores
#' `(d - 2p) / sqrt(2p(1 - p))` (plink2 `center` / `variance-standardize`),
#' with `p` from `freq` or the observed A1 frequency; variants with
#' `p(1 - p) = 0` are dropped under `"standardize"`. All chromosomes are
#' scored as diploid (0..2), so restrict to autosomes to match plink2.
#'
#' @param weights data.frame with a variant ID column, an A1 column and one or
#'   more numeric weight columns (one per trait); one row per variant.
#' @param bedfile,pgenfile,sample_vec,threads As in [geno_open()].
#' @param id_col,a1_col Names of the ID and A1 columns of `weights`.
#' @param score_cols Weight columns; `NULL` takes every other numeric column.
#' @param impute `"mean"`, `"median"` or `"none"` (see Details).
#' @param freq Optional data.frame of A1 frequencies with columns `ID`, `A1`,
#'   `FREQ` (`A1` may be either allele; e.g. a plink2 `.afreq` with `ALT`,
#'   `ALT_FREQS` renamed). Variants without a frequency fall back to the
#'   observed one.
#' @param transform `"none"`, `"center"` or `"standardize"`.
#' @param dup_ids `"error"` (default) or `"skip"` (plink2 `ignore-dup-ids`)
#'   for weight IDs that occur more than once in the file.
#' @return A list with `sum` (n x q, plink2 `<name>_SUM`), `avg` (`sum /
#'   allele_ct`, plink2 `<name>_AVG`), `allele_ct` (plink2 `ALLELE_CT`),
#'   `dosage_sum` (A1 dosage sum, plink2 `NAMED_ALLELE_DOSAGE_SUM`), `sample`,
#'   `variants` (scored variants: ID, A1, file index, `flip` = A1 is the
#'   file's other allele, `p` = A1 frequency used) and `skipped` (ID, A1,
#'   reason).
#' @export
geno_score <- function(weights, bedfile = NULL, pgenfile = NULL, id_col = "ID", a1_col = "A1",
                       score_cols = NULL, impute = c("mean", "median", "none"), freq = NULL,
                       transform = c("none", "center", "standardize"), sample_vec = NULL,
                       dup_ids = c("error", "skip"), threads = 4L) {
  impute <- match.arg(impute)
  transform <- match.arg(transform)
  dup_ids <- match.arg(dup_ids)
  if (is.null(bedfile) == is.null(pgenfile)) stop("Give exactly one of bedfile and pgenfile.", call. = FALSE)
  w <- data.table::as.data.table(weights)
  if (!all(c(id_col, a1_col) %in% names(w))) stop("weights lacks ", id_col, " or ", a1_col, ".", call. = FALSE)
  if (is.null(score_cols)) {
    score_cols <- setdiff(names(w)[vapply(w, is.numeric, logical(1))], c(id_col, a1_col))
  }
  if (!length(score_cols)) stop("No numeric weight columns in weights.", call. = FALSE)
  ids <- as.character(w[[id_col]])
  wa1 <- toupper(as.character(w[[a1_col]]))
  if (anyDuplicated(ids)) stop("Variant IDs repeated in weights, e.g. ", ids[anyDuplicated(ids)], call. = FALSE)
  W <- as.matrix(w[, score_cols, with = FALSE]) + 0
  if (any(!is.finite(W))) stop("Weights must be finite.", call. = FALSE)

  vi <- .geno_score_variants(bedfile, pgenfile)
  skip <- rep(NA_character_, length(ids))
  idx <- match(ids, vi$snp)
  skip[is.na(idx)] <- "not in file"
  dup <- !is.na(idx) & ids %in% vi$snp[duplicated(vi$snp)]
  if (any(dup) && dup_ids == "error") {
    stop("Variant ID ", ids[which(dup)[1L]], " occurs more than once in the file; ",
         "use dup_ids = \"skip\" (plink2 ignore-dup-ids).", call. = FALSE)
  }
  skip[is.na(skip) & dup] <- "duplicate ID in file"
  multi <- is.na(skip) & grepl(",", vi$a1[idx], fixed = TRUE)
  skip[multi] <- "multiallelic"
  ca1 <- toupper(vi$a1[idx]); ca2 <- toupper(vi$a2[idx])
  flip <- wa1 == ca2
  skip[is.na(skip) & !(wa1 == ca1 | flip)] <- "A1 matches neither allele"

  # A1 frequency used for imputation / transforms, on the weight's A1.
  p <- rep(NA_real_, length(ids))
  if (!is.null(freq)) {
    fq <- data.table::as.data.table(freq)
    if (!all(c("ID", "A1", "FREQ") %in% names(fq))) stop("freq needs columns ID, A1, FREQ.", call. = FALSE)
    k <- match(ids, as.character(fq$ID))
    fa <- toupper(as.character(fq$A1))[k]
    fp <- as.numeric(fq$FREQ)[k]
    p <- ifelse(fa == wa1, fp, ifelse(!is.na(k), 1 - fp, NA_real_))
    p[!is.na(k) & !(fa %in% c(ca1, ca2))] <- NA_real_
  }
  ok <- is.na(skip)
  if (!any(ok)) stop("No weight variant matched the file.", call. = FALSE)
  j <- idx[ok]

  m <- length(j)
  open <- function(const, imp) {
    geno_open(bedfile = bedfile, pgenfile = pgenfile, snp_vec = j, sample_vec = sample_vec,
              impute = imp, scale = FALSE, threads = threads,
              fill = if (is.null(const)) NULL else list(const = const))
  }
  # Counted-allele fill: 2 p_counted where a frequency is given (NaN = built-in mean).
  pc <- ifelse(flip[ok], 1 - p[ok], p[ok])
  X0 <- X1 <- NULL
  if (impute == "mean") {
    X <- open(if (any(is.finite(pc))) ifelse(is.finite(pc), 2 * pc, NaN) else NULL, "mean")
    p_obs <- X$center / 2  # observed-mean fill keeps the observed mean
  } else {
    if (impute == "none" || (transform != "none" && !all(is.finite(pc)))) {
      # Missing calls filled with 0 and with 1: their difference is the missing indicator.
      X0 <- open(rep(0, m), "mean")
      X1 <- open(rep(1, m), "mean")
      nmiss_v <- (X1$center - X0$center) * X0$n
      p_obs <- X0$center * X0$n / pmax(X0$n - nmiss_v, 1) / 2
    } else {
      p_obs <- rep(NA_real_, m)
    }
    X <- if (impute == "none") X0 else open(NULL, "median")
  }
  n <- X$n
  pa1 <- ifelse(is.finite(p[ok]), p[ok], ifelse(flip[ok], 1 - p_obs, p_obs))

  # A1 dosage d = c + s x (c = 2 for a flipped variant, s = -1); scored value a (d - off).
  s <- ifelse(flip[ok], -1, 1); c0 <- 2 * flip[ok]
  a <- rep(1, m); off <- rep(0, m); keep <- rep(TRUE, m)
  if (transform != "none") {
    off <- 2 * pa1
    if (transform == "standardize") {
      v <- 2 * pa1 * (1 - pa1)
      keep <- is.finite(v) & v > 0
      a <- ifelse(keep, 1 / sqrt(pmax(v, 1e-300)), 0)
      skip[which(ok)[!keep]] <- "zero variance (standardize)"
    }
  }
  Wa <- W[ok, , drop = FALSE] * a
  q <- ncol(Wa)
  B <- cbind(Wa * s, s * keep)                       # coefficients of x
  C <- cbind(Wa * (c0 - off), c0 * keep)             # constants per call
  if (impute == "none") {
    # Observed calls only: X0 B + (1 - M) C with M = X1 - X0.
    one <- matrix(keep + 0, m, 1)
    MC <- geno_multiply(X1, cbind(C, one), threads = threads) - geno_multiply(X0, cbind(C, one), threads = threads)
    R <- sweep(geno_multiply(X0, B, threads = threads), 2L, colSums(C), "+") - MC[, seq_len(q + 1L), drop = FALSE]
    allele_ct <- 2 * (sum(keep) - round(MC[, q + 2L]))
  } else {
    R <- sweep(geno_multiply(X, B, threads = threads), 2L, colSums(C), "+")
    allele_ct <- rep(2 * sum(keep), n)
  }
  S <- R[, seq_len(q), drop = FALSE]
  dimnames(S) <- list(X$sample, score_cols)
  used <- which(ok)[keep]
  list(sum = S, avg = S / allele_ct, allele_ct = allele_ct, dosage_sum = as.numeric(R[, q + 1L]),
       sample = X$sample,
       variants = data.table::data.table(ID = ids[used], A1 = wa1[used], idx = idx[used],
                                         flip = flip[used], p = pa1[keep]),
       skipped = data.table::data.table(ID = ids[!is.na(skip)], A1 = wa1[!is.na(skip)],
                                        reason = skip[!is.na(skip)]),
       impute = impute, transform = transform)
}

# Variant IDs and alleles of a BED or PGEN, without opening the genotypes.
.geno_score_variants <- function(bedfile, pgenfile) {
  if (!is.null(bedfile)) {
    base <- sub("\\.bed$", "", .geno_bed_path(bedfile))
    bim <- data.table::fread(paste0(base, ".bim"), header = FALSE, colClasses = "character")
    return(list(snp = bim[[2L]], a1 = bim[[5L]], a2 = bim[[6L]]))
  }
  pv <- data.table::fread(.geno_pgen_files(pgenfile)$pvar, skip = "#CHROM", colClasses = "character")
  list(snp = pv[["ID"]], a1 = pv[["ALT"]], a2 = pv[["REF"]])
}
