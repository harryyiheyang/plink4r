#' Polygenic scores from a PGEN, with plink2 or LD-anchored filling of missing calls
#'
#' Computes `S = D W` for every sample, with `D` the dosage of each weight's A1 allele and `W` one column
#' of weights per trait, as `plink2 --score` with `--score-col-nums` does. The function is a composition of the
#' genotype layer: the weight variants are opened a chunk at a time with [geno_open()] (hard calls, dosages and
#' missing calls alike; a dosage is held to 0.01 as everywhere in the package), their missing calls are filled as
#' by [genotype_fill()], and the chunk adds [geno_multiply()] of its columns to the scores. Memory is the n x q
#' score matrix plus one chunk of variants (or one LD block).
#'
#' Matching follows plink2: a weight's `A1` must equal the variant's ALT or REF (`2 - x` is used for REF);
#' no strand flipping. Variants whose ID is absent, whose `A1` matches neither allele, that are
#' multiallelic, all missing, or whose ID occurs more than once in the .pvar (`dup_ids = "skip"`, plink2
#' `ignore-dup-ids`) are listed in `skipped`.
#'
#' Missing calls (`impute`):
#' * `"anchor"` (default): the fill of [genotype_fill()]. With an A1 frequency `p` in `weights[[freq_col]]`
#'   the missing calls of a variant take the constant `c = (2 n p - sum_obs) / n_miss` (clipped to [0, 2]),
#'   so observed + filled allele count equals `2 n p`; without one, `c` is the observed mean, plink2's default
#'   mean imputation. With `ld`, a variant of the LD panel that has missing calls takes the LD fill of
#'   [genotype_fill()] instead (per-cell LD conditional mean, stage 2 anchored to `p` or to the panel's A1Freq):
#'   the panel SNPs of its LD block that are in the file are opened and condition the fill. Variants outside the
#'   panel, or without an anchor frequency inside (0, 1), keep the constant.
#' * `"none"`: plink2 `no-mean-imputation`; a missing call adds nothing and lowers `allele_ct` by 2.
#'
#' A fill is stored on the grid of the genotype store (1/16384 for a file of hard calls, 0.01 for a file with
#' dosages), so a filled cell can differ from plink2's by half a step.
#'
#' `transform = "center"` scores `d - 2p` and `"standardize"` `(d - 2p) / sqrt(2p(1 - p))` (plink2 `center`,
#' `variance-standardize`), with `p` the A1 frequency from `freq_col` or else the observed one; variants with
#' `p(1 - p) = 0` are skipped under `"standardize"`. All chromosomes are scored as diploid, so use autosomes
#' to compare with plink2.
#'
#' @param pfile PLINK2 prefix (.pgen/.pvar/.psam; .pvar may be .pvar.zst).
#' @param weights data.frame: variant ID, A1, one or more numeric weight columns, optionally `freq_col`.
#' @param id_col,a1_col Names of the ID and A1 columns of `weights`.
#' @param score_cols Weight columns; `NULL` takes every numeric column other than `freq_col`.
#' @param freq_col Column of `weights` with the A1 frequency used as anchor (plink2 `--read-freq`); ignored
#'   when absent; `NA` entries fall back to the observed frequency.
#' @param impute `"anchor"` or `"none"`.
#' @param ld `NULL`, or `list(dir = the LD panel folder, ref = optional table in place of its
#'   snpinfo.parquet)` as for [genotype_fill()]: the panel is aligned to the .pvar by chromosome, position and
#'   allele pair, in either orientation; only panel SNPs present in the .pvar enter the imputation.
#' @param transform `"none"`, `"center"` or `"standardize"`.
#' @param sample IIDs to score (`NULL` = all .psam rows), in the row order wanted.
#' @param dup_ids `"error"` or `"skip"`.
#' @param para named list overriding [geno_score_para()].
#' @param threads OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return list(sum = n x q scores (plink2 `<name>_SUM`), avg = sum / allele_ct (`<name>_AVG`),
#'   allele_ct (`ALLELE_CT`), sample, variants = data.table(ID, A1, flip (A1 is REF), fill ("mean",
#'   "anchor", "ld", "none"), n_miss, p (A1 frequency used)), skipped = data.table(ID, A1, reason),
#'   ld_blocks = data.table(block, m (panel SNPs opened), scored, cells (LD-filled cells of the scored
#'   variants), sec), simd (instruction set picked at run time for the SIMD kernels: avx512, avx2 or generic;
#'   environment variable PLINK4R_SIMD caps it), resources, time).
#' @export
geno_score <- function(pfile, weights, id_col = "ID", a1_col = "A1", score_cols = NULL, freq_col = "A1Freq",
                       impute = c("anchor", "none"), ld = NULL, transform = c("none", "center", "standardize"),
                       sample = NULL, dup_ids = c("error", "skip"), para = list(), threads = 4L) {
  t0 <- proc.time()[["elapsed"]]
  impute <- match.arg(impute); transform <- match.arg(transform); dup_ids <- match.arg(dup_ids)
  P <- geno_score_para()
  bad <- setdiff(names(para), names(P))
  if (length(bad)) stop("unknown para: ", paste(bad, collapse = ", "), call. = FALSE)
  P[names(para)] <- para
  threads <- .threads(threads)
  use_ld <- !is.null(ld) && impute == "anchor"
  w <- data.table::as.data.table(weights)
  if (!all(c(id_col, a1_col) %in% names(w))) stop("weights lacks ", id_col, " or ", a1_col, ".", call. = FALSE)
  if (is.null(score_cols))
    score_cols <- setdiff(names(w)[vapply(w, is.numeric, logical(1))], c(id_col, a1_col, freq_col))
  if (!length(score_cols)) stop("No numeric weight columns in weights.", call. = FALSE)
  ids <- as.character(w[[id_col]])
  wa1 <- toupper(as.character(w[[a1_col]]))
  if (anyDuplicated(ids)) stop("Variant IDs repeated in weights, e.g. ", ids[anyDuplicated(ids)], call. = FALSE)
  W <- as.matrix(w[, score_cols, with = FALSE]) + 0
  if (any(!is.finite(W))) stop("Weights must be finite.", call. = FALSE)
  q <- ncol(W)
  pfreq <- if (freq_col %in% names(w)) as.numeric(w[[freq_col]]) else rep(NA_real_, length(ids))
  pfreq[!(pfreq >= 0 & pfreq <= 1)] <- NA_real_

  # ---- variants: the weight's A1 is the file's counted allele (ALT) or the other one (flip)
  gf <- geno_file(pgenfile = pfile)
  v <- gf$variants
  srow <- if (is.null(sample)) NULL else
    .geno_index(as.character(sample), gf$samples$iid, paste(gf$samples$fid, gf$samples$iid, sep = "_"), "sample")
  n <- if (is.null(srow)) nrow(gf$samples) else length(srow)
  skip <- rep(NA_character_, length(ids))
  idx <- match(ids, v$id)
  skip[is.na(idx)] <- "not in pvar"
  dup <- !is.na(idx) & ids %in% v$id[duplicated(v$id)]
  if (any(dup) && dup_ids == "error")
    stop("Variant ID ", ids[which(dup)[1L]], " occurs more than once in the .pvar; use dup_ids = \"skip\" ",
         "(plink2 ignore-dup-ids).", call. = FALSE)
  skip[is.na(skip) & dup] <- "duplicate ID in pvar"
  skip[is.na(skip) & v$allele_ct[idx] > 2L] <- "multiallelic"
  flip <- wa1 == v$a2[idx]
  skip[is.na(skip) & !(wa1 == v$a1[idx] | flip)] <- "A1 matches neither allele"
  ok <- which(is.na(skip))
  if (!length(ok)) stop("No weight variant matched the .pvar.", call. = FALSE)
  af_given <- ifelse(flip, 1 - pfreq, pfreq)  # the given anchor, as a frequency of the counted allele

  # ---- LD candidates: weight variants of the panel with an anchor inside (0, 1), by LD block
  blk <- rep(NA_integer_, length(ids))
  if (use_ld) {
    al <- .ld_align(.ld_panel(ld$ref, ld$dir), v[, .(CHROM = chr, POS = pos, REF = a2, ALT = a1, row, allele_ct)])
    k <- match(idx[ok] - 1L, al$vi)
    anc <- ifelse(is.finite(af_given[ok]), af_given[ok], ifelse(al$sign[k] > 0, al$A1Freq[k], 1 - al$A1Freq[k]))
    cand <- !is.na(k) & is.finite(anc) & anc > 0 & anc < 1
    blk[ok[cand]] <- al$Block[k[cand]]
  }
  fill_para <- c(P[intersect(names(P), names(genotype_fill_para()))], list(threads = threads))

  # ---- groups of weight variants, each opened as one object: chunks in file order, then the candidates of each LD block
  ord <- ok[order(idx[ok])]
  in_ld <- !is.na(blk[ord])
  plain <- ord[!in_ld]
  groups <- c(unname(split(plain, (seq_along(plain) - 1L) %/% as.integer(P$chunk))), unname(split(ord[in_ld], blk[ord[in_ld]])))
  gblock <- c(rep(NA_integer_, length(groups) - length(unique(blk[ord[in_ld]]))), sort(unique(blk[ord[in_ld]])))
  S <- matrix(0, n, q); allele_ct <- numeric(n)
  n_miss <- p1v <- rep(NA_real_, length(ids)); fillv <- rep(NA_character_, length(ids)); status <- rep(NA_integer_, length(ids))
  ld_info <- list(); t_ld <- 0
  pb <- if (isTRUE(P$progress)) utils::txtProgressBar(0, length(groups), style = 3) else NULL
  for (gi in seq_along(groups)) {
    j <- groups[[gi]]
    X <- geno_open(file = gf, snp_vec = idx[j], sample_vec = srow, scale = FALSE, threads = threads)
    fl <- flip[j]; nm <- X$n_miss
    p1 <- ifelse(is.finite(pfreq[j]), pfreq[j], ifelse(fl, 1 - X$af, X$af))
    st <- ifelse(nm == n, 1L, ifelse(transform == "standardize" & !(p1 > 0 & p1 < 1), 2L, 0L))
    sc <- st == 0L
    off <- if (transform == "none") 0 else 2 * p1
    a <- if (transform == "standardize") 1 / sqrt(2 * p1 * (1 - p1)) else 1
    # a cell holding x (count of the file's counted allele) adds a w (cc + s x - off) = B x + C
    B <- W[j, , drop = FALSE] * (ifelse(fl, -1, 1) * a)
    C <- W[j, , drop = FALSE] * (a * (ifelse(fl, 2, 0) - off))
    B[!sc, ] <- 0; C[!sc, ] <- 0
    ct <- 2 * sum(sc)
    lab <- ifelse(is.finite(pfreq[j]), "anchor", "mean")
    if (impute == "none") {  # a missing cell holds 0, and its C is taken back
      X <- .geno_refill(X, const = rep(0, X$p), inplace = TRUE, impute = "custom", threads = threads)
      Mm <- geno_miss_mx_cpp(X$ptr, cbind(C, sc + 0), threads)
      Sg <- geno_multiply(X, B, threads) - Mm[, seq_len(q), drop = FALSE]
      ct <- ct - 2 * Mm[, q + 1L]
      lab[] <- "none"
    } else if (!is.na(gblock[gi]) && any(sc & nm > 0)) {  # the LD block of these variants conditions their missing cells
      tb <- proc.time()[["elapsed"]]
      wr <- .ld_rows(al, gblock[gi])
      Xb <- geno_open(file = gf, snp_vec = wr$vi + 1L, sample_vec = srow, scale = FALSE, threads = threads)
      cj <- match(idx[j] - 1L, wr$vi)
      afb <- rep(NA_real_, Xb$p); afb[cj] <- af_given[j]
      Xb <- genotype_fill(Xb, af = afb, ld = list(dir = ld$dir, ref = wr[, intersect(.panel_cols, names(wr)), with = FALSE]),
                          para = fill_para)
      Bb <- matrix(0, Xb$p, q); Bb[cj, ] <- B
      Sg <- geno_multiply(Xb, Bb, threads)
      lab[sc & nm > 0] <- "ld"
      dt <- proc.time()[["elapsed"]] - tb; t_ld <- t_ld + dt
      ld_info[[length(ld_info) + 1L]] <- data.table::data.table(block = gblock[gi], m = Xb$p, scored = sum(sc),
                                                                cells = sum(nm[sc]), sec = dt)
      rm(Xb)
    } else {
      X <- .geno_refill(X, af = af_given[j], inplace = TRUE, impute = "af", threads = threads)
      Sg <- geno_multiply(X, B, threads)
    }
    S <- S + Sg + rep(colSums(C), each = n)
    allele_ct <- allele_ct + ct
    n_miss[j] <- nm; p1v[j] <- p1; fillv[j] <- lab; status[j] <- st
    rm(X)
    if (!is.null(pb)) utils::setTxtProgressBar(pb, gi)
  }
  if (!is.null(pb)) close(pb)
  skip[ok[status[ok] == 1L]] <- "all calls missing"
  skip[ok[status[ok] == 2L]] <- "zero variance (standardize)"
  ok <- ok[status[ok] == 0L]
  if (!length(ok)) stop("No variant left to score.", call. = FALSE)
  rn <- if (is.null(srow)) gf$samples$iid else gf$samples$iid[srow]
  dimnames(S) <- list(rn, score_cols)
  list(sum = S, avg = S / allele_ct, allele_ct = allele_ct, sample = rn,
       variants = data.table::data.table(ID = ids[ok], A1 = wa1[ok], flip = flip[ok], fill = fillv[ok],
                                         n_miss = n_miss[ok], p = p1v[ok]),
       skipped = data.table::data.table(ID = ids[!is.na(skip)], A1 = wa1[!is.na(skip)], reason = skip[!is.na(skip)]),
       ld_blocks = if (length(ld_info)) data.table::rbindlist(ld_info),
       simd = simd_isa_cpp(), resources = c(threads = threads, mem_gb = .mem_gb(P$mem_gb, threads)),
       time = c(ld = t_ld, total = proc.time()[["elapsed"]] - t0))
}

#' Defaults of [geno_score()]
#'
#' chunk 4096 weight variants opened at a time (memory: about n / 4 bytes per variant of hard calls, more with
#' dosages); progress: a bar over the chunks and LD blocks; and, for the LD fill, the values of
#' [genotype_fill_para()] with the same meaning and defaults: delta, anchor_min_miss, flank, thresh, mem_gb,
#' dosage_step.
#' @export
geno_score_para <- function() {
  c(list(chunk = 4096L, progress = interactive()),
    genotype_fill_para()[c("delta", "anchor_min_miss", "flank", "thresh", "mem_gb", "dosage_step")])
}
