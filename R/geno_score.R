#' Polygenic scores from a PGEN, with plink2 or LD-anchored filling of missing calls
#'
#' Computes `S = D W` for every sample, with `D` the dosage of each weight's A1 allele and `W` one column
#' of weights per trait, as `plink2 --score` with `--score-col-nums` does. The genotype matrix is never
#' formed: variants are streamed through a per-variant lookup table (code 0/1/2/missing -> value added to
#' each trait), in parallel over sample ranges, so memory is the n x q score matrix plus one chunk of
#' variants.
#'
#' Matching follows plink2: a weight's `A1` must equal the variant's ALT or REF (`2 - x` is used for REF);
#' no strand flipping. Variants whose ID is absent, whose `A1` matches neither allele, that are
#' multiallelic, all missing, or whose ID occurs more than once in the .pvar (`dup_ids = "skip"`, plink2
#' `ignore-dup-ids`) are listed in `skipped`.
#'
#' Missing calls (`impute`):
#' * `"anchor"` (default): scheme 3 of [genotype_fill()]. With an A1 frequency `p` in `weights[[freq_col]]`
#'   the missing calls of a variant take the constant `c = (2 n p - sum_obs) / n_miss` (clipped to [0, 2]),
#'   so observed + filled allele count equals `2 n p`; without one, `c` is the observed mean, which is
#'   exactly plink2's default mean imputation. With `ld`, variants in the eigen.bin panel are filled by
#'   scheme 4 instead (per-cell LD conditional mean, stage 2 anchored to `p` or to `AF_PanUKB`), one LD
#'   block at a time through the same kernel as [genotype_impute()]; variants outside the panel keep
#'   scheme 3.
#' * `"none"`: plink2 `no-mean-imputation`; a missing call adds nothing and lowers `allele_ct` by 2.
#'
#' `transform = "center"` scores `d - 2p` and `"standardize"` `(d - 2p) / sqrt(2p(1 - p))` (plink2 `center`,
#' `variance-standardize`), with `p` the A1 frequency from `freq_col` or else the observed one; variants with
#' `p(1 - p) = 0` are skipped under `"standardize"`. All chromosomes are scored as diploid, so use autosomes
#' to compare with plink2. Hard calls only: a PGEN with dosages is refused.
#'
#' @param pfile PLINK2 prefix (.pgen/.pvar/.psam; .pvar may be .pvar.zst).
#' @param weights data.frame: variant ID, A1, one or more numeric weight columns, optionally `freq_col`.
#' @param id_col,a1_col Names of the ID and A1 columns of `weights`.
#' @param score_cols Weight columns; `NULL` takes every numeric column other than `freq_col`.
#' @param freq_col Column of `weights` with the A1 frequency used as anchor (plink2 `--read-freq`); ignored
#'   when absent; `NA` entries fall back to the observed frequency.
#' @param impute `"anchor"` or `"none"`.
#' @param ld `NULL`, or `list(dir = folder of block<N>.eigen.bin and the panel SNP .parquet, ref = optional
#'   table in place of that parquet)` as for
#'   [genotype_impute()]: the panel is aligned to the .pvar by ID and alleles, else by chromosome, position and
#'   allele pair, in either orientation; only panel SNPs present in the .pvar enter the imputation.
#' @param transform `"none"`, `"center"` or `"standardize"`.
#' @param sample IIDs to score (`NULL` = all .psam rows), in the row order wanted.
#' @param dup_ids `"error"` or `"skip"`.
#' @param para named list overriding [geno_score_para()].
#' @param threads OpenMP threads; 0 takes them from [plink4r_resources()].
#' @return list(sum = n x q scores (plink2 `<name>_SUM`), avg = sum / allele_ct (`<name>_AVG`),
#'   allele_ct (`ALLELE_CT`), sample, variants = data.table(ID, A1, flip (A1 is REF), fill ("mean",
#'   "anchor", "ld", "none"), n_miss, p (A1 frequency used)), skipped = data.table(ID, A1, reason),
#'   ld_blocks, simd (instruction set picked at run time for the LD kernel: avx512, avx2 or scalar;
#'   environment variable PLINK4R_SIMD caps it), time).
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
  threads <- .threads(threads); P$mem_gb <- .mem_gb(P$mem_gb, threads)
  if (!is.null(ld)) .check_delta(P$delta)
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
  pfreq <- if (freq_col %in% names(w)) as.numeric(w[[freq_col]]) else rep(NA_real_, length(ids))
  pfreq[!(pfreq >= 0 & pfreq <= 1)] <- NA_real_

  # ---- variants
  pvar <- .read_pvar5(pfile)
  act <- pvar$allele_ct
  psam <- data.table::fread(paste0(pfile, ".psam"), colClasses = "character")
  iid <- psam[[intersect(c("#IID", "IID"), names(psam))[1L]]]
  si <- if (is.null(sample)) seq_along(iid) else match(as.character(sample), iid)
  if (anyNA(si)) stop(sum(is.na(si)), " samples not in .psam", call. = FALSE)
  n <- length(si)
  pgen <- normalizePath(paste0(pfile, ".pgen"))

  skip <- rep(NA_character_, length(ids))
  idx <- match(ids, pvar$ID)
  skip[is.na(idx)] <- "not in pvar"
  dup <- !is.na(idx) & ids %in% pvar$ID[duplicated(pvar$ID)]
  if (any(dup) && dup_ids == "error")
    stop("Variant ID ", ids[which(dup)[1L]], " occurs more than once in the .pvar; use dup_ids = \"skip\" ",
         "(plink2 ignore-dup-ids).", call. = FALSE)
  skip[is.na(skip) & dup] <- "duplicate ID in pvar"
  skip[is.na(skip) & act[idx] > 2L] <- "multiallelic"
  flip <- wa1 == toupper(pvar$REF[idx])
  skip[is.na(skip) & !(wa1 == toupper(pvar$ALT[idx]) | flip)] <- "A1 matches neither allele"
  ok <- which(is.na(skip))
  if (!length(ok)) stop("No weight variant matched the .pvar.", call. = FALSE)

  st <- pgen_stats_cpp(pgen, length(iid), act, as.integer(idx[ok] - 1L), as.integer(si - 1L), threads)
  if (isTRUE(st$dosage)) stop("geno_score scores hard calls; this PGEN holds dosages.", call. = FALSE)
  n_obs <- st$n_obs; sum_obs <- st$sum; n_miss <- n - n_obs
  allmiss <- n_obs == 0
  skip[ok[allmiss]] <- "all calls missing"
  # A1 frequency: given, else observed.
  p_obs_alt <- sum_obs / (2 * pmax(n_obs, 1))
  p1 <- ifelse(is.finite(pfreq[ok]), pfreq[ok], ifelse(flip[ok], 1 - p_obs_alt, p_obs_alt))
  if (transform == "standardize") skip[ok[is.na(skip[ok]) & !(p1 > 0 & p1 < 1)]] <- "zero variance (standardize)"
  keep <- is.na(skip[ok])
  ok <- ok[keep]; n_obs <- n_obs[keep]; sum_obs <- sum_obs[keep]; n_miss <- n_miss[keep]; p1 <- p1[keep]
  m <- length(ok)
  if (!m) stop("No variant left to score.", call. = FALSE)
  fl <- flip[ok]
  q <- ncol(W)

  # ---- fill of a missing call, in ALT units, and which variants go through the LD panel
  fill <- rep(if (impute == "none") "none" else "mean", m)
  p_alt_anchor <- ifelse(fl, 1 - pfreq[ok], pfreq[ok])
  f_alt <- sum_obs / pmax(n_obs, 1)
  if (impute == "anchor") {
    a <- is.finite(p_alt_anchor) & n_miss > 0
    f_alt[a] <- pmin(2, pmax(0, (2 * n * p_alt_anchor[a] - sum_obs[a]) / n_miss[a]))
    fill[is.finite(p_alt_anchor)] <- "anchor"
  }
  in_ld <- rep(FALSE, m)
  if (!is.null(ld) && impute == "anchor") {
    al <- .ld_align(.ld_panel(ld$ref, ld$dir), pvar)  # panel rows in the PGEN, with vi and sign
    al[, `:=`(arow = .I, anchor = ifelse(sign > 0, AF_PanUKB, 1 - AF_PanUKB))]  # anchor: PGEN-ALT frequency
    k <- match(idx[ok] - 1L, al$vi)
    anc <- ifelse(is.finite(p_alt_anchor), p_alt_anchor, al$anchor[k])
    in_ld <- !is.na(k) & is.finite(anc) & anc > 0 & anc < 1 & n_miss > 0
    fill[in_ld] <- "ld"
  }

  # ---- per-variant tables: value a (d - off) w for A1 dosage d = c0 + s x, x = ALT count of code 0/1/2
  s <- ifelse(fl, -1, 1); c0 <- 2 * fl
  a <- rep(1, m); off <- rep(0, m)
  if (transform != "none") {
    off <- 2 * p1
    if (transform == "standardize") a <- 1 / sqrt(2 * p1 * (1 - p1))
  }
  Wa <- W[ok, , drop = FALSE] * a
  tab <- array(0, c(q, 4L, m))
  for (code in 0:2) tab[, code + 1L, ] <- t(Wa * (c0 + s * code - off))
  d_fill <- c0 + s * f_alt
  tab[, 4L, ] <- t(Wa * ifelse(fill == "none", 0, ifelse(in_ld, c0 - off, d_fill - off)))
  ct <- matrix(2L, 4L, m)
  ct[4L, fill == "none"] <- 0L
  gemm <- if (q >= 2L) 2L else 0L  # one trait: lookup table; several: float32 product per chunk
  # One progress bar: a step per chromosome of the stream (scheme 3 / plink / none), then per LD block.
  chr <- pvar$CHROM[idx[ok]]
  chrs <- unique(chr)
  ld_blocks_todo <- if (any(in_ld)) sort(unique(al$block[k[in_ld]])) else integer(0)
  pb <- if (isTRUE(P$progress)) utils::txtProgressBar(0, length(chrs) + length(ld_blocks_todo), style = 3) else NULL
  step <- 0L
  tick <- function() { step <<- step + 1L; if (!is.null(pb)) utils::setTxtProgressBar(pb, step) }
  S <- matrix(0, n, q); allele_ct <- integer(n)
  for (cc in chrs) {
    sel <- which(chr == cc)
    res <- score_stream_cpp(pgen, length(iid), act, as.integer(idx[ok][sel] - 1L), as.integer(si - 1L),
                            as.numeric(tab[, , sel, drop = FALSE]), as.integer(ct[, sel, drop = FALSE]), q, threads,
                            as.integer(P$chunk), gemm)
    S <- S + res$score; allele_ct <- allele_ct + res$allele_ct
    tick()
  }
  t_stream <- proc.time()[["elapsed"]] - t0

  # ---- LD panel variants: S += (s a w) g over their filled cells, one block at a time
  ld_info <- NULL
  if (any(in_ld)) {
    kk <- k[in_ld]
    ua <- is.finite(p_alt_anchor[in_ld])
    al$anchor[kk[ua]] <- p_alt_anchor[in_ld][ua]
    Bw <- matrix(0, nrow(al), q)
    Bw[kk, ] <- Wa[in_ld, , drop = FALSE] * s[in_ld]
    blocks <- ld_blocks_todo
    ld_info <- vector("list", length(blocks))
    for (b in seq_along(blocks)) {
      tb <- proc.time()[["elapsed"]]
      wr <- .ld_rows(al, blocks[b])
      wr <- wr[anchor > 0 & anchor < 1]
      ef <- file.path(ld$dir, paste0("block", blocks[b], ".eigen.bin"))
      Bb <- Bw[wr$arow, , drop = FALSE]
      out <- score_ld_block_cpp(pgen, length(iid), act, as.integer(wr$vi), as.integer(si - 1L), normalizePath(ef),
                                wr$anchor, as.integer(wr$sign), P$thresh, P$delta, P$anchor_min_miss, Bb, threads,
                                if (is.null(P$nb_people)) 0L else as.integer(P$nb_people), P$mem_gb,
                                as.integer(wr$ld_index))
      S <- S + out$score
      tick()
      ld_info[[b]] <- data.table::data.table(block = blocks[b], m = nrow(wr), scored = sum(rowSums(Bb != 0) > 0),
                                             cells = out$cells, persons_solved = out$persons_solved,
                                             sec = proc.time()[["elapsed"]] - tb)
    }
    ld_info <- data.table::rbindlist(ld_info)
  }
  if (!is.null(pb)) close(pb)
  dimnames(S) <- list(iid[si], score_cols)
  list(sum = S, avg = S / allele_ct, allele_ct = allele_ct, sample = iid[si],
       variants = data.table::data.table(ID = ids[ok], A1 = wa1[ok], flip = fl, fill = fill, n_miss = n_miss, p = p1),
       skipped = data.table::data.table(ID = ids[!is.na(skip)], A1 = wa1[!is.na(skip)], reason = skip[!is.na(skip)]),
       ld_blocks = ld_info, simd = simd_isa_cpp(), resources = c(threads = threads, mem_gb = P$mem_gb),
       time = c(stream = t_stream, total = proc.time()[["elapsed"]] - t0))
}

#' Defaults of [geno_score()]
#'
#' chunk 512 variants read per streaming step; progress: one bar over chromosomes then LD blocks; for the LD path (as [genotype_impute()]): thresh 0 = all
#' stored eigenpairs, delta 0.01 ridge of stage 1, anchor_min_miss 0.01 (stage 2 at or above this missing
#' rate), nb_people NULL (people per batch from mem_gb), mem_gb 2.
#' @export
geno_score_para <- function() {
  list(chunk = 512L, progress = interactive(), thresh = 0, delta = 0.01, anchor_min_miss = 0.01, nb_people = NULL, mem_gb = 2)
}
