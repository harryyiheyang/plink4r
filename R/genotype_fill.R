#' Locate SNPs in a PGEN by CHR / BP / A1 / A2 and compute the fills of their missing genotypes
#'
#' One call gives everything a fine-mapping genotype reader needs to fill missing cells
#' (SuSiE4I calls it for `X = list(pfile, snp, ld, para)`):
#' 1. SNPs of `snp` are matched to the .pvar by chromosome, position and the allele pair
#'    ({A1, A2} = {ALT, REF}, either orientation; values are ALT dosages).
#' 2. Per-SNP observed count and ALT-dosage sum over `sample` (dosage-aware) -> missing rate and observed AF;
#'    SNPs with missing rate > `miss_max` or observed MAF < `maf` are dropped.
#' 3. Fill of a missing cell, by what is given:
#'    - default: observed mean dosage;
#'    - `A1Freq` in `snp` (scheme 3): c_j = (2 n p_j - sum_obs g) / m_j clipped to [0, 2], p_j the ALT frequency;
#'    - `ld` (scheme 4): per-cell conditional mean from the eigen.bin panel ([genotype_impute()]: stage 1
#'      z_M = R_MO (R_OO + delta I)^-1 z_O on the SNPs +- `flank` panel rows, stage 2 anchor shift when the
#'      missing rate >= `anchor_min_miss`), anchored to `A1Freq` when given, else to the panel AF.
#'      The panel is aligned to the PGEN by ID and alleles, else by chromosome, position and allele pair, in
#'      either orientation (a swapped REF/ALT enters with a negative sign); SNPs absent from the panel keep the
#'      per-SNP constant.
#'
#' @param pfile PLINK2 prefix (.pgen/.pvar/.psam); hard calls, dosages, or both, with missing cells.
#' @param snp data.frame with CHR, BP37 or BP38 (the build of the .pvar POS), A1, A2; optional A1Freq.
#' @param ld NULL, or list(dir = folder of block<N>.eigen.bin and the panel SNP .parquet, ref = optional
#'   parquet path or data.frame in place of that table; see [genotype_impute()]).
#' @param sample IIDs (NULL = all .psam rows), in the row order wanted.
#' @param para named list overriding the defaults of [genotype_fill_para()].
#' @return list(snp = one row per kept SNP: ID (pvar), idx (1-based pvar row), CHR, POS, REF, ALT, A1, A2,
#'   a1_is_alt, n_miss, miss_rate, af_alt_obs, fill_const, in_ld; cells = data.table(ID, IID, g) per-cell fills
#'   (scheme 4 only, ALT dosage); dropped = SNPs not used with the reason; method; para; time)
#' @export
genotype_fill <- function(pfile, snp, ld = NULL, sample = NULL, para = list()) {
  t0 <- proc.time()[["elapsed"]]
  P <- genotype_fill_para()
  bad <- setdiff(names(para), names(P))
  if (length(bad)) stop("unknown para: ", paste(bad, collapse = ", "))
  P[names(para)] <- para
  P$threads <- .threads(P$threads); P$mem_gb <- .mem_gb(P$mem_gb, P$threads)
  if (!is.null(ld)) .check_delta(P$delta)
  snp <- data.table::copy(data.table::as.data.table(snp))
  bpcol <- if (!is.null(P$build)) paste0("BP", P$build) else intersect(c("BP37", "BP38"), names(snp))
  if (length(bpcol) != 1L) stop("`snp` needs exactly one of BP37 / BP38, or set para$build")
  miss_col <- setdiff(c("CHR", bpcol, "A1", "A2"), names(snp))
  if (length(miss_col)) stop("`snp` lacks ", paste(miss_col, collapse = ", "))
  chr_std <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

  pvar <- .read_pvar5(pfile)
  act <- pvar$allele_ct

  q <- data.table::data.table(k = seq_len(nrow(snp)), CHROM = chr_std(snp$CHR), POS = as.integer(snp[[bpcol]]),
                              A1 = toupper(snp$A1), A2 = toupper(snp$A2),
                              A1Freq = if ("A1Freq" %in% names(snp)) as.numeric(snp$A1Freq) else NA_real_)
  m <- merge(q, pvar[, .(CHROM, POS, ID, REF, ALT, row)], by = c("CHROM", "POS"), allow.cartesian = TRUE)
  m <- m[(A1 == ALT & A2 == REF) | (A1 == REF & A2 == ALT)]
  m[, a1_is_alt := A1 == ALT]
  m <- m[!duplicated(k)]
  dropped <- data.table::data.table(k = setdiff(q$k, m$k), reason = "not in pvar (CHR/BP/alleles)")
  data.table::setorder(m, k)

  psam <- data.table::fread(paste0(pfile, ".psam"), colClasses = "character")
  iid <- psam[[intersect(c("IID", "#IID"), names(psam))[1]]]
  si <- if (is.null(sample)) seq_along(iid) else match(as.character(sample), iid)
  if (anyNA(si)) stop(sum(is.na(si)), " samples not in .psam")
  n <- length(si)
  pg <- normalizePath(paste0(pfile, ".pgen"))
  st <- pgen_stats_cpp(pg, length(iid), act, as.integer(m$row - 1L), as.integer(si - 1L), as.integer(P$threads))
  m[, `:=`(n_miss = n - st$n_obs, miss_rate = 1 - st$n_obs / n, af_alt_obs = st$sum / (2 * pmax(st$n_obs, 1)), sum_obs = st$sum)]
  r_miss <- m$miss_rate > P$miss_max
  r_maf <- !r_miss & pmin(m$af_alt_obs, 1 - m$af_alt_obs) < P$maf
  dropped <- rbind(dropped, data.table::data.table(k = m$k[r_miss], reason = paste0("miss_rate > ", P$miss_max)),
                   data.table::data.table(k = m$k[r_maf], reason = paste0("MAF < ", P$maf)))
  m <- m[!(r_miss | r_maf)]
  if (!nrow(m)) stop("no SNP left after matching and QC")

  m[, p_alt := data.table::fifelse(a1_is_alt, A1Freq, 1 - A1Freq)]
  m[, fill_const := data.table::fifelse(is.finite(p_alt) & n_miss > 0,
                                        pmin(2, pmax(0, (2 * n * p_alt - sum_obs) / pmax(n_miss, 1))),
                                        sum_obs / pmax(n - n_miss, 1))]
  method <- if (all(is.finite(m$p_alt))) "scheme3" else if (any(is.finite(m$p_alt))) "scheme3_partial" else "mean"

  cells <- data.table::data.table(ID = character(0), IID = character(0), g = numeric(0))
  m[, in_ld := FALSE]
  if (!is.null(ld)) {
    pan <- .ld_panel(ld$ref, ld$dir)
    al <- .ld_align(pan, pvar)
    m[, in_ld := (row - 1L) %in% al$vi]
    if (any(m$in_ld)) {
      ua <- m[in_ld == TRUE & is.finite(p_alt)]
      gi <- genotype_impute(pfile, pan, ld$dir, snps = m$ID[m$in_ld], flank = P$flank, sample = iid[si],
                            delta = P$delta, anchor_min_miss = P$anchor_min_miss, thresh = P$thresh,
                            threads = P$threads, mem_gb = P$mem_gb, af = stats::setNames(ua$p_alt, ua$ID))
      cells <- gi$fill[, .(ID = snp, IID, g)]
      method <- paste0("scheme4+", method)
    }
  }
  dropped <- merge(dropped, snp[, `:=`(k = .I)][, c("k", "CHR", bpcol, "A1", "A2"), with = FALSE], by = "k")
  list(snp = m[, .(ID, idx = row, CHR = CHROM, POS, REF, ALT, A1, A2, a1_is_alt, n_miss, miss_rate, af_alt_obs,
                   fill_const, in_ld)],
       cells = cells, dropped = dropped[, !"k"], method = method, para = P,
       time = proc.time()[["elapsed"]] - t0)
}

#' Defaults of [genotype_fill()]
#'
#' build NULL = the one BP column present; maf 0.01 (observed); miss_max 1 = keep all (drop SNPs whose
#' missing rate exceeds it); delta 0.01 ridge of stage 1; anchor_min_miss 0.01 (stage 2 only at or above
#' this missing rate); flank 2500 panel rows each side; thresh 0 = all stored eigenpairs; threads; mem_gb
#' memory budget for the people batches of stage 1.
#' @export
genotype_fill_para <- function() {
  list(build = NULL, maf = 0.01, miss_max = 1, delta = 0.01, anchor_min_miss = 0.01, flank = 2500L,
       thresh = 0, threads = 4L, mem_gb = 2)
}
