#' Impute hard-call missing genotypes from eigen.bin LD (whole block or a SNP window)
#'
#' Stage 1 (OpenMP over people): z_M = R_MO (R_OO + delta I)^{-1} z_O per person on the working SNP set,
#' standardized with the anchor AF (AF_PanUKB). Whole block: R = U diag(lam) U' via C = U W U',
#' W = lam/(lam+delta) (delta = 0: generalized inverse). Window: the working set is the requested SNPs plus
#' `flank` SNPs on each side inside the same LD block (blocks are independent in the panel, so a region
#' spanning blocks is imputed block by block); C = I - delta (R_WW + delta I)^{-1}.
#' Stage 2: for SNPs with missing rate >= anchor_min_miss, shift the filled cells so that observed + filled
#' allele count = 2 n p, then clip to [0, 2]; below the cutoff stage 1 is kept. anchor_min_miss < 0 disables it.
#' @param pfile PLINK2 prefix (.pgen/.pvar/.psam); variant IDs must match `ref$ID`.
#' @param ref SNP reference parquet path or data.frame with ID, block, ld_index (0-based eigen row), AF_PanUKB.
#' @param ld_dir directory with block<N>.eigen.bin.
#' @param block block number for whole-block imputation (used when `snps` is NULL).
#' @param snps IDs to impute (e.g. a fine-mapping region); NULL = whole `block`.
#' @param flank SNPs added on each side of `snps` within each block (eigen-row order).
#' @param sample IIDs to impute (NULL = all .psam rows).
#' @param thresh eigenvalue variance kept; <= 0 keeps every stored pair.
#' @param nb_people people per batch (NULL = largest batch whose S_i, rhs and SNP lists fit in `mem_gb`).
#' @param nbr_frac,nbr_max rhs restriction by C_lj^2 (whole block only; loses accuracy, NULL = off).
#' @return list(fill = data.table(IID, snp, g1, g), snp = per-SNP n_miss / c / shift, info, time)
#' @export
genotype_impute <- function(pfile, ref, ld_dir, block = NULL, snps = NULL, flank = 2500L, sample = NULL,
                            delta = 0.01, anchor_min_miss = 0.01, thresh = 0, threads = 4L, nb_people = NULL,
                            mem_gb = 2, nbr_frac = NULL, nbr_max = NULL) {
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  r <- if (is.data.frame(ref)) data.table::as.data.table(ref) else
    data.table::as.data.table(arrow::read_parquet(ref, col_select = c("ID", "block", "ld_index", "AF_PanUKB")))
  pvar <- .read_pvar(pfile)
  psam <- data.table::fread(paste0(pfile, ".psam"))
  iid <- as.character(psam[[intersect(c("IID", "#IID"), names(psam))[1]]])
  si <- if (is.null(sample)) seq_along(iid) - 1L else match(as.character(sample), iid) - 1L
  if (anyNA(si)) stop(sum(is.na(si)), " samples not in .psam")
  act <- as.integer(2L + lengths(regmatches(pvar$ALT, gregexpr(",", pvar$ALT))))
  if (is.null(snps)) {
    if (is.null(block)) stop("give `block` or `snps`")
    jobs <- list(list(block = block, lo = NULL, hi = NULL, keep = NULL))
  } else {
    rs <- r[r$ID %in% snps]
    if (nrow(rs) < length(unique(snps))) warning(length(unique(snps)) - nrow(rs), " of `snps` not in `ref`")
    jobs <- lapply(sort(unique(rs$block)), function(b) {
      li <- rs$ld_index[rs$block == b]
      list(block = b, lo = min(li) - flank, hi = max(li) + flank, keep = rs$ID[rs$block == b])
    })
  }
  res <- lapply(jobs, function(jb) {
    blk <- jb$block
    rb <- r[r$block == blk]
    ef <- file.path(ld_dir, paste0("block", blk, ".eigen.bin"))
    con <- file(ef, "rb"); m <- readBin(con, "integer", 1L, size = 4L); close(con)
    snp <- rep(NA_character_, m); af <- rep(NA_real_, m)
    snp[rb$ld_index + 1L] <- rb$ID; af[rb$ld_index + 1L] <- rb$AF_PanUKB
    rows <- if (is.null(jb$lo)) integer(0) else seq.int(max(0L, jb$lo), min(m - 1L, jb$hi))
    if (length(rows)) { snp <- snp[rows + 1L]; af <- af[rows + 1L] }
    vi <- match(snp, pvar$ID) - 1L
    vi[is.na(vi)] <- -1L
    if (sum(vi >= 0) == 0) stop("no SNP of block ", blk, " found in ", pfile, ".pvar")
    kr <- if (is.null(jb$keep)) rep(TRUE, length(snp)) else snp %in% jb$keep
    af2 <- af; af2[vi < 0] <- 0.5
    out <- genotype_impute_cpp(normalizePath(paste0(pfile, ".pgen")), length(iid), act, as.integer(vi), as.integer(si),
                               normalizePath(ef), af2, thresh, delta, anchor_min_miss, kr, as.integer(threads),
                               if (is.null(nb_people)) 0L else as.integer(nb_people), mem_gb,
                               if (is.null(nbr_frac)) 0 else nbr_frac, if (is.null(nbr_max)) 0L else as.integer(nbr_max),
                               as.integer(rows))
    f <- data.table::as.data.table(out$fill)
    f[, `:=`(IID = iid[si[i] + 1L], snp = snp[j])]
    s <- data.table::as.data.table(out$snp)
    s[, snp := snp[j]]
    list(fill = f[, .(IID, snp, g1, g)], snp = s[!is.na(snp) & n_miss > 0, .(snp, n_miss, c_raw, c, mean_stage1, shift)],
         info = c(list(block = blk, window = length(rows), threads = threads, mem_gb = mem_gb), out$info), time = out$time)
  })
  if (length(res) == 1L) return(res[[1]])
  list(fill = data.table::rbindlist(lapply(res, `[[`, "fill")), snp = data.table::rbindlist(lapply(res, `[[`, "snp")),
       info = lapply(res, `[[`, "info"), time = do.call(rbind, lapply(res, `[[`, "time")))
}
