#' Impute hard-call missing genotypes from eigen.bin LD (whole block or a SNP window)
#'
#' Works on the panel SNPs that are in the PGEN (the intersection S; panel SNPs absent from the PGEN are not
#' carried). The panel need not share the PGEN's REF/ALT orientation, IDs or build: SNPs are aligned by ID and
#' alleles, else by chromosome, position and allele pair, and a swapped pair enters with a negative sign.
#' Stage 1 (OpenMP over people): z_M = R_MO (R_OO + delta I)^{-1} z_O per person on S, standardized with the
#' anchor AF (AF_PanUKB in PGEN orientation, or `af`), from C_S = I - delta (R_SS + delta I)^{-1}, built by
#' Woodbury in the k eigenpairs or directly in |S|, whichever is smaller; delta must be > 0.
#' Window: S is the panel rows of the requested SNPs plus `flank` eigen rows on each side inside the same LD
#' block (blocks are independent in the panel, so a region spanning blocks is imputed block by block).
#' Stage 2: for SNPs with missing rate >= anchor_min_miss, shift the filled cells so that observed + filled
#' allele count = 2 n p, then clip to [0, 2]; below the cutoff stage 1 is kept. anchor_min_miss < 0 disables it.
#' A `geno` object (hard calls) may be given in place of the PGEN: its columns, counted allele (ALT = `a1`)
#' and missing cells take the place of the .pvar and the file's missing calls, every row is imputed (subset
#' first with [geno_subset()]), and BED-based objects work the same way. Dosage objects are not accepted.
#' @param pfile PLINK2 prefix (.pgen/.pvar/.psam), or a `geno` object.
#' @param ref panel SNP table; NULL (default) reads the .parquet in `ld_dir`. A parquet path or data.frame with
#'   ID, block, ld_index (0-based eigen row), AF_PanUKB (frequency of alt) and, for orientation and position
#'   matching, chr, ref, alt and bp37 and/or bp38 (one column: the PGEN must be on that build; both: the build
#'   placing more SNPs is used). Without ref/alt the panel is matched by ID alone in the PGEN orientation.
#' @param ld_dir directory with block<N>.eigen.bin and the panel's SNP .parquet.
#' @param block block number for whole-block imputation (used when `snps` is NULL).
#' @param snps IDs to impute (panel or PGEN IDs, e.g. a fine-mapping region); NULL = whole `block`.
#' @param flank eigen rows added on each side of `snps` within each block.
#' @param sample IIDs to impute (NULL = all .psam rows).
#' @param thresh eigenvalue variance kept; <= 0 keeps every stored pair.
#' @param nb_people people per batch (NULL = largest batch whose S_i, rhs and SNP lists fit in `mem_gb`).
#' @param mem_gb memory budget of the people batches; 0 takes it from [plink4r_resources()].
#' @param af optional anchor frequencies of the PGEN ALT allele, named by PGEN ID; they replace AF_PanUKB.
#' @return list(fill = data.table(IID, snp, g1, g), snp = per-SNP n_miss / c / shift, info, time); `snp` is
#'   the PGEN ID and g the ALT dosage. For a `geno` object `snp` is the column name, g the dosage of `a1`,
#'   and `fill` also has `i` (row) and `j` (column) of the object.
#' @export
genotype_impute <- function(pfile, ref = NULL, ld_dir, block = NULL, snps = NULL, flank = 2500L, sample = NULL,
                            delta = 0.01, anchor_min_miss = 0.01, thresh = 0, threads = 4L, nb_people = NULL,
                            mem_gb = 2, af = NULL) {
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  .check_delta(delta)
  obj <- inherits(pfile, "geno")
  if (obj) {
    if (!is.null(sample)) stop("with a geno object, sample must be NULL (use geno_subset() first)", call. = FALSE)
    if (pfile$dosage) stop("LD imputation on an object needs hard calls; for a PGEN with dosages give the PGEN prefix.", call. = FALSE)
    X <- pfile
    pvar <- .geno_vartab(X)  # appended columns have allele_ct 0: never aligned to the panel
    iid <- X$sample
    si <- seq_len(X$n) - 1L
    what <- "the geno object"
  } else {
    pvar <- .read_pvar5(pfile)
    iid <- .psam_ids(pfile)
    si <- .sample_rows(sample, iid) - 1L
    what <- paste0(pfile, ".pvar")
  }
  if (is.null(snps)) {
    if (is.null(block)) stop("give `block` or `snps`")
    al <- .ld_align(.ld_panel(ref, ld_dir, blocks = block), pvar)
    jobs <- list(list(block = block, lo = NULL, hi = NULL, keep = NULL))
  } else {
    al <- .ld_align(.ld_panel(ref, ld_dir), pvar)  # whole panel: `snps` may be panel or PGEN IDs
    rs <- al[al$ID %in% snps | pvar$ID[al$vi + 1L] %in% snps]
    if (nrow(rs) < length(unique(snps))) warning(length(unique(snps)) - nrow(rs), " of `snps` not in both `ref` and the PGEN")
    if (!nrow(rs)) stop("none of `snps` is in both `ref` and ", what)
    jobs <- lapply(sort(unique(rs$block)), function(b) {
      li <- rs$ld_index[rs$block == b]
      list(block = b, lo = min(li) - flank, hi = max(li) + flank, keep = rs$vi[rs$block == b])
    })
  }
  if (!is.null(af)) {  # user anchors are PGEN-ALT frequencies; the panel column is in panel orientation
    pa <- as.numeric(af[pvar$ID[al$vi + 1L]])
    ok <- is.finite(pa)
    al$AF_PanUKB[ok] <- ifelse(al$sign[ok] > 0, pa[ok], 1 - pa[ok])
  }
  res <- lapply(jobs, function(jb) {
    blk <- jb$block
    ef <- file.path(ld_dir, paste0("block", blk, ".eigen.bin"))
    w <- .ld_rows(al, blk, jb$lo, jb$hi)
    w <- w[AF_PanUKB > 0 & AF_PanUKB < 1]
    if (!nrow(w)) stop("no SNP of block ", blk, " found in ", what)
    snp <- pvar$ID[w$vi + 1L]
    kr <- if (is.null(jb$keep)) rep(TRUE, nrow(w)) else w$vi %in% jb$keep
    p_alt <- ifelse(w$sign > 0, w$AF_PanUKB, 1 - w$AF_PanUKB)
    nbp <- if (is.null(nb_people)) 0L else as.integer(nb_people)
    out <- if (obj) {
      genotype_impute_geno_cpp(X$ptr, as.integer(pvar$row[w$vi + 1L] - 1L), normalizePath(ef), p_alt,
                               as.integer(w$sign), thresh, delta, anchor_min_miss, kr, as.integer(threads), nbp,
                               mem_gb, as.integer(w$ld_index))
    } else {
      genotype_impute_cpp(normalizePath(paste0(pfile, ".pgen")), length(iid), pvar$allele_ct, as.integer(w$vi),
                          as.integer(si), normalizePath(ef), p_alt, as.integer(w$sign), thresh, delta,
                          anchor_min_miss, kr, as.integer(threads), nbp, mem_gb, as.integer(w$ld_index))
    }
    f <- data.table::as.data.table(out$fill)
    col <- if (obj) pvar$row[w$vi + 1L][f$j] else NULL
    f[, `:=`(IID = iid[si[i] + 1L], snp = snp[j])]
    f <- if (obj) f[, .(IID, snp, g1, g, i, j = col)] else f[, .(IID, snp, g1, g)]
    s <- data.table::as.data.table(out$snp)
    s[, snp := snp[j]]
    list(fill = f, snp = s[n_miss > 0, .(snp, n_miss, c_raw, c, mean_stage1, shift)],
         info = c(list(block = blk, window = if (is.null(jb$lo)) 0L else nrow(w), working_rows = nrow(w),
                       flipped = sum(w$sign < 0), threads = threads, mem_gb = mem_gb), out$info),
         time = out$time)
  })
  if (length(res) == 1L) return(res[[1]])
  list(fill = data.table::rbindlist(lapply(res, `[[`, "fill")), snp = data.table::rbindlist(lapply(res, `[[`, "snp")),
       info = lapply(res, `[[`, "info"), time = do.call(rbind, lapply(res, `[[`, "time")))
}
