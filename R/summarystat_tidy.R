.read_table <- function(x) {
  if (is.data.frame(x)) return(data.table::as.data.table(x))
  if (grepl("\\.parquet$", x, ignore.case = TRUE)) return(data.table::as.data.table(arrow::read_parquet(x)))
  data.table::fread(x)
}

#' Align GWAS summary statistics to the SNP table of an LD panel and recompute Z
#'
#' Required columns (exact names): CHR, BP37 or BP38 (as `build`), A1, A2, A1Freq, BETA, SE, N.
#' Input may be txt/tsv/csv (optionally .gz), parquet, or a data.frame.
#' The GWAS must have SNPs on every chromosome 1 to 22 (an error otherwise); X is optional: with it the panel's X
#' SNPs are aligned too, without it the X part of the panel is left out. CHR is read as the chromosome code of
#' [geno_file()] (X, chrX and 23 are the same).
#' Matching is by chromosome, position and alleles against the panel; a panel SNP without a position on `build`
#' is matched to nothing. Rows whose A1/A2 are the panel's A2/A1 are converted (BETA -> -BETA,
#' A1Freq -> 1 - A1Freq). Z is recomputed from t = BETA/SE with df = N - 2.
#' QC as SBayesRC::tidy: |A1Freq - panel A1Freq| <= freq_thresh and N within mean +- N_sd_range SD
#' (rate2pq is not applied: it depends on the scale of BETA).
#' An optional column `r2` (imputation quality; 1 = typed) is carried through, so [sumstat_impute()] uses only
#' rows with r2 == 1 as typed and re-imputes the rest; without it every row is typed (r2 = 1).
#' @param gwas file path or data.frame.
#' @param ref the LD panel: its folder (`snpinfo.parquet` is read), a parquet path, or a data.frame with CHR,
#'   BP37 and/or BP38, A1, A2, A1Freq, Block, Index.
#' @param build 37 or 38; the panel must carry that build's position column.
#' @param output optional output path (.parquet or text).
#' @return data.table, one row per aligned panel SNP, by Block and Index: CHR, BP37, BP38, A1, A2 (the panel's),
#'   A1Freq (of the GWAS, for the panel's A1), Z, N, r2, Block, Index.
#' @export
summarystat_tidy <- function(gwas, ref, build = 37, freq_thresh = 0.2, N_sd_range = 3, output = NULL) {
  g <- .read_table(gwas)
  bp <- paste0("BP", build)
  need <- c("CHR", bp, "A1", "A2", "A1Freq", "BETA", "SE", "N")
  miss <- setdiff(need, names(g))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  has_r2 <- "r2" %in% names(g)
  g <- g[, c(need, if (has_r2) "r2"), with = FALSE]
  if (!has_r2) g[, r2 := 1]
  data.table::setnames(g, bp, "BP")
  g[, `:=`(CHR = .chr_code(CHR), BP = as.integer(BP), A1 = .by_unique(as.character(A1), toupper),
           A2 = .by_unique(as.character(A2), toupper), r2 = as.numeric(r2))]
  for (v in c("A1Freq", "BETA", "SE", "N")) data.table::set(g, j = v, value = as.numeric(g[[v]]))
  message(nrow(g), " SNPs in summary data")
  g <- g[!is.na(CHR) & !is.na(BP) & is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(A1Freq) & is.finite(N) & N > 2 &
           is.finite(r2)]
  if (has_r2) message(sum(g$r2 < 1), " rows with r2 < 1 (imputed): kept, flagged for re-imputation")
  message(nrow(g), " valid SNPs")
  .need_autosomes(g$CHR, "the GWAS")

  r <- .ld_panel(ref)
  if (!bp %in% names(r)) stop("the panel has no ", bp, " column: it holds only ",
                              paste(intersect(c("BP37", "BP38"), names(r)), collapse = ", "), call. = FALSE)
  r <- .with_bp(r)[CHR %in% unique(g$CHR)]  # the chromosomes of the GWAS: X (23) only when the GWAS has it
  rk <- data.table::data.table(pi = seq_len(nrow(r)), CHR = r$CHR, BP = r[[bp]], A1 = r$A1, A2 = r$A2)[!is.na(BP)]
  same <- rk[g, on = .(CHR, BP, A1, A2), nomatch = NULL]
  swap <- rk[g, on = .(CHR, BP, A1 = A2, A2 = A1), nomatch = NULL]
  swap[, `:=`(BETA = -BETA, A1Freq = 1 - A1Freq)]
  message(nrow(same), " SNPs aligned as is; ", nrow(swap), " SNPs with A1/A2 swapped (converted)")
  x <- rbind(same[, .(pi, A1Freq, BETA, SE, N, r2)], swap[, .(pi, A1Freq, BETA, SE, N, r2)])
  dup <- unique(x$pi[duplicated(x$pi)])
  if (length(dup)) { x <- x[!pi %in% dup]; message(length(dup), " panel SNPs matched more than one GWAS row; removed") }
  x <- x[abs(A1Freq - r$A1Freq[pi]) <= freq_thresh]
  message(nrow(x), " SNPs passed allele-frequency check (threshold ", freq_thresh, ")")
  mN <- mean(x$N); sN <- stats::sd(x$N)
  x <- x[N >= mN - N_sd_range * sN & N <= mN + N_sd_range * sN]
  message(nrow(x), " SNPs with N within mean +- ", N_sd_range, " SD (mean ", round(mN), ", SD ", round(sN), ")")
  t <- x$BETA / x$SE
  z <- -sign(t) * stats::qnorm(stats::pt(-abs(t), df = x$N - 2, log.p = TRUE), log.p = TRUE)
  p <- r[x$pi]
  out <- data.table::data.table(CHR = p$CHR, BP37 = p$BP37, BP38 = p$BP38, A1 = p$A1, A2 = p$A2, A1Freq = x$A1Freq,
                                Z = z, N = x$N, r2 = x$r2, Block = p$Block, Index = p$Index)
  data.table::setorder(out, Block, Index)
  message(nrow(out), " SNPs after tidy; ", nrow(r), " panel SNPs on the chromosomes of the GWAS (",
          formatC(100 * (1 - nrow(out) / nrow(r)), format = "f", digits = 2), "% missing)")
  if (!is.null(output)) {
    if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
  }
  out[]
}
