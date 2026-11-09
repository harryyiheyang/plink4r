.read_table <- function(x) {
  if (is.data.frame(x)) return(data.table::as.data.table(x))
  if (grepl("\\.parquet$", x, ignore.case = TRUE)) return(data.table::as.data.table(arrow::read_parquet(x)))
  data.table::fread(x)
}

#' Align GWAS summary statistics to the plink4r SNP reference and recompute Z
#'
#' Required columns (exact names): CHR, BP37 or BP38 (as `build`), A1, A2, A1Freq, BETA, SE, N.
#' Input may be txt/tsv/csv (optionally .gz), parquet, or a data.frame.
#' Matching is by CHR + BP + alleles against `ref`; rows whose A1/A2 are the reference ALT/REF swapped are
#' converted (BETA -> -BETA, A1Freq -> 1 - A1Freq). Z is recomputed from t = BETA/SE with df = N - 2.
#' QC as SBayesRC::tidy: |A1Freq - AF_ref| <= freq_thresh and N within mean +- N_sd_range SD
#' (rate2pq is not applied: it depends on the scale of BETA).
#' @param gwas file path or data.frame.
#' @param ref reference parquet (columns ID, chr, bp37, bp38, ref, alt, AF_PanUKB, block, ld_index) or data.frame.
#' @param build 37 or 38.
#' @param output optional output path (.parquet or text).
#' @return data.table: ID, chr, bp37, bp38, A1 (= ref alt), A2 (= ref ref), freq, Z, N, block, ld_index.
#' @export
summarystat_tidy <- function(gwas, ref, build = 37, freq_thresh = 0.2, N_sd_range = 3, output = NULL) {
  g <- .read_table(gwas)
  bp <- paste0("BP", build)
  need <- c("CHR", bp, "A1", "A2", "A1Freq", "BETA", "SE", "N")
  miss <- setdiff(need, names(g))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
  g <- g[, ..need]
  data.table::setnames(g, bp, "BP")
  for (v in c("CHR", "BP")) data.table::set(g, j = v, value = as.integer(g[[v]]))
  for (v in c("A1Freq", "BETA", "SE", "N")) data.table::set(g, j = v, value = as.numeric(g[[v]]))
  g[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]
  message(nrow(g), " SNPs in summary data")
  g <- g[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(A1Freq) & is.finite(N) & N > 2]
  message(nrow(g), " valid SNPs")

  r <- .read_table(ref)
  r <- r[, .(ID, chr, bp37, bp38, ref, alt, AF_PanUKB, block, ld_index)]
  r[, BP := if (build == 37) bp37 else bp38]
  same <- merge(g, r, by.x = c("CHR", "BP", "A1", "A2"), by.y = c("chr", "BP", "alt", "ref"))
  swap <- merge(g, r, by.x = c("CHR", "BP", "A1", "A2"), by.y = c("chr", "BP", "ref", "alt"))
  swap[, `:=`(BETA = -BETA, A1Freq = 1 - A1Freq)]
  message(nrow(same), " SNPs aligned as is; ", nrow(swap), " SNPs with A1/A2 swapped (converted)")
  x <- rbind(same[, .(ID, chr = CHR, bp37, bp38, freq = A1Freq, BETA, SE, N, AF_PanUKB, block, ld_index)],
             swap[, .(ID, chr = CHR, bp37, bp38, freq = A1Freq, BETA, SE, N, AF_PanUKB, block, ld_index)])
  dup <- x[duplicated(ID), unique(ID)]
  if (length(dup)) { x <- x[!ID %in% dup]; message(length(dup), " reference SNPs matched more than one GWAS row; removed") }
  x <- x[abs(freq - AF_PanUKB) <= freq_thresh]
  message(nrow(x), " SNPs passed allele-frequency check (threshold ", freq_thresh, ")")
  mN <- mean(x$N); sN <- stats::sd(x$N)
  x <- x[N >= mN - N_sd_range * sN & N <= mN + N_sd_range * sN]
  message(nrow(x), " SNPs with N within mean +- ", N_sd_range, " SD (mean ", round(mN), ", SD ", round(sN), ")")
  t <- x$BETA / x$SE
  x[, Z := -sign(t) * stats::qnorm(stats::pt(-abs(t), df = N - 2, log.p = TRUE), log.p = TRUE)]
  r2 <- r[, .(ID, A1 = alt, A2 = ref)]
  x <- r2[x, on = "ID"]
  x <- x[, .(ID, chr, bp37, bp38, A1, A2, freq, Z, N, block, ld_index)]
  data.table::setorder(x, block, ld_index)
  message(nrow(x), " SNPs after tidy; ", nrow(r), " SNPs in reference (", formatC(100 * (1 - nrow(x) / nrow(r)), format = "f", digits = 2), "% missing)")
  if (!is.null(output)) {
    if (grepl("\\.parquet$", output)) arrow::write_parquet(x, output) else data.table::fwrite(x, output, sep = "\t")
  }
  x[]
}
