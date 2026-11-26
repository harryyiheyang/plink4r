# usage: Rscript run_var.R <ld_dir> <snp_table> <sim_dir>  -> sim_dir/<method>.tsv (ID, BETA_STD, BETA_STD_<annot>)
#   and sim_dir/<method>_h2.tsv (annotation, h2_annot = (sum_j A_jk) sigma2_k, h2_shapley = b^(k)' R b for P / PL)
# Methods: V viprs; A viprs_annot (the variational selection step, moment M-step); P partprs (A, then short Gibbs + RB);
# PL partprs started at ldsc_strat(). The intercept (residual PRS) is on every SNP.
# env: BINS (LD-bin file r2 > 0.98 for the snp table; made with ld_bins() when missing), SKIP (annotation columns of
#      the LD-bin suppression, comma separated; "" none, unset every 0/1 column), METHODS (default all),
#      LDSCORE (parquet / tsv of precomputed stratified LD scores by ID for the ldsc_strat() start: a base column
#      L2 / L_base / base / baseL2 and one column per annotation named L_<annot>, L2_<annot>, <annot>L2 or <annot>;
#      unset: ld_score()).
a <- commandArgs(TRUE); suppressPackageStartupMessages({library(plink4r); library(data.table)})
tab <- fread(a[2])
rf <- tab[, .(CHR = chr, BP37 = bp37, A1, A2, A1Freq = freq, Block = block, Index = ld_index)]
ids <- tab[, .(ID, Block = block, Index = ld_index)]
ma <- fread(file.path(a[3], "trait.ma"))
ss <- ids[ma[, .(ID = SNP, Z = b / se, N)], on = "ID"][rf, on = .(Block, Index), nomatch = NULL][, .(A1, A2, Z, N, Block, Index)]
an <- fread(file.path(a[3], "annot.tsv"))[, ID := NULL]
bins <- Sys.getenv("BINS", file.path(dirname(a[2]), "ldbin0.98.bin"))
if (!file.exists(bins)) ld_bins(a[1], r2 = 0.98, ref = rf, output = bins, threads = 4)
skip <- if (Sys.getenv("SKIP", "\001") == "\001") NULL else setdiff(strsplit(Sys.getenv("SKIP"), ",")[[1]], "")
meth <- strsplit(Sys.getenv("METHODS", "V,A,P,PL"), ",")[[1]]
tm <- function(e) { t0 <- Sys.time(); x <- e; list(x, as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
ls <- NULL
if ("PL" %in% meth) {
  t0 <- Sys.time()
  acols <- setdiff(names(an), c("Block", "Index"))
  if (nzchar(Sys.getenv("LDSCORE"))) {
    f <- Sys.getenv("LDSCORE")
    l <- if (grepl("parquet$", f)) as.data.table(arrow::read_parquet(f)) else fread(f)
    pick <- function(cands) { h <- intersect(cands, names(l)); if (!length(h)) NA_character_ else h[1] }
    base <- pick(c("L2", "L_base", "base", "baseL2", "base_L2", "baseline", "baselineL2"))
    lc <- vapply(acols, function(x) pick(c(paste0("L_", x), paste0("L2_", x), paste0(x, "L2"), paste0(x, "_L2"), x)), "")
    if (is.na(base) || anyNA(lc)) stop("LDSCORE: no column for ", paste(c(if (is.na(base)) "base L2", acols[is.na(lc)]),
                                       collapse = ", "), "; its columns: ", paste(names(l), collapse = " "))
    l <- l[, c("ID", base, lc), with = FALSE]; setnames(l, c("ID", "L2", paste0("L2_", acols)))
    L <- l[ids, on = "ID", nomatch = NULL][rf[, .(Block, Index, A1, A2)], on = .(Block, Index), nomatch = NULL]
    L <- an[L, on = .(Block, Index)]
    setnames(L, acols, paste0("A_", acols)); L[, ID := NULL]
    cat("LDSCORE: ", nrow(L), " SNPs with precomputed LD scores\n")
  } else L <- ld_score(a[1], annot = an[rf[, .(Block, Index, CHR, BP37)], on = .(Block, Index)][, -c("Block", "Index")],
                       ref = rf, build = 37, threads = 4)
  ls <- ldsc_strat(ss, L)
  cat("ldsc_strat", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
}
fit <- function(k) switch(k,
  V = viprs(ss, a[1], ref = rf, threads = 4),
  A = viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip,
                  control = viprs_annot_control(exact_em = FALSE), threads = 4),
  P = partprs(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip, threads = 4),
  PL = partprs(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip, sigma2_init = ls, threads = 4))
for (k in meth) {
  f <- tm(fit(k)); r <- ids[f[[1]], on = .(Block, Index)]
  bc <- grep("^BETA_STD", names(r), value = TRUE)
  fwrite(r[, c("ID", bc), with = FALSE], file.path(a[3], paste0(k, ".tsv")), sep = "\t")
  ha <- attr(f[[1]], "h2_annot")
  if (!is.null(ha)) {
    hs <- attr(f[[1]], "h2_shapley")
    fwrite(data.table(annotation = names(ha), h2_annot = ha, h2_shapley = if (is.null(hs)) NA_real_ else hs[names(ha)]),
           file.path(a[3], paste0(k, "_h2.tsv")), sep = "\t")
  }
  if (!is.null(attr(f[[1]], "trace"))) fwrite(attr(f[[1]], "trace"), file.path(a[3], paste0(k, "_trace.tsv")), sep = "\t")
  cat(k, "time", round(f[[2]], 1), "pi", signif(attr(f[[1]], "hyp")$pi, 3), "h2", round(attr(f[[1]], "hyp")$h2, 3),
      "suppressed", if (is.null(attr(f[[1]], "n_suppressed"))) 0 else attr(f[[1]], "n_suppressed"), "\n")
}
