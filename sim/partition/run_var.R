# usage: Rscript run_var.R <ld_dir> <snp_table> <sim_dir>  -> sim_dir/<method>.tsv (ID, BETA_STD, BETA_STD_<annot>)
#   and sim_dir/<method>_h2.tsv (annotation, h2_shapley = Shapley value of the column in beta' R beta, its sd,
#   var_shapley = b_k' R b, h2_share = model share)
# Methods: V viprs; A1 viprs_annot, one Gaussian slab; AM 6:3:1 slab mixture (final Gibbs fit); AMV AM with the VI
# refit as the final fit;
# AML AM started at ldsc_strat(). Intercept (sigma2_0) on every SNP.
# env: SHARED (directory of annot.tsv, default sim_dir), GIBBS_DRAWS (draws of the final Gibbs fit of A1, AM, AML; default 300),
#      PI_BY (criterion of the pi grid: elbo (default) or pseudo_r), CUTOFFS (eigen cutoffs, default 0.995,0.999),
#      GRID_EXTEND (pi grid points above the official top, default 4), COLS (annotation columns given to the model, comma or +
#      separated; unset all), BINS (LD-bin file r2 > 0.98 for the snp table; made with ld_bins() when missing), SKIP (annotation columns of
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
an <- fread(file.path(Sys.getenv("SHARED", a[3]), "annot.tsv"))[, ID := NULL]
if (nzchar(Sys.getenv("COLS"))) an <- an[, c("Block", "Index", strsplit(Sys.getenv("COLS"), "[,+]")[[1]]), with = FALSE]
bins <- Sys.getenv("BINS", file.path(dirname(a[2]), "ldbin0.98.bin"))
if (!file.exists(bins)) ld_bins(a[1], r2 = 0.98, ref = rf, output = bins, threads = 4)
skip <- if (Sys.getenv("SKIP", "\001") == "\001") NULL else setdiff(strsplit(Sys.getenv("SKIP"), ",")[[1]], "")
meth <- strsplit(Sys.getenv("METHODS", "V,A1,AM,AMV,AML"), ",")[[1]]
tm <- function(e) { t0 <- Sys.time(); x <- e; list(x, as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
ls <- NULL
if ("AML" %in% meth) {
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
gd <- as.integer(Sys.getenv("GIBBS_DRAWS", "300")); pb <- Sys.getenv("PI_BY", "elbo")
ge <- as.integer(Sys.getenv("GRID_EXTEND", "4"))
co <- if (nzchar(Sys.getenv("CUTOFFS"))) as.numeric(strsplit(Sys.getenv("CUTOFFS"), ",")[[1]]) else c(0.995, 0.999)
fit <- function(k) switch(k,
  V = viprs(ss, a[1], ref = rf, threads = 4),
  A1 = viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip,
                   control = viprs_annot_control(mix_prop = 1, mix_var = 1, gibbs_draws = gd, pi_by = pb, cutoffs = co, grid_extend = ge), threads = 4),
  AM = viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip,
                   control = viprs_annot_control(gibbs_draws = gd, pi_by = pb, cutoffs = co, grid_extend = ge), threads = 4),
  AMV = viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip,
                    control = viprs_annot_control(gibbs_draws = 0, pi_by = pb, cutoffs = co, grid_extend = ge), threads = 4),
  AML = viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip,
                    sigma2_init = ls, control = viprs_annot_control(gibbs_draws = gd, pi_by = pb, cutoffs = co, grid_extend = ge), threads = 4))
for (k in meth) {
  f <- tm(fit(k)); r <- ids[f[[1]], on = .(Block, Index)]
  bc <- grep("^BETA_STD", names(r), value = TRUE)
  fwrite(r[, c("ID", bc), with = FALSE], file.path(a[3], paste0(k, ".tsv")), sep = "\t")
  hp <- attr(f[[1]], "h2_shapley")
  if (!is.null(hp)) fwrite(data.table(annotation = names(hp), h2_shapley = hp, h2_shapley_sd = attr(f[[1]], "h2_shapley_sd"),
                                      var_shapley = attr(f[[1]], "var_shapley"), h2_share = attr(f[[1]], "h2_share")),
                           file.path(a[3], paste0(k, "_h2.tsv")), sep = "\t")
  cat(k, "time", round(f[[2]], 1), "pi", signif(attr(f[[1]], "hyp")$pi, 3), "h2", round(attr(f[[1]], "hyp")$h2, 3),
      "suppressed", if (is.null(attr(f[[1]], "n_suppressed"))) 0 else attr(f[[1]], "n_suppressed"), "\n")
  g <- attr(f[[1]], "grid")
  if (!is.null(g)) { cat(k, "grid (cutoff", attr(f[[1]], "cutoff"), "chosen)\n"); options(width = 200); print(as.data.frame(g), digits = 4, row.names = FALSE) }
}
