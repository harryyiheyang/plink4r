# usage: Rscript run_var.R <ld_dir> <snp_table> <sim_dir>  -> sim_dir/<method>.tsv (ID, BETA_STD, BETA_STD_<annot>)
#   and sim_dir/<method>_h2.tsv (annotation, h2_shapley = Shapley value of the column in beta' R beta, h2_share =
#   model share)
# Methods (one path each: VI chooses the cutoff and pi, the final fit is Gibbs + RB; intercept sigma2_0 on every SNP):
#   ORACLE viprs_annot with sigma2_k and pi fixed at the oracle of gen_var.py (oracle.tsv: the design columns; every
#   other column at 0);
#   MASH viprs_annot (MASH prior on sigma2_k) started at ldsc_strat(); MASH_PI the same with pi fixed at the truth of
#   oracle.tsv (no pi grid); BASE viprs_annot with the intercept only (no annotation, no suppression): the
#   baseline the annotated fits must not fall below. Eigen cutoff: CUTOFFS, default 1 (every component of
#   eigen.bin, the panel the data are simulated from: no cutoff choice).
# env: SHARED (directory of annot.tsv, default sim_dir), GIBBS_DRAWS (draws of the final Gibbs fit; default 300),
#      CUTOFFS (eigen cutoffs, default 1),
#      GRID_EXTEND (pi grid points above the official top, default 4), GRID_MAX_ITER (iterations per grid point, default
#      1000: at 200 the 57-column fits stop unconverged for pi >= 0.015), COLS (annotation columns given to MASH, comma or +
#      separated; unset all), BINS (LD-bin file r2 > 0.98 for the snp table; made with ld_bins() when missing), SKIP (annotation columns of
#      the LD-bin suppression, comma separated; "" none, unset every column, as gen_var.py with BINS), METHODS (default all),
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
an0 <- an   # every column (ORACLE takes the design columns from it)
if (nzchar(Sys.getenv("COLS"))) an <- an[, c("Block", "Index", strsplit(Sys.getenv("COLS"), "[,+]")[[1]]), with = FALSE]
bins <- Sys.getenv("BINS", file.path(dirname(a[2]), "ldbin0.98.bin"))
if (!file.exists(bins)) ld_bins(a[1], r2 = 0.98, ref = rf, output = bins, threads = 4)
skip <- if (Sys.getenv("SKIP", "\001") == "\001") NULL else setdiff(strsplit(Sys.getenv("SKIP"), ",")[[1]], "")
meth <- strsplit(Sys.getenv("METHODS", "BASE,ORACLE,MASH,MASH_PI"), ",")[[1]]
tm <- function(e) { t0 <- Sys.time(); x <- e; list(x, as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
ls <- NULL
if (any(c("MASH", "MASH_PI") %in% meth)) {
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
gd <- as.integer(Sys.getenv("GIBBS_DRAWS", "300"))
ge <- as.integer(Sys.getenv("GRID_EXTEND", "4")); gm <- as.integer(Sys.getenv("GRID_MAX_ITER", "1000"))
co <- if (nzchar(Sys.getenv("CUTOFFS"))) as.numeric(strsplit(Sys.getenv("CUTOFFS"), ",")[[1]]) else 1
ctl <- function(...) viprs_annot_control(gibbs_draws = gd, cutoffs = co, grid_extend = ge, grid_max_iter = gm, ...)
pt <- function() { p <- fread(file.path(a[3], "oracle.tsv"))$pi[1]; list(pi_grid = p, cut_pi = p) }
fit <- function(k) switch(k,
  ORACLE = {
    # every column (the same LD-bin suppression as the data), sigma2 at the oracle: the design columns as in oracle.tsv,
    # the others 0
    o <- fread(file.path(a[3], "oracle.tsv")); ac0 <- setdiff(names(an0), c("Block", "Index"))
    s2 <- setNames(rep(0, length(ac0) + 1), c("intercept", ac0)); s2[o$annotation] <- o$sigma2
    viprs_annot(ss, a[1], annot = an0, ref = rf, bin_file = bins, skip_annot = skip, sigma2_init = s2,
                control = ctl(fix_sigma2 = TRUE, start_scale = 1, start_zero_share = 0,
                              pi_grid = o$pi[1], cut_pi = o$pi[1]),
                threads = 4)
  },
  BASE = viprs_annot(ss, a[1], annot = an[, .(Block, Index)], ref = rf, bin_file = bins, skip_annot = character(0),
                     control = ctl(), threads = 4),
  MASH = viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip, sigma2_init = ls,
                     control = ctl(), threads = 4),
  MASH_PI = viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bins, skip_annot = skip, sigma2_init = ls,
                        control = do.call(ctl, pt()), threads = 4),
  stop("unknown method ", k))
for (k in meth) {
  f <- tm(fit(k)); r <- ids[f[[1]], on = .(Block, Index)]
  bc <- grep("^BETA_STD", names(r), value = TRUE)
  fwrite(r[, c("ID", bc), with = FALSE], file.path(a[3], paste0(k, ".tsv")), sep = "\t")
  hp <- attr(f[[1]], "h2_shapley")
  if (!is.null(hp)) fwrite(data.table(annotation = names(hp), h2_shapley = hp, h2_share = attr(f[[1]], "h2_share")),
                           file.path(a[3], paste0(k, "_h2.tsv")), sep = "\t")
  cat(k, "time", round(f[[2]], 1), "pi", signif(attr(f[[1]], "hyp")$pi, 3), "h2", round(attr(f[[1]], "hyp")$h2, 3),
      "suppressed", if (is.null(attr(f[[1]], "n_suppressed"))) 0 else attr(f[[1]], "n_suppressed"), "\n")
  ms_ <- attr(f[[1]], "mash")
  if (!is.null(ms_) && length(ms_$mean_atom) > 1 && !all(ms_$q == 0))
    cat(k, "mash sigma2_0", signif(ms_$sigma2_0, 3), "\n  E_q[a_k]", paste(round(ms_$mean_atom[-1], 2), collapse = " "),
        "\n  q(a_k=0)", paste(round(ms_$q[-1, 1], 2), collapse = " "), "\n  pzero_gibbs", paste(round(ms_$pzero_gibbs[-1], 2), collapse = " "), "\n")
  g <- attr(f[[1]], "grid")
  if (!is.null(g)) { cat(k, "grid (cutoff", attr(f[[1]], "cutoff"), "chosen)\n"); options(width = 200); print(as.data.frame(g), digits = 4, row.names = FALSE) }
}
