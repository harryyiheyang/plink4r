# usage: Rscript run_ideal.R <ld_dir> <snp_table> <sim_dir>  -> sim_dir/{V,A1,AM,AML}.tsv
# V viprs; A1 viprs_annot with one Gaussian slab; AM viprs_annot with the 6:3:1 slab mixture; AML AM started at
# ldsc_strat(). Residual annotation only on unannotated SNPs (RESID_ALL=1: on every SNP).
a <- commandArgs(TRUE); suppressPackageStartupMessages({library(plink4r); library(data.table)})
RA <- Sys.getenv("RESID_ALL", "0") == "1"
tab <- fread(a[2])
rf <- tab[, .(CHR = chr, BP37 = bp37, A1, A2, A1Freq = freq, Block = block, Index = ld_index)]
ids <- tab[, .(ID, Block = block, Index = ld_index)]
ma <- fread(file.path(a[3], "trait.ma"))
ss <- ids[ma[, .(ID = SNP, Z = b / se, N)], on = "ID"][rf, on = .(Block, Index), nomatch = NULL][, .(A1, A2, Z, N, Block, Index)]
an <- fread(file.path(a[3], "annot.tsv"))[, ID := NULL]
bin0 <- file.path(a[3], "ldbin0.bin"); if (!file.exists(bin0)) writeBin(rep(-1L, nrow(rf)), bin0, size = 2)
tm <- function(e) { t0 <- Sys.time(); x <- e; list(x, as.numeric(difftime(Sys.time(), t0, units = "secs"))) }
t0 <- Sys.time()
L <- ld_score(a[1], annot = an[rf[, .(Block, Index, CHR, BP37)], on = .(Block, Index)][, -c("Block", "Index")],
              ref = rf, build = 37, threads = 4)
ls <- ldsc_strat(ss, L)
cat("ldsc_strat", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s; tau",
    paste(ls$Parameter_summary$Parameter, signif(ls$Parameter_summary$Estimate, 3), collapse = " "), "\n")
one <- viprs_annot_control(mix_prop = 1, mix_var = 1)
fits <- list(
  A1 = tm(viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bin0, resid_all = RA, control = one, threads = 4)),
  AM = tm(viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bin0, resid_all = RA, threads = 4)),
  AML = tm(viprs_annot(ss, a[1], annot = an, ref = rf, bin_file = bin0, resid_all = RA, sigma2_init = ls, threads = 4)),
  V = tm(viprs(ss, a[1], ref = rf, threads = 4)))
for (k in names(fits)) {
  r <- ids[fits[[k]][[1]], on = .(Block, Index)]
  bc <- grep("^BETA_STD", names(r), value = TRUE)
  fwrite(r[, c("ID", bc), with = FALSE], file.path(a[3], paste0(k, ".tsv")), sep = "\t")
  sh <- attr(fits[[k]][[1]], "h2_share")
  cat(k, "time", round(fits[[k]][[2]], 1), "pi", signif(attr(fits[[k]][[1]], "hyp")$pi, 3),
      if (!is.null(sh)) paste("share", paste(names(sh), round(sh, 3), collapse = " ")), "\n")
}
