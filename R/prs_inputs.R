# Inputs shared by the summary-statistics PRS methods (viprs(), viprs_annot(), sbayesrc()): the GWAS rows on the
# eigen panel.
#
# sumstat: summarystat_tidy() / sumstat_impute() output (or a file of it), aligned to this panel: A1, A2 (the
# panel's), Z, N, Block, Index (r2 optional: rows with r2 <= 0 carry no data and are dropped). Returns x (by Block,
# Index: CHR BP37 BP38 A1 A2 A1Freq Z N Block Index, A1Freq = the panel's), blocks (ids with at least one SNP), keep
# (0-based eigen rows per block), n_panel (panel SNPs) and start (first row of each of `blocks` in the panel order
# Block, Index; per-SNP panel files are indexed by it).
.prs_inputs <- function(sumstat, ref, ld_dir) {
  s <- .read_table(sumstat)
  miss <- setdiff(c("A1", "A2", "Z", "N", "Block", "Index"), names(s))
  if (length(miss)) stop("sumstat lacks ", paste(miss, collapse = ", "), ": give the output of summarystat_tidy()",
                         call. = FALSE)
  s <- s[is.finite(Z) & is.finite(N) & N > 1]
  if ("r2" %in% names(s)) s <- s[is.finite(r2) & r2 > 0]
  s[, `:=`(Block = as.integer(Block), Index = as.integer(Index))]
  r <- .with_bp(.ld_panel(ref, ld_dir))
  data.table::setorder(r, Block, Index)
  i <- r[s, on = .(Block, Index), which = TRUE]
  if (anyNA(i) || anyDuplicated(i) || !isTRUE(all(r$A1[i] == as.character(s$A1) & r$A2[i] == as.character(s$A2))))
    stop("sumstat does not belong to this panel (its Block, Index, A1, A2 must name panel SNPs, each once): ",
         "make it with summarystat_tidy() and this panel", call. = FALSE)
  x <- r[i, .(CHR, BP37, BP38, A1, A2, A1Freq, Block, Index)]
  x[, `:=`(Z = as.numeric(s$Z), N = as.numeric(s$N))]
  data.table::setorder(x, Block, Index)
  if (!nrow(x)) stop("no sumstat row on the panel", call. = FALSE)
  message(nrow(x), " SNPs on the panel")
  bl <- x[, .(keep = list(Index)), by = Block]
  pn <- r[, .N, keyby = Block]
  st <- c(0, cumsum(as.numeric(pn$N)))[match(bl$Block, pn$Block)]
  list(x = x, blocks = bl$Block, keep = bl$keep, n_panel = nrow(r), start = st)
}

# Annotation matrix of the fitted SNPs x (.prs_inputs()$x): annot keyed by Block and Index, or by CHR and a position
# column (BP37 or BP38; with both, BP37), one numeric column per annotation; SNPs without a row get 0.
.prs_annot <- function(annot, x) {
  an <- .read_table(annot)
  key <- if (all(c("Block", "Index") %in% names(an))) c("Block", "Index") else {
    bp <- intersect(c("BP37", "BP38"), names(an))[1L]
    if (!"CHR" %in% names(an) || is.na(bp)) stop("annot needs Block and Index, or CHR and BP37 / BP38", call. = FALSE)
    c("CHR", bp)
  }
  cols <- setdiff(names(an), c("CHR", "BP37", "BP38", "A1", "A2", "Block", "Index", "SNP", "ID"))
  if (key[1L] == "CHR") an[, CHR := .chr_code(CHR)]
  for (k in key[key != "CHR"]) data.table::set(an, j = k, value = as.integer(an[[k]]))
  an <- unique(an, by = key)
  i <- an[x, on = key, which = TRUE]
  A <- if (length(cols)) as.matrix(an[i, cols, with = FALSE]) else matrix(0, nrow(x), 0)
  storage.mode(A) <- "double"
  A[is.na(A)] <- 0
  colnames(A) <- cols
  message(length(cols), " annotation columns; ", sum(is.na(i)), " SNPs without annotation (set to 0)")
  A
}

.prs_write <- function(out, output) {
  if (is.null(output)) return(invisible())
  if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
}
