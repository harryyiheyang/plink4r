# Inputs shared by the summary-statistics PRS methods (viprs(), sbayesrc()): the GWAS rows matched to the eigen panel.
#
# sumstat: summarystat_tidy() / summarystat_impute() output (or a file of it) with ID, A1, A2, Z, N (r2 optional:
# rows with r2 <= 0 carry no data, e.g. imputed rows in blocks without typed SNPs, and are dropped). Rows are matched
# to the panel by ID; A1 / A2 equal to the panel alt / ref keep Z, the swapped pair flips it, other rows are dropped.
# Returns x (panel order: block, ld_index; columns ID chr bp37 bp38 A1 A2 freq Z N block ld_index, A1 = panel alt,
# freq = AF_PanUKB), blocks (ids with at least one SNP), keep (0-based eigen rows per block), n_panel (panel SNPs) and
# start (first row of each of `blocks` in the panel order block, ld_index; per-SNP panel files are indexed by it).
.prs_inputs <- function(sumstat, ref, ld_dir) {
  s <- .read_table(sumstat)
  need <- c("ID", "A1", "A2", "Z", "N")
  miss <- setdiff(need, names(s))
  if (length(miss)) stop("sumstat lacks ", paste(miss, collapse = ", "), call. = FALSE)
  s <- s[is.finite(Z) & is.finite(N) & N > 1]
  if ("r2" %in% names(s)) s <- s[is.finite(r2) & r2 > 0]
  s <- s[!ID %in% ID[duplicated(ID)], .(ID, sA1 = toupper(A1), sA2 = toupper(A2), Z = as.numeric(Z), N = as.numeric(N))]
  r <- .with_bp(.ld_panel(ref, ld_dir))
  x <- r[s, on = "ID", nomatch = NULL]
  sw <- x$sA1 == x$ref & x$sA2 == x$alt
  ok <- (x$sA1 == x$alt & x$sA2 == x$ref) | sw
  x[sw, Z := -Z]
  if (sum(!ok)) message(sum(!ok), " rows with alleles that differ from the panel dropped")
  x <- x[ok, .(ID, chr, bp37, bp38, A1 = alt, A2 = ref, freq = AF_PanUKB, Z, N, block = as.integer(block),
               ld_index = as.integer(ld_index))]
  data.table::setorder(x, block, ld_index)
  if (anyDuplicated(x, by = c("block", "ld_index"))) stop("panel rows matched more than once", call. = FALSE)
  if (!nrow(x)) stop("no sumstat row matches the panel", call. = FALSE)
  message(nrow(x), " SNPs matched to the panel (", sum(sw), " with A1 / A2 swapped)")
  bl <- x[, .(keep = list(ld_index)), by = block]
  pn <- r[, .N, keyby = .(block = as.integer(block))]
  st <- c(0, cumsum(as.numeric(pn$N)))[match(bl$block, pn$block)]
  list(x = x, blocks = bl$block, keep = bl$keep, n_panel = nrow(r), start = st)
}

.prs_write <- function(out, output) {
  if (is.null(output)) return(invisible())
  if (grepl("\\.parquet$", output)) arrow::write_parquet(out, output) else data.table::fwrite(out, output, sep = "\t")
}
