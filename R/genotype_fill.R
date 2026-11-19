# Filling the missing genotype cells, group by group.
#
# Layers (each function computes one thing and calls only its own layer or the ones below):
#   C++ units       fill_value_cpp (the per-column value: anchor / observed mean), geno_fill_cpp (write per-column
#                   values), genotype_fill_ld_jobs_cpp (LD conditional means of several jobs on one eigen.bin, read
#                   once, written into the store)
#   geno_fill()     pointer in, pointer out: fills an open geno object group by group
#   genotype_fill() file in, pointer out: geno_open() + geno_fill(); the same layer as geno_read()

#' Fill the missing cells of an open geno object, group by group
#'
#' `fill` is a list of groups (one group may be given as it is). A group is `list(snp = , ld = )`:
#' * `snp`: data.frame with `CHR`, `BP`, `A1`, `A2` and optionally `AF` (frequency of `A1`). Its rows are
#'   matched to the columns of `X` by chromosome, position and the allele pair in either orientation (IDs are
#'   not used; `BP` must be on the build of `X`, nothing is converted). Rows not found in `X` are listed in the
#'   result and otherwise ignored.
#' * `ld`: `NULL`, or the folder of an eigen LD panel (`block<N>.eigen.bin` and the panel SNP `.parquet`), or
#'   `list(dir = folder, ref = panel SNP table)`.
#'
#' The SNPs of a group are filled from the SNPs of that group alone:
#' * anchor: `AF` when given, else the observed frequency of the column (then the fill is the observed mean);
#' * without `ld` every missing cell of a column takes the anchored value `clip((2 n af - S) / m, 0, 2)`
#'   (n rows, S observed sum, m missing cells);
#' * with `ld` every missing cell takes its conditional mean given the person's observed calls at the other
#'   SNPs of the group in the same LD block (the LD of the group's rows of that block's eigen.bin), shifted to
#'   the anchor when the column's missing rate is at least `anchor_min_miss`; a SNP of the group that is not in
#'   the panel takes the anchored value. A group that spans several LD blocks is filled block by block (the
#'   panel has no LD between blocks).
#'
#' Columns of `X` in no group are left as they are. Groups are applied in the order given; a column in two
#' groups ends with the fill of the later one. Each eigen.bin is read once, whatever the number of groups on
#' it; the files are processed one after the other, each with all the threads.
#'
#' @param X A `geno` object (hard calls for groups with `ld`).
#' @param fill A group, or a list of groups (see above).
#' @param para Named list overriding [genotype_fill_para()].
#' @return A `geno` object (new external pointer; `X` is unchanged) with `impute = "fill"` and
#'   `fill = list(snp = one row per matched SNP: group, row (of the group's table), col, snp, a1_is_counted,
#'   af (anchor used, of the counted allele), af_source ("AF" / "observed"), n_miss, method ("ld" / "anchor");
#'   dropped = rows not found in X; ld = per-job tables and timings of the LD kernel)`.
#' @export
geno_fill <- function(X, fill, para = list()) {
  if (!inherits(X, "geno")) stop("X must be a geno object.", call. = FALSE)
  P <- .fill_para(para)
  groups <- .fill_groups(fill, need_file = FALSE)
  n <- X$n; p <- X$p
  info <- geno_info_cpp(X$ptr)
  m <- info$n_miss
  vt <- list(chr = X$chr, pos = X$pos, a1 = X$a1, a2 = X$a2)
  f <- rep(NaN, p)  # NaN = keep the current fill
  rep_snp <- vector("list", length(groups)); rep_drop <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    s <- groups[[k]]$snp
    mt <- .match_alleles(vt, s$CHR, s$BP, s$A1, s$A2)
    hit <- which(!is.na(mt$j) & !duplicated(mt$j))
    cols <- mt$j[hit]
    af <- rep(NA_real_, p)
    if (!is.null(s$AF)) af[cols] <- ifelse(mt$flip[hit], 1 - s$AF[hit], s$AF[hit])
    f[cols] <- fill_value_cpp(n, info$obs_sum, m, af, numeric(0), TRUE)[cols]
    given <- is.finite(af[cols])
    af[cols[!given]] <- X$af[cols[!given]]
    groups[[k]]$cols <- cols; groups[[k]]$af <- af
    rep_snp[[k]] <- data.table::data.table(group = k, row = hit, col = cols, snp = X$snp[cols], a1_is_counted = !mt$flip[hit],
                                           af = af[cols], af_source = ifelse(given, "AF", "observed"), n_miss = m[cols],
                                           method = "anchor")
    rep_drop[[k]] <- data.table::data.table(group = rep(k, sum(is.na(mt$j))), row = which(is.na(mt$j)))
  }
  f[m == 0] <- NaN
  ptr <- geno_fill_cpp(X$ptr, f, integer(0), integer(0), numeric(0), FALSE, P$threads)
  rep_snp <- data.table::rbindlist(rep_snp)
  ld_out <- .fill_ld(X, ptr, groups, m, P)
  if (length(ld_out)) {
    done <- data.table::rbindlist(lapply(ld_out, function(o) data.table::data.table(group = o$group, col = o$cols)))
    rep_snp[done, method := "ld", on = c("group", "col")]
  }
  out <- .geno_new(ptr, .geno_cols(X), list(sample = X$sample, fid = X$fid), impute = "fill", scale = X$scale,
                   threads = X$threads, file = X$file, appended = X$appended)
  out$fill <- list(snp = rep_snp[], dropped = data.table::rbindlist(rep_drop), ld = ld_out)
  out
}

# The LD step of geno_fill(): for every panel folder, the jobs (group x LD block) are listed first, then each
# eigen.bin is read once and its jobs are solved one after the other, written straight into the store ptr.
# A job's working rows are the SNPs of its group in that block (in X and in the panel, anchor strictly inside
# (0, 1)); a job without a missing cell is not run. Returns one entry per job run.
.fill_ld <- function(X, ptr, groups, m, P) {
  has_ld <- which(!vapply(groups, function(g) is.null(g$ld), NA))
  if (!length(has_ld)) return(list())
  if (X$dosage) stop("an LD fill needs hard calls: X holds dosages.", call. = FALSE)
  .check_delta(P$delta)
  mem_gb <- .mem_gb(P$mem_gb, P$threads)
  dirs <- vapply(groups[has_ld], function(g) g$ld$dir, "")
  out <- list()
  for (d in unique(dirs)) {
    gs <- has_ld[dirs == d]
    al <- .ld_align(.ld_panel(groups[[gs[1L]]]$ld$ref, d), .geno_vartab(X))
    jobs <- list()  # by block, in group order
    for (k in gs) {
      w <- al[which((al$vi + 1L) %in% groups[[k]]$cols)]
      w[, afx := groups[[k]]$af[vi + 1L]]
      w <- w[which(is.finite(afx) & afx > 0 & afx < 1)]
      data.table::setorder(w, block, ld_index)
      w <- w[!duplicated(w[, .(block, ld_index)])]
      for (b in unique(w$block)) {
        wb <- w[which(w$block == b)]
        if (!any(m[wb$vi + 1L] > 0)) next
        jobs[[as.character(b)]] <- c(jobs[[as.character(b)]],
                                     list(list(group = k, var_idx = as.integer(wb$vi), af = wb$afx,
                                               sign = as.integer(wb$sign), rows = as.integer(wb$ld_index))))
      }
    }
    for (b in names(jobs)) {
      ef <- normalizePath(file.path(d, paste0("block", b, ".eigen.bin")))
      res <- genotype_fill_ld_jobs_cpp(ptr, ef, jobs[[b]], P$thresh, P$delta, P$anchor_min_miss, P$threads, 0L, mem_gb)
      for (q in seq_along(res)) {
        jb <- jobs[[b]][[q]]
        out[[length(out) + 1L]] <- list(group = jb$group, dir = d, block = as.integer(b), cols = jb$var_idx + 1L,
                                        snp = res[[q]]$snp, info = res[[q]]$info, time = res[[q]]$time)
      }
    }
  }
  out
}

#' Open a BED / PGEN and fill its missing cells, group by group
#'
#' One call from the file to a filled geno object: the variant table of the file is read once, the SNPs of all
#' groups are opened once ([geno_open()], missing cells at their observed mean) and filled by [geno_fill()].
#' A group is `list(bedfile = or pgenfile = , snp = , ld = )` with `snp` and `ld` as in [geno_fill()]; all
#' groups name the same file. The columns of the result are the SNPs of the groups found in the file, in the
#' order of the groups. E.g. `list(list(pgenfile = f, snp = snp1, ld = dir), list(pgenfile = f, snp = snp2))`
#' fills the SNPs of `snp1` from each other through the LD panel and the SNPs of `snp2` by their anchor.
#'
#' @param fill A group, or a list of groups.
#' @param sample Rows to open, as in `sample_vec` of [geno_open()]; `NULL` = all, in file order.
#' @param scale Whether products of the result are of the standardized matrix (see [geno_open()]).
#' @param para Named list overriding [genotype_fill_para()].
#' @return A `geno` object as from [geno_fill()].
#' @export
genotype_fill <- function(fill, sample = NULL, scale = TRUE, para = list()) {
  P <- .fill_para(para)
  groups <- .fill_groups(fill, need_file = TRUE)
  file <- unique(lapply(groups, function(g) g[c("type", "path")]))
  if (length(file) != 1L) stop("all groups must name the same file.", call. = FALSE)
  gf <- if (file[[1L]]$type == "bed") geno_file(bedfile = file[[1L]]$path) else geno_file(pgenfile = file[[1L]]$path)
  rows <- unlist(lapply(groups, function(g) .match_alleles(gf$variants, g$snp$CHR, g$snp$BP, g$snp$A1, g$snp$A2)$j))
  rows <- unique(rows[!is.na(rows)])
  if (!length(rows)) stop("no SNP of `fill` is in the file (by CHR, BP and alleles).", call. = FALSE)
  X <- geno_open(file = gf, snp_vec = rows, sample_vec = sample, impute = "mean", scale = scale, threads = P$threads)
  geno_fill(X, lapply(groups, function(g) g[c("snp", "ld")]), para)
}

# Groups of geno_fill() / genotype_fill() checked and brought to one form: list(snp = data.table(CHR, BP, A1,
# A2[, AF]), ld = NULL or list(dir, ref)[, type, path]).
.fill_groups <- function(fill, need_file) {
  if (!is.list(fill) || is.data.frame(fill) || !length(fill)) stop("fill must be a group or a list of groups.", call. = FALSE)
  if (is.data.frame(fill$snp)) fill <- list(fill)
  lapply(seq_along(fill), function(k) {
    g <- fill[[k]]
    what <- paste0("group ", k)
    if (!is.list(g) || !is.data.frame(g$snp)) stop(what, ": `snp` must be a data.frame.", call. = FALSE)
    bad <- setdiff(names(g), c("bedfile", "pgenfile", "snp", "ld"))
    if (length(bad) || is.null(names(g)) || any(!nzchar(names(g))))
      stop(what, ": a group is list(", if (need_file) "bedfile = or pgenfile = , ", "snp = , ld = ).", call. = FALSE)
    s <- data.table::as.data.table(g$snp)
    lack <- setdiff(c("CHR", "BP", "A1", "A2"), names(s))
    if (length(lack)) stop(what, ": `snp` lacks ", paste(lack, collapse = ", "), ".", call. = FALSE)
    s <- s[, intersect(c("CHR", "BP", "A1", "A2", "AF"), names(s)), with = FALSE]
    if (!is.null(s$AF)) {
      s[, AF := as.numeric(AF)]
      if (any(s$AF < 0 | s$AF > 1, na.rm = TRUE)) stop(what, ": AF must lie in [0, 1].", call. = FALSE)
    }
    ld <- g$ld
    if (is.character(ld) && length(ld) == 1L) ld <- list(dir = ld, ref = NULL)
    if (!is.null(ld)) {
      if (!is.list(ld) || !is.character(ld$dir) || length(ld$dir) != 1L || !dir.exists(ld$dir))
        stop(what, ": `ld` must be the folder of an eigen LD panel, or list(dir = , ref = ).", call. = FALSE)
      ld <- list(dir = normalizePath(ld$dir, winslash = "/"), ref = ld$ref)
    }
    out <- list(snp = s, ld = ld)
    if (need_file) {
      if (is.null(g$bedfile) == is.null(g$pgenfile)) stop(what, ": give exactly one of bedfile and pgenfile.", call. = FALSE)
      out$type <- if (is.null(g$bedfile)) "pgen" else "bed"
      out$path <- if (is.null(g$bedfile)) g$pgenfile else g$bedfile
    } else if (!is.null(g$bedfile) || !is.null(g$pgenfile)) {
      stop(what, ": X is already open; a file belongs to genotype_fill().", call. = FALSE)
    }
    out
  })
}

.fill_para <- function(para) {
  P <- genotype_fill_para()
  bad <- setdiff(names(para), names(P))
  if (length(bad)) stop("unknown para: ", paste(bad, collapse = ", "), call. = FALSE)
  P[names(para)] <- para
  P$threads <- .threads(P$threads)
  P
}

# Per-column values f written into the missing cells of X (NaN = keep the current fill); a new object, or X's
# own store when inplace (only for an object nothing else references).
.geno_fill_value <- function(X, f, impute, inplace = FALSE, threads = X$threads) {
  ptr <- geno_fill_cpp(X$ptr, f, integer(0), integer(0), numeric(0), inplace, .threads(threads))
  .geno_new(ptr, .geno_cols(X), list(sample = X$sample, fid = X$fid), impute = impute, scale = X$scale,
            threads = X$threads, file = X$file, appended = X$appended)
}

#' Defaults of [geno_fill()] and [genotype_fill()]
#'
#' delta 0.01: ridge of the LD conditional mean; anchor_min_miss 0.01: an LD-filled column is shifted to its
#' anchor only at or above this missing rate (< 0 = never); thresh 0: all stored eigenpairs of the panel (else
#' the share of the eigenvalue sum to keep); threads; mem_gb: memory budget of the people batches of the LD
#' step (0 takes threads and memory from [plink4r_resources()]).
#' @export
genotype_fill_para <- function() {
  list(delta = 0.01, anchor_min_miss = 0.01, thresh = 0, threads = 4L, mem_gb = 2)
}
