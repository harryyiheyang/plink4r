# Filling the missing genotype cells, group by group.
#
# Layers (each function computes one thing and calls only its own layer or the ones below):
#   C++ units       fill_value_cpp (the per-column value: anchor / observed mean), geno_fill_cpp (write per-column
#                   values), genotype_fill_ld_jobs_cpp (LD conditional means of several jobs on one eigen.bin, read
#                   once, written into the stores)
#   geno_fill()     pointers in, pointers out: fills open geno objects, one group per object
#   genotype_fill() file in, pointer out: geno_open() + geno_fill(); the same layer as geno_read()

#' Fill the missing cells of open geno objects, one group per object
#'
#' `X` is a geno object and `fill` its group, or `X` is a list of geno objects and `fill` the list of their
#' groups (group k fills object k). A group is `list(snp = , ld = )`:
#' * `snp`: data.frame with `CHR`, `BP`, `A1`, `A2` and optionally `AF` (frequency of `A1`). Its rows are
#'   matched to the columns of the object by chromosome, position and the allele pair in either orientation
#'   (IDs are not used; `BP` must be on the build of the object, nothing is converted). Rows not found are
#'   listed in the result and otherwise ignored.
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
#'   panel has no LD between blocks); the block of every SNP is in the result.
#'
#' Columns of an object that are not in its group are left as they are. The jobs (object x LD block) are
#' listed first; each eigen.bin is then read once, whatever the number of objects on it, and the files are
#' processed one after the other, each with all the threads.
#'
#' @param X A `geno` object, or a list of them (hard calls for groups with `ld`).
#' @param fill The group of `X`, or the list of the groups of the objects of `X`.
#' @param para Named list overriding [genotype_fill_para()].
#' @return A `geno` object, or a list of them named as `fill` (`X1`, `X2`, ... when unnamed): new external
#'   pointers, the inputs are unchanged. Each has `impute = "fill"` and `fill = list(snp = one row per
#'   matched SNP: row (of the group's table), col, snp, a1_is_counted, af (anchor used, of the counted
#'   allele), af_source ("AF" / "observed"), n_miss, method ("ld" / "anchor"), block (LD block of an "ld"
#'   SNP); dropped = rows of the table not found in the object; ld = per-block tables and timings of the
#'   LD kernel)`.
#' @export
geno_fill <- function(X, fill, para = list()) {
  one <- inherits(X, "geno")
  if (one) X <- list(X)
  if (!is.list(X) || !length(X) || !all(vapply(X, inherits, NA, "geno"))) stop("X must be a geno object or a list of them.", call. = FALSE)
  P <- .fill_para(para)
  groups <- .fill_groups(fill, need_file = FALSE)
  if (length(groups) != length(X)) stop("one group per geno object: ", length(X), " object(s), ", length(groups), " group(s).", call. = FALSE)
  K <- length(X)
  ptr <- vector("list", K); rep_snp <- vector("list", K); rep_drop <- vector("list", K)
  for (k in seq_len(K)) {
    Xk <- X[[k]]; s <- groups[[k]]$snp
    info <- geno_info_cpp(Xk$ptr)
    m <- info$n_miss
    mt <- .match_alleles(list(chr = Xk$chr, pos = Xk$pos, a1 = Xk$a1, a2 = Xk$a2), s$CHR, s$BP, s$A1, s$A2)
    hit <- which(!is.na(mt$j) & !duplicated(mt$j))
    cols <- mt$j[hit]
    af <- rep(NA_real_, Xk$p)
    if (!is.null(s$AF)) af[cols] <- ifelse(mt$flip[hit], 1 - s$AF[hit], s$AF[hit])
    f <- rep(NaN, Xk$p)  # NaN = keep the current fill
    f[cols] <- fill_value_cpp(Xk$n, info$obs_sum, m, af, numeric(0), TRUE)[cols]
    f[m == 0] <- NaN
    given <- is.finite(af[cols])
    af[cols[!given]] <- Xk$af[cols[!given]]
    groups[[k]]$cols <- cols; groups[[k]]$af <- af; groups[[k]]$n_miss <- m
    ptr[[k]] <- geno_fill_cpp(Xk$ptr, f, FALSE, P$threads)
    rep_snp[[k]] <- data.table::data.table(row = hit, col = cols, snp = Xk$snp[cols], a1_is_counted = !mt$flip[hit],
                                           af = af[cols], af_source = ifelse(given, "AF", "observed"), n_miss = m[cols],
                                           method = "anchor", block = NA_integer_)
    rep_drop[[k]] <- data.table::data.table(row = which(is.na(mt$j)))
  }
  ld_out <- .fill_ld(X, ptr, groups, P)
  out <- lapply(seq_len(K), function(k) {
    Xk <- X[[k]]
    mine <- ld_out[vapply(ld_out, function(o) o$object == k, NA)]
    for (o in mine) rep_snp[[k]][match(o$cols, col), `:=`(method = "ld", block = o$block)]
    Y <- .geno_new(ptr[[k]], .geno_cols(Xk), list(sample = Xk$sample, fid = Xk$fid), impute = "fill", scale = Xk$scale,
                   threads = Xk$threads, file = Xk$file, appended = Xk$appended)
    Y$fill <- list(snp = rep_snp[[k]][], dropped = rep_drop[[k]], ld = mine)
    Y
  })
  if (one) return(out[[1L]])
  names(out) <- if (is.null(names(groups))) paste0("X", seq_len(K)) else names(groups)
  out
}

# The LD step of geno_fill(): for every panel folder, the jobs (object x LD block) are listed first, then each
# eigen.bin is read once and its jobs are solved one after the other, written straight into the stores ptr.
# A job's working rows are the SNPs of the object's group in that block (in the object and in the panel,
# anchor strictly inside (0, 1)); a job without a missing cell is not listed, so a block without one is never
# read. Returns one entry per job run.
.fill_ld <- function(X, ptr, groups, P) {
  # an object whose group has no missing cell is skipped here, before any panel table or eigen.bin is read
  has_ld <- which(vapply(groups, function(g) !is.null(g$ld) && any(g$n_miss[g$cols] > 0), NA))
  if (!length(has_ld)) return(list())
  if (any(vapply(X[has_ld], function(x) x$dosage, NA))) stop("an LD fill needs hard calls: the object holds dosages.", call. = FALSE)
  .check_delta(P$delta)
  mem_gb <- .mem_gb(P$mem_gb, P$threads)
  dirs <- vapply(groups[has_ld], function(g) g$ld$dir, "")
  out <- list()
  for (d in unique(dirs)) {
    gs <- has_ld[dirs == d]
    pan <- .ld_panel(groups[[gs[1L]]]$ld$ref, d)
    jobs <- list()  # by block, in object order
    for (k in gs) {
      al <- .ld_align(pan, .geno_vartab(X[[k]]))
      w <- al[which((al$vi + 1L) %in% groups[[k]]$cols)]
      w[, afx := groups[[k]]$af[vi + 1L]]
      w <- w[which(is.finite(afx) & afx > 0 & afx < 1)]
      data.table::setorder(w, block, ld_index)
      w <- w[!duplicated(w[, .(block, ld_index)])]
      for (b in unique(w$block)) {
        wb <- w[which(w$block == b)]
        if (!any(groups[[k]]$n_miss[wb$vi + 1L] > 0)) next
        jobs[[as.character(b)]] <- c(jobs[[as.character(b)]],
                                     list(list(object = k, ptr = ptr[[k]], var_idx = as.integer(wb$vi), af = wb$afx,
                                               sign = as.integer(wb$sign), rows = as.integer(wb$ld_index))))
      }
    }
    for (b in names(jobs)) {
      ef <- normalizePath(file.path(d, paste0("block", b, ".eigen.bin")))
      res <- genotype_fill_ld_jobs_cpp(ef, jobs[[b]], P$thresh, P$delta, P$anchor_min_miss, P$threads, 0L, mem_gb)
      for (q in seq_along(res)) {
        jb <- jobs[[b]][[q]]
        out[[length(out) + 1L]] <- list(object = jb$object, dir = d, block = as.integer(b), cols = jb$var_idx + 1L,
                                        snp = res[[q]]$snp, info = res[[q]]$info, time = res[[q]]$time)
      }
    }
  }
  out
}

#' Open BED / PGEN files and fill their missing cells: one geno object per group
#'
#' A group is `list(bedfile = or pgenfile = , snp = , ld = )` with `snp` and `ld` as in [geno_fill()]. Each
#' group is opened as its own geno object ([geno_open()]: the SNPs of its table found in its file, in file
#' order, missing cells at their observed mean) and filled from its own SNPs by [geno_fill()]; a file named
#' by several groups has its variant and sample tables read once. E.g.
#' `list(X1 = list(pgenfile = f, snp = snp1, ld = dir), X2 = list(pgenfile = f, snp = snp2))` returns `X1`,
#' the SNPs of `snp1` filled from each other through the LD panel, and `X2`, the SNPs of `snp2` at their
#' anchored value; [geno_merge()] puts them side by side.
#'
#' @param fill A group, or a list of groups.
#' @param sample Rows to open, as in `sample_vec` of [geno_open()]; `NULL` = all, in file order.
#' @param scale Whether products of the results are of the standardized matrix (see [geno_open()]).
#' @param para Named list overriding [genotype_fill_para()].
#' @return As [geno_fill()]: a `geno` object for one group, else a list of them, one per group.
#' @export
genotype_fill <- function(fill, sample = NULL, scale = TRUE, para = list()) {
  P <- .fill_para(para)
  one <- is.list(fill) && is.data.frame(fill$snp)
  groups <- .fill_groups(fill, need_file = TRUE)
  key <- vapply(groups, function(g) paste(g$type, g$path), "")
  gfs <- lapply(split(seq_along(groups), key), function(i) {
    g <- groups[[i[1L]]]
    if (g$type == "bed") geno_file(bedfile = g$path) else geno_file(pgenfile = g$path)
  })
  X <- lapply(seq_along(groups), function(k) {
    g <- groups[[k]]; gf <- gfs[[key[k]]]
    rows <- .match_alleles(gf$variants, g$snp$CHR, g$snp$BP, g$snp$A1, g$snp$A2)$j
    rows <- sort(unique(rows[!is.na(rows)]))
    if (!length(rows)) stop("group ", k, ": no SNP of `snp` is in the file (by CHR, BP and alleles).", call. = FALSE)
    geno_open(file = gf, snp_vec = rows, sample_vec = sample, impute = "mean", scale = scale, threads = P$threads)
  })
  gl <- lapply(groups, function(g) list(snp = g$snp, ld = g$ld))
  if (one) geno_fill(X[[1L]], gl[[1L]], para) else geno_fill(X, gl, para)
}

# Groups of geno_fill() / genotype_fill() checked and brought to one form: list(snp = data.table(CHR, BP, A1,
# A2[, AF]), ld = NULL or list(dir, ref)[, type, path]).
.fill_groups <- function(fill, need_file) {
  if (!is.list(fill) || is.data.frame(fill) || !length(fill)) stop("fill must be a group or a list of groups.", call. = FALSE)
  if (is.data.frame(fill$snp)) fill <- list(fill)
  out <- lapply(seq_along(fill), function(k) {
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
  names(out) <- names(fill)
  out
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
  ptr <- geno_fill_cpp(X$ptr, f, inplace, .threads(threads))
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
