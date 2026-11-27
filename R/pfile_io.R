# A .pvar.zst as plain text, for the compiled reader of the variant table (geno_file()).

# A zstd-compressed text expanded into a temporary file (R's zstd connection, R >= 4.5, else zstdcat); returns the
# path of that file, which the caller removes.
.zst_expand <- function(path) {
  tmp <- tempfile(fileext = ".pvar")
  if (!exists("zstdfile", baseenv())) {
    if (system(paste("zstdcat", shQuote(path), ">", shQuote(tmp))) != 0L) {
      unlink(tmp)
      stop("Could not expand ", path, " (zstdcat).", call. = FALSE)
    }
    return(tmp)
  }
  con <- get("zstdfile", baseenv())(path, "rb")
  out <- file(tmp, "wb")
  on.exit({ close(out); close(con) })
  repeat {
    b <- readBin(con, "raw", 64L * 1024L^2)
    if (!length(b)) break
    writeBin(b, out)
  }
  tmp
}

# f(x) for an elementwise string function f, computed once per distinct value.
.by_unique <- function(x, f) {
  u <- unique(x)
  fu <- f(u)
  if (identical(fu, u)) x else fu[match(x, u)]
}
