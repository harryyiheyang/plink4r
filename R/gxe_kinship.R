# The kinship route of gxe() (and its Cox model): the null model is the one of gwas() with the same arguments
# (.reml_null: the same model, the same kinship parameter, fitted once), and the variants are streamed against that
# null-model object (gxe_kin_stream_cpp). The exposures and interaction covariates are columns of its design; with
# V the covariance of the (working) response under the null model, P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1 and
# Z = [g, g e_1, ..., g e_{S-1}], a variant gets c = Z' P y and A = Z' P Z (S x S): the association test of g is
# c_1, A_11 (that of gwas()), and beta = A^-1 c with covariance A^-1 gives the main effect, the interactions and
# their Wald tests (generalised least squares; one step from the null model for the logistic and Cox models).
.gxe_kinship <- function(formula, data, exposure, out, int_covar, gf, kinship, family, id_col, center,
                         threads, mem_gb, block, para) {
  t0 <- proc.time()[["elapsed"]]
  P <- .reml_para(para)
  threads <- .threads(threads); mem_gb <- .mem_gb(mem_gb, threads)
  D <- .mm_design(formula, data, family, id_col, gf$samples$iid)
  lab <- gsub("`", "", D$labels)
  nt <- setdiff(gsub("`", "", c(exposure, int_covar)), lab)
  if (length(nt))
    stop("Not a term of the formula: ", paste(nt, collapse = ", "), ". Its terms: ", paste(lab, collapse = ", "), ".",
         call. = FALSE)
  cols_of <- function(v) which(D$assign %in% match(gsub("`", "", v), lab))
  ie <- cols_of(exposure); ii <- cols_of(int_covar)
  if (length(intersect(ie, ii))) stop("A term is given both as exposure and as int_covar.", call. = FALSE)
  E <- D$X[, c(ie, ii), drop = FALSE]
  if (any(apply(E, 2L, stats::var) == 0))
    stop("An exposure or interaction covariate is constant among the analysed samples.", call. = FALSE)
  # centring changes the parametrisation only: the null model is not touched
  if (center) E <- sweep(E, 2L, colMeans(E))
  terms <- colnames(E); ne <- length(ie)
  nl <- .reml_null(D, kinship, family, P, proj = TRUE)
  t1 <- proc.time()[["elapsed"]]
  res <- gxe_kin_stream_cpp(gf$ptr, out, as.integer(D$srow - 1L), nl$null, E, ne,
                            threads, mem_gb, as.integer(block))
  t2 <- proc.time()[["elapsed"]]
  invisible(c(list(out = out, variants = res$variants, bytes = res$bytes),
              .kin_value(nl, family, res),
              list(exposure_terms = terms[seq_len(ne)], int_covar_terms = terms[-seq_len(ne)],
                   time = c(null = t1 - t0, stream = t2 - t1, total = t2 - t0))))
}
