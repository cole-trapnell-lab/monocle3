#' Preprocess a cds to prepare for trajectory inference
#'
#' @description Most analyses (including trajectory inference, and clustering)
#' in Monocle3, require various normalization and preprocessing steps.
#' \code{preprocess_cds} executes and stores these preprocessing steps.
#'
#' Specifically, depending on the options selected, \code{preprocess_cds} first
#' normalizes the data by log and size factor to address depth differences, or
#' by size factor only. Next, \code{preprocess_cds} calculates a lower
#' dimensional space that will be used as the input for further dimensionality
#' reduction like tSNE and UMAP.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param method a string specifying the initial dimension method to use,
#'   currently either PCA or LSI. For LSI (latent semantic indexing), it
#'   converts the (sparse) expression matrix into tf-idf matrix and then
#'   performs SVD to decompose the gene expression / cells into certain
#'   modules / topics. Default is "PCA".
#' @param pca_method a string specifying the PCA algorithm to use, which can be
#'  irlba, rsvd, or svdr. irlba is the implicitly restarted Lanczos
#' bidiagonalization algorithm, rsvd is the randomized SVD method in the
#' rsvd package, and svdr is the randomized SVD method in the irlba package.
#' Default is "irlba".
#' @param num_dim the dimensionality of the reduced space.
#' @param norm_method Determines how to transform expression values prior to
#'   reducing dimensionality. Options are "log", "size_only", and "none".
#'   Default is "log". Users should only use "none" if they are confident that
#'   their data is already normalized.
#' @param use_genes NULL or a list of gene IDs. If a list of gene IDs, only
#'   this subset of genes is used for dimensionality reduction. Default is
#'   NULL.
#' @param residual_model_formula_str NULL or a string model formula specifying
#'   any effects to subtract from the data before dimensionality reduction.
#'   Uses a linear model to subtract effects. For non-linear effects, use
#'   alignment_group. Default is NULL.
#' @param alignment_group String specifying a column of colData to use for
#'  aligning groups of cells. The column specified must be a factor.
#'  Alignment can be used to subtract batch effects in a non-linear way.
#'  For correcting continuous effects, use residual_model_formula_str.
#'  Default is NULL.
#' @param pseudo_count NULL or the amount to increase expression values before
#'   normalization and dimensionality reduction. If NULL (default), a
#'   pseudo_count of 1 is added for log normalization and 0 is added for size
#'   factor only normalization.
#' @param scaling When this argument is set to TRUE (default), it will scale
#'   each gene before running trajectory reconstruction. Relevant for
#'   method = PCA only.
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param ... additional arguments to pass to limma::lmFit if
#'   residual_model_formula is not NULL
#' @return an updated cell_data_set object
#' @export
preprocess_cds <- function(cds, method = c('PCA', "LSI"),
                           pca_method = c( 'irlba', 'rsvd', 'svdr'),
                           num_dim=50,
                           norm_method = c("log", "size_only", "none"),
                           use_genes = NULL,
                           residual_model_formula_str=NULL,
                           alignment_group=NULL,
                           pseudo_count=NULL,
                           scaling = TRUE,
                           verbose=FALSE,
                           ...) {

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'PCA' or 'LSI'")
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(pca_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "pca_method must be one of 'irlba' or 'rsvd'")
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "norm_method must be one of 'log', 'size_only' or 'none'")
  assertthat::assert_that(assertthat::is.count(num_dim))
  if(!is.null(use_genes)) {
    assertthat::assert_that(is.character(use_genes))
    assertthat::assert_that(all(use_genes %in% row.names(rowData(cds))),
                            msg = paste("use_genes must be NULL, or all must",
                            "be present in the row.names of rowData(cds)"))
  }
  assertthat::assert_that(!is.null(size_factors(cds)),
             msg = paste("You must call estimate_size_factors before calling",
                         "preprocess_cds."))
  assertthat::assert_that(sum(is.na(size_factors(cds))) == 0,
                          msg = paste("One or more cells has a size factor of",
                                      "NA."))

  method <- match.arg(method)
  pca_method <- match.arg(pca_method)
  norm_method <- match.arg(norm_method)

  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)
  FM <- normalize_expr_data(cds, norm_method, pseudo_count)

  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }

  if (!is.null(use_genes)) {
    FM <- FM[use_genes, ]
  }

  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  if(method == 'PCA') {
    if (verbose) message("Remove noise by PCA ...")

  irlba_res <- switch(pca_method,
      irlba = sparse_prcomp_irlba(Matrix::t(FM),
                                  n = min(num_dim,min(dim(FM)) - 1),
                                  center = scaling, scale. = scaling),
      rsvd = sparse_prcomp_rsvd(Matrix::t(FM),
                                  n = min(num_dim,min(dim(FM)) - 1),
                                  center = scaling, scale. = scaling)
      svdr = sparse_prcomp_svdr(Matrix::t(FM),
                                  n = min(num_dim,min(dim(FM)) - 1),
                                  center = scaling, scale. = scaling)
    )
    preproc_res <- irlba_res$x
    row.names(preproc_res) <- colnames(cds)

    irlba_rotation <- irlba_res$rotation
    row.names(irlba_rotation) <- rownames(FM)
    cds@preprocess_aux$gene_loadings <- irlba_rotation
    cds@preprocess_aux$prop_var_expl <- irlba_res$sdev^2 / sum(irlba_res$sdev^2)

  } else if(method == "LSI") {

    preproc_res <- tfidf(FM)
    irlba_res <- irlba::irlba(Matrix::t(preproc_res),
                              nv = min(num_dim,min(dim(FM)) - 1))

    preproc_res <- irlba_res$u %*% diag(irlba_res$d)
    row.names(preproc_res) <- colnames(cds)

    irlba_rotation = irlba_res$v
    row.names(irlba_rotation) = rownames(FM)
    cds@preprocess_aux$gene_loadings = irlba_rotation

  }

  row.names(preproc_res) <- colnames(cds)

  reducedDims(cds)[[method]] <- as.matrix(preproc_res)
  cds@preprocess_aux$beta = NULL

  cds
}


# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(cds,
                                norm_method = c("log", "size_only", "none"),
                                pseudo_count = NULL) {
  norm_method <- match.arg(norm_method)

  FM <- SingleCellExperiment::counts(cds)

  # If we're going to be using log, and the user hasn't given us a
  # pseudocount set it to 1 by default.
  if (is.null(pseudo_count)){
    if(norm_method == "log")
      pseudo_count <- 1
    else
      pseudo_count <- 0
  }

  if (norm_method == "log") {
    # If we are using log, normalize by size factor before log-transforming

    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))

    if (pseudo_count != 1 || is_sparse_matrix(SingleCellExperiment::counts(cds)) == FALSE){
      FM <- FM + pseudo_count
      FM <- log2(FM)
    } else {
      FM@x = log2(FM@x + 1)
    }

  } else if (norm_method == "size_only") {
    FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
    FM <- FM + pseudo_count
  }
  return (FM)
}

# Andrew's tfidf
tfidf <- function(count_matrix, frequencies=TRUE, log_scale_tf=TRUE,
                  scale_factor=100000, block_size=2000e6) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf <- Matrix::t(Matrix::t(count_matrix) / Matrix::colSums(count_matrix))
  } else {
    # "raw count" method
    tf <- count_matrix
  }

  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }

  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(count_matrix) / Matrix::rowSums(count_matrix > 0))

  # Try to just to the multiplication and fall back on delayed array
  # TODO hopefully this actually falls back and not get jobs killed in SGE
  tf_idf_counts = tryCatch({
    tf_idf_counts = tf * idf
    tf_idf_counts
  }, error = function(e) {
    print(paste("TF*IDF multiplication too large for in-memory, falling back",
                "on DelayedArray."))
    options(DelayedArray.block.size=block_size)
    DelayedArray:::set_verbose_block_processing(TRUE)

    tf = DelayedArray::DelayedArray(tf)
    idf = as.matrix(idf)

    tf_idf_counts = tf * idf
    tf_idf_counts
  })

  rownames(tf_idf_counts) = rownames(count_matrix)
  colnames(tf_idf_counts) = colnames(count_matrix)
  tf_idf_counts = methods::as(tf_idf_counts, "sparseMatrix")
  return(tf_idf_counts)
}
