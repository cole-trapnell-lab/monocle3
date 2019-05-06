#' Project a cell_data_set object into a lower dimensional PCA (or ISI) space
#' after normalize the data
#'
#'
#'### Data will still be size_factor_normalized! ###
#' @description For most analysis (including trajectory inference, clustering)
#' in Monocle 3, it requires us to to start from a low dimensional PCA space.
#' \code{preprocess_cds} will be used to first project a cell_data_set object
#' into a lower dimensional PCA space before we apply clustering with community
#' detection algorithm or other non-linear dimension reduction method, for
#' example UMAP, tSNE, etc.  While tSNE is especially suitable for visualizing
#' clustering results, comparing to UMAP, the global distance in tSNE space is
#' not meaningful. UMAP can either be used for visualizing clustering result or
#' as a general non-linear dimension reduction method. SimplePPT, DDRTree and
#' L1-graph are two complementary trajectory inference method where the first
#' one is very great at learning a tree structure but the later is general and
#' can learn any arbitrary graph structure. Both methods can be applied to the
#' UMAP space.
#'
#' @details
#' In Monocle3, we overhauled the code from Monocle2 so that a standard
#' Monocle3 workflow works as following:
#' 1. run \code{preprocess_cds} to project a cell_data_set object into a lower
#' dimensional PCA space after normalizing the data
#' 2. run \code{reduce_dimension} to further project the PCA space into much
#' lower dimension space with non-linear dimension reduction techniques,
#' including tSNE, UMAP.
#' 3. run \code{cluster_cells} to partition cells into different graphs based
#' on a similar approach proposed by Alex Wolf and colleagues. We then
#' reconstruct the trajectory in each partition with the \code{learn_graph}
#' function.
#' 4. run \code{learn_graph} to reconstruct developmental trajectory with
#' reversed graph embedding algorithms. In monocle 3, we enabled the capability
#' to learn multiple disjointed trajectory with either tree or loop structure,
#' etc.
#'
#' Prior to reducing the dimensionality of the data, it usually helps to
#' normalize it so that highly expressed or highly variable genes don't
#' dominate the computation. \code{reduce_dimension()} automatically transforms
#' the data in one of several ways depending on the \code{expression_family} of
#' the cell_data_set object. If the expression_family is \code{"negbinomial"}
#' or \code{"negbinomial.size"}, the data are variance-stabilized. If the
#' expression_family is \code{"Tobit"}, the data are adjusted by adding a
#' pseudocount (of 1 by default) and then log-transformed. If you don't want
#' any transformation at all, set \code{norm_method} to "none" and
#' \code{pseudo_count} to 0. This maybe useful for single-cell qPCR data, or
#' data you've already transformed yourself in some way.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param method the initial dimension method to use, current either PCA or
#'   LSI. For LSI (latent semantic indexing), it converts the (sparse)
#'   expression matrix into tf-idf (term-frequency-inverse document frequency
#'   which increases proportionally to the gene expression value appears in the
#'   cell and is offset by the frequency of the gene in the entire dataset,
#'   which helps to adjust for the fact that some gene appear more frequently
#'   across cell in general.) matrix and then performs SVD to decompose the
#'   gene expression / cells into certain modules / topics. This method can be
#'   used to find associated gene modules and cell clusters at the same time.
#'   It removes noise in the data and thus makes the UMAP result even better.
#' @param num_dim the dimensionality of the reduced space. Ignored if
#'   method = "none".
#' @param norm_method Determines how to transform expression values prior to
#'   reducing dimensionality
#' @param residual_model_formula_str A model formula specifying the effects to
#'   subtract from the data before clustering.
#' @param pseudo_count amount to increase expression values before
#'   dimensionality reduction. If NULL (default), pseudo_count of 1 is added
#'   for log normalization and 0 is added for no normalization.
#' @param scaling When this argument is set to TRUE (default), it will scale
#'   each gene before running trajectory reconstruction. Relevant for
#'   method = PCA only.
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param ... additional arguments to pass to the dimensionality reduction
#'   function
#' @return an updated cell_data_set object
#' @export
preprocess_cds <- function(cds, method = c('PCA', "LSI"),
                           num_dim=50,
                           norm_method = c("log", "size_only"),
                           use_genes = NULL,
                           residual_model_formula_str=NULL,
                           pseudo_count=NULL,
                           scaling = TRUE,
                           verbose=FALSE,
                           ...) {
  extra_arguments <- list(...)
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "method must be one of 'PCA' or 'LSI'")
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(norm_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "norm_method must be one of 'log' or 'size_only'")
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
  norm_method <- match.arg(norm_method)

  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)
  FM <- normalize_expr_data(cds, norm_method, pseudo_count)

  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }

  # If the user has selected a subset of genes for use in ordering the cells
  # via set_ordering_filter(), subset the expression matrix.
  # TO DO
  if (!is.null(use_genes)) {
    FM <- FM[use_genes, ]
  }

  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  if(method == 'PCA') {
    if (verbose) message("Remove noise by PCA ...")

    irlba_res <- sparse_prcomp_irlba(Matrix::t(FM),
                                     n = min(num_dim,min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    preproc_res <- irlba_res$x
    row.names(preproc_res) <- colnames(cds)

  } else if(method == "LSI") {
    preproc_res <- tfidf(FM)
    do_svd <- function(tf_idf_counts, dims=50) {
      pca.results = irlba::irlba(Matrix::t(tf_idf_counts), nv=dims)
      final_result = pca.results$u %*% diag(pca.results$d)
      rownames(final_result) = colnames(tf_idf_counts)
      colnames(final_result) = paste0('C_', 1:dims)
      return(final_result)
    }
    preproc_res <- do_svd(preproc_res, num_dim)
  }

  row.names(preproc_res) <- colnames(cds)

  if (!is.null(residual_model_formula_str)) {
    if (verbose) message("Removing batch effects")
    X.model_mat <- Matrix::sparse.model.matrix(stats::as.formula(residual_model_formula_str),
                                       data = colData(cds),
                                       drop.unused.levels = TRUE)

    fit <- limma::lmFit(Matrix::t(preproc_res), X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    preproc_res <- Matrix::t(as.matrix(Matrix::t(preproc_res)) -
                                 beta %*% Matrix::t(X.model_mat[, -1]))
  }

  reducedDims(cds)[[method]] <- as.matrix(preproc_res)

  cds
}


# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(cds,
                                norm_method = c("log", "size_only"),
                                pseudo_count = NULL) {
  norm_method <- match.arg(norm_method)

  FM <- counts(cds)

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

    if (pseudo_count != 1 || is_sparse_matrix(counts(cds)) == FALSE){
      FM <- FM + pseudo_count
      FM <- log2(FM)
    } else {
      FM@x = log2(FM@x + 1)
    }

  } else if (norm_method == "size_only"){
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
  idf = log(1 + ncol(count_matrix) / Matrix::rowSums(count_matrix))

  # Try to just to the multiplication and fall back on delayed array
  # TODO hopefully this actually falls back and not get jobs killed in SGE
  tf_idf_counts = tryCatch({
    tf_idf_counts = tf * idf
    tf_idf_counts
  }, error = function(e) {
    print("TF*IDF multiplication too large for in-memory, falling back on DelayedArray.")
    options(DelayedArray.block.size=block_size)
    DelayedArray:::set_verbose_block_processing(TRUE)

    tf = DelayedArray(tf)
    idf = as.matrix(idf)

    tf_idf_counts = tf * idf
    tf_idf_counts
  })

  rownames(tf_idf_counts) = rownames(count_matrix)
  colnames(tf_idf_counts) = colnames(count_matrix)
  tf_idf_counts = as(tf_idf_counts, "sparseMatrix")
  return(tf_idf_counts)
}
