#' Project a cell_data_set object into a lower dimensional PCA (or ISI) space
#' after normalize the data
#'
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
#' 3. run \code{partition_cells} to partition cells into different graphs based
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
#' \code{pseudo_expr} to 0. This maybe useful for single-cell qPCR data, or
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
#' @param num_dim the dimensionality of the reduced space
#' @param norm_method Determines how to transform expression values prior to
#'   reducing dimensionality
#' @param residualModelFormulaStr A model formula specifying the effects to
#'   subtract from the data before clustering.
#' @param pseudo_expr amount to increase expression values before
#'   dimensionality reduction
#' @param scaling When this argument is set to TRUE (default), it will scale
#'   each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality
#'   reduction
#' @param ... additional arguments to pass to the dimensionality reduction
#'   function
#' @return an updated cell_data_set object
#' @export
preprocess_cds <- function(cds, method = c('PCA', 'none'),
                           num_dim=50,
                           norm_method = c("log", "none"),
                           residualModelFormulaStr=NULL,
                           pseudo_expr=1,
                           scaling = TRUE,
                           verbose=FALSE,
                           ...) {
  extra_arguments <- list(...)
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)
  #ensure results from RNG sensitive algorithms are the same on all calls
  set.seed(2016)
  FM <- normalize_expr_data(cds, norm_method, pseudo_expr)

  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }

  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  cds@aux_ordering_data$normalize_expr_data <- FM

  if(method == 'PCA') {
    if (verbose) message("Remove noise by PCA ...")

    irlba_res <- sparse_prcomp_irlba(Matrix::t(FM),
                                     n = min(num_dim,min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    irlba_pca_res <- irlba_res$x
    row.names(irlba_pca_res) <- colnames(cds)

  } else if(method == 'none') {
    irlba_pca_res <- Matrix::t(FM)
  } else {
    stop('unknown preprocessing method, stop!')
  }
  row.names(irlba_pca_res) <- colnames(cds)

  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose) message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(stats::as.formula(residualModelFormulaStr),
                                       data = pData(cds),
                                       drop.unused.levels = TRUE)

    fit <- limma::lmFit(t(irlba_pca_res), X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    irlba_pca_res <- Matrix::t(as.matrix(Matrix::t(irlba_pca_res)) -
                                 beta %*% Matrix::t(X.model_mat[, -1]))
  } else {
    X.model_mat <- NULL
  }

  cds@normalized_data_projection <- as.matrix(irlba_pca_res)

  cds
}


# Helper function to normalize the expression data prior to dimensionality
# reduction
normalize_expr_data <- function(cds,
                                norm_method = c("log", "none"),
                                pseudo_expr = 1) {
  FM <- exprs(cds)
  use_for_ordering <- NULL
  # If the user has selected a subset of genes for use in ordering the cells
  # via set_ordering_filter(), subset the expression matrix.
  if (is.null(fData(cds)$use_for_ordering) == FALSE &&
      nrow(subset(fData(cds), use_for_ordering == TRUE)) > 0) {
    FM <- FM[fData(cds)$use_for_ordering, ]
  }

  norm_method <- match.arg(norm_method)
  if (cds@expression_family %in% c("negbinomial", "negbinomial.size")) {

    # If we're going to be using log, and the user hasn't given us a
    # pseudocount set it to 1 by default.
    if (is.null(pseudo_expr)){
      if(norm_method == "log")
        pseudo_expr = 1
      else
        pseudo_expr = 0
    }

    check_size_factors(cds)

    if (norm_method == "log") {
      # If we are using log, normalize by size factor before log-transforming

      FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))

      if(is.null(pseudo_expr))
        pseudo_expr <- 1
      if (pseudo_expr != 1 || is_sparse_matrix(exprs(cds)) == FALSE){
        FM <- FM + pseudo_expr
        FM <- log2(FM)
      } else {
        FM@x = log2(FM@x + 1)
      }

    } else if (norm_method == "none"){
      # If we are using log, normalize by size factor before log-transforming
      FM <- Matrix::t(Matrix::t(FM)/size_factors(cds))
      FM <- FM + pseudo_expr
    }
  }else if (cds@expression_family == "binomialff") {
    if (norm_method == "none"){
      #If this is binomial data, transform expression values into TF-IDF scores.
      ncounts <- FM > 0
      ncounts[ncounts != 0] <- 1
      FM <- ncounts * log(1 + ncol(ncounts)/rowSums(ncounts))
    } else {
      stop(paste("Error: the only normalization method supported with",
                 "binomial data is 'none'"))
    }
  } else if (cds@expression_family == "Tobit") {
    FM <- FM + pseudo_expr
    if (norm_method == "none"){

    } else if (norm_method == "log"){
      FM <- log2(FM)
    } else {
      stop(paste("Error: the only normalization methods supported with",
                 "Tobit-distributed (e.g. FPKM/TPM) data are 'log'",
                 "(recommended) or 'none'"))
    }
  } else if (cds@expression_family == "uninormal") {
    if (norm_method == "none"){
      FM <- FM + pseudo_expr
    } else {
      stop(paste("Error: the only normalization method supported with",
                 "gaussian data is 'none'"))
    }
  }
  return (FM)
}

