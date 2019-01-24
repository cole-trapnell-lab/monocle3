#' Project a cell_data_set object into a lower dimensional PCA (or ISI) space after normalize the data
#'
#' @description For most analysis (including trajectory inference, clustering) in Monocle 3, it requires us to to start from a
#' low dimensional PCA space. preprocessCDS will be used to first project a cell_data_set object into a lower dimensional PCA space
#' before we apply clustering with community detection algorithm or other non-linear dimension reduction method, for example
#' UMAP, tSNE, DDRTree, L1-graph, etc.  While tSNE is especially suitable for visualizing clustering results, comparing
#' to UMAP, the global distance in tSNE space is not meaningful. UMAP can either be used for visualizing clustering result or as a general
#' non-linear dimension reduction method.
#' SimplePPT, DDRTree and L1-graph are two complementary trajectory inference method where the first one is very great at learning a tree structure
#' but the later is general and can learn any arbitrary graph structure. Both methods can be applied to the UMAP space.
#'
#' @details
#' In Monocle 3, we overhauled the code from Monocle2 so that a standard Monocle 3 workingflow works as following:
#' 1. run \code{preprocessCDS} to project a cell_data_set object into a lower dimensional PCA space after
#' normalize the data
#' 2. run \code{reduceDimension} to further project the PCA space into much lower dimension space with non-linear
#' dimension reduction techniques, including tSNE, UMAP.
#' 3. run \code{smoothEmbedding} (optional) to smooth noisy embedding from 2 to facilitate visualization and learning
#' of the graph structure.
#' 4. run \code{partitionCells} to partition cells into different graphs based on a similar approach proposed by Alex Wolf and colleagues.
#' We then reconstruct the trajectory in each partition with the \code{learnGraph} function.
#' 5. run \code{learnGraph} to reconstruct developmental trajectory with reversed graph embedding algorithms. In monocle 3, we enabled the
#' the capability to learn multiple disjointed trajectory with either tree or loop structure, etc.
#'
#' Prior to reducing the dimensionality of the data, it usually helps
#' to normalize it so that highly expressed or highly variable genes don't
#' dominate the computation. \code{reduceDimension()} automatically transforms
#' the data in one of several ways depending on the \code{expressionFamily} of
#' the cell_data_set object. If the expressionFamily is \code{"negbinomial"} or \code{"negbinomial.siz"e}, the
#' data are variance-stabilized. If the expressionFamily is \code{"tobit"}, the data
#' are adjusted by adding a pseudocount (of 1 by default) and then log-transformed.
#' If you don't want any transformation at all, set norm_method to "none" and
#' pseudo_expr to 0. This maybe useful for single-cell qPCR data, or data you've
#' already transformed yourself in some way.

#' @param cds the cell_data_set upon which to perform this operation
#' @param method the initial dimension method to use, current either PCA or LSI. For LSI (latent semantic indexing),
#' it converts the (sparse) expression matrix into tf-idf (term-frequency-inverse document frequency
#' which increases proportionally to the gene expression value appears in the cell and is offset by the frequency of
#' the gene in the entire dataset, which helps to adjust for the fact that some gene appear more frequently across cell in general.) matrix and then performs a
#' SVD to decompose the gene expression / cells into certain modules / topics. This method can be used to find associated gene modules
#  and cell clusters at the same time. It removes noise in the data and thus makes the UMAP result even better.
#' @param num_dim the dimensionality of the reduced space
#' @param norm_method Determines how to transform expression values prior to reducing dimensionality
#' @param residualModelFormulaStr A model formula specifying the effects to subtract from the data before clustering.
#' @param pseudo_expr amount to increase expression values before dimensionality reduction
#' @param relative_expr When this argument is set to TRUE (default), we intend to convert the expression into a relative expression.
#' @param scaling When this argument is set to TRUE (default), it will scale each gene before running trajectory reconstruction.
#' @param verbose Whether to emit verbose output during dimensionality reduction
#' @param ... additional arguments to pass to the dimensionality reduction function
#' @return an updated cell_data_set object
#' @export
preprocess_cds <- function(cds, method = c('PCA', 'none'), #, 'LSI' , 'NMF'
                          num_dim=50,
                          norm_method = c("log", "vstExprs", "none"),
                          residualModelFormulaStr=NULL,
                          pseudo_expr=1,
                          relative_expr=TRUE,
                          scaling = TRUE,
                          verbose=FALSE,
                          ...) {
  extra_arguments <- list(...)
  method <- match.arg(method)
  norm_method <- match.arg(norm_method)
  set.seed(2016) #ensure results from RNG sensitive algorithms are the same on all calls

  FM <- normalize_expr_data(cds, norm_method, pseudo_expr, relative_expr)

  if (nrow(FM) == 0) {
    stop("Error: all rows have standard deviation zero")
  }

  fm_rowsums = Matrix::rowSums(FM)
  FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

  cds@auxOrderingData$normalize_expr_data <- FM

  if(method == 'PCA') {
    if (verbose)
      message("Remove noise by PCA ...")

    irlba_res <- sparse_prcomp_irlba(t(FM), n = min(num_dim, min(dim(FM)) - 1),
                                     center = scaling, scale. = scaling)
    irlba_pca_res <- irlba_res$x
    row.names(irlba_pca_res) <- colnames(cds)

  } else if(method == 'none') {
    irlba_pca_res <- t(FM)
  } else {
    stop('unknown preprocessing method, stop!')
  }
  row.names(irlba_pca_res) <- colnames(cds)

  if (is.null(residualModelFormulaStr) == FALSE) {
    if (verbose)
      message("Removing batch effects")
    X.model_mat <- sparse.model.matrix(as.formula(residualModelFormulaStr),
                                       data = pData(cds), drop.unused.levels = TRUE)

    fit <- limma::lmFit(t(irlba_pca_res), X.model_mat, ...)
    beta <- fit$coefficients[, -1, drop = FALSE]
    beta[is.na(beta)] <- 0
    irlba_pca_res <- t(as.matrix(t(irlba_pca_res)) - beta %*% t(X.model_mat[, -1]))
  }else{
    X.model_mat <- NULL
  }

  cds@normalized_data_projection <- as.matrix(irlba_pca_res)

  cds
}
