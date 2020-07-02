#' Cluster genes into modules that are co-expressed across cells.
#'
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param reduction_method The dimensionality reduction method used to generate the lower dimensional space in which genes will be clustered. Currently only UMAP is supported.
#' @param max_components The number of dimensions in which to cluster genes into modules.
#' @param umap.metric Metric used by UMAP for measuring similarity between genes .
#' @param umap.min_dist Minimum distance parameter passed to UMAP.
#' @param umap.n_neighbors Number of nearest neighbors used by UMAP.
#' @param umap.fast_sgd Whether to allow UMAP to perform fast stochastic gradient descent. Defaults to TRUE. Setting FALSE will result in slower, but deterministic behavior (if cores=1).
#' @param umap.nn_method The method used for nearest neighbor network construction during UMAP.
#' @param k number of kNN used in creating the k nearest neighbor graph for Louvain clustering. The number of kNN is related to the resolution of the clustering result, bigger number of kNN gives low resolution and vice versa. Default to be 20
#' @param leiden_iter Integer number of iterations used for Leiden clustering. The clustering result with the largest modularity score is used as the final clustering result.  Default to be 1.
#' @param partition_qval Significance threshold used in Louvain community graph partitioning.
#' @param weight A logic argument to determine whether or not we will use
#'   Jaccard coefficient for two nearest neighbors (based on the overlapping of
#'   their kNN) as the weight used for Louvain clustering. Default to be FALSE.
#' @param resolution Resolution parameter passed to Louvain. Can be a list. If
#'   so, this method will evaluate modularity at each resolution and use the
#'   one with the highest value.
#' @param random_seed  the seed used by the random number generator in Leiden.
#' @param cores number of cores computer should use to execute function
#' @param verbose Whether or not verbose output is printed.
#' @param ... Additional arguments passed to UMAP and Louvain analysis.
#'
#' @return A dataframe with genes and the modules to which they are assigned.
#'
#' @export
find_gene_modules <- function(cds,
                          reduction_method = c("UMAP"),
                          max_components = 2,
                          umap.metric = "cosine",
                          umap.min_dist = 0.1,
                          umap.n_neighbors = 15L,
                          umap.fast_sgd = FALSE,
                          umap.nn_method = "annoy",
                          k = 20,
                          leiden_iter = 1,
                          partition_qval = 0.05,
                          weight = FALSE,
                          resolution = NULL,
                          random_seed = 0L,
                          cores=1,
                          verbose = F,
                          ...) {
  method = 'leiden'
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'UMAP', 'PCA' or 'tSNE'")

  reduction_method <- match.arg(reduction_method)

  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(is.character(reduction_method))
  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(is.logical(weight))
  assertthat::assert_that(assertthat::is.count(leiden_iter))
  ## TO DO what is resolution?
  assertthat::assert_that(is.numeric(partition_qval))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before running cluster_cells"))

  preprocess_mat <- cds@preprocess_aux$gene_loadings
  if (is.null(cds@preprocess_aux$beta) == FALSE){
    preprocess_mat <- sweep( preprocess_mat, 2, cds@preprocess_aux$beta[,1], '*')
  }
  preprocess_mat <- preprocess_mat[intersect(rownames(cds), row.names(preprocess_mat)),]

  # uwot::umap uses a random number generator
  if( random_seed != 0L )
    set.seed( random_seed )

  umap_res = uwot::umap(as.matrix(preprocess_mat),
                        n_components = max_components,
                        metric = umap.metric,
                        min_dist = umap.min_dist,
                        n_neighbors = umap.n_neighbors,
                        fast_sgd = umap.fast_sgd,
                        n_threads=cores,
                        verbose=verbose,
                        nn_method= umap.nn_method,
                        ...)

  row.names(umap_res) <- row.names(preprocess_mat)
  colnames(umap_res) <- paste0('dim_', 1:ncol(umap_res))
  reduced_dim_res <- umap_res

  if(verbose)
    message("Running leiden clustering algorithm ...")

  cluster_result <- leiden_clustering(data = reduced_dim_res,
                                    pd = rowData(cds)[
                                      row.names(reduced_dim_res),,drop=FALSE],
                                    k = k,
                                    weight = weight,
                                    num_iter = leiden_iter,
                                    resolution_parameter = resolution,
                                    random_seed = random_seed,
                                    verbose = verbose, ...)

  cluster_graph_res <- compute_partitions(cluster_result$g,
                                          cluster_result$optim_res,
                                          partition_qval, verbose)
  partitions <-
    igraph::components(cluster_graph_res$cluster_g)$membership[
      cluster_result$optim_res$membership]
  names(partitions) <- row.names(reduced_dim_res)
  partitions <- as.factor(partitions)

  gene_module_df <- tibble::tibble(id = row.names(preprocess_mat),
                                   module = factor(
                                     igraph::membership(cluster_result$optim_res)),
                                   supermodule = partitions)
  gene_module_df <- tibble::as_tibble(cbind(gene_module_df, umap_res))

  return(gene_module_df)
}

#' A function to aggregate columns within a matrix.
#' @keywords internal
my.aggregate.Matrix = function (x, groupings = NULL, form = NULL, fun = "sum", ...)
{
  if (!methods::is(x, "Matrix"))
    x <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
  if (fun == "count")
    x <- x != 0
  groupings2 <- data.frame(A=as.factor(groupings))
  if (is.null(form))
    form <- stats::as.formula("~0+.")
  form <- stats::as.formula(form)
  mapping <- Matrix.utils::dMcast(groupings2, form)
  colnames(mapping) <- substring(colnames(mapping), 2)
  result <- Matrix::t(mapping) %*% x
  if (fun == "mean")
    result <- result/as.numeric(table(groupings)[rownames(result)])
  attr(result, "crosswalk") <- grr::extract(groupings, match(rownames(result),
                                                             groupings2$A))
  return(result)
}

#' Creates a matrix with aggregated expression values for arbitrary groups of
#' genes
#'
#' @param cds The cell_data_set on which this function operates
#' @param gene_group_df A dataframe in which the first column contains gene ids
#'   or short gene names and the second contains groups. If NULL, genes are not
#'   grouped.
#' @param cell_group_df A dataframe in which the first column contains cell ids
#'   and the second contains groups. If NULL, cells are not grouped.
#' @param norm_method How to transform gene expression values before
#'   aggregating them. If "log", a pseudocount is added. If "size_only", values
#'   are divided by cell size factors prior to aggregation.
#' @param pseudocount Value to add to expression prior to log transformation
#'   and aggregation.
#' @param scale_agg_values Whether to center and scale aggregated groups of
#'   genes.
#' @param max_agg_value If scale_agg_values is TRUE, the maximum value the
#'   resulting Z scores can take. Higher values are capped at this threshold.
#' @param min_agg_value If scale_agg_values is TRUE, the minimum value the
#'   resulting Z scores can take. Lower values are capped at this threshold.
#' @param exclude.na Logical indicating whether or not to exclude NA values
#'   from the aggregated matrix.
#'
#' @return A matrix of dimension NxM, where N is the number of gene groups and
#'   M is the number of cell groups.
#' @export
aggregate_gene_expression <- function(cds,
                                      gene_group_df = NULL,
                                      cell_group_df = NULL,
                                      norm_method=c("log", "binary",
                                                    "size_only"),
                                      pseudocount=1,
                                      scale_agg_values=TRUE,
                                      max_agg_value=3,
                                      min_agg_value=-3,
                                      exclude.na=TRUE){
  if (is.null(gene_group_df) && is.null(cell_group_df))
    stop("Error: one of either gene_group_df or cell_group_df must not be NULL")
  agg_mat <- normalized_counts(cds, norm_method=norm_method,
                               pseudocount=pseudocount)
  if (is.null(gene_group_df) == FALSE){
    gene_group_df <- as.data.frame(gene_group_df)
    gene_group_df <- gene_group_df[gene_group_df[,1] %in%
                                     fData(cds)$gene_short_name |
                                     gene_group_df[,1] %in%
                                     row.names(fData(cds)),,drop=FALSE]

    # Convert gene short names to rownames if necessary. The more
    # straightforward single call to recode took much longer.
    # Thanks to Christopher Johnstone who posted this on github.
    short_name_mask <- gene_group_df[[1]] %in% fData(cds)$gene_short_name
    if (any(short_name_mask)) {
      geneids <- as.character(gene_group_df[[1]])
      geneids[short_name_mask] <- row.names(fData(cds))[match(
                  geneids[short_name_mask], fData(cds)$gene_short_name)]
      gene_group_df[[1]] <- geneids
    }

    # gene_group_df = gene_group_df[row.names(fData(cds)),]

    # FIXME: this should allow genes to be part of multiple groups. group_by
    # over the second column with a call to colSum should do it.
    agg_mat = as.matrix(my.aggregate.Matrix(agg_mat[gene_group_df[,1],],
                                            as.factor(gene_group_df[,2]),
                                            fun="sum"))
	if (scale_agg_values){
      agg_mat <- t(scale(t(agg_mat)))
      agg_mat[agg_mat < min_agg_value] <- min_agg_value
      agg_mat[agg_mat > max_agg_value] <- max_agg_value
    }
  }

  if (is.null(cell_group_df) == FALSE){

    cell_group_df <- as.data.frame(cell_group_df)
    cell_group_df <- cell_group_df[cell_group_df[,1] %in% row.names(pData(cds)),,
                                  drop=FALSE]
    agg_mat <- agg_mat[,cell_group_df[,1]]
    agg_mat <- my.aggregate.Matrix(Matrix::t(agg_mat),
                                  as.factor(cell_group_df[,2]),
                                  fun="mean")
    agg_mat <- Matrix::t(agg_mat)
  }

  if (exclude.na){
    agg_mat <- agg_mat[rownames(agg_mat) != "NA", colnames(agg_mat) != "NA",drop=FALSE]
  }
  return(agg_mat)
}
