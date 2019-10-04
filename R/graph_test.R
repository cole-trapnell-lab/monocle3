#' Test genes for differential expression based on the low dimensional
#' embedding and the principal graph
#'
#' @description We are often interested in finding genes that are
#' differentially expressed across a single-cell trajectory. Monocle3
#' introduces a new approach for finding such genes that draws on a powerful
#' technique in spatial correlation analysis, the Moran’s I test. Moran’s I is
#' a measure of multi-directional and multi-dimensional spatial
#' autocorrelation. The statistic tells you whether cells at nearby positions
#' on a trajectory will have similar (or dissimilar) expression levels for the
#' gene being tested. Although both Pearson correlation and Moran’s I ranges
#' from -1 to 1, the interpretation of Moran’s I is slightly different: +1
#' means that nearby cells will have perfectly similar expression; 0 represents
#' no correlation, and -1 means that neighboring cells will be
#' *anti-correlated*.
#'
#' @param cds a cell_data_set object upon which to perform this operation
#' @param neighbor_graph String indicating what neighbor graph to use.
#'   "principal_graph" and "knn" are supported. Default is "knn", but
#'   "principal_graph" is recommended for trajectory analysis.
#' @param reduction_method character, the method used to reduce dimension.
#'   Currently only supported for "UMAP".
#' @param k Number of nearest neighbors used for building the kNN graph which
#'   is passed to knn2nb function during the Moran's I (Geary's C) test
#'   procedure.
#' @param method a character string specifying the method (currently only
#'   'Moran_I' is supported) for detecting significant genes showing
#'   correlation along the principal graph embedded in the low dimensional
#'   space.
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of greater (default), less or two.sided.
#' @param expression_family a character string specifying the expression family
#'   function used for the test.
#' @param cores the number of cores to be used while testing each gene for
#'   differential expression.
#' @param verbose Whether to show spatial test (Moran's I) errors and warnings.
#'   Only valid for cores = 1.
#' @return a data frame containing the p values and q-values from the Moran's I
#'   test on the parallel arrays of models.
#' @seealso \code{\link[sdmoran]{moran.test}} \code{\link[sdmoran]{geary.test}}
#' @export
graph_test <- function(cds,
                       neighbor_graph = c("knn", "principal_graph"),
                       reduction_method = "UMAP",
                       k = 25,
                       method = c('Moran_I'),
                       alternative = 'greater',
                       expression_family="quasipoisson",
                       cores=1,
                       verbose=FALSE) {
  neighbor_graph <- match.arg(neighbor_graph)
  lw <- calculateLW(cds,
                    k = k,
                    verbose = verbose,
                    neighbor_graph = neighbor_graph,
                    reduction_method = reduction_method)

  if(verbose) {
    message("Performing Moran's I test: ...")
  }
  exprs_mat <- SingleCellExperiment::counts(cds)[, attr(lw, "region.id"), drop=FALSE]
  sz <- size_factors(cds)[attr(lw, "region.id")]

  wc <- sdmoran::spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)
  test_res <- pbmcapply::pbmclapply(row.names(exprs_mat),
                                    FUN = function(x, sz, alternative,
                                                   method, expression_family) {
    exprs_val <- exprs_mat[x, ]

    if (expression_family %in% c("uninormal", "binomialff")){
      exprs_val <- exprs_val
    }else{
      exprs_val <- log10(exprs_val / sz + 0.1)
    }

    test_res <- tryCatch({
      if(method == "Moran_I") {
        mt <- suppressWarnings(sdmoran::moran.test.alt(exprs_val, lw, wc, alternative = alternative))
        data.frame(status = 'OK', p_value = mt$p.value,
                   morans_test_statistic = mt$statistic,
                   morans_I = mt$estimate[["Moran I statistic"]])
      } else if(method == 'Geary_C') {
        gt <- suppressWarnings(sdmoran::geary.test.alt(exprs_val, lw, wc, alternative = alternative))
        data.frame(status = 'OK', p_value = gt$p.value,
                   geary_test_statistic = gt$statistic,
                   geary_C = gt$estimate[["Geary C statistic"]])
      }
    },
    error = function(e) {
      data.frame(status = 'FAIL', p_value = NA, morans_test_statistic = NA,
                 morans_I = NA)
    })
  }, sz = sz, alternative = alternative, method = method,
  expression_family = expression_family, mc.cores=cores,
  ignore.interactive = TRUE)

  if(verbose) {
    message("returning results: ...")
  }

  test_res <- do.call(rbind.data.frame, test_res)
  row.names(test_res) <- row.names(cds)
  test_res <- merge(test_res, rowData(cds), by="row.names")
  #remove the first column and set the row names to the first column
  row.names(test_res) <- test_res[, 1]
  test_res[, 1] <- NULL
  test_res$q_value <- 1
  test_res$q_value[which(test_res$status == 'OK')] <-
    stats::p.adjust(subset(test_res, status == 'OK')[, 'p_value'], method="BH")
  test_res$status = as.character(test_res$status)
  # make sure gene name ordering in the DEG test result is the same as the CDS
  test_res[row.names(cds), ]
}


#' Function to calculate the neighbors list with spatial weights for the chosen
#' coding scheme from a cell dataset object
#'
#' @description This function first retrieves the association from each cell to
#' any principal points, then builds a kNN graph for all cells
#' and removes edges that connected between groups that disconnected in the
#' corresponding principal graph and finally uses this kNN graph to calculate a
#' global Moran's I and get the p-value
#' @param cds The cell_data_set object where the neighbors list is calculated
#'   from
#' @param  k The maximum number of nearest neighbors to compute
#' @param verbose A logic flag that determines whether or not to print
#' execution details
#' @keywords internal
#'
calculateLW <- function(cds,
                        k,
                        neighbor_graph,
                        reduction_method,
                        verbose = FALSE
) {
  if(verbose) {
    message("retrieve the matrices for Moran's I test...")
  }
  knn_res <- NULL
  principal_g <- NULL

  cell_coords <- reducedDims(cds)[[reduction_method]]
  if (neighbor_graph == "knn") {
    knn_res <- RANN::nn2(cell_coords, cell_coords,
                         min(k + 1, nrow(cell_coords)),
                         searchtype = "standard")[[1]]
  } else if(neighbor_graph == "principal_graph") {
    pr_graph_node_coords <- cds@principal_graph_aux[[reduction_method]]$dp_mst
    principal_g <-
      igraph::get.adjacency(
        cds@principal_graph[[reduction_method]])[colnames(pr_graph_node_coords),
                                                 colnames(pr_graph_node_coords)]
  }

  exprs_mat <- exprs(cds)
  if(neighbor_graph == "knn") {
    if(is.null(knn_res)) {
      knn_res <- RANN::nn2(cell_coords, cell_coords,
                           min(k + 1, nrow(cell_coords)),
                           searchtype = "standard")[[1]]
    }
    links <- jaccard_coeff(knn_res[, -1], F)
    links <- links[links[, 1] > 0, ]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
    knn_res_graph <- igraph::graph.data.frame(relations, directed = T)

      knn_list <- lapply(1:nrow(knn_res), function(x) knn_res[x, -1])
      region_id_names <- colnames(cds)

      id_map <- 1:ncol(cds)
      names(id_map) <- id_map

      points_selected <- 1:nrow(knn_res)

    knn_list <- lapply(points_selected,
                       function(x) id_map[as.character(knn_res[x, -1])])
  }
  else if (neighbor_graph == "principal_graph") {
    # mapping from each cell to the principal points
    cell2pp_map <-
      cds@principal_graph_aux[[
        reduction_method]]$pr_graph_cell_proj_closest_vertex
    if(is.null(cell2pp_map)) {
      stop(paste("Error: projection matrix for each cell to principal",
                 "points doesn't exist, you may need to rerun learn_graph"))
    }

    # This cds object might be a subset of the one on which ordering was
    # performed, so we may need to subset the nearest vertex and low-dim
    # coordinate matrices:
    cell2pp_map <-  cell2pp_map[row.names(cell2pp_map) %in%
                                  row.names(colData(cds)),, drop=FALSE]
    cell2pp_map <- cell2pp_map[colnames(cds), ]

    if(verbose) {
      message("Identify connecting principal point pairs ...")
    }
    # an alternative approach to make the kNN graph based on the principal
    # graph
    knn_res <- RANN::nn2(cell_coords, cell_coords,
                         min(k + 1, nrow(cell_coords)),
                         searchtype = "standard")[[1]]
    # convert the matrix of knn graph from the cell IDs into a matrix of
    # principal points IDs
    # kNN_res_pp_map <- matrix(cell2pp_map[knn_res], ncol = k + 1, byrow = F)

    # kNN can be built within group of cells corresponding to each principal
    # points
    principal_g_tmp <- principal_g
    diag(principal_g_tmp) <- 1 # so set diagnol as 1
    cell_membership <- as.factor(cell2pp_map)
    uniq_member <- sort(unique(cell_membership))

    membership_matrix <- Matrix::sparse.model.matrix( ~ cell_membership + 0)
    colnames(membership_matrix) <- levels(uniq_member)
    # sparse matrix multiplication for calculating the feasible space
    feasible_space <- membership_matrix %*%
      Matrix::tcrossprod(principal_g_tmp[as.numeric(levels(uniq_member)),
                                         as.numeric(levels(uniq_member))],
                         membership_matrix)

    links <- jaccard_coeff(knn_res[, -1], F)
    links <- links[links[, 1] > 0, ]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
    knn_res_graph <- igraph::graph.data.frame(relations, directed = T)

    # remove edges across cells belong to two disconnected principal points
    tmp_a <- igraph::get.adjacency(knn_res_graph)
    block_size <- 10000
    num_blocks = ceiling(nrow(tmp_a) / block_size)
    if(verbose) {
      message('start calculating valid kNN graph ...')
    }

    tmp <- NULL

    for (j in 1:num_blocks){
      if (j < num_blocks){
        block_a <- tmp_a[((((j-1) * block_size)+1):(j*block_size)), ]
        block_b <- feasible_space[((((j-1) * block_size)+1):(j*block_size)), ]
      }else{
        block_a <- tmp_a[((((j-1) * block_size)+1):(nrow(tmp_a))), ]
        block_b <- feasible_space[((((j-1) * block_size)+1):(nrow(tmp_a))), ]
      }

      cur_tmp <- block_a * block_b

      if(is.null(tmp)) {
        tmp <- cur_tmp
      } else {
        tmp <- Matrix::rBind(tmp, cur_tmp)
      }
    }

    #close(pb_feasible_knn)
    if(verbose) {
      message('Calculating valid kNN graph, done ...')
    }

      region_id_names <- colnames(cds)

      id_map <- 1:ncol(cds)
      names(id_map) <- id_map

    knn_list <-
      slam::rowapply_simple_triplet_matrix(slam::as.simple_triplet_matrix(tmp),
                                           function(x) {
                                             res <- which(as.numeric(x) > 0)
                                             if(length(res) == 0)
                                               res <- 0L
                                             res
                                           })
  } else {
    stop("Error: unrecognized neighbor_graph option")
  }
  # create the lw list for moran.test
  names(knn_list) <- id_map[names(knn_list)]
  class(knn_list) <- "nb"
  attr(knn_list, "region.id") <- region_id_names
  attr(knn_list, "call") <- match.call()
  # attr(knn_list, "type") <- "queen"
  lw <- sdmoran::nb2listw(knn_list, zero.policy = TRUE)
  lw
}

