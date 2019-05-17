#' @title Test genes for differential expression based on the low dimensional embedding and the principal graph
#'
#' @description We are often interested in finding genes that are differentially expressed across a single-cell trajectory.
#' Monocle 3 introduces a new approach for finding such genes that draws on a powerful technique in spatial
#' correlation analysis, the Moran’s I test. Moran’s I is a measure of multi-directional and multi-dimensional
#' spatial autocorrelation. The statistic tells you whether cells at nearby positions on a trajectory will have
#' similar (or dissimilar) expression levels for the gene being tested. Although both Pearson correlation and
#' Moran’s I ranges from -1 to 1, the interpretation of Moran’s I is slightly different: +1 means that nearby
#' cells will have perfectly similar expression; 0 represents no correlation,
#' and -1 means that neighboring cells will be *anti-correlated*.
#'
#' @param cds a cell_data_set object upon which to perform this operation
# #' @param landmark_num Number of landmark cells selected for performing aggregate Moran's I test, default is NULL (no landmark selection and all cells are used)
#' @param k Number of nearest neighbors used for building the kNN graph which is passed to knn2nb function during the Moran's I (Geary's C) test procedure.
#' @param method a character string specifying the method (the default 'Moran_I' or 'Geary_C') for detecting significant genes showing correlated genes along the principal graph embedded in the low dimensional space.
#' @param alternative a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.
#' @param cores the number of cores to be used while testing each gene for differential expression.
#' @param interactive Whether or not to allow the user to choose a point or region in the scene, then to only identify genes spatially correlated for those selected cells.
#' @param verbose Whether to show spatial test (Moran's I or Gearys' C test) errors and warnings. Only valid for cores = 1.
#' @return a data frame containing the p values and q-values from the Moran's I test on the parallel arrays of models.
#' @seealso \code{\link[spdep]{moran.test}} \code{\link[spdep]{geary.test}}
#' @export
graph_test <- function(cds,
                                 neighbor_graph = c("principal_graph", "knn"),
                                 reduction_method = "UMAP",
                                 k = 25,
                                 method = c('Moran_I'),
                                 alternative = 'greater',
                                 expression_family="quasipoisson",
                                 cores=1,
                                 interactive = FALSE,
                                 verbose=FALSE) {
  lw <- calculateLW(cds,
                    k = k,
                    interactive = interactive,
                    verbose = verbose,
                    neighbor_graph = neighbor_graph,
                    reduction_method = reduction_method)

  if(verbose) {
    message("Performing Moran's I test: ...")
  }
  exprs_mat <- counts(cds)[, attr(lw, "region.id"), drop=FALSE]
  sz <- size_factors(cds)[attr(lw, "region.id")]

  wc <- spdep::spweights.constants(lw, zero.policy = TRUE, adjust.n = TRUE)
  test_res <- pbmcapply::pbmclapply(row.names(exprs_mat), FUN = function(x, sz, alternative, method, expression_family) {
    exprs_val <- exprs_mat[x, ]

    if (expression_family %in% c("uninormal", "binomialff")){
      exprs_val <- exprs_val
    }else{
        exprs_val <- log10(exprs_val / sz + 0.1)
    }

    test_res <- tryCatch({
      if(method == "Moran_I") {
        mt <- my.moran.test(exprs_val, lw, wc, alternative = alternative)
        data.frame(status = 'OK', p_value = mt$p.value, morans_test_statistic = mt$statistic, morans_I = mt$estimate[["Moran I statistic"]])
      } else if(method == 'Geary_C') {
        gt <- my.geary.test(exprs_val, lw, wc, alternative = alternative)
        data.frame(status = 'OK', p_value = gt$p.value, geary_test_statistic = gt$statistic, geary_C = gt$estimate[["Geary C statistic"]])
      }
    },
    error = function(e) {
      data.frame(status = 'FAIL', p_value = NA, morans_test_statistic = NA, morans_I = NA)
    })
  }, sz = sz, alternative = alternative, method = method, expression_family = expression_family, mc.cores=cores, ignore.interactive = TRUE)

  if(verbose) {
    message("returning results: ...")
  }

  test_res <- do.call(rbind.data.frame, test_res)
  row.names(test_res) <- row.names(cds)
  test_res <- merge(test_res, rowData(cds), by="row.names")
  row.names(test_res) <- test_res[, 1] #remove the first column and set the row names to the first column
  test_res[, 1] <- NULL
  test_res$q_value <- 1
  test_res$q_value[which(test_res$status == 'OK')] <- stats::p.adjust(subset(test_res, status == 'OK')[, 'p_value'], method="BH")
  test_res$status = as.character(test_res$status)
  test_res[row.names(cds), ] # make sure gene name ordering in the DEG test result is the same as the CDS
}

my.moran.test <- function (x, listw, wc, alternative = "greater", randomisation = TRUE) {
  zero.policy = TRUE
  adjust.n = TRUE
  na.action = stats::na.fail
  drop.EI2 = FALSE
  xname <- deparse(substitute(x))
  wname <- deparse(substitute(listw))
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")

  S02 <- wc$S0 * wc$S0
  res <- spdep::moran(x, listw, wc$n, wc$S0, zero.policy = zero.policy,
                      NAOK = NAOK)
  I <- res$I
  K <- res$K

  EI <- (-1)/wc$n1
  if (randomisation) {
    VI <- wc$n * (wc$S1 * (wc$nn - 3 * wc$n + 3) - wc$n *
                    wc$S2 + 3 * S02)
    tmp <- K * (wc$S1 * (wc$nn - wc$n) - 2 * wc$n * wc$S2 +
                  6 * S02)
    if (tmp > VI)
      warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
    VI <- (VI - tmp)/(wc$n1 * wc$n2 * wc$n3 * S02)
    if (!drop.EI2)
      VI <- (VI - EI^2)
    if (VI < 0)
      warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  }
  else {
    VI <- (wc$nn * wc$S1 - wc$n * wc$S2 + 3 * S02)/(S02 *
                                                      (wc$nn - 1))
    if (!drop.EI2)
      VI <- (VI - EI^2)
    if (VI < 0)
      warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  }
  ZI <- (I - EI)/sqrt(VI)
  statistic <- ZI
  names(statistic) <- "Moran I statistic standard deviate"
  if (alternative == "two.sided")
    PrI <- 2 * stats::pnorm(abs(ZI), lower.tail = FALSE)
  else if (alternative == "greater")
    PrI <- stats::pnorm(ZI, lower.tail = FALSE)
  else PrI <- stats::pnorm(ZI)
  if (!is.finite(PrI) || PrI < 0 || PrI > 1)
    warning("Out-of-range p-value: reconsider test arguments")
  vec <- c(I, EI, VI)
  names(vec) <- c("Moran I statistic", "Expectation", "Variance")
  method <- paste("Moran I test under", ifelse(randomisation,
                                               "randomisation", "normality"))

  res <- list(statistic = statistic, p.value = PrI, estimate = vec)
  if (!is.null(na.act))
    attr(res, "na.action") <- na.act
  class(res) <- "htest"
  res
}

my.geary.test <- function (x, listw, wc, randomisation = TRUE, alternative = "greater")
{
  zero.policy = TRUE
  adjust.n = TRUE
  spChk = NULL
  alternative <- match.arg(alternative, c("less", "greater",
                                          "two.sided"))
  if (!inherits(listw, "listw"))
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (!is.numeric(x))
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  if (any(is.na(x)))
    stop("NA in X")
  n <- length(listw$neighbours)
  if (n != length(x))
    stop("objects of different length")
  if (is.null(spChk))
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw))
    stop("Check of data and weights ID integrity failed")
  S02 <- wc$S0 * wc$S0
  res <- spdep::geary(x, listw, wc$n, wc$n1, wc$S0, zero.policy)
  C <- res$C
  if (is.na(C))
    stop("NAs generated in geary - check zero.policy")
  K <- res$K
  EC <- 1
  if (randomisation) {
    VC <- (wc$n1 * wc$S1 * (wc$nn - 3 * n + 3 - K * wc$n1))
    VC <- VC - ((1/4) * (wc$n1 * wc$S2 * (wc$nn + 3 * n -
                                            6 - K * (wc$nn - n + 2))))
    VC <- VC + (S02 * (wc$nn - 3 - K * (wc$n1^2)))
    VC <- VC/(n * wc$n2 * wc$n3 * S02)
  }
  else {
    VC <- ((2 * wc$S1 + wc$S2) * wc$n1 - 4 * S02)/(2 * (n +
                                                          1) * S02)
  }
  ZC <- (EC - C)/sqrt(VC)
  statistic <- ZC
  names(statistic) <- "Geary C statistic standard deviate"
  PrC <- NA
  if (is.finite(ZC)) {
    if (alternative == "two.sided")
      PrC <- 2 * stats::pnorm(abs(ZC), lower.tail = FALSE)
    else if (alternative == "greater")
      PrC <- stats::pnorm(ZC, lower.tail = FALSE)
    else PrC <- stats::pnorm(ZC)
    if (!is.finite(PrC) || PrC < 0 || PrC > 1)
      warning("Out-of-range p-value: reconsider test arguments")
  }
  vec <- c(C, EC, VC)
  names(vec) <- c("Geary C statistic", "Expectation", "Variance")
  method <- paste("Geary C test under", ifelse(randomisation,
                                               "randomisation", "normality"))
  data.name <- paste(deparse(substitute(x)), "\nweights:",
                     deparse(substitute(listw)), "\n")
  res <- list(statistic = statistic, p.value = PrC, estimate = vec,
              alternative = ifelse(alternative == "two.sided", alternative,
                                   paste("Expectation", alternative, "than statistic")),
              method = method, data.name = data.name)
  class(res) <- "htest"
  res
}

#' Function to calculate the neighbors list with spatial weights for the chosen
#' coding scheme from a cell dataset object
#'
#' @description This function first retrieves the association from each cell to
#' any principal points, then builds a kNN graph for all cells
#' and removes edges that connected between groups that disconnected in the
#' corresponding principal graph and finally uses this kNN graph to calculate a
#' global Moran's I and get the p-value
#' @param cds The cell_data_set object where the neighbors list is calculated from
#' @param  k The maximum number of nearest neighbors to compute
#' @param verbose A logic flag that determines whether or not to print
#' execution details
#' @param interactive Whether or not to allow the user to choose a point or
#' region in the scene, then to only identify genes spatially correlated for
#' those selected cells.
#' @keywords internal
#'
calculateLW <- function(cds,
                        k,
                        neighbor_graph,
                        reduction_method,
                        interactive = FALSE,
                        verbose = FALSE
                        ) {
  if(verbose) {
    message("retrieve the matrices for Moran's I test...")
  }
  knn_res <- NULL
  principal_g <- NULL

  cell_coords <- reducedDims(cds)$UMAP
  if (neighbor_graph == "knn") {
    knn_res <- RANN::nn2(cell_coords, cell_coords, min(k + 1, nrow(cell_coords)), searchtype = "standard")[[1]]
  } else if(neighbor_graph == "principal_graph") {
    pr_graph_node_coords <- cds@principal_graph_aux[[reduction_method]]$dp_mst
    principal_g <-  igraph::get.adjacency(cds@principal_graph[[reduction_method]])[colnames(pr_graph_node_coords), colnames(pr_graph_node_coords)]
  }

  exprs_mat <- counts(cds)
  if(neighbor_graph == "knn") {
    if(is.null(knn_res)) {
      cell_coords <- reducedDims(cds)$UMAP
      knn_res <- RANN::nn2(cell_coords, cell_coords, min(k + 1, nrow(cell_coords)), searchtype = "standard")[[1]]
    }
    links <- jaccard_coeff(knn_res[, -1], F)
    links <- links[links[, 1] > 0, ]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
    knn_res_graph <- igraph::graph.data.frame(relations, directed = T)

    if(interactive) {
      if('Marked' %in% names(colData(cds))) {
        points_selected <- which(colData(cds)$Marked)
      } else {
        cat("Left click or drag multiple points to select a group of cells\n")
        cat("Press <Esc> in the rgl screen to exit \n")
        points_selected <- sort(select_cells(cds))
      }
      knn_list <- lapply(points_selected, function(x) intersect(knn_res[x, -1], points_selected))
      knn_res_graph <- knn_res_graph[points_selected, points_selected]
      region_id_names <- colnames(cds)[points_selected]

      id_map <- 1:length(points_selected)
      names(id_map) <- points_selected
    } else {
      knn_list <- lapply(1:nrow(knn_res), function(x) knn_res[x, -1])
      region_id_names <- colnames(cds)

      id_map <- 1:ncol(cds)
      names(id_map) <- id_map

      points_selected <- 1:nrow(knn_res)
    }

    knn_list <- lapply(points_selected, function(x) id_map[as.character(knn_res[x, -1])])
  }
  else if (neighbor_graph == "principal_graph") {
    cell2pp_map <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex # mapping from each cell to the principal points
    if(is.null(cell2pp_map)) {
      stop("Error: projection matrix for each cell to principal points doesn't exist, you may need to rerun learnGraph")
    }

    # This cds object might be a subset of the one on which ordering was performed,
    # so we may need to subset the nearest vertex and low-dim coordinate matrices:
    cell2pp_map <-  cell2pp_map[row.names(cell2pp_map) %in% row.names(colData(cds)),, drop=FALSE]
    cell2pp_map <- cell2pp_map[colnames(cds), ]

    if(verbose) {
      message("Identify connecting principal point pairs ...")
    }
    # an alternative approach to make the kNN graph based on the principal graph
    knn_res <- RANN::nn2(cell_coords, cell_coords, min(k + 1, nrow(cell_coords)), searchtype = "standard")[[1]]
    # kNN_res_pp_map <- matrix(cell2pp_map[knn_res], ncol = k + 1, byrow = F) # convert the matrix of knn graph from the cell IDs into a matrix of principal points IDs

    principal_g_tmp <- principal_g # kNN can be built within group of cells corresponding to each principal points
    diag(principal_g_tmp) <- 1 # so set diagnol as 1
    cell_membership <- as.factor(cell2pp_map)
    uniq_member <- sort(unique(cell_membership))

    membership_matrix <- Matrix::sparse.model.matrix( ~ cell_membership + 0)
    colnames(membership_matrix) <- levels(uniq_member)
    # sparse matrix multiplication for calculating the feasible space
    feasible_space <- membership_matrix %*% Matrix::tcrossprod(principal_g_tmp[as.numeric(levels(uniq_member)), as.numeric(levels(uniq_member))], membership_matrix)

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

    #pb_feasible_knn <- utils::txtProgressBar(max = num_blocks, file = "", style = 3, min = 0)
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
        tmp <- rBind(tmp, cur_tmp)
      }
      #utils::setTxtProgressBar(pb = pb_feasible_knn, value = pb_feasible_knn$getVal() + 1)
    }

    #close(pb_feasible_knn)
    if(verbose) {
      message('Calculating valid kNN graph, done ...')
    }

    if(interactive) {
      if('Marked' %in% names(colData(cds))) {
        points_selected <- which(colData(cds)$Marked)
      } else {
        cat("Left click or drog multiple points to select a group of cells\n")
        cat("Press <Esc> in the rgl screen to exit \n")
        points_selected <- sort(select_cells(cds, return_index = T))
      }
      tmp <- tmp[points_selected, points_selected]
      region_id_names <- colnames(cds)[points_selected]

      id_map <- 1:length(points_selected)
      names(id_map) <- points_selected
    } else {
      region_id_names <- colnames(cds)

      id_map <- 1:ncol(cds)
      names(id_map) <- id_map
    }

    knn_list <- slam::rowapply_simple_triplet_matrix(slam::as.simple_triplet_matrix(tmp), function(x) {
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
  lw <- spdep::nb2listw(knn_list, zero.policy = TRUE)
  lw
}

