# These functions are due to Maddy Duran.

# Purpose: return the most frequently occurring cell label.
# Parameters:
#   x  list of cell labels.
which_mode <- function(x, top_threshold=0.5, top_over_second_threshold=1.5) {
  # Make a contigency table of cell labels sorted from
  # most to least frequently occurring.
  ta <- sort(table(x), decreasing=TRUE)
  tam <- dplyr::first(ta)
  freq <- tam/length(x)
  top_to_second <- dplyr::nth(ta, 1)/dplyr::nth(ta, 2)
  if (freq > top_threshold)
    mod <- names(ta)[ta == tam]
  else if (top_to_second >= top_over_second_threshold)
    mod <- names(ta)[ta == tam]
  else
    mod <- NA_character_
  return(mod)
}


# Purpose: get cell labels of nearest neighbor cells for discrete values.
# Called by: transfer_cell_labels
# Calls: which_mode
get_nn_cell_label <- function(query_data, query_search, ref_colData_df, cell_label_type, top_threshold=0.5, top_over_second_threshold=1.5) {
  # Loop through the query cells.
  query_nns <- sapply(seq(1, nrow(query_data)), function(i) {
    # Get labels of neighboring cells in reference space.
    ref_neighbors <- query_search[['idx']][i,]
    # Get corresponding cell labels.
    ref_labels <- ref_colData_df[ref_neighbors, cell_label_type]
    # Find modal cell label over a threshold.
    top_label <- which_mode(ref_labels, top_threshold=top_threshold, top_over_second_threshold=top_over_second_threshold)
  })
}


# Purpose: get means from nearest neighbor cells for continuous values.
get_nn_means <- function(query_data, query_search, ref_colData_df, cell_label_type) {
  query_nns <- sapply(seq(1, nrow(query_data)), function(i) {
    # Get labels of neighboring cells in reference space.
    ref_neighbors <- query_search[['idx']][i,]
    # Get corresponding reference cell label.
    ref_labels <- ref_colData_df[ref_neighbors, cell_label_type]
    # Find the modal cell label.
    top_label <- mean(ref_labels, na.rm=TRUE)
  })
}


# Purpose: transfer labels from reference to query.
# Calls: get_nn_cell_label or get_nn_means
# Notes:
#   o  on entry, the cds contains the query data set that was
#      projected into the reference space
#   o  the reference-related reduce_dim_aux (model/nn) information
#      is loaded in place of the query version (is this
#      necessary? Wasn't it loaded in order to project the
#      query cells? Maybe the nn indexes (in the query cds) were
#      recalculated using the query cells?)
#   o  the reference annoy index is searched for (reference cell)
#      nn's to the query data set cells
#   o  test the ref_colData for the required columns
#   o  the cell_label_type value must be in the colnames(colData(cds)).
transfer_cell_labels <- function(cds, reduction_method=c('UMAP', 'PCA'), in_model_dir, ref_colData, cell_label_type, nn_method=c('annoy'), k=10, search_k=100 * k, cores=1, top_threshold=0.5, top_over_second_threshold=1.5) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste0('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be 'UMAP' or 'PCA'")
  assertthat::assert_that(cell_label_type %in% colnames(ref_colData_df),
                          msg=paste0('cell_label_type \'', cell_label_type, '\' is not in the ref_colData'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be 'annoy'")
  nn_method <- match.arg(nn_method)

  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(assertthat::is.count(search_k))

  assertthat::assert_that(is.double(top_threshold),
                          msg=paste0('top_threshold value is not numeric'))
  assertthat::assert_that(is.double(top_over_second_threshold),
                          msg=paste0('top_over_second_threshold value is not numeric'))
  reduction_method <- match.arg(reduction_method)

  ref_colData_df <- ref_colData
  if(!is.data.frame(ref_colData_df)) {
    ref_colData_df <- as.data.frame(ref_colData_df)
  }

  # Are the cell_label_type values discrete?
  label_data_are_discrete <- !is.double(colData(ref_colData_df[[cell_label_type]])[1])

  # Load the reference projection models and nn indexes
  # into the query cds.
  cds <- load_transform_models(cds, in_model_dir)
  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]]),
                          msg=paste0("Reduction Method '", reduction_method, "' is not in the",
                                    "loaded model object."))
  colData(cds)[[cell_label_type]] <- NULL



  # Search the reference reduction_method space for nearest neighbors
  # to the query cells.
  # The cds_reduced_dims contains the query cell coordinates
  # after projection into the reference space.
  # The cds@reduce_dim_aux[[reduction_method]] contains the reduction_method
  # coordinates for the reference data set, which were
  # loaded using load_transform_models() above.
  cds_res <- search_nn_index(cds, reduction_method=reduction_method, nn_method, k=k, search_k=search_k, cores=cores)
 
  # Get the best reference cell label for the query cells.
  cds_reduced_dims <- reducedDims(cds)[[reduction_method]]
  if(label_data_are_discrete) {
    nn_labels <- get_nn_cell_label(cds_reduced_dims, 
                                   cds_res, 
                                   ref_colData_df, 
                                   cell_label_type=cell_label_type,
                                   top_threshold=top_threshold,
                                   top_over_second_threshold=top_over_second_threshold)
  } else {
    nn_labels <- get_nn_means(cds_reduced_dims,
                              cds_res,
                              ref_colData_df,
                              cell_label_type=cell_label_type)
  }
  
  colData(cds)[[cell_label_type]] <- nn_labels

  return(cds)
}


# fill in NAs ----------------------------------------------------------------

edit_cell_label <- function(curr_label, other_labels, top_threshold=0.5, top_over_second_threshold=1.5) {
  # Change NAs to strings.
  curr_label <- replace_na(curr_label, "NA")
  other_labels <- replace_na(other_labels, "NA")
  top_label <- which_mode(other_labels, top_threshold, top_over_second_threshold)
  
  # Switch back if necessary.
  top_label <- gsub("NA", NA_character_, top_label)
  
  return(top_label)
}


edit_query_cell_labels <- function(preproc_res, query_colData, query_nn_index,
                                   cell_label_type,
                                   nn_method='annoy', k=10, search_k=100 * k,
                                   top_threshold=0.5, top_over_second_threshold=1.5, ...) {


  query_search <- uwot:::annoy_search(X=preproc_res, k=k+1, ann=query_nn_index, search_k=search_k, ...)
  query_nns <- sapply(seq(1, nrow(query_search[['idx']])), function(i) {
    # Get neighbors in reference space.
    query_neighbors <- query_search[['idx']][i,]
    # Get corresponding reference cell label.
    query_labels <- query_colData[query_neighbors, cell_label_type]
    curr_label <- query_labels[1]
    other_labels <- query_labels[2:k]
    top_label <- edit_cell_label(curr_label=curr_label,
                                 other_labels=other_labels,
                                 top_threshold=top_threshold,
                                 top_over_second_threshold=top_over_second_threshold)
  })
}


# Purpose: replace missing cell labels.
# Notes:
#   the cell_label_type value must be in the colnames(colData(cds)).
fix_missing_cell_labels <- function(cds, reduction_method=c('UMAP', 'PCA'), out_notna_model_dir=NULL,
                                    cell_label_type=NULL,
                                    nn_method=c('annoy'),  nn_metric=c('euclidean', 'cosine', 'manhattan', 'hamming'),
                                    n_trees=50, k=10, search_k=100 * k,
                                    top_threshold=0.5, top_over_second_threshold=1.5, ...) {
  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be 'UMAP' or 'PCA'")
  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]]),
                          msg=paste0("Reduction method '", reduction_method, "' is not in the cds."))
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(cell_label_type %in% colnames(colData(cds)),
                          msg=paste0('cell_label_type \'', cell_label_type, '\' is not in the cds colData'))

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_method must be 'annoy'")
  nn_method <- match.arg(nn_method)

  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(nn_metric) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "nn_metric must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'")

  nn_metric <- match.arg(nn_metric)

  assertthat::assert_that(assertthat::is.count(n_trees))
  
  na_cds <- cds[, is.na(colData(cds)[[cell_label_type]])]
  
  # Build on rest of cds, where there is a label.
  notna_cds <- cds[, !is.na(colData(cds)[[cell_label_type]])]
  notna_cds <- build_nn_index(notna_cds, 
                              reduction_method=reduction_method,
                              nn_method=nn_method,
                              nn_metric=nn_metric,
                              n_trees=n_trees)
  
  notna_nn_index <- notna_cds@reduce_dim_aux[[reduction_method]][['nn_index']][['annoy_index']]
  notna_colData <- as.data.frame(colData(notna_cds))

  if(!is.null(out_notna_model_dir)) {
    save_transform_models(notna_cds, out_notna_model_dir)
  }
  
  new_cell_labels <- edit_query_cell_labels(preproc_res=reducedDims(na_cds)[[reduction_method]],
                                            query_colData=notna_colData, 
                                            query_nn_index=notna_nn_index, 
                                            cell_label_type=cell_label_type,
                                            nn_method=nn_method,
                                            k=k,
                                            search_k=search_k,
                                            top_threshold=top_threshold,
                                            top_over_second_threshold=top_over_second_threshold,
                                            ...)
  
  return(new_cell_labels)
}

