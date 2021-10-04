# These functions are due to Maddy Duran.

# Purpose: return the most frequently occurring cell label.
# Parameters:
#   x  list of cell labels.
# which_mode <- function(x, top_threshold=0.5, top_two_ratio_threshold=1.5) {
#   # Make a contigency table of cell labels sorted from
#   # most to least frequently occurring.
#   ta <- sort(table(x), decreasing=TRUE)
#   tam <- dplyr::first(ta)
#   freq <- tam / length(x)
#   top_to_second <- dplyr::nth(ta, 1) / dplyr::nth(ta, 2)
#   if(freq > top_threshold)
#     mod <- names(ta)[ta == tam]
#   else
#   if(top_to_second >= top_two_ratio_threshold)
#     mod <- names(ta)[ta == tam]
#   else
#     mod <- NA_character_
#   return(mod)
# }


# Purpose: return the most frequently occurring cell label.
# Parameters:
#   x  list of cell labels.
which_mode <- function(x, top_threshold=0.5, top_two_ratio_threshold=1.5) {
  # Make a contigency table of cell labels sorted from
  # most to least frequently occurring.
  ta <- sort(table(x), decreasing=TRUE)
  tam <- ta[1]
  freq <- tam / length(x)
  top_to_second <- ta[1] / ta[2]
  if(freq > top_threshold)
    mod <- names(tam)[1]
  else
  if(top_to_second >= top_two_ratio_threshold)
    mod <- names(tam)[1]
  else
    mod <- NA_character_
  return(mod)
}


# Purpose: get cell labels of nearest neighbor cells for discrete values.
# Called by: transfer_cell_labels
# Calls: which_mode
get_nn_cell_label <- function(query_data, query_search, ref_colData, transfer_cell_label, top_threshold=0.5, top_two_ratio_threshold=1.5) {
  # Loop through the query cells.
  query_nns <- sapply(seq(1, nrow(query_data)), function(i) {
    # Get labels of neighboring cells in reference space.
    ref_neighbors <- query_search[['nn.idx']][i,]
    # Get corresponding cell labels.
    ref_labels <- ref_colData[ref_neighbors, transfer_cell_label]
    # Find modal cell label over a threshold.
    top_label <- which_mode(x=ref_labels, top_threshold=top_threshold, top_two_ratio_threshold=top_two_ratio_threshold)
  })
}


# Purpose: get means from nearest neighbor cells for continuous values.
get_nn_means <- function(query_data, query_search, ref_colData, transfer_cell_label) {
  query_nns <- sapply(seq(1, nrow(query_data)), function(i) {
    # Get labels of neighboring cells in reference space.
    ref_neighbors <- query_search[['nn.idx']][i,]
    # Get corresponding reference cell label.
    ref_labels <- ref_colData[ref_neighbors, transfer_cell_label]
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
#   o  the transfer_cell_label value must be in the colnames(colData(cds)).
#   o  keep k=10 or so for confidence in the label transfer
transfer_cell_labels <- function(cds,
                                 reduction_method=c('UMAP', 'PCA', 'LSI'),
                                 in_model_dir,
                                 ref_colData,
                                 transfer_cell_label,
                                 k=10,
                                 nn_control,
                                 top_threshold=0.5,
                                 top_two_ratio_threshold=1.5,
                                 verbose=FALSE) {

  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste0('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be 'UMAP', 'PCA', or 'LSI'")

  assertthat::assert_that(assertthat::is.count(k))

  assertthat::assert_that(is.double(top_threshold),
                          msg=paste0('top_threshold value is not numeric'))
  assertthat::assert_that(is.double(top_two_ratio_threshold),
                          msg=paste0('top_two_ratio_threshold value is not numeric'))
  reduction_method <- match.arg(reduction_method)

  if(!is.data.frame(ref_colData)) {
    ref_colData <- as.data.frame(ref_colData)
  }

  assertthat::assert_that(transfer_cell_label %in% colnames(ref_colData),
                          msg=paste0('transfer_cell_label \'', transfer_cell_label, '\' is not in the ref_colData'))

  nn_control <- set_nn_control(nn_control=nn_control, k=k, method_default='annoy', verbose=verbose)

  # Are the transfer_cell_label values discrete?
  label_data_are_discrete <- !is.double(ref_colData[[transfer_cell_label]][1])

  # Load the reference projection models and nn indexes
  # into the query cds.
  cds <- load_transform_models(cds=cds, directory_path=in_model_dir)
  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]]),
                          msg=paste0("Reduction Method '", reduction_method, "' is not in the",
                                    "loaded model object."))
  colData(cds)[[transfer_cell_label]] <- NULL

  # Search the reference reduction_method space for nearest neighbors
  # to the query cells.
  # The cds_reduced_dims contains the query cell coordinates
  # after projection into the reference space.
  # The cds@reduce_dim_aux[[reduction_method]] contains the reduction_method
  # coordinates for the reference data set, which were
  # loaded using load_transform_models() above.
  cds_reduced_dims <- reducedDims(cds)[[reduction_method]]
  cds_nn_index <- get_cds_nn_index(cds=cds, reduction_method=reduction_method, nn_control=nn_control, verbose=verbose)
  cds_res <- search_nn_index(query_matrix=cds_reduced_dims, nn_index=cds_nn_index,
                             k=k, nn_control=nn_control, verbose=verbose)
 
  # Get the best reference cell label for the query cells.
  if(label_data_are_discrete) {
    nn_labels <- get_nn_cell_label(query_data=cds_reduced_dims, 
                                   query_search=cds_res, 
                                   ref_colData=ref_colData, 
                                   transfer_cell_label=transfer_cell_label,
                                   top_threshold=top_threshold,
                                   top_two_ratio_threshold=top_two_ratio_threshold)
  } else {
    nn_labels <- get_nn_means(query_data=cds_reduced_dims,
                              query_search=cds_res,
                              ref_colData=ref_colData,
                              transfer_cell_label=transfer_cell_label)
  }
  
  colData(cds)[[transfer_cell_label]] <- nn_labels

  return(cds)
}


# fill in NAs ----------------------------------------------------------------

edit_cell_label <- function(curr_label, other_labels, top_threshold=0.5, top_two_ratio_threshold=1.5) {
  # Change NAs to strings.
  curr_label <- tidyr::replace_na(curr_label, "NA")
  other_labels <- tidyr::replace_na(other_labels, "NA")
  top_label <- which_mode(x=other_labels, top_threshold=top_threshold, top_two_ratio_threshold=top_two_ratio_threshold)
  
  # Switch back if necessary.
  top_label <- gsub("NA", NA_character_, top_label)
  
  return(top_label)
}


edit_query_cell_labels <- function(preproc_res,
                                   query_colData,
                                   query_nn_index,
                                   transfer_cell_label,
                                   k=10,
                                   nn_control=nn_control,
                                   top_threshold=0.5,
                                   top_two_ratio_threshold=1.5,
                                   verbose=FALSE) {

  query_search <- search_nn_index(query_matrix=preproc_res,
                                  nn_index=query_nn_index,
                                  k=k+1,
                                  nn_control=nn_control,
                                  verbose=verbose)

  query_nns <- sapply(seq(1, nrow(query_search[['nn.idx']])), function(i) {
    # Get neighbors in reference space.
    query_neighbors <- query_search[['nn.idx']][i,]
    # Get corresponding reference cell label.
    query_labels <- query_colData[query_neighbors, transfer_cell_label]
    curr_label <- query_labels[1]
    other_labels <- query_labels[2:k]
    top_label <- edit_cell_label(curr_label=curr_label,
                                 other_labels=other_labels,
                                 top_threshold=top_threshold,
                                 top_two_ratio_threshold=top_two_ratio_threshold)
  })
}


# Purpose: replace missing cell labels.
# Notes:
#   the transfer_cell_label value must be in the colnames(colData(cds)).
#' export
fix_missing_cell_labels <- function(cds,
                                    reduction_method=c('UMAP', 'PCA', 'LSI'),
                                    out_notna_model_dir=NULL,
                                    transfer_cell_label=NULL,
                                    k=10,
                                    nn_control=nn_control,
                                    top_threshold=0.5,
                                    top_two_ratio_threshold=1.5,
                                    verbose=FALSE) {

  assertthat::assert_that(class(cds) == 'cell_data_set',
                          msg=paste('cds parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be 'UMAP', 'PCA', or 'LSI'")
  assertthat::assert_that(!is.null(cds@reduce_dim_aux[[reduction_method]]),
                          msg=paste0("Reduction method '", reduction_method, "' is not in the cds."))
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(transfer_cell_label %in% colnames(colData(cds)),
                          msg=paste0('transfer_cell_label \'', transfer_cell_label, '\' is not in the cds colData'))

  nn_control <- set_nn_control(nn_control=nn_control, k=k, method_default='annoy', verbose=verbose)

  # Partition cds.
  notna_cds <- cds[, !is.na(colData(cds)[[transfer_cell_label]])]
  
  # Build index on not NA cds, where there is a label.
  notna_nn_index <- make_nn_index(subject_matrix=reducedDims(notna_cds)[[reduction_method]],
                                  nn_control=nn_control,
                                  verbose=verbose)
  notna_colData <- as.data.frame(colData(notna_cds))

  if(!is.null(out_notna_model_dir)) {
    save_transform_models(cds=notna_cds, directory_path=out_notna_model_dir)
  }
  
  na_cds <- cds[, is.na(colData(cds)[[transfer_cell_label]])]
  new_cell_labels <- edit_query_cell_labels(preproc_res=reducedDims(na_cds)[[reduction_method]],
                                            query_colData=notna_colData, 
                                            query_nn_index=notna_nn_index, 
                                            transfer_cell_label=transfer_cell_label,
                                            k=k,
                                            nn_control=nn_control,
                                            top_threshold=top_threshold,
                                            top_two_ratio_threshold=top_two_ratio_threshold,
                                            verbose=verbose)
  
  return(new_cell_labels)
}

