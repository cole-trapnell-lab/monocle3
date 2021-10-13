# These functions are due to Maddy Duran.

# Purpose: return the most frequently occurring cell label.
# Parameters:
#   x  list of cell labels.
which_mode <- function(x, top_frac_threshold=0.5, top_next_ratio_threshold=1.5) {
  # Make a contigency table of cell labels sorted from
  # most to least frequently occurring.
  ta <- sort(table(x), decreasing=TRUE)
  tam <- ta[1]
  freq <- tam / length(x)
  top_to_second <- ta[1] / ta[2]
  if(freq > top_frac_threshold)
    mod <- names(tam)[1]
  else
  if(top_to_second >= top_next_ratio_threshold)
    mod <- names(tam)[1]
  else
    mod <- NA_character_
  return(mod)
}


# Purpose: get cell labels of nearest neighbor cells for discrete values.
# Called by: transfer_cell_labels
# Calls: which_mode
get_nn_cell_label <- function(query_data, query_search, ref_col_data, ref_column_name, top_frac_threshold=0.5, top_next_ratio_threshold=1.5) {
  # Loop through the query cells.
  query_nns <- sapply(seq(1, nrow(query_data)), function(i) {
    # Get labels of neighboring cells in reference space.
    ref_neighbors <- query_search[['nn.idx']][i,]
    # Get corresponding cell labels.
    ref_labels <- ref_col_data[ref_neighbors, ref_column_name]
    # Find modal cell label over a threshold.
    top_label <- which_mode(x=ref_labels, top_frac_threshold=top_frac_threshold, top_next_ratio_threshold=top_next_ratio_threshold)
  })
}


# Purpose: get means from nearest neighbor cells for continuous values.
get_nn_means <- function(query_data, query_search, ref_col_data, ref_column_name) {
  query_nns <- sapply(seq(1, nrow(query_data)), function(i) {
    # Get labels of neighboring cells in reference space.
    ref_neighbors <- query_search[['nn.idx']][i,]
    # Get corresponding reference cell label.
    ref_labels <- ref_col_data[ref_neighbors, ref_column_name]
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
#   o  test the ref_col_data for the required columns
#   o  the ref_column_name value must be in the colnames(colData(cds)).
#   o  keep k=10 or so for confidence in the label transfer
#

#' Transfer cell column data from a reference to a query
#' cell_data_set.
#'
#' @description For each cell in a query cell_data_set,
#' transfer_cell_labels finds the most similar cell data
#' in a reference cell_data_set and copies it to the query.
#'
#' transfer_cell_labels() requires a nearest neighbor index
#' made from a reference reduced dimension matrix, the reference
#' cell data to transfer, and a query cds. The index can be made
#' from UMAP coordinates using the build_nn_index=TRUE option in
#' the reduce_dimensions(..., build_nn_index=TRUE) function, for
#' example. The query cds must have been processed with the
#' preprocess_transform() and reduce_dimension_transform()
#' functions using the models created when the reference cds
#' was processed, rather than with preprocess_cds() and
#' reduce_dimension(). The models are made when the
#' reference cds is processed and must be saved to disk at
#' that time using save_transform_models(). The
#' load_transform_models() function loads the models into the
#' query cds where they can be used by preprocess_transform()
#' and reduce_dimension(). The cells in the reference and query
#' cdses must be similar in the sense that they map to similar
#' reduced dimension coordinates. For each cell in the query cds,
#' transfer_cell_labels() finds the k nearest neighbor reference
#' cells and, if at least top_frac_threshold fraction of them
#' have the same ref_column_name value, it copies the value to
#' the query_column_name column in the query cds. If not, it
#' checks whether the ratio of the most frequent to the second
#' most frequent value is at least top_next_ratio_threshold, in
#' which case it copies it. Otherwise, it sets the query value
#' to NA.
#'
#' Note: Monocle3 does not have an align_transform() function to
#' apply align_cds()-related transforms at this time. If your
#' data sets require batch correction, you need to co-embed them.
#'
#' @param cds_query the cell_data_set upon which to perform this operation
#' @param reduction_method a string specifying the reduced dimension
#'   matrix to use for the label transfer. These are "PCA", 
#'   "LSI", and "UMAP". Default is "UMAP".
#' @param ref_col_data the reference cds colData data frame, which
#'   is obtained using the colData(cds_ref) function.
#' @param ref_column_name a string giving the name of the reference
#'   cds column with the value to copy to the query cds.
#' @param query_column_name a string giving the name of the
#'   query cds column where you want the cell stored. The default
#'   is ref_column_name.
#' @param transform_models_dir a string giving the name of the transform
#'   model directory to load into the query cds. If it is NULL, use
#'   the transform models in the query cds, which requires that the
#'   reference models were loaded into the query cds before the
#'   call to transfer_cell_labels(). The default is NULL.
#' @param k an integer giving the number of reference nearest
#'   neighbors to find. This value must be large enough to find
#'   meaningful column value fractions. The default is 10.
#' @param nn_control a list of parameters used to make the nearest
#'   neighbors index. See the set_nn_control help for additional
#'   details. The default is to use the global nn_control list.
#' @param top_frac_threshold a numeric value giving the minimum
#'   value of the top fraction of reference values required
#'   for transfering the reference value to the query. The top
#'   fraction is the fraction of the k neighbors with the most
#'   frequent value. The default is 0.5.
#' @param top_next_ratio_threshold a numeric value giving the
#'   minimum value ratio of the counts of the most frequent to the
#'   second most frequent reference values. The default is 1.5.
#' @param verbose a boolean controlling verbose output.
#'
#' @return an updated cell_data_set object
#' @export
transfer_cell_labels <- function(cds_query,
                                 reduction_method=c('UMAP', 'PCA', 'LSI'),
                                 ref_col_data,
                                 ref_column_name,
                                 query_column_name=ref_column_name,
                                 transform_models_dir=NULL,
                                 k=10,
                                 nn_control=list(),
                                 top_frac_threshold=0.5,
                                 top_next_ratio_threshold=1.5,
                                 verbose=FALSE) {

  assertthat::assert_that(class(cds_query) == 'cell_data_set',
                          msg=paste0('cds_query parameter is not a cell_data_set'))
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be 'UMAP', 'PCA', or 'LSI'")

  assertthat::assert_that(assertthat::is.count(k))

  assertthat::assert_that(is.double(top_frac_threshold),
                          msg=paste0('top_frac_threshold value is not numeric'))
  assertthat::assert_that(is.double(top_next_ratio_threshold),
                          msg=paste0('top_next_ratio_threshold value is not numeric'))
  reduction_method <- match.arg(reduction_method)

  if(!is.data.frame(ref_col_data)) {
    ref_col_data <- as.data.frame(ref_col_data)
  }

  assertthat::assert_that(ref_column_name %in% colnames(ref_col_data),
                          msg=paste0('ref_column_name \'', ref_column_name, '\' is not in the ref_col_ata'))

  nn_control <- set_nn_control(mode=2,
                               nn_control=nn_control,
                               k=k,
                               nn_control_default=get_global_variable('nn_control_2'),
                               verbose=verbose)

  # Are the ref_column_name values discrete?
  label_data_are_discrete <- !is.double(ref_col_data[[ref_column_name]][1])

  # Load the reference projection models and nn indexes
  # into the query cds.
  if(!is.null(transform_models_dir)) {
    cds_query <- load_transform_models(cds=cds_query, directory_path=transform_models_dir)
  }

  assertthat::assert_that(!is.null(cds_query@reduce_dim_aux[[reduction_method]]),
                          msg=paste0("Reduction Method '", reduction_method, "' is not in the",
                                    "loaded model object."))


  # Search the reference reduction_method space for nearest neighbors
  # to the query cells.
  # The cds_reduced_dims contains the query cell coordinates
  # after projection into the reference space.
  # The cds@reduce_dim_aux[[reduction_method]] contains the reduction_method
  # coordinates for the reference data set, which were
  # loaded using load_transform_models() above.
  cds_reduced_dims <- reducedDims(cds_query)[[reduction_method]]
  cds_nn_index <- get_cds_nn_index(cds=cds_query, reduction_method=reduction_method, nn_control=nn_control, verbose=verbose)
  cds_res <- search_nn_index(query_matrix=cds_reduced_dims, nn_index=cds_nn_index,
                             k=k, nn_control=nn_control, verbose=verbose)
 
  # Get the best reference cell label for the query cells.
  if(label_data_are_discrete) {
    nn_labels <- get_nn_cell_label(query_data=cds_reduced_dims, 
                                   query_search=cds_res, 
                                   ref_col_data=ref_col_data, 
                                   ref_column_name=ref_column_name,
                                   top_frac_threshold=top_frac_threshold,
                                   top_next_ratio_threshold=top_next_ratio_threshold)
  } else {
    nn_labels <- get_nn_means(query_data=cds_reduced_dims,
                              query_search=cds_res,
                              ref_col_data=ref_col_data,
                              ref_column_name=ref_column_name)
  }
  
  colData(cds_query)[[query_column_name]] <- NULL
  colData(cds_query)[[query_column_name]] <- nn_labels

  return(cds_query)
}


# fill in NAs ----------------------------------------------------------------

edit_cell_label <- function(curr_label, other_labels, top_frac_threshold=0.5, top_next_ratio_threshold=1.5) {
  # Change NAs to strings.
  curr_label <- tidyr::replace_na(curr_label, "NA")
  other_labels <- tidyr::replace_na(other_labels, "NA")
  top_label <- which_mode(x=other_labels, top_frac_threshold=top_frac_threshold, top_next_ratio_threshold=top_next_ratio_threshold)
  
  # Switch back if necessary.
  top_label <- gsub("NA", NA_character_, top_label)
  
  return(top_label)
}


edit_query_cell_labels <- function(preproc_res,
                                   query_col_data,
                                   query_nn_index,
                                   column_name,
                                   k=10,
                                   nn_control=nn_control,
                                   top_frac_threshold=0.5,
                                   top_next_ratio_threshold=1.5,
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
    query_labels <- query_col_data[query_neighbors, column_name]
    curr_label <- query_labels[1]
    other_labels <- query_labels[2:k]
    top_label <- edit_cell_label(curr_label=curr_label,
                                 other_labels=other_labels,
                                 top_frac_threshold=top_frac_threshold,
                                 top_next_ratio_threshold=top_next_ratio_threshold)
  })
}


# Purpose: replace missing cell labels.
# Notes:
#   the column_name value must be in the colnames(colData(cds)).

#' Try to replace NA cell column data values.
#'
#' @description Try to replace NA values left in a query cds by
#' transfer_cell_labels().
#'
#' fix_missing_cell_labels() uses non-NA cell data values in the query cds
#' to replace NAs in nearby cells. It partitions the cells into a set with
#' NA and a set with non-NA column data values. It makes a nearest
#' neighbor index using cells with non-NA values, and for each cell
#' with NA, it tries to find an acceptable non-NA column data value as
#' follows. If at least top_frac_threshold fraction of them have the same
#' value, it replaces the NA with it. If not, it checks whether the ratio
#' of the most frequent to the second most frequent values is at least
#' top_next_ratio_threshold, in which case it copies the most frequent
#' value. Otherwise, it leaves the NA.
#'
#' @param cds_query the cell_data_set upon which to perform this operation
#' @param reduction_method a string specifying the reduced dimension
#'   matrix to use for the label transfer. These are "PCA",
#'   "LSI", and "UMAP". Default is "UMAP".
#' @param column_name a string giving the name of the query
#'   cds column to fix.
#' @param out_notna_models_dir a string giving the name of the transform
#'   model directory to which to save the not-NA transform models,
#'   which includes the nearest neighbor index.
#' @param k an integer giving the number of reference nearest
#'   neighbors to find. This value must be large enough to find
#'   meaningful column value fractions. The default is 10.
#' @param nn_control a list of parameters used to make the nearest
#'   neighbors index. See the set_nn_control help for additional
#'   details. The default is to use the global nn_control list.
#' @param top_frac_threshold a numeric value giving the minimum
#'   value of the top fraction of reference values required
#'   for transfering the reference value to the query. The top
#'   fraction is the fraction of the k neighbors with the most
#'   frequent value. The default is 0.5.
#' @param top_next_ratio_threshold a numeric value giving the
#'   minimum value ratio of the counts of the most frequent to the
#'   second most frequent reference values. The default is 1.5.
#' @param verbose a boolean controlling verbose output.
#'
#' @return an updated cell_data_set object
#' @export
fix_missing_cell_labels <- function(cds,
                                    reduction_method=c('UMAP', 'PCA', 'LSI'),
                                    column_name,
                                    out_notna_models_dir=NULL,
                                    k=10,
                                    nn_control=list(),
                                    top_frac_threshold=0.5,
                                    top_next_ratio_threshold=1.5,
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
  assertthat::assert_that(column_name %in% colnames(colData(cds)),
                          msg=paste0('column_name \'', column_name, '\' is not in the cds colData'))

  nn_control <- set_nn_control(mode=3,
                               nn_control=nn_control,
                               k=k,
                               nn_control_default=get_global_variable('nn_control_2'),
                               verbose=verbose)

  # Partition cds.
  notna_cds <- cds[, !is.na(colData(cds)[[column_name]])]
  
  # Build index on not NA cds, where there is a label.
  notna_nn_index <- make_nn_index(subject_matrix=reducedDims(notna_cds)[[reduction_method]],
                                  nn_control=nn_control,
                                  verbose=verbose)
  notna_col_data <- as.data.frame(colData(notna_cds))

  if(!is.null(out_notna_models_dir)) {
    save_transform_models(cds=notna_cds, directory_path=out_notna_models_dir)
  }
  
  na_cds <- cds[, is.na(colData(cds)[[column_name]])]
  new_cell_labels <- edit_query_cell_labels(preproc_res=reducedDims(na_cds)[[reduction_method]],
                                            query_col_data=notna_col_data, 
                                            query_nn_index=notna_nn_index, 
                                            column_name=column_name,
                                            k=k,
                                            nn_control=nn_control,
                                            top_frac_threshold=top_frac_threshold,
                                            top_next_ratio_threshold=top_next_ratio_threshold,
                                            verbose=verbose)
  
  na_label_index <- which(is.na(colData(cds)[[column_name]]))
  colData(cds)[[column_name]][na_label_index] <- new_cell_labels

  return(cds)
}

