# Test whether a matrix is one of our supported sparse matrices
is_sparse_matrix <- function(x){
  any(class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix", "CsparseMatrix"))
}

#' Function to calculate size factors for single-cell RNA-seq data
#'
#' @param cds The cell_data_set
#' @param round_exprs A logic flag to determine whether or not the expression
#'   value should be rounded
#' @param method A string to specify the size factor calculation approach.
#'   Options are "mean-geometric-mean-total" (default),
#'   "mean-geometric-mean-log-total".
#'
#' @return Updated cell_data_set object with a new colData column called
#'   'Size_Factor'.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     colData(cds)[['Size_Factor']] <- NULL
#'     cds <- estimate_size_factors(cds)
#'   }
#'
#' @export
estimate_size_factors <- function(cds,
                                  round_exprs=TRUE,
                                  method=c("mean-geometric-mean-total",
                                           'mean-geometric-mean-log-total'))
{
  method <- match.arg(method)
  if(is(SingleCellExperiment::counts(cds), 'IterableMatrix')) {
    if(any(BPCells::colSums(SingleCellExperiment::counts(cds)) == 0)) {
      warning("Your CDS object contains cells with zero reads. ",
                    "This causes size factor calculation to fail. Please remove ",
                    "the zero read cells using ",
                    "cds <- cds[,Matrix::colSums(exprs(cds)) != 0] and then ",
                    "run cds <- estimate_size_factors(cds)")
      return(cds)
    }
  }
  else {
    if(any(Matrix::colSums(SingleCellExperiment::counts(cds)) == 0)) {
      warning("Your CDS object contains cells with zero reads. ",
                    "This causes size factor calculation to fail. Please remove ",
                    "the zero read cells using ",
                    "cds <- cds[,Matrix::colSums(exprs(cds)) != 0] and then ",
                    "run cds <- estimate_size_factors(cds)")
      return(cds)
    }
  }

  if (is_sparse_matrix(SingleCellExperiment::counts(cds))){
    size_factors(cds) <- estimate_sf_sparse(SingleCellExperiment::counts(cds),
                                            round_exprs=round_exprs,
                                            method=method)
  } else if(is(SingleCellExperiment::counts(cds), 'IterableMatrix')) {
    size_factors(cds) <- estimate_sf_bpcells(SingleCellExperiment::counts(cds),
                                           round_exprs=round_exprs,
                                           method=method)
  }else{
    size_factors(cds) <- estimate_sf_dense(SingleCellExperiment::counts(cds),
                                           round_exprs=round_exprs,
                                           method=method)
  }
  return(cds)
}

# Estimate size factors for each column, given a sparseMatrix from the Matrix
# package
estimate_sf_sparse <- function(counts,
                               round_exprs=TRUE,
                               method="mean-geometric-mean-total"){
  if (round_exprs)
    counts <- round(counts)

  if(method == 'mean-geometric-mean-total') {
    cell_total <- Matrix::colSums(counts)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }else if(method == 'mean-geometric-mean-log-total') {
    cell_total <- Matrix::colSums(counts)
    sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

# Estimate size factors for each column, given an IterableMatrix (or
# derived class) from the BPCells package.
estimate_sf_bpcells <- function(counts,
                               round_exprs=TRUE,
                               method="mean-geometric-mean-total"){
  if (round_exprs)
    counts <- round(counts)

  if(method == 'mean-geometric-mean-total') {
    cell_total <- BPCells::colSums(counts)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }else if(method == 'mean-geometric-mean-log-total') {
    cell_total <- BPCells::colSums(counts)
    sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

# Estimate size factors for each column, given a matrix
estimate_sf_dense <- function(counts,
                              round_exprs=TRUE,
                              method="mean-geometric-mean-total"){

  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  if(method == "mean-geometric-mean-log-total") {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
  } else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

sparse_apply <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
  if (convert_to_dense){
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...)
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...)
      }, FUN, ...)
    }
  }else{
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...)
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...)
      }, FUN, ...)
    }
  }

  return(res)

}

# 20230424 bge: my test suggests that I can subset a BPCells matrix using
#               the subset operator with the form x[i, , drop = FALSE]
#' @noRd
split_rows <- function (x, ncl) {
  lapply(parallel::splitIndices(nrow(x), ncl),
         function(i) x[i, , drop = FALSE])
}

# 20230424 bge: my test suggests that I can subset a BPCells matrix using
#               the subset operator with the form x[i, , drop = FALSE]
#' @noRd
split_cols <- function (x, ncl) {
  lapply(parallel::splitIndices(ncol(x), ncl),
         function(i) x[, i, drop = FALSE])
}

#' @noRd
sparse_par_r_apply <- function (cl, x, FUN, convert_to_dense, ...) {
  par_res <- do.call(c, parallel::clusterApply(cl = cl,
                                               x = split_rows(x,
                                                              length(cl)),
                                               fun = sparse_apply, MARGIN = 1L,
                                               FUN = FUN,
                                               convert_to_dense=convert_to_dense, ...),
                     quote = TRUE)
  names(par_res) <- row.names(x)
  par_res
}

#' @noRd
sparse_par_c_apply <- function (cl = NULL, x, FUN, convert_to_dense, ...) {
  par_res <- do.call(c, parallel::clusterApply(cl = cl,
                                               x = split_cols(x,
                                                              length(cl)),
                                               fun = sparse_apply, MARGIN = 2L,
                                               FUN = FUN,
                                               convert_to_dense=convert_to_dense, ...),
                     quote = TRUE)
  names(par_res) <- colnames(x)
  par_res
}


#' Multicore apply-like function for cell_data_set
#'
#' mc_es_apply computes the row-wise or column-wise results of FUN, just like
#' esApply. Variables in colData from cds are available in FUN.
#'
#' @param cds A cell_data_set object.
#' @param MARGIN The margin to apply to, either 1 for rows (samples) or 2 for
#'   columns (features).
#' @param FUN Any function.
#' @param required_packages A list of packages FUN will need. Failing to
#'   provide packages needed by FUN will generate errors in worker threads.
#' @param convert_to_dense Whether to force conversion of a sparse matrix to a
#'   dense one before calling FUN.
#' @param reduction_method character, the method used to reduce dimension.
#'   Default "UMAP".
#' @param ... Additional parameters for FUN.
#' @param cores The number of cores to use for evaluation.
#'
#' @importFrom Biobase multiassign
#' @return The result of with(colData(cds) apply(counts(cds)), MARGIN, FUN, ...))
mc_es_apply <- function(cds, MARGIN, FUN, required_packages, cores=1,
                        convert_to_dense=TRUE,
                        reduction_method="UMAP", ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  coldata_df = as.data.frame(colData(cds))
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
  }, error = function(e) {} )

  tryCatch({
    coldata_df$pseudotime = pseudotime(cds)
  }, error = function(e) {} )

  Biobase::multiassign(names(as.data.frame(coldata_df)),
                       as.data.frame(coldata_df), envir=e1)
  environment(FUN) <- e1


  platform <- Sys.info()[['sysname']]

  # Temporarily disable OpenMP threading in functions to be run in parallel
  old_omp_num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
  if (is.na(old_omp_num_threads)){
    old_omp_num_threads = 1
  }
  RhpcBLASctl::omp_set_num_threads(1)

  # Temporarily set the number of threads the BLAS library can use to be 1
  old_blas_num_threads = as.numeric(Sys.getenv("OPENBLAS_NUM_THREADS"))
  if (is.na(old_omp_num_threads)){
    old_blas_num_threads = 1
  }
  RhpcBLASctl::blas_set_num_threads(1)

  # Note: use outfile argument to makeCluster for debugging
  if (platform == "Windows")
    cl <- parallel::makeCluster(cores)
  if (platform %in% c("Linux", "Darwin"))
    cl <- parallel::makeCluster(cores, type="FORK")

  cleanup <- function(){
    parallel::stopCluster(cl)
    RhpcBLASctl::omp_set_num_threads(old_omp_num_threads)
    RhpcBLASctl::blas_set_num_threads(old_blas_num_threads)
  }
  on.exit(cleanup)

  if (is.null(required_packages) == FALSE){
    parallel::clusterCall(cl, function(pkgs) {
      options(conflicts.policy =
                list(error = FALSE,
                     warn = FALSE,
                     generics.ok = TRUE,
                     can.mask = c("base", "methods", "utils",
                                  "grDevices", "graphics",
                                  "stats"),
                     depends.ok = TRUE))
      for (req in pkgs) {
        suppressMessages(library(req, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE,
                verbose=FALSE))
      }
    }, required_packages)
  }

  #
  # 20230424 bge: These two following calls using counts(cds) appear to work based on my simple test.
  #
  if (MARGIN == 1){
    suppressWarnings(res <- sparse_par_r_apply(cl=cl, x=SingleCellExperiment::counts(cds), FUN=FUN,
                                               convert_to_dense=convert_to_dense, ...))
  }else{
    suppressWarnings(res <- sparse_par_c_apply(cl=cl, x=SingleCellExperiment::counts(cds), FUN=FUN,
                                               convert_to_dense=convert_to_dense, ...))
  }

  res
}

#' @importFrom Biobase multiassign
smart_es_apply <- function(cds, MARGIN, FUN, convert_to_dense,
                           reduction_method="UMAP", ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  coldata_df = as.data.frame(colData(cds))
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    coldata_df$pseudotime = pseudotime(cds)
  }, error = function(e) {} )
  Biobase::multiassign(names(as.data.frame(coldata_df)),
                       as.data.frame(coldata_df), envir=e1)
  environment(FUN) <- e1

  if (is_sparse_matrix(SingleCellExperiment::counts(cds))){
    res <- sparse_apply(SingleCellExperiment::counts(cds), MARGIN, FUN, convert_to_dense, ...)
  } else {
    res <- pbapply::pbapply(SingleCellExperiment::counts(cds), MARGIN, FUN, ...)
  }

  if (MARGIN == 1)
  {
    names(res) <- row.names(cds)
  }else{
    names(res) <- colnames(cds)
  }

  res
}


#' Detects genes above minimum threshold.
#'
#' @description For each gene in a cell_data_set object, detect_genes counts
#' how many cells are expressed above a minimum threshold. In addition, for
#' each cell, detect_genes counts the number of genes above this threshold that
#' are detectable. Results are added as columns num_cells_expressed and
#' num_genes_expressed in the rowData and colData tables respectively.
#'
#' @param cds Input cell_data_set object.
#' @param min_expr Numeric indicating expression threshold
#' @return Updated cell_data_set object
#' @export
#' @examples
#' \dontrun{
#'    cds <- detect_genes(cds, min_expr=0.1)
#' }
detect_genes <- function(cds, min_expr=0){
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(is.numeric(min_expr))

  mat_bin <- counts(cds) > min_expr
  rowData(cds)$num_cells_expressed <- Matrix::rowSums(mat_bin > min_expr)
  colData(cds)$num_genes_expressed <- Matrix::colSums(mat_bin > min_expr)

  cds
}

#' Return a size-factor normalized and (optionally) log-transformed expression
#' matrix
#'
#' @param cds A CDS object to calculate normalized expression matrix from.
#' @param norm_method String indicating the normalization method. Options are
#'   "log" (Default), "binary" and "size_only".
#' @param pseudocount A pseudocount to add before log transformation. Ignored
#'   if norm_method is not "log". Default is 1.
#' @return Size-factor normalized, and optionally log-transformed, expression
#'   matrix.
#'
#' @examples
#'   \donttest{
#'     cds <- load_a549()
#'     normalized_matrix <- normalized_counts(cds)
#'   }
#'
#' @export
normalized_counts <- function(cds,
                              norm_method=c("log", "binary", "size_only"),
                              pseudocount=1) {
  norm_method <- match.arg(norm_method)

  norm_mat <- SingleCellExperiment::counts(cds)

  if (norm_method == "binary"){
    if(is(norm_mat, 'IterableMatrix')) {
      norm_mat <- BPCells::binarize(norm_mat, threshold=0, strict_inequality=TRUE)
#      stop('normalized_counts: binary normalization is unimplemented at this time...check back soon')
    }
    else {
      # The '+ 0' coerces the matrix to type numeric. It's possible
      # to use 'as.numeric(norm_mat > 0)' but the matrix
      # attributes disappear...
      norm_mat <- (norm_mat > 0) + 0
      if (is_sparse_matrix(norm_mat)) {
        norm_mat = methods::as(norm_mat, "dgCMatrix")
      }
    }
  }
  else {
    assertthat::assert_that(!is.null(size_factors(cds)))
    if(is(norm_mat, 'IterableMatrix')) {
      if(norm_method == 'log' && pseudocount != 1) {
        stop('normalized_counts: pseudocount must be 1 for sparse expression matrices and norm_method log')
      }
      norm_mat <- BPCells::t(BPCells::t(norm_mat) / size_factors(cds))
      if(norm_method == 'log' && pseudocount == 1) {
        norm_mat <- log1p(norm_mat) / log(10)
      }
    }
    else
    if (is_sparse_matrix(norm_mat)){
      norm_mat <- norm_mat
      norm_mat@x = norm_mat@x / rep.int(size_factors(cds), diff(norm_mat@p))
      if (norm_method == "log"){
        if (pseudocount == 1){
          norm_mat@x = log10(norm_mat@x + pseudocount)
        }else{
          stop("Pseudocount must equal 1 with sparse expression matrices")
        }
      }
    }else{
      norm_mat = Matrix::t(Matrix::t(norm_mat) / size_factors(cds))
      if (norm_method == "log"){
#          norm_mat@x <- log10(norm_mat + pseudocount)
          norm_mat <- log10(norm_mat + pseudocount)
      }
    }
  }

  return(norm_mat)
}



#' Combine a list of cell_data_set objects
#'
#' This function will combine a list of cell_data_set objects into a new
#' cell_data_set object.
#'
#' @details If any of the counts matrices is BPCells class, the combined
#'   counts matrix will be BPCells class.
#'
#' @param cds_list List of cds objects to be combined.
#' @param keep_all_genes Logical indicating what to do if there is a mismatch
#'   in the gene sets of the CDSs. If TRUE, all genes are kept and cells from
#'   CDSs missing a given gene will be filled in with zeroes. If FALSE, only
#'   the genes in common among all of the CDSs will be kept. Default is TRUE.
#' @param cell_names_unique Logical indicating whether all of the cell IDs
#'   across all of the CDSs are unique. If FALSE, the CDS name is appended to
#'   each cell ID to prevent collisions. These cell IDs are used as count matrix
#'   column names and colData(cds) row names. Cell names stored in other
#'   cds locations are not modified so you will need to modify them manually
#'   for consistency. Default is FALSE.
#' @param sample_col_name A string to be the column name for the colData column
#'   that indicates which original cds the cell derives from. Default is
#'   "sample".
#' @param keep_reduced_dims Logical indicating whether to keep the reduced
#'   dimension matrices. Do not keep the reduced dimensions unless you know
#'   that the reduced dimensions are the same in each CDS. This is true for
#'   projected data sets, for example. Default is FALSE.
#'
#' @return A combined cell_data_set object.
#' @export
#'
combine_cds <- function(cds_list,
                        keep_all_genes = TRUE,
                        cell_names_unique = FALSE,
                        sample_col_name = "sample",
                        keep_reduced_dims = FALSE,
                        matrix_control = list(),
                        diagnostics=FALSE) {

  if(diagnostics) { message('Diag: here 1') }

  assertthat::assert_that(is.list(cds_list),
                          msg=paste("cds_list must be a list."))

  assertthat::assert_that(all(sapply(cds_list, class) == "cell_data_set"),
                          msg=paste("All members of cds_list must be",
                                    "cell_data_set class."))
  assertthat::assert_that(is.character(sample_col_name))

  if (sample_col_name == "sample" &
      any(sapply(cds_list, function(cds) "sample" %in% names(colData(cds))))) {
    warning("By default, the combine_cds function adds a column called ",
                   "'sample' which indicates which initial cds a cell came ",
                   "from. One or more of your input cds objects contains a ",
                   "'sample' column, which will be overwritten. We recommend ",
                   "you rename this column or provide an alternative column ",
                   "name using the 'sample_col_name' parameter.")
  }
  assertthat::assert_that(!any(sapply(cds_list, function(cds)
    sum(is.na(names(colData(cds)))) != 0)),
                          msg = paste0("One of the input CDS' has a colData ",
                                       "column name that is NA, please ",
                                       "remove or rename that column before ",
                                       "proceeding."))
  assertthat::assert_that(!any(sapply(cds_list, function(cds)
    sum(is.na(names(rowData(cds)))) != 0)),
    msg = paste0("One of the input CDS' has a rowData ",
                 "column name that is NA, please ",
                 "remove or rename that column before ",
                 "proceeding."))

  if(diagnostics) { message('Diag: here 2') }

  num_cells <- sapply(cds_list, ncol)
  if(sum(num_cells == 0) != 0) {
    message("Some CDS' have no cells, these will be skipped.")
    cds_list <- cds_list[num_cells != 0]
  }
  if(length(cds_list) == 1) return(cds_list[[1]])

  assertthat::assert_that(is.logical(keep_all_genes))
  assertthat::assert_that(is.logical(cell_names_unique))

  list_named <- TRUE
  if(is.null(names(cds_list))) {
    list_named <- FALSE
  }

  if(diagnostics) { message('Diag: here 3') }

  if(!is.null(matrix_control[['matrix_class']]) &&
     matrix_control[['matrix_class']] == 'BPCells') {
    bpcells_matrix_flag <- TRUE
  }
  else {
    bpcells_matrix_flag <- FALSE
    # Are any of the count matrices BPCells class?
    for(i in seq(1, length(cds_list), 1)) {
      if(is(counts(cds_list[[i]]), 'IterableMatrix')) {
        bpcells_matrix_flag <- TRUE
        break
      }
    }
  }

  if(diagnostics) { message('Diag: here 4') }

  check_matrix_control(matrix_control=matrix_control, control_type='counts', check_conditional=FALSE)
  if(bpcells_matrix_flag ||
     (!is.null(matrix_control[['matrix_class']]) && matrix_control[['matrix_class']] == 'BPCells')) {
    matrix_control_default <- get_global_variable('matrix_control_bpcells_counts')
  }
  else {
    matrix_control_default <- get_global_variable('matrix_control_csparsematrix_counts')
  }
  matrix_control <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='unrestricted')

  if(diagnostics) { message('Diag: here 5') }

  exprs_list <- list()
  fd_list <- list()
  pd_list <- list()
  gene_list <- c()
  overlap_list <- c(row.names(fData(cds_list[[1]])))
  pdata_cols <- c()
  fdata_cols <- c()
  all_cells <- c()

  if(diagnostics) { message('Diag: here 6') }

  for(cds in cds_list) {
    # Make a vector of gene names either all names or
    # only names in common to all CDSes.
    gene_list <-  c(gene_list, row.names(fData(cds)))
    overlap_list <- intersect(overlap_list, row.names(fData(cds)))
    if (!keep_all_genes) {
      gene_list <- overlap_list
    }

  if(diagnostics) { message('Diag: here 7') }

    # Make concatenated vectors of column (header) names of the cells
    # and features, and of cell names.
    pdata_cols <- c(pdata_cols, names(pData(cds)))
    fdata_cols <- c(fdata_cols, names(fData(cds)))
    all_cells <- c(all_cells, row.names(pData(cds)))
  }

  if(diagnostics) { message('Diag: here 8') }

  # Remove duplicate gene names, feature and cell column
  # names, and cell names.
  gene_list <- unique(gene_list)
  if(length(overlap_list) == 0) {
    if (keep_all_genes) {
      warning("No genes are shared amongst all the CDS objects.")
    } else {
      stop("No genes are shared amongst all the CDS objects. To generate ",
                 "a combined CDS with all genes, use keep_all_genes = TRUE")
    }
  }

  if(diagnostics) { message('Diag: here 9') }

  pdata_cols <- unique(pdata_cols)
  fdata_cols <- unique(fdata_cols)
  if (sum(duplicated(all_cells)) != 0 & cell_names_unique) {
    stop("Cell names are not unique across CDSs - cell_names_unique ",
               "must be FALSE.")
  }
  all_cells <- unique(all_cells)

  if(diagnostics) { message('Diag: here 10') }

  # Give all CDSes the same set of rows (amongst other things),
  # looping through the CDSes.
  for(i in seq(1, length(cds_list), 1)) {
    pd <- as.data.frame(pData(cds_list[[i]]))

  if(diagnostics) { message('Diag: here 11') }

    # Counts matrix rows of genes common to the CDSes examined
    # up to this pass through the loop.
    exp <- exprs(cds_list[[i]])
    if(bpcells_matrix_flag && !is(exp, 'IterableMatrix')) {
      exp <- as(exp, 'IterableMatrix')
    }
    exp <- exp[intersect(row.names(exp), gene_list),, drop=FALSE]

  if(diagnostics) { message('Diag: here 12') }

    # Make cell names distinct, if necessary, assign cell names to
    # pd, the sample names to a column in pd, and cell names to
    # the counts matrix columns.
    if (!cell_names_unique) {
      if(list_named) {
        row.names(pd) <- paste(row.names(pd), names(cds_list)[[i]], sep="_")
        pd[,sample_col_name] <- names(cds_list)[[i]]
      } else {
        row.names(pd) <- paste(row.names(pd), i, sep="_")
        pd[,sample_col_name] <- i
      }
      colnames(exp) <- row.names(pd)
    } else {
      if(list_named) {
        pd[,sample_col_name] <- names(cds_list)[[i]]
      } else {
        pd[,sample_col_name] <- i
      }
    }

  if(diagnostics) { message('Diag: here 13') }

    # Initialize new entries in pd to 'NA'.
    not_in <- pdata_cols[!pdata_cols %in% names(pd)]
    for (n in not_in) {
      pd[,n] <- NA
    }

  if(diagnostics) { message('Diag: here 14') }

    # Select feature data frame rows that are common to
    # fd and gene_list.
    fd <- as.data.frame(fData(cds_list[[i]]))
    fd <- fd[intersect(row.names(fd), gene_list),, drop=FALSE]

    # Make a vector of gene names that are that are not in
    # fd.
    not_in <- fdata_cols[!fdata_cols %in% names(fd)]
    for(col in names(fd)) {
      if(methods::is(fd[,col], "factor")) {
        fd[,col] <- as.character(fd[,col])
      }
    }
    for (n in not_in) {
      fd[,n] <- NA
    }
    not_in_g <- gene_list[!gene_list %in% row.names(fd)]

  if(diagnostics) { message('Diag: here 15') }

    # Make an empty matrix (and fd data frame) with the rows
    # that need to be added to the counts matrix for cds_list[[i]],
    # and append it to the accumulating counts matrix.
    if (length(not_in_g) > 0) {

  if(diagnostics) { message('Diag: here 16') }

      not_in_g_df <- as.data.frame(matrix(NA, nrow = length(not_in_g), ncol=ncol(fd)))
      row.names(not_in_g_df) <- not_in_g
      names(not_in_g_df) <- names(fd)
      fd <- rbind(fd, not_in_g_df)      # bge rbind

      extra_rows <- Matrix::Matrix(0, ncol=ncol(exp),
                                   sparse=TRUE,
                                   nrow=length(not_in_g))
      row.names(extra_rows) <- not_in_g
      colnames(extra_rows) <- colnames(exp)

  if(diagnostics) { message('Diag: here 17') }

      # Append additional rows.
      if(bpcells_matrix_flag) {
        exp <- rbind2(exp, as(extra_rows, 'IterableMatrix'))
      }
      else {
        exp <- rbind(exp, extra_rows)
      }
#      exp <- exp

  if(diagnostics) { message('Diag: here 18') }

    }

  if(diagnostics) { message('Diag: here 19') }

    # Gather matrices and data frames into lists.
    exprs_list[[i]] <- exp[gene_list, , drop=FALSE]
    fd_list[[i]] <- fd[gene_list, , drop=FALSE]
    pd_list[[i]] <- pd
  }

  if(diagnostics) { message('Diag: here 20') }

  all_fd <- array(NA,dim(fd_list[[1]]),dimnames(fd_list[[1]]))

  if(diagnostics) { message('Diag: here 21') }

  for (fd in fd_list) {
    for (j in colnames(fd)) {
      col_info <- all_fd[,j]
      col_info[is.na(col_info)] <- fd[is.na(col_info),j]
      col_info[col_info != fd[,j]] <- "conf"
      all_fd[,j] <- col_info
    }
  }

  if(diagnostics) { message('Diag: here 22') }

  confs <- sum(all_fd == "conf", na.rm=TRUE)

  if (confs > 0) {
   warning("When combining rowData, conflicting values were found - ",
                  "conflicts will be labelled 'conf' in the combined cds ",
                  "to prevent conflicts, either change conflicting values to ",
                  "match, or rename columns from different cds' to be unique.")
  }
  #all_fd <- do.call(cbind, fd_list)
  all_fd <- all_fd[,fdata_cols, drop=FALSE]

  if(diagnostics) { message('Diag: here 23') }

  # Build the final, comprehensive cell data frame and counts matrix.
  all_pd <- do.call(rbind, pd_list)    # bge rbind
  if(bpcells_matrix_flag) {
    num_exprs <- length(exprs_list)
    all_exp <- exprs_list[[1]]
    for(i in seq(2, num_exprs, 1)) {
      all_exp <- cbind2(all_exp, exprs_list[[i]])
    }
  }
  else {
    all_exp <- do.call(cbind, exprs_list)
  }

  if(diagnostics) { message('Diag: here 24') }

  # Filter counts matrix by fd and pd names.
  all_exp <- all_exp[row.names(all_fd), row.names(all_pd), drop=FALSE]

  if(diagnostics) { message('Diag: here 25') }

  # Make a BPCells count matrix, if necessary.
  if(bpcells_matrix_flag) {
    all_exp <- set_matrix_class(mat=all_exp, matrix_control=matrix_control)
  }

  if(diagnostics) { message('Diag: here 26') }

  # Make a combined CDS from all_exp, all_pd, and all_fd.
  new_cds <- new_cell_data_set(all_exp, cell_metadata = all_pd, gene_metadata = all_fd)

  if(diagnostics) { message('Diag: here 27') }

  # Add in preprocessing results.
  if(keep_reduced_dims) {

  if(diagnostics) { message('Diag: here 28') }

    # Find intersection of reduced dim names, for example, 'PCA', 'UMAP', 'Aligned'.
    reduced_dim_names <- names(reducedDims(cds_list[[1]]))
    for(i in seq(2, length(cds_list), 1)) {
      reduced_dim_names <- intersect(reduced_dim_names, names(reducedDims(cds_list[[i]])))
    }

  if(diagnostics) { message('Diag: here 29') }

#    for(red_dim in names(SingleCellExperiment::reducedDims(cds_list[[1]]))) {
    for(red_dim in reduced_dim_names) {
      reduced_dims_list <- list()
      for(j in seq(1, length(cds_list), 1)) {
        reduced_dims_list[[j]] <- SingleCellExperiment::reducedDims(cds_list[[j]])[[red_dim]]
      }
      SingleCellExperiment::reducedDims(new_cds, withDimnames=FALSE)[[red_dim]] <- do.call(rbind, reduced_dims_list, quote=FALSE)
      # The following should not happen; the accessor appears to ensure the
      # correct row order.
      if(!identical(rownames(SingleCellExperiment::reducedDims(new_cds)[[red_dim]]), rownames(all_pd))) {
        stop('Mis-ordered reduced matrix rows.')
      }
    }

  if(diagnostics) { message('Diag: here 30') }

  }

  # Add a BPCells row-major order matrix to assays
  # for BPCells count matrices.
  if(bpcells_matrix_flag) {
    new_cds <- set_cds_row_order_matrix(new_cds)
  }

  if(diagnostics) { message('Diag: here 31') }

  matrix_id <-  get_unique_id(counts(cds))
  new_cds <- initialize_counts_metadata(new_cds) 
  new_cds <- set_counts_identity(new_cds, 'combin_cds', matrix_id)

  if(diagnostics) { message('Diag: here 32') }

  new_cds
}


#' Clear CDS slots
#'
#' Function to clear all CDS slots besides colData, rowData and expression data.
#'
#' @param cds cell_data_set to be cleared
#'
#' @return A cell_data_set with only expression, rowData and colData present.
#' @export
#'
clear_cds_slots <- function(cds) {
  cds@reduce_dim_aux <- S4Vectors::SimpleList()
  cds@principal_graph_aux <- S4Vectors::SimpleList()
  cds@principal_graph <- S4Vectors::SimpleList()
  cds@clusters <- S4Vectors::SimpleList()
  SingleCellExperiment::reducedDims(cds) <- S4Vectors::SimpleList()
  cds
}


add_citation <- function(cds, citation_key) {
  citation_map <- list(
    UMAP = c("UMAP", "McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. Preprint at https://arxiv.org/abs/1802.03426 (2018)."),
    MNN_correct = c("MNN Correct", "Haghverdi, L. et. al. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat. Biotechnol. 36, 421-427 (2018). https://doi.org/10.1038/nbt.4091"),
    partitions = c("paritioning", c("Levine, J. H., et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.047",
                                  "Wolf, F. A. et. al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. Genome Biol. 20, 59 (2019). https://doi.org/10.1186/s13059-019-1663-x")),
    clusters = c("clustering", "Levine, J. H. et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.047"),
    leiden = c("leiden", "Traag, V.A., Waltman, L. & van Eck, N.J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reportsvolume 9, Article number: 5233 (2019). https://doi.org/10.1038/s41598-019-41695-z" )
  )
  if (is.null(S4Vectors::metadata(cds)$citations) | citation_key == "Monocle") {
    S4Vectors::metadata(cds)$citations <- data.frame(method = c("Monocle", "Monocle", "Monocle"),
                                          citations = c("Trapnell C. et. al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat. Biotechnol. 32, 381-386 (2014). https://doi.org/10.1038/nbt.2859",
                                                        "Qiu, X. et. al. Reversed graph embedding resolves complex single-cell trajectories. Nat. Methods 14, 979-982 (2017). https://doi.org/10.1038/nmeth.4402",
                                                        "Cao, J. et. al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496-502 (2019). https://doi.org/10.1038/s41586-019-0969-x"))
  }
  S4Vectors::metadata(cds)$citations <- rbind(S4Vectors::metadata(cds)$citations,
                                   data.frame(method = citation_map[[citation_key]][1],
                                              citations = citation_map[[citation_key]][2]))
  cds
}

#' Access citations for methods used during analysis.
#'
#' @param cds The cds object to access citations from.
#'
#' @return A data frame with the methods used and the papers to be cited.
#' @export
#'
#' @examples {
#'   \dontrun{
#'      get_citations(cds)
#'   }
#' }
get_citations <- function(cds) {
  message("Your analysis used methods from the following recent work. ",
                "Please cite them wherever you are presenting your analyses.")
  if(is.null(S4Vectors::metadata(cds)$citations)) {
    cds <- add_citation(cds, "Monocle")
  }
  S4Vectors::metadata(cds)$citations
}


# Make a unique identifier string.
get_unique_id <- function(object=NULL) {
  #
  # I don't have a way to calculate a checksum
  # without creating an in memory matrix copy
  # so skip if this is a BPCells matrix.
  if(is(object, 'IterableMatrix')) {
    return('BPcells matrix')
  }

  if(!is.null(object)) {
    object_dim <- dim(object)
    if(!is(object, 'IterableMatrix')) {
      object_checksum <- digest::digest(object)
    }
    else {
      object_checksum <- digest::digest(as(object, 'dgCMatrix'))
    }
    if(!is.null(object_dim))
      object_id <- list(checksum=object_checksum, dim=object_dim)
    else
      object_id <- list(checksum=object_checksum, dim=length(object))
  }
  else {
    id_count <- get_global_variable('id_count')
    rtime <- as.numeric(Sys.time()) * 100000 + id_count
    object_id <- openssl::md5(as.character(rtime))
    id_count <- id_count + 1
    set_global_variable('id_count', id_count)
  }

  return(object_id)
}


# Get current time stamp.
get_time_stamp <- function() {
  time_stamp <- format(Sys.time(), "%Y%m%d:%H%M%S" )
  return(time_stamp)
}


# Manage parallel processing for matrix multiplication.
matrix_multiply_multicore <- function(mat_a, mat_b, cores=1L) {
  assertthat::assert_that(is.matrix(mat_a) || is_sparse_matrix(mat_a) || methods::is(mat_a, 'DelayedMatrix'),
    msg=paste0('mat_a must be either a matrix or a sparse matrix'))
  assertthat::assert_that(is.matrix(mat_b) || is_sparse_matrix(mat_b) || methods::is(mat_b, 'DelayedMatrix'),
    msg=paste0('mat_b must be either a matrix or a sparse matrix'))

  if(cores > 1) {
    omp_num_threads <- get_global_variable('omp_num_threads')
    blas_num_threads <- get_global_variable('blas_num_threads')
  
    RhpcBLASctl::omp_set_num_threads(1L)
    RhpcBLASctl::blas_set_num_threads(1L)
  
    DelayedArray::setAutoBPPARAM(BPPARAM=BiocParallel::MulticoreParam(workers=as.integer(cores)))
  
    mat_c <- mat_a %*% mat_b
  
    DelayedArray::setAutoBPPARAM(BPPARAM=BiocParallel::SerialParam())
  
    RhpcBLASctl::omp_set_num_threads(as.integer(omp_num_threads))
    RhpcBLASctl::blas_set_num_threads(as.integer(blas_num_threads))
  } else {
    mat_c <- mat_a %*% mat_b
  }

  return(mat_c)
}


# Return the call stack as a character vector.
get_call_stack <- function () {
  cv<-as.vector(sys.calls())
  lcv <- length(cv)
  n <- lcv - 1

  ocv <- vector()
  for(i in seq(1,n,1)) {
    elem <- stringr::str_split(as.character(cv[i]), '[(]', n=2)[[1]][[1]]
    ocv <- c(ocv, elem)
  }

  return(ocv)
}


# Return the call stack as a character string.
get_call_stack_as_string <- function() {
  cs <- get_call_stack()
  scs <- ''
  for(i in seq(length(cs)-1)) {
    csep <- ifelse(i == 1, '', ' => ')
    scs <- sprintf("%s%s%s()", scs, csep, cs[[i]])
  }

  return(scs)
}


# Return the name of object as a string
object_name_to_string <- function( object ) {
  str <- deparse(substitute(object))

  return( str )
}


stop_no_noise <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


# number of tasks per block for multicore processing
# task vector limits
#   vbeg <- c(0,cumsum(tasks_per_block(11,3)))[1:3]+1
#   vend <- cumsum(tasks_per_block(11,3))
tasks_per_block <- function(ntask=NULL, nblock=NULL) {
  tasks_block <- rep(trunc(ntask/nblock), nblock)
  remain <- ntask %% nblock
  if(remain)
    for(i in seq(remain)) 
      tasks_block[i] <- tasks_block[i] + 1
  return(tasks_block)
}


#
# Initialize Monocle3 timer.
#
tick <- function(msg="") {
  set_global_variable('monocle3_timer_t0', Sys.time())
  set_global_variable('monocle3_timer_msg', msg)
}

#
# Return time elapsed since call to tick.
#
tock <- function() {
  t1 <- Sys.time()
  t0 <- get_global_variable('monocle3_timer_t0')
  msg <- get_global_variable('monocle3_timer_msg')
  if(length(msg) > 0) {
    message(sprintf('%s %.2f seconds.',msg, t1 - t0))
  }
  else {
    return(t1 - t0)
  }
}
