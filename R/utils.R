# Test whether a matrix is one of our supported sparse matrices
is_sparse_matrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix")
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
#' @export
estimate_size_factors <- function(cds,
                                  round_exprs=TRUE,
                                  method=c("mean-geometric-mean-total",
                                           'mean-geometric-mean-log-total'))
{
  method <- match.arg(method)
  if(any(Matrix::colSums(SingleCellExperiment::counts(cds)) == 0)) {
    warning(paste("Your CDS object contains cells with zero reads.",
                  "This causes size factor calculation to fail. Please remove",
                  "the zero read cells using",
                  "cds <- cds[,Matrix::colSums(exprs(cds)) != 0] and then",
                  "run cds <- estimate_size_factors(cds)"))
    return(cds)
  }
  if (is_sparse_matrix(SingleCellExperiment::counts(cds))){
    size_factors(cds) <- estimate_sf_sparse(SingleCellExperiment::counts(cds),
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

#' @keywords internal
split_rows <- function (x, ncl) {
  lapply(parallel::splitIndices(nrow(x), ncl),
         function(i) x[i, , drop = FALSE])
}

#' @keywords internal
split_cols <- function (x, ncl) {
  lapply(parallel::splitIndices(ncol(x), ncl),
         function(i) x[, i, drop = FALSE])
}

#' @keywords internal
sparse_par_r_apply <- function (cl, x, FUN, convert_to_dense, ...) {
  par_res <- do.call(c, BiocGenerics::clusterApply(cl = cl,
                                                   x = split_rows(x,
                                                                  length(cl)),
                                     fun = sparse_apply, MARGIN = 1L,
                                     FUN = FUN,
                                     convert_to_dense=convert_to_dense, ...),
                     quote = TRUE)
  names(par_res) <- row.names(x)
  par_res
}

#' @keywords internal
sparse_par_c_apply <- function (cl = NULL, x, FUN, convert_to_dense, ...) {
  par_res <- do.call(c, BiocGenerics::clusterApply(cl = cl,
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
    BiocGenerics::clusterCall(cl, function(pkgs) {
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

  if (MARGIN == 1){
    suppressWarnings(res <- sparse_par_r_apply(cl, SingleCellExperiment::counts(cds), FUN,
                                               convert_to_dense, ...))
  }else{
    suppressWarnings(res <- sparse_par_c_apply(cl, SingleCellExperiment::counts(cds), FUN,
                                               convert_to_dense, ...))
  }

  res
}

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


#' Principal Components Analysis
#'
#' Efficient computation of a truncated principal components analysis of a
#' given data matrix using an implicitly restarted Lanczos method from the
#' \code{\link{irlba}} package.
#'
#' @param x a numeric or complex matrix (or data frame) which provides the data
#'   for the principal components analysis.
#' @param retx a logical value indicating whether the rotated variables should
#'   be returned.
#' @param center a logical value indicating whether the variables should be
#'   shifted to be zero centered. Alternately, a centering vector of length
#'   equal the number of columns of \code{x} can be supplied.
#' @param scale. a logical value indicating whether the variables should be
#'   scaled to have unit variance before the analysis takes place. The default
#'   is \code{FALSE} for consistency with S, but scaling is often advisable.
#'   Alternatively, a vector of length equal the number of columns of \code{x}
#'   can be supplied.
#'
#'   The value of \code{scale} determines how column scaling is performed
#'   (after centering). If \code{scale} is a numeric vector with length equal
#'   to the number of columns of \code{x}, then each column of \code{x} is
#'   divided by the corresponding value from \code{scale}. If \code{scale} is
#'   \code{TRUE} then scaling is done by dividing the (centered) columns of
#'   \code{x} by their standard deviations if \code{center=TRUE}, and the root
#'   mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#'   See \code{\link{scale}} for more details.
#' @param n integer number of principal component vectors to return, must be
#'   less than \code{min(dim(x))}.
#' @param ... additional arguments passed to \code{\link{irlba}}.
#'
#' @return
#' A list with class "prcomp" containing the following components:
#' \itemize{
#'    \item{sdev} {the standard deviations of the principal components (i.e.,
#'      the square roots of the eigenvalues of the covariance/correlation
#'      matrix, though the calculation is actually done with the singular
#'      values of the data matrix).}
#'   \item{rotation} {the matrix of variable loadings (i.e., a matrix whose
#'     columns contain the eigenvectors).}
#'   \item {x} {if \code{retx} is \code{TRUE} the value of the rotated data
#'     (the centered (and scaled if requested) data multiplied by the
#'     \code{rotation} matrix) is returned. Hence, \code{cov(x)} is the
#'     diagonal matrix \code{diag(sdev^2)}.}
#'   \item{center, scale} {the centering and scaling used, or \code{FALSE}.}
#' }
#'
#' @note
#' The signs of the columns of the rotation matrix are arbitrary, and so may
#' differ between different programs for PCA, and even between different builds
#' of R.
#'
#' NOTE DIFFERENCES WITH THE DEFAULT \code{\link{prcomp}} FUNCTION! The
#' \code{tol} truncation argument found in \code{prcomp} is not supported. In
#' place of the truncation tolerance in the original function, the
#' \code{prcomp_irlba}  function has the argument \code{n} explicitly giving
#' the number of principal components to return. A warning is generated if the
#' argument \code{tol} is used, which is interpreted differently between the
#' two functions.
#'
#' @examples
#' \dontrun{set.seed(1)
#' x  <- matrix(rnorm(200), nrow=20)
#' p1 <- prcomp_irlba(x, n=3)
#' summary(p1)
#'
#' # Compare with
#' p2 <- prcomp(x, tol=0.7)
#' summary(p2)}
#'
#' @seealso \code{\link{prcomp}}
sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE,
                                scale. = FALSE, ...)
{
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`
            function to control that algorithm's convergence tolerance. See
            `?prcomp_irlba` for help.")
  orig_x <- x
  if (class(x) != "DelayedMatrix")
    x = DelayedArray::DelayedArray(x)

  args <- list(A=orig_x, nv=n)
  if (is.logical(center))
  {
    if (center) args$center <- DelayedMatrixStats::colMeans2(x)
  } else args$center <- center
  if (is.logical(scale.))
  {
    if (is.numeric(args$center))
    {
      scale. <- sqrt(DelayedMatrixStats::colVars(x))
      if (ans$scale) ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale. ^ 2)
    } else
    {
      if (ans$scale)
      {
        scale. <-
          sqrt(DelayedMatrixStats::colSums2(x ^ 2) / (max(1, nrow(x) - 1L)))
        ans$totalvar <-
          sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) /
                     (nrow(x) - 1L)))
      } else
      {
        ans$totalvar <-
          sum(DelayedMatrixStats::colSums2(x ^ 2) / (nrow(x) - 1L))
      }
    }
    if (ans$scale) args$scale <- scale.
  } else
  {
    args$scale <- scale.
    ans$totalvar <-
      sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) /
                 (nrow(x) - 1L)))
  }
  if (!missing(...)) args <- c(args, list(...))

  s <- do.call(irlba::irlba, args=args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  ans$svd_scale <- args$scale
  if (retx)
  {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN=`*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
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
#' cds <- detect_genes(cds, min_expr=0.1)
#' }
detect_genes <- function(cds, min_expr=0){
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(is.numeric(min_expr))

  rowData(cds)$num_cells_expressed <- Matrix::rowSums(SingleCellExperiment::counts(cds) > min_expr)
  colData(cds)$num_genes_expressed <- Matrix::colSums(SingleCellExperiment::counts(cds) > min_expr)

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
#'
#' @export
normalized_counts <- function(cds,
                              norm_method=c("log", "binary", "size_only"),
                              pseudocount=1){
  norm_method = match.arg(norm_method)
  norm_mat = SingleCellExperiment::counts(cds)
  if (norm_method == "binary"){
    norm_mat = norm_mat > 0
    if (is_sparse_matrix(norm_mat)){
      norm_mat = methods::as(norm_mat, "dgCMatrix")
    }
  }
  else {
    if (is_sparse_matrix(norm_mat)){
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
          norm_mat@x <- log10(norm_mat + pseudocount)
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
#'
#' @return A combined cell_data_set object.
#' @export
#'
combine_cds <- function(cds_list,
                        keep_all_genes = TRUE,
                        cell_names_unique = FALSE,
                        sample_col_name = "sample",
                        keep_reduced_dims = FALSE) {

  assertthat::assert_that(is.list(cds_list),
                          msg=paste("cds_list must be a list."))

  assertthat::assert_that(all(sapply(cds_list, class) == "cell_data_set"),
                          msg=paste("All members of cds_list must be",
                                    "cell_data_set class."))
  assertthat::assert_that(is.character(sample_col_name))

  if (sample_col_name == "sample" &
      any(sapply(cds_list, function(cds) "sample" %in% names(colData(cds))))) {
    warning(paste0("By default, the combine_cds function adds a column called ",
                   "'sample' which indicates which initial cds a cell came ",
                   "from. One or more of your input cds objects contains a ",
                   "'sample' column, which will be overwritten. We recommend ",
                   "you rename this column or provide an alternative column ",
                   "name using the 'sample_col_name' parameter."))
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

  exprs_list <- list()
  fd_list <- list()
  pd_list <- list()
  gene_list <- c()
  overlap_list <- c(row.names(fData(cds_list[[1]])))
  pdata_cols <- c()
  fdata_cols <- c()
  all_cells <- c()

  for(cds in cds_list) {
    gene_list <-  c(gene_list, row.names(fData(cds)))
    overlap_list <- intersect(overlap_list, row.names(fData(cds)))
    if (!keep_all_genes) {
      gene_list <- overlap_list
    }

    pdata_cols <- c(pdata_cols, names(pData(cds)))
    fdata_cols <- c(fdata_cols, names(fData(cds)))
    all_cells <- c(all_cells, row.names(pData(cds)))
  }

  gene_list <- unique(gene_list)
  if(length(overlap_list) == 0) {
    if (keep_all_genes) {
      warning(paste("No genes are shared amongst all the CDS objects."))
    } else {
      stop(paste("No genes are shared amongst all the CDS objects. To generate",
                 "a combined CDS with all genes, use keep_all_genes = TRUE"))
    }
  }
  pdata_cols <- unique(pdata_cols)
  fdata_cols <- unique(fdata_cols)
  if (sum(duplicated(all_cells)) != 0 & cell_names_unique) {
    stop(paste("Cell names are not unique across CDSs - cell_names_unique",
               "must be FALSE."))
  }
  all_cells <- unique(all_cells)
  for(i in 1:length(cds_list)) {
    pd <- as.data.frame(pData(cds_list[[i]]))
    exp <- exprs(cds_list[[i]])
    exp <- exp[intersect(row.names(exp), gene_list),, drop=FALSE]
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
    not_in <- pdata_cols[!pdata_cols %in% names(pd)]
    for (n in not_in) {
      pd[,n] <- NA
    }

    fd <- as.data.frame(fData(cds_list[[i]]))
    fd <- fd[intersect(row.names(fd), gene_list),, drop=FALSE]
    not_in <- fdata_cols[!fdata_cols %in% names(fd)]
    for(col in names(fd)) {
      if(class(fd[,col]) == "factor") {
        fd[,col] <- as.character(fd[,col])
      }
    }
    for (n in not_in) {
      fd[,n] <- NA
    }
    not_in_g <- gene_list[!gene_list %in% row.names(fd)]

    if (length(not_in_g) > 0) {
      not_in_g_df <- as.data.frame(matrix(NA, nrow = length(not_in_g), ncol=ncol(fd)))
      row.names(not_in_g_df) <- not_in_g
      names(not_in_g_df) <- names(fd)
      fd <- rbind(fd, not_in_g_df)

      extra_rows <- Matrix::Matrix(0, ncol=ncol(exp),
                                   sparse=TRUE,
                                   nrow=length(not_in_g))
      row.names(extra_rows) <- not_in_g
      colnames(extra_rows) <- colnames(exp)
      exp <- rbind(exp, extra_rows)
      exp <- exp
    }


    exprs_list[[i]] <- exp[gene_list, , drop=FALSE]
    fd_list[[i]] <- fd[gene_list, , drop=FALSE]
    pd_list[[i]] <- pd

  }
  all_fd <- array(NA,dim(fd_list[[1]]),dimnames(fd_list[[1]]))

  for (fd in fd_list) {
    for (j in colnames(fd)) {
      col_info <- all_fd[,j]
      col_info[is.na(col_info)] <- fd[is.na(col_info),j]
      col_info[col_info != fd[,j]] <- "conf"
      all_fd[,j] <- col_info
    }
  }

  confs <- sum(all_fd == "conf", na.rm=TRUE)

  if (confs > 0) {
   warning(paste0("When combining rowData, conflicting values were found - ",
                  "conflicts will be labelled 'conf' in the combined cds ",
                  "to prevent conflicts, either change conflicting values to ",
                  "match, or rename columns from different cds' to be unique."))
  }
  #all_fd <- do.call(cbind, fd_list)
  all_fd <- all_fd[,fdata_cols, drop=FALSE]

  all_pd <- do.call(rbind, pd_list)
  all_exp <- do.call(cbind, exprs_list)

  all_exp <- all_exp[row.names(all_fd), row.names(all_pd), drop=FALSE]

  new_cds <- new_cell_data_set(all_exp, cell_metadata = all_pd, gene_metadata = all_fd)

  if(keep_reduced_dims) {
    for(red_dim in names(reducedDims(cds_list[[1]]))) {
      reduced_dims_list <- list()
      for(j in 1:length(cds_list)) {
        reduced_dims_list[[j]] <- reducedDims(cds_list[[j]])[[red_dim]]
      }
      reducedDims(new_cds)[[red_dim]] <- do.call(rbind, reduced_dims_list, quote=FALSE)
      # The following should not happen; the accessor appears to ensure the
      # correct row order.
      if(!identical(rownames(reducedDims(new_cds)[[red_dim]]), rownames(all_pd))) {
        stop('Mis-ordered reduced matrix rows.')
      }
    }
  }
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
  cds@reduce_dim_aux <- SimpleList()
  cds@principal_graph_aux <- SimpleList()
  cds@principal_graph <- SimpleList()
  cds@clusters <- SimpleList()
  reducedDims(cds) <- SimpleList()
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
  if (is.null(metadata(cds)$citations) | citation_key == "Monocle") {
    metadata(cds)$citations <- data.frame(method = c("Monocle", "Monocle", "Monocle"),
                                          citations = c("Trapnell C. et. al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat. Biotechnol. 32, 381-386 (2014). https://doi.org/10.1038/nbt.2859",
                                                        "Qiu, X. et. al. Reversed graph embedding resolves complex single-cell trajectories. Nat. Methods 14, 979-982 (2017). https://doi.org/10.1038/nmeth.4402",
                                                        "Cao, J. et. al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496-502 (2019). https://doi.org/10.1038/s41586-019-0969-x"))
  }
  metadata(cds)$citations <- rbind(metadata(cds)$citations,
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
#' \dontrun{
#' get_citations(cds)
#' }
#' }
get_citations <- function(cds) {
  message(paste("Your analysis used methods from the following recent work.",
                "Please cite them wherever you are presenting your analyses."))
  if(is.null(metadata(cds)$citations)) {
    cds <- add_citation(cds, "Monocle")
  }
  metadata(cds)$citations
}


initialize_counts_metadata <- function(cds) {
  if(is.null(int_metadata(cds))) {
    int_metadata(cds) <- list()
  }
  int_metadata(cds)[['counts_metadata']] <- list()
  return(cds)
}


set_counts_identity <- function(cds, matrix_type, matrix_id) {
  int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']] <- matrix_type
  int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']] <- matrix_id
  return(cds)
}


get_counts_identity <- function(cds) {
  if(!is.null(int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']])) {
    matrix_id <- int_metadata(cds)[['counts_metadata']][['identity']][['matrix_id']]
  } else {
    matrix_id <- 'none'
  }
  if(!is.null(int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']])) {
    matrix_type <- int_metadata(cds)[['counts_metadata']][['identity']][['matrix_type']]
  } else {
    matrix_type <- 'matrix:counts'
  }
  return(list(matrix_id=matrix_id, matrix_type=matrix_type))
}


# Note: int_metadata(cds) requires a list
initialize_reduce_dim_metadata <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  if(is.null(int_metadata(cds))) {
    int_metadata(cds) <- list()
  }
  if(is.null(int_metadata(cds)[['reduce_dim_metadata']])) {
    int_metadata(cds)[['reduce_dim_metadata']] <- list()
  }
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]] <- list()
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']] <- list()
  return(cds)
}


set_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                           matrix_type, matrix_id,
                                           prev_matrix_type, prev_matrix_id,
                                           model_type, model_id) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_type']] <- matrix_type
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_id']] <- matrix_id
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_type']] <- prev_matrix_type
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_id']] <- prev_matrix_id
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_type']] <- model_type
  int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_id']] <- model_id
  return(cds)
}


get_reduce_dim_matrix_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  return(list(matrix_type=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_type']],
              matrix_id=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['matrix_id']],
              prev_matrix_type=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_type']],
              prev_matrix_id=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['prev_matrix_id']],
              model_type=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_type']],
              model_id=int_metadata(cds)[['reduce_dim_metadata']][[reduction_method]][['identity']][['model_id']]))
}


initialize_reduce_dim_aux_model <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  if(is.null(cds@reduce_dim_aux[[reduction_method]])) {
    cds@reduce_dim_aux[[reduction_method]] <- SimpleList()
  }
  cds@reduce_dim_aux[[reduction_method]][['model']] <- SimpleList()
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']] <- SimpleList()
  return(cds)
}


set_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'),
                                          model_type, model_id,
                                          prev_model_type, prev_model_id,
                                          model_path='none') {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_type']] <- model_type
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_id']] <- model_id
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_type']] <- prev_model_type
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_id']] <- prev_model_id
  cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- model_path
  return(cds)
}


get_reduce_dim_model_identity <- function(cds, reduction_method=c('PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP')) {
  assertthat::assert_that(
    tryCatch(expr = ifelse(match.arg(reduction_method) == "",TRUE, TRUE),
             error = function(e) FALSE),
    msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
  reduction_method <- match.arg(reduction_method)
  return(list(model_type=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_type']],
              model_id=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_id']],
              prev_model_type=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_type']],
              prev_model_id=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['prev_model_id']],
              model_path=cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']]))
}


set_model_identity_path <- function(cds, reduction_method, model_path='none') {
  if(!is.null(cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']])) {
    cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- model_path
  } else {
    cds@reduce_dim_aux[[reduction_method]][['model']][['identity']][['model_path']] <- 'none'
  }
  return(cds)
}


# Make a unique identifier string.
..global_unique_id_counter.. <- 1
get_unique_id <- function() {
  ..global_unique_id_counter.. <- ..global_unique_id_counter.. + 1
  rtime <- as.numeric(Sys.time())*100000 + ..global_unique_id_counter..
  id_hash <- openssl::md5(as.character(rtime))
  return(id_hash)
}


# Get current time stamp.
get_time_stamp <- function() {
  time_stamp <- format(Sys.time(), "%Y%m%d:%H%M%S" )
  return(time_stamp)
}


