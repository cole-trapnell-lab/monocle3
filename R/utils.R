# Check if class is a sparseMatrix from Matrix package
is_sparse_matrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix")
}

#' Function to calculate the size factor for the single-cell RNA-seq data
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
estimate_size_factors <- function(cds, locfunc = stats::median,
                                  round_exprs=TRUE,
                                  method=c("mean-geometric-mean-total",
                                           'mean-geometric-mean-log-total'))
{
  method <- match.arg(method)
  if (is_sparse_matrix(counts(cds))){
    size_factors(cds) <- estimate_sf_sparse(counts(cds),
                                            round_exprs=round_exprs,
                                            method=method)
  }else{
    size_factors(cds) <- estimate_sf_dense(counts(cds),
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
  lapply(parallel::splitIndices(nrow(x), ncl), function(i) x[i, , drop = FALSE])
}

#' @keywords internal
split_cols <- function (x, ncl) {
  lapply(parallel::splitIndices(ncol(x), ncl), function(i) x[, i, drop = FALSE])
}

#' @keywords internal
sparse_par_r_apply <- function (cl, x, FUN, convert_to_dense, ...)
{
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
sparse_par_c_apply <- function (cl = NULL, x, FUN, convert_to_dense, ...)
{
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
#' mcesApply computes the row-wise or column-wise results of FUN, just like
#' esApply. Variables in colData from cds are available in FUN.
#'
#' @param cds a cell_data_set object
#' @param MARGIN The margin to apply to, either 1 for rows (samples) or 2 for columns (features)
#' @param FUN Any function
#' @param required_packages A list of packages FUN will need. Failing to provide packages needed by FUN will generate errors in worker threads.
#' @param convert_to_dense Whether to force conversion a sparse matrix to a dense one before calling FUN
#' @param ... Additional parameters for FUN
#' @param cores The number of cores to use for evaluation
#'
#' @return The result of with(colData(cds) apply(counts(cds)), MARGIN, FUN, ...))
mc_es_apply <- function(cds, MARGIN, FUN, required_packages, cores=1, convert_to_dense=TRUE, reduction_method="UMAP", ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  coldata_df = as.data.frame(colData(cds))
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
  }, error = function(e) {} )

  coldata_df$pseudotime = pseudotime(cds) # Once we resolve issue #69.
  Biobase::multiassign(names(as.data.frame(coldata_df)), as.data.frame(coldata_df), envir=e1)
  environment(FUN) <- e1

  # Note: use outfile argument to makeCluster for debugging
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
      for (req in pkgs) {
        library(req, character.only=TRUE, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE)
      }
    }, required_packages)
  }

  if (MARGIN == 1){
    suppressWarnings(res <- sparse_par_r_apply(cl, counts(cds), FUN, convert_to_dense, ...))
  }else{
    suppressWarnings(res <- sparse_par_c_apply(cl, counts(cds), FUN, convert_to_dense, ...))
  }

  res
}

smart_es_apply <- function(cds, MARGIN, FUN, convert_to_dense, reduction_method="UMAP", ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  coldata_df = as.data.frame(colData(cds))
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
  }, error = function(e) {} )
  coldata_df$pseudotime = pseudotime(cds) # Once we resolve issue #69.
  Biobase::multiassign(names(as.data.frame(coldata_df)), as.data.frame(coldata_df), envir=e1)
  environment(FUN) <- e1

  if (is_sparse_matrix(counts(cds))){
    res <- sparse_apply(counts(cds), MARGIN, FUN, convert_to_dense, ...)
  } else {
    res <- pbapply::pbapply(counts(cds), MARGIN, FUN, ...)
  }

  if (MARGIN == 1)
  {
    names(res) <- row.names(cds)
  }else{
    names(res) <- colnames(cds)
  }

  res
}



#' Build a cell_data_set from the data stored in inst/extdata directory.
#' @export
load_a549 <- function(){
  small_a549_colData_df <- readRDS(system.file("extdata",
                                               "small_a549_dex_pdata.rda",
                                               package = "monocle3"))
  small_a549_rowData_df <- readRDS(system.file("extdata",
                                               "small_a549_dex_fdata.rda",
                                               package = "monocle3"))
  small_a549_exprs <- readRDS(system.file("extdata",
                                          "small_a549_dex_exprs.rda",
                                          package = "monocle3"))
  small_a549_exprs <- small_a549_exprs[,row.names(small_a549_colData_df)]

  cds <- new_cell_data_set(expression_data = small_a549_exprs,
                           cell_metadata = small_a549_colData_df,
                           gene_metadata = small_a549_rowData_df)
  cds
}

#' Build a cell_data_set from the data stored in inst/extdata directory.
#' @keywords internal
load_worm_embryo <- function(){
  expression_matrix = readRDS(url("http://jpacker-data.s3.amazonaws.com/for-cole/for-cole.ciliated.amphid.neurons.exprs.rds"))
  cell_metadata = readRDS(url("http://jpacker-data.s3.amazonaws.com/for-cole/for-cole.ciliated.amphid.neurons.pData.rds"))
  gene_annotation = readRDS(url("http://jpacker-data.s3.amazonaws.com/for-cole/for-cole.ciliated.amphid.neurons.fData.rds"))
  gene_annotation$use_for_ordering = NULL

  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  cds <- estimate_size_factors(cds)
  cds
}


# Test whether a matrix is one of our supported sparse matrices
is_sparse_matrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix")
}



#' Principal Components Analysis
#'
#' Efficient computation of a truncated principal components analysis of a given data matrix
#' using an implicitly restarted Lanczos method from the \code{\link{irlba}} package.
#'
#' @param x a numeric or complex matrix (or data frame) which provides
#'          the data for the principal components analysis.
#' @param retx a logical value indicating whether the rotated variables should be returned.
#' @param center a logical value indicating whether the variables should be
#'          shifted to be zero centered. Alternately, a centering vector of length
#'          equal the number of columns of \code{x} can be supplied.
#' @param scale. a logical value indicating whether the variables should be
#'          scaled to have unit variance before the analysis takes place.
#'          The default is \code{FALSE} for consistency with S, but scaling is often advisable.
#'          Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.
#'
#'          The value of \code{scale} determines how column scaling is performed
#'          (after centering).  If \code{scale} is a numeric vector with length
#'          equal to the number of columns of \code{x}, then each column of \code{x} is
#'          divided by the corresponding value from \code{scale}.  If \code{scale} is
#'          \code{TRUE} then scaling is done by dividing the (centered) columns of
#'          \code{x} by their standard deviations if \code{center=TRUE}, and the
#'          root mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#'          See \code{\link{scale}} for more details.
#' @param n integer number of principal component vectors to return, must be less than
#' \code{min(dim(x))}.
#' @param ... additional arguments passed to \code{\link{irlba}}.
#'
#' @return
#' A list with class "prcomp" containing the following components:
#' \itemize{
#'    \item{sdev} {the standard deviations of the principal components (i.e.,
#'          the square roots of the eigenvalues of the
#'          covariance/correlation matrix, though the calculation is
#'          actually done with the singular values of the data matrix).}
#'   \item{rotation} {the matrix of variable loadings (i.e., a matrix whose columns
#'          contain the eigenvectors).}
#'   \item {x} {if \code{retx} is \code{TRUE} the value of the rotated data (the centred
#'          (and scaled if requested) data multiplied by the \code{rotation}
#'         matrix) is returned.  Hence, \code{cov(x)} is the diagonal matrix
#'          \code{diag(sdev^2)}.}
#'   \item{center, scale} {the centering and scaling used, or \code{FALSE}.}
#' }
#'
#' @note
#' The signs of the columns of the rotation matrix are arbitrary, and
#' so may differ between different programs for PCA, and even between
#' different builds of R.
#'
#' NOTE DIFFERENCES WITH THE DEFAULT \code{\link{prcomp}} FUNCTION!
#' The \code{tol} truncation argument found in \code{prcomp} is not supported.
#' In place of the truncation tolerance in the original function, the
#' \code{prcomp_irlba}  function has the argument \code{n} explicitly giving the
#' number of principal components to return. A warning is generated if the
#' argument \code{tol} is used, which is interpreted differently between
#' the two functions.
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
#'
#' @seealso \code{\link{prcomp}}
sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, ...)
{
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba` function to
            control that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
  # Try to convert to a matrix...
  #if (!is.matrix(x)) x <- as.matrix(x)
  orig_x = x
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
      scale. = sqrt(DelayedMatrixStats::colVars(x))
      if (ans$scale) ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale. ^ 2)
    } else
    {
      if (ans$scale)
      {
        scale. = sqrt(DelayedMatrixStats::colSums2(x ^ 2) / (max(1, nrow(x) - 1L)))
        ans$totalvar = sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) / (nrow(x) - 1L)))
      } else
      {
        ans$totalvar = sum(DelayedMatrixStats::colSums2(x ^ 2) / (nrow(x) - 1L))
      }
    }
    if (ans$scale) args$scale <- scale.
  } else
  {
    args$scale <- scale.
    ans$totalvar = sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) / (nrow(x) - 1L)))
  }
  if (!missing(...)) args <- c(args, list(...))

  s <- do.call(irlba::irlba, args=args)
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
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
#' @description For each feature in a cell_data_set object, detect_genes counts
#' how many cells are expressed above a minimum threshold. In addition, for
#' each cell, detect_genes counts the number of genes above this threshold that
#' are detectable. Results are added as columns num_cells_expressed and
#' num_genes_expressed in the rowData and colData tables respectively.
#'
#' @param cds Input cell_data_set object.
#' @param min_expr Expression threshold
#' @return Updated cell_data_set object
#' @export
#' @examples
#' \dontrun{
#' cds <- detect_genes(cds, min_expr=0.1)
#' }
detect_genes <- function(cds, min_expr=0){
  assertthat::assert_that(is(cds, "cell_data_set"))
    assertthat::assert_that(is.numeric(min_expr))

  rowData(cds)$num_cells_expressed <- Matrix::rowSums(counts(cds) > min_expr)
  colData(cds)$num_genes_expressed <- Matrix::colSums(counts(cds) > min_expr)

  cds
}

#' Return a size-factor normalized and (optionally) log-transformed expression matrix
#' @export
normalized_counts <- function(cds, norm_method=c("log", "binary", "size_only"), pseudocount=1){
  norm_method = match.arg(norm_method)
  norm_mat = counts(cds)
  if (norm_method == "binary"){
    norm_mat = norm_mat > 0
    if (is_sparse_matrix(norm_mat)){
      norm_mat = as(norm_mat, "dgCMatrix")
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
        norm_mat@x = log10(norm_mat + pseudocount)
      }
    }
  }
  return(norm_mat)
}
