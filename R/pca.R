#
# pca_control elements
#   matrix_class  'CsparseMatrix' or 'BPCells'  default: 'CsparseMatrix'
#   matrix_mode   'mem' or 'dir'  default: 'mem'
#   matrix_type   'float' or 'double'
#   matrix_path   default: NULL
#   matrix_buffer_size: <integer> default: 8192L

# Usage
#   matrix_class: default: 'CsparseMatrix'
#     matrix_mode: default: 'mem'
#   matrix_class: 'BPCells'
#     matrix_mode: 'mem'  default: 'dir'
#       matrix_type: 'float', 'double'
#     matrix_mode: 'dir'
#       matrix_type: 'float', 'double' default: 'double'
#       matrix_path: <path to directory or file> default: 'NULL' -> temporary directory in pwd
#       matrix_buffer_size: <integer> default: 8192L

select_pca_parameter_value <- function(parameter, pca_control, pca_control_default, default_value) {
  if(!is.null(pca_control[[parameter]])) {
    return(pca_control[[parameter]])
  }
  else
  if(!is.null(pca_control_default[[parameter]])) {
    return(pca_control_default[[parameter]])
  }
  return(default_value)
}


#' Verify and set the pca_control parameter list.
#'
#' @description Verifies and sets the list of parameter values
#'   that is used when storing matrices that are used for PCA
#'   and LSI. To see the default values,
#'   call "set_pca_control(pca_control=list(show_values=TRUE))".
#'   "show_values=TRUE" can be used in functions with the
#'   pca_control list parameter in which case the function will
#'   show the pca_control values to be used and then stop.
#' @param pca_control Input control list.
#' @return pca_control Output control list.
#'
#' @section pca_control parameters:
#' \describe{
#'   \item{matrix_class}{Specifies the matrix class to use for
#'      matrix storage. The acceptable values are "CsparseMatrix"
#'      and "BPCells".}
#'   \item{matrix_type}{Specifies whether to store the matrix
#'      single precision "floats" (matrix_type="float") or
#'      double precision "doubles" (matrix_type="double").
#'      "matrix_type" is used only for BPCells class matrices.}
#'   \item{matrix_mode}{Specifies whether to store the BPCells
#'      class matrix in memory (matrix_mode="mem") or on
#'      disk (matrix_mode="dir"). "matrix_mode" is used only
#'      for BPCells class matrices.}
#'   \item{matrix_path}{Specifies the disk directory where the
#'      BPCells on-disk matrix data are stored. The default is a
#'      directory, with a randomized name, in the directory
#'      where R is running. "matrix_path" is used only for
#'      BPCells class matrices with matrix_mode="dir".}
#'   \item{matrix_compress}{Specifies whether to use bit-packing
#'      compression to store BPCells matrix values. This
#'      reduces somewhat storage requirements but increases
#'      run time. "matrix_compress" is used only for BPCells
#'      class matrices.}
#'   \item{matrix_buffer_size}{Specifies how many items of
#'      data to buffer in memory before flushing to disk. This
#'      is used for matrix_class="BPCells" with matrix_mode="dir".}
#' @export
set_pca_control <- function(pca_control=list()) {

  assertthat::assert_that(methods::is(pca_control, "list"))

  allowed_control_parameters <- c('matrix_class',
                                  'matrix_mode',
                                  'matrix_type',
                                  'matrix_path',
                                  'matrix_buffer_size',
                                  'show_values')

  allowed_matrix_class <- c('CsparseMatrix', 'BPCells')


  assertthat::assert_that(methods::is(pca_control, "list"))

  if(is.null(pca_control[['matrix_class']])) {
    pca_control[['matrix_class']] <- 'CsparseMatrix'
  }
  else
  if(!(pca_control[['matrix_class']] %in% allowed_matrix_class)) {
    stop('  matrix_class value must be "CsparseMatrix" or "BPCells"')
  }

  if(pca_control[['matrix_class']] == 'CsparseMatrix') {
    pca_control_default <- get_global_variable('pca_control_csparsematrix')
  }
  else
  if(pca_control[['matrix_class']] == 'BPCells') {
    pca_control_default <- get_global_variable('pca_control_bpcells')
  }

  assertthat::assert_that(all(names(pca_control) %in% allowed_control_parameters),
                          msg = "set_pca_control: unknown variable in pca_control")

  assertthat::assert_that(all(names(pca_control_default) %in% allowed_control_parameters),
                          msg = "set_pca_control: unknown variable in pca_control_default")

  #
  # Last resort fall-back parameter values.
  #
  default_matrix_class <- 'CsparseMatrix'
  default_matrix_mode <- 'mem'
  default_matrix_type <- 'double'
  default_matrix_path <- NULL
  default_matrix_buffer_size <- 8192L

  error_string = list()

  pca_control_out = list()
  pca_control_out[['matrix_class']] <- select_pca_parameter_value('matrix_class', pca_control, pca_control_default, default_matrix_class)

  if(pca_control_out[['matrix_class']] == 'CsparseMatrix') {
    pca_control_out[['matrix_mode']] <- select_pca_parameter_value('matrix_mode', pca_control, pca_control_default, default_matrix_mode)
    if(!(pca_control_out[['matrix_mode']] %in% c('mem'))) {
      error_string <- list(error_string, paste0('  ',
                           pca_control_out[['matrix_mode']],
                           ' is not valid for matrix_class CsparseMatrix.'))
      stop(error_string)
    }
  }
  else
  if(pca_control[['matrix_class']] == 'BPCells') {
    pca_control_out[['matrix_mode']] <- select_pca_parameter_value('matrix_mode', pca_control, pca_control_default, default_matrix_mode)
    if(!(pca_control_out[['matrix_mode']] %in% c('mem', 'dir'))) {
      error_string <- list(error_string, paste0('  ',
                           pca_control_out[['matrix_mode']],
                           ' is not a valid matrix_mode for matrix_class "BPCells".'))
      stop(error_string)
    }
    if(pca_control_out[['matrix_mode']] == 'mem') {
      pca_control_out[['matrix_type']] <- select_pca_parameter_value('matrix_type', pca_control, pca_control_default, default_matrix_type)
      if(!(pca_control_out[['matrix_type']] %in% c('float', 'double'))) {
        error_string <- list(error_string, paste0('  ',
                             pca_control_out[['matrix_type']],
                             ' is not a valid matrix_type for matrix_class "BPCells".'))
      }
      if(length(error_string) > 0) {
        stop(error_string)
      }
    }
    else
    if(pca_control_out[['matrix_mode']] == 'dir') {
      pca_control_out[['matrix_type']] <- select_pca_parameter_value('matrix_type', pca_control, pca_control_default, default_matrix_type)
      if(!(pca_control_out[['matrix_type']] %in% c('float', 'double'))) {
        error_string <- list(error_string, paste0('  ',
                             pca_control_out[['matrix_type']],
                             ' is not a valid matrix_type for matrix_class "BPCells".'))
      }
      pca_control_out[['matrix_path']] <- select_pca_parameter_value('matrix_path', pca_control, pca_control_default, default_matrix_path)
      if(!is.null(pca_control_out[['matrix_path']]) &&
         !is.character(pca_control_out[['matrix_path']])) {
        error_string <- list(error_string, paste0('  ',
                             pca_control_out[['matrix_path']],
                             ' is not a valid matrix_path.'))

      }
      pca_control_out[['matrix_buffer_size']] <- select_pca_parameter_value('matrix_buffer_size', pca_control, pca_control_default, default_matrix_buffer_size)
      if(!is.integer(pca_control_out[['matrix_buffer_size']])) {
        error_string <- list(error_string, paste0('  ',
                             pca_control_out[['matrix_buffer_size']],
                             ' is not an integer.'))
      }
      if(length(error_string) > 0) {
        stop(error_string)
      }
    }
  }

  #
  # Set BPCells out-of-core directory name, if necessary.
  #
  if(pca_control_out[['matrix_class']] == 'BPCells' && is.null(pca_control_out[['matrix_path']])) {
    if(pca_control_out[['matrix_mode']] == 'dir') {
      pca_control_out[['matrix_path']] <- tempfile(pattern='monocle.bpcells.', tmpdir='.', fileext='.tmp')[[1]]
    }
  }

  #
  # Display pca_control list if pca_control[['show_values']] <- TRUE.
  #
  if(!is.null(pca_control[['show_values']]) && pca_control[['show_values']] == TRUE)
  {
    report_pca_control('  pca_control: ', pca_control=pca_control_out)
    stop_no_noise()
  }

  return(pca_control_out)
}


# Report pca_control list values.
report_pca_control <- function(label=NULL, pca_control) {
  indent <- ''
  if(!is.null(label)) {
    indent <- '  '
  }

  message(ifelse(!is.null(label), label, ''))

  message(indent, '  matrix_class: ', ifelse(!is.null(pca_control[['matrix_class']]), pca_control[['matrix_class']], as.character(NA)))

  if(pca_control[['matrix_class']] == 'CsparseMatrix') {
    message(indent, '  matrix_mode: ', ifelse(!is.null(pca_control[['matrix_mode']]), pca_control[['matrix_mode']], as.character(NA)))
  }
  else
  if(pca_control[['matrix_class']] == 'BPCells') {
    if(pca_control[['matrix_mode']] == 'mem') {
      message(indent, '  matrix_type: ', ifelse(!is.null(pca_control[['matrix_type']]), pca_control[['matrix_type']], as.character(NA)))
    }
    else
    if(pca_control[['matrix_mode']] == 'dir') {
      message(indent, '  matrix_type: ', ifelse(!is.null(pca_control[['matrix_type']]), pca_control[['matrix_type']], as.character(NA)))
      message(indent, '  matrix_path: ', ifelse(!is.null(pca_control[['matrix_path']]), pca_control[['matrix_path']], as.character(NA)))
      message(indent, '  matrix_buffer_size: ', ifelse(!is.null(pca_control[['matrix_buffer_size']]), pca_control[['matrix_buffer_size']], as.character(NA)))
    }
    else {
      stop('report_nn_control: unsupported pca class/mode/...\'', pca_control[['method']], '\'')
    }
  }
  else {
    stop('report_nn_control: unsupported pca class/mode/...\'', pca_control[['method']], '\'')
  }
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
#' \dontrun{
#'   set.seed(1)
#'   x  <- matrix(rnorm(200), nrow=20)
#'   p1 <- irlba::prcomp_irlba(x, n=3)
#'   summary(p1)
#'
#'   # Compare with
#'   p2 <- prcomp(x, tol=0.7)
#'   summary(p2)}
#'
#' @seealso \code{\link{prcomp}}
sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE,
                                scale. = FALSE, verbose = FALSE, ...)
{
  if(verbose) {
    message('pca: sparse_prcomp_irlba: matrix class: ', class(x))
  }

  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`
            function to control that algorithm's convergence tolerance. See
            `?prcomp_irlba` for help.")
  orig_x <- x
  if (!methods::is(x, "DelayedMatrix")) {
    x = DelayedArray::DelayedArray(x)
  }

  args <- list(A=orig_x, nv=n)
  if (is.logical(center))
  {
    if (center) args$center <- DelayedMatrixStats::colMeans2(x)
  }
  else args$center <- center
  if (is.logical(scale.))
  {
    if (is.numeric(args$center))
    {
      scale. <- sqrt(DelayedMatrixStats::colVars(x))
      if (ans$scale) ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale. ^ 2)
    }
    else {
      if (ans$scale)
      {
        scale. <-
          sqrt(DelayedMatrixStats::colSums2(x ^ 2) / (max(1, nrow(x) - 1L)))
        ans$totalvar <-
          sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) /
                     (nrow(x) - 1L)))
      }
      else {
        ans$totalvar <-
          sum(DelayedMatrixStats::colSums2(x ^ 2) / (nrow(x) - 1L))
      }
    }
    if (ans$scale) args$scale <- scale.
  }
  else {
    args$scale <- scale.
    ans$totalvar <-
      sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.) ^ 2) /
                 (nrow(x) - 1L)))
  }
  if (!missing(...)) args <- c(args, list(...))


  if(verbose) {
    message('start irlba: ', Sys.time())
  }
  s <- do.call(irlba::irlba, args=args)
  if(verbose) {
    message('end irlba: ', Sys.time())
  }

  if(verbose) {
    message('singular values (head)')
    message(paste(head(s$d), collapse=' '))
  }

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

  if(verbose) {
    message("umat: ", paste(dim(s$u), collapse=" "))
    message("vtmat: ", paste(dim(s$v), collapse=" "))
  }

  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}


#' Principal Components Analysis on BPCells IterableMatrix
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
#' \dontrun{
#'   set.seed(1)
#'   x  <- matrix(rnorm(200), nrow=20)
#'   p1 <- irlba::prcomp_irlba(x, n=3)
#'   summary(p1)
#'
#'   # Compare with
#'   p2 <- prcomp(x, tol=0.7)
#'   summary(p2)}
#'
#' @seealso \code{\link{prcomp}}
bpcells_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE,
                                 scale. = FALSE, pca_control=list(),
                                 verbose = FALSE, ...)
{
  if(verbose) {
    message('pca: bpcells_prcomp_irlba: "x" matrix class: ', class(bpcells_find_base_matrix(x)))
    message('pca: bpcells_prcomp_irlba: "x" matrix info')
    message(bpcells_base_matrix_info(x))
  }


  assertthat::assert_that(pca_control[['matrix_class']] == 'BPCells')

  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`
            function to control that algorithm's convergence tolerance. See
            `?prcomp_irlba` for help.")

  if(pca_control[['matrix_mode']] == 'mem') {
    x_commit <- BPCells::write_matrix_memory(x, compress=FALSE)
  }
  else
  if(pca_control[['matrix_mode']] == 'dir') {
    matrix_path <- paste0(pca_control[['matrix_path']], '.2')
    buffer_size <- pca_control[['matrix_buffer_size']]
    x_commit <- BPCells::write_matrix_dir(mat=x, dir=matrix_path, compress=FALSE, buffer_size=buffer_size, overwrite=FALSE)
    push_matrix_path(matrix_path, 'bpcells_dir')
  }

  stats <- BPCells::matrix_stats(matrix = x_commit, row_stats = 'none', col_stats = 'variance')
  center <- stats[['col_stats']]['mean',]
  scale <- sqrt(stats[['col_stats']]['variance',])

  # BPCells:::linear_operator() is meant to reduce irlba run time.
  args <- list(A = BPCells:::linear_operator(x_commit), nv = n, center = center, scale = scale)

  if (!missing(...)) args <- c(args, list(...))

  if(verbose) {
    message('start time: ', Sys.time())
  }
  s <- do.call(irlba::irlba, args=args)
  if(verbose) {
    message('end time: ', Sys.time())
  }

  if(pca_control[['matrix_mode']] == 'dir') {
    unlink(matrix_path, recursive=TRUE)
    rm(x_commit)
  }

  if(verbose) {
    message('singular values (head)')
    message(paste(head(s$d), collapse=' '))
  }

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

  if(verbose) {
    message("umat: ", paste(dim(s$u), collapse=" "))
    message("vtmat: ", paste(dim(s$v), collapse=" "))
  }

  class(ans) <- c("irlba_prcomp", "prcomp")

  ans
}


make_pca_matrix <- function(FM, pca_control=list()) {
  matrix_class <- pca_control[['matrix_class']]
  matrix_path <- NULL
  if(matrix_class == 'CsparseMatrix') {
    if(!is_sparse_matrix(FM)) {
      FM <- as(FM, 'dgCMatrix')
    }
  }
  else
  if(matrix_class == 'BPCells') {
    # Force running queued operations by 'writing' a matrix because
    # we want to not add subsequent operations to the mat queue.
    matrix_mode <- pca_control[['matrix_mode']]
    matrix_type <- pca_control[['matrix_type']]
    if(matrix_mode == 'mem') {
      FM <- BPCells::write_matrix_memory(mat=BPCells::convert_matrix_type(FM, type=matrix_type), compress=FALSE)
    } 
    else
    if(matrix_mode == 'dir') {
      matrix_path <- pca_control[['matrix_path']]
      assertthat::assert_that(!is.null(matrix_path))
      matrix_buffer_size <- pca_control[['matrix_buffer_size']]
      if(!is(FM, 'IterableMatrix') && !is(FM, 'CsparseMatrix')) {
        FM <- as(FM, 'CsparseMatrix')
      }
      FM <- BPCells::write_matrix_dir(mat=BPCells::convert_matrix_type(FM, type=matrix_type), dir=matrix_path, compress=FALSE, buffer_size=matrix_buffer_size, overwrite=FALSE)
      push_matrix_path(matrix_path, 'bpcells_dir')
    }
  }

  return(list(mat=FM, matrix_path=matrix_path))
}


