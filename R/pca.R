svd_rebuild_matrix <- function(u, s, v, filename) {
  # m = u * s * vt
  mat <- u %*% diag(s, nrow=ncol(u), ncol=nrow(t(v))) %*% t(v)
  write.table(as.numeric(mat), file=filename)
}


#
# Notes on setting the PCA matrix control using a relatively
# complex approach: do not require input_matrix_control[['matrix_class']] but allow edits
# Possible conditions:
#   1  input_matrix_control has no values set (OK!)
#        o  want input_matrix_info with missing values from default for matrix_class/pca (probably no missing values)
#
#   2  input_matrix_control has matrix_class and other values set and input_matrix_control[['matrix_class']] == input_matrix_info[['matrix_class']]  (same as 4)
#        o  want input_matrix_info with edits by input_matrix_control
#             example: input_matrix_control[['matrix_path'] replaces input_matrix_info[['matrix_path']
#
#   3  input_matrix_control has matrix_class and other values set and input_matrix_control[['matrix_class']] != input_matrix_info[['matrix_class']] (OK? is there any reason to use input_matrix_info values?)
#        o  want input_matrix_control with missing values from default for matrix_class/pca
#
#   4  input_matrix_control does not have matrix_class set but other values are set  (same as 2)
#        o  want input_matrix_info with edits from input_matrix_control
#             example: input_matrix_control[['matrix_path'] replaces input_matrix_info[['matrix_path']
#
#
# Is any matrix_control command line parameter value set? (Check the following summary carefully
# before proceeding with it.)
#   o  yes: is input_matrix_control[['matrix_class']] set?
#     o  is input_matrix_control[['matrix_class']] the same as the input matrix matrix_class
#       o  yes: use matrix_control:         input_matrix_control  (2)  (*** wrong ***)
#                   matrix_control_default: input matrix info (edit matrix_path as required)
#       o  no:  use matrix_control:         input_matrix_control  (3)
#                   matrix_control_default: default matrix_control for input_matrix_control[['matrix_class']] (global environment value)
#     o  no: use matrix_control: input_matrix_control  (4)  (*** wrong ***)
#                matrix_control_default: default matrix_control for input_matrix_control[['matrix_class']] (global environment value)
#   o  no: use matrix_control: input matrix info (edit matrix_path as required)  (1)
#              matrix_control_default: default matrix_control for input_matrix_control[['matrix_class']] (global environment value)
#
#  Notes:
#    o  use the simple approach because it's substantially more provable, maintainable, and documentable.
#       Also, users cannot set matrix_control in preprocess_cds() and preprocess_transform().
#    o  report the resulting matrix_control after calling set_matrix_control_pca(), when verbose=TRUE
#


#
# Use the simple approach (relatively simple): require either input_matrix_control
# is zero length list or input_matrix_control[['matrix_class']] is set with all
# desired non-default values.
#
#  Notes:
#    o  if input_matrix_control is not set, use input_matrix_info with modifications required for pca and matrix_path
#    o  if input_matrix_control is set, use matrix_control=input_matrix_control and matrix_control_default for matrix_class/pca.
#       In this case input_matrix_control[['matrix_class']] must be set; otherwise, it's invalid and an error results. Missing
#       values in input_matrix_control are taken from the default.
#    o  watch for issues with matrix_path
#
set_matrix_control_pca <- function(mat, matrix_control=list(), verbose=FALSE) {

  assertthat::assert_that(monocle3:::is_matrix(mat),
                          msg=paste0('set_matrix_control_pca: input matrix parameter object is not a matrix'))
  assertthat::assert_that(is.list(matrix_control),
                          msg=paste0('set_matrix_control_pca: input matrix_control parameter object is not a list'))

  # Get matrix_control values of input matrix.
  matrix_info <- tryCatch(monocle3:::get_matrix_info(mat),
                          error = function(c) { stop(paste0(trimws(c),
                                                           '\n* error in set_matrix_control_pca')) })
  if(length(matrix_control) > 0) {
    # Use matrix_control=matrix_control
    assertthat::assert_that(!is.null(matrix_control[['matrix_class']]),
                            msg=paste0('set_matrix_control_pca: matrix_control[[\'matrix_class\']] is not set'))
#    if(matrix_control[['matrix_class']] == 'BPCells')
#      matrix_control_default <- get_global_variable('matrix_control_bpcells_pca')
#    else
#      matrix_control_default <- get_global_variable('matrix_control_csparsematrix_pca')

    matrix_control_default <- tryCatch(set_matrix_control_default(matrix_control, 'pca'),
                                       error = function(c) { stop(paste0(trimws(c), '\n* error in set_matrix_control_pca')) })

    matrix_control_res <- set_matrix_control(matrix_control=matrix_control, matrix_control_default=matrix_control_default, control_type='pca')
  }
  else {
    # Use matrix_control=matrix_info
    if(matrix_info[['matrix_class']] == 'BPCells') {
      # Trim off name of input matrix directory.
      tmp_matrix_info <- matrix_info
      tmp_matrix_info[['matrix_path']] <- dirname(tmp_matrix_info[['matrix_path']])
    }
    else {
      tmp_matrix_info <- matrix_info
    }

    matrix_control_default <- tryCatch(set_matrix_control_default(matrix_info, 'pca'),
                                error = function(c) { stop(paste0('\n* error in set_matrix_control_pca')) })

    matrix_control_res <- set_matrix_control(matrix_control=tmp_matrix_info, matrix_control_default=matrix_control_default, control_type='pca')
  }

  if(verbose) {
    message('set_matrix_control_pca: matrix_control_res:')
    show_matrix_control(matrix_control_res)
  }

  return(matrix_control_res)
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
    message(paste0('matrix_info:\n', show_matrix_info(matrix_info=get_matrix_info(mat=x), indent='  ')), appendLF=FALSE)
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

  # Diagnostic test.
  #svd_rebuild_matrix(s$u, s$d, s$v, 'dgcmatrix.vec')

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
                                 scale. = FALSE, verbose = FALSE, ...)
{
  if(verbose) {
    message('pca: bpcells_prcomp_irlba: matrix class: ', class(x))
    message(paste0('matrix_info:\n', show_matrix_info(matrix_info=get_matrix_info(mat=x), indent='  ')), appendLF=FALSE)
  }

  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`
            function to control that algorithm's convergence tolerance. See
            `?prcomp_irlba` for help.")

  # Use the same matrix_control for the 'x_commit' matrix as used for the
  # input matrix 'x'.
  matrix_control_res <- set_matrix_control_pca(mat=x, verbose=verbose)
  x_commit <- set_matrix_class(mat=x, matrix_control=matrix_control_res)

  stats <- BPCells::matrix_stats(matrix = x_commit, row_stats = 'none', col_stats = 'variance')
  center <- stats[['col_stats']]['mean',]
  scale <- sqrt(stats[['col_stats']]['variance',])

  if(verbose) {
    message('pca: bpcells_prcomp_irlba: x_commit:')
    message(show_matrix_info(matrix_info=get_matrix_info(mat=x_commit), indent='  '), appendLF=FALSE)
  }

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

  rm_bpcells_dir(mat=x_commit)

  # Diagnostic test.
  #message('bpcells svd')
  #svd_rebuild_matrix(s$u, s$d, s$v, 'bpcmatrix.vec')

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

