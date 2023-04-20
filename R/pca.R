set_pca_matrix_control <- function(mat, matrix_control=list()) {
  matrix_info <- get_matrix_info(mat)
  if(!is.null(matrix_control[['matrix_class']])) {
    if(matrix_control[['matrix_class']] == 'dgCMatrix') {
       matrix_control_default <- get_global_variable('pca_control_csparsematrix')
    }
    else
    if(matrix_control[['matrix_class']] == 'BPCells') {
      matrix_control_default <- get_global_variable('pca_control_bpcells')
    }
    else {
      stop('unrecognized matrix_control[[\'matrix_class\']] value')
    }
    matrix_control_res <- set_matrix_control(matrix_control=matrix_control,
                                             matrix_control_default=matrix_control_default,
                                             control_type='pca')
  }
  else
  if(matrix_info[['matrix_class']] == 'r_dense_matrix' ||
     matrix_info[['matrix_class']] == 'dgCMatrix' ||
     matrix_info[['matrix_class']] == 'dgTMatrix') {
    matrix_control_res = list()
    matrix_control_res[['matrix_class']] <- 'dgCMatrix'
  }
  else
  if(matrix_info[['matrix_class']] == 'BPCells') {
    matrix_control_res <- matrix_info
    if(matrix_control_res[['matrix_type']] == 'uint32_t') {
      matrix_control_res[['matrix_type']] <- 'double'
    }
    if(matrix_control_res[['matrix_compress']] == TRUE) {
      matrix_control_res[['matrix_compress']] <- FALSE
    }
    if(matrix_control_res[['matrix_mode']] == 'dir') {
      matrix_control_res[['matrix_path']] <- dirname(matrix_control_res[['matrix_path']])
    }
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
                                 scale. = FALSE, matrix_control=list(),
                                 verbose = FALSE, ...)
{
  if(verbose) {
    message(show_matrix_info(get_matrix_info(mat=x), indent='  '), appendLF=FALSE)
  }

  message('bpcells_prcomp_irlba: matrix_control matrix_class: ', matrix_control[['matrix_class']])

  if(!is.null(matrix_control[['class']]) && matrix_control[['matrix_class']] != 'BPCells') {
    stop('bpcells_prcomp_irlba: matrix_control[[\'matrix_class\']] must be \'BPCells\'')
  }

  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  if ("tol" %in% a)
    warning("The `tol` truncation argument from `prcomp` is not supported by
            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`
            function to control that algorithm's convergence tolerance. See
            `?prcomp_irlba` for help.")

  message('\n==== bpcells_prcomp_irlba ====')
  matrix_control <- get_matrix_info(x)
  if(matrix_control[['matrix_class']] == 'BPCells' &&
     matrix_control[['matrix_mode']] == 'dir') {
    matrix_control[['matrix_path']] <- dirname(matrix_control[['matrix_path']])
  }
  x_commit <- set_matrix_class(mat=x, matrix_control=matrix_control)

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

  rm_bpcells_dirs(x_commit)
  rm(x_commit)

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

