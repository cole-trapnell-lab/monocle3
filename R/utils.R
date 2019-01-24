# Check if class is a sparseMatrix from Matrix package
is_sparse_matrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix")
}

#' Function to calculate the size factor for the single-cell RNA-seq data
#'
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#'
estimate_sf_matrix <- function(counts, locfunc = stats::median, round_exprs=TRUE,  method="mean-geometric-mean-total")
{
  if (is_sparse_matrix(counts)){
    estimate_sf_sparse(counts, locfunc = locfunc, round_exprs=round_exprs, method=method)
  }else{
    estimate_sf_dense(counts, locfunc = locfunc, round_exprs=round_exprs,  method=method)
  }

}


# Estimate size factors for each column, given a sparseMatrix from the Matrix
# package
estimate_sf_sparse <- function(counts,
                               locfunc = median,
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
estimate_sf_dense <- function(counts, locfunc = median, round_exprs=TRUE, method="mean-geometric-mean-total"){

  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  if (method == "weighted-median"){
    log_medians <- apply(CM, 1, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowMeans(log(CM))

    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    row_median <- apply(CM, 1, median)
    sfs <- apply(Matrix::t(Matrix::t(CM) - row_median), 2, median)
  }else if(method == 'mode'){
    sfs <- estimate_t(CM)
  }else if(method == 'geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
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

#' @importFrom parallel splitIndices
#' @keywords internal
split_rows <- function (x, ncl) {
  lapply(splitIndices(nrow(x), ncl), function(i) x[i, , drop = FALSE])
}

#' @importFrom parallel splitIndices
#' @keywords internal
split_cols <- function (x, ncl) {
  lapply(splitIndices(ncol(x), ncl), function(i) x[, i, drop = FALSE])
}

#' @importFrom BiocGenerics clusterApply
#' @keywords internal
sparse_par_r_apply <- function (cl, x, FUN, convert_to_dense, ...)
{
  par_res <- do.call(c, clusterApply(cl = cl, x = split_rows(x, length(cl)),
                                     fun = sparse_apply, MARGIN = 1L, FUN = FUN, convert_to_dense=convert_to_dense, ...), quote = TRUE)
  names(par_res) <- row.names(x)
  par_res
}

#' @importFrom BiocGenerics clusterApply
#' @keywords internal
sparse_par_c_apply <- function (cl = NULL, x, FUN, convert_to_dense, ...)
{
  par_res <- do.call(c, clusterApply(cl = cl, x = split_cols(x, length(cl)),
                                     fun = sparse_apply, MARGIN = 2L, FUN = FUN, convert_to_dense=convert_to_dense, ...), quote = TRUE)
  names(par_res) <- colnames(x)
  par_res
}


#' Multicore apply-like function for CellDataSet
#'
#' mcesApply computes the row-wise or column-wise results of FUN, just like esApply.
#' Variables in pData from X are available in FUN.
#'
#' @param X a CellDataSet object
#' @param MARGIN The margin to apply to, either 1 for rows (samples) or 2 for columns (features)
#' @param FUN Any function
#' @param required_packages A list of packages FUN will need. Failing to provide packages needed by FUN will generate errors in worker threads.
#' @param convert_to_dense Whether to force conversion a sparse matrix to a dense one before calling FUN
#' @param ... Additional parameters for FUN
#' @param cores The number of cores to use for evaluation
#'
#' @return The result of with(pData(X) apply(exprs(X)), MARGIN, FUN, ...))
#' @importFrom parallel makeCluster stopCluster
#' @importFrom BiocGenerics clusterCall parRapply parCapply
#' @importFrom Biobase pData exprs multiassign
#' @export
mc_es_apply <- function(X, MARGIN, FUN, required_packages, cores=1, convert_to_dense=TRUE, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1

  # Note: use outfile argument to makeCluster for debugging
  platform <- Sys.info()[['sysname']]
  if (platform == "Windows")
    cl <- makeCluster(cores)
  if (platform %in% c("Linux", "Darwin"))
    cl <- makeCluster(cores, type="FORK")

  cleanup <- function(){
    stopCluster(cl)
  }
  on.exit(cleanup)

  if (is.null(required_packages) == FALSE){
    clusterCall(cl, function(pkgs) {
      for (req in pkgs) {
        library(req, character.only=TRUE)
      }
    }, required_packages)
  }
  #clusterExport(cl, ls(e1), e1)
  #force(exprs(X))
  if (MARGIN == 1){
    suppressWarnings(res <- sparse_par_r_apply(cl, exprs(X), FUN, convert_to_dense, ...))
  }else{
    suppressWarnings(res <- sparse_par_c_apply(cl, exprs(X), FUN, convert_to_dense, ...))
  }

  res
}

#' @importFrom Biobase multiassign
#' @importFrom pbapply pbapply
smart_es_apply <- function(X, MARGIN, FUN, convert_to_dense, ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1

  if (isSparseMatrix(exprs(X))){
    res <- sparse_apply(exprs(X), MARGIN, FUN, convert_to_dense, ...)
  }else{
    res <- pbapply(exprs(X), MARGIN, FUN, ...)
  }

  if (MARGIN == 1)
  {
    names(res) <- row.names(X)
  }else{
    names(res) <- colnames(X)
  }

  res
}



#' Build a CellDataSet from the data stored in inst/extdata directory.
#' @importFrom Biobase pData pData<- exprs fData
#' @export
load_a549 <- function(){

  baseLoc <- system.file(package="monocle3")
  #baseLoc <- './inst'
  extPath <- file.path(baseLoc, "extdata")
  small_a549_pdata_df = readRDS(file.path(extPath, "small_a549_dex_pdata.rda"))
  small_a549_fdata_df = readRDS(file.path(extPath, "small_a549_dex_fdata.rda"))
  small_a549_exprs = readRDS(file.path(extPath, "small_a549_dex_exprs.rda"))

  small_a549_exprs <- small_a549_exprs[,row.names(small_a549_pdata_df)]

  pd <- new("AnnotatedDataFrame", data = small_a549_pdata_df)
  fd <- new("AnnotatedDataFrame", data = small_a549_fdata_df)

  # Now, make a new CellDataSet using the RNA counts
  cds <- new_cell_data_set(small_a549_exprs,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily="negbinomial")
  pData(cds)$Size_Factor = small_a549_pdata_df$Size_Factor

  cds
}

# Test whether a matrix is one of our supported sparse matrices
#' @import Matrix
isSparseMatrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix")
}




