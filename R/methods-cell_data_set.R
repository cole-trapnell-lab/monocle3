#' Methods for the cell_data_set class
#' @name cell_data_set-methods
#' @docType methods
#' @rdname cell_data_set-methods
#' @param object The cell_data_set object
setValidity( "cell_data_set", function( object ) {
  TRUE
} )

#' @rdname cell_data_set-methods
#' @aliases cell_data_set,ANY,ANY,ANY-method
#' @param x the cell_data_set object
#' @param i index (or name) to extract or replace
#' @param j index (or name) to extract or replace
#' @param ... extra arguments passed to method
#' @param drop If TRUE the result is coerced to the lowest possible dimension
#'   (see the examples). This only works for extracting elements, not for the
#'   replacement.
#' @docType methods
#' @rdname extract-methods
setMethod("[", "cell_data_set", function(x, i, j, ..., drop = FALSE) {
  if (missing(drop))
    drop <- FALSE
  if (missing(i) && missing(j)) {
    if (!missing(...))
      stop("specify genes or samples to subset; use '",
           substitute(x), "$", names(list(...))[[1]],
           "' to access phenoData variables")
    return(x)
  }
  if (!isVersioned(x) || !isCurrent(x)["eSet"])
    x <- updateObject(x)
  if (!missing(j)) {
    phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
    protocolData(x) <- protocolData(x)[j,, ..., drop = drop]
  }
  if (!missing(i))
    featureData(x) <- featureData(x)[i,,..., drop=drop]
  ## assayData; implemented here to avoid function call
  orig <- assayData(x)
  storage.mode <-  Biobase:::assayDataStorageMode(orig)
  assayData(x) <-
    switch(storage.mode,
           environment =,
           lockedEnvironment = {
             aData <- new.env(parent=emptyenv())
             if (missing(i))  {                   # j must be present
               for(nm in ls(orig)) {
                 aData[[nm]] <- orig[[nm]][, j, ..., drop = drop]
               }
             } else if (missing(j)) { # j may or may not be present
               for(nm in ls(orig)) {
                 aData[[nm]] <- orig[[nm]][i,, ..., drop = drop]
               }
             }  else {
               for(nm in ls(orig)) {
                 aData[[nm]] <- orig[[nm]][i, j, ..., drop = drop]
               }
             }
             if ("lockedEnvironment" == storage.mode) assayDataEnvLock(aData)
             aData
           },
           list = {
             if (missing(i))                     # j must be present
               lapply(orig, function(obj) obj[, j, ..., drop = drop])
             else {                              # j may or may not be present
               if (missing(j))
                 lapply(orig, function(obj) obj[i,, ..., drop = drop])
               else
                 lapply(orig, function(obj) obj[i, j, ..., drop = drop])
             }
           })
  x@auxOrderingData = as.environment(as.list(x@auxOrderingData, all.names=TRUE))
  x@auxClusteringData = as.environment(as.list(x@auxClusteringData,
                                               all.names=TRUE))
  x
})


#' @rdname cell_data_set-methods
#' @aliases cell_data_set,ANY,ANY-method
setMethod("sizeFactors", signature(object="cell_data_set"), function(object) {
  sf <- pData(object)$Size_Factor
  names( sf ) <- colnames( exprs(object) )
  sf})

#' @rdname cell_data_set-methods
#' @aliases cell_data_set,ANY,ANY-method
#' @param value A vector of size factors, with length equal to the cells in
#'   object
setReplaceMethod("sizeFactors",
                 signature(object="cell_data_set", value="numeric"),
                 setSizeFactors <- function(object, value) {
  pData(object)$Size_Factor <- value
  validObject( object )
  object
})


#' @rdname cell_data_set-methods
#' @param locfunc A function applied to the geometric-mean-scaled expression
#'   values to derive the size factor.
#' @aliases cell_data_set,ANY,ANY-method
setMethod("estimateSizeFactors",
          signature(object="cell_data_set"),
          function( object, locfunc=median, ... )
          {
            sizeFactors(object) <- estimate_sf_matrix(exprs(object),
                                                                locfunc=locfunc,
                                                                ...)
            object
          })

#' @rdname cell_data_set-methods
#' @param modelFormulaStr A model formula, passed as a string, specifying how
#'   to group the cells prior to estimated dispersion. The default groups all
#'   cells together.
#' @param relative_expr Whether to transform expression into relative values
#' @param min_cells_detected Only include genes detected above
#'   lowerDetectionLimit in at least this many cells in the dispersion
#'   calculation
#' @param remove_outliers Whether to remove outliers (using Cook's distance)
#'   when estimating dispersions
#' @param cores The number of cores to use for computing dispersions
#' @aliases cell_data_set,ANY,ANY-method
setMethod("estimateDispersions",
          signature(object="cell_data_set"),
          function(object, modelFormulaStr="~ 1", relative_expr=TRUE,
                   min_cells_detected=1, remove_outliers=TRUE, cores=1,...)
          {
            dispModelName="blind"
            stopifnot( is( object, "cell_data_set" ) )

            if(!(identical("negbinomial.size",
                           object@expressionFamily) ||
                 identical("negbinomial", object@expressionFamily))) {
              stop("Error: estimateDispersions only works, and is only needed, when you're using a cell_data_set with a negbinomial or negbinomial.size expression family")
            }

            if( any( is.na( sizeFactors(object) ) ) )
              stop( "NAs found in size factors. Have you called 'estimateSizeFactors'?" )

            if( length(list(...)) != 0 )
              warning( "in estimateDispersions: Ignoring extra argument(s)." )

            # Remove results from previous fits
            object@dispFitInfo = new.env( hash=TRUE )

            dfi <- estimateDispersionsForcell_data_set(object,
                                                     modelFormulaStr,
                                                     relative_expr,
                                                     min_cells_detected,
                                                     remove_outliers,
                                                     cores)
            object@dispFitInfo[[dispModelName]] <- dfi

            validObject( object )
            object
          })

checkSizeFactors <- function(cds)
{
  if (cds@expressionFamily %in% c("negbinomial", "negbinomial.size"))
  {
    if (is.null(sizeFactors(cds))){
      stop("Error: you must call estimateSizeFactors() before calling this function.")
    }
    if (sum(is.na(sizeFactors(cds))) > 0){
      stop("Error: one or more cells has a size factor of NA.")
    }
  }
}



