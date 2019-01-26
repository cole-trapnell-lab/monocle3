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
#' @import Matrix
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
  x@aux_ordering_data = as.environment(as.list(x@aux_ordering_data, all.names=TRUE))
  x@aux_clustering_data = as.environment(as.list(x@aux_clustering_data,
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
#'   lower_detection_limit in at least this many cells in the dispersion
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
                           object@expression_family) ||
                 identical("negbinomial", object@expression_family))) {
              stop(paste("Error: estimateDispersions only works, and is only",
                         "needed, when you're using a cell_data_set with a",
                         "negbinomial or negbinomial.size expression family"))
            }

            if( any( is.na( sizeFactors(object) ) ) )
              stop(paste("NAs found in size factors. Have you called",
                         "'estimateSizeFactors'?"))

            if( length(list(...)) != 0 )
              warning( "in estimateDispersions: Ignoring extra argument(s).")

            # Remove results from previous fits
            object@disp_fit_info = new.env( hash=TRUE )

            dfi <- estimate_dispersions_cds(object,
                                            modelFormulaStr,
                                            relative_expr,
                                            min_cells_detected,
                                            remove_outliers,
                                            cores)
            object@disp_fit_info[[dispModelName]] <- dfi

            validObject( object )
            object
          })

checkSizeFactors <- function(cds) {
  if (cds@expression_family %in% c("negbinomial", "negbinomial.size"))
  {
    if (is.null(sizeFactors(cds))){
      stop(paste("Error: you must call estimateSizeFactors() before calling",
                 "this function."))
    }
    if (sum(is.na(sizeFactors(cds))) > 0){
      stop("Error: one or more cells has a size factor of NA.")
    }
  }
}


#' Retrieves the coordinates of each cell in the reduced-dimensionality space generated by calls to
#' reduceDimension.
#'
#' Reducing the dimensionality of the expression data is a core step in the Monocle
#' workflow. After you call reduceDimension(), this function will return the new
#' coordinates of your cells in the reduced space.
#' @param cds A cell_data_set object.
#' @return A matrix, where rows are cell coordinates and columns correspond to dimensions of the
#' reduced space.
#' @export
#' @examples
#' \dontrun{
#' S <- reducedDimS(HSMM)
#' }
reducedDimS <- function( cds ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimS
}

#' Set embedding coordinates of each cell in a cell_data_set
#'
#' This function sets the coordinates of each cell in a new
#' (reduced-dimensionality) space. Not intended to be called directly.
#'
#' @param cds A cell_data_set object.
#' @param value A matrix of coordinates specifying each cell's position in the reduced-dimensionality space.
#' @return An update cell_data_set object
#' @examples
#' \dontrun{
#' cds <- reducedDimS(S)
#' }
`reducedDimS<-` <- function( cds, value ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimS <- value
  validObject( cds )
  cds
}

#' Get the whitened expression values for a cell_data_set
#'
#' Retrieves the expression values for each cell (as a matrix) after whitening
#' during dimensionality reduction.
#'
#' @param cds A cell_data_set object.
#' @return A matrix, where each row is a set of whitened expression values for a feature and columns are cells.
#' @export
#' @examples
#' \dontrun{
#' W <- reducedDimW(HSMM)
#' }
reducedDimW <- function( cds ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimW
}

#' @title Sets the the whitening matrix during independent component analysis.
#'
#' @description Sets the the whitening matrix during independent component analysis.
#'
#' @param cds A cell_data_set object.
#' @param value a numeric matrix
#' @return A matrix, where each row is a set of whitened expression values for a feature and columns are cells.
#' @docType methods
#' @examples
#' \dontrun{
#' cds <- reducedDimK(K)
#' }
`reducedDimK<-` <- function( cds, value ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimK <- value
  validObject( cds )
  cds
}

#' Retrieves the the whitening matrix during independent component analysis.
#'
#' @param cds A cell_data_set object.
#' @return A matrix, where each row is a set of whitened expression values for a feature and columns are cells.
#' @docType methods
#' @export
#' @examples
#' \dontrun{
#' K <- reducedDimW(HSMM)
#' }
reducedDimK <- function( cds ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimK
}

#' Sets the whitened expression values for each cell prior to independent component analysis. Not intended to be called directly.
#'
#' @param cds A cell_data_set object.
#' @param value A whitened expression data matrix
#' @return An updated cell_data_set object
#' @examples
#' \dontrun{
#' #' cds <- reducedDimA(A)
#' }
`reducedDimW<-` <- function( cds, value ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimW <- value
  validObject( cds )
  cds
}

#' Get the weights needed to lift cells back to high dimensional expression space.
#'
#' Retrieves the weights that transform the cells' coordinates in the reduced
#' dimension space back to the full (whitened) space.
#'
#' @param cds A cell_data_set object.
#' @return A matrix that when multiplied by a reduced-dimension set of coordinates for the cell_data_set,
#' recovers a matrix in the full (whitened) space
#' @export
#' @examples
#' \dontrun{
#' A <- reducedDimA(HSMM)
#' }
reducedDimA <- function( cds ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimA
}

#' Get the weights needed to lift cells back to high dimensional expression space.
#'
#' Sets the weights transform the cells' coordinates in the reduced dimension
#' space back to the full (whitened) space.
#'
#' @param cds A cell_data_set object.
#' @param value A whitened expression data matrix
#' @return An updated cell_data_set object
#' @export
#' @examples
#' \dontrun{
#' cds <- reducedDimA(A)
#' }
`reducedDimA<-` <- function( cds, value ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@reducedDimA <- value
  validObject( cds )
  cds
}

#' Retrieves the minimum spanning tree generated by Monocle during cell ordering.
#'
#' Retrieves the minimum spanning tree (MST) that Monocle constructs during order_cells().
#' This MST is mostly used in plot_spanning_tree to help assess the accuracy
#' of Monocle\'s ordering.
#' @param cds expression data matrix for an experiment
#' @return An igraph object representing the cell_data_set's minimum spanning tree.
#' @export
#' @examples
#' \dontrun{
#' T <- principal_graph(HSMM)
#' }
principal_graph <- function( cds ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@principal_graph
}

#' Set the minimum spanning tree generated by Monocle during cell ordering.
#'
#' Sets the minimum spanning tree used by Monocle during cell ordering. Not intended to be called directly.
#'
#' @param cds A cell_data_set object.
#' @param value an igraph object describing the minimum spanning tree.
#' @return An updated cell_data_set object
#' @export
#' @examples
#' \dontrun{
#' cds <- principal_graph(T)
#' }
`principal_graph<-` <- function( cds, value ) {
  stopifnot( is( cds, "cell_data_set" ) )
  cds@principal_graph <- value
  validObject( cds )
  cds
}

