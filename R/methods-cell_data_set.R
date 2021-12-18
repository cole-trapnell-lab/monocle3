#' Methods for the cell_data_set class
#' @name cell_data_set-methods
#' @docType methods
#' @rdname cell_data_set-methods
#' @param object The cell_data_set object
setValidity( "cell_data_set", function( object ) {
  TRUE
} )

#' Set the size factor values in the cell_data_set
#'
#' @param cds A cell_data_set object.
#' @param value the size factor values.
#' @return An updated cell_data_set object
`size_factors<-` <- function( cds, value ) {
  stopifnot( methods::is( cds, "cell_data_set" ) )
  colData(cds)$Size_Factor <- value
  methods::validObject( cds )
  cds
}

#' Get the size factors from a cds object.
#'
#' A wrapper around \code{colData(cds)$Size_Factor}
#'
#' @param cds A cell_data_set object.
#' @return An updated cell_data_set object
#' @export
#' @examples
#'   cds <- load_a549()
#'   size_factors(cds)
size_factors <- function( cds ) {
  stopifnot( methods::is( cds, "cell_data_set" ) )
  sf <- colData(cds)$Size_Factor
  names( sf ) <- colnames( SingleCellExperiment::counts(cds) )
  sf
}

