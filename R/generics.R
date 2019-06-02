#' Generic to extract pseudotime from CDS object
#'
#' @param x CDS object
#' @param reduction_method Reduced dimension to extract pseudotime for.
#'
#' @export
setGeneric("pseudotime", function(x, reduction_method = "UMAP")
  standardGeneric("pseudotime"))

#' Method to extract pseudotime from CDS object
#' @param cell_data_set CDS object
#'
#' @export
setMethod("pseudotime", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@principal_graph_aux[[reduction_method]]$pseudotime
            if (is.null(value)) {
              stop(paste0("No pseudotime calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "order_cells with reduction_method = ",
                          reduction_method, "."))
            }
            return(value)
          })

#' Generic to extract clusters from CDS object
#'
#' @param x CDS object
#' @param reduction_method Reduced dimension to extract clusters for.
#' @export
setGeneric("clusters", function(x, reduction_method = "UMAP")
  standardGeneric("clusters"))

#' Method to extract clusters from CDS object
#' @param cell_data_set CDS object
#'
#' @export
setMethod("clusters", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@clusters[[reduction_method]]$clusters
            if (is.null(value)) {
              stop(paste0("No clusters calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "cluster_cells with reduction_method = ",
                          reduction_method, "."))
            }
            return(value)
          })

#' Generic to extract partitions from CDS object
#'
#' @param x CDS object
#' @param reduction_method Reduced dimension to extract clusters for.
#'
#' @export
setGeneric("partitions", function(x, reduction_method = "UMAP")
  standardGeneric("partitions"))

#' Method to extract partitions from CDS object
#' @param cell_data_set CDS object
#'
#' @export
setMethod("partitions", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@clusters[[reduction_method]]$partitions
            if (is.null(value)) {
              stop(paste0("No partitions calculated for reduction_method = ",
                          reduction_method, ". Please first run ",
                          "cluster_cells with reduction_method = ",
                          reduction_method, "."))
            }
            return(value)
          })

#' @export
setGeneric("principal_graph", function(x) standardGeneric("principal_graph"))

#' @export
setGeneric("principal_graph<-", function(x, value)
  standardGeneric("principal_graph<-"))

#' @export
setMethod("principal_graph", "cell_data_set", function(x) {
  value <- x@principal_graph
  return(value)
})

#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("principal_graph", "cell_data_set", function(x, value) {
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@principal_graph <- value
  validObject(x)
  return(x)
})


#' @export
setGeneric("principal_graph_aux", function(x)
  standardGeneric("principal_graph_aux"))

#' @export
setGeneric("principal_graph_aux<-", function(x, value)
  standardGeneric("principal_graph_aux<-"))

#' @export
setMethod("principal_graph_aux", "cell_data_set", function(x) {
  value <- x@principal_graph_aux
  return(value)
})

#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("principal_graph_aux", "cell_data_set", function(x, value) {
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@principal_graph_aux <- value
  validObject(x)
  return(x)
})


# Set of wrappers for easy transition from monocle.

#' @export
setGeneric("exprs", function(x) standardGeneric("exprs"))

#' @export
setMethod("exprs", "cell_data_set", function(x) {
  value <- assays(x)$counts
  return(value)
})

#' @export
setGeneric("pData", function(x) standardGeneric("pData"))

#' @export
setGeneric("pData<-", function(x, value) standardGeneric("pData<-"))

#' @export
setMethod("pData", "cell_data_set", function(x) {
  value <- colData(x)
  return(value)
})

#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("pData", "cell_data_set", function(x, value) {
  colData(x) <- value
  validObject(x)
  return(x)
})

#' @export
setGeneric("fData", function(x) standardGeneric("fData"))

#' @export
setGeneric("fData<-", function(x, value) standardGeneric("fData<-"))

#' @export
setMethod("fData", "cell_data_set", function(x) {
  value <- rowData(x)
  return(value)
})

#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("fData", "cell_data_set", function(x, value) {
  rowData(x) <- value
  validObject(x)
  return(x)
})
