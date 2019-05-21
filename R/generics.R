#' Generic to extract clusters from CDS object
#'
#' @param x CDS object
#' @param reduction_method Reduced dimension to extract clusters for.
#'
#' @export
setGeneric("get_clusters", function(x, reduction_method = "UMAP")
  standardGeneric("get_clusters"))

#' Method to extract clusters from CDS object
#' @param cell_data_set CDS object
#'
#' @export
setMethod("get_clusters", "cell_data_set",
          function(x, reduction_method = "UMAP") {
  value <- x@clusters[[reduction_method]]$clusters
  return(value)
})

#' @export
setGeneric("clusters", function(x) standardGeneric("clusters"))

#' @export
setGeneric("clusters<-", function(x, value) standardGeneric("clusters<-"))

#' @export
setGeneric("clusterNames", function(x) standardGeneric("clusterNames"))

#' @export
setGeneric("clusterNames<-", function(x, value)
  standardGeneric("clusterNames<-"))

# Getter/setter functions for clusters
# Copied from reduceDims functions from SingleCellExperiment

#' @export
setMethod("clusters", "cell_data_set", function(x) {
  value <- x@clusters
  return(value)
})

#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("clusters", "cell_data_set", function(x, value) {
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@clusters <- value
  validObject(x)
  return(x)
})

#' @export
setMethod("clusterNames", "cell_data_set", function(x) {
  names(clusters(x))
})


#' @export
setReplaceMethod("clusterNames", c("cell_data_set", "character"),
                 function(x, value) {
  out <- clusters(x)
  names(out) <- value
  x@clusters <- out
  return(x)
})

#' Generic to extract partitions from CDS object
#'
#' @param x CDS object
#' @param reduction_method Reduced dimension to extract clusters for.
#'
#' @export
setGeneric("get_partitions", function(x, reduction_method = "UMAP")
  standardGeneric("get_partitions"))

#' Method to extract partitions from CDS object
#' @param cell_data_set CDS object
#'
#' @export
setMethod("get_partitions", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@clusters[[reduction_method]]$partitions
            return(value)
          })

#' @export
setGeneric("partitions", function(x) standardGeneric("partitions"))

#' @export
setGeneric("principal_graph", function(x)
  standardGeneric("principal_graph"))

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

#' @export
setGeneric("exprs", function(x) standardGeneric("exprs"))

# Getter/setter wrappers for exprs
# Copied from reduceDims functions from SingleCellExperiment

#' @export
setMethod("exprs", "cell_data_set", function(x) {
  value <- assays(x)$counts
  return(value)
})


#' @export
setGeneric("pData", function(x) standardGeneric("pData"))

#' @export
setGeneric("pData<-", function(x, value) standardGeneric("pData<-"))

# Getter/setter wrappers for pdata
# Copied from reduceDims functions from SingleCellExperiment

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

# Getter/setter wrappers for pdata
# Copied from reduceDims functions from SingleCellExperiment

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
