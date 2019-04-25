#' @export
setGeneric("partitions", function(x, ...) standardGeneric("partitions"))

#' @export
setGeneric("partitions<-", function(x, value) standardGeneric("partitions<-"))

#' @export
setGeneric("partitionNames", function(x) standardGeneric("partitionNames"))

#' @export
setGeneric("partitionNames<-", function(x, value) standardGeneric("partitionNames<-"))

# Getter/setter functions for partitions
# Copied from reduceDims functions from SingleCellExperiment

#' @export
setMethod("partitions", "cell_data_set", function(x) {
  value <- x@partitions
  return(value)
})

#' @export
#' @importClassesFrom S4Vectors List
setReplaceMethod("partitions", "cell_data_set", function(x, value) {
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  x@partitions <- value
  validObject(x)
  return(x)
})

#' @export
setMethod("partitionNames", "cell_data_set", function(x) {
  names(partitions(x))
})


#' @export
setReplaceMethod("partitionNames", c("cell_data_set", "character"), function(x, value) {
  out <- partitions(x)
  names(out) <- value
  x@partitions <- out
  return(x)
})






#' @export
setGeneric("principal_graph", function(x, ...) standardGeneric("principal_graph"))

#' @export
setGeneric("principal_graph<-", function(x, value) standardGeneric("principal_graph<-"))

# Getter/setter functions for partitions
# Copied from reduceDims functions from SingleCellExperiment

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
setGeneric("principal_graph_aux", function(x, ...) standardGeneric("principal_graph_aux"))

#' @export
setGeneric("principal_graph_aux<-", function(x, value) standardGeneric("principal_graph_aux<-"))


# Getter/setter functions for partitions
# Copied from reduceDims functions from SingleCellExperiment

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
setGeneric("exprs", function(x, ...) standardGeneric("exprs"))

# Getter/setter wrappers for exprs
# Copied from reduceDims functions from SingleCellExperiment

#' @export
setMethod("exprs", "cell_data_set", function(x) {
  value <- assays(x)$counts
  return(value)
})

#' @export
setGeneric("pData", function(x, ...) standardGeneric("pData"))

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
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  colData(x) <- value
  validObject(x)
  return(x)
})

#' @export
setGeneric("fData", function(x, ...) standardGeneric("fData"))

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
  value <- methods::as(value, "List")
  if (is.null(names(value))) {
    names(value) <- character(length(value))
  }
  rowData(x) <- value
  validObject(x)
  return(x)
})
