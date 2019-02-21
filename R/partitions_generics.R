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
