# Rescued from archived grr 0.9.5 package.
#
#'Extract/return parts of objects
#'
#'Alternative to built-in \code{\link{Extract}} or \code{[}.  Allows for
#'extraction operations that are ambivalent to the data type of the object. For
#'example, \code{extract(x,i)} will work on lists, vectors, data frames,
#'matrices, etc.
#'
#'Extraction is 2-100x faster on data frames than with the built in operation -
#'but does not preserve row names.
#'
#'@param x object from which to extract elements
#'@param i,j indices specifying elements to extract.  Can be \code{numeric},
#'  \code{character}, or \code{logical} vectors.
#'@export
#'@examples
#'#Typically about twice as fast on normal subselections
#'orders<-data.frame(orderNum=1:1e5,
#'  sku=sample(1e3, 1e5, TRUE),
#'  customer=sample(1e4,1e5,TRUE))
#'a<-sample(1e5,1e4)
#' system.time(b<-orders[a,])
#'system.time(c<-extract(orders,a))
#'rownames(b)<-NULL
#'rownames(c)<-NULL
#'identical(b,c)
#'
#'#Speedup increases to 50-100x with oversampling
#'a<-sample(1e5,1e6,TRUE)
#' system.time(b<-orders[a,])
#'system.time(c<-extract(orders,a))
#'rownames(b)<-NULL
#'rownames(c)<-NULL
#'identical(b,c)
#'
#'#Can create function calls that work for multiple data types
#'alist<-as.list(1:50)
#'avector<-1:50
#'extract(alist,1:5)
#'extract(avector,1:5)
#'extract(orders,1:5)#'
#'
#'\dontrun{
#'orders<-data.frame(orderNum=as.character(sample(1e5, 1e6, TRUE)),
#'  sku=sample(1e3, 1e6, TRUE),
#'  customer=sample(1e4,1e6,TRUE))
#'system.time(a<-sample(1e6,1e7,TRUE))
#' system.time(b<-orders[a,])
#'system.time(c<-extract(orders,a))
#'}
grr_extract<-function(x,i=NULL,j=NULL)
{
  if(is.null(dim(x)))
  {
     x<-x[i]
    return(x)
  }
  else
    if(!is.null(j))
      x<-x[,j]
  if(!is.null(i))
  {
  if(is.data.frame(x))
    x<-as.data.frame(lapply(x,function (a) a[i]))
  else
    x<-x[i,]
  }
  return(x)
}

