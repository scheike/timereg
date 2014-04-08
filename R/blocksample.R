##' Sample blockwise from clustered data
##' 
##' @title Block sampling
##' @param data Data frame
##' @param idvar Column defining the clusters
##' @param size Size of samples
##' @param replace Logical indicating wether to sample with replacement
##' @param \dots additional arguments to lower level functions
##' @return \code{data.frame}
##' @author Klaus K. Holst
##' @keywords models utilities
##' @export
##' @examples
##' 
##' dd <- data.frame(x=rnorm(6), z=rnorm(6), id=c(4,4,10,10,5,5), v=rnorm(6))
##' blocksample(dd,idvar="id",ids=c(4,5))
##' blocksample(dd,size=10)
##'
blocksample <- function(data, size, idvar="id", replace=TRUE, ids, ...) {
  if (is.matrix(data))    
      data <- as.data.frame(data)
  if (length(idvar)==1 && is.character(idvar)) {
      idvar <- data[,idvar]
  }
  ii <-  cluster.index(idvar)
  id0 <- idvar[ii$firstclust+1]
  ids <- sample(seq(ii$uniqueclust), size=ifelse(missing(size),length(idvar),size),replace=replace)
  idx <- na.omit(as.vector(t(ii$idclustmat[ids,])))+1
  len <- unlist(lapply(ids,function(x) rep(x,ii$cluster.size[x])))
  
  return(data[idx,])
}
