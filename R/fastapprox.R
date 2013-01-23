##' Fast approximation
##'
##' @title Fast approximation
##' @param time Original ordered time points
##' @param new.time New time points 
##' @param y Optional function values (default is simply the identity) evaluated at the original time points \code{time}
##' @param ... Optional additional arguments
##' @author Klaus K. Holst
##' @return List with the elements
##' 
##' \item{approx}{Approximated function values (\code{y1}) evaluated in \code{t2}}
##' 
##' \item{pos}{The positions in \code{t1} of the elements of \code{t2} starting from index 0}
##' 
##' @export
fast.approx <- function(time,new.time,y=time,...) {
  if (is.matrix(time)) {
    y <- time[,-1,drop=TRUE]; time <- time[,1,drop=TRUE]
  }    
  arglist <- list("FastApprox",
                  a=time,
                  t=new.time,
                  z=y,
                  DUP=FALSE,PACKAGE="mets")
  res <- do.call(".Call",arglist)
  return(res)
}
