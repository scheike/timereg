##' Fast approximation
##'
##' @title Fast approximation
##' @param t1 Original ordered time points
##' @param t2 New time points
##' @param y1 Optional function values (default is simply the identity) corresponding to the original time points \code{t1}
##' @param ... Optional additional arguments
##' @author Klaus K. Holst
##' @return List with the elements
##' 
##' \item{approx}{Approximated function values (\code{y1}) evaluated in \code{t2}}
##' 
##' \item{pos}{The positions in \code{t1} of the elements of \code{t2}}
##' 
##' @export
fastapprox <- function(t1,t2,y1=t1,...) {
  if (is.matrix(t1)) {
    y <- t1[,-1]; t1 <- t1[,1]
  }    
  arglist <- list(name="FastApprox",
                  a=t1,
                  t=y1,
                  z=t2,
                  DUP=FALSE)##,PACKAGE="bptwin")  
  res <- do.call(".Call",arglist)
  return(res)
}
