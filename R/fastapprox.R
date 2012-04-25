##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Fast approximation
##' @param t1 
##' @param t2 
##' @param y1 
##' @author Klaus Kähler Holst
##' @export
fastapprox <- function(t1,t2,y1) {
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
