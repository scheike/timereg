##' Fast approximation
##'
##' @title Fast approximation
##' @param time Original ordered time points
##' @param new.time New time points 
##' @param equal If TRUE a list is returned with additional element
##' @param ... Optional additional arguments
##' @author Klaus K. Holst
##' @export
fast.approx <- function(time,new.time,equal=FALSE,...) {
  arglist <- list("FastApprox",
                  time=time,
                  newtime=new.time,
                  equal=equal,
                  DUP=FALSE,PACKAGE="mets")
  res <- do.call(".Call",arglist)
  return(res)
}
