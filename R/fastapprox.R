##' Fast approximation
##'
##' @title Fast approximation
##' @param time Original ordered time points
##' @param new.time New time points 
##' @param equal If TRUE a list is returned with additional element
##' @param right If FALSE the closest element in time is chosen,
##' otherwise the closest value above new.time is returned.
##' @param ... Optional additional arguments
##' @author Klaus K. Holst
##' @examples
##' id <- c(1,1,2,2,7,7,10,10)
##' fast.approx(unique(id),id)
##' 
##' t <- 0:6
##' n <- c(-1,0,0.1,0.9,1,1.1,1.2,6,6.5)
##' fast.approx(t,n,equal=TRUE)
##' @export
fast.approx <- function(time,new.time,equal=FALSE,right=FALSE,...) {
    ## if (sort) {
    ##     ord <- order(time)
    ##     time <- time[ord]
    ## }
    arglist <- list("FastApprox",
                    time=sort(time),
                    newtime=new.time,
                    equal=equal,
                    right=right,
                    PACKAGE="mets")
    res <- do.call(".Call",arglist)
    return(res)
}
