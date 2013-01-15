##' Fast pattern
##'
##' @title Fast pattern
##' @param x Matrix (binary) of patterns. Optionally if \code{y} is
##' also passed as argument, then the pattern matrix is defined as the
##' elements agreeing in the two matrices.
##' @param y Optional matrix argument with same dimensions as
##' \code{x} (see above)
##' @param ... Optional additional arguments
##' @author Klaus K. Holst
##' @export
fast.pattern <- function(x,y,...) {
  if (missing(y)) y <- NULL
   .Call("FastPattern",x,y,DUP=FALSE)
}
 
