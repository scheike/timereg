##' @export
fast.cluster <- function(x,...) {
    arglist <- list("FastCluster",
                    time=as.integer(x),
                    PACKAGE="mets")
    res <- do.call(".Call",arglist)
    return(as.vector(res))
}
