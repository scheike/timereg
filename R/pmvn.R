##' @export
pbvn <- function(upper,rho,sigma) {
    if (!missing(sigma)) {
        rho <- cov2cor(sigma)[1,2]
        upper <- upper/diag(sigma)^0.5
    }
    arglist <- list("bvncdf",
                    a=upper[1],
                    b=upper[2],
                    r=rho,
                    DUP=FALSE,PACKAGE="mets")
    res <- do.call(".Call",arglist)
    return(res)
}

##' @export
pmvn <- function(lower,upper,mu=rep(0,ncol(sigma)),sigma,notcor=TRUE) {
    if (missing(sigma)) stop("Specify variance matrix 'sigma'")
    if (missing(lower)) {
        if (missing(upper)) stop("Lower or upper integration bounds needed")
        lower <- upper; lower[] <- -Inf
    }
    p <- ncol(rbind(mu))
    if (missing(upper)) {
        upper <- lower; upper[] <- Inf
    }
    if (ncol(rbind(sigma))!=p)
        stop("Incompatible dimensions of mean and variance")
    if (ncol(rbind(lower))!=p || ncol(rbind(upper))!=p)
        stop("Incompatible integration bounds")    
    arglist <- list("pmvn",
                    lower=rbind(lower),
                    upper=rbind(upper),
                    mu=rbind(mu),
                    sigma=rbind(sigma),
                    notcor=as.logical(notcor[1]))
    res <- do.call(".Call",arglist)
    return(as.vector(res))
}
