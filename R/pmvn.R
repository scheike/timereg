pmvn <- function(upper,rho,sigma) {
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
