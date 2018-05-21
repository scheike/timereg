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
                    PACKAGE="mets")
    res <- do.call(".Call",arglist)
    return(res)
}


##' Multivariate normal distribution function
##'
##' Multivariate normal distribution function
##' @aliases pmvn pbvn loglikMVN scoreMVN dmvn rmvn
##' @export
##' @examples
##' lower <- rbind(c(0,-Inf),c(-Inf,0))
##' upper <- rbind(c(Inf,0),c(0,Inf))
##' mu <- rbind(c(1,1),c(-1,1))
##' sigma <- diag(2)+1
##' pmvn(lower=lower,upper=upper,mu=mu,sigma=sigma)
##' @param lower lower limits
##' @param upper upper limits
##' @param mu mean vector
##' @param sigma variance matrix or vector of correlation coefficients
##' @param cor if TRUE sigma is treated as standardized (correlation matrix)
pmvn <- function(lower,upper,mu,sigma,cor=FALSE) {
    if (missing(sigma)) stop("Specify variance matrix 'sigma'")
    if (missing(lower)) {
        if (missing(upper)) stop("Lower or upper integration bounds needed")
        lower <- upper; lower[] <- -Inf
    }
    p <- ncol(rbind(lower))
    if (missing(upper)) {
        upper <- lower; upper[] <- Inf
    }
    if (missing(mu)) mu <- rep(0,p)
    sigma <- rbind(sigma)
    ncor <- p*(p-1)/2
    if (ncol(sigma)!=p && ncol(sigma)!=ncor)
        stop("Incompatible dimensions of mean and variance")    
    if (ncol(rbind(lower))!=p || ncol(rbind(upper))!=p)
        stop("Incompatible integration bounds")    
    arglist <- list("pmvn0",
                   lower=rbind(lower),
                   upper=rbind(upper),
                   mu=rbind(mu),
                   sigma=rbind(sigma),
                   cor=as.logical(cor[1]),
                   PACKAGE="mets")
    res <- do.call(".Call",arglist)
    return(as.vector(res))
}


introotpn <- function(p) {
    ## Find integer root of x^2-x-2*p=0
    n <- 0.5*(1+sqrt(1+8*p))
    if (floor(n)!=n) n <- NA
    return(n)
}

##' @export
rmvn <- function(n,mu,sigma,rho,...) {
    if (!missing(rho)) {
        if (is.vector(rho)) rho <- rbind(rho)
        if (missing(mu)) {
            p <- introotpn(NCOL(rho))
            mu <- rep(0,p)
        }
        return (.Call("_mets_rmvn",
                 n=as.integer(n),
                 mu=rbind(mu),
                 rho=rbind(rho)))
    }
    if (!missing(mu) && missing(sigma)) sigma <- diag(nrow=length(mu))
    if (missing(sigma)) sigma <- matrix(1)
    if (is.vector(sigma)) sigma <- diag(sigma,ncol=length(sigma))
    if (missing(mu)) mu <- rep(0,ncol(sigma))    
    PP <- with(svd(sigma), v%*%diag(sqrt(d),ncol=length(d))%*%t(u))
    res <- matrix(rnorm(ncol(sigma)*n),ncol=ncol(sigma))%*%PP
    if (NROW(mu)==nrow(res) && NCOL(mu)==ncol(res)) return(res+mu)
    return(res+cbind(rep(1,n))%*%mu)
}

##' @export
dmvn <- function(x,mu,sigma,rho,log=FALSE,nan.zero=TRUE,...) {    
    if (!missing(rho)) {
        if (is.vector(rho)) rho <- rbind(rho)
        if (is.vector(x)) x <- rbind(x)
        if (missing(mu)) {
            p <- NCOL(x)
            mu <- rep(0,p)
        }
        res <- .Call("_mets_dmvn",
                    u=x,
                    mu=rbind(mu),
                    rho=rho)
        if (!log) res <- exp(res)
        return(res)
    }
    if (!missing(mu) && missing(sigma)) sigma <- diag(nrow=length(mu))
    if (missing(sigma)) sigma <- matrix(1)
    if (is.vector(sigma)) sigma <- diag(sigma,ncol=length(sigma))
    if (missing(mu)) mu <- rep(0,ncol(sigma))
    
    if (length(sigma)==1) {
        k <- 1
        isigma <- structure(cbind(1/sigma),det=as.vector(sigma))

    } else {
        k <- ncol(sigma)
        isigma <- Inverse(sigma)
    }
    if (!missing(mu)) {
        if (NROW(mu)==NROW(x) && NCOL(mu)==NCOL(x)) {
            x <- x-mu
        } else {
            x <- t(t(x)-mu)
        }
    }
    logval <- -0.5*(base::log(2*base::pi)*k+
                    base::log(attributes(isigma)$det)+
                    rowSums((x%*%isigma)*x))
    if (nan.zero) logval[is.nan(logval)] <- -Inf
    if (log) return(logval)
    return(exp(logval))
}
