##' @export
or2prob <- function(OR,marg) {
    p1 <- marg[1]; p2 <- marg[2]
    if (OR==1) {
        PP <- outer(c(p1,1-p1),c(p2,1-p2))
    } else {
        b <- 1+(p1+p2)*(OR-1)
        a <- (1-OR)
        ac <- -OR*(1-OR)*p1*p2
        d <- sqrt(b^2-4*ac)
        ##p11b <- (-b-d)/(2*a)
        p11 <- (-b+d)/(2*a)
        PP <- c(p11,p1-p11,p2-p11)
        PP <- c(PP,1-sum(PP))
    }
    structure(matrix(PP,2),marg=c(p1,p2))
}


##' @export
tetrachoric <- function(P,OR,approx=0,...) {
    if (!missing(OR)) {
        ## Assuming P[1],P[2] is the marginals
        P <- or2prob(OR,P)
        p1 <- attributes(P)$marg[1]
        p2 <- attributes(P)$marg[2]
    } else {
        ## Assuming P contains the joint probabilities
        if (is.vector(P)) {
            if (length(P)==3) P <- c(P,1-sum(P))
            P <- matrix(P,2)
        }
        if (!all.equal(sum(P),1)) stop("Not a probability matrix")
        p1 <- colSums(P)[1]
        p2 <- rowSums(P)[1]
    }
    if (approx>0) {
        k <- (1-abs(p1-p2)/5 - (.5-min(p1,p2))^2)/2
        if (missing(OR)) OR <- prod(diag(P))/prod(revdiag(P))
        return(cos(pi/(1+OR^k)))
    }    
    q1 <- qnorm(p1)
    q2 <- qnorm(p2)
    lo <- rbind(c(0,0),c(0,-Inf),c(-Inf,0),c(-Inf,-Inf))
    hi <- rbind(c(Inf,Inf),c(Inf,0),c(0,Inf),c(0,0))
    mu <- cbind(q1,q2)%x%cbind(rep(1,4))
    obj <- function(r) {
        Pr <- pmvn(lower=lo,upper=hi,mu=mu,sigma=r,cor=TRUE)
        return(mean(abs(P-Pr)^2))
        ##(P[1,1]-pmvn(lower=c(0,0),mu=c(q1,q2),sigma=r,cor=TRUE))^2
    }
    optimize(obj,interval=c(-1,1))$minimum
}


assoc <- function(x,...) {
    N <- sum(x)
    P <- x
    if (N!=1) P <- P/N
    p1 <- colSums(P)[1]
    p2 <- rowSums(P)[1]
    OR <- prod(diag(P))/prod(revdiag(P))
    rho <- tetrachoric(P)
    list(P=P, OR=OR, rho=rho)
}
