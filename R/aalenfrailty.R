##' Aalen additive hazards model with (gamma) frailty
##'
##' Aalen frailty model
##' @title Aalen frailty model
##' @param time Time variable
##' @param status Status variable (0,1)
##' @param X Covariate design matrix
##' @param id cluster variable
##' @param theta list of thetas (returns score evaluated here), or
##' starting point for optimization (defaults to magic number 0.1)
##' @param ... Additional arguments to lower level functions
##' @return Parameter estimates
##' @author Klaus K. Holst
##' @export
##' @examples
##' dd <- simAalenFrailty(5000)
##' f <- ~1##+x
##' X <- model.matrix(f,dd) ## design matrix for non-parametric terms
##' system.time(out<-aalen(update(f,Surv(time,status)~.),dd,n.sim=0,robust=0))
##' dix <- which(dd$status==1)
##' t1 <- system.time(bb <- .Call("Bhat",as.integer(dd$status),
##'                               X,theta,as.integer(dd$id),NULL,NULL,-1,
##'                               package="mets"))
##' spec <- 1
##' ##plot(out,spec=spec)
##' plot(dd$time[dix],bb$B2[,spec],col="red",type="s",
##'      ylim=c(0,max(dd$time)*c(beta0,beta)[spec]))
##' abline(a=0,b=c(beta0,beta)[spec])
##' ##'
##' 
##' \dontrun{thetas <- seq(0.1,2,length.out=10)
##' Us <- unlist(aalenfrailty(dd$time,dd$status,X,dd$id,as.list(thetas)))
##' plot(thetas,Us,type="l",ylim=c(-.5,1)); abline(h=0,lty=2); abline(v=theta,lty=2)
##' op <- aalenfrailty(dd$time,dd$status,X,dd$id)
##' op}
aalenfrailty <- function(time,status,X,id,theta,B=NULL,...) {  
  dix <- which(status==1)
  cc <- cluster.index(id)
  ncluster <- length(cc$clusters)
  U <- function(theta,indiv=FALSE) {
    if (is.null(B)) {
        BB <- .Call("Bhat",as.integer(status),X,theta,as.integer(cc$clusters),cc$idclust,as.integer(cc$cluster.size))$B2
    } else {
        BB <- B*time[dix]
    }
    Hij0 <- apply(X[dix,,drop=FALSE]*BB,1,sum)
    Hij <- Cpred(cbind(time[dix],Hij0),time)[,2,drop=FALSE]
##    if (is.na(Hij[1])) browser()
    res <- .Call("Uhat",as.integer(status),Hij,theta,cc$idclust,as.integer(cc$cluster.size))
    if (!indiv) res <- mean(res,na.rm=TRUE)
    res
  }
  if (missing(theta)) theta <- 0.1
  if (is.list(theta)) return(lapply(theta,function(x) U(x,...)))

  op <- nlminb(theta,function(x) U(x)^2)
  uu <- U(op$par,TRUE)
  du <- numDeriv::grad(U,op$par)
  return(list(theta=op$par, sd=(mean(uu^2)/du^2/ncluster)^0.5))
}
                  

##' Simulate observations from Aalen Frailty model with Gamma
##' distributed frailty and constant intensity.
##'
##' @title Simulate from the Aalen Frailty model
##' @param K Number of clusters
##' @param n Number of observations in each cluster
##' @param eta 1/variance
##' @param beta Effect (log hazard ratio) of covariate
##' @param stoptime Stopping time
##' @param left Left truncation
##' @param pairleft pairwise (1) left truncation or individual (0)
##' @author Klaus K. Holst
##' @export
simAalenFrailty <- function(n=5e3,theta=0.3,K=2,beta0=1.5,beta=1,cens=1.5,...) {
    ## beta0 (constant baseline intensity)
    ## beta (covariate effect)
    Z <- rep(rgamma(n/K,1/theta,1/theta),each=K) ## Frailty, mean 1, var theta
    id <- rep(seq(n/K),each=K) ## Cluster indicator
    Ai <- function(u) (u/Z)/(beta0+beta*x) 
    x <- rbinom(n,1,0.5) ## Binary covariate
    uu <- -log(runif(n)) ##rexp(n,1)   
    dat <- data.frame(time=Ai(uu), x=x, status=1, id=id, Z=Z)
    if (cens==0) cens <- Inf else cens <- -log(runif(n))/cens
    dat$status <- (dat$time<=cens)*1
    dat$time <- apply(cbind(dat$time,cens),1,min)
    dat <- dat[order(dat$time),] ## order after event/censoring tim
    return(dat)
}
