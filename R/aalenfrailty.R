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
##' n <- 5e3
##' theta <- 0.3
##' K <- 2
##' Z <- rep(rgamma(n/K,1/theta,1/theta),each=K)
##' id <- rep(seq(n/K),each=K)
##' beta0 <- 1.5; beta <- 0##-0.5 #-.15 
##' Ai <- function(u) (u/Z)/(beta0+beta*x)
##' x <- rbinom(n,1,0.5)
##' uu <- -log(runif(n)) ##rexp(n,1)
##' dd <- data.frame(time=Ai(uu), x=x, status=1, id=id, Z=Z)
##' cens <- runif(n,10,20)
##' dd$status <- (dd$time<=cens)*1
##' dd$time <- apply(cbind(dd$time,cens),1,min)
##' dd <- dd[order(dd$time),] ## order after event/censoring tim
##' f <- ~1##+x
##' X <- model.matrix(f,dd) ## design matrix for non-parametric terms
##' system.time(out<-aalen(update(f,Surv(time,status)~.),dd,n.sim=0,robust=0))
##' dix <- which(dd$status==1)
##' t1 <- system.time(bb <- .Call("Bhat",dd$status,X,theta,dd$id,NULL,NULL,-1,package="mets"))
##' spec <- 1
##' ##plot(out,spec=spec)
##' plot(dd$time[dix],bb$B2[,spec],col="red",type="s",ylim=c(0,max(dd$time)*c(beta0,beta)[spec]))
##' abline(a=0,b=c(beta0,beta)[spec])
##'
##' 
##' thetas <- seq(0.1,2,length.out=10)
##' Us <- unlist(aalenfrailty(dd$time,dd$status,X,dd$id,as.list(thetas)))
##' plot(thetas,Us,type="l",ylim=c(-.5,1)); abline(h=0,lty=2); abline(v=theta,lty=2)
##' op <- aalenfrailty(dd$time,dd$status,X,dd$id)
##' op
aalenfrailty <- function(time,status,X,id,theta,
                         ...) {  
  dix <- which(status==1)
  fB <- fastapprox(time[dix],time)
  cc <- cluster.index(id)
  ncluster <- length(cc$clusters)
  U <- function(theta,indiv=FALSE) {
    B <- .Call("Bhat",status,X,theta,cc$clusters,cc$idclust,cc$cluster.size)$B2
    Ba <- B[fB$pos+1,,drop=FALSE]
    Hij <- as.vector(X*Ba)
    res <- .Call("Uhat",status,Hij,theta,cc$idclust,cc$cluster.size)
    if (!indiv) res <- mean(res)
    res
  }
  if (missing(theta)) theta <- 0.1
  if (is.list(theta)) return(lapply(theta,function(x) U(x,...)))

  op <- nlminb(theta,function(x) U(x)^2)
  uu <- U(op$par,TRUE)
  du <- grad(U,op$par)
  return(list(theta=op$par, sd=(mean(uu^2)/du^2/ncluster)^0.5))
}
                  
  
