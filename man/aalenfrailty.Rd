\name{aalenfrailty}
\alias{aalenfrailty}
\title{Aalen frailty model}
\usage{
  aalenfrailty(time, status, X, id, theta, B = NULL, ...)
}
\arguments{
  \item{time}{Time variable}

  \item{status}{Status variable (0,1)}

  \item{X}{Covariate design matrix}

  \item{id}{cluster variable}

  \item{theta}{list of thetas (returns score evaluated
  here), or starting point for optimization (defaults to
  magic number 0.1)}

  \item{B}{(optional) Cumulative coefficients (update theta
  by fixing B)}

  \item{...}{Additional arguments to lower level functions}
}
\value{
  Parameter estimates
}
\description{
  Additive hazards model with (gamma) frailty
}
\details{
  Aalen frailty model
}
\examples{
dd <- simAalenFrailty(5000)
f <- ~1##+x
X <- model.matrix(f,dd) ## design matrix for non-parametric terms
system.time(out<-aalen(update(f,Surv(time,status)~.),dd,n.sim=0,robust=0))
dix <- which(dd$status==1)
t1 <- system.time(bb <- .Call("Bhat",as.integer(dd$status),
                              X,0.2,as.integer(dd$id),NULL,NULL,-1,
                              package="mets"))
spec <- 1
##plot(out,spec=spec)
## plot(dd$time[dix],bb$B2[,spec],col="red",type="s",
##      ylim=c(0,max(dd$time)*c(beta0,beta)[spec]))
## abline(a=0,b=c(beta0,beta)[spec])
##'

\dontrun{thetas <- seq(0.1,2,length.out=10)
Us <- unlist(aalenfrailty(dd$time,dd$status,X,dd$id,as.list(thetas)))
plot(thetas,Us,type="l",ylim=c(-.5,1)); abline(h=0,lty=2); abline(v=theta,lty=2)
op <- aalenfrailty(dd$time,dd$status,X,dd$id)
op}
}
\author{
  Klaus K. Holst
}

