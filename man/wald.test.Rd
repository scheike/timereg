\name{wald.test}
\alias{wald.test}
\title{Makes wald test}
\description{Makes wald test on timereg object, either by contrast matrix or
testing components to 0. Can also specify the regression coefficients and 
the variance matrix. 
}
\usage{
wald.test(object=NULL,coef=NULL,Sigma=NULL,contrast,coef.null=NULL,null=NULL)
}
\arguments{
  \item{object}{timereg object}
  \item{coef}{estimates from some model}
  \item{Sigma}{variance of  estimates}
  \item{contrast}{contrast matrix for testing}
  \item{coef.null}{which indeces to test to 0}
  \item{null}{mean of null, 0 by default}
}
\author{Thomas Scheike}
\examples{
data(sTRACE)
# Fits Cox model 
out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
prop(vf)+prop(chf)+prop(diabetes),data=sTRACE,n.sim=0)

wald.test(out,coef.null=c(1,2,3))
### test age=sex   vf=chf
wald.test(out,contrast=rbind(c(1,-1,0,0,0),c(0,0,1,-1,0)))

### now same with direct specifation of estimates and variance
wald.test(coef=out$gamma,Sigma=out$var.gamma,coef.null=c(1,2,3))
wald.test(coef=out$gamma,Sigma=out$robvar.gamma,coef.null=c(1,2,3))
### test age=sex   vf=chf
wald.test(coef=out$gamma,Sigma=out$var.gamma,
	  contrast=rbind(c(1,-1,0,0,0),c(0,0,1,-1,0)))
}
\keyword{survival}
