\name{summary.cum.residuals}
\alias{summary.cum.residuals}
\title{Prints summary statistics for goodness-of-fit tests based on cumulative residuals}
\description{
Computes p-values for extreme behaviour relative to the model of 
various cumulative residual processes.}
\usage{\method{summary}{cum.residuals}(object,digits=3,...) }
\arguments{
\item{object}{output from the cum.residuals() function.}
\item{digits}{number of digits in printouts.}
\item{...}{unused arguments - for S3 compatibility}
}
\author{Thomas Scheike}
\examples{
library(survival)
data(sTRACE)
# Fits Aalen model and returns residuals
out<-aalen(Surv(time,status==9)~age+sex+diabetes+chf+vf,
sTRACE,max.time=7,n.sim=0,residuals=1)

# constructs and simulates cumulative residuals versus age groups
# and versus covariates of model
out.mg<-cum.residuals(out,sTRACE,
modelmatrix=model.matrix(~-1+factor(cut(age,4)),sTRACE),cum.resid=1,n.sim=100)

summary(out.mg)
}
\keyword{survival}
