\name{pe.sasieni}
\alias{pe.sasieni}
\alias{summary.pe-sasieni}
\title{Fits Proportional excess hazards model with fixed offsets}
\description{
Fits proportional excess hazards model. The Sasieni proportional
excess risk model. 

The models are written using the survival modelling given in the 
survival package.
}
\usage{
pe.sasieni(formula=formula(data),data=sys.parent(),
id=NULL,start.time=0,max.time=NULL,offsets=0,Nit=50,detail=0,n.sim=500)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a `~' operator, and
the terms on the right.  The response must be a survival object as
returned by the `Surv' function.}
\item{data}{a data.frame with the variables.}
\item{id}{gives the number of individuals.}
\item{start.time}{starting time for considered time-period.}
\item{max.time}{stopping considered time-period if different from 0. 
Estimates thus computed from [0,max.time] if max.time>0. Default is max of data.}
\item{offsets}{fixed offsets giving the mortality.}
\item{Nit}{number of itterations.}
\item{detail}{if detail is one, prints iteration details.}
\item{n.sim}{number of simulations, 0 for no simulations.}
}
\details{
The program assumes that there are no ties, and if such are present 
random noise is added to break the ties.
}
\value{Returns an object of type "pe.sasieni". 
With the following arguments:
\item{cum}{baseline of Cox model excess risk.}
\item{var.cum}{pointwise variance estimates for estimated cumulatives.}
\item{gamma}{estimate of relative risk terms of model.}
\item{var.gamma}{variance estimates for gamma.}
\item{Ut}{score process for Cox part of model.}
\item{D2linv}{The inverse of the second derivative.}
\item{score}{final score}
\item{test.Prop}{re-sampled absolute supremum values.}
\item{pval.Prop}{p-value based on resampling.}
}
\references{
Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer Verlag (2006).

Sasieni, P.D., Proportional excess hazards, Biometrika (1996),
127--41. 

Cortese, G. and Scheike, T.H., Dynamic regression hazards models
for relative survival (2007), submitted.  
}
\author{Thomas Scheike}
\examples{
data(mela.pop)
out<-pe.sasieni(Surv(start,stop,status==1)~age+sex,mela.pop,
id=1:205,Nit=10,max.time=7,offsets=mela.pop$rate,detail=0,n.sim=100)
summary(out)

ul<-out$cum[,2]+1.96*out$var.cum[,2]^.5
ll<-out$cum[,2]-1.96*out$var.cum[,2]^.5
plot(out$cum,type="s",ylim=range(ul,ll))
lines(out$cum[,1],ul,type="s"); lines(out$cum[,1],ll,type="s")
# see also prop.excess function
}
\keyword{survival}
