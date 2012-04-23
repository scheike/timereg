\name{summary.cor}
\alias{summary.cor}
\title{Summary for dependence models for competing risks} 
\description{ Computes concordance and casewise concordance for dependence models for competing risks 
models of the type cor.cif, rr.cif or or.cif for the given cumulative incidences and the different dependence
measures in the object..
}
\usage{
\method{summary}{cor}(object,marg.cif=NULL,marg.cif2=NULL,digits=3,...)
}
\arguments{
\item{object}{object from cor.cif rr.cif or or.cif for dependence between competing risks data for two causes.}
\item{marg.cif}{a number that gives the cumulative incidence in one time point for which concordance and 
	casewise concordance are computed.}
\item{marg.cif2}{the cumulative incidence for cause 2 for concordance and 
	casewise concordance are computed. Default is that it is the same as marg.cif.}
\item{digits}{digits in output.}
\item{...}{ekstra arguments.}
}
\value{prints summary for dependence model. 
\item{casewise}{gives casewise concordance that is, probability of cause 2 (related to cif2) given that cause 1 (related to cif1)
	has occured.}
\item{concordance}{gives concordance that is, probability of cause 2 (related to cif2) and cause 1 (related to cif1).}
\item{cif1}{cumulative incidence for cause1.}
\item{cif2}{cumulative incidence for cause1.}
}
\references{
Cross odds ratio Modelling of dependence for
Multivariate Competing Risks Data, Scheike and Sun (2010), work in progress.

A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
Scheike, Zhang, Sun, Jensen (2010), Biometrika. 
}
\author{Thomas Scheike}
\examples{
data(multcif) # simulated data 
multcif$cause[multcif$cause==0] <- 2

times=seq(0.05,3,by=0.1) # to speed up computations use only these time-points
add<-comp.risk(Surv(time,status>0)~const(X)+cluster(id),data=multcif,
               multcif$cause,n.sim=0,times=times)
###
out1<-cor.cif(add,data=multcif,cause1=1,cause2=1,theta=log(2+1))
summary(out1)

pad <- predict(add,X=1,Z=0,se=0,uniform=0)$P1
summary(out1,marg.cif=pad)

}
\keyword{survival}
