\name{concordance}
\alias{concordance}
\title{Computes concordance and casewise concordance for dependence models for competing risks} 
\description{
Computes concordance and casewise concordance for dependence models for competing risks 
models of the type cor.cif, rr.cif or or.cif.
}
\usage{
concordance(object,cif1,cif2=NULL,messages=TRUE,model=NULL,coefs=NULL,...)
}
\arguments{
\item{object}{object from cor.cif rr.cif or or.cif for dependence between competing risks data for two causes.}
\item{cif1}{vector of cumulative incidences for different time-points for which concordance and 
	casewise concordance are computed for cause1.}
\item{cif2}{vector of cumulative incidences for different time-points for which concordance and 
	casewise concordance are computed for cause2. Default is the same as cause 1. Must be computed in the same time-points as cif1.}
\item{messages}{messages?.}
\item{model}{specifies which model to use }
\item{coefs}{coefficientsf for something }
\item{...}{ekstra arguments.}
}
\value{returns concordance and casewise concordance for dependence models for competing risks. List for each dependence
	measure in object, that contains
\item{concordance}{gives concordance that is, probability of cause 2 (related to cif2) and cause 1 (related to cif1).}
\item{casewise.concordance}{gives casewise concordance that is, probability 
	of cause 2 (related to cif2) given that cause 1 (related to cif1) has occured.}
\item{cif1}{cumulative incidence for cause1.}
\item{cif2}{cumulative incidence for cause1.}
}
\references{
Estimating heritability for cause specific mortality based on twin data,
Scheike, Holst and Hjelmborg (2011), submitted.

Cross odds ratio Modelling of dependence for Multivariate Competing Risks Data, 
Scheike and Sun (2010), submitted.

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

out1<-cor.cif(add,data=multcif,cause1=1,cause2=1,theta=log(2+1))
summary(out1)

pad <- predict(add,X=1,Z=0,se=0,uniform=0)$P1[,2] 
summary(out1,marg.cif=pad)
}
\keyword{survival}
