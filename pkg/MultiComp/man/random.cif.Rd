\name{random.cif}
\alias{random.cif}
\title{Random effects model for competing risks data} 
\description{
Fits a random effects  model describing the dependence in the cumulative 
incidence curves for subjects within a cluster.  Given the gamma distributed
random effects it is assumed that the cumulative incidence curves are indpendent, and
that the marginal cumulative incidence curves are on additive form
\deqn{
P(T \leq t, cause=1 | x,z) = P_1(t,x,z) = 1- exp( -x^T A(t) - t z^T \beta)
}

We allow a regression structure for the random effects variances that may depend on
cluster covariates.

}
\usage{
random.cif(cif,data=sys.parent(),cause,causeS=1,cens.code=0,
cens.model="KM",Nit=40,detail=0,clusters=NULL,theta=NULL,
theta.des=NULL,step=1,same.cens=FALSE,entry=NULL,trunkp=1,notaylor=0)
}
\arguments{
\item{cif}{a model object from the comp.risk function with the 
marginal cumulative incidence of cause2, i.e., the event that is conditioned on, and whose
odds the comparision is made with respect to}
\item{data}{a data.frame with the variables.}
\item{cause}{specifies the causes  related to the death
	times, the value cens.code is the censoring value.}
\item{causeS}{specificies the cause of interst.}
\item{cens.code}{specificies the code for the censoring.}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{clusters}{specifies the cluster structure.}
\item{theta}{specifies starting values for the cross-odds-ratio parameters of the model.}
\item{theta.des}{specifies a regression design for the cross-odds-ratio parameters.}
\item{step}{specifies the step size for the Newton-Raphson algorithm.}
\item{same.cens}{if true then censoring within clusters are assumed to be the same variable, default is independent censoring.}
\item{entry}{entry-age in case of delayed entry. Then two causes must be given.} 
\item{trunkp}{gives probability of survival for delayed entry, and related to entry-ages given above.}
\item{notaylor}{if TRUE then the standard errors ignores the uncertainty 
from the marginal cumulative incidence. This will speed things up considerably.} 
}
\value{returns an object of type 'cor'. With the following arguments:
\item{theta}{estimate of proportional odds parameters of model.}
\item{var.theta}{variance for gamma.  }
\item{hess}{the derivative of the used score.}
\item{score}{scores at final stage.}
\item{score}{scores at final stage.}
\item{theta.iid}{matrix of iid decomposition of parametric effects.}
}
\references{
A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
Scheike, Zhang, Sun, Jensen (2010), Biometrika. 

Cross odds ratio Modelling of dependence for
Multivariate Competing Risks Data, Scheike and Sun (2012), work in progress.
}
\author{Thomas Scheike}
\examples{
data(np)

add<-comp.risk(Surv(time,cause>0)~const(country)+cluster(id),data=np,np$cause,n.sim=0)

out1<-random.cif(add,data=np,causeS=1,cens.model="KM",Nit=15,detail=0,theta=1,theta.des=NULL,step=1.0)
summary(out1)

theta.des<-model.matrix(~-1+factor(zyg),data=np)

out2<-random.cif(add,data=np,causeS=1,cens.model="KM",Nit=15,detail=0,theta=1,theta.des=theta.des,step=1.0)
summary(out2)
}
\keyword{survival}
