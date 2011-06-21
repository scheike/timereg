\name{rr.cif}
\alias{rr.cif}
\title{Ratio regression for competing risks} 
\description{
Fits a parametric model for the log-ratio 
of the joint probabilty of two events compared to the similar 
probability under independence. 
We consider the 
\deqn{
rr_{1,i;2,k}(t) = h(\theta,t) = 
\log ( \frac{P(T_1 \leq t,cause_1=i , T_2 \leq t,cause_2=k)}{
P(T_1 \leq t,cause_1=i) P(T_2 \leq t,cause_2=k)} ) 
}
the joint probabilty divided by the product.

\eqn{h()} may be specified as an R-function of the parameters, 
see example below, but the default is that it is simply \eqn{\theta^T z}.
}
\usage{
rr.cif(cif,data=sys.parent(),cause,times=NULL,
cause1=1,cause2=1,cens.code=0,cif2=NULL,cens.model="KM",Nit=40,detail=0,
clusters=NULL,theta=NULL,theta.des=NULL,parfunc=NULL,dparfunc=NULL,
step=1,sym=0,colnames=NULL,dimpar=NULL,weights=NULL,
notaylor=0,same.cens=FALSE,stab.cens=FALSE,entry=NULL,estimator=1,trunkp=1,
admin.cens=NULL)
}
\arguments{
\item{cif}{a model object from the comp.risk function with the 
marginal cumulative incidence of cause2, i.e., the event that is conditioned on, and whose
odds the comparision is made with respect to}
\item{data}{a data.frame with the variables.}
\item{cause}{specifies the causes  related to the death
	times, the value cens.code is the censoring value.}
\item{times}{time-vector that specifies the times used for the estimating euqations for the cross-odds-ratio estimation.}
\item{cause1}{specificies the cause considered.}
\item{cause2}{specificies the cause that is conditioned on.}
\item{cif2}{marginal model for cause 2, must be specified when cause1 and cause2 differs}
\item{cens.code}{specificies the code for the censoring.}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{clusters}{specifies the cluster structure.}
\item{theta}{specifies starting values for the cross-odds-ratio parameters of the model.}
\item{theta.des}{specifies a regression design for the cross-odds-ratio parameters.}
\item{parfunc}{R-function to specify any relationship for the parameters of the cross-odds-ratio.}
\item{dparfunc}{R-function to specify the derivative of the parfunc function with respect to
the parameters.}
\item{step}{specifies the step size for the Newton-Raphson algorithm.}
\item{sym}{specifies if symmetry is used in the model, for sym=0 uses only pairs i<j for concordance.}
\item{colnames}{specifies possible colnames to go with the parameters.}
\item{dimpar}{number of parameters in case where functions specify the model.}
\item{weights}{weights for estimating equations.}
\item{notaylor}{if TRUE then the standard errors ignores the uncertainty from the marginal
cumulative incidence. This will speed things up considerably.} 
\item{same.cens}{if true then censoring within clusters are assumed to be the same variable, default is independent censoring.}
\item{stab.cens}{when the censoring weights are 0 in some places, it may
be useful to use modified censoring weights Gc(min(T,s)) in the score 
equations, and to restrict attention to those s where Gc() is positive.}
\item{entry}{entry-age in case of delayed entry. Then two causes must be given.} 
\item{estimator}{gives a design to fit models for each strata (cause combinations).} 
\item{trunkp}{gives probability of survival for delayed entry, and related to entry-ages given above.}
\item{admin.cens}{censoring times for administrative censoring.}
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
Cross odds ratio Modelling of dependence for
Multivariate Competing Risks Data, Scheike and Sun (2010), work in progress.

A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
Scheike, Zhang, Sun, Jensen (2010), Biometrika. 
}
\author{Thomas Scheike}
\examples{
library(MultiComp)
data(multcif) # simulated data with cross-odds-ratio 3
multcif$cause[multcif$cause==0] <- 2

times12<-sort(multcif$time[multcif$status \%in\% c(1,2)])[-(1:10)]

add1<-comp.risk(Surv(time,status>0)~const(X),data=multcif,times=times12,
                causeS=1,multcif$cause,clusters=multcif$id,n.sim=0)

add2<-comp.risk(Surv(time,status>0)~const(X),data=multcif,times=times12,
                causeS=2,multcif$cause,clusters=multcif$id,n.sim=0)

out12<-rr.cif(add1,data=multcif,multcif$cause,cif2=add2,cens.code=0,cause1=1,cause2=2,
               cens.model="KM",Nit=10,clusters=multcif$id)
summary(out12)
}
\keyword{survival}
