\name{random.cif2cause}
\alias{random.cif2cause}
\title{Random effects model for competing risks data} 
\description{
Fits a random effects  model describing the dependence in the cumulative 
incidence curves for subjects within a cluster.  Given the gamma distributed
random effects it is assumed that the cumulative incidence curves are indpendent, and
that the marginal cumulative incidence curves are on additive form
for two different causes: 
\deqn{ P(T_1 \leq t, cause=i | x,z) = P_1(t,x,z) = 1- exp( -x_i^T A_i(t) + - (t z_i)^T \beta_i) }
and 
\deqn{ P(T_2 \leq t, cause=k | x,z) = P_1(t,x,z) = 1- exp( -x_k^T A_k(t) + - (t z_k)^T \beta_k) }
These models are specified by fitting two separate marginals models in the same
time-points (using the times argument of the comp.risk function).

We allow a regression structure for the random effects variances that may depend on
cluster covariates.
}
\usage{
random.cif2cause(cif,data=sys.parent(),cause,
cause1=1,cause2=1,cens.code=0,cif2=NULL, 
cens.model="KM",Nit=40,detail=0,clusters=NULL,theta=NULL,
theta.des=NULL,step=1,squarepar=0,same.cens=FALSE)
}
\arguments{
\item{cif}{a model object from the comp.risk function with the 
marginal cumulative incidence of cause1.}
\item{data}{a data.frame with the variables.}
\item{cause}{specifies the causes  related to the death
	times, the value cens.code is the censoring value.}
\item{cause1}{specificies one cause of interest.}
\item{cause2}{specificies one cause of interest.}
\item{cens.code}{specificies the code for the censoring.}
\item{cif2}{a model object from the comp.risk function with the 
marginal cumulative incidence of cause2.}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{clusters}{specifies the cluster structure.}
\item{theta}{specifies starting values for the cross-odds-ratio parameters of the model.}
\item{theta.des}{specifies a regression design for the cross-odds-ratio parameters.}
\item{step}{specifies the step size for the Newton-Raphson algorithm.}
\item{squarepar}{uses squared parameters to make sure that the random effect becomes positive.}
\item{same.cens}{if true then censoring within clusters are assumed to be the same variable,
default is independent censoring.}
}
\value{returns an object of type 'cor'. With the following arguments:
\item{theta}{estimate of proportional odds parameters of model.}
\item{var.theta}{variance for gamma.  }
\item{hess}{the derivative of the used score.}
\item{score}{scores at final stage.}
}
\references{
A Semiparametric Random Effects Model for Multivariate Competing Risks Data,
Scheike, Zhang, Sun, Jensen (2009), under review. 

Cross odds ratio Modelling of dependence for
Multivariate Competing Risks Data, Scheike and Sun (2009), work in progress.
}
\author{Thomas Scheike}
\examples{
data(multcif)
multcif$cause[multcif$cause==0] <- 2

times<-sort(multcif$time[multcif$status \%in\% c(1,2)])

add1<-comp.risk(Surv(time,status>0)~const(X)+cluster(id),data=multcif,causeS=1,
               multcif$cause,n.sim=0,times=times)

add2<-comp.risk(Surv(time,status>0)~const(X)+cluster(id),data=multcif,causeS=2,
               multcif$cause,n.sim=0,times=times)

out1<-random.cif2cause(add1,data=multcif,cause1=1,cause2=2,cif2=add2,
                       cens.model="KM",theta=1,theta.des=NULL,step=1.0)
summary(out1)
}
\keyword{survival}
