\name{two.stage}
\alias{two.stage}
\title{Fit Clayton-Oakes-Glidden Two-Stage model}
\description{
Fit Clayton-Oakes-Glidden Two-Stage model with 
Cox-Aalen marginals and regression on the variance parameters.

The model specifikatin allows a regression structure on 
the variance of the random effects, such it is allowed to depend on 
covariates fixed within clusters 
\deqn{ \theta_{k} =  Q_{k}^T \nu }. This is particularly useful to
model jointly different groups and  to compare their variances. 

Fits an Cox-Aalen survival model. 
Time dependent variables and counting process data
(multiple events per subject) are not possible !

The marginal baselines are on the Cox-Aalen form
\deqn{ \lambda_{ki}(t) = Y_{ki}(t) 
( X_{ki}^T(t) \alpha(t) )  \exp(Z_{ki}^T \beta ) 
}

The model thus contains the Cox's regression model and the 
additive hazards model as special cases. (see cox.aalen function
for more on this).

The modelling formula uses the standard survival modelling given in the
\bold{survival} package. Only for right censored survival data. 
}
\usage{
two.stage(margsurv,data=sys.parent(),Nit=60, detail=0,
start.time=0,max.time=NULL,id=NULL,clusters=NULL,robust=1,
theta=NULL,theta.des=NULL,var.link=0,step=0.5,notaylor=0,se.clusters=NULL)
}
\arguments{
\item{margsurv}{fit of marginal survival cox.aalen model with residuals=2, and
resample.iid=1 to get fully correct standard errors. See notaylor below.} 
\item{data}{ a data.frame with the variables.}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{end of observation period where estimates are computed. 
Estimates thus computed from [start.time, max.time]. Default is max of data.}
\item{id}{For timevarying covariates the variable must 
associate each record with the id of a subject.}
\item{clusters}{cluster variable for computation of robust variances.}
\item{robust}{if 0 then totally omits computation of standard errors.}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details is printed during iterations, if 1 details are given.}
\item{theta}{starting values for the frailty variance (default=0.1).}
\item{theta.des}{design for regression for variances. The defauls is NULL that 
is equivalent to just one theta and the design with only a baseline.}
\item{var.link}{default "0" is that the regression design on the 
variances is without a link, and "1" uses the link function exp. }
\item{step}{step size for Newton-Raphson.}
\item{notaylor}{if 1 then ignores variation due to survival model, this is quicker
and then resample.iid=0 and residuals=0 is ok for marginal survival model that then is much quicker.}
\item{se.clusters}{cluster variable for sandwich estimator of variance.}
}
\details{
The data for a subject is presented as multiple rows or 'observations',
each of which applies to an interval of observation (start, stop]. 
For counting process data with the )start,stop] notation is used
the 'id' variable is needed to identify the records for each subject.
Only one record per subject is allowed in the current
implementation for the estimation of theta. 
The program assumes that there are no ties, and if such are present 
random noise is added to break the ties.

Left truncation is dealt with. Here the key assumption is that the maginals are
correctly estimated and that we have a common truncation time within each cluster. 
}
\value{
returns an object of type "two.stage". With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.}
\item{var.cum}{the martingale based pointwise variance estimates.}
\item{robvar.cum}{robust pointwise variances estimates.}
\item{gamma}{estimate of parametric components of model.}
\item{var.gamma}{variance for gamma.}
\item{robvar.gamma}{robust variance for gamma.}
\item{D2linv}{inverse of the derivative of the score function from marginal model.}
\item{score}{value of score for final estimates.}
\item{theta}{estimate of Gamma variance for frailty.}
\item{var.theta}{estimate of variance of theta.}
\item{SthetaInv}{inverse of derivative of score of theta.}
\item{theta.score}{score for theta parameters.}
}
\references{
Glidden (2000), A Two-Stage estimator of the dependence parameter
for the Clayton Oakes model.

Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer (2006).}
\author{Thomas Scheike}
\examples{
library(timereg)
data(diabetes)
# Marginal Cox model  with treat as covariate
marg <- cox.aalen(Surv(time,status)~prop(treat)+prop(adult)+
	  cluster(id),data=diabetes,resample.iid=1)
fit<-two.stage(marg,data=diabetes,theta=1.0,Nit=40)
summary(fit)

# using coxph and giving clusters, but SE wittout cox uncetainty
margph <- coxph(Surv(time,status)~treat,data=diabetes)
fit<-two.stage(margph,data=diabetes,theta=1.0,Nit=40,clusters=diabetes$id)


# Stratification after adult 
theta.des<-model.matrix(~-1+factor(adult),diabetes);
des.t<-model.matrix(~-1+factor(treat),diabetes);
design.treat<-cbind(des.t[,-1]*(diabetes$adult==1),
                    des.t[,-1]*(diabetes$adult==2))

# test for common baselines included here 
marg1<-cox.aalen(Surv(time,status)~-1+factor(adult)+prop(design.treat)+cluster(id),
 data=diabetes,resample.iid=1,Nit=50)

fit.s<-two.stage(marg1,data=diabetes,Nit=40,theta=1,theta.des=theta.des)
summary(fit.s)

# with common baselines  and common treatment effect (although test reject this)
fit.s2<-two.stage(marg,data=diabetes,Nit=40,theta=1,theta.des=theta.des)
summary(fit.s2)

# test for same variance among the two strata
theta.des<-model.matrix(~factor(adult),diabetes);
fit.s3<-two.stage(marg,data=diabetes,Nit=40,theta=1,theta.des=theta.des)
summary(fit.s3)

# to fit model without covariates, use beta.fixed=1 and prop or aalen function
marg <- aalen(Surv(time,status)~+1+cluster(id),
	 data=diabetes,resample.iid=1,n.sim=0)
fita<-two.stage(marg,data=diabetes,theta=0.95,detail=0)
summary(fita)

# same model but se's without variation from marginal model to speed up computations
marg <- aalen(Surv(time,status) ~+1+cluster(id),data=diabetes,
	      resample.iid=0,n.sim=0)
fit<-two.stage(marg,data=diabetes,theta=0.95,detail=0)
summary(fit)

# same model but se's now with fewer time-points for approx of iid decomp of marginal 
# model to speed up computations
marg <- cox.aalen(Surv(time,status) ~+prop(treat)+cluster(id),data=diabetes,
	      resample.iid=1,n.sim=0,max.timepoint.sim=5,beta.fixed=1,beta=0)
fit<-two.stage(marg,data=diabetes,theta=0.95,detail=0)
summary(fit)
}
\keyword{survival}
