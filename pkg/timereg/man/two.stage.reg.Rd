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
\bold{survival} package.
}
\usage{
two.stage(formula=formula(data),data=sys.parent(),beta=0,Nit=60,
detail=0,start.time=0,max.time=NULL,id=NULL,clusters=NULL,
robust=1,rate.sim=1,beta.fixed=0,theta=NULL,theta.des=NULL,var.link=0,step=1)
}
\arguments{
\item{formula}{
  a formula object with the response on the left of a '~' operator, and
  the independent terms on the right as regressors. The response must be
  a survival object as returned by the `Surv' function. Terms with a
  proportional effect are specified by the wrapper prop(), and cluster
  variables (for computing robust variances) by the wrapper cluster().
}
\item{data}{ a data.frame with the variables.}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{end of observation period where estimates are computed. 
Estimates thus computed from [start.time, max.time]. Default is max of data.}
\item{robust}{to compute robust variances and construct processes for
resampling. May be set to 0 to save memory.}
\item{id}{For timevarying covariates the variable must 
associate each record with the id of a subject.}
\item{clusters}{cluster variable for computation of robust variances.}
\item{beta}{starting value for relative risk estimates}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details is printed during iterations, if 1 details are given.}
\item{rate.sim}{rate.sim=1 such that resampling of residuals is based on estimated martingales and thus valid in rate case, rate.sim=0 means that resampling is based on counting processes and thus only valid in intensity case.}
\item{beta.fixed}{option for keeping beta in the Cox-Aalen model fixed.}
\item{theta}{starting values for the frailty variance (default=0.1).}
\item{theta.des}{design for regression for variances. The defauls is NULL that 
is equivalent to just one theta and the design with only a baseline.}
\item{var.link}{default "0" is that the regression design on the variances 
is without a link, and "1" uses the link function exp. }
\item{step}{step size for Newton-Raphson.}
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
}
\value{
returns an object of type "two.stage". With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.}
\item{var.cum}{the martingale based pointwise variance estimates.  }
\item{robvar.cum}{robust pointwise variances estimates.}
\item{gamma}{estimate of parametric components of model.}
\item{var.gamma}{variance for gamma.}
\item{robvar.gamma}{robust variance for gamma.}
\item{D2linv}{inverse of the derivative of the score function.}
\item{score}{value of score for final estimates.}
\item{theta}{estimate of Gamma variance for frailty.}
\item{var.theta}{estimate of variance of theta.}
\item{S.theta}{estimate of derivative of score of theta.}
\item{theta.score}{score for theta parameters.}
}
\references{
Glidden (2000), A Two-Stage estimator of the dependence parameter
for the Clayton Oakes model.

Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer (2006).}
\author{Thomas Scheike}
\examples{
library(survival)
data(diabetes)
# Marginal Cox model  with treat as covariate
fit<-two.stage(Surv(time,status) ~ prop(treat) + cluster(id),
data=diabetes,Nit=40,theta=1)
summary(fit)

# Stratification after adult 
theta.des<-model.matrix(~-1+factor(adult),diabetes);
des.t<-model.matrix(~-1+factor(treat),diabetes);
design.treat<-cbind(des.t[,-1]*(diabetes$adult==1),
                    des.t[,-1]*(diabetes$adult==2))
fit.s<-two.stage(Surv(time,status) ~
-1+factor(adult)+prop(design.treat)+cluster(id),
data=diabetes,Nit=40,theta=1,theta.des=theta.des)
summary(fit.s)

# test for common baselines 
fit.s1<-cox.aalen(Surv(time,status) ~ 
factor(adult)+prop(design.treat)+cluster(id),data=diabetes)
summary(fit.s1)

# with common baselines  and common treatment effect (although test reject this)
fit.s2<-two.stage(Surv(time,status) ~+1+prop(treat) + cluster(id),
data=diabetes,Nit=40,theta=1,theta.des=theta.des)
summary(fit.s2)

# test for same variance among the two strata
theta.des<-model.matrix(~factor(adult),diabetes);
fit.s3<-two.stage(Surv(time,status) ~+1+prop(treat)+cluster(id),
data=diabetes,Nit=40,theta=1,theta.des=theta.des)
summary(fit.s3)

# to fit model without covariates, beta.fixed=1, but still need prop term !
fit<-two.stage(Surv(time,status) ~ prop(treat) + cluster(id),
data=diabetes,theta=0.95,detail=0,beta.fixed=1)
summary(fit)
}
\keyword{survival}
