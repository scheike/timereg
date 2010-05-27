\name{cox.aalen}
\alias{cox.aalen}
\title{Fit Cox-Aalen survival model}
\description{
Fits an Cox-Aalen survival model. 
Time dependent variables and counting process data
(multiple events per subject) are possible.

\deqn{ \lambda_{i}(t) = Y_i(t) 
( X_{i}^T(t) \alpha(t) )  \exp(Z_{i}^T \beta ) 
}

The model thus contains the Cox's regression  model as special case. 

Resampling is used for computing p-values for tests of
time-varying effects.
Test for proportionality is considered by
considering the score processes for the proportional effects of model.

The modelling formula uses the standard survival modelling given in the
\bold{survival} package.
}
\usage{
cox.aalen(formula=formula(data),data=sys.parent(),beta=0,Nit=10,detail=0,
start.time=0,max.time=NULL,id=NULL,clusters=NULL,n.sim=500,residuals=0,
robust=1,weighted.test=0,covariance=0,resample.iid=0,weights=NULL,rate.sim=1,
beta.fixed=0)
}
\arguments{
\item{formula}{
  a formula object with the response on the left of a '~' operator, and
  the independent terms on the right as regressors. The response must be
  a survival object as returned by the `Surv' function. Terms with a
  proportional effect are specified by the wrapper prop(), and cluster
  variables (for computing robust variances) by the wrapper cluster().
}
\item{data}{
a data.frame with the variables.
}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{end of observation period where estimates are computed. 
Estimates thus computed from [start.time, max.time]. Default is max of data.}
\item{robust}{to compute robust variances and construct processes for
resampling. May be set to 0 to save memory.}
\item{id}{For timevarying covariates the variable must 
associate each record with the id of a subject.}
\item{clusters}{cluster variable for computation of robust 
variances.}
\item{n.sim}{number of simulations in resampling.}
\item{weighted.test}{to compute a variance weighted version of the 
test-processes used for testing
time-varying effects.}
\item{residuals}{to returns residuals that can be used for
model validation in the function cum.residuals. Estimated martingale increments (dM)
and corresponding time vector (time).
}
\item{covariance}{to compute covariance estimates for
nonparametric terms rather than just the variances.}
\item{resample.iid}{to return i.i.d. representation for
nonparametric and parametric terms.}
\item{beta}{starting value for relative risk estimates}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details is printed during iterations, if 1 details are given.}
\item{weights}{weights for weighted analysis.}
\item{rate.sim}{rate.sim=1 such that resampling of residuals is based on estimated
martingales and thus valid in rate case, rate.sim=0 means that resampling is based 
on counting processes and thus only valid in intensity case.}
\item{beta.fixed}{option for computing score process for fixed relative risk parameter}
}
\details{
The data for a subject is presented as multiple rows or 'observations',
each of which applies to an interval of observation (start, stop]. 
For counting process data with the )start,stop] notation is used
the 'id' variable is needed to identify the records for each subject.
The program assumes that there are no ties, and if such are present 
random noise is added to break the ties.
}
\value{
returns an object of type "cox.aalen". With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.
}
\item{var.cum}{the martingale based pointwise variance estimates.  }
\item{robvar.cum}{robust pointwise variances estimates.  }
\item{gamma}{estimate of parametric components of model.  }
\item{var.gamma}{variance for gamma.  }
\item{robvar.gamma}{robust variance for gamma.  }
\item{residuals}{list with residuals.}
\item{obs.testBeq0}{observed absolute value of supremum of 
cumulative components scaled with the variance.}
\item{pval.testBeq0}{p-value for covariate effects based on supremum test.}
\item{sim.testBeq0}{resampled supremum values.}
\item{obs.testBeqC}{observed absolute value of supremum of 
difference between observed cumulative
process and estimate under null of constant effect.}
\item{pval.testBeqC}{p-value based on resampling.}
\item{sim.testBeqC}{resampled supremum values.}
\item{obs.testBeqC.is}{observed integrated squared differences between 
observed cumulative and estimate under null of constant effect.}
\item{pval.testBeqC.is}{p-value based on resampling.}
\item{sim.testBeqC.is}{resampled supremum values.}
\item{conf.band}{resampling based constant to construct robust 95\% uniform 
confidence bands. }
\item{test.procBeqC}{observed test-process of 
difference between observed cumulative 
process and estimate under null of constant effect over time. 
}
\item{sim.test.procBeqC}{list of 50 random realizations of
test-processes under null based on resampling.}
\item{covariance}{covariances for nonparametric terms of model.}
\item{B.iid}{Resample processes for nonparametric terms of model.}
\item{gamma.iid}{Resample processes for parametric terms of model.}
\item{loglike}{approximate log-likelihood for model, similar to
Cox's partial likelihood. Only computed when robust=1.}
\item{D2linv}{inverse of the derivative of the score function.}
\item{score}{value of score for final estimates.}
\item{test.procProp}{observed score process for proportional part of model.}
\item{var.score}{variance of score process (optional variation estimator
for beta.fixed=1 and robust estimator otherwise).}
\item{pval.Prop}{p-value based on resampling.}
\item{sim.supProp}{re-sampled absolute supremum values.}
\item{sim.test.procProp}{list of 50 random realizations of 
test-processes for proportionality under the model based on resampling.}
}
\references{
Martinussen and Scheike, Dynamic Regression Models for Survival Data, 
Springer (2006).  }
\author{Thomas Scheike}
\examples{
library(survival)
data(sTRACE)
# Fits Cox model 
out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
prop(vf)+prop(chf)+prop(diabetes),sTRACE,max.time=7,n.sim=100)

# makes Lin, Wei, Ying test for proportionality
summary(out)
par(mfrow=c(2,3))
plot(out,score=1) 

# Fits Cox-Aalen model 
out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
vf+chf+prop(diabetes),sTRACE,max.time=7,n.sim=100)

summary(out)
par(mfrow=c(2,3))
plot(out)
}
\keyword{survival}
