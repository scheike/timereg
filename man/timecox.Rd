\name{timecox}
\alias{timecox}
\title{Fit Cox model with partly timevarying effects.
}
\description{
Fits proportional hazards model with some effects time-varying and 
some effects constant. 
Time dependent variables and counting process data
(multiple events per subject) are possible.

Resampling is used for computing p-values for tests of
timevarying effects.

The modelling formula uses the standard survival modelling given in the
\bold{survival} package.
}
\usage{
timecox(formula=formula(data),data=sys.parent(),
start.time=0,max.time=NULL,id=NULL,clusters=NULL,n.sim=1000,
residuals=0,robust=1,Nit=20,bandwidth=0.5,method="basic",
weighted.test=0,degree=1,covariance=0)
}
\arguments{
\item{formula}{
  a formula object with the response on the left of a '~' operator, and
  the independent terms on the right as regressors. The response must be
  a survival object as returned by the `Surv' function. Time-invariant
  regressors are specified by the wrapper const(), and cluster variables
  (for computing robust variances) by the wrapper cluster().
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
model validation in the function cum.residuals
}
\item{covariance}{to compute covariance estimates for
nonparametric terms rather than just the variances.}
\item{Nit}{number of iterations for score equations. }
\item{bandwidth}{bandwidth for local iterations. Default is 50 \% 
of the range of the considered observation period. 
}
\item{method}{Method for estimation. This refers to different parametrisations
of the baseline of the model. Options are "basic" where the baseline
is written as \eqn{\lambda_0(t) = \exp(\alpha_0(t))}
or the "breslow" version where the baseline is parametrised as 
\eqn{\lambda_0(t)}.
}
\item{degree}{gives the degree of the local linear smoothing, that is
local smoothing. Possible values are 1 or 2.
}
}
\details{
The data for a subject is presented as multiple rows or 'observations',
each of which applies to an interval of observation (start, stop]. 
When counting process data with the )start,stop] notation is used
the 'id' variable is needed to identify the records for each subject.
The program assumes that there are no ties, and if such are present 
random noise is added to break the ties.
}

\value{
Returns an object of type "timecox". With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.}
\item{var.cum}{the martingale based pointwise variance estimates.  }
\item{robvar.cum}{robust pointwise variances estimates.  }
\item{gamma}{estimate of parametric components of model.  }
\item{var.gamma}{variance for gamma.  }
\item{robvar.gamma}{robust variance for gamma.  }
\item{residuals}{list with residuals. Estimated martingale increments (dM) and
corresponding time vector (time).}
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
\item{schoenfeld.residuals}{Schoenfeld residuals are returned 
for "breslow" parametrisation.}
}
\references{
Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer (2006).
}
\author{Thomas Scheike}
\examples{
data(sTRACE)
# Fits time-varying Cox model 
out<-timecox(Surv(time/365,status==9)~age+sex+diabetes+chf+vf,
data=sTRACE,max.time=7,n.sim=100)

summary(out)
par(mfrow=c(2,3))
plot(out)
par(mfrow=c(2,3))
plot(out,score=TRUE)

# Fits semi-parametric time-varying Cox model
out<-timecox(Surv(time/365,status==9)~const(age)+const(sex)+
const(diabetes)+chf+vf,data=sTRACE,max.time=7,n.sim=100)

summary(out)
par(mfrow=c(2,3))
plot(out)
}
\keyword{survival}
