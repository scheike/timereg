\name{prop.odds.subdist}
\alias{prop.odds.subdist}
\title{Fit Semiparametric Proportional 0dds Model for the competing risks subdistribution}
\description{
Fits a semiparametric  proportional odds model:
\deqn{
logit(F_1(t;X,Z)) = log( A(t)) + \beta^T Z 
}
where A(t) is increasing but otherwise unspecified. Model is 
fitted by maximising the modified partial likelihood. 
A goodness-of-fit test by
considering the score functions is also computed by resampling methods. 

An alternative way of writing the model :
\deqn{
F_1(t;X,Z)  = \frac{ \exp(  \beta^T Z )}{
 (A(t)) + \exp(  \beta^T Z)  }
}
such that \eqn{\beta} is the log-odds-ratio of cause 1 before time t,
and \eqn{A(t)} is the odds-ratio.

The modelling formula uses the standard survival modelling given in the
\bold{survival} package.
}
\usage{
prop.odds.subdist(formula,data=sys.parent(),cause=1,beta=NULL,
Nit=10,detail=0,start.time=0,max.time=NULL,id=NULL,n.sim=500,weighted.test=0,
profile=1,sym=0,cens.model="KM",clusters=NULL,max.clust=1000,baselinevar=1,weights=NULL,
cens.weights=NULL)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a '~' operator, and
the terms on the right.  The response must be an object as
returned by the `Event' function.
}
j\item{data}{a data.frame with the variables.}
\item{cause}{cause indicator for competing risks.}
\item{beta}{starting value for relative risk estimates}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details is printed during iterations, if 1 details are given.}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{ end of observation period where estimates are computed. 
       	Estimates thus computed from [start.time, max.time]. 
       This is very useful to obtain stable estimates, especially for the baseline. Default is max of data.}
\item{id}{For timevarying covariates the variable must associate each record with the id of a subject.}
\item{n.sim}{number of simulations in resampling.}
\item{weighted.test}{to compute a variance weighted version of the test-processes used for 
	testing time-varying effects.}
\item{profile}{use profile version of score equations.}
\item{sym}{to use symmetrized second derivative in the case of the estimating equation approach (profile=0). 
This may improve the numerical performance.  }
\item{cens.model}{specifies censoring model. So far only Kaplan-Meier "KM".}
\item{clusters}{to compute cluster based standard errors.}
\item{max.clust}{number of maximum clusters to be used, to save time in iid decomposition.}
\item{baselinevar}{set to 0 to save time on computations.}
\item{weights}{additional weights.}
\item{cens.weights}{specify censoring weights related to the observations.}
}
\details{
The data for a subject is presented as multiple rows or "observations",
each of which applies to an interval of observation (start, stop].
The program essentially assumes no ties, and if such are present 
a little random noise is added to break the ties.
}
\value{
returns an object of type 'cox.aalen'. With the following arguments:

\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.}
\item{var.cum}{the martingale based pointwise variance estimates.  }
\item{robvar.cum}{robust pointwise variances estimates.  }
\item{gamma}{estimate of proportional odds parameters of model.}
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
process and estimate under null of constant effect over time.}
\item{loglike}{modified partial likelihood, pseudo profile likelihood 
for regression parameters.}
\item{D2linv}{inverse of the derivative of the score function.}
\item{score}{value of score for final estimates.}
\item{test.procProp}{observed score process for proportional odds 
regression effects.}
\item{pval.Prop}{p-value based on resampling.}
\item{sim.supProp}{re-sampled supremum values.}
\item{sim.test.procProp}{list of 50 random realizations of
test-processes for constant proportional odds under the model 
based on resampling.}
}
\references{
Eriksson, Li, Zhang and Scheike (2014), The proportional odds cumulative incidence model for competing risks,
Biometrics, to appear.

Scheike, A flexible semiparametric transformation model for
survival data, Lifetime Data Anal. (2007). 

Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer (2006).
}
\author{Thomas Scheike}
\examples{
library(timereg)
data(bmt)
# Fits Proportional odds model 
out <- prop.odds.subdist(Event(time,cause)~platelet+age+tcell,data=bmt,
 cause=1,cens.model="KM",detail=0,n.sim=1000)
summary(out) 
par(mfrow=c(2,3))
plot(out,sim.ci=2); 
plot(out,score=1) 

# simple predict function without confidence calculations 
pout <- predictpropodds(out,X=model.matrix(~platelet+age+tcell,data=bmt)[,-1])
matplot(pout$time,pout$pred,type="l")

# predict function with confidence intervals
pout2 <- predict(out,Z=c(1,0,1))
plot(pout2,col=2)
pout1 <- predictpropodds(out,X=c(1,0,1))
lines(pout1$time,pout1$pred,type="l")

# Fits Proportional odds model with stratified baseline, does not work yet!
###out <- Gprop.odds.subdist(Surv(time,cause==1)~-1+factor(platelet)+
###prop(age)+prop(tcell),data=bmt,cause=bmt$cause,
###cens.code=0,cens.model="KM",causeS=1,detail=0,n.sim=1000)
###summary(out) 
###par(mfrow=c(2,3))
###plot(out,sim.ci=2); 
###plot(out,score=1) 
}
\keyword{survival}
