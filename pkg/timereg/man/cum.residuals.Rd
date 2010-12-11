\name{cum.residuals}
\alias{cum.residuals}
\title{Model validation based on cumulative residuals}
\description{Computes cumulative residuals and
	approximative p-values based on resampling techniques.}
\usage{
cum.residuals(object,data=sys.parent(),modelmatrix=0,cum.resid=1,
n.sim=500,weighted.test=1,start.design=1)
	}
\arguments{
\item{object}{an object of class 'aalen', 'timecox', 'cox.aalen' where
the residuals are returned ('residuals=1')}
\item{data}{data frame based on which residuals are computed.}
\item{modelmatrix}{specifies a grouping of the data that 
is used for cumulating residuals. Must have same size as data and be
ordered in the same way.}
\item{n.sim}{number of simulations in resampling.}
\item{weighted.test}{to compute a variance weighted version of the
test-processes used for testing constant effects of covariates.}
\item{cum.resid}{to compute residuals versus each of
the continuous covariates in the model.}
\item{start.design}{if '1' the groupings specified in modelmatrix
changes over time, i.e. in the case with time-dependent 
covariates.}
}
\value{
returns an object of type "cum.residuals" with the following arguments:
\item{cum}{cumulative residuals versus time for the groups specified by
modelmatrix.
}
\item{var.cum}{the martingale based pointwise variance estimates.}
\item{robvar.cum}{robust pointwise variances estimates of cumulatives.}
\item{obs.testBeq0}{observed absolute value of supremum of 
cumulative components scaled with the variance.}
\item{pval.testBeq0}{p-value covariate effects based on supremum test.}
\item{sim.testBeq0}{resampled supremum value.}
\item{conf.band}{resampling based constant to construct robust 95\% 
uniform confidence bands for cumulative residuals.}
\item{obs.test}{absolute value of supremum of observed test-process.}
\item{pval.test}{p-value for supremum test statistic.}
\item{sim.test}{resampled absolute value of supremum cumulative residuals.}
\item{proc.cumz}{observed cumulative residuals versus all continuous
covariates of model.}
\item{sim.test.proccumz}{list of 50 random realizations of 
test-processes under model for all continuous covariates.}
}
\references{
Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer (2006).
}
\author{Thomas Scheike}
\examples{
library(survival)
data(sTRACE)
# Fits Aalen model and returns residuals
fit<-aalen(Surv(time,status==9)~age+sex+diabetes+chf+vf,
sTRACE,max.time=7,n.sim=0,residuals=1)

# constructs and simulates cumulative residuals versus age groups
fit.mg<-cum.residuals(fit,sTRACE,n.sim=100,
modelmatrix=model.matrix(~-1+factor(cut(age,4)),sTRACE))

par(mfrow=c(1,4))
# cumulative residuals with confidence intervals
plot(fit.mg);
# cumulative residuals versus processes under model
plot(fit.mg,score=1); 
summary(fit.mg)

# cumulative residuals vs. covariates Lin, Wei, Ying style 
fit.mg<-cum.residuals(fit,sTRACE,cum.resid=1,n.sim=100)

par(mfrow=c(2,4))
plot(fit.mg,score=2)
summary(fit.mg)
}
\keyword{survival}
