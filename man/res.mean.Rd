\name{res.mean}
\alias{res.mean}
\title{Residual mean life (restricted)} 
\description{
Fits a semiparametric model for the residual life (estimator=1):
\deqn{
	E( \min(Y,\tau) -t | Y>=t) = h_1(  g(t,x,z) )
}
or cause specific years lost of Andersen (2012) (estimator=3)
\deqn{
E( \tau- \min(Y_j,\tau) | Y>=0) = \int_0^t (1-F_j(s)) ds = h_2(  g(t,x,z) )
}
where \eqn{Y_j = \sum_j Y I(\epsilon=j) + \infty * I(\epsilon=0)}
or (estimator=2)
\deqn{
E( \tau- \min(Y_j,\tau) | Y<\tau, \epsilon=j) = h_3(  g(t,x,z) ) = h_2(g(t,x,z)) F_j(\tau,x,z)
}
where 
\eqn{F_j(s,x,z) = P(Y<\tau, \epsilon=j | x,z )}
for a known link-function \eqn{h()}
and known prediction-function \eqn{g(t,x,z)}

Uses the IPCW for the score equations based on 
\deqn{
w(t) \Delta(\tau)/P(\Delta(\tau)=1| T,\epsilon,X,Z) ( Y(t) - h_1(t,X,Z))
}
and where 
\eqn{\Delta(\tau)} is the at-risk indicator given data and requires a IPCW model.
}
\usage{
res.mean(formula,data=sys.parent(),cause=1,restricted=NULL,times=NULL,Nit=50,
clusters=NULL,gamma=0,n.sim=0,weighted=0,model="additive",detail=0,interval=0.01,
resample.iid=1,cens.model="KM",cens.formula=NULL,time.pow=NULL,time.pow.test=NULL,
silent=1,conv=1e-6,estimator=1,cens.weights=NULL,conservative=1,weights=NULL)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a '~'
  operator, and the terms on the right. The response must be a survival
  object as returned by the `Event' function. The status indicator is not
  important here. Time-invariant regressors are specified by the wrapper
  const(), and cluster variables (for computing robust variances) by the
  wrapper cluster().
}
\item{data}{a data.frame with the variables.}
\item{cause}{For competing risk models specificies which cause we consider.}
\item{restricted}{gives a possible restriction times for means.}
\item{times}{specifies the times at which the estimator is
	considered. Defaults to all the times where an event of interest occurs, with the first 
 10 percent or max 20 jump points removed for numerical stability in simulations.}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{clusters}{specifies cluster structure, for backwards compability.}
\item{gamma}{starting value for constant effects.}
\item{n.sim}{number of simulations in resampling.}
\item{weighted}{Not implemented. To compute a variance weighted version of the test-processes used for 
testing time-varying effects.}
\item{model}{"additive", "prop"ortional.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{interval}{specifies that we only consider timepoints where the Kaplan-Meier of the censoring distribution
is larger than this value.}
\item{resample.iid}{to return the iid decomposition, that can be used to construct confidence bands for predictions}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox" or
"aalen" model for further flexibility.}
\item{cens.formula}{specifies the regression terms used for the regression model for chosen regression model. 
When cens.model is specified, the default
is to use the same design as specified for the competing risks model. "KM","cox","aalen","weights". "weights" are user
specified weights given is cens.weight argument. 
}
\item{time.pow}{specifies that the power at which the time-arguments  is transformed, for each of the arguments of the const() terms, default is 1 for the additive model and 0 for the proportional model.}
\item{time.pow.test}{specifies that the power the time-arguments is 
transformed for each of the arguments of the non-const() terms. This is 
relevant for testing if a coefficient function is consistent with the 
specified form A_l(t)=beta_l t^time.pow.test(l). Default is 1 for the 
additive model and 0 for the proportional model.}
\item{silent}{if 0 information on convergence problems due to non-invertible derviates of scores are printed.}
\item{conv}{gives convergence criterie in terms of sum of absolute change of parameters of model}
\item{estimator}{specifies what that is estimated.}
\item{cens.weights}{censoring weights for estimating equations.}
\item{conservative}{for slightly conservative standard errors.}
\item{weights}{weights for estimating equations.}
}
\details{
Since timereg version 1.8.4. the response must be specified with the
\code{\link{Event}} function instead of the \code{\link{Surv}}
function and the arguments. 
}
\value{returns an object of type 'comprisk'. With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.}
\item{var.cum}{pointwise variances estimates.  }
\item{gamma}{estimate of proportional odds parameters of model.}
\item{var.gamma}{variance for gamma.  }
\item{score}{sum of absolute value of scores.}
\item{gamma2}{estimate of constant effects based on the
non-parametric estimate. Used for testing of constant effects.}
\item{obs.testBeq0}{observed absolute value of supremum of 
cumulative components scaled with the variance.}
\item{pval.testBeq0}{p-value for covariate effects based on supremum test.}
\item{obs.testBeqC}{observed absolute value of supremum of difference between observed cumulative process 
and estimate under null of constant effect.}
\item{pval.testBeqC}{p-value based on resampling.}
\item{obs.testBeqC.is}{observed integrated squared differences between 
observed cumulative and estimate under null of constant effect.}
\item{pval.testBeqC.is}{p-value based on resampling.}
\item{conf.band}{resampling based constant to construct 95\% uniform confidence bands.}
\item{B.iid}{list of iid decomposition of non-parametric effects.}
\item{gamma.iid}{matrix of iid decomposition of parametric effects.}
\item{test.procBeqC}{observed test process for testing of time-varying effects}
\item{sim.test.procBeqC}{50 resample processes for for testing of time-varying effects}
\item{conv}{information on convergence for time points used for estimation.}
}
\references{
Andersen (2013), Decomposition of number of years lost according to causes of death, Statistics in Medicine,  5278-5285.

Scheike, and Cortese (2015), Regression Modelling of Cause Specific Years Lost,  

Scheike, Cortese and Holmboe (2015), Regression Modelling of Restricted Residual Mean with Delayed Entry,  
}
\author{Thomas Scheike}
\examples{
data(bmt); 
tau <- 100 

### residual restricted mean life
out<-res.mean(Event(time,cause>=1)~factor(tcell)+factor(platelet),data=bmt,cause=1,
	      times=0,restricted=tau,n.sim=0,model="additive",estimator=1); 
summary(out)

out<-res.mean(Event(time,cause>=1)~factor(tcell)+factor(platelet),data=bmt,cause=1,
	      times=seq(0,90,5),restricted=tau,n.sim=0,model="additive",estimator=1); 
par(mfrow=c(1,3))
plot(out)

### restricted years lost given death
out21<-res.mean(Event(time,cause)~factor(tcell)+factor(platelet),data=bmt,cause=1,
	      times=0,restricted=tau,n.sim=0,model="additive",estimator=2); 
summary(out21)
out22<-res.mean(Event(time,cause)~factor(tcell)+factor(platelet),data=bmt,cause=2,
	      times=0,restricted=tau,n.sim=0,model="additive",estimator=2); 
summary(out22)


### total restricted years lost 
out31<-res.mean(Event(time,cause)~factor(tcell)+factor(platelet),data=bmt,cause=1,
	      times=0,restricted=tau,n.sim=0,model="additive",estimator=3); 
summary(out31)
out32<-res.mean(Event(time,cause)~factor(tcell)+factor(platelet),data=bmt,cause=2,
	      times=0,restricted=tau,n.sim=0,model="additive",estimator=3); 
summary(out32)


### delayed entry 
nn <- nrow(bmt)
entrytime <- rbinom(nn,1,0.5)*(bmt$time*runif(nn))
bmt$entrytime <- entrytime

bmtw <- prep.comp.risk(bmt,times=tau,time="time",entrytime="entrytime",cause="cause")

out<-res.mean(Event(time,cause>=1)~factor(tcell)+factor(platelet),data=bmtw,cause=1,
	      times=0,restricted=tau,n.sim=0,model="additive",estimator=1,
              cens.model="weights",weights=bmt$cw,cens.weights=1/bmtw$weights); 
summary(out)
}
\keyword{survival}
