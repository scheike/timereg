\name{comp.risk}
\alias{comp.risk}
\title{Competings Risks Regression} 
\description{
Fits a semiparametric model for the cause-specific quantities :
\deqn{
P(T < t, cause=1 | x,z) = P_1(t,x,z) = h(  g(t,x,z) )
}
for a known link-function \eqn{h()}
and known prediction-function \eqn{g(t,x,z)}
for the probability of dying from cause 1 in a situation with 
competing causes of death. 

We consider the following models :  
1) the additive model where
\eqn{h(x)=1-\exp(-x)} and 
\deqn{
g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta
}
2)  
the proportional setting that includes the Fine & Gray (FG) "prop" model and
some extensions where 
\eqn{h(x)=1-\exp(-\exp(x))} and 
\deqn{
g(t,x,z) = \exp(x^T A(t) + (diag(t^p) z)^T \beta)
}
The FG model is obtained when \eqn{x=1}. 

3) a "logistic" model where 
\eqn{h(x)=\exp(x)/( 1+\exp(x))} and 
\deqn{
g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta
}

4) the relative cumulative incidence function "rcif" model where
\eqn{h(x)=\exp(x)} and 
\deqn{
g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta
}) 
Where p by default is 1 for the additive model and 0 for 
the other models.
In general p may be powers of the same length as z. 
}
\usage{
comp.risk(formula,data=sys.parent(),cause,times=NULL,Nit=50,
clusters=NULL,est=NULL,fix.gamma=0,gamma=0,n.sim=500,weighted=0,model="additive",
causeS=1,cens.code=0,detail=0,interval=0.01,resample.iid=1,
cens.model="KM",time.pow=NULL,time.pow.test=NULL,silent=1,conv=1e-6,
weights=NULL,max.clust=1000,n.times=50,first.time.p=0.05,
trunc.p=NULL,entry.time=NULL,cens.weight=NULL,admin.cens=NULL)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a '~'
  operator, and the terms on the right. The response must be a survival
  object as returned by the `Surv' function. The status indicator is not
  important here. Time-invariant regressors are specified by the wrapper
  const(), and cluster variables (for computing robust variances) by the
  wrapper cluster().
}
\item{data}{a data.frame with the variables.}
\item{cause}{specifies the causes  related to the death
	times, the value 0 is the censoring value.}
\item{times}{specifies the times at which the estimator is
	considered. The default is cause "1" jump times, with the first 
 10 percent or max 20 jump points removed for numerical stability in simulations.}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{clusters}{specifies cluster structure, for backwards compability.}
\item{est}{possible starting value for nonparametric component of model.}
\item{fix.gamma}{to keep gamma fixed, possibly at 0.}
\item{gamma}{starting value for constant effects.}
\item{n.sim}{number of simulations in resampling.}
\item{weighted}{Not implemented. To compute a variance weighted version of the test-processes used for testing time-varying effects.}
\item{model}{"additive", "prop"ortional, "rcif", or "logistic".}
\item{causeS}{specificies which cause we consider.}
\item{cens.code}{specificies the code for the censoring.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{interval}{specifies that we only consider timepoints where the Kaplan-Meier of the censoring distribution  is larger than this value.}
\item{resample.iid}{to return the iid decomposition, that can be used to construct confidence bands for predictions}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"}
\item{time.pow}{specifies that the power at which the time-arguments  is transformed, for each of the arguments of the const() terms, default is 1 for the additive model and 0 for the proportional model.}
\item{time.pow.test}{specifies that the power the time-arguments is 
transformed for each of the arguments of the non-const() terms. This is 
relevant for testing if a coefficient function is consistent with the 
specified form A_l(t)=beta_l t^time.pow.test(l). Default is 1 for the 
additive model and 0 for the proportional model.}
\item{silent}{if 0 information on convergence problems due to non-invertible derviates of scores are printed.}
\item{conv}{gives convergence criterie in terms of sum of absolute change of parameters of model}
\item{weights}{weights for estimating equations.}
\item{max.clust}{sets the total number of i.i.d. terms in i.i.d. decompostition. This can limit the amount of memory used by coarsening 
the clusters. When NULL then all clusters are used. 
Default is 1000 to save memory and time. } 
\item{first.time.p}{first point for estimation is pth percentile of causeS jump times.}
\item{n.times}{only uses 50 points for estimation, if NULL then uses all points, subject to p.start condition.}
\item{trunc.p}{truncation weight for delayed entry, P(T > entry.time | Z_i), typically Cox model.}
\item{entry.time}{entry times}
\item{cens.weight}{censoring weights can be given here rather than calculated using the KM, cox or aalen models.}
\item{admin.cens}{censoring times for the administrative censoring}
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
\item{obs.testBeqC}{observed absolute value of supremum of difference between observed cumulative process and estimate under null of constant effect.}
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
Scheike, Zhang and Gerds (2008), Predicting cumulativeincidence probability by direct binomial regression,Biometrika, 95, 205-220.           

Scheike and Zhang (2007), Flexible competing risks regression modelling and goodness of fit, LIDA, 14, 464-483.

Martinussen and Scheike (2006), Dynamic regression models for survival data, Springer.
}
\author{Thomas Scheike}
\examples{
library(timereg)
data(bmt); 

add<-comp.risk(Surv(time,cause>0)~platelet+age+tcell,data=bmt,
bmt$cause,causeS=1,resample.iid=1,n.sim=100)
summary(add)

par(mfrow=c(2,4))
plot(add); plot(add,score=1)

ndata<-data.frame(platelet=c(1,0,0),age=c(0,1,0),tcell=c(0,0,1))
par(mfrow=c(2,3))
out<-predict(add,ndata,uniform=1,n.sim=100)
par(mfrow=c(2,2))
plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=1)

## fits additive model with some constant effects 
add.sem<-comp.risk(Surv(time,cause>0)~
const(platelet)+const(age)+const(tcell),data=bmt,
bmt$cause,causeS=1,resample.iid=1,n.sim=100)
summary(add.sem)

out<-predict(add.sem,ndata,uniform=1,n.sim=100)
par(mfrow=c(2,2))
plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=0)

## Fine & Gray model 
fg<-comp.risk(Surv(time,cause>0)~
const(platelet)+const(age)+const(tcell),data=bmt,
bmt$cause,causeS=1,resample.iid=1,model="prop",n.sim=100)
summary(fg)

out<-predict(fg,ndata,uniform=1,n.sim=100)

par(mfrow=c(2,2))
plot(out,multiple=1,uniform=0,col=1:3,lty=1,se=0)

## extended model with time-varying effects
fg.npar<-comp.risk(Surv(time,cause>0)~platelet+age+const(tcell),
data=bmt,bmt$cause,causeS=1,resample.iid=1,model="prop",n.sim=100)
summary(fg.npar); 

out<-predict(fg.npar,ndata,uniform=1,n.sim=100)
head(out$P1[,1:5]); head(out$se.P1[,1:5])

par(mfrow=c(2,2))
plot(out,multiple=1,uniform=0,col=1:3,lty=1,se=0)

## Fine & Gray model with alternative parametrization
fg<-comp.risk(Surv(time,cause>0)~
const(platelet)+const(age)+const(tcell),data=bmt,
bmt$cause,causeS=1,resample.iid=1,model="fg",n.sim=100)
summary(fg)

## Left truncated cumulative incidence estimation 
tt <- runif(408)*bmt$time*rbinom(408,1,0.5)*0.5
## truncation weights  computed first, by Cox model
fgt<- cox.aalen(Surv(tt,time,cause>0)~prop(platelet)+prop(age)+prop(tcell),data=bmt,n.sim=0,robust=0,silent=0)
cumt <- Cpred(fgt$cum,tt)[,2]  
psurv <- exp(-cumt); 


fgt<-comp.risk(Surv(time,cause>0)~ const(platelet)+const(age)+const(tcell),data=bmt,
bmt$cause,causeS=1,model="fg",n.sim=0,trunc.p=psurv,entry.time=tt)
summary(fgt)
pfgt <- predict(fgt,X=1,Z=c(0,0,0))
pfgt <- predict(fgt,X=1)
plot(pfgt)

## Left truncated cumulative incidence estimation 
## without covariates but using dummy covariate, equivalent with Aalen-Johansen estimator
fgt<- aalen(Surv(tt,time,cause>0)~+1,data=bmt,n.sim=0,robust=0)
cumt <- Cpred(fgt$cum,tt)[,2]
psurv <- exp(-cumt); 

pcif<-comp.risk(Surv(time,cause>0)~const(platelet),data=bmt,
bmt$cause,fix.gamma=1,causeS=1,model="fg",n.sim=0,trunc.p=psurv,entry.time=tt)
pout <- predict(pcif,X=1)
}
\keyword{survival}
