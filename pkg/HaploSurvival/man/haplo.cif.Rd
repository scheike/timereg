\name{haplo.cif}
\alias{haplo.cif}
\title{Competings Risks Regression for haplotype effects} 
\description{
Fits a semiparametric model for the cause-specific quantities :
\deqn{
P(T \leq t, cause=1 | x) = P_1(t,x) = 1 - \exp(- g(t,x(h)) )
}
for the probability of dying from cause 1 in a situation with 
competing causes of death where  
\eqn{x(h)} is a design vector that depends on the unobserved haplotype pair. We only
observe the related genotype. The haplotype pair is integer coded using the ordered 
haplotypes (uniqueHaploNames from the geno.setup() function). 

The model is considered in two situations : 
\deqn{
g(t,x,z,h) = x(h)^T A(t) + (diag(t^p) z(h))^T \beta
}
the additive situation (with additive subdistribution hazard) and 
the proportional setting that includes the Fine & Gray (FG) model and
some extensions
\deqn{
g(t,x,z,h) = \exp(x(h)^T A(t) + (diag(t^p) z(h))^T \beta)
}
The FG model is obtained when \eqn{x=1}. Where p is 1 for the additive
model and 0 for the proportional model. In general p may be
powers of the same length as z. 

The marginal model, given genotype g and covariates x,z are then given as 
\deqn{
P(T \leq t, cause=1 | x,z,g) =  \sum_h P_1(t,x,z,h) P_\theta(H=h|G=g) 
}
and this is the model that is fitted to the data. 

The call of the function must contain information about the MLE for
the underlying haplotype frequencies under HWE equilibrium 
\deqn{ P(H=(h_k,h_j)) = \pi_k pi_j }
with K different haplotypes, where we parametrize the frequencies as
\deqn{ 
\pi_k= \frac{ \exp(\theta_k) }{1+\theta_1+...+\theta_{K-1}+0}
}
We allow a regression design on the haplotype parameters to reduce the 
dimensionality 
\deqn{ \theta = G \alpha }

Resampling is used for computing p-values for tests of time-varying effects.
}
\usage{
haplo.cif(formula,data=sys.parent(),
cause,times,designfuncX,designfuncZ,Nit=50,match=FALSE,
clusters=NULL,gamma=0,n.sim=500,weighted=0,model="additive",
causeS=1,cens.code=0,detail=0,interval=0.01,resample.iid=1,
cens.model="KM",time.pow=0,fix.haplofreq=0,haplo.freq=NULL,
alpha.iid=NULL,geno.setup=NULL,fit.haplofreq=NULL,
design.test=0,covnamesX=NULL,covnamesZ=NULL)
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
	considered. This is typically all cause "1" jump times.}
\item{designfuncX}{For nonparametric model R-function(x,h) that gives the X(h) 
design given x,z, and for semiparametric model 
R-function(x,z,h) that gives the X(h) design given x,z and h. See ?design for more details on this.}
\item{designfuncZ}{R-function(x,z,h) that gives the Z design given x,z and h.}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{match}{if match is true it the matched survival model, see haplomatch.cif()}
\item{clusters}{specifies cluster structure, for backwards compability.}
\item{gamma}{starting value for constant effects.}
\item{n.sim}{number of simulations in resampling.}
\item{weighted}{Not implemented. To compute a variance weighted version of the test-processes used for testing time-varying effects.}
\item{model}{"additive" or "prop"ortional.}
\item{causeS}{specificies which cause we consider.}
\item{cens.code}{specificies the code for the censoring.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{interval}{specifies that we only consider timepoints where the Kaplan-Meier of the censoring distribution  is larger than this value.}
\item{resample.iid}{to return the iid decomposition, that can be used to construct confidence bands for predictions}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"}
\item{time.pow}{specifies that the power at which the time-arguments  is transformed, for each of the arguments of the const() terms, default is 1 for the additive model and 0 for the proportional model.}
\item{fix.haplofreq}{is 1 when haplofrequencies are considered as known. Default is 0}
\item{haplo.freq}{the known haplofrequencies, or estimated MLE.}
\item{alpha.iid}{iid decomposition of haplofrequency paramters for uncertainty
related to haplotypefrequencies.}
\item{geno.setup}{analysis of genotype by geno.setup() function}
\item{fit.haplofreq}{MLE fit of haplofrequency model using HWE.}
\item{design.test}{prints out the design matrices  so one can check that the design is
correctly specified in the functions, designfuncX and designfuncX.}
\item{covnamesX}{names for output related to additive part of model.}
\item{covnamesZ}{names for output related to proportional part of model.}
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
}
\references{
Scheike and  Zhang and Gerds (2007), Predicting cumulative
incidence probability by direct binomial regression,
Biometrika, 2008. 

Scheike and Zhang (2007), Competing risks with missing covariates: Effect of matching
haplotypes on BMT patients, work in progress. 
}
\author{Thomas Scheike}
\examples{
library(timereg)
data(bmt); 
times<-bmt$time[bmt$cause==1]; 

# simulated genotype data 
Gbmt<-matrix(rbinom(4*nrow(bmt),1,0.5),nrow(bmt),4)
setup<-geno.setup(Gbmt)

# fits MLE under HWE 
out<-haplo.freqs(Gbmt,geno.setup=setup) 

designX<-function(x,h) { 
hh<-(h[1]==0)+(h[2]==0) # count occurrences of haplotype "0"
y<-c(x,hh)
return(y)
}

# nonparametric additive model 
add<-haplo.cif(Surv(time,cause>0)~platelet+age+tcell,bmt,
bmt$cause,times[-c(1:5)],designX,designZ,geno.setup=setup,Nit=40,
causeS=1,fit.haplofreq=out)
summary(add)

par(mfrow=c(2,4))
plot(add)

pout<-predict(add,X=rbind(
c(1,0,0,0,1), # with 1 haplotype 0
c(1,0,0,0,0)  # without haplotype 0
))

par(mfrow=c(1,3))
plot(pout,multiple=0,col=1:2,lty=1:2)
plot(pout,multiple=1,uniform=0,col=1:2,lty=1,se=0)

# fits semiparmatric model, first specifies design vectors

designZ<-function(x,z,h) {
h<-round(h); 
hh<-(h[1]==0)+(h[2]==0) # count occurrences of haplotype "0"
y<-c(z,hh) 
return(y)
}

designX<-function(x,z,h) { 
hh<-(h[1]==0)+(h[2]==0) # count occurrences of haplotype "0"
y<-c(x,hh)
return(x) # returns x so no nonpar haplotype effect
}

sadd<-haplo.cif(Surv(time,cause>0)~const(platelet)+const(age)+const(tcell),bmt,
bmt$cause,times[-c(1:5)],designX,designZ,geno.setup=setup,Nit=40,
causeS=1,fit.haplofreq=out)
summary(sadd)

par(mfrow=c(2,2))
plot(sadd)

pout<- predict(sadd, X=rbind(c(1),c(1)),
Z=rbind( c(0,0,0,1), # with 1 haplotype 0
         c(0,0,0,0))  # without haplotype 0
)

plot(pout,multiple=1,uniform=0,col=1:2,lty=1,se=0)

plot(pout,multiple=0,uniform=1,col=1:2,lty=1:2,se=1)
}
\keyword{survival}
