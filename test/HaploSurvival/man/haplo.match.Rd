\name{haplo.match}
\alias{haplo.match}
\title{Fit Cox-Aalen survival model for haplo-match effects for BMT data with 
donor and patient with matching genotypes.}
\description{
Fits an Cox-Aalen survival model. 
Time dependent variables and counting process data
(multiple events per subject) are possible.

\deqn{ \lambda_{i}(t) = Y_i(t) 
( X_{i}(hd,hp)^T(t) \alpha(t) )  \exp(Z_{i}(hd,hp)^T \beta ) 
}
where \eqn{X_i(hd,hp)} and \eqn{Z_i(hd,hp)} are the design vectors
specified by the design functions given in the call. The function
is designed for bone-marrow-transplant data where a donor and patient
are matched on genotype and we wish to consider the effect of their
haplotype pairs.

Resampling is used for computing p-values for tests of time-varying effects.
Test for proportionality is considered by
considering the score processes for the proportional effects of model.

The modelling formula uses the standard survival modelling given in the
\bold{survival} package.

The haplotype model assumes Hardy-Weinberg equilibrium 
\deqn{ P(H=(h_k,h_j)) = \pi_k pi_j }
with K different haplotypes, where we parametrize the frequencies as
\deqn{ 
\pi_k= \frac{ \exp(\theta_k) }{1+\theta_1+...+\theta_{K-1}+0}
}
We allow a regression design on the haplotype parameters to reduce the 
dimensionality 
\deqn{ \theta = G \alpha }
}
\usage{
haplo.match(formula = formula(data),data=sys.parent(),
designfuncX,designfuncZ, beta = 0, Nit = 10, detail = 0, 
start.time = 0, max.time = NULL, id = NULL, n.sim = 500, 
geno.type = NULL, geno.setup=NULL, haplo.freq = NULL,
fix.beta = 0, fix.haplofreq = 0, two.stage = 0, weighted.test = 0,
step = 1, lev.marq=1, min.lev.marq=0,
haplo.design=NULL,haplo.baseline=NULL,alpha=NULL,resample.iid=1,
covnamesX=NULL,covnamesZ=NULL)
}
\arguments{
\item{formula}{
  a formula object with the response on the left of a '~' operator, and
  the independent terms on the right as regressors. The response must be
  a survival object as returned by the `Surv' function. Terms with a
  proportional effect are specified by the wrapper prop(). 
  Haplotype effects specified using the two designfunctions. 
}
\item{data}{a data.frame with the variables.}
\item{designfuncX}{R-function(x,z,h) that gives the X design given x,z and h.
See ?design for more details on this.}
\item{designfuncZ}{R-function(x,z,h) that gives the Z design given x,z and h.}
\item{beta}{starting value for relative risk estimates}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details is printed during iterations, if 1 details are given.}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{end of observation period where estimates are computed. 
Estimates thus computed from [start.time, max.time]. Default is max of data.}
\item{id}{For timevarying covariates the variable must 
associate each record with the id of a subject.}
\item{n.sim}{number of simulations in resampling.}
\item{geno.type}{geno types associated with subjects.}
\item{geno.setup}{genotype setup based on genotype function.}
\item{haplo.freq}{haplo type frequencies for starting values or for fixed values.}
\item{fix.beta}{option for keeping beta fixed.}
\item{fix.haplofreq}{option for keeping haplotype frequencies fixed.}
\item{two.stage}{option to specify that genotype likelihood is based only on genotype likelihood, if two.stage=0 score equations for haplotype parameters also uses information from survival data.}
\item{weighted.test}{to compute a variance weighted version of the 
test-processes used for testing time-varying effects.}
\item{step}{step size for iteration steps.}
\item{lev.marq}{starting value for Levenberg-Marquart algorithm.}
\item{min.lev.marq}{stopping value for Levenberg-Marquart algorithm after initial steps, min.lev.marq=0 is stopping value for Levenberg-Marquart algorithm after initial steps, min.lev.marq=0 is
Newton-Raphson steps, that are needed to get correct standard errors.}
\item{haplo.design}{design for haplo parameters, default is diagonal of size K-1}
\item{haplo.baseline}{gives the name of the geno.setup that should be used as baseline.}
\item{alpha}{starting value of parameters for  haplo.design}
\item{resample.iid}{set to 1 get i.i.d. decomposition of baseline and  beta.}
\item{covnamesX}{names for output related to additive part of model.}
\item{covnamesZ}{names for output related to proportional part of model.}
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
computed within the estimation interval. }
\item{var.cum}{the martingale based pointwise variance estimates.  }
\item{robvar.cum}{robust pointwise variances estimates.  }
\item{gamma}{estimate of parametric components of model.  }
\item{var.gamma}{variance for gamma.}
\item{robvar.gamma}{robust variance for gamma.  }
\item{haplo.pars}{estimates for haplotype frequency parameters.}
\item{haplo.freqs}{estimates for haplotype frequencies based on parameters.}
\item{haplo.alpha}{estimates for haplotype frequencies based on haplo.design.}
\item{freq.names}{names of haplotypes.}
\item{var.haplo.alpha}{variance matrix for alpha}
\item{robvar.haplo.alpha}{variance matrix for alpha}
\item{D2linv.all}{inverse of derivative of score for all parameters.}
\item{D2linv.gamma}{inverse of derivative of score for relative risk parameters.}
\item{D2linv.haplo.alpha}{inverse of derivative of score for alpha (haplo-frequency) parameters.}
\item{score}{score of last iteration.}
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
\item{conf.band}{resampling based constant to construct robust 95\% uniform confidence bands. }
\item{test.procBeqC}{observed test-process of difference between observed cumulative 
process and estimate under null of constant effect over time.  }
\item{sim.test.procBeqC}{list of 50 random realizations of
test-processes under null based on resampling.}
\item{covariance}{covariances for nonparametric terms of model.}
\item{B.iid}{Resample processes for nonparametric terms of model.}
\item{gamma.iid}{Resample processes for parametric terms of model.}
\item{D2linv}{inverse of the derivative of the score function.}
\item{score}{value of score for final estimates.}
\item{test.procProp}{observed score process for proportional part of model.}
\item{var.score}{variance of score process (optional variation estimator
for beta.fixed=1 and robust estimator otherwise).}
\item{pval.Prop}{p-value based on resampling.}
\item{sim.supProp}{re-sampled absolute supremum values.}
\item{sim.test.procProp}{list of 50 random realizations of 
test-processes for proportionality under the model based on resampling.}
\item{haplo.design}{design that relates haplotype parameters to alpha parameters.}
}
\references{
Scheike, Martinussen and Silver,  Haplotype effects for survival data, Submitted.

Scheike, Martinussen and Zhang,  Missing data for the additive risk model: Estimating haplotype effects, Submitted.
}
\author{Thomas Scheike}
\examples{
### Simulating data 
n<-400; 
haplo.types<-rbind(c(0,0),c(1,0),c(0,1),c(1,1))
### randomly selects haplotypes
pairs<-matrix(sample(1:4,n*2,replace=TRUE),n,2) 

gp<-matrix(haplo.types[pairs,],n,4) ### observed genotype for patient 
                                    ### haplotype-pair 1=g[c(1,3),]
                                    ### haplotype-pair 2=g[c(2,4),]
gs<-geno.setup(gp); 

gd<-sim.donor(gp,gs,haplo.types);### draws donor haplotype pairs such
                                 ### genotype match, and returns 0/1 for match 

Z<-matrix(rnorm(n),n,1); Z<-cbind(Z,Z*0.5+1+rnorm(n))
X<-matrix(rbinom(n,1,0.5),n,1); X<-cbind(1,X);

time<-exp(Z \%*\% c(-0.3,0.3) + gd * 0.5 ) * rexp(n)/ X \%*\% c(0.2,0.5)
status<-(time<3); time[status==0]<-3
simdata=data.frame(time=time,status=status,x1=X[,1],x2=X[,2],
             z1=Z[,1],z2=Z[,2],zh=gd)

### specifies Cox-Aalen model 
### design for baseline, does not depend on haplotype
designX<-function(x,z,hd,hp) 
\{
return(x)
\}

### design for proportional part of model 
designZ<-function(x,z,hd,hp) 
\{
hd<-round(hd); hp<-round(hp); 
zh<-(hd[1]==hp[1])*(hp[2]==hd[2]) ### 1 if match and 0 otherwise
y<-c(z,zh)
return(y)
\}
designZ(c(1,3),c(1,1),c(0,0),c(0,1)); 
designZ(c(1,3),c(1,1),c(1,0),c(1,0))

hapfit<-haplo.freqs(g,geno.setup=gs)

### MLE starting value of haplotype frequencies 
out<-haplo.match(Surv(time,status)~x2+prop(z1)+prop(z2),
data=simdata,designfuncX=designX,designfuncZ=designZ,geno.setup=gs,
two.stage=1,alpha=hapfit$haplo.alpha) 

summary(out);
plot(out,sim.ci=1) ### plots baseline with confidence bands
plot(out,score=1)  ### score process test for proportionality
out$haplo.freq     ### haplotype frequencies in model 

### predicts survival with Z=0 and with 0 and 2 number of "0,0"
pout<-predict(out,X=rbind(c(1,1),c(1,1)),Z=rbind(c(0,0,0),c(0,0,2)))
plot(pout)                 ### the two survival curves
plot(pout,uniform=1,se=1)  ### with 95 % confidence bands and intervals

######################################################################
########## structured haplotype model ################################
######################################################################
X<-matrix(1,3,1); 
hapfit<-haplo.freqs(g,geno.setup=gs,haplo.design=X)

out<-haplo.match(Surv(time,status)~x2+prop(z1)+prop(z2),data=simdata,
designfuncX=designX,designfuncZ=designZ,geno.setup=gs,
two.stage=1,haplo.design=X,alpha=hapfit$haplo.alpha)
summary(out);
out$haplo.freq ### haplotype frequencies in model 
par(mfrow=c(1,2))
plot(out,sim.ci=1) 

#######################################################################
########### Fits additive hazard model ################################
#######################################################################

### design for baseline, does depend on haplotype
designXa<-function(x,z,hd,hp) 
\{ 
hd<-round(hd); hp<-round(hp); 
zh<-(hd[1]==hp[1])*(hp[2]==hd[2]) ### 1 if match and 0 otherwise
y<-c(x,zh) 
return(y)
\}

### returns proportional part of model
designZa<-function(x,z,hd,hp) 
\{ 
return(z)
\}

### Nit=1 because we start with the correct MLE haplotype estimates 

out<-haplo.match(Surv(time,status)~x2,data=simdata,
designfuncX=designXa,designfuncZ=designZa,geno.setup=gs,Nit=1,
fix.beta=1,two.stage=1,haplo.design=X,alpha=hapfit$haplo.alpha)
summary(out);
par(mfrow=c(1,3)); 
plot(out,sim.ci=1) ### plots baseline with confidence bands

### predicts survival with Z=0 and with 0 and 2 number of "0,0"
pout<-predict(out,X=rbind(c(1,1,0),c(1,1,2)),Z=rbind(0,0))
par(mfrow=c(1,3))
plot(pout,multiple=1,se=0,uniform=0,col=2:3) 
plot(pout,uniform=1,se=1,transparency=1)   
}
\keyword{survival}
