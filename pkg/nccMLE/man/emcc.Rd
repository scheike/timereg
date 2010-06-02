\name{em.ncc}
\alias{em.ncc}
\title{Fits Cox's regression model for nested case-control data and 
case-cohort data}
\description{
Z covariates observed only for cases and controls 

X observed for all subjects in the study such as sex and age 

The distribution of Z given X needs is specified non-parametrically
on a number of strata specified by X through the stratum variable. 

Intensity model 
\deqn{
\lambda_{X,Z}(t)) = \lambda_S(X) \exp( \beta^T Z )
}
The baseline stratification based on \eqn{S(X)}  is specified
by baseline.stratum.

When fullobsX is given the baseline may be also on Cox form 
\deqn{
\lambda_{X}(t)) = \lambda_{S(X_1)} \exp( \beta_X^T X_2 )
}
where a stratum is used for some covariates and other covariates
are modelled as regression effects.  This model is specified by
fullobsX=\eqn{X_2}, and \eqn{S(X_1)}= baseline.stratum

When stratZ==1 the baseline is stratified by the Z components as 
specified by baseline.stratum=S(Z)
\deqn{
\lambda_{X,Z}(t)) = \lambda_S(Z) \exp( \beta^T X )
}

}
\usage{
em.ncc(Zobs,Tobs,status,tau,nno,Nit=1,beta=0,detail=0,
betait=2,em.dif=1,stratum=0,baseline=0,emvar=0,
fullobsZ=0,baseline.stratum=0,stratZ=0)
}
\arguments{
\item{Zobs}{covariates observed for cases and controls, fullobsX is
also given it means that the first columns of Zobs is the Z from
the above intensity and the second part is the covariates X.}
\item{Tobs}{observed survival times for cases and controls.}
\item{status}{censoring status for observed cases and controls.}
\item{tau}{end of observation period. If tau has length equal to 
the number of subjects that are not cases and controls it is individual 
censoring times.}
\item{nno}{number of subjects in the cohort that are not cases or
controls. }
\item{Nit}{number of itterations for EM algorithm.}
\item{beta}{starting value for regression parameter.}
\item{detail}{prints out iteration details.}
\item{betait}{number of itterations for Cox regression score in
EM, should be at least 2.}
\item{em.dif}{constant used for EM aided differentation, default is
(em.dif/number of subjects in cohort).}
\item{stratum}{case specific strata for number of strata for
specification of the distribution of Z given X.}
\item{baseline}{starting value for baseline estimates.}
\item{emvar}{computes EM based variance by EM aided differentation. }
\item{fullobsZ}{The fully observed covariates if these are needed
for a possible regression model for the baseline}
\item{baseline.stratum}{stratum for the baseline, defined as the
grouping given by S(X), a vector of length ntot=nno+nrow(Zobs).}
\item{stratZ}{stratifies after the Z covariate that is only
observed for cases and controls, stratification given by
baseline.stratum}
}
\value{
returns an object of type "em.ncc". With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.
}
\item{baseline}{the baseline estimates.}
\item{beta}{estimate of parametric components of model.}
\item{var.beta}{variance for beta.}
\item{p}{the distribution of covariates.}
\item{delta}{something about convergence.}
\item{em.dif}{parameter used for EM aided differentation.}
\item{covz}{the covariates realtes do p.}
\item{konv}{something about convergence.}
}
\references{
Scheike and Juul, Biostatistics

Scheike and Maritnussen, SJS 
}
\author{Thomas Scheike}
\examples{
###library(nccMLE)

# ud<-simNCC(3728,2,rate=1/0.004,win=1) realistic simulation

ud<-simNCC(500,2,rate=1/0.04,win=1)
Zobs<-ud$Xobs; Tobs<-ud$Tobs; status<-ud$status;
nno<-ud$nno;  tau<-15;

###################################################################
# Fits model with one common baseline and Z independent of X
# \lambda_0(t) \exp( Z^T \beta)
###################################################################
out<-em.ncc(Zobs,Tobs,status,tau,nno,Nit=100,beta=c(-.2,.2),
betait=20,detail=0,emvar=1)
print(c(out$beta,diag(out$var)^.5)); # estimates related to Z
plot(out$baseline,type="l")          # baseline estimates

#### using numDeriv for second derivative
out<-em.ncc(Zobs,Tobs,status,tau,nno,Nit=100,beta=c(-.2,.2),
betait=20,detail=0,emvar=2)
print(c(out$beta,diag(out$var)^.5)); # estimates related to Z


####################################################################
## slow way of fitting the same model, when censoring not the same 
####################################################################
out1<-em.ncc(Zobs,Tobs,status,rep(tau,nno),nno,Nit=100,beta=c(-.2,.2),
betait=20,detail=0,emvar=1)
print(c(out1$beta,diag(out1$var)^.5));     # estimates related to Z
lines(out$baseline,lty=2  )                # baseline estimates
#out$p                                     # distribution of Z

###################################################################
# Fits model : \lambda_S(X)(t) \exp( Z^T \beta)
#  Z | S(X) is equivalent gives conditional distribution of Z | X
###################################################################
strat<-rbinom(nrow(Zobs),1,0.5)
stratnoobs<-rbinom(500-nrow(Zobs),1,0.5)

out.strat<-em.ncc(Zobs,Tobs,status,tau,nno,Nit=100,beta=0,betait=20,
detail=0,stratum=c(strat,stratnoobs),emvar=1,
baseline.stratum=c(strat,stratnoobs))

print(cbind(out.strat$beta,diag(out.strat$var)^.5));
plot(out.strat$baseline[,1:2],type="l")     # baseline estimates
lines(out.strat$baseline[,c(1,3)],type="l") # baseline estimates
out.strat$p                                 # distribution of Z | S(X)
}
\keyword{survival}
