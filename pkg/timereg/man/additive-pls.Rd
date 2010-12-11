\name{additive.pls}
\alias{additive.pls}
\alias{additive.plsR}
\alias{pls.surv.cv}
\title{Fits PLS for additive hazards model}
\description{
Fits the partial least squares estimator for the additive risk
model.
Time dependent variables and counting process data 
(multiple events per subject) are possible.

The modelling formula uses the standard survival modelling given in the 
\bold{survival} package.

Covariates \eqn{Z_1,...,Z_p}, fixed covariates \eqn{F(t)}. 
Algorithm : 
1) 
For pls components \eqn{X_1,...,X_K},  fits 
\deqn{
\lambda_0(t) + \alpha^T(t) F(t) + 
\sum_{j=1}^K X_j(t) \gamma_j(t) + Z_i(t) \beta_i^{K+1}
}
for \eqn{i=1,...,p}
2) compute new pls components \eqn{X_{K+1} = \sum \beta_i^{K+1} Z_i(t)}
   and iterate. 
}
\usage{
additive.pls(formula = formula(data), data = sys.parent(), 
start.time=0,max.time=NULL,id=NULL, pls.dim=1, 
silent=1)
}
\arguments{
\item{formula}{
a formula object with the response on the left of a '~' 
operator, and the independent terms on the right as regressors. 
The response must be a 
survival object as returned by the `Surv' function.

The const terms are kept as fixed covariates that are not involved in the 
pls variable reduction. This may be covariates that are know to be of clinical 
importance. 
}
\item{data}{a data.frame with the variables.}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{end of observation period where estimates are computed. 
Estimates thus computed from [start.time, max.time]. Default is max of data.}
\item{id}{For timevarying covariates the variable must 
associate each record with the id of a subject.}
\item{pls.dim}{number of pls components}
\item{silent}{set to 1 to avoid printing of warnings for non-inverible
design-matrices for different timepoints, default is 0.}
}
\details{
const() specifies the \eqn{F(t)} covariates in the above model. 
}
\value{
returns an object with the following arguments:
\item{baseline}{baseline of the semparametric additive risk
model}
\item{pls.comp}{the pls components, i.e. the covariates multipled on the beta
coefficients.}
\item{beta}{risk regression coefficients related to the pls components.}
\item{beta.pls}{regression coefficients that defines the pls componets.}
\item{tbeta.pls}{the combined regression coefficients from the pls components and the regression coefficients, these leads to risk predictions when applied to
new covariates.}
}
\references{
Martinussen and Scheike, Dynamic Regression Models for
Survival Data, Springer (2006).

Martinussen and Scheike, The Aalen additive hazards model
with high-dimensional regressors, (2009), Lifetime Data Anal.
}
\author{Thomas Scheike}
\examples{
data(mypbc)
pbc<-mypbc
pbc$time<-pbc$time+runif(418)*0.1; pbc$time<-pbc$time/365
pbc<-subset(pbc,complete.cases(pbc));
covs<-as.matrix(pbc[,-c(1:3,6)])
covs<-cbind(covs[,c(1:6,16)],log(covs[,7:15]))
covs<-scale(covs);

### 5 PLS components for the 16 covariates 
### based on Gui-Li approach with constant effects
out<-additive.pls(Surv(time,status>=1)~covs,pbc,max.time=7,pls.dim=5)

### survival predictions
par(mfrow=c(1,3))
pout=predict.pls(out,Z=covs)
matplot(pout$times,t(pout$surv[1:20,]),type="l",ylim=c(0,1),xlab="time",ylab="survival")

### cross-validated number of pls components

outcv<-pls.surv.cv(Surv(time,status>=1)~covs,pbc,max.time=7,pls.dims=0:2)
plot(outcv$pls.dims,outcv$cv,type="l")

out2<-additive.pls(Surv(time,status>=1)~covs,pbc,max.time=7,
pls.dim=outcv$cv.dim)

pout=predict.pls(out2, Z= covs)
matplot(pout$times,t(pout$surv[1:20,]),type="l",ylim=c(0,1),
xlab="time",ylab="survival")

### 2 PLS components for the 3 covariates 
### age and sex are fixed covariates 
data(sTRACE)
out<-additive.pls(Surv(time,status==9)~const(age)+const(sex)+
vf+diabetes+chf,sTRACE,max.time=7,pls.dim=2,silent=1)
}
\keyword{survival}
