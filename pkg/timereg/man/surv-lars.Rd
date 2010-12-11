\name{surv.lars}
\alias{mylars}
\alias{surv.lars}
\alias{surv.lars.cv}
\alias{mycoef.lars}
\alias{mypredict.lars}
\title{Fits LASSO model for  additive hazards model by Lars
algorithm}
\description{
Fits the LASSO estimator for the additive risk
model based on the least squares fitting criterion

\deqn{
L(\beta,S,s) = \beta^T S \beta - 2 \beta^T s
}
where \eqn{S=\int Z H Z dt} and  \eqn{s=\int Z H dN}.

This is equivalent to an appropriate normal least squares
problem on least squares data. 

Much quicker implementaion in ahaz package that is recommended. 
}
\usage{
surv.lars(S, s, n, l1.weights = NULL, ...)
}
\arguments{
\item{S}{the S matrix defined above.}
\item{s}{the s vector defined above.}
\item{n}{number of subjects, used in lars to decide dimension.}
\item{l1.weights}{NOT WORKING yet. specifies the weights for the L1 penalty.}
\item{...}{unused arguments - for S3 compatibility}
}
\details{
Modified version of standard lars program in the LARS package.
Essentially all scalings are removed from the lars program. 

The matrices S and S can be computed using the aalen() function. 
}
\value{
see lars. 
}
\references{
Martinussen and Scheike, Model selection for the the additive
risks hazards model,  Scandinavian Journal of Statistics, 2008, pages . 

Martinussen and Scheike, Dynamic Regression Models for
Survival Data, Springer (2006).

LARS, Efron et al. 
}
\author{Thomas Scheike}
\examples{
## makes data for pbc complete case
data(mypbc)
pbc<-mypbc
pbc$time<-pbc$time+runif(418)*0.1; pbc$time<-pbc$time/365
pbc<-subset(pbc,complete.cases(pbc)); 
covs<-as.matrix(pbc[,-c(1:3,6)])
covs<-cbind(covs[,c(1:6,16)],log(covs[,7:15]))
covs<-scale(covs); 

## computes the matrices needed for the least squares 
## criterion 

out<-aalen(Surv(time,status>=1)~const(covs),data=pbc,robust=0)
S<-out$intZHZ; s<-out$intZHdN
n<-nrow(pbc)

## lasso for survival data 
if (require(lars)) {
fit<-surv.lars(S,s,n)
plot(fit)

cv<-surv.lars.cv(Surv(time,status>=1)~ const(covs),data=pbc)
beta<-mypredict.lars(fit,cv$cv.frac,type="coefficients",mode = "fraction")$coef
c(beta)

## fitting survival model with these coefficients
out<-aalen(Surv(time,status>=1)~ const(covs),data=pbc,robust=0,gamma=beta)
pout<-predict(out,Z=covs[1:20,],uniform=0,se=0)
plot(pout,multiple=1,se=0,uniform=0,main="Lars survival predictions"); 

full<-fit$beta[nrow(fit$beta),] # least squares solution
c(full)
}
}
\keyword{survival}
