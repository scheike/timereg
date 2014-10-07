\name{predict.timereg}
\alias{predict.timereg}
\alias{predict.aalen}
\alias{predict.comprisk}
\alias{predict.cox.aalen}
\title{Predictions for Survival and Competings Risks Regression for timereg} 
\description{
Make predictions based on the survival models 
(Aalen and Cox-Aalen) and the competing risks
models for the cumulative incidence function
(comp.risk). Computes confidence intervals and
confidence bands based on resampling. 
}
\usage{
\method{predict}{timereg}(object,newdata=NULL,X=NULL,times=NULL,Z=NULL,
n.sim=500,uniform=TRUE,se=TRUE,alpha=0.05,resample.iid=0,...)
}
\arguments{
  \item{object}{an object belonging to one of the following
    classes: comprisk, aalen or cox.aalen}
  \item{newdata}{specifies the data at which the
    predictions are wanted.}
  \item{X}{alternative to newdata, specifies the
    nonparametric components for predictions.}
  \item{Z}{alternative to newdata, specifies the
    parametric components of the model for predictions.}
  \item{times}{times in which predictions are computed, default is all time-points for baseline}
  \item{n.sim}{number of simulations in resampling.}
  \item{uniform}{computes resampling based uniform
    confidence bands.}
  \item{se}{computes pointwise standard errors}
  \item{alpha}{specificies the significance levelwhich cause we consider.}
  \item{resample.iid}{set to 1 to return iid decomposition of estimates, 3-dim matrix (predictions x times x subjects)}
  \item{...}{unused arguments - for S3 compatability}
}
\value{
\item{time}{vector of time points where the predictions are computed.}
\item{unif.band}{resampling based constant to construct 95\% uniform 
confidence bands.}
\item{model}{specifies what model that was fitted.}
\item{alpha}{specifies the significance level for the confidence intervals.
This relates directly to the constant given in unif.band.}
\item{newdata}{specifies the newdata given in the call.}
\item{RR}{gives relative risk terms for Cox-type models.}
\item{call}{gives call for predict funtion.}
\item{initial.call}{gives call for underlying object used for predictions.}
\item{P1}{gives cumulative inicidence predictions for competing risks models.
Predictions given in 
matrix form with different subjects in different rows.} 
\item{S0}{gives survival predictions for survival models. 
Predictions given in 
matrix form with different subjects in different rows.} 
\item{se.P1}{pointwise standard errors for predictions of P1.}
\item{se.S0}{pointwise standard errors for predictions of S0.}
}
\references{
Scheike, Zhang and Gerds (2008), Predicting cumulative
incidence probability by direct binomial regression,
Biometrika, 95, 205-220. 

Scheike and Zhang (2007), Flexible competing risks regression modelling and goodness of fit, LIDA, 14, 464-483 . 

Martinussen and Scheike (2006), Dynamic regression
models for survival data, Springer. 
}
\author{Thomas Scheike, Jeremy Silver}
\examples{
library(timereg)
data(bmt); 

add<-comp.risk(Event(time,cause)~platelet+age+tcell,data=bmt,cause=1)

ndata<-data.frame(platelet=c(1,0,0),age=c(0,1,0),tcell=c(0,0,1))
out<-predict(add,newdata=ndata,uniform=1,n.sim=1000)
par(mfrow=c(2,2))
plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=1)
# see comp.risk for further examples. 

add<-comp.risk(Event(time,cause)~factor(tcell),data=bmt,cause=1)
summary(add)
out<-predict(add,newdata=ndata,uniform=1,n.sim=1000)
plot(out,multiple=1,uniform=1,col=1:3,lty=1,se=1)

\donttest{
## SURVIVAL predictions aalen function
data(sTRACE)
out<-aalen(Surv(time,status==9)~sex+ diabetes+chf+vf,
data=sTRACE,max.time=7,n.sim=0,resample.iid=1)

pout<-predict(out,X=rbind(c(1,0,0,0,0),rep(1,5)))
head(pout$S0[,1:5]); head(pout$se.S0[,1:5])
par(mfrow=c(2,2))
plot(pout,multiple=1,se=0,uniform=0,col=1:2,lty=1:2)
plot(pout,multiple=0,se=1,uniform=1,col=1:2)

out<-aalen(Surv(time,status==9)~const(age)+const(sex)+
const(diabetes)+chf+vf,
data=sTRACE,max.time=7,n.sim=0,resample.iid=1)

pout<-predict(out,X=rbind(c(1,0,0),c(1,1,0)),
Z=rbind(c(55,0,1),c(60,1,1)))
head(pout$S0[,1:5]); head(pout$se.S0[,1:5])
par(mfrow=c(2,2))
plot(pout,multiple=1,se=0,uniform=0,col=1:2,lty=1:2)
plot(pout,multiple=0,se=1,uniform=1,col=1:2)

pout<-predict(out,uniform=0,se=0,newdata=sTRACE[1:10,]) 
plot(pout,multiple=1,se=0,uniform=0)

out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
prop(diabetes)+chf+vf,
data=sTRACE,max.time=7,n.sim=0,resample.iid=1)

pout<-predict(out,X=rbind(c(1,0,0),c(1,1,0)),Z=rbind(c(55,0,1),c(60,1,1)))
head(pout$S0[,1:5]); head(pout$se.S0[,1:5])
par(mfrow=c(2,2))
plot(pout,multiple=1,se=0,uniform=0,col=1:2,lty=1:2)
plot(pout,multiple=0,se=1,uniform=1,col=1:2)

pout<-predict(out,uniform=0,se=0,newdata=sTRACE[1:10,]) 
plot(pout,multiple=1,se=0,uniform=0)
}
}
\keyword{survival}
