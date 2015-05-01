\name{prep.comp.risk}
\alias{prep.comp.risk}
\title{Set up weights for delayed-entry competing risks data for comp.risk function}
\description{
Computes the weights of Geskus (2011) modified to the setting of the comp.risk function.
The returned weights are 
\eqn{1/(H(min(T_i,tau))*G_c(min(T_i,tau)))} and tau is the max of the times argument, 
here \eqn{H} is the estimator of the truncation distribution and \eqn{G_c} is the right censoring distribution.
}
\usage{
prep.comp.risk(data,times=NULL,entrytime=NULL,time="time",cause="cause",
	  	    strata=NULL,nocens.out=TRUE,cens.formula=NULL)
}
\arguments{
\item{data}{ data frame for comp.risk.}
\item{times}{ times for estimating equations.}
\item{entrytime}{name of delayed entry variable, if not given computes right-censoring 
case.}
\item{time}{name of survival time variable.}
\item{cause}{name of cause indicator}
\item{strata}{strata variable to obtain stratified weights.}
\item{nocens.out}{returns only uncensored part of data-frame}
\item{cens.formula}{censoring model formula for Cox models for the truncation and censoring 
model.}
}
\value{
Returns an object. With the following arguments:
\item{dataw}{a data.frame with weights.}

The function wants to make two new variables "weights" and "cw" so if these already are in
the data frame it tries to add an "_" in the names. 
}
\references{
Geskus (2011), Cause-Specific Cumulative Incidence Estimation and the Fine and Gray Model Under
Both Left Truncation and Right Censoring, Biometrics (2011), pp 39-49.  

Shen (2011), Proportional subdistribution hazards regression for left-truncated competing 
risks data, Journal of Nonparametric Statistics (2011), 23, 885-895
}
\author{Thomas Scheike}
\examples{
data(bmt)
nn <- nrow(bmt)
entrytime <- rbinom(nn,1,0.5)*(bmt$time*runif(nn))
bmt$entrytime <- entrytime
times <- seq(5,70,by=1)

bmtw <- prep.comp.risk(times,data=bmt,time="time",entrytime="entrytime",cause="cause")

out <- comp.risk(Event(time,cause)~const(tcell)+const(platelet),data=bmtw,
		 cause=1,cens.weight=bmtw$cw,weights=bmtw$weights,times=times,n.sim=0)
summary(out)

outn <- comp.risk(Event(time,cause)~const(tcell)+const(platelet),data=bmt,
		  cause=1,times=times,n.sim=0)
summary(outn)

out <- comp.risk(Event(time,cause)~+1,data=bmtw,cause=1,model="rcif2",
		 cens.weight=bmtw$cw,weights=bmtw$weights,times=c(5,30,70),n.sim=0)
out$cum
### same as
###out <- prodlim(Hist(entry=entrytime,time,cause)~+1,data=bmtw)
###summary(out,cause="1",times=c(5,30,70))

out <- comp.risk(Event(time,cause)~+1,data=bmtw,cause=1,model="rcif2",
		 times=c(5,30,70),n.sim=0)
out$cum
### same as
###out <- prodlim(Hist(time,cause)~+1,data=bmtw)
###summary(out,cause="1",times=c(5,30,70))




}
\keyword{survival}
