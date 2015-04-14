\name{prep.comp.risk}
\alias{prep.comp.risk}
\title{Set up weights for delayed-entry competing risks data for comp.risk function}
\description{
Computes the weights of Geskus (2011) modified to the setting of the comp.risk function
}
\usage{
prep.comp.risk(data,times=NULL,entrytime=NULL,time="time",cause="cause",
	       strata=NULL,nocens.out=TRUE)
}
\arguments{
\item{data}{ data frame for comp.risk.}
\item{times}{ times for estimating equations.}
\item{entrytime}{name of delayed entry variable, if not given computes right-censoring case.}
\item{time}{name of survival time variable.}
\item{cause}{name of cause indicator}
\item{strata}{strata variable to obtain stratified weights.}
\item{nocens.out}{returns only uncensored part of data-frame}
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
}
\author{Thomas Scheike}
\examples{
data(bmt)
nn <- nrow(bmt)
entrytime <- rbinom(nn,1,0.5)*(bmt$time*runif(nn))
bmt$entrytime <- entrytime
times <- seq(5,70,by=1)

bmtw <- prep.comp.risk(times,data=bmt,time="time",entrytime="entrytime",cause="cause")

out <- comp.risk(Event(time,cause)~const(tcell)
	 +const(platelet),data=bmtw,cause=1,
 cens.weight=bmtw$cw,weights=bmtw$weights,times=times,n.sim=0)
summary(out)

outn <- comp.risk(Event(time,cause)~const(tcell)
     +const(platelet),data=bmt,cause=1,times=times,n.sim=0)
summary(outn)
}
\keyword{survival}
