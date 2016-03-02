\name{event.split}
\alias{event.split}
\title{EventSplit (SurvSplit).}
\description{ contstructes start stop formulation of event time data after
a variable in the data.set, otherwise similar to SurvSplit of the survival package.
Can also split after random time given in data frame.
}
\usage{
event.split(data, time="time",status="status",cuts="cuts",
	     name.start="start",id=NULL,name.id="id",
			entry=NULL,cens.code=0,num="num")
}
\arguments{
  \item{data}{data to be split }
  \item{time}{time variable.}
  \item{status}{status variable.}
  \item{cuts}{contrast matrix for testing}
  \item{name.start}{name of start variable in new data.}
  \item{id}{name of id variable.}
  \item{name.id}{name of id variable in new data.}
  \item{entry}{entry time for delayed entry.}
  \item{cens.code}{code for the censoring.}
  \item{num}{numbering of events to same subject.}
}
\author{Thomas Scheike}
\examples{
set.seed(1)
d <- data.frame(event=round(5*runif(5),2),
		start=round(5*runif(5),2),
		time=1:5,
		status=rbinom(5,1,0.5),
		x=1:5)
d$fcut <- 2
d
event.split(d,cuts="event",num="nume")
### like survSplit
event.split(d,cuts="fcut",num="numf")

### successive splitting 
de <- event.split(d,cuts="event",num="nume")
de
def <- event.split(de,cuts="fcut",entry="start",id="id",num="numf")
###
### successive splitting 
df <- event.split(d,cuts="fcut",num="numf")
df
dfe <- event.split(df,cuts="event",entry="start",id="id",num="nume")
dfe <- dfe[,names(def)]
dfe
def
dfe==def
}
\keyword{survival}
