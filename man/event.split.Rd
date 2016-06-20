\name{event.split}
\alias{event.split}
\title{EventSplit (SurvSplit).}
\description{ contstructs start stop formulation of event time data after
a variable in the data.set. Similar to SurvSplit of the survival package but
can also split after random time given in data frame.
}
\usage{
event.split(data, time="time",status="status",cuts="cuts",
                  name.id="id", name.start="start", cens.code=0, order.id=TRUE,
		  time.group=TRUE)
}
\arguments{
  \item{data}{data to be split }
  \item{time}{time variable.}
  \item{status}{status variable.}
  \item{cuts}{cuts variable or numeric cut (only one value)}
  \item{name.start}{name of start variable in data, start can also be numeric "0"}
  \item{name.id}{name of id variable.}
  \item{cens.code}{code for the censoring.}
  \item{order.id}{order data after id and start.}
  \item{time.group}{make variable "before"."cut" that keeps track of wether start,stop is before (1) or after cut (0).}
}
\author{Thomas Scheike}
\examples{
set.seed(1)
d <- data.frame(event=round(5*runif(5),2),start=1:5,time=2*1:5,
		status=rbinom(5,1,0.5),x=1:5)
d

d0 <- event.split(d,cuts="event",name.start=0)
d0

dd <- event.split(d,cuts="event")
dd
ddd <- event.split(dd,cuts=3.5)
ddd
event.split(ddd,cuts=5.5)

### successive cutting for many values 
dd <- d
for  (cuts in seq(2,3,by=0.3)) dd <- event.split(dd,cuts=cuts)
dd

}
\keyword{survival}
