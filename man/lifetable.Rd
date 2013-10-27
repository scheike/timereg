\name{lifetable}
\alias{lifetable}
\title{Life table}
\usage{
  lifetable(time, status, entry, strata = list(),
    breaks = c(), confint = FALSE)
}
\arguments{
  \item{time}{Time (or matrix with columns time,status or
  entry,time,status)}

  \item{status}{status (event=TRUE,censoring=FALSE)}

  \item{entry}{Entry time}

  \item{strata}{Strata}

  \item{breaks}{Time intervals}

  \item{confint}{If TRUE 95% confidence limits are
  calculated}
}
\description{
  Creat simple life table
}
\examples{
library(timereg)
data(TRACE)
with(TRACE, lifetable(cbind(time,status==9),breaks=c(0.2,0.5),confint=TRUE))

d <- with(TRACE, lifetable(time,status==9,strata=list(sex=sex,vf=vf),breaks=c(0.2,0.5)))
summary(glm(events ~ offset(log(atrisk))+interval+sex+vf,data=d,poisson))
}
\author{
  Klaus K. Holst
}

