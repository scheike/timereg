\name{ClaytonOakes}
\alias{ClaytonOakes}
\title{Clayton-Oakes model with piece-wise constant hazards}
\usage{
  ClaytonOakes(formula, data = parent.frame(), id,
    var.formula = ~1, cuts = NULL, type = "co", start,
    control = list(), ...)
}
\arguments{
  \item{formula}{}

  \item{data}{}

  \item{id}{}

  \item{var.formula}{}

  \item{cuts}{}

  \item{type}{}

  \item{start}{}

  \item{control}{}

  \item{...}{}
}
\description{
  .. content for description (no empty lines) ..
}
\details{
  .. content for details ..
}
\examples{
d <- simClaytonOakes(2000,4,2,1,stoptime=2,left=0.5)
e <- ClaytonOakes(Surv(lefttime,time,status)~x1+id(~1,cluster),cuts=c(0,0.5,1,2),data=subset(d,!truncated))
e
plot(e,add=FALSE)
}
\author{
  Klaus K. Holst
}

