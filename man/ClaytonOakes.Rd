\name{ClaytonOakes}
\alias{ClaytonOakes}
\title{Clayton-Oakes model with piece-wise constant hazards}
\usage{
  ClaytonOakes(formula, data = parent.frame(), id,
    var.formula = ~1, cuts = NULL, type = "co", start,
    control = list(), ...)
}
\arguments{
  \item{formula}{formula specifying the marginal
  proportional (piecewise constant) hazard structure with
  the right-hand-side being a survival object (Surv)
  specifying the entry time (optional), the follow-up time,
  and event/censoring status at follow-up. The clustering
  can be specified using the special function \code{id}
  (see example below).}

  \item{data}{Data frame}

  \item{id}{Variable defining the clustering (if not given
  in the formula)}

  \item{var.formula}{Formula specifying the variance
  component structure (if not given via the id special
  function in the formula)}

  \item{cuts}{Cut points defining the piecewise constant
  hazard}

  \item{type}{Type of estimation (Clayton-Oakes or
  conditional frailty model)}

  \item{start}{Optional starting values}

  \item{control}{Control parameters to the optimization
  routine}

  \item{...}{Additional arguments}
}
\description{
  .. content for description (no empty lines) ..
}
\details{
  .. content for details ..
}
\examples{
d <- simClaytonOakes(2000,4,2,1,stoptime=2,left=0.5)
d$z <- d$x1<0.5
e <- ClaytonOakes(Surv(lefttime,time,status)~x1+id(~1+z,cluster),cuts=c(0,0.5,1,2),data=subset(d,!truncated))
e
plot(e,add=FALSE)
}
\author{
  Klaus K. Holst
}

