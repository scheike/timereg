\name{ClaytonOakes}
\alias{ClaytonOakes}
\title{Clayton-Oakes model with piece-wise constant hazards}
\usage{
ClaytonOakes(formula, data = parent.frame(), cluster, var.formula = ~1,
  cuts = NULL, type = "piecewise", start, control = list(),
  var.invlink = exp, ...)
}
\arguments{
  \item{formula}{formula specifying the marginal
  proportional (piecewise constant) hazard structure with
  the right-hand-side being a survival object (Surv)
  specifying the entry time (optional), the follow-up time,
  and event/censoring status at follow-up. The clustering
  can be specified using the special function
  \code{cluster} (see example below).}

  \item{data}{Data frame}

  \item{cluster}{Variable defining the clustering (if not
  given in the formula)}

  \item{var.formula}{Formula specifying the variance
  component structure (if not given via the cluster special
  function in the formula) using a linear model with
  log-link.}

  \item{cuts}{Cut points defining the piecewise constant
  hazard}

  \item{type}{when equal to \code{two.stage}, the
  Clayton-Oakes-Glidden estimator will be calculated via
  the \code{timereg} package}

  \item{start}{Optional starting values}

  \item{control}{Control parameters to the optimization
  routine}

  \item{var.invlink}{Inverse link function for variance
  structure model}

  \item{...}{Additional arguments}
}
\description{
Clayton-Oakes frailty model
}
\examples{
set.seed(1)
d <- subset(simClaytonOakes(500,4,2,1,stoptime=2,left=2),!truncated)
e <- ClaytonOakes(Surv(lefttime,time,status)~x1+cluster(~1,cluster),
                  cuts=c(0,0.5,1,2),data=d)
e


d2 <- simClaytonOakes(500,4,2,1,stoptime=2,left=0)
d2$z <- rep(1,nrow(d2)); d2$z[d2$cluster\%in\%sample(d2$cluster,100)] <- 0
## Marginal=Cox Proportional Hazards model:
ts <- ClaytonOakes(Surv(time,status)~prop(x1)+cluster(~1,cluster),
                   data=d2,type="two.stage")
## Marginal=Aalens additive model:
ts2 <- ClaytonOakes(Surv(time,status)~x1+cluster(~1,cluster),
                    data=d2,type="two.stage")
## Marginal=Piecewise constant:
e2 <- ClaytonOakes(Surv(time,status)~x1+cluster(~-1+factor(z),cluster),
                   cuts=c(0,0.5,1,2),data=d2)
e2
plot(ts)
plot(e2,add=TRUE)

e3 <- ClaytonOakes(Surv(time,status)~x1+cluster(~1,cluster),cuts=c(0,0.5,1,2),
                   data=d,var.invlink=identity)
e3
}
\author{
Klaus K. Holst
}

