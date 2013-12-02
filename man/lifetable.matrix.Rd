\name{lifetable.matrix}
\alias{lifetable}
\alias{lifetable.formula}
\alias{lifetable.matrix}
\title{Life table}
\usage{
  \method{lifetable}{matrix}(x, strata = list(), breaks =
    c(), confint = FALSE, ...)

  \method{lifetable}{formula}(x, data=parent.frame(),
    breaks = c(), confint = FALSE, ...)
}
\arguments{
  \item{x}{time formula (Surv) or matrix/data.frame with
  columns time,status or entry,exit,status}

  \item{strata}{Strata}

  \item{data}{data.frame}

  \item{breaks}{Time intervals}

  \item{confint}{If TRUE 95% confidence limits are
  calculated}

  \item{...}{Additional arguments to lower level functions}
}
\description{
  Create simple life table
}
\examples{
library(timereg)
data(TRACE)
\donttest{
    lifetable(Surv(time,status==9)~sex+I(cut(wmi,c(-Inf,1,1.5,Inf))),
              data=TRACE,breaks=c(0.2,0.5),confint=TRUE)
}

d <- with(TRACE,lifetable(Surv(time,status==9)~sex+vf,breaks=c(0.2,0.5)))
summary(glm(events ~ offset(log(atrisk))+interval*vf + sex*vf,
            data=d,poisson))
}
\author{
  Klaus K. Holst
}

