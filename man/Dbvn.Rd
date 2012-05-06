\name{Dbvn}
\alias{Dbvn}
\title{Derivatives of the bivariate normal cumulative distribution function}
\usage{
  Dbvn(p, design = function(p, ...) {
        return(list(mu = cbind(p[1], p[1]), dmu = cbind(1, 1), S = matrix(c(p[2], p[3], p[3], p[4]), ncol = 2), dS = rbind(c(1, 0, 0, 0), c(0, 1, 1, 0), c(0, 0, 0, 1))))
    }, Y = cbind(0, 0))
}
\arguments{
  \item{p}{}

  \item{design}{}

  \item{Y}{}
}
\description{
  .. content for description (no empty lines) ..
}
\details{
  .. content for details ..
}
\author{
  Klaus K. Holst
}

