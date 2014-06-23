\name{twinsim}
\alias{twinsim}
\title{Simulate twin data}
\usage{
twinsim(nMZ = 100, nDZ = nMZ, b1 = c(), b2 = c(), mu = 0,
  acde = c(1, 1, 0, 1), randomslope = NULL, threshold = 0, cens = FALSE,
  wide = FALSE, ...)
}
\arguments{
  \item{nMZ}{Number of monozygotic twin pairs}

  \item{nDZ}{Number of dizygotic twin pairs}

  \item{b1}{Effect of covariates (labelled x1,x2,...) of
  type 1. One distinct covariate value for each
  twin/individual.}

  \item{b2}{Effect of covariates (labelled g1,g2,...) of
  type 2. One covariate value for each twin pair.}

  \item{mu}{Intercept parameter.}

  \item{acde}{Variance of random effects (in the order
  A,C,D,E)}

  \item{randomslope}{Logical indicating wether to include
  random slopes of the variance components w.r.t.
  x1,x2,...}

  \item{threshold}{Treshold used to define binary outcome
  y0}

  \item{cens}{Logical variable indicating whether to censor
  outcome}

  \item{wide}{Logical indicating if wide data format should
  be returned}

  \item{...}{Additional arguments parsed on to lower-level
  functions}
}
\description{
Simulate twin data from a linear normal ACE/ADE/AE model.
}
\author{
Klaus K. Holst
}
\seealso{
\code{\link{twinlm}}
}
\keyword{models}
\keyword{regression}

