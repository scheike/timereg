\name{simAalenFrailty}
\alias{simAalenFrailty}
\title{Simulate from the Aalen Frailty model}
\usage{
simAalenFrailty(n = 5000, theta = 0.3, K = 2, beta0 = 1.5, beta = 1,
  cens = 1.5, cuts = 0, ...)
}
\arguments{
  \item{n}{Number of observations in each cluster}

  \item{theta}{Dependence paramter (variance of frailty)}

  \item{K}{Number of clusters}

  \item{beta0}{Baseline (intercept)}

  \item{beta}{Effect (log hazard ratio) of covariate}

  \item{cens}{Censoring rate}

  \item{cuts}{time cuts}

  \item{...}{Additional arguments}
}
\description{
Simulate observations from Aalen Frailty model with Gamma
distributed frailty and constant intensity.
}
\author{
Klaus K. Holst
}

