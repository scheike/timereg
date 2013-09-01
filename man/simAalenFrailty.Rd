\name{simAalenFrailty}
\alias{simAalenFrailty}
\title{Simulate from the Aalen Frailty model}
\usage{
  simAalenFrailty(n = 5000, theta = 0.3, K = 2,
    beta0 = 1.5, beta = 1, cens = 1.5, ...)
}
\arguments{
  \item{K}{Number of clusters}

  \item{n}{Number of observations in each cluster}

  \item{eta}{1/variance}

  \item{beta}{Effect (log hazard ratio) of covariate}

  \item{stoptime}{Stopping time}

  \item{left}{Left truncation}

  \item{pairleft}{pairwise (1) left truncation or
  individual (0)}
}
\description{
  Simulate observations from Aalen Frailty model with Gamma
  distributed frailty and constant intensity.
}
\author{
  Klaus K. Holst
}

