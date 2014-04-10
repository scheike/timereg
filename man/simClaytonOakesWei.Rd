\name{simClaytonOakesWei}
\alias{simClaytonOakesWei}
\title{Simulate from the Clayton-Oakes frailty model}
\usage{
simClaytonOakesWei(K, n, eta, beta, stoptime, weiscale = 1, weishape = 2,
  left = 0, pairleft = 0)
}
\arguments{
  \item{K}{Number of clusters}

  \item{n}{Number of observations in each cluster}

  \item{eta}{1/variance}

  \item{beta}{Effect (log hazard ratio) of covariate}

  \item{stoptime}{Stopping time}

  \item{weiscale}{weibull scale parameter}

  \item{weishape}{weibull shape parameter}

  \item{left}{Left truncation}

  \item{pairleft}{pairwise (1) left truncation or
  individual (0)}
}
\description{
Simulate observations from the Clayton-Oakes copula model
with Weibull type baseline and Cox marginals.
}
\author{
Klaus K. Holst
}

