\name{ipw}
\alias{ipw}
\title{Inverse Probability of Censoring Weights}
\usage{
  ipw(formula, data, cluster, samecens = FALSE,
    obsonly = TRUE, weightname = "w", pairs = FALSE)
}
\arguments{
  \item{formula}{Formula specifying the censoring model}

  \item{data}{data frame}

  \item{cluster}{clustering variable}

  \item{samecens}{For clustered data, should same censoring
  be assumed (bivariate probability calculated as mininum
  of the marginal probabilities)}

  \item{obsonly}{Return data with uncensored observations
  only}

  \item{weightname}{Name of weight variable in the new
  data.frame}

  \item{pairs}{For paired data (e.g. twins) only the
  complete pairs are returned (With pairs=TRUE)}
}
\description{
  Calculates Inverse Probability of Censoring Weights
  (IPCW) and adds them to a data.frame
}
\author{
  Klaus K. Holst
}

