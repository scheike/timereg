\name{bicomprisk}
\alias{bicomprisk}
\title{Estimation of concordance in bivariate competing risks data}
\usage{
  bicomprisk(formula, data, cause = c(1, 1), cens = 0,
    causes, indiv, strata = NULL, id, num, prodlim = FALSE,
    messages = TRUE, model, return.data = 0, uniform = 0,
    robust = 0, ...)
}
\arguments{
  \item{formula}{Formula with left-hand-side being a
  \code{Hist} object (see example below) and the
  left-hand-side specying the covariate structure}

  \item{data}{Data frame}

  \item{cause}{Causes (default (1,1)) for which to estimate
  the bivariate cumulative incidence}

  \item{cens}{The censoring code}

  \item{causes}{causes}

  \item{indiv}{indiv}

  \item{strata}{Strata}

  \item{id}{Clustering variable}

  \item{num}{num}

  \item{prodlim}{prodlim}

  \item{messages}{Control amount of output}

  \item{model}{Type of competing risk model}

  \item{return.data}{Should data be returned (skipping
  modeling)}

  \item{uniform}{uniform}

  \item{robust}{robust}

  \item{...}{Additional arguments to lower level functions}
}
\description{
  Estimation of concordance in bivariate competing risks
  data
}
\author{
  Thomas Scheike, Klaus K. Holst
}

