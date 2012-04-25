\name{biprobit}
\alias{biprobit}
\title{Bivariate Probit model}
\usage{
  biprobit(formula, data, id, time, strata = NULL,
    eqmarg = TRUE, indep = FALSE, weight = NULL,
    biweight = function(x) 1/min(x), samecens = TRUE,
    randomeffect = FALSE, vcov = "robust",
    pairsonly = FALSE,
    allmarg = samecens & !is.null(weight),
    control = list(trace = 0, method = "qausi"),
    bound = FALSE, messages = 1, p, ...)
}
\arguments{
  \item{formula}{}

  \item{data}{}

  \item{id}{}

  \item{time}{}

  \item{strata}{}

  \item{eqmarg}{}

  \item{indep}{}

  \item{weight}{}

  \item{biweight}{}

  \item{samecens}{}

  \item{randomeffect}{}

  \item{vcov}{}

  \item{pairsonly}{}

  \item{allmarg}{}

  \item{control}{}

  \item{bound}{}

  \item{messages}{}

  \item{p}{}

  \item{...}{}
}
\value{
  \code{biprobit}
}
\description{
  Bivariate Probit model
}
\details{
  .. content for \details{} ..
}
\author{
  Klaus K. Holst
}
\keyword{regression}

