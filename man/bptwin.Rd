\name{bptwin}
\alias{bptwin}
\title{Liability model for twin data}
\usage{
  bptwin(formula, data, id, zyg, DZ, DZos, weight = NULL,
    biweight = function(x) 1/min(x), strata = NULL,
    messages = 1, control = list(trace = 0), type = "ace",
    eqmean = TRUE, param = 0, pairsonly = FALSE,
    stderr = TRUE, robustvar = TRUE, p, indiv = FALSE,
    constrain, samecens = TRUE,
    allmarg = samecens & !is.null(weight), bound = FALSE,
    debug = FALSE, ...)
}
\arguments{
  \item{formula}{}

  \item{data}{}

  \item{id}{}

  \item{zyg}{}

  \item{DZ}{}

  \item{DZos}{}

  \item{weight}{}

  \item{biweight}{}

  \item{strata}{}

  \item{messages}{}

  \item{control}{}

  \item{type}{}

  \item{eqmean}{}

  \item{param}{}

  \item{pairsonly}{}

  \item{stderr}{}

  \item{robustvar}{}

  \item{p}{}

  \item{indiv}{}

  \item{constrain}{}

  \item{samecens}{}

  \item{allmarg}{}

  \item{bound}{}

  \item{debug}{}

  \item{...}{}
}
\description{
  .. content for \description{} (no empty lines) ..
}
\details{
  .. content for \details{} ..
}
\author{
  Klaus Kähler Holst
}

