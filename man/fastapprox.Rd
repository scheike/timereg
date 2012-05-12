\name{fastapprox}
\alias{fastapprox}
\title{Fast approximation}
\usage{
  fastapprox(t1, t2, y1 = t1, ...)
}
\arguments{
  \item{t1}{Original ordered time points}

  \item{t2}{New time points}

  \item{y1}{Optional function values (default is simply the
  identity) corresponding to the original time points
  \code{t1}}

  \item{...}{Optional additional arguments}
}
\value{
  List with the elements

  \item{approx}{Approximated function values (\code{y1})
  evaluated in \code{t2}}

  \item{pos}{The positions in \code{t1} of the elements of
  \code{t2}}
}
\description{
  Fast approximation
}
\author{
  Klaus K. Holst
}

