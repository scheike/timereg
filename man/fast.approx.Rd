\name{fast.approx}
\alias{fast.approx}
\title{Fast approximation}
\usage{
  fast.approx(time, new.time, y = time, ...)
}
\arguments{
  \item{time}{Original ordered time points}

  \item{new.time}{New time points}

  \item{y}{Optional function values (default is simply the
  identity) evaluated at the original time points
  \code{time}}

  \item{...}{Optional additional arguments}
}
\value{
  List with the elements

  \item{approx}{Approximated function values (\code{y1})
  evaluated in \code{t2}}

  \item{pos}{The positions in \code{t1} of the elements of
  \code{t2} starting from index 0}
}
\description{
  Fast approximation
}
\author{
  Klaus K. Holst
}

