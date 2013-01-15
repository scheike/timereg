\name{fast.reshape}
\alias{fast.reshape}
\title{Fast reshape}
\usage{
  fast.reshape(data, id, varying, num, sep = "", ...)
}
\arguments{
  \item{data}{data.frame or matrix}

  \item{id}{id-variable. If omitted then reshape from wide
  to long.}

  \item{varying}{Vector of prefix-names of the time varying
  variables. Optional for long->wide reshaping.}

  \item{num}{Optional number/time variable}

  \item{sep}{String seperating prefix-name with
  number/time}

  \item{...}{Optional additional arguments to the
  \code{reshape} function used in the wide->long reshape.}
}
\description{
  Simple reshape/tranpose of data
}
\examples{
m <- lvm(c(y1,y2,y3,y4)~x)
d <- sim(m,1e1)

dd <- fast.reshape(d,var="y")
d1 <- fast.reshape(dd,"id")

## From wide-format
d1 <- fast.reshape(dd,"id")
d2 <- fast.reshape(dd,"id",var="y")
d3 <- fast.reshape(dd,"id",var="y",num="time")

## From long-format
fast.reshape(d,var="y",idvar="a",timevar="b")
fast.reshape(d,var=list(c("y1","y2","y3","y4")),idvar="a",timevar="b")
}
\author{
  Thomas Scheike, Klaus K. Holst
}

