\name{fast.reshape}
\alias{fast.reshape}
\title{Fast reshape}
\usage{
  fast.reshape(data, id, varying, num, sep = "",
    drop = NULL, idname = "id", numname = "num", ...)
}
\arguments{
  \item{data}{data.frame or matrix}

  \item{id}{id-variable. If omitted then reshape
  Wide->Long.}

  \item{varying}{Vector of prefix-names of the time varying
  variables. Optional for Long->Wide reshaping.}

  \item{num}{Optional number/time variable}

  \item{sep}{String seperating prefix-name with
  number/time}

  \item{drop}{Vector of column names to drop}

  \item{idname}{Name of id-variable (Wide->Long)}

  \item{numname}{Name of number-variable (Wide->Long)}

  \item{...}{Optional additional arguments}
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
d3 <- fast.reshape(dd,"id",var="y",num="num")

d4 <- fast.reshape(data.matrix(dd),"id",var="y")

## From long-format
fast.reshape(d,var="y",idvar="a",timevar="b")
fast.reshape(d,var=list(c("y1","y2","y3","y4")),idvar="a",timevar="b")

data(prt)
head(fast.reshape(prt,"id",var="cancer"))
}
\author{
  Thomas Scheike, Klaus K. Holst
}

