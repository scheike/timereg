\name{fast.reshape}
\alias{fast.reshape}
\title{Fast reshape
Simple reshape/tranpose of data}
\usage{
  fast.reshape(data, id, varying, num, sep = "", keep,
    idname = "id", numname = "num", factors.keep = TRUE,
    ...)
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

  \item{keep}{Vector of column names to keep}

  \item{idname}{Name of id-variable (Wide->Long)}

  \item{numname}{Name of number-variable (Wide->Long)}

  \item{factors.keep}{If false all factors are converted to
  integers}

  \item{...}{Optional additional arguments}
}
\description{
  Fast reshape Simple reshape/tranpose of data
}
\examples{
m <- lvm(c(y1,y2,y3,y4)~x)
d <- sim(m,5)
fast.reshape(fast.reshape(d,var="y"),id="id")

##### From wide-format
(dd <- fast.reshape(d,var="y"))
## Same with explicit setting new id and number variable/column names and seperator "" (default) and dropping x
fast.reshape(d,var="y",idname="a",timevar="b",sep="",keep=c())
## Same with 'reshape' list-syntax
fast.reshape(d,var=list(c("y1","y2","y3","y4")))

##### From long-format
fast.reshape(dd,id="id")
## Restrict set up within-cluster varying variables
fast.reshape(dd,id="id",var="y")
fast.reshape(dd,id="id",var="y",keep="x",sep=".")

#####
x <- data.frame(id=c(5,5,6,6,7),y=1:5,x=1:5,tv=c(1,2,2,1,2))
(xw <- fast.reshape(x,id="id"))
(xl <- fast.reshape(xw,varying=c("y","x"),idname="id2",keep=c()))
(xl <- fast.reshape(xw,varying=c("y","x","tv")))
(xw2 <- fast.reshape(xl,id="id",num="num"))
fast.reshape(xw2,varying=c("y","x"),idname="id")


##### Family cluster example
d <- mets:::sim.bin.fam(3)
d
dd <- fast.reshape(d,var="y")
dd

##### Prostate cancer example
data(prt)
head(prtw <- fast.reshape(prt,"id",var="cancer"))
ftable(cancer1~cancer2,data=prtw)
rm(prtw)
}
\author{
  Thomas Scheike, Klaus K. Holst
}

