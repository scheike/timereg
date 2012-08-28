\name{probandwise}
\alias{probandwise}
\title{Estimates the probandwise concordance based on Concordance and marginal estimate}
\usage{
  probandwise(conc, marg, cause.prodlim = 1)
}
\arguments{
  \item{conc}{Concordance}

  \item{marg}{Marginal estimate}

  \item{cause.prodlim}{specififes which cause that should
  be used for marginal cif based on prodlim}
}
\description{
  .. content for description (no empty lines) ..
}
\examples{
data(prt);

### marginal cumulative incidence of prostate cancer
outm <- prodlim(Hist(time,status)~+1,data=prt)

times <- 60:100
cifmz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="MZ"))
cifdz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="DZ"))

### concordance for MZ and DZ twins
cc <- bicomprisk(Hist(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2),prodlim=TRUE)
cdz <- cc$model$"DZ"
cmz <- cc$model$"MZ"

cdz <- probandwise(cdz,outm,cause.prodlim=2)
cmz <- probandwise(cmz,outm,cause.prodlim=2)

plot(cmz,ci=NULL,ylim=c(0,0.5),xlim=c(60,100),legend=TRUE,col=c(3,2,1))
par(new=TRUE)
plot(cdz,ci=NULL,ylim=c(0,0.5),xlim=c(60,100),legend=TRUE)
summary(cdz)
summary(cmz)
}
\author{
  Thomas Scheike
}

