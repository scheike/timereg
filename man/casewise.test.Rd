\name{casewise.test}
\alias{casewise.test}
\title{Estimates the casewise concordance based on Concordance and marginal estimate using timereg and performs test for independence}
\usage{
  casewise.test(conc, marg, test = "no-test")
}
\arguments{
  \item{conc}{Concordance}

  \item{marg}{Marginal estimate}

  \item{test}{Type of test for independence assumption.
  "conc" makes test on concordance scale and "case" means a
  test on the casewise concordance}
}
\description{
  .. content for description (no empty lines) ..
}
\details{
  Uses cluster based conservative standard errors for
  marginal
}
\examples{
data(prt);

prt <- prt[which(prt$id \%in\% sample(unique(prt$id),10000)),]
### marginal cumulative incidence of prostate cancer
times <- seq(60,100,by=2)
outm <- comp.risk(Surv(time,status==0)~+1,data=prt,prt$status,causeS=2,times=times)

cifmz <- predict(outm,X=1,uniform=0,resample.iid=1)
cifdz <- predict(outm,X=1,uniform=0,resample.iid=1)

### concordance for MZ and DZ twins
cc <- bicomprisk(Hist(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2))
cdz <- cc$model$"DZ"
cmz <- cc$model$"MZ"

### To compute casewise cluster argument must be passed on, here with a max of 100 to limit comp-time
outm <-comp.risk(Surv(time,status==0)~+1,data=prt,prt$status,causeS=2,times=times,max.clust=100)
cifmz <-predict(outm,X=1,uniform=0,resample.iid=1)
cc <-bicomprisk(Hist(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2),se.clusters=outm$clusters)
cdz <- cc$model$"DZ"
cmz <- cc$model$"MZ"

cdz <- casewise.test(cdz,cifmz,test="case") ## test based on casewise
cmz <- casewise.test(cmz,cifmz,test="conc") ## test based on concordance:w

plot(cmz,ylim=c(0,0.7),xlim=c(60,100))
par(new=TRUE)
plot(cdz,ylim=c(0,0.7),xlim=c(60,100))
}
\author{
  Thomas Scheike
}

