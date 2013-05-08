\name{fast.cox}
\alias{fast.cox}
\title{Fast Cox PH regression}
\usage{
  fast.cox(formula, data, ...)
}
\arguments{
  \item{formula}{formula with 'Surv' outcome (see
  \code{coxph})}

  \item{data}{data frame}

  \item{...}{Additional arguments to lower level funtions}
}
\description{
  Fast Cox PH regression
}
\examples{
simcox <- function(n=1000, seed=1, beta=c(1,1), entry=TRUE) {
  if (!is.null(seed))
    set.seed(seed)
  library(lava)
  m <- lvm()
  regression(m,T~X1+X2) <- beta
  distribution(m,~T+C) <- coxWeibull.lvm(scale=1/100)
  distribution(m,~entry) <- coxWeibull.lvm(scale=1/10)
  m <- eventTime(m,time~min(T,C=0),"status")
  d <- sim(m,n);
  if (!entry) d$entry <- 0
  else d <- subset(d, time>entry,select=-c(T,C))
  return(d)
}

n <- 1e3;
d <- simcox(n); d$id <- seq(nrow(d)); d$group <- factor(rbinom(nrow(d),1,0.5))
fast.cox(Surv(entry,time,status)~X1*X2+strata(group)+cluster(id),data=d)

(m1 <- fast.cox(Surv(entry,time,status)~X1+X2,data=d))
(m2 <- coxph(Surv(entry,time,status)~X1+X2+cluster(id),data=d))
(coef(m3 <-cox.aalen(Surv(entry,time,status)~prop(X1)+prop(X2),data=d)))

(m1b <- fast.cox(Surv(entry,time,status)~X1+X2+strata(group),data=d))
(m2b <- coxph(Surv(entry,time,status)~X1+X2+cluster(id)+strata(group),data=d))
(coef(m3b <-cox.aalen(Surv(entry,time,status)~-1+group+prop(X1)+prop(X2),data=d)))
}
\author{
  Klaus K. Holst
}

