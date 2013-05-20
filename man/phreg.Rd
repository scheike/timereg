\name{phreg}
\alias{phreg}
\title{Fast Cox PH regression}
\usage{
  phreg(formula, data, ...)
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
n <- 1e3;
d <- mets:::simCox(n); d$id <- seq(nrow(d)); d$group <- factor(rbinom(nrow(d),1,0.5))

(m1 <- phreg(Surv(entry,time,status)~X1+X2,data=d))
(m2 <- coxph(Surv(entry,time,status)~X1+X2+cluster(id),data=d))
(coef(m3 <-cox.aalen(Surv(entry,time,status)~prop(X1)+prop(X2),data=d)))

\dontrun{
(m1b <- phreg(Surv(entry,time,status)~X1+X2+strata(group),data=d))
(m2b <- coxph(Surv(entry,time,status)~X1+X2+cluster(id)+strata(group),data=d))
(coef(m3b <-cox.aalen(Surv(entry,time,status)~-1+group+prop(X1)+prop(X2),data=d)))
}

m <- phreg(Surv(entry,time,status)~X1*X2+strata(group)+cluster(id),data=d)
m
plot(m,ylim=c(0,5))
}
\author{
  Klaus K. Holst
}

