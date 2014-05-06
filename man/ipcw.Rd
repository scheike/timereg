\name{ipcw}
\alias{ipcw}
\title{Inverse Probability of Censoring Weights}
\usage{
  ipcw(formula, data, cluster, samecens = FALSE,
    obsonly = TRUE, weightname = "w", iid = TRUE,
    cens.model = "aalen", pairs = FALSE, ...)
}
\arguments{
  \item{formula}{Formula specifying the censoring model}

  \item{data}{data frame}

  \item{cluster}{clustering variable}

  \item{samecens}{For clustered data, should same censoring
  be assumed (bivariate probability calculated as mininum
  of the marginal probabilities)}

  \item{obsonly}{Return data with uncensored observations
  only}

  \item{weightname}{Name of weight variable in the new
  data.frame}

  \item{cens.model}{Censoring model (default Aalens
  additive model)}

  \item{pairs}{For paired data (e.g. twins) only the
  complete pairs are returned (With pairs=TRUE)}

  \item{...}{Additional arguments to censoring model}
}
\description{
  Calculates Inverse Probability of Censoring Weights
  (IPCW) and adds them to a data.frame
}
\examples{
data(prt)
prtw <- ipcw(Surv(time,status==0)~country, data=prt[sample(nrow(prt),5000),],
            cluster="id",weightname="w")
plot(0,type="n",xlim=range(prtw$time),ylim=c(0,1),xlab="Age",ylab="Probability")
count <- 0
for (l in unique(prtw$country)) {
    count <- count+1
    prtw <- prtw[order(prtw$time),]
    with(subset(prtw,country==l),
         lines(time,w,col=count,lwd=2))
}
legend("topright",legend=unique(prtw$country),col=1:4,pch=-1,lty=1)
}
\author{
  Klaus K. Holst
}

