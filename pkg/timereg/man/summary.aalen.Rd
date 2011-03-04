\name{summary.aalen}
\alias{summary.aalen}
\alias{summary.cox.aalen}
\alias{summary.prop.excess}
\alias{summary.timecox}
\alias{summary.dynreg}
\title{Prints summary statistics}
\description{
Computes p-values for test of significance for nonparametric
terms of model, p-values for test of constant effects based on
both supremum and integrated squared difference.

Returns parameter estimates and their standard errors. 
}
\usage{
\method{summary}{aalen}(object,digits=3,...)
}
\arguments{
\item{object}{an aalen object.}
\item{digits}{number of digits in printouts.}
\item{...}{unused arguments - for S3 compatibility}
}
\references{
Martinussen and Scheike, 
}
\author{Thomas Scheike}
\examples{
### see help(aalen)
}
\keyword{survival}
