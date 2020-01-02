\name{bmt}
\alias{bmt}
\non_function{}
\title{The Bone Marrow Transplant Data}
\description{
  Bone marrow transplant data with 408 rows and 5 columns.
}
\format{
The data has 408 rows and 5 columns.
\describe{
\item{cause}{a numeric vector code. 
Survival status. 1: dead from treatment related causes, 
2: relapse , 0: censored.}
\item{time}{ a numeric vector. Survival time.  }
\item{platelet}{a numeric vector code. Plalelet 
1: more than 100 x \eqn{10^9} per L, 0: less.}
\item{tcell}{a numeric vector. T-cell depleted BMT 1:yes,
0:no.}
\item{age}{a numeric vector code. Age of patient, scaled and
centered ((age-35)/15).}
}
}
\source{
Simulated data  
}
\references{
NN 
}
\examples{
data(bmt)
names(bmt)
}
\keyword{datasets}
