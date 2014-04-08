\name{blocksample}
\alias{blocksample}
\title{Block sampling}
\usage{
  blocksample(data, size, idvar = "id", replace = TRUE,
    ids, ...)
}
\arguments{
  \item{data}{Data frame}

  \item{idvar}{Column defining the clusters}

  \item{size}{Size of samples}

  \item{replace}{Logical indicating wether to sample with
  replacement}

  \item{\dots}{additional arguments to lower level
  functions}
}
\value{
  \code{data.frame}
}
\description{
  Sample blockwise from clustered data
}
\examples{
dd <- data.frame(x=rnorm(6), z=rnorm(6), id=c(4,4,10,10,5,5), v=rnorm(6))
blocksample(dd,idvar="id",ids=c(4,5))
blocksample(dd,size=10)
}
\author{
  Klaus K. Holst
}
\keyword{models}
\keyword{utilities}

