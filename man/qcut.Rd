\name{qcut}
\alias{qcut}
\title{Cut a variable}
\description{
Calls the cut function to cut variables on data frame. 
}
\usage{
qcut(x,cuts=4,breaks=NULL,...)
}
\arguments{

  \item{x}{ variable to cut }
  \item{cuts}{ number of groups, 4 gives quartiles }
  \item{breaks}{ can also give breaks }
  \item{...}{ other argument for cut function of R}
}
\author{Thomas Scheike}
\examples{
data(sTRACE)
gx <- qcut(sTRACE$age)
table(gx)
}
\keyword{survival}
