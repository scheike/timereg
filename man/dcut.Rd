\name{qcut}
\title{Cut a variables}
\description{
Calls the cut function to cut variables on data frame and attach them to data-frame that is returned. 
}
\usage{
qcut(x,cuts=4,breaks=NULL,...)
}
\author{Thomas Scheike}
\examples{
###library(timereg)
data(sTRACE)

gx <- qcut(sTRACE$age)
head(gx)
}
\keyword{survival}
