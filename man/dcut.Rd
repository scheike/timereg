\name{dcut}
\alias{dcut}
\alias{qcut}
\alias{drm}
\title{Cuts variables and data frame}
\description{
Calls the cut function to cut variables on data frame and attach them to data-frame that is returned. 
}
\usage{
dcut(x,data=NULL,cuts=4,breaks=NULL,...)
}
\arguments{
  \item{x}{name of variable, or fomula, or names of variables on data frame.}
  \item{data}{if x is formula or names for data frame then data frame is needed.}
  \item{cuts}{vector of number of groups, 4 is default and gives quartiles}
  \item{breaks}{possible breaks for cutting.}
}
\author{Thomas Scheike}
\examples{
###library(timereg)
data(sTRACE)
sTRACE$age2 <- sTRACE$age^2
sTRACE$age3 <- sTRACE$age^3
names(sTRACE)

mm <- dcut(~age+wmi,data=sTRACE)
head(mm)

mm <- dcut(c("age","wmi"),data=sTRACE)
head(mm)

mm <- dcut(c("age","wmi"),data=sTRACE,cuts=c(2,4))
head(mm)

gx <- dcut(sTRACE$age)
head(gx)

gx <- qcut(sTRACE$age)
head(gx)

## Removes all cuts  variables with these names wildcards
mm <- drm(c("*.2","*.4"),mm)
head(mm)

## wildcards, for age, age2, age4 and wmi 
mm <- dcut(c("a*","*m*"),data=mm)
head(mm)
}
\keyword{survival}
