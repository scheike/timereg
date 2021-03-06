% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-pc-hazard.r
\name{invsubdist}
\alias{invsubdist}
\title{Finds inverse of piecwise linear sub-distribution}
\usage{
invsubdist(F1, u, entry = NULL, cond = 1, ptrunc = NULL)
}
\arguments{
\item{F1}{matrix with x, and F1(x)}

\item{u}{points for which to compute inverse}

\item{entry}{possible delayed entry points}

\item{cond}{1 indcates that we draw given that this subdistribution is used, so draws from F1(t)/F1(tau)}

\item{ptrunc}{possible trunction weigth for delayed entry, if NULL then uses ptrunc=1-F1(entry)}
}
\description{
Finds inverse of piecwise linear sub-distribution to be used for simulation of subdistributions
}
\examples{
F1 <- cbind(c(0,5,8,10),c(0,0.1,0.3,0.9))
plot(F1,type="l")
u <- runif(100)
Fiu <- invsubdist(F1,u,cond=0)
points(Fiu$time,u,pch="x")

F1cond <- F1
F1cond[,2] <- F1cond[,2]/0.9
plot(F1cond,type="l")
u <- runif(100)
Ficond <- invsubdist(F1cond,u,cond=0)
points(Ficond$time,u,pch="-")
Fiu <- invsubdist(F1,u,cond=1)
points(Fiu$time,u,pch="x")


## simulation of distribution with delayed entry starting at 3
par(mfrow=c(1,1))
F1 <- cbind(c(0,5,8,10),c(0,0.5,0.6,0.9))
F1
plot(F1,ylim=c(0,1),type="l")

n <- 100000
entry <- c(rep(3,10000),runif(n)*7+3)
###entry <- rep(3,n)
u <- runif(n+10000)
###
Fiu <- invsubdist(F1,u,cond=0,entry=entry)
###
# library(prodlim)
# pp <- prodlim(Hist(time,status,entry=entry)~+1,data=Fiu)
# plot(pp,xlim=c(3,10))
###
entry <- 3
###
F1entry <- subdist(F1,entry)[,2]
ptrunc <- 1-F1entry
###
F1entry5 <- F1
F1entry5[,1] <- F1entry5[,1]-entry
F1entry5[,2] <- (F1entry5[,2]-F1entry)/ptrunc
pos <- F1entry5[,1]>=0
F1entry5 <- rbind(c(0,0),F1entry5[pos,])
#
# lines(entry+F1entry5[,1],1-F1entry5[,2],col=2)

##############################################################
## Simulations of two cumulative incidence functions with truncation 
##############################################################
par(mfrow=c(1,1))
F1 <- cbind(c(0,5,8,10),c(0,0.5,0.6,0.9)*0.3)
F2 <- cbind(c(0,5,8,10),c(0,0.5,0.6,0.9)*0.5)
plot(F1,ylim=c(0,1),type="l")
lines(F2,col=2)
entry1 <- 3
###
F1entry <- subdist(F1,entry1)[,2]
F2entry <- subdist(F2,entry1)[,2]
ptrunc <- 1-F1entry-F2entry
###
F1e <- F1
F1e[,1] <- F1e[,1]-entry1
F1e[,2] <- (F1e[,2]-F1entry)/ptrunc
pos <- F1e[,1]>=0
F1e <- rbind(c(0,0),F1e[pos,])
F2e <- F2
F2e[,1] <- F2e[,1]-entry1
F2e[,2] <- (F2e[,2]-F2entry)/ptrunc
pos <- F2e[,1]>=0
F2e <- rbind(c(0,0),F2e[pos,])
#
# truncated identifiable version 
lines(entry1+F1e[,1],F1e[,2],col=1)
lines(entry1+F2e[,1],F2e[,2],col=2)

n <- 10000
entry <- c(rep(entry1,10000),runif(n)*(10-entry1)+entry1)
u <- runif(n+10000)
###
F1entry <- subdist(F1,entry)[,2]
F2entry <- subdist(F2,entry)[,2]
ptrunc <- 1-( F1entry+F2entry)
Fiu1 <- invsubdist(F1,u,cond=1,entry=entry,ptrunc=ptrunc)
Fiu2 <- invsubdist(F1,u,cond=1,entry=entry,ptrunc=ptrunc)
###
ptot <- (tail(F1[,2],1)+tail(F2[,2],1)-F1entry-F2entry)/(ptrunc)
rt <- rbinom(n+10000,1,ptot)
p1 <- ((tail(F1[,2],1)-F1entry)/ptrunc)
p2 <- ((tail(F2[,2],1)-F2entry)/ptrunc)
rb <- rbinom(n+10000,1,p1/ptot)
cause=ifelse(rb==1,1,2)
time=ifelse(cause==1,Fiu1$time,Fiu2$time)
cause <- rt*cause
time[cause==0] <- 10
### simulated data, now checking that things are working 
# pp <- prodlim(Hist(time,cause,entry=entry)~+1)
# plot(pp,xlim=c(entry1,10),cause=1)
# plot(pp,xlim=c(entry1,10),cause=2,add=TRUE)
###
# lines(entry1+F1e[,1],F1e[,2],col=2)
# lines(entry1+F2e[,1],F2e[,2],col=2)

}
\author{
Thomas Scheike
}
\keyword{survival}
