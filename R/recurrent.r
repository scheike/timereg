#' Estimates marginal mean of recurrent events 
#' 
#' Fitting two aalen models for death and recurent events these are
#' combined to prducte the estimator 
#' \deqn{ \int_0^t  S(u) dR(u) } the mean number of recurrent events, here
#' \deqn{ S(u) }  is the probability of survival, and 
#' \deqn{ dR(u) }  is the probability of an event among survivors. 
#' 
#' IID versions used for Ghosh & Lin (2000) variance. See also mets package for 
#' quick version of this for large data mets:::recurrent.marginal, these two 
#' version should give the same when there are no ties. 
#' 
#' @param recurrent aalen model for recurrent events
#' @param death     aalen model for recurrent events
#' @author Thomas Scheike
#' @references 
#' Ghosh and Lin (2002) Nonparametric Analysis of Recurrent events and death, 
#'                      Biometrics, 554--562.
#' @keywords survival
#' @examples
#' \donttest{
#' ### get some data using mets simulaitons 
#' library(mets)
#' data(base1cumhaz)
#' data(base4cumhaz)
#' data(drcumhaz)
#' dr <- drcumhaz
#' base1 <- base1cumhaz
#' base4 <- base4cumhaz
#' rr <- simRecurrent(100,base1,death.cumhaz=dr)
#' rr$x <- rnorm(nrow(rr)) 
#' rr$strata <- floor((rr$id-0.01)/50)
#' drename(rr) <- start+stop~entry+time
#' 
#' ar <- aalen(Surv(start,stop,status)~+1+cluster(id),data=rr,resample.iid=1
#'                                                      ,max.clust=NULL)
#' ad <- aalen(Surv(start,stop,death)~+1+cluster(id),data=rr,resample.iid=1,
#'                                                      ,max.clust=NULL)
#' mm <- recurrent.marginal.mean(ar,ad)
#' with(mm,plot(times,mu,type="s"))
#' with(mm,lines(times,mu+1.96*se.mu,type="s",lty=2))
#' with(mm,lines(times,mu-1.96*se.mu,type="s",lty=2))
#' }
##' @export
recurrent.marginal.mean <- function(recurrent,death) 
{# {{{
  axr <- recurrent
  adr <- death
  St <- exp(-adr$cum[,2])
  ### construct iid at same time points , combined jump-points 
  timesr <- axr$cum[,1]
  timesd <- adr$cum[,1]
  times <- c(timesr[-1],timesd[-1])
  or    <- order(times)
  times <- times[or]
  keepr <- order(or)[1:length(timesr[-1])]

  rid <- sindex.prodlim(timesd,times,strict=FALSE)
  rir <- sindex.prodlim(timesr,times,strict=FALSE)
  Stt <- St[rid]
  ###
  ariid <- axr$cum[rir,2]
  mu <- cumsum(Stt*diff(c(0,ariid)))
  nc <- length(axr$B.iid)
  muiid <-  matrix(0,length(times),nc)
###  muiid1 <- matrix(0,length(times),nc)
###  muiid2 <- matrix(0,length(times),nc)
###  muiid3 <- matrix(0,length(times),nc)
  for (i in  1:nc) {
          mriid <- axr$B.iid[[i]]
	  mdiid <- adr$B.iid[[i]]
	  mriid <- mriid[rir] 
	  mdiid <- mdiid[rid]
	  dmridd <- diff(c(0,mriid))
	  dmdidd <- diff(c(0,mdiid))
	  muiid[,i]  <- cumsum(Stt*dmridd)-mu*cumsum(dmdidd)+cumsum(mu*dmdidd)
###	  muiid1[,i] <- cumsum(Stt*dmridd)
###	  muiid2[,i] <- -mu*cumsum(dmdidd)
###	  muiid3[,i] <- cumsum(mu*dmdidd)
  }
  var1 <- apply(muiid^2,1,sum)
###  varl1 <- apply(muiid1^2,1,sum)
###  varl2 <- apply(muiid2^2,1,sum)
###  varl3 <- apply(muiid3^2,1,sum)
######
###  cov12 <- apply(muiid1*muiid2,1,sum)
###  cov13 <- apply(muiid1*muiid3,1,sum)
###  cov23 <- apply(muiid2*muiid3,1,sum)
  out=list(times=times[keepr],mu=mu[keepr],var.mu=var1[keepr],se.mu=var1[keepr]^.5,
	   St=St,Stt=Stt[keepr])
###	   covs=cbind(cov12,cov13,cov23)[keepr,],
###	   vari=cbind(varl1,varl2,varl3)[keepr,])
}# }}}

#' Estimates marginal mean of recurrent events  based on two cox models 
#' 
#' Fitting two Cox models for death and recurent events these are
#' combined to prducte the estimator 
#' \deqn{ \int_0^t  S(u|x=0) dR(u|x=0) }{} the mean number of recurrent events, here
#' \deqn{ S(u|x=0) }{}  is the probability of survival, and 
#' \deqn{ dR(u|x=0) }{}  is the probability of an event among survivors. 
#' For now the estimator is based on the two-baselines so \deqn{x=0}{}, but covariates
#' can be rescaled to look at different x's and extensions possible. 
#' 
#' IID versions along the lines of Ghosh & Lin (2000) variance. See also mets package for 
#' quick version of this for large data. 
#' IID versions used for Ghosh & Lin (2000) variance. See also mets package for 
#' quick version of this for large data mets:::recurrent.marginal, these two 
#' version should give the same when there are now ties. 
#' 
#' @param recurrent aalen model for recurrent events
#' @param death     cox.aalen (cox) model for death events
#' @author Thomas Scheike
#' @references 
#' Ghosh and Lin (2002) Nonparametric Analysis of Recurrent events and death, 
#'                      Biometrics, 554--562.
#' @keywords survival
#' @examples
#' \donttest{
#' ### do not test because iid slow  and uses data from mets
#' library(mets)
#' data(base1cumhaz)
#' data(base4cumhaz)
#' data(drcumhaz)
#' dr <- drcumhaz
#' base1 <- base1cumhaz
#' base4 <- base4cumhaz
#' rr <- simRecurrent(100,base1,death.cumhaz=dr)
#' rr$x <- rnorm(nrow(rr)) 
#' rr$strata <- floor((rr$id-0.01)/50)
#' drename(rr) <- start+stop~entry+time
#'
#' ar <- cox.aalen(Surv(start,stop,status)~+1+prop(x)+cluster(id),data=rr,
#'                    resample.iid=1,,max.clust=NULL,max.timepoint.sim=NULL)
#' ad <- cox.aalen(Surv(start,stop,death)~+1+prop(x)+cluster(id),data=rr,
#'                    resample.iid=1,,max.clust=NULL,max.timepoint.sim=NULL)
#' mm <- recurrent.marginal.coxmean(ar,ad)
#' with(mm,plot(times,mu,type="s"))
#' with(mm,lines(times,mu+1.96*se.mu,type="s",lty=2))
#' with(mm,lines(times,mu-1.96*se.mu,type="s",lty=2))
#' }
##' @export
recurrent.marginal.coxmean <- function(recurrent,death) 
{# {{{
out <- recurrent.marginal.mean(recurrent,death)
return(out)
}# }}}

