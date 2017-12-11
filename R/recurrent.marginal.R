##' Fast recurrent marginal mean when death is possible
##'
##' Fast Marginal means of recurrent events. Using the Lin and Ghosh (2000) 
##' standard errors.  
##' Fitting two models for death and recurent events these are
##' combined to prducte the estimator 
##' \deqn{ \int_0^t  S(u|x=0) dR(u|x=0) } the mean number of recurrent events, here
##' \deqn{ S(u|x=0) }  is the probability of survival for the baseline group, and 
##' \deqn{ dR(u|x=0) }  is the hazard rate of an event among survivors for the baseline. 
##' Here \deqn{ S(u|x=0) }  is estimated by \deqn{ exp(-\Lambda_d(u|x=0) }  with 
##'  \deqn{\Lambda_d(u|x=0) } being the cumulative baseline for death.
##' 
##' Assumes no ties so break ties before running with for example strata, use tie.breaker function. 
##' 
##' @param recurrent phreg object with recurrent events
##' @param death     phreg object with deaths
##' @param fixbeta   to force the estimation of standard errors to think of regression coefficients as known/fixed
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases recurrent.marginal 
##' @examples
##' 
##' data(simrecurrent)
##' ##  to fit non-parametric models with just a baseline 
##' xr <- phreg(Surv(start,stop,status)~cluster(id),data=simd)
##' dr <- phreg(Surv(start,stop,death)~cluster(id),data=simd)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE)
##' ### robust standard errors 
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' basehazplot.phreg(rxr,se=TRUE,robust=TRUE,add=TRUE,col=4)
##' 
##' ## marginal mean of expected number of recurrent events 
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=2)
##' 
##' ########################################################################
##' ###   with strata     ##################################################
##' ########################################################################
##' xr <- phreg(Surv(start,stop,status)~strata(x.V1)+cluster(id),data=simd)
##' dr <- phreg(Surv(start,stop,death)~strata(x.V1)+cluster(id),data=simd)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE)
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' basehazplot.phreg(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)
##'
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=1:2)
##'
##' ########################################################################
##' ###   cox case        ##################################################
##' ########################################################################
##' xr <- phreg(Surv(start,stop,status)~x.V1+x.V2+cluster(id),data=simd)
##' dr <- phreg(Surv(start,stop,death)~x.V1+x.V2+cluster(id),data=simd)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE)
##' rxr <-   robust.phreg(xr)
##' basehazplot.phreg(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)
##'
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=1:2)
##' 
##' \donttest{
##' ### do not test to avoid dependence on frailtypack for data 
##' library(frailtypack)
##' data(readmission)
##' r <-  readmission
##' ### breaking ties with respect to both event types 
##' r <-  transform(r,jump=(event==1) | (death==1))
##' r <- tie.breaker(r,stop="t.stop",start="t.start",status="jump")
##'
##' ## to fit non-parametric models with just a baseline 
##' xr <- phreg(Surv(t.start,t.stop,event)~cluster(id),data=r)
##' dr <- phreg(Surv(t.start,t.stop,death)~cluster(id),data=r)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE,xlim=c(0,2000))
##' ### robust standard errors 
##' rxr <-   robust.phreg(xr)
##' basehazplot.phreg(rxr,se=TRUE,xlim=c(0,2000),robust=TRUE,add=TRUE,col=4)
##' title(main="recurrent")
##' 
##' ## marginal mean of expected number of recurrent events 
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=2)
##' title(main="marginal mean ")
##' 
##' ########################################################################
##' ###   with strata     ##################################################
##' ########################################################################
##' xr <- phreg(Surv(t.start,t.stop,event)~strata(chemo)+cluster(id),data=r)
##' dr <- phreg(Surv(t.start,t.stop,death)~strata(chemo)+cluster(id),data=r)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE,xlim=c(0,2000),ylim=c(0,3))
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' basehazplot.phreg(rxr,se=TRUE,xlim=c(0,2000),robust=TRUE,add=TRUE,col=1:2)
##' title(main="recurrent")
##'
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=1:2,xlim=c(0,2000),ylim=c(0,3))
##' title(main="marginal mean ")
##'
##' ########################################################################
##' ###   cox case        ##################################################
##' ########################################################################
##' xr <- phreg(Surv(t.start,t.stop,event)~chemo+charlson+cluster(id),data=r)
##' dr <- phreg(Surv(t.start,t.stop,death)~charlson+dukes+cluster(id),data=r)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE,xlim=c(0,2000),ylim=c(0,3))
##' rxr <-   robust.phreg(xr)
##' basehazplot.phreg(rxr,se=TRUE,xlim=c(0,2000),robust=TRUE,add=TRUE,col=1:2)
##' title(main="recurrent")
##'
##' ## for covariates equal to 0
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=1:2,xlim=c(0,2000),ylim=c(0,3))
##' title(main="marginal mean ")
##' 
##' }
##' 
##' ########################################################################
##' ###   CIF             ##################################################
##' ########################################################################
##' ### use of function to compute cumulative incidence (cif) with robust standard errors
##' data(bmt)
##' bmt$id <- 1:nrow(bmt)
##' xr  <- phreg(Surv(time,cause==1)~cluster(id),data=bmt)
##' dr  <- phreg(Surv(time,cause!=0)~cluster(id),data=bmt)
##' 
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="cumulative incidence")
##' 
##' @export
recurrent.marginal <- function(recurrent,death,fixbeta=NULL,...)
{# {{{
  xr <- recurrent
  dr <- death 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if (is.null(xr$opt) | is.null(xr$coef)) fixbeta<- 1 else fixbeta <- 0

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  Gt <- exp(-dr$cumhaz[,2])
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  xx$time[xx$jumps+1]
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
  ###
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazR <-  cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
###  mm <- cbind(cumhazR,cumhazDR,St,xx$time,x$cox.prep$time)
###  head(cbind(cumhazR,cumhazDR,St,xx$time,x$cox.prep$time),500)
  mu <- cumhazDR[,2]
# }}}

  ### robust standard errors 
  ### 1. sum_k ( int_0^t S(s)/S_0^r(s) dM_k.^r(s) )^2
# {{{
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
  U[xx$jumps+1,] <- x$U
  ###    
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumS0i2 <- cumsumstrata(St*S0i2,xx$strata,xx$nstrata)
  if (fixbeta==0) {
  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  Ht <- apply(St*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  HtR <- Ht
  }
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
  else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  ### also weights 
  w <- c(xx$weights)
  xxx <- w*(St*S0i-rr*c(cumS0i2))
  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
  ss <- c(revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2^2)
  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2)
  varA1 <- c(ssf+ss+2*covv)
  cumS0i2R <- cumS0i2 
  xxxR <- xxx
  rrR <- rr

  if (fixbeta==0) {# {{{
      invhess <- -solve(x$hessian)
      MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-EdLam0)*rr*c(xx$weights)
      UU <- apply(MGt,2,sumstrata,id,max(id)+1)
      betaiidR <- UU %*% invhess
###     betaiidR <- iid(x)
     vbeta <- crossprod(betaiidR)
     varbetat <-   rowSums((Ht %*% vbeta)*Ht)
     ### writing each beta for all individuals 
     betakt <- betaiidR[id+1,,drop=FALSE]
     ###
     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- c(covk2)*c(cumS0i2)
     ###
     varA1 <- varA1+varbetat-2*apply((covk1-covk2)*Ht,1,sum)
  }# }}}
# }}}
  ### 2. mu(t)^2 * sum_k ( int_0^t 1/S_0^d(s) dM_k.^d(s) )^2
# {{{
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  if (fixbeta==0) {
	  Z <- xx$X
	  U <- E <- matrix(0,nrow(xx$X),x$p)
	  E[xx$jumps+1,] <- x$E
	  U[xx$jumps+1,] <- x$U
	  Ht <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  }
  ###    
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumS0i2 <- cumsumstrata(S0i2,xx$strata,xx$nstrata)
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
  else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  ### also weights 
  w <- c(xx$weights)
  xxx1 <- w*(S0i-rr*c(cumS0i2))
  ssf <- cumsumidstratasum(xxx1,id,mid,xx$strata,xx$nstrata)$sumsquare
  ss <- c(revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2^2)
  covv <- covfridstrata(xxx1,w*rr,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2)
  varA2 <- mu^2*c(ssf+ss+2*covv)
  cumS0i2D1 <- cumS0i2 
  xxxD1 <- xxx1
  rrD1 <- rr
  if (fixbeta==0) {# {{{
      HtD1 <- mu*Ht
      invhess <- -solve(x$hessian)
      MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-Ht)*rr*c(xx$weights)
      UU <- apply(MGt,2,sumstrata,id,max(id)+1)
      betaiidD <- UU %*% invhess
###     betaiidD1 <- iid(x)
     vbetaD <- crossprod(betaiidD)
     varbetat <-   rowSums((HtD1 %*% vbetaD)*HtD1)
     ### writing each beta for all individuals 
     betakt <- betaiidD[id+1,,drop=FALSE]
     ###
     covk1 <- apply(xxx1*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- c(covk2)*c(cumS0i2D1)
     ###
     varA2 <- varA2+varbetat-2*apply((covk1-covk2)*HtD1*mu,1,sum)
  }# }}}

# }}}
  ### 3. sum_k( int_0^t mu(s) /S_0^d(s) dM_k.^d(s))^2
# {{{
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  if (fixbeta==0) {
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
  }
###  U[xx$jumps+1,] <- x$U
  xr$cox.prep$time  - xx$time
  xr$cox.prep$time[xr$cox.prep$jumps+1]
  xr$cox.prep$time[dr$cox.prep$jumps+1]
  ###    
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumS0i2 <- cumsumstrata(mu*S0i2,xx$strata,xx$nstrata)

  if (fixbeta==0) {
  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  Ht <- apply(mu*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  HtD <- Ht
  }
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
  else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  ### also weights 
  w <- c(xx$weights)
  xxx <- w*(mu*S0i-rr*c(cumS0i2))
  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
  ss <- c(revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2^2)
  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2)
  varA3 <- c(ssf+ss+2*covv)
  cumS0i2D <- cumS0i2
  rrD <- rr
  xxxD <- xxx
 
  if (fixbeta==0) {# {{{
     invhess <- -solve(x$hessian)
     varbetat <-   rowSums((Ht %*% vbetaD)*Ht)
     ### writing each beta for all individuals 
     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- c(covk2)*c(cumS0i2)
     ###
     varA3 <- varA3+varbetat-2*apply((covk1-covk2)*Ht,1,sum)
  }# }}}
  varA <-  varA1+varA2+varA3 
# }}}

 ### covariances between different terms  13 23  12 12
# {{{
 cov23 <- c(cumsumidstratasumCov(xxxD,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)*mu
 cov232 <- c(revcumsumidstratasum(w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2D*mu*cumS0i2D1)
### cov233 <- covfridstrataCov(xxxD,w*rrD,xxxD1,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1*mu)
 cov233 <- covfridstrata(xxxD,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1*mu)
 cov234 <- covfridstrata(xxxD1,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D*mu)
 cov23A <- -c(cov23+(cov232+cov233+cov234))
###
 cov12 <-  c(cumsumidstratasumCov(xxxR,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)*mu
 cov122 <- c(revcumsumidstratasumCov(w*rrR,w*rrD1,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2R*mu*cumS0i2D1)
 cov123 <- covfridstrataCov(xxxR,w*rrR,xxxD,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1*mu)
 cov124 <- mu*covfridstrataCov(xxxD1,w*rrD1,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2R)
 cov12A <- -c(cov12+(cov122+cov123+cov124))
###
 cov13 <-  c(cumsumidstratasumCov(xxxR,xxxD,id,mid,xx$strata,xx$nstrata)$sumsquare)
 cov132 <- c(revcumsumidstratasumCov(w*rrR,w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2R*cumS0i2D)
 cov133 <- covfridstrataCov(xxxR,w*rrR,xxxD,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D)
 cov134 <- covfridstrataCov(xxxD,w*rrD,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2R)
###
 cov13A <- c(cov13+(cov132+cov133+cov134))
 varA <- varA+2*cov12A+2*cov23A+2*cov13A
# }}}
 

 cov12aa <- cov13aa <- cov23aa <- 0
 if (fixbeta==0) {
### covariances between different terms and  beta's 
# {{{
 covbetaRD <- t(betaiidR) %*% betaiidD
 DHt <- HtD1-HtD
### covbeta1.23 <-   -2*rowSums((HtR %*% covbetaRD)*DHt)
 covbeta1.12 <-   -2*rowSums((HtR %*% covbetaRD)*HtD1)
 covbeta1.13 <-   2*rowSums((HtR %*% covbetaRD)*HtD)
 covbetaD.23 <-   -2*rowSums((HtD %*% vbetaD)*(HtD1))

 ### D versus betaD from two terms  cov23 wrt beta 
 betakt <- betaiidD[id+1,,drop=FALSE]
 covk1 <-apply(xxxD1*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2 <-apply(w*rrD1*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2D1)
 covD1.D <- 2*apply((covk1-covk2)*mu*HtD,1,sum)
 ###
 covk1 <-apply(xxxD*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2 <-apply(w*rrD*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2D)
 covD.D1 <- 2*apply((covk1-covk2)*HtD1,1,sum)
 cov23aa <- covbetaD.23  + covD1.D + covD.D1

 ### cov12 wrt betaD and betaR
 betakt <- betaiidD[id+1,,drop=FALSE]
 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2R)
 covRD12 <- 2*apply((covk1-covk2)*HtD1,1,sum)
###
 betakt <- betaiidR[id+1,,drop=FALSE]
 covk1 <-apply(xxxD1*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2 <-apply(w*rrD1*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2D1)
 covRD21 <- 2*apply((covk1-covk2)*mu*HtR,1,sum)
 cov12aa <- covbeta1.12 + covRD12+covRD21


 ### cov13 wrt betaD and betaR
 betakt <- betaiidD[id+1,,drop=FALSE]
 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2R)
 covRD13 <- -2*apply((covk1-covk2)*HtD,1,sum)
###
 betakt <- betaiidR[id+1,,drop=FALSE]
 covk1 <-apply(xxxD*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2 <-apply(w*rrD*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2D)
 covRD31 <- -2*apply((covk1-covk2)*HtR,1,sum)
 cov13aa <- covbeta1.13 + covRD13+covRD31

 varA <-varA+ cov23aa+cov13aa+cov12aa
# }}}
 }

###
 varrs <- data.frame(mu=mu,cumhaz=mu,se.mu=varA^.5,time=xr$time,se.cumhaz=varA^.5,strata=xr$strata,St=St)
###     varA1=varA1,varA2=varA2,varA3=varA3,cov12=2*cov12A+cov12aa,cov13=2*cov13A+cov13aa, 
###     cov23=2*cov23A+cov23aa); 
 varrs <- varrs[c(xr$cox.prep$jumps)+1,]

 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,se.mu=varrs$se.mu,times=varrs$time,
     St=varrs$St,
     cumhaz=cbind(varrs$time,varrs$mu),se.cumhaz=cbind(varrs$time,varrs$se.mu),
     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),strata.name=xr$strata.name)
###  vari=varrs[,c("varA1","varA2","varA3")],covs=varrs[,c("cov12","cov13","cov23")])
 return(out)
}# }}}

##' @export
tie.breaker <- function(data,stop="time",start="entry",status="status",id=NULL,ddt=NULL)
 {# {{{
   if (!is.null(id)) id <- data[,id]
   ord <- 1:nrow(data)
   stat <- data[,status]
   time <- data[,stop]
   time1 <- data[stat==1,stop]
   ties <- duplicated(c(time1))
   nties <- sum(ties)
   ordties <- ord[stat==1][ties]
   if (is.null(ddt)) ddt <- min(abs(diff(data[,stop])))*0.5
   time[ordties] <- time[ordties]+runif(nties)*ddt

   data[ordties,stop] <- time[ordties]
   ties <- (ord %in% ordties)
   if (!is.null(id)) {
   lagties <- dlag(ties)
   ### also move next start time if id the same 
   change.start <- lagties==TRUE & id==dlag(id)
   change.start[is.na(change.start)] <- FALSE
   ocs <- ord[change.start]
   data[ocs,start] <- data[ocs-1,stop]
   data[,"tiebreaker"] <- FALSE
   data[ocs,"tiebreaker"] <- TRUE
   }
   
   return(data)
 } # }}}

##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Must give hazard of death and recurrent events.  Possible with two
##' event types and their dependence can be specified but the two recurrent events need
##' to share random effect.
##' combined to prducte the estimator 
##'
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param cumhaz.death cumulative hazard of death 
##' @param cumhaz2  cumulative hazard of recurrent events  of type 2
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param dhaz rate for death hazard if it is extended to time-range of first event 
##' @param haz2 rate of second cause  if it is extended to time-range of first event 
##' @param dependence  =0 independence, =1 all share same random effect with variance var.z
##'                    =2 random effect exp(normal) with correlation structure from cor.mat,
##'                    first random effect is z1 and shared for a possible second cause,
##'                    second random effect is for death 
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##'
##' ## getting some rates to mimick 
##' ########################################
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##'
##'  cor.mat <- corM <- rbind(c(1.0, 0.6, 0.9), c(0.6, 1.0, 0.5), c(0.9, 0.5, 1.0))
##'
##'  ######################################################################
##'  ### simulating simple model that mimicks data 
##'  ######################################################################
##'  rr <- sim.recurrent(10,base1,death.cumhaz=dr)
##'  dlist(rr,.~id,n=0)
##'
##'  rr <- sim.recurrent(1000,base1,death.cumhaz=dr)
##'  par(mfrow=c(1,3))
##'  showfitsim(causes=1)
##'
##' ######################################################################
##' ### simulating simple model that mimicks data 
##' ### now with two event types and second type has same rate as death rate
##' ######################################################################
##'
##'  rr <- sim.recurrent(1000,base1,cumhaz2=base4,death.cumhaz=dr)
##'  dtable(rr,~death+status)
##'  par(mfrow=c(2,2))
##'  showfitsim(causes=2)
##'
##' ######################################################################
##' ### simulating simple model 
##' ### random effect for all causes (Z shared for death and recurrent) 
##' ######################################################################
##'
##'  rr <- sim.recurrent(1000,base1,
##'         death.cumhaz=dr,dependence=1,var.gamma=0.4)
##'  ### marginals do fit after input after integrating out
##'  par(mfrow=c(2,2))
##'  showfitsim(causes=1)
##'
##' @export
sim.recurrent <- function(n,cumhaz,death.cumhaz=NULL,cumhaz2=NULL,
			    gap.time=FALSE,
			    max.recurrent=100,dhaz=NULL,haz2=NULL,
			    dependence=0,var.z=2,cor.mat=NULL,...) 
  {# {{{

  if (dependence==0) { z1 <- z2 <- zd <- rep(1,n) } else if (dependence==1) {
###	      zz <- rgamma(n,1/var.gamma[1])*var.gamma[1]
	      zz <- exp(rnorm(n,1)*var.z[1]^.5)
	      var(zz)
	      z1 <- zz; z2 <- zz; zd <- zz
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(2*n),n,2)
	      z <- (z%*% chol(covv))
	      z1 <- exp(z[,1]); zd <- exp(z[,3])
	      apply(exp(z),2,mean); cov(exp(z))
      }

  cumhaz <- rbind(c(0,0),cumhaz)
  ll <- nrow(cumhaz)
  max.time <- tail(cumhaz[,1],1)

  ## extend cumulative for cause 2 to full range of cause 1
  if (!is.null(cumhaz2)) {# {{{
      cumhaz2 <- rbind(c(0,0),cumhaz2)
	      if (is.null(haz2)) 
		      haz2 <- tail(cumhaz2[,2],1)/tail(cumhaz2[,1],1)
	      ### linear extrapolation of mortality using given dhaz or 
      if (tail(cumhaz2[,1],1)<max.time) {
	      cumhaz2 <- rbind(cumhaz2,c(max.time,haz2*max.time)) 
       }
  }# }}}

  ## extend cumulative for death to full range  of cause 1
  if (!is.null(death.cumhaz)) {# {{{
      cumhazd <- rbind(c(0,0),death.cumhaz) 
     if (tail(death.cumhaz[,1],1)<max.time) {
	      ### linear extrapolation of mortality using given dhaz or 
	      if (is.null(dhaz)) 
	        dhaz <- tail(death.cumhaz[,2],1)/tail(death.cumhaz[,1],1)
              cumhazd <- rbind(cumhazd,c(max.time,dhaz*max.time)) 
       } 
  }# }}}

  ## sum two cumulatives to get combined events 
  if (!is.null(cumhaz2)) {# {{{
	  times <- sort(unique(c(cumhaz[,1],cumhaz2[,1])))
          cumhaz1t <- approx(cumhaz[,1],cumhaz[,2],times,rule=2)$y
          cumhaz2t <- approx(cumhaz2[,1],cumhaz2[,2],times,rule=2)$y
	  cumhaz1 <- cumhaz
	  cumhaz <- cbind(times,cumhaz1t+cumhaz2t)
  }# }}}

### recurrent first time
  tall <- pc.hazard(cumhaz,z1)
  tall$id <- 1:n
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- pc.hazard(cumhazd,zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
  } else { tall$dtime <- max.time; tall$fdeath <- 0; cumhazd <- NULL }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=1,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
###  gemsim <- as.data.frame(matrix(0,max.recurrent*n,ncol(tall)))
###  names(gemsim) <- names(tall)
###  gemsim[1:n,] <- tall; nrr <- n
  i <- 1; 
  while (any(tt$time<tt$dtime) & i < max.recurrent) {
	  i <- i+1
	  still <- subset(tt,time<dtime)
          tt <- pc.hazard(cumhaz,rr=z1[still$id],entry=still$time)
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath))
	  tt <- dtransform(tt,death=1,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
###	  gemsim[(nrr+1):(nrr+nt),] <- tt
	  tall <- rbind(tall,tt)
###	  nrr <- nrr+nt
  }
###  tall <- gemsim[1:nrr,]
  dsort(tall) <- ~id+entry+time

  ### cause 2 is there then decide if jump is 1 or 2
  if (!is.null(cumhaz2)) {# {{{
      haz1 <- apply(cumhaz1,2,diff)
      haz1 <- haz1[,2]/haz1[,1]
      ## hazard2 at times 
      haz2 <- apply(cumhaz2,2,diff)
      haz2 <- haz2[,2]/haz2[,1]
      ll1 <- nrow(cumhaz1)
      ll2 <- nrow(cumhaz2)
      haz1t <- Cpred(cbind(cumhaz1[-ll1,1],haz1),tall$time)[,2]
      haz2t <- Cpred(cbind(cumhaz2[-ll2,1],haz2),tall$time)[,2]
      p2t <- haz2t/(haz1t+haz2t)
      tall$p2t <- p2t
      tall$status <- (1+rbinom(nrow(tall),1,p2t))*(tall$status>=1)
  }# }}}

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2

  return(tall)
  }# }}}


##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Must give hazard of death and two recurrent events.  Possible with two
##' event types and their dependence can be specified but the two recurrent events need
##' to share random effect. Based on drawing the from cumhaz and cumhaz2 and taking the first event rather
##' the cumulative and then distributing it out. Key advantage of this is that there is 
##' more flexibility wrt random effects 
##'
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param cumhaz2  cumulative hazard of recurrent events  of type 2
##' @param cumhaz.death cumulative hazard of death 
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param dhaz rate for death hazard if it is extended to time-range of first event 
##' @param haz2 rate of second cause  if it is extended to time-range of first event 
##' @param dependence  =0 independence, =1 all share same random effect with variance var.z
##'                    =2 random effect exp(normal) with correlation structure from cor.mat,
##'                    first random effect is z1 and for N1
##'                    second random effect is z2 and for N2
##'                    third random effect is for death 
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##'
##' ########################################
##' ## getting some rates to mimick 
##' ########################################
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##'
##'  cor.mat <- corM <- rbind(c(1.0, 0.6, 0.9), c(0.6, 1.0, 0.5), c(0.9, 0.5, 1.0))
##' 
##' ######################################################################
##' ### simulating simple model that mimicks data 
##' ### now with two event types and second type has same rate as death rate
##' ######################################################################
##'
##' rr <- sim.recurrentII(1000,base1,dr,death.cumhaz=base4)
##' dtable(rr,~death+status)
##' par(mfrow=c(2,2))
##' showfitsim(causes=2)
##'
##' @export
sim.recurrentII <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,
		    gap.time=FALSE,max.recurrent=100,dhaz=NULL,haz2=NULL,
		    dependence=0,var.z=0.22,cor.mat=NULL,...) 
  {# {{{

  if (dependence==0) { z <- z1 <- z2 <- zd <- rep(1,n) 
     } else if (dependence==1) {
###	      zz <- rgamma(n,1/var.gamma[1])*var.gamma[1]
	      z <- exp(rnorm(n,1)*var.z[1]^.5)
	      z1 <- z; z2 <- z; zd <- z
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      }

  cumhaz <- rbind(c(0,0),cumhaz)
  ll <- nrow(cumhaz)
  max.time <- tail(cumhaz[,1],1)

  ## extend cumulative for cause 2 to full range of cause 1
  if (!is.null(cumhaz2)) {# {{{
      cumhaz2 <- rbind(c(0,0),cumhaz2)
	      if (is.null(haz2)) 
		      haz2 <- tail(cumhaz2[,2],1)/tail(cumhaz2[,1],1)
	      ### linear extrapolation of mortality using given dhaz or 
      if (tail(cumhaz2[,1],1)<max.time) {
	      cumhaz2 <- rbind(cumhaz2,c(max.time,haz2*max.time)) 
       }
  }# }}}


  ## extend cumulative for death to full range  of cause 1
  if (!is.null(death.cumhaz)) {# {{{
     cumhazd <- rbind(c(0,0),death.cumhaz)
     if (tail(death.cumhaz[,1],1)<max.time) {
	      ### linear extrapolation of mortality using given dhaz or 
	      if (is.null(dhaz)) 
                 dhaz <- tail(death.cumhaz[,2],1)/tail(death.cumhaz[,1],1)
              cumhazd <- rbind(cumhazd,c(max.time,dhaz*max.time)) 
       }
  }# }}}


### recurrent first time
  tall1 <- pc.hazard(cumhaz,z1)
  tall2 <- pc.hazard(cumhaz2,z2)
  tall <- tall1 
  tall$status <- ifelse(tall1$time<tall2$time,tall1$status,2*tall2$status)
  tall$time <- ifelse(tall1$time<tall2$time,tall1$time,tall2$time)
  tall$id <- 1:n
  tall$rr2 <- tall2$rr
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- pc.hazard(cumhazd,zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
  } else { tall$dtime <- max.time; tall$fdeath <- 0; cumhazd <- NULL }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=1,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  ### setting aside memory 
  tt1 <- tt2 <- tt
###  gemsim <- as.data.frame(matrix(0,max.recurrent*n,ncol(tall)))
###  names(gemsim) <- names(tall)
###  gemsim[1:n,] <- tall; nrr <- n
  i <- 1; 
  while (any(tt$time<tt$dtime) & i < max.recurrent) {
	  i <- i+1
	  still <- subset(tt,time<dtime)
	  nn <- nrow(still)
          tt1 <- pc.hazard(cumhaz,rr=z1[still$id],entry=still$time)
          tt2 <- pc.hazard(cumhaz2,rr=z2[still$id],entry=still$time)
	  tt <- tt1
###          drename(tt1,paste(names(tt1),"1",sep="")) <- ~.
###          drename(tt2,paste(names(tt2),"2",sep="")) <- ~.
          tt$status <- ifelse(tt1$time<=tt2$time,tt1$status,2*tt2$status)
          tt$time <-   ifelse(tt1$time<=tt2$time,tt1$time,tt2$time)
	  tt$rr2 <- tt2$rr
###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath))
	  tt <- dtransform(tt,death=1,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
###	  gemsim[(nrr+1):(nrr+nt),] <- tt
	  tall <- rbind(tall,tt[1:nn,])
###	  nrr <- nrr+nt
  }
###  tall <- gemsim[1:nrr,]
  dsort(tall) <- ~id+entry+time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"z") <- z

  return(tall)
  }# }}}

##' @export
count.history <- function(data,status="status",id="id",types=1:2,names.count="Count")
{# {{{
stat <- data[,status]

clusters <- data[,id]
if (is.numeric(clusters)) {
      clusters <- fast.approx(unique(clusters), clusters) - 1
      max.clust <- max(clusters)
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

data[,"lbnr__id"] <- cumsumstrata(rep(1,nrow(data)),clusters,max.clust+1) 
for (i in types)  {
data[,paste(names.count,i,sep="")] <- 
   cumsumidstratasum((stat==i),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
}

return(data)
}# }}}

##' @export
showfitsim <- function(causes=2) 
{# {{{
  drr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)
  basehazplot.phreg(drr,ylim=c(0,8))
  lines(dr,col=2)
###
  xrr <- phreg(Surv(entry,time,status==1)~cluster(id),data=rr)
  basehazplot.phreg(xrr,add=TRUE)
###  basehazplot.phreg(xrr)
  lines(base1,col=2)
  if (causes>=2) {
	  xrr2 <- phreg(Surv(entry,time,status==2)~cluster(id),data=rr)
	  basehazplot.phreg(xrr2,add=TRUE)
	  lines(base4,col=2)
  }
  meanr1 <-   recurrent.marginal(xrr,drr)
  basehazplot.phreg(meanr1,se=TRUE)
  if (causes>=2) {
	  meanr2 <-   recurrent.marginal(xrr2,drr)
	  basehazplot.phreg(meanr2,se=TRUE,add=TRUE,col=2)
  }
}# }}}


