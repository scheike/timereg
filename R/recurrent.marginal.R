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
##' ## no.opt to fit non-parametric models with just a baseline 
##' ## covariate just given to make cox call  possible 
##' xr <- phreg(Surv(start,stop,status)~x.V1+x.V2+cluster(id),data=simd,no.opt=TRUE)
##' dr <- phreg(Surv(start,stop,death)~x.V1+x.V2+cluster(id),data=simd,no.opt=TRUE)
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
##' xr <- phreg(Surv(start,stop,status)~strata(x.V1)+x.V2+cluster(id),data=simd,no.opt=TRUE)
##' dr <- phreg(Surv(start,stop,death)~strata(x.V1)+x.V2+cluster(id),data=simd,no.opt=TRUE)
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
##'
##' ## no.opt to fit non-parametric models with just a baseline 
##' ## covariate just given to make cox call  possible 
##' xr <- phreg(Surv(t.start,t.stop,event)~charlson+cluster(id),data=r,no.opt=TRUE)
##' dr <- phreg(Surv(t.start,t.stop,death)~charlson+cluster(id),data=r,no.opt=TRUE)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE,xlim=c(0,2000))
##' ### robust standard errors 
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' basehazplot.phreg(rxr,se=TRUE,xlim=c(0,2000),robust=TRUE,add=TRUE,col=4)
##' 
##' ## marginal mean of expected number of recurrent events 
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=2)
##' 
##' ########################################################################
##' ###   with strata     ##################################################
##' ########################################################################
##' xr <- phreg(Surv(t.start,t.stop,event)~strata(chemo)+charlson+cluster(id),data=r,no.opt=TRUE)
##' dr <- phreg(Surv(t.start,t.stop,death)~strata(chemo)+charlson+cluster(id),data=r,no.opt=TRUE)
##' par(mfrow=c(1,3))
##' basehazplot.phreg(dr,se=TRUE)
##' title(main="death")
##' basehazplot.phreg(xr,se=TRUE,xlim=c(0,2000),ylim=c(0,3))
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' basehazplot.phreg(rxr,se=TRUE,xlim=c(0,2000),robust=TRUE,add=TRUE,col=1:2)
##'
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=1:2,xlim=c(0,2000),ylim=c(0,3))
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
##'
##' out <- recurrent.marginal(xr,dr)
##' basehazplot.phreg(out,se=TRUE,ylab="marginal mean",col=1:2,xlim=c(0,2000),ylim=c(0,3))
##' 
##' }
##' 
##' ########################################################################
##' ###   CIF             ##################################################
##' ########################################################################
##' ### use of function to compute cumulative incidence (cif) with robust standard errors
##' bmt$id <- 1:nrow(bmt)
##' xr  <- phreg(Surv(time,cause==1)~age+cluster(id),data=bmt,no.opt=TRUE)
##' dr  <- phreg(Surv(time,cause!=0)~age+cluster(id),data=bmt,no.opt=TRUE)
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
  if (is.null(xr$opt)) fixbeta<- 1 else fixbeta <- 0

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
  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  Ht <- apply(St*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
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
  HtR <- Ht

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
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
  U[xx$jumps+1,] <- x$U
  ###    
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumS0i2 <- cumsumstrata(S0i2,xx$strata,xx$nstrata)
  Ht <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
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
  HtD1 <- mu*Ht
  if (fixbeta==0) {# {{{
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
  Z <- xx$X
  U <- E <- matrix(0,nrow(xx$X),x$p)
  E[xx$jumps+1,] <- x$E
###  U[xx$jumps+1,] <- x$U
  xr$cox.prep$time  - xx$time
  xr$cox.prep$time[xr$cox.prep$jumps+1]
  xr$cox.prep$time[dr$cox.prep$jumps+1]
  ###    
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumS0i2 <- cumsumstrata(mu*S0i2,xx$strata,xx$nstrata)

  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  Ht <- apply(mu*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
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
  HtD <- Ht
 
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


