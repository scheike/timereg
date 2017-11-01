##' Fast recurrent marginal mean when death is possible
##'
##' Fast Marginal means of recurrent events. Using the Lin and Ghosh (2000) 
##' standard errors.  
##' @param recurrent phreg object with recurrent events
##' @param death     phreg object with deaths
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases recurrent.marginal 
##' @examples
##' \donttest{
##' ### do not test to avoid dependence on frailtypack for data 
##' library(frailtypack)
##' data(readmission)
##' r <-  readmission
##'
##' ## no.opt to fit non-parametric models with just a baseline 
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
##' }
##' @export
recurrent.marginal <- function(recurrent,death,fixbeta=1,...)
{# {{{
  dr <- death 
  xr <- recurrent

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
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
  ###
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazR <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
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
     betaiidR <- iid(x)
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
     betaiidD1 <- iid(x)
     vbeta <- crossprod(betaiidD1)
     varbetat <-   rowSums((HtD1 %*% vbeta)*HtD1)
     ### writing each beta for all individuals 
     betakt <- betaiidD1[id+1,,drop=FALSE]
     ###
     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- c(covk2)*c(cumS0i2)
     ###
     varA2 <- varA2+varbetat-2*apply((covk1-covk2)*Ht,1,sum)
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
  U[xx$jumps+1,] <- x$U
  ###    
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumS0i2 <- cumsumstrata(mu*S0i2,xx$strata,xx$nstrata)
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
     betaiidD <- iid(x)
     vbetaD <- crossprod(betaiidD)
     varbetat <-   rowSums((Ht %*% vbetaD)*Ht)
     ### writing each beta for all individuals 
     betakt <- betaiidD[id+1,,drop=FALSE]
     ###
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
 cov233 <- covfridstrataCov(xxxD,w*rrD,xxxD1,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1)
 cov234 <- covfridstrataCov(xxxD1,w*rrD1,xxxD,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D*mu)
 cov23A <- -c(cov23+(cov232+cov233+cov234))
### cov23A <- cov23A[sort(c(xr$jumps,x$jumps))]
### timeo <- xr$time[c(xr$jumps,x$jumps),]
### timeo <- xr$time[sort(c(xr$jumps,x$jumps)),]
### cbind(cov23A,cov23i)
 cov12 <-  c(cumsumidstratasumCov(xxxR,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)*mu
 cov122 <- c(revcumsumidstratasumCov(w*rrR,w*rrD1,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(mu*cumS0i2R*cumS0i2D1)
 cov123 <- covfridstrataCov(xxxR,w*rrR,xxxD,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1*mu)
 cov124 <- covfridstrataCov(xxxD1,w*rrD1,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2R)
 cov12A <- -c(cov12+(cov122+cov123+cov124))
### cov12A <- cov12A[sort(c(xr$jumps,x$jumps))]
### cbind(cov12A,cov12i)
 cov13 <-  c(cumsumidstratasumCov(xxxR,xxxD,id,mid,xx$strata,xx$nstrata)$sumsquare)
 cov132 <- c(revcumsumidstratasumCov(w*rrR,w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2R*cumS0i2D)
 cov133 <- covfridstrataCov(xxxR,w*rrR,xxxD,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D)
 cov134 <- covfridstrataCov(xxxD,w*rrD,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2R)
###
 cov13A <- c(cov13+(cov132+cov133+cov134))
### cov13A <- cov13A[sort(c(xr$jumps,x$jumps))]
### cbind(cov13A,cov13i)
###
 varA <- varA+2*cov12A+2*cov23A+2*cov13A
# }}}
 if (fixbeta==0) {
### covariances between different terms  beta's 
# {{{
 covbetaRD <- t(betaiiidR) %*% betaiidD
 DHt <- HtD-HtD1*mu
 covbeta1.23 <-   2*rowSums((HtR %*% covbetaRD)*DHt)
 covbetaD.23 <-   2*rowSums((HtD %*% covbetaD)*(HtD1*mu))
 ### R, vesus betaD
 betakt <- betaiidD[id+1,,drop=FALSE]
 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2R)
 covR.23 <- -2*apply((covk1-covk2)*DHt,1,sum)
 ### D versus betaR 
 betakt <- betaiidR[id+1,,drop=FALSE]
 covk1<-apply(xxxD*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2<-apply(w*rrD*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2D)
 covD.1 <- -2*apply((covk1-covk2)*HtD,1,sum)
 covk1<-apply(xxxD1*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
 covk2<-apply(w*rrD1*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
 covk2 <- c(covk2)*c(cumS0i2D1)*mu
 covD.2 <- -2*apply((covk1-covk2)*HtD1,1,sum)

 varA <- varA+ covR.23 + covbeta1.23+covbetaD.23  + covD.1 + covD.2
# }}}
 }


###
 varrs <- data.frame(mu=mu,cumhaz=mu,se.mu=varA^.5,time=xr$time,se.cumhaz=varA^.5,strata=xr$strata)
###		   varA0=varA0, varA1=varA1,varA2=varA2,varA3,
###                   cov12=cov12A, cov23=cov23A,cov13=cov13A,jumps=1:length(mu))
 varrs <- varrs[sort(c(xr$jumps,x$jumps)),]
 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,se.mu=varrs$se.mu,times=varrs$time,
	     cumhaz=cbind(varrs$time,varrs$mu),se.cumhaz=cbind(varrs$time,varrs$se.mu),
	     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),strata.name=xr$strata.name)
 return(out)
}# }}}


