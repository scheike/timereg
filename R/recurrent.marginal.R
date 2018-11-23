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
##' Assumes no ties in the sense that jump times needs to be unique, this is particularly so for the stratified version.
##' 
##' @param recurrent phreg object with recurrent events
##' @param death     phreg object with deaths
##' @param fixbeta   to force the estimation of standard errors to think of regression coefficients as known/fixed
##' @param km  if true then uses Kaplan-Meier for death, otherwise exp(- Nelson-Aalen ) 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' 
##' @references 
##'             Ghosh and Lin (2002) Nonparametric Analysis of Recurrent events and death, 
##'             Biometrics, 554--562.
##' @examples
##'
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##' rr <- simRecurrent(1000,base1,death.cumhaz=dr)
##' rr$x <- rnorm(nrow(rr)) 
##' rr$strata <- floor((rr$id-0.01)/500)
##' 
##' ##  to fit non-parametric models with just a baseline 
##' xr <- phreg(Surv(entry,time,status)~cluster(id),data=rr)
##' dr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)
##' par(mfrow=c(1,3))
##' bplot(dr,se=TRUE)
##' title(main="death")
##' bplot(xr,se=TRUE)
##' ### robust standard errors 
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=4)
##' 
##' ## marginal mean of expected number of recurrent events 
##' out <- recurrentMarginal(xr,dr)
##' bplot(out,se=TRUE,ylab="marginal mean",col=2)
##' 
##' ########################################################################
##' ###   with strata     ##################################################
##' ########################################################################
##' xr <- phreg(Surv(entry,time,status)~strata(strata)+cluster(id),data=rr)
##' dr <- phreg(Surv(entry,time,death)~strata(strata)+cluster(id),data=rr)
##' par(mfrow=c(1,3))
##' bplot(dr,se=TRUE)
##' title(main="death")
##' bplot(xr,se=TRUE)
##' rxr <-   robust.phreg(xr,fixbeta=1)
##' bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)
##'
##' out <- recurrentMarginal(xr,dr)
##' bplot(out,se=TRUE,ylab="marginal mean",col=1:2)
##'
##' ########################################################################
##' ###   cox case        ##################################################
##' ########################################################################
##' xr <- phreg(Surv(entry,time,status)~x+cluster(id),data=rr)
##' dr <- phreg(Surv(entry,time,death)~x+cluster(id),data=rr)
##' par(mfrow=c(1,3))
##' bplot(dr,se=TRUE)
##' title(main="death")
##' bplot(xr,se=TRUE)
##' rxr <-   robust.phreg(xr)
##' bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)
##'
##' out <- recurrentMarginal(xr,dr)
##' bplot(out,se=TRUE,ylab="marginal mean",col=1:2)
##' 
##' ########################################################################
##' ###   CIF  #############################################################
##' ########################################################################
##' ### use of function to compute cumulative incidence (cif) with robust standard errors
##'  data(bmt)
##'  bmt$id <- 1:nrow(bmt)
##'  xr  <- phreg(Surv(time,cause==1)~cluster(id),data=bmt)
##'  dr  <- phreg(Surv(time,cause!=0)~cluster(id),data=bmt)
##' 
##'  out <- recurrentMarginal(xr,dr,km=TRUE)
##'  bplot(out,se=TRUE,ylab="cumulative incidence")
##' 
##' @export
##' @aliases recurrentMarginal tie.breaker  recmarg
recurrentMarginal <- function(recurrent,death,fixbeta=NULL,km=TRUE,...)
{# {{{
  xr <- recurrent
  dr <- death 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if (is.null(xr$opt) | is.null(xr$coef)) fixbeta<- 1 else fixbeta <- 0

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## survival at t- to also work in competing risks situation
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
    St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazR <-  cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
  mu <- cumhazDR[,2]
# }}}

  ### robust standard errors 
  ### 1. sum_k ( int_0^t S(s)/S_0^r(s) dM_k.^r(s) )^2
 resIM1 <-  squareintHdM(xr,ft=St,fixbeta=fixbeta)
 ### 2. mu(t)^2 * sum_k ( int_0^t 1/S_0^d(s) dM_k.^d(s) )^2
 resIM2 <-  squareintHdM(dr,ft=NULL,fixbeta=fixbeta)
  ### 3. sum_k( int_0^t mu(s) /S_0^d(s) dM_k.^d(s))^2
 resIM3 <-  squareintHdM(dr,ft=mu,fixbeta=fixbeta)

 varA <-  resIM1$varInt+mu^2*resIM2$varInt+resIM3$varInt 

## covariances between different terms  13 23  12 12
## to allow different strata for xr and dr, but still nested strata
 if ((xr$nstrata>1 & dr$nstrata==1)) {
    cM1M3 <- covIntH1dM1IntH2dM2(resIM1,resIM3,fixbeta=fixbeta,mu=NULL)
    cM1M2 <- covIntH1dM1IntH2dM2(resIM1,resIM2,fixbeta=fixbeta,mu=mu)
 } else  {
    cM1M3 <- covIntH1dM1IntH2dM2(resIM3,resIM1,fixbeta=fixbeta,mu=NULL)
    cM1M2 <- covIntH1dM1IntH2dM2(resIM2,resIM1,fixbeta=fixbeta,mu=mu)
 }
 cM2M3 <- covIntH1dM1IntH2dM2(resIM2,resIM3,fixbeta=fixbeta,mu=mu)

 varA <- varA+2*cM1M3$cov12A-2*cM1M2$cov12A-2*cM2M3$cov12A 
### varA <- varA-2*cM1M3$cov12A+2*cM1M2$cov12A+2*cM2M3$cov12A 

 cov12aa <- cov13aa <- cov23aa <- 0

 if (fixbeta==0) {
    varA <-varA + cM2M3$covbeta - cM1M3$covbeta + cM1M2$covbeta 
 }

 varrs <- data.frame(mu=mu,cumhaz=mu,se.mu=varA^.5,time=xr$time,
		     se.cumhaz=varA^.5,strata=xr$strata,St=St)
 varrs <- varrs[c(xr$cox.prep$jumps)+1,]

 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,se.mu=varrs$se.mu,times=varrs$time,
     St=varrs$St,
     cumhaz=cbind(varrs$time,varrs$mu),se.cumhaz=cbind(varrs$time,varrs$se.mu),
     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
     strata.name=xr$strata.name,strata.level=recurrent$strata.level)
 return(out)
}# }}}

###recurrentMarginal <- function(recurrent,death,fixbeta=NULL,km=FALSE,...)
###{# {{{
###  xr <- recurrent
###  dr <- death 
###
###  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
###  if (is.null(fixbeta)) 
###  if (is.null(xr$opt) | is.null(xr$coef)) fixbeta<- 1 else fixbeta <- 0
###
###  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
###  # {{{
###  strat <- dr$strata[dr$jumps]
###  ###
###  x <- dr
###  xx <- x$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
###  S0i[xx$jumps+1] <-  1/x$S0
###  S0i2[xx$jumps+1] <- 1/x$S0^2
###  ## survival at t- to also work in competing risks situation
###  if (!km) { 
###    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
###    St      <- exp(-cumhazD)
###  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
###  x <- xr
###  xx <- x$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
###  S0i[xx$jumps+1] <-  1/x$S0
###  S0i2[xx$jumps+1] <- 1/x$S0^2
###  cumhazR <-  cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
###  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
###  mu <- cumhazDR[,2]
#### }}}
###
###  test <- 0
###  ### robust standard errors 
###  ### 1. sum_k ( int_0^t S(s)/S_0^r(s) dM_k.^r(s) )^2
### if (test==1) {
#### {{{
###  x <- xr
###  xx <- x$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
###  S0i[xx$jumps+1] <-  1/x$S0
###  S0i2[xx$jumps+1] <- 1/x$S0^2
###  Z <- xx$X
###  U <- E <- matrix(0,nrow(xx$X),x$p)
###  E[xx$jumps+1,] <- x$E
###  U[xx$jumps+1,] <- x$U
###  ###    
###  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
###  cumS0i2 <- cumsumstrata(St*S0i2,xx$strata,xx$nstrata)
###  if (fixbeta==0) {
###	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
###	  Ht <- apply(St*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
###	  HtR <- Ht
###  }
###  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
###  else rr <- c(xx$sign*exp(xx$offset))
###  id <-   xx$id
###  mid <- max(id)+1
###  ### also weights 
###  w <- c(xx$weights)
###  xxx <- w*(St*S0i-rr*c(cumS0i2))
###  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
###  ss <- c(revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2^2)
###  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2)
###  varA1 <- c(ssf+ss+2*covv)
###  cumS0i2R <- cumS0i2; xxxR <- xxx; rrR <- rr
###
###  if (fixbeta==0) {# {{{
###      invhess <- -solve(x$hessian)
###      MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-EdLam0)*rr*c(xx$weights)
###      UU <- apply(MGt,2,sumstrata,id,max(id)+1)
###      betaiidR <- UU %*% invhess
######     betaiidR <- iid(x)
###     vbeta <- crossprod(betaiidR)
###     varbetat <-   rowSums((Ht %*% vbeta)*Ht)
###     ### writing each beta for all individuals 
###     betakt <- betaiidR[id+1,,drop=FALSE]
###     ###
###     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###     covk2 <- c(covk2)*c(cumS0i2)
###     ###
###     varA1 <- varA1+varbetat-2*apply((covk1-covk2)*Ht,1,sum)
###  }# }}}
#### }}}
### }
### resIM1 <-  squareintHdM(xr,ft=St,fixbeta=fixbeta,...)
###
###
### ### 2. mu(t)^2 * sum_k ( int_0^t 1/S_0^d(s) dM_k.^d(s) )^2
### resIM2 <-  squareintHdM(dr,ft=NULL,fixbeta=fixbeta,...)
###
### if (test==1) {
#### {{{
###  x <- dr
###  xx <- x$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
###  S0i[xx$jumps+1] <-  1/x$S0
###  S0i2[xx$jumps+1] <- 1/x$S0^2
###  if (fixbeta==0) {
###	  Z <- xx$X
###	  U <- E <- matrix(0,nrow(xx$X),x$p)
###	  E[xx$jumps+1,] <- x$E
###	  U[xx$jumps+1,] <- x$U
###	  Ht <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
###  }
###  ###    
###  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
###  cumS0i2 <- cumsumstrata(S0i2,xx$strata,xx$nstrata)
###  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
###  else rr <- c(xx$sign*exp(xx$offset))
###  id <-   xx$id
###  mid <- max(id)+1
###  ### also weights 
###  w <- c(xx$weights)
###  xxx1 <- w*(S0i-rr*c(cumS0i2))
###  ssf <- cumsumidstratasum(xxx1,id,mid,xx$strata,xx$nstrata)$sumsquare
###  ss <- c(revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2^2)
###  covv <- covfridstrata(xxx1,w*rr,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2)
###  varA2 <- mu^2*c(ssf+ss+2*covv)
###  cumS0i2D1 <- cumS0i2 
###  xxxD1 <- xxx1
###  rrD1 <- rr
###  if (fixbeta==0) {# {{{
###      HtD1 <- mu*Ht
###      invhess <- -solve(x$hessian)
###      MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-Ht)*rr*c(xx$weights)
###      UU <- apply(MGt,2,sumstrata,id,max(id)+1)
###      betaiidD <- UU %*% invhess
######     betaiidD1 <- iid(x)
###     vbetaD <- crossprod(betaiidD)
###     varbetat <-   rowSums((HtD1 %*% vbetaD)*HtD1)
###     ### writing each beta for all individuals 
###     betakt <- betaiidD[id+1,,drop=FALSE]
###     ###
###     covk1 <- apply(xxx1*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###     covk2 <- c(covk2)*c(cumS0i2D1)
###     ###
###     varA2 <- varA2+varbetat-2*apply((covk1-covk2)*HtD1*mu,1,sum)
###  }# }}}
###
#### }}}
### }
###
###  ### 3. sum_k( int_0^t mu(s) /S_0^d(s) dM_k.^d(s))^2
### if (test==1) {
#### {{{
###  x <- dr
###  xx <- x$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
###  S0i[xx$jumps+1] <-  1/x$S0
###  S0i2[xx$jumps+1] <- 1/x$S0^2
###  if (fixbeta==0) {
###	  Z <- xx$X
###	  U <- E <- matrix(0,nrow(xx$X),x$p)
###	  E[xx$jumps+1,] <- x$E
###  }
######  U[xx$jumps+1,] <- x$U
###  xr$cox.prep$time  - xx$time
###  xr$cox.prep$time[xr$cox.prep$jumps+1]
###  xr$cox.prep$time[dr$cox.prep$jumps+1]
###  ###    
###  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
###  cumS0i2 <- cumsumstrata(mu*S0i2,xx$strata,xx$nstrata)
###
###  if (fixbeta==0) {
###	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
###	  Ht <- apply(mu*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
###	  HtD <- Ht
###  }
###  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
###  else rr <- c(xx$sign*exp(xx$offset))
###  id <-   xx$id
###  mid <- max(id)+1
###  ### also weights 
###  w <- c(xx$weights)
###  xxx <- w*(mu*S0i-rr*c(cumS0i2))
###  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
###  ss <- c(revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2^2)
###  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2)
###  varA3 <- c(ssf+ss+2*covv)
###  cumS0i2D <- cumS0i2
###  rrD <- rr
###  xxxD <- xxx
### 
###  if (fixbeta==0) {# {{{
###     invhess <- -solve(x$hessian)
###     varbetat <-   rowSums((Ht %*% vbetaD)*Ht)
###     ### writing each beta for all individuals 
###     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###     covk2 <- c(covk2)*c(cumS0i2)
###     ###
###     varA3 <- varA3+varbetat-2*apply((covk1-covk2)*Ht,1,sum)
###  }# }}}
###  varA <-  varA1+varA2+varA3 
#### }}}
### }
### resIM3 <-  squareintHdM(dr,ft=mu,fixbeta=fixbeta,...)
### varA <-  resIM1$varInt+mu^2*resIM2$varInt+resIM3$varInt 
###
### if (test==1) {# {{{
###	 print("=====var ================="); 
###	 print(summary(varA1))
###	 print(summary(varA2)); 
###	 print(summary(varA3)); 
###	 print("--------------------------"); 
###	 print(summary(resIM1$varInt))
###	 print(summary(mu^2*resIM2$varInt))
###	 print(summary(resIM3$varInt))
###	 print("======================"); 
### }# }}}
###
### ### covariances between different terms  13 23  12 12
### if (test==1) {
### # {{{
### cov23<-c(cumsumidstratasumCov(xxxD,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)
### cov232<-c(revcumsumidstratasum(w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2D*cumS0i2D1)
###### cov232<-c(revcumsumidstratasumCov(w*rrD,w*w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2D*cumS0i2D1)
### cov233 <- covfridstrata(xxxD,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1)
###### cov233 <- covfridstrataCov(xxxD,w*rrD1,xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1)
### cov234 <- covfridstrata(xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D)
###### cov234 <- covfridstrataCov(xxxD1,w*rrD1,xxxD,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D)
### cov23A <- -c(cov23+(cov232+cov233+cov234))*mu
###### print(summary(cov23))
###### print(summary(cov232))
###### print(summary(cov233))
###### print(summary(cov234))
###### print(" ====================== 23 slut ")
######
###
### cov12 <-  c(cumsumidstratasumCov(xxxR,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)
### cov122 <- c(revcumsumidstratasumCov(w*rrR,w*rrD1,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2R*cumS0i2D1)
### cov123 <- covfridstrataCov(xxxR,w*rrR,xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1)
### cov124 <- covfridstrataCov(xxxD1,w*rrD1,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2R)
### cov12A <- -c(cov12+cov122+cov123+cov124)*mu
### print("________________12____________________"); 
### print(summary(c(xxxR))); print(summary(c(xxxD1))); 
### print(summary(c(rrD1))); print(summary(c(rrR)))
### print(summary(c(cumS0i2R))); print(summary(c(cumS0i2D1)))
### print(summary(cov12)); 
### print(summary(cov122)); 
### print(summary(c(cov123))); 
### print(summary(c(cov124))); 
### print("______________________________________"); 
###
######
### cov13 <-  c(cumsumidstratasumCov(xxxR,xxxD,id,mid,xx$strata,xx$nstrata)$sumsquare)
### cov132 <- c(revcumsumidstratasumCov(w*rrR,w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2R*cumS0i2D)
### cov133 <- covfridstrataCov(xxxR,w*rrR,xxxD,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D)
### cov134 <- covfridstrataCov(xxxD,w*rrD,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2R)
######
### cov13A <- c(cov13+(cov132+cov133+cov134))
### varAg <- varA+2*cov12A+2*cov23A+2*cov13A
#### }}}
### }
###
###
###cM1M3 <- covIntH1dM1IntH2dM2(resIM3,resIM1,fixbeta=fixbeta,mu=NULL)
###cM1M2 <- covIntH1dM1IntH2dM2(resIM2,resIM1,fixbeta=fixbeta,mu=mu)
###cM2M3 <- covIntH1dM1IntH2dM2(resIM2,resIM3,fixbeta=fixbeta,mu=mu)
###
######print(lapply(cM1M3,summary))
######cM1M3 <- covIntH1dM1IntH2dM2(resIM1,resIM3,fixbeta=fixbeta,mu=NULL)
######print(lapply(cM1M3,summary))
######print("____________________________")
######print(lapply(cM1M2,summary))
######cM1M2 <- covIntH1dM1IntH2dM2(resIM1,resIM2,fixbeta=fixbeta,mu=mu)
######print(lapply(cM1M2,summary))
######print("____________________________")
######
######print(lapply(cM2M3,summary))
######cM2M3 <- covIntH1dM1IntH2dM2(resIM2,resIM3,fixbeta=fixbeta,mu=mu)
######print(lapply(cM2M3,summary))
######print("____________________________")
###
###
### if (test==1) {# {{{
###	 print("=======cov AAA==============="); 
###	 print(summary(cov13A))
###	 print(summary(cov12A))
###	 print(summary(cov23A))
###	 print("_______________________________"); 
###	 print(summary(cM1M3$cov12A))
###	 print(summary(-1*cM1M2$cov12A))
###	 print(summary(-1*cM2M3$cov12A))
###	 print("======================"); 
### }# }}}
### varA <- varA+2*cM1M3$cov12A-2*cM1M2$cov12A-2*cM2M3$cov12A 
###
### if (test==1) {# {{{
### print(" var Ag"); 
### print(summary(varAg)); 
### print(" var A"); 
### print(summary(varA)); 
### }# }}}
###
### cov12aa <- cov13aa <- cov23aa <- 0
### if (fixbeta==0 & test==1 ) {
### ### covariances between different terms and  beta's 
###	# {{{
###	 covbetaRD <- t(betaiidR) %*% betaiidD
###	 DHt <- HtD1-HtD
###	### covbeta1.23 <-   -2*rowSums((HtR %*% covbetaRD)*DHt)
###	 covbeta1.12 <-   -2*rowSums((HtR %*% covbetaRD)*HtD1)
###	 covbeta1.13 <-   2*rowSums((HtR %*% covbetaRD)*HtD)
###	 covbetaD.23 <-   -2*rowSums((HtD %*% vbetaD)*(HtD1))
###
###	 ### D versus betaD from two terms  cov23 wrt beta 
###	 betakt <- betaiidD[id+1,,drop=FALSE]
###	 covk1 <-apply(xxxD1*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###	 covk2 <-apply(w*rrD1*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###	 covk2 <- c(covk2)*c(cumS0i2D1)
###	 covD1.D <- 2*apply((covk1-covk2)*mu*HtD,1,sum)
###	 ###
###	 covk1 <-apply(xxxD*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###	 covk2 <-apply(w*rrD*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###	 covk2 <- c(covk2)*c(cumS0i2D)
###	 covD.D1 <- 2*apply((covk1-covk2)*HtD1,1,sum)
###	 cov23aa <- covbetaD.23  + covD1.D + covD.D1
###
###	 ### cov12 wrt betaD and betaR
###	 betakt <- betaiidD[id+1,,drop=FALSE]
###	 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###	 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###	 covk2 <- c(covk2)*c(cumS0i2R)
###	 covRD12 <- 2*apply((covk1-covk2)*HtD1,1,sum)
###	###
###	 betakt <- betaiidR[id+1,,drop=FALSE]
###	 covk1 <-apply(xxxD1*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###	 covk2 <-apply(w*rrD1*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###	 covk2 <- c(covk2)*c(cumS0i2D1)
###	 covRD21 <- 2*apply((covk1-covk2)*mu*HtR,1,sum)
###	 cov12aa <- covbeta1.12 + covRD12+covRD21
######         print("--------12------")
######	 print(summary(covbeta1.12))
######	 print(summary(covRD12))
######	 print(summary(covRD21))
######	 print("--------------")
###
###
###	 ### cov13 wrt betaD and betaR
###	 betakt <- betaiidD[id+1,,drop=FALSE]
###	 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###	 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###	 covk2 <- c(covk2)*c(cumS0i2R)
###	 covRD13 <- -2*apply((covk1-covk2)*HtD,1,sum)
###	###
###	 betakt <- betaiidR[id+1,,drop=FALSE]
###	 covk1 <-apply(xxxD*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
###	 covk2 <-apply(w*rrD*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
###	 covk2 <- c(covk2)*c(cumS0i2D)
###	 covRD31 <- -2*apply((covk1-covk2)*HtR,1,sum)
###	 cov13aa <- covbeta1.13 + covRD13+covRD31
### 	 varAg <-varA+ cov23aa+cov13aa+cov12aa
#### }}}
### }
###
### if (test==1) {
### print("=========cov beta============"); 
###	 print(summary(cov13aa))
###	 print(summary(cov12aa))
###	 print(summary(cov23aa))
### print("_______________________________"); 
### print(summary(-1*cM1M3$covbeta)); 
### print(summary(cM1M2$covbeta))
### print(summary(cM2M3$covbeta));
### print("======================"); 
### }
###
### if (fixbeta==0) {
###    varA <-varA+ cM2M3$covbeta - cM1M3$covbeta + cM1M2$covbeta 
### if (test==1) {
### print(summary(varAg))
### print(summary(varA))
### }
### }
###
### varrs <- data.frame(mu=mu,cumhaz=mu,se.mu=varA^.5,time=xr$time,
###		     se.cumhaz=varA^.5,strata=xr$strata,St=St)
### varrs <- varrs[c(xr$cox.prep$jumps)+1,]
###
### ### to use basehazplot.phreg
### ### making output such that basehazplot can work also
### out <- list(mu=varrs$mu,se.mu=varrs$se.mu,times=varrs$time,
###     St=varrs$St,
###     cumhaz=cbind(varrs$time,varrs$mu),se.cumhaz=cbind(varrs$time,varrs$se.mu),
###     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
###     strata.name=xr$strata.name,strata.level=recurrent$strata.level)
### return(out)
###}# }}}

##' @export
recurrentMarginalgam <- function(recurrent,death,fixbeta=NULL,km=TRUE,...)
{# {{{
  xr <- recurrent
  dr <- death 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if (is.null(xr$opt) | is.null(xr$coef)) fixbeta<- 1 else fixbeta <- 0

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  ## survival at t- to also work in competing risks situation
  if (!km) { 
    cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
    St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazR <-  cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
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
  cumS0i2R <- cumS0i2; xxxR <- xxx; rrR <- rr

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
 cov23<-c(cumsumidstratasumCov(xxxD,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)
 cov232<-c(revcumsumidstratasum(w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2D*cumS0i2D1)
### cov232<-c(revcumsumidstratasumCov(w*rrD,w*w*rrD,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2D*cumS0i2D1)
 cov233 <- covfridstrata(xxxD,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1)
### cov233 <- covfridstrataCov(xxxD,w*rrD1,xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1)
 cov234 <- covfridstrata(xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D)
### cov234 <- covfridstrataCov(xxxD1,w*rrD1,xxxD,w*rrD,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D)
 cov23A <- -c(cov23+(cov232+cov233+cov234))*mu
### print(summary(cov23))
### print(summary(cov232))
### print(summary(cov233))
### print(summary(cov234))
### print(" ====================== 23 slut ")
###

 cov12 <-  c(cumsumidstratasumCov(xxxR,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)
 cov122 <- c(revcumsumidstratasumCov(w*rrR,w*rrD1,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2R*cumS0i2D1)
 cov123 <- covfridstrataCov(xxxR,w*rrR,xxxD,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2D1)
 cov124 <- covfridstrataCov(xxxD1,w*rrD1,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2R)
 cov12A <- -c(cov12+(cov122+cov123+cov124))*mu

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
###	 DHt <- HtD1-HtD
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
	 cov23aa <- (covbetaD.23  + covD1.D + covD.D1)
###	 print("D1D_________________")
###	 print(summary(covbetaD.23))
###	 print(summary(covD1.D))
###	 print(summary(covD.D1))

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
###	 print("RD12_________________")
###	 print(summary(covbeta1.12))
###	 print(summary(covRD12))
###	 print(summary(covRD21))


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
###	 print("RD13_________________")
###	 print(summary(covbeta1.13))
###	 print(summary(covRD13))
###	 print(summary(covRD31))

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
     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),strata.name=xr$strata.name,
     strata.level=recurrent$strata.level)
###  vari=varrs[,c("varA1","varA2","varA3")],covs=varrs[,c("cov12","cov13","cov23")])
 return(out)
}# }}}

##' @export
recmarg <- function(recurrent,death,Xr=NULL,Xd=NULL,km=TRUE,...)
{# {{{
  xr <- recurrent
  dr <- death 

  if (!is.null(Xr)) rr <- exp(sum(xr$coef * Xr)) else rr <- 1
  if (!is.null(Xd)) rrd <- exp(sum(dr$coef * Xd)) else rrd <- 1

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  strat <- dr$strata[dr$jumps]
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD*rrd)
  } else St <- exp(rrd*c(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  ###
  x <- xr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  cumhazDR <- cbind(xx$time,cumsumstrata(St*S0i,xx$strata,xx$nstrata))
  mu <- rr*cumhazDR[,2]
# }}}

 varrs <- data.frame(mu=mu,time=xr$time,strata=xr$strata,St=St)
 varrs  <-  varrs[c(xr$cox.prep$jumps)+1,]

 ### to use basehazplot.phreg
 ### making output such that basehazplot can work also
 out <- list(mu=varrs$mu,time=varrs$time,
	     St=varrs$St,cumhaz=cbind(varrs$time,varrs$mu),
             strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
	     strata.name=xr$strata.name,strata.level=recurrent$strata.level)
 return(out)
}# }}}

##' @export
squareintHdM <- function(phreg,ft=NULL,fixbeta=NULL,...)
{# {{{
###  sum_k ( int_0^t f(s)/S_0^r(s) dM_k.^r(s) )^2
###  strata "r" from object and "k" id from cluster 
  if (class(phreg)!="phreg") stop("Must be phreg object\n"); 

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if (is.null(phreg$opt) | is.null(phreg$coef)) fixbeta<- 1 else fixbeta <- 0

  x <- phreg
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
  if (is.null(ft))  ft <- rep(1,length(xx$time))
  cumS0i2 <- c(cumsumstrata(ft*S0i2,xx$strata,xx$nstrata))
  if (fixbeta==0) {
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  Ht <- apply(ft*E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  } else Ht <- NULL
  if (fixbeta==0) rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
  else rr <- c(xx$sign*exp(xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  ### also weights 
  w <- c(xx$weights)
  xxx <- (ft*S0i-rr*cumS0i2)
  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
  ss <-  revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare*cumS0i2^2
  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*cumS0i2
  varA1 <- c(ssf+ss-2*covv)

  vbeta <- betaiidR <- NULL
  if (fixbeta==0) {# {{{
     invhess <- -solve(x$hessian)
     MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-EdLam0)*rr*c(xx$weights)
     UU <- apply(MGt,2,sumstrata,id,mid)
     betaiidR <- UU %*% invhess
     vbeta <- crossprod(betaiidR)
     varbetat <-   rowSums((Ht %*% vbeta)*Ht)
     ### writing each beta for all individuals 
     betakt <- betaiidR[id+1,,drop=FALSE]
     ###
     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- c(covk2)*cumS0i2
     covv <- covk1-covk2
     ###
     varA1 <- varA1+varbetat-2*apply(covv*Ht,1,sum)
  }# }}}

  return(list(xx=xx,Ht=Ht,varInt=varA1,xxx=xxx,rr=rr,
	      cumhaz=cumhaz,cumS0i2=cumS0i2,mid=mid,id=id,
	      betaiid=betaiidR,vbeta=vbeta,covv=covv))
} # }}}

##' @export
covIntH1dM1IntH2dM2 <- function(square1,square2,fixbeta=1,mu=NULL)
{# {{{

 ### strata and id same for two objects 
 xx <- square1$xx; xx2 <- square2$xx
 xxxR <- square1$xxx;     xxxD1 <- square2$xxx
 rrR  <- square1$rr;       rrD1 <- square2$rr
 id   <- id1 <- square1$id; id2 <- square2$id
 mid  <- square1$mid; w <- c(xx$weights)

 if (is.null(mu)) mu <- rep(1,length(xx$strata))

 cov12 <-  c(cumsumidstratasumCov(xxxR,xxxD1,id,mid,xx$strata,xx$nstrata)$sumsquare)
 cov122 <- c(revcumsumidstratasumCov(w*rrR,w*rrD1,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(square1$cumS0i2*square2$cumS0i2)
 cov123 <- covfridstrataCov(xxxR,w*rrR,xxxD1,w*rrD1,id,mid,xx$strata,xx$nstrata)$covs*c(square2$cumS0i2)
 cov124 <- covfridstrataCov(xxxD1,w*rrD1,xxxR,w*rrR,id,mid,xx$strata,xx$nstrata)$covs*c(square1$cumS0i2)
 cov12A <- c(cov12+cov122+cov123+cov124)
 
 test <- 0
 if (test==1) {# {{{
	 print("________cov cov ___________________________"); 
	 print(summary(c(xxxR))); print(summary(c(xxxD1))); 
	 print(summary(c(rrD1))); print(summary(c(rrR)))
	 print(summary(c(square1$cumS0i2))); print(summary(c(square2$cumS0i2)));
	 print("-----------"); 
	 print(summary(cov12)); 
	 print(summary(cov122)); 
	 print(summary(c(cov123))); 
	 print(summary(c(cov124))); 
	 print("______________________________________"); 
	 jumps <- c(square1$xx$jumps,square2$xx$jumps)+1
	 print(summary(jumps))
         print(summary(cov12[jumps])); 
	 print(summary(cov122[jumps])); 
	 print(summary(c(cov123[jumps]))); 
	 print(summary(c(cov124[jumps]))); 
 }# }}}

 cov12aa <- 0
 if (fixbeta==0) {
 ### covariances between different terms and  beta's 
 # {{{
	 betaiidR <- square1$betaiid; betaiidD <- square2$betaiid
	 HtR <- square1$Ht; HtD <- square2$Ht
	 covbetaRD <- t(betaiidR) %*% betaiidD
	 covbeta <-   -1*rowSums((HtR %*% covbetaRD)*HtD)
###	 print(summary(covbeta))

	 ### cov12 wrt betaD and betaR
	 betakt <- betaiidD[id1+1,,drop=FALSE]
	 covk1 <-apply(xxxR*betakt,2,cumsumidstratasum,id1,mid,xx$strata,xx$nstrata,type="sum")
	 covk2 <-apply(w*rrR*betakt,2,revcumsumidstratasum,id1,mid,xx$strata,xx$nstrata,type="lagsum")
	 covk2 <- c(covk2)*c(square1$cumS0i2)
	 covRD12 <- apply((covk1-covk2)*HtD,1,sum)
	###
	 betakt <- betaiidR[id2+1,,drop=FALSE]
	 covk1 <-apply(xxxD1*betakt,2,cumsumidstratasum,id2,mid,xx2$strata,xx$nstrata,type="sum")
	 covk2 <-apply(w*rrD1*betakt,2,revcumsumidstratasum,id2,mid,xx2$strata,xx$nstrata,type="lagsum")
	 covk2 <- c(covk2)*c(square2$cumS0i2)
	 covRD21 <- apply((covk1-covk2)*HtR,1,sum)
	 cov12aa <- 2*(covbeta + covRD12+covRD21)
	 test <- 0
	 if (test==1) {
		 print("--------------")
		 print(summary(2*mu*covbeta))
		 print(summary(2*mu*covRD12))
		 print(summary(2*mu*covRD21))
		 print("--------------")
	 }
 } # }}}

 cov12 <- (cov12A-cov12aa)*mu

 return(list(cov=cov12,cov12A=cov12A*mu,covbeta=cov12aa*mu))
} # }}}

##' @export
tie.breaker <- function(data,stop="time",start="entry",status="status",id=NULL,ddt=NULL,exit.unique=TRUE)
{# {{{

   if (!is.null(id)) id <- data[,id]
   ord <- 1:nrow(data)
   stat <- data[,status]
   time <- data[,stop]
   dupexit <- duplicated(time)
   time1 <- data[stat==1,stop]
   time0 <- data[stat!=1,stop]
   lt0 <- length(time0)
   ddp <- duplicated(c(time0,time1))
   if (exit.unique) ties <-ddp[(lt0+1):nrow(data)] else ties <- duplicated(c(time1))
   nties <- sum(ties)
   ordties <- ord[stat==1][ties]
   if (is.null(ddt)) {
	   abd <- abs(diff(data[,stop]))
	   abd <- min(abd[abd>0])
	   ddt <- abd*0.5
   }
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
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param death.cumhaz cumulative hazard of death 
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
##'  ######################################################################
##'  ### simulating simple model that mimicks data 
##'  ######################################################################
##'  rr <- simRecurrent(5,base1,death.cumhaz=dr)
##'  dlist(rr,.~id,n=0)
##'
##'  rr <- simRecurrent(1000,base1,death.cumhaz=dr)
##'  par(mfrow=c(1,3))
##'  showfitsim(causes=1,rr,dr,base1,base1)
##'
##' ######################################################################
##' ### simulating simple model that mimicks data 
##' ### now with two event types and second type has same rate as death rate
##' ######################################################################
##'
##'  rr <- simRecurrent(1000,base1,death.cumhaz=dr,cumhaz2=base4)
##'  dtable(rr,~death+status)
##'  par(mfrow=c(2,2))
##'  showfitsim(causes=2,rr,dr,base1,base4)
##'
##' ######################################################################
##' ### simulating simple model 
##' ### random effect for all causes (Z shared for death and recurrent) 
##' ######################################################################
##'
##'  rr <- simRecurrent(1000,base1,
##'         death.cumhaz=dr,dependence=1,var.gamma=0.4)
##'  ### marginals do fit after input after integrating out
##'  par(mfrow=c(2,2))
##'  showfitsim(causes=1,rr,dr,base1,base1)
##'
##' @export
##' @aliases simRecurrent showfitsim  simRecurrentGamma covIntH1dM1IntH2dM2 recurrentMarginalgam squareintHdM addCums
simRecurrent <- function(n,cumhaz,death.cumhaz=NULL,cumhaz2=NULL,
			    gap.time=FALSE,
			    max.recurrent=100,dhaz=NULL,haz2=NULL,
			    dependence=0,var.z=2,cor.mat=NULL,...) 
  {# {{{
  dtime <- NULL ## to avoid R-check 

  if (dependence==0) { z1 <- z2 <- zd <- rep(1,n) } else if (dependence==1) {# {{{
###	      zz <- rgamma(n,1/var.gamma[1])*var.gamma[1]
	      zz <- exp(rnorm(n,1)*var.z[1]^.5)
	      var(zz)
	      z1 <- zz; z2 <- zz; zd <- zz
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- (z%*% chol(covv))
	      z1 <- exp(z[,1]); zd <- exp(z[,3])
	      apply(exp(z),2,mean); cov(exp(z))
      } # }}}

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
  tall <- timereg::pc.hazard(cumhaz,z1)
  tall$id <- 1:n
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- timereg::pc.hazard(cumhazd,zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
  } else { tall$dtime <- max.time; tall$fdeath <- 0; cumhazd <- NULL }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=1,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  nrr <- n
  i <- 1; 
  while (any(tt$time<tt$dtime) & i < max.recurrent) {
	  i <- i+1
	  still <- subset(tt,time<dtime)
          tt <- timereg::pc.hazard(cumhaz,z1[still$id],entry=still$time)
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=1,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt,row.names=NULL)
	  nrr <- nrr+nt
  }
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
      haz1t <- timereg::Cpred(cbind(cumhaz1[-ll1,1],haz1),tall$time)[,2]
      haz2t <- timereg::Cpred(cbind(cumhaz2[-ll2,1],haz2),tall$time)[,2]
      p2t <- haz2t/(haz1t+haz2t)
      tall$p2t <- p2t
      tall$status <- (1+rbinom(nrow(tall),1,p2t))*(tall$status>=1)
  }# }}}

  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2

  return(tall)
  }# }}}


lin.approx <- function(x2,xfx,x=1) {# {{{
   ### x=1   gives  f(x2) 
   ### x=-1  gives  f^-1(x2) 
   breaks <- xfx[,x]
   fx     <- xfx[,-x]
   ri <- sindex.prodlim(breaks,x2)
   rrr <- (x2-breaks[ri])/(breaks[ri+1]-breaks[ri])
   res <- rrr*(fx[ri+1]-fx[ri])+fx[ri]
   res[is.na(res)] <- tail(fx,1)
   return(res)
}#

##' @export
addCums <- function(cumB,cumA,max=5)
{# {{{
 times <- sort(unique(c(cumB[,1],cumA[,1])))
 times <- times[times<max]
 cumBjx <- lin.approx(times,cumB,x=1)
 cumAjx <- lin.approx(times,cumA,x=1)
 cumBA <- cumBjx+cumAjx
 return(cbind(times,cumBA))
}# }}}

##' @export
simRecurrentGamma <- function(n,haz=0.5,death.haz=0.1,haz2=0.1,max.recurrent=100,var.z=2,times=5000) 
{# {{{

  dtime <- NULL ## to avoid R-check 

  max.time <- times
  cumhaz1 <- rbind(c(0,0),c(times,times*haz))
  cumhaz2 <- rbind(c(0,0),c(times,times*haz2))
  death.cumhaz <- rbind(c(0,0),c(times,death.haz))
  z <- rgamma(1/var.z)*var.z

  cumhaz <- cbind(times,cumhaz1+cumhaz2)

### recurrent first time
  tall <- timereg::pc.hazard(cumhaz,z)
  tall$id <- 1:n
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- timereg::pc.hazard(cumhazd,n)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
  } else { tall$dtime <- max.time; tall$fdeath <- 0; cumhazd <- NULL }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=1,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tt <- tall
  i <- 1; 
  while (any(tt$time<tt$dtime) & i < max.recurrent) {
	  i <- i+1
	  still <- subset(tt,time<dtime)
          tt <- timereg::pc.hazard(cumhaz,z[still$id],entry=still$time)
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=1,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt,row.names=NULL)
  }
  dsort(tall) <- ~id+entry+time

  ### cause 2 is there then decide if jump is 1 or 2
  if (!is.null(haz2)) {# {{{
      p2t <- haz2/(haz+haz2)
      tall$p2t <- p2t
      tall$status <- (1+rbinom(nrow(tall),1,p2t))*(tall$status>=1)
  }# }}}

  tall$start <- tall$entry
  tall$stop  <- tall$time
  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  ### haz*haz2*(var.z+1)

  return(tall)
}# }}}

##' Simulation of recurrent events data based on cumulative hazards II 
##'
##' Simulation of recurrent events data based on cumulative hazards 
##'
##' Must give hazard of death and two recurrent events.  Possible with two
##' event types and their dependence can be specified but the two recurrent events need
##' to share random effect. Based on drawing the from cumhaz and cumhaz2 and 
##' taking the first event rather
##' the cumulative and then distributing it out. Key advantage of this is that 
##' there is  more flexibility wrt random effects 
##'
##' @param n number of id's 
##' @param cumhaz  cumulative hazard of recurrent events 
##' @param cumhaz2  cumulative hazard of recurrent events  of type 2
##' @param death.cumhaz cumulative hazard of death 
##' @param gap.time if true simulates gap-times with specified cumulative hazard
##' @param max.recurrent limits number recurrent events to 100
##' @param dhaz rate for death hazard if it is extended to time-range of first event 
##' @param haz2 rate of second cause  if it is extended to time-range of first event 
##' @param dependence 0:independence; 1:all share same random effect with variance var.z; 2:random effect exp(normal) with correlation structure from cor.mat; 3:additive gamma distributed random effects, z1= (z11+ z12)/2 such that mean is 1 , z2= (z11^cor.mat(1,2)+ z13)/2, z3= (z12^(cor.mat(2,3)+z13^cor.mat(1,3))/2, with z11 z12 z13 are gamma with mean and variance 1 , first random effect is z1 and for N1 second random effect is z2 and for N2 third random effect is for death  
##' @param var.z variance of random effects 
##' @param cor.mat correlation matrix for var.z variance of random effects 
##' @param cens rate of censoring exponential distribution
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
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
##' rr <- simRecurrentII(1000,base1,base4,death.cumhaz=dr)
##' dtable(rr,~death+status)
##' par(mfrow=c(2,2))
##' showfitsim(causes=2,rr,dr,base1,base4)
##'
##' @aliases simRecurrentIII 
##' @export
simRecurrentII <- function(n,cumhaz,cumhaz2,death.cumhaz=NULL,
		    gap.time=FALSE,max.recurrent=100,dhaz=NULL,haz2=NULL,
		    dependence=0,var.z=0.22,cor.mat=NULL,cens=NULL,...) 
  {# {{{

  fdeath <- dtime <- NULL # to avoid R-check 

  if (dependence==0) { z <- z1 <- z2 <- zd <- rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
###	      z <- exp(rnorm(n,1)*var.z[1]^.5)
	      z1 <- z; z2 <- z; zd <- z
	      if (!is.null(cor.mat)) { zd <- rep(1,n); }
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
###	      print(summary(z))
###	      print(cor(z))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
###	      print(summary(z))
###	      print(cor(z))
      } else stop("dependence 0-3"); # }}}

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
  tall1 <- timereg::pc.hazard(cumhaz,z1)
  tall2 <- timereg::pc.hazard(cumhaz2,z2)
  tall <- tall1 
  tall$status <- ifelse(tall1$time<tall2$time,tall1$status,2*tall2$status)
  tall$time <- ifelse(tall1$time<tall2$time,tall1$time,tall2$time)
  tall$id <- 1:n
  tall$rr2 <- tall2$rr
### death time simulated
  if (!is.null(death.cumhaz)) {
	  timed   <- timereg::pc.hazard(cumhazd,zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
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
          tt1 <- timereg::pc.hazard(cumhaz,z1[still$id],entry=still$time)
          tt2 <- timereg::pc.hazard(cumhaz2,z2[still$id],entry=still$time)
	  tt <- tt1
###          drename(tt1,paste(names(tt1),"1",sep="")) <- ~.
###          drename(tt2,paste(names(tt2),"2",sep="")) <- ~.
          tt$status <- ifelse(tt1$time<=tt2$time,tt1$status,2*tt2$status)
          tt$time <-   ifelse(tt1$time<=tt2$time,tt1$time,tt2$time)
	  tt$rr2 <- tt2$rr
###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  nt <- nrow(tt)
###	  gemsim[(nrr+1):(nrr+nt),] <- tt
	  tall <- rbind(tall,tt[1:nn,],row.names=NULL)
###	  nrr <- nrr+nt
  }
###  tall <- gemsim[1:nrr,]
  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"z") <- z

  return(tall)
  }# }}}

##' @export
showfitsim <- function(causes=2,rr,dr,base1,base4) 
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
  meanr1 <-   recurrentMarginal(xrr,drr)
  basehazplot.phreg(meanr1,se=TRUE)
  if (causes>=2) {
	  meanr2 <-   recurrentMarginal(xrr2,drr)
	  basehazplot.phreg(meanr2,se=TRUE,add=TRUE,col=2)
  }
}# }}}

##' @export
simRecurrentIII <- function(n,cumhaz,cumhaz2,beta1,beta2,death.cumhaz=NULL,
		    gap.time=FALSE,max.recurrent=100,dhaz=NULL,haz2=NULL,
		    dependence=0,var.z=0.22,cor.mat=NULL,cens=NULL,...) 
{# {{{

  fdeath <- dtime <- NULL # to avoid R-check 

  if (dependence==0) { z <- z1 <- z2 <- zd <- rep(1,n) # {{{
     } else if (dependence==1) {
	      z <- rgamma(n,1/var.z[1])*var.z[1]
###	      z <- exp(rnorm(n,1)*var.z[1]^.5)
	      z1 <- z; z2 <- z; zd <- z
	      if (!is.null(cor.mat)) { zd <- rep(1,n); }
      } else if (dependence==2) {
              stdevs <- var.z^.5
              b <- stdevs %*% t(stdevs)  
              covv  <- b * cor.mat  
	      z <- matrix(rnorm(3*n),n,3)
	      z <- exp(z%*% chol(covv))
###	      print(summary(z))
###	      print(cor(z))
	      z1 <- z[,1]; z2 <- z[,2]; zd <- z[,3]; 
      } else if (dependence==3) {
	      z <- matrix(rgamma(3*n,1),n,3)
              z1 <- (z[,1]^cor.mat[1,1]+z[,2]^cor.mat[1,2]+z[,3]^cor.mat[1,3])
              z2 <- (z[,1]^cor.mat[2,1]+z[,2]^cor.mat[2,2]+z[,3]^cor.mat[2,3])
              zd <- (z[,1]^cor.mat[3,1]+z[,2]^cor.mat[3,2]+z[,3]^cor.mat[3,3])
	      z <- cbind(z1,z2,zd)
###	      print(summary(z))
###	      print(cor(z))
      } else stop("dependence 0-3"); # }}}

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
  tall1 <- timereg::pc.hazard(cumhaz,z1)
  tall2 <- timereg::pc.hazard(cumhaz2,z2)
  tall <- tall1 
  tall$status <- ifelse(tall1$time<tall2$time,tall1$status,2*tall2$status)
  tall$time <- ifelse(tall1$time<tall2$time,tall1$time,tall2$time)
  tall$id <- 1:n
  tall$rr2 <- tall2$rr
### death time simulated
  if (!is.null(death.cumhaz)) {# {{{
	  timed   <- timereg::pc.hazard(cumhazd,zd)
	  tall$dtime <- timed$time
	  tall$fdeath <- timed$status
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  } else { 
	  tall$dtime <- max.time; 
	  tall$fdeath <- 0; 
	  cumhazd <- NULL 
	  if (!is.null(cens)) { 
             ctime <- rexp(n)/cens
	     tall$fdeath[tall$dtime>ctime] <- 0; 
	     tall$dtime[tall$dtime>ctime] <- ctime[tall$dtime>ctime] 
	  }
  }# }}}

### fixing the first time to event
  tall$death <- 0
  tall <- dtransform(tall,death=fdeath,time>dtime)
  tall <- dtransform(tall,status=0,time>dtime)
  tall <- dtransform(tall,time=dtime,time>dtime)
  tall <- count.history(tall,lag=FALSE)
  ddrop(tall) <- ~lbnr__id
  dtable(tall,~"Count*")
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
	  rr1 <- exp(beta1*still$Count1)*z2[still$id]
	  rr2 <- exp(beta2*still$Count2)*z1[still$id]
###	  print(dtable(still,~"Count*"))
          tt1 <- timereg::pc.hazard(cumhaz,rr2,entry=still$time)
          tt2 <- timereg::pc.hazard(cumhaz2,rr1,entry=still$time)
	  tt <- tt1
          tt$status <- ifelse(tt1$time<=tt2$time,tt1$status,2*tt2$status)
          tt$time <-   ifelse(tt1$time<=tt2$time,tt1$time,tt2$time)
	  tt$rr2 <- tt2$rr
           ###
	  tt <- cbind(tt,dkeep(still,~id+dtime+death+fdeath),row.names=NULL)
	  tt <- dtransform(tt,death=fdeath,time>dtime)
	  tt <- dtransform(tt,status=0,time>dtime)
	  tt <- dtransform(tt,time=dtime,time>dtime)
	  tt$Count1 <- still$Count1+(tt$status==1)
	  tt$Count2 <- still$Count2+(tt$status==2)
	  nt <- nrow(tt)
	  tall <- rbind(tall,tt[1:nn,],row.names=NULL)
  }
  dsort(tall) <- ~id+entry+time
  tall$start <- tall$entry
  tall$stop  <- tall$time

  attr(tall,"death.cumhaz") <- cumhazd
  attr(tall,"cumhaz") <- cumhaz
  attr(tall,"cumhaz2") <- cumhaz2
  attr(tall,"z") <- z

  return(tall)
  }# }}}


##' Counts the number of previous events of two types for recurrent events processes
##'
##' Counts the number of previous events of two types for recurrent events processes
##'
##' @param data data-frame
##' @param status name of status 
##' @param id  id 
##' @param types types of the events (code) related to status
##' @param names.count name of Counts, for example Count1 Count2 when types=c(1,2)
##' @param lag if true counts previously observed, and if lag=FALSE counts up to know
##' @author Thomas Scheike
##' @examples
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
##' ######################################################################
##' ### simulating simple model that mimicks data 
##' ### now with two event types and second type has same rate as death rate
##' ######################################################################
##'
##' rr <- simRecurrentII(1000,base1,dr,death.cumhaz=base4)
##' rr <-  count.history(rr)
##' dtable(rr,~"Count*"+status,level=1)
##'
##' @export
count.history <- function(data,status="status",id="id",types=1:2,names.count="Count",lag=TRUE)
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
if (lag==TRUE)
data[,paste(names.count,i,sep="")] <- 
   cumsumidstratasum((stat==i),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
   else 
data[,paste(names.count,i,sep="")] <- 
   cumsumidstratasum((stat==i),rep(0,nrow(data)),1,clusters,max.clust+1)$sum 
}

return(data)
}# }}}


##' Estimation of probability of more that k events for recurrent events process
##'
##' Estimation of probability of more that k events for recurrent events process
##' where there is terminal event, based on this also estimate of variance of recurrent events. The estimator is based on cumulative incidence of exceeding "k" events.
##' In contrast the probability of exceeding k events can also be computed as a 
##' counting process integral, and this is implemented in prob.exceedRecurrent
##'
##' @param data data-frame
##' @param type type of evnent (code) related to status
##' @param status name of status 
##' @param death  name of death indicator 
##' @param start start stop call of Hist() of prodlim 
##' @param stop start stop call of Hist() of prodlim 
##' @param id  id 
##' @param times time at which to get probabilites P(N1(t) >= n)
##' @param exceed n's for which which to compute probabilites P(N1(t) >= n)
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @references 
##'             Scheike, Eriksson, Tribler (2018) 
##'             The mean, variance and correlation for bivariate recurrent events
##'             with a terminal event,  work in progress
##'
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
##' cor.mat <- corM <- rbind(c(1.0, 0.6, 0.9), c(0.6, 1.0, 0.5), c(0.9, 0.5, 1.0))
##' rr <- simRecurrent(1000,base1,cumhaz2=base4,death.cumhaz=dr)
##' rr <-  count.history(rr)
##' dtable(rr,~death+status)
##' 
##' oo <- prob.exceedRecurrent(rr,1)
##' bplot(oo)
##' 
##' par(mfrow=c(1,2))
##' with(oo,plot(time,mu,col=2,type="l"))
##' ###
##' with(oo,plot(time,varN,type="l"))
##' 
##' 
##' ### Bivariate probability of exceeding 
##' oo <- prob.exceedBiRecurrent(rr,1,2,exceed1=c(1,5,10),exceed2=c(1,2,3))
##' with(oo, matplot(time,pe1e2,type="s"))
##' nc <- ncol(oo$pe1e2)
##' legend("topleft",legend=colnames(oo$pe1e2),lty=1:nc,col=1:nc)
##' 
##' 
##' \donttest{
##' ### do not test to avoid dependence on prodlim 
##' ### now estimation based on cumualative incidence, but do not test to avoid dependence on prodlim 
##' library(prodlim)
##' pp <- prob.exceed.recurrent(rr,1,status="status",death="death",start="entry",stop="time",id="id")
##' with(pp, matplot(times,prob,type="s"))
##' ###
##' with(pp, matlines(times,se.lower,type="s"))
##' with(pp, matlines(times,se.upper,type="s"))
##' }
##' @export
##' @aliases prob.exceedRecurrent prob.exceedBiRecurrent
prob.exceed.recurrent <- function(data,type,status="status",death="death",
 start="start",stop="stop",id="id",times=NULL,exceed=NULL)
{# {{{
### setting up data 
stat <- data[,status]
dd   <- data[,death]
tstop <- data[,stop]
tstart <- data[,start]
clusters <- data[,id]

if (sum(stat==type)==0) stop("none of these events")

if (is.numeric(clusters)) {
      clusters <- fast.approx(unique(clusters), clusters) - 1
      max.clust <- max(clusters)
}
else {
     max.clust <- length(unique(clusters))
     clusters <- as.integer(factor(clusters, labels = seq(max.clust))) - 1
}

 count <- cumsumstrata((stat==type),clusters,max.clust+1)
### count  <- cumsumidstratasum((stat==type),rep(0,nrow(data)),1,clusters,max.clust+1)$lagsum 
 mc <- max(count)+1
 idcount <- clusters*mc + count
 idcount <- cumsumstrata(rep(1,length(idcount)),idcount,mc*(max.clust+1))

if (is.null(times)) times <- sort(unique(tstop[stat==1]))
if (is.null(exceed)) exceed <- sort(unique(count))

###form <- as.formula(paste("Event(entry=",start,",",stop,",statN)~+1",sep=""))
form <- as.formula(paste("Hist(entry=",start,",",stop,",statN)~+1",sep=""))

probs.orig <- se.probs <- probs <- matrix(0,length(times),length(exceed))
se.lower <-  matrix(0,length(times),length(exceed))
se.upper <-  matrix(0,length(times),length(exceed))
i <- 1
for (n1 in exceed[-1]) {# {{{
	i <- i+1
	### first time that get to n1
	keep <- (count<n1 ) | (count==n1 & idcount==1)
	### status, censoring, get to n1, or die
        statN <- rep(0,nrow(data))
	statN[count==n1] <- 1
	statN[dd==1] <- 2
	statN <- statN[keep]
	pN1 <-  suppressWarnings(prodlim::prodlim(form,data=data[keep,]))

	if (sum(statN)==0) {
		se.lower[,i] <- se.upper[,i] <- se.probs[,i] <- probs[,i] <- rep(0,length(times)) } else  {
		mps  <- suppressWarnings(summary(pN1,times=times,cause=1)$table$"1")
		probs.orig[,i] <- ps <- mps[,5]
		mm <- which.max(ps)
		probs[,i] <- ps
		probs[is.na(ps),i] <- ps[mm]
		se.probs[,i] <- mps[,6]
		se.probs[is.na(ps),i] <- se.probs[mm,i]
		se.lower[,i] <- mps[,7] 
		se.lower[is.na(ps),i] <- se.lower[mm,i]
		se.upper[,i] <- mps[,8]
		se.upper[is.na(ps),i] <- se.upper[mm,i]
	}
	if (i==2) { probs[,1]    <- 1-probs[,2]; 
                    se.probs[,1] <- se.probs[,2]; 
                    se.lower[,1] <- 1-se.lower[,2]; 
                    se.upper[,1] <- 1-se.upper[,2]; 
	}
}# }}}

dp <- -t(apply(cbind(probs[,-1],0),1,diff))
meanN <- apply(probs[,-1,drop=FALSE],1,sum)
meanN2 <- apply(t(exceed[-1]^2 * t(dp)),1,sum)
 
colnames(probs) <- c(paste("N=",exceed[1],sep=""),paste("exceed>=",exceed[-1],sep=""))
colnames(se.probs) <- c(paste("N=",exceed[1],sep=""),paste("exceed>=",exceed[-1],sep=""))

return(list(time=times,times=times,prob=probs,se.prob=se.probs,meanN=meanN,probs.orig=probs.orig[,-1],
	    se.lower=se.lower,se.upper=se.upper,meanN2=meanN2,varN=meanN2-meanN^2,exceed=exceed[-1]))
}# }}}

##' @export
prob.exceedRecurrent <- function(data,type1,km=TRUE,status="status",death="death",
                start="start",stop="stop",id="id",names.count="Count",...)
{# {{{

formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~cluster(",id,")",sep=""))
###
form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~strata(",names.count,type1,")+cluster(",id,")",sep=""))

dr <- phreg(formdr,data=data)
base1   <- phreg(form1,data=data)
base1.2 <- phreg(form1C,data=data)

cc <- base1$cox.prep
risk <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###### risk stratified after count 1
cc <- base1.2$cox.prep
risk1 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
pstrata <- risk1/risk
pstrata[risk1==0] <- 0

### marginal int_0^t G(s) P(N1(t-)==k|D>t) \lambda_{1,N1=k}(s) ds 
### strata og count skal passe sammen
  # {{{
  strat <- dr$strata[dr$jumps]
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- base1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  mu <- c(cumsumstrata(St*S0i,xx$strata,xx$nstrata))
###
  x <- base1.2
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  xstrata <- xx$strata
  vals1 <- sort(unique(data[,paste("Count",type1,sep="")]))
  valjumps <- vals1[xx$strata+1]
  fk <- (valjumps+1)^2-valjumps^2
  EN2 <- c(cumsumstrata(fk*St*pstrata*S0i,rep(0,lss),1))
###  pstrata <- matdoubleindex(nrisk,1:nrow(nrisk),xstrata+1)
  pcumhaz <- cbind(xx$time,
              cumsumstrata(pstrata*St*S0i,xx$strata,xx$nstrata))
# }}}
  EN2     <- EN2[xx$jumps+1]
  cumhaz <- pcumhaz[xx$jumps+1,]
  mu     <- mu[xx$jumps+1]

  exceed.name <- paste("Exceed>=",vals1+1,sep="")

  out=list(cumhaz=cumhaz,time=cumhaz[,1],
	varN=EN2-mu^2,mu=mu,
	nstrata=base1.2$nstrata,strata=base1.2$strata[xx$jumps+1],
	jumps=1:nrow(cumhaz),
	strat.cox.name=base1.2$strata.name,
	strat.cox.level=base1.2$strata.level,exceed=vals1+1,
        strata.name=exceed.name,strata.level=exceed.name)

### use recurrentMarginal estimator til dette via strata i base1 
### strata og count skal passe sammen
### se beregning via recurrent marginal function

###  base1$cox.prep$strata <- base1.2$cox.prep$strata
###  base1$cox.prep$nstrata <- base1.2$cox.prep$nstrata
###  base1$nstrata <- base1.2$cox.prep$nstrata
###  base1$strata <- base1.2$strata
###  base1$strata.name <- base1.2$strata.name
###  base1$strata.level <- base1.2$strata.level

###  mm <- recurrentMarginal(base1,dr,km=km,...)
###  out=c(mm,list(varN=EN2-mu^2))

  return(out)
}# }}}

##' @export
prob.exceedBiRecurrent <- function(data,type1,type2,km=TRUE,status="status",death="death",
      start="start",stop="stop",id="id",names.count="Count",exceed1=NULL,exceed2=NULL)
{# {{{

formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))

form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~cluster(",id,")",sep=""))
form2 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~cluster(",id,")",sep=""))
###
###form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))
###form2C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~
###  strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))

form2Ccc <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~
  ",names.count,type1,"+",names.count,type2,"+","
  strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))
form1Ccc <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~
  ",names.count,type1,"+",names.count,type2,"+","
  strata(",names.count,type1,",",names.count,type2,")+cluster(",id,")",sep=""))


### stratified and with counts in covariate matrix 
bb2.12 <- phreg(form2Ccc,data=data,no.opt=TRUE)
bb1.12 <- phreg(form1Ccc,data=data,no.opt=TRUE)

dr <- phreg(formdr,data=data)
base1   <- phreg(form1,data=data)
base2   <- phreg(form2,data=data)

cc <- base1$cox.prep
risk1 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###### risk stratified after count 1 og count2
cc <- bb1.12$cox.prep
risk1.12 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
pstrata1 <- risk1.12/risk1
pstrata1[1] <- 0

cc <- base2$cox.prep
risk2 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###### risk stratified after count 1 og count2
cc <- bb2.12$cox.prep
risk2.12 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
pstrata2 <- risk2.12/risk2
pstrata2[1] <- 0

### marginal int_0^t G(s) P(N1(t-)==k|D>t) \lambda_{1,N1=k}(s) ds 
### strata og count skal passe sammen

  # {{{
  strat <- dr$strata[dr$jumps]
  Gt <- exp(-dr$cumhaz[,2])
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  if (!km) { 
     cumhazD <- c(cumsumstratasum(S0i,xx$strata,xx$nstrata)$lagsum)
     St      <- exp(-cumhazD)
  } else St <- c(exp(cumsumstratasum(log(1-S0i),xx$strata,xx$nstrata)$lagsum))
  x <- base1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  mu <- c(cumsumstrata(St*S0i,rep(0,lss),1))
###

  x <- bb1.12
  xx1 <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx1$jumps+1] <-  1/x$S0
  dcumhaz1 <- cbind(xx1$time,pstrata1*St*S0i)
###              cumsumstrata(pstrata1*St*S0i,xx1$strata,xx1$nstrata))
  x <- bb2.12
  xx2 <- x$cox.prep
  lss <- length(xx2$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx2$jumps+1] <-  1/x$S0
  dcumhaz2 <- cbind(xx2$time,pstrata2*St*S0i)
###   cumsumstrata(pstrata2*St*S0i,xx2$strata,xx2$nstrata))

  n1 <- length(xx1$jumps)
  n2 <- length(xx2$jumps)
  ojumps <- order(c(xx1$jumps,xx2$jumps))
  jumps <- sort(c(xx1$jumps,xx2$jumps))
  times <- xx1$time[jumps+1]

  dcumhaz1 <-  dcumhaz1[jumps+1,] 
  dcumhaz2 <-  dcumhaz2[jumps+1,] 
  x1 <- xx1$X[jumps+1,]
  x2 <- xx2$X[jumps+1,]

###  dcumhaz1 <- dcumhaz1[xx1$jumps+1,]
###  dcumhaz2 <- dcumhaz2[xx2$jumps+1,]
###  x1 <- xx1$X[xx1$jumps+1,]
###  x2 <- xx2$X[xx2$jumps+1,]
# }}}

  if (is.null(exceed1)) exceed1 <- 1:max(x1[,1])
  if (is.null(exceed2)) exceed2 <- 1:max(x1[,2])

  pe1e2 <- matrix(0,n1+n2,length(exceed1)*length(exceed2))
  m <- 0; nn <- c()
  for (i in exceed1) 
  for (j in exceed2)  {
	  m <- m+1
	  strat1 <- (x1[,2]>=j)*(x1[,1]==(i-1))
	  strat2 <- (x2[,1]>=i)*(x2[,2]==(j-1))
	  pe1e2[,m] <- cumsum(strat1*dcumhaz1[,2]) + cumsum(strat2*dcumhaz2[,2])
	  nn <- c(nn,paste("N_1(t)>=",i,",N_2(t)>=",j,sep="")) 
  }

  colnames(pe1e2) <- nn

  out=list(time=times,pe1e2=pe1e2,
###        pcumhaz1=pcumhaz1,pcumhaz2=pcumhaz2,
	x1=x1,x2=x2)

###	jumps1=1:nrow(pcumhaz1),
###	nstrata=bb1.12$nstrata,
###	strata1=bb1.12$strata[xx1$jumps+1],
###        strata.name1=bb1.12$strata.name,
###	strata.level1=bb1.12$strata.level)

  return(out)
}# }}}

##' Estimation of covariance for bivariate recurrent events with terminal event
##'
##' Estimation of probability of more that k events for recurrent events process
##' where there is terminal event 
##'
##' @param data data-frame
##' @param type1 type of first event (code) related to status
##' @param type2 type of second event (code) related to status
##' @param status name of status 
##' @param death  name of death indicator 
##' @param start start stop call of Hist() of prodlim 
##' @param stop start stop call of Hist() of prodlim 
##' @param id  id 
##' @param names.count name of count for number of previous event of different types, here generated by count.history()
##' @author Thomas Scheike
##' @references 
##'             Scheike, Eriksson, Tribler (2018) 
##'             The mean, variance and correlation for bivariate recurrent events
##'             with a terminal event,  work in progress
##'
##' @examples
##'
##' ########################################
##' ## getting some data to work on 
##' ########################################
##' data(base1cumhaz)
##' data(base4cumhaz)
##' data(drcumhaz)
##' dr <- drcumhaz
##' base1 <- base1cumhaz
##' base4 <- base4cumhaz
##' rr <- simRecurrent(1000,base1,cumhaz2=base4,death.cumhaz=dr)
##' rr <- count.history(rr)
##' rr$strata <- 1
##' dtable(rr,~death+status)
##' 
##' covrp <- covarianceRecurrent(rr,1,2,status="status",death="death",
##'                         start="entry",stop="time",id="id",names.count="Count")
##' par(mfrow=c(1,3)) 
##' plot(covrp)
##' 
##' ### with strata, each strata in matrix column, provides basis for fast Bootstrap
##' covrpS <- covarianceRecurrentS(rr,1,2,status="status",death="death",
##'         start="entry",stop="time",strata="strata",id="id",names.count="Count")
##' 
##' @aliases plot.covariace.recurrent  covarianceRecurrentS Bootcovariancerecurrence BootcovariancerecurrenceS 
##' @export
covarianceRecurrent <- function(data,type1,type2,status="status",death="death",
                   	 start="start",stop="stop",id="id",names.count="Count")
{# {{{

formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~ cluster(",id,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~cluster(",id,")",sep=""))
form2 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~cluster(",id,")",sep=""))
###
form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~strata(",names.count,type2,")+cluster(",id,")",sep=""))
form2C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~strata(",names.count,type1,")+cluster(",id,")",sep=""))

dr <- phreg(formdr,data=data)
###basehazplot.phreg(dr)
###
base1   <- phreg(form1,data=data)
base1.2 <- phreg(form1C,data=data)
###
base2   <- phreg(form2,data=data)
base2.1 <- phreg(form2C,data=data)

###marginal.mean1 <- recurrentMarginal(base1,dr)
###marginal.mean2 <- recurrentMarginal(base2,dr)

marginal.mean1 <- recmarg(base1,dr)
marginal.mean2 <- recmarg(base2,dr)

cc <- base2$cox.prep
risk <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###### risk stratified after count 1
cc <- base2.1$cox.prep
risk1 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
ssshed1 <- risk1/risk
ssshed1[1,1] <- 1
sshed1   <- list(cumhaz=cbind(cc$time,ssshed1),
	      strata=cc$strata,nstrata=cc$nstrata,
	      jumps=1:length(cc$time),
	      strata.name=paste("prob",type1,sep=""),
	      strata.level=base2.1$strata.level)
###
riskstrata <- .Call("riskstrataR",cc$sign,cc$strata,cc$nstrata)$risk
nrisk <- apply(riskstrata,2,revcumsumstrata,rep(0,nrow(riskstrata)),1)
ntot <- apply(nrisk,1,sum)
vals1 <- sort(unique(data[,paste("Count",type1,sep="")]))
mean1risk <- apply(t(nrisk)*vals1,2,sum)/ntot
mean1risk[1] <- 0

cc <- base1$cox.prep
S0 <- rep(0,length(cc$strata))
risk <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
###
cc <- base1.2$cox.prep
S0 <- rep(0,length(cc$strata))
risk2 <- revcumsumstrata(cc$sign,cc$strata,cc$nstrata)
ssshed2 <- risk2/risk
ssshed2[1,1] <- 1
###
sshed2  <- list(cumhaz=cbind(cc$time,ssshed2),
	      strata=cc$strata,nstrata=cc$nstrata,
              jumps=1:length(cc$time),
	      strata.name=paste("prob",type2,sep=""),
	      strata.level=base1.2$strata.level)
###
riskstrata <- .Call("riskstrataR",cc$sign,cc$strata,cc$nstrata)$risk
nrisk <- apply(riskstrata,2,revcumsumstrata,rep(0,nrow(riskstrata)),1)
ntot <- apply(nrisk,1,sum)
vals2 <- sort(unique(data[,paste("Count",type2,sep="")]))
mean2risk <- apply(t(nrisk)*vals2,2,sum)/ntot
mean2risk[1] <- 0

mu1 <-timereg::Cpred(rbind(c(0,0),marginal.mean1$cumhaz),cc$time)[,2]
mu2 <-timereg::Cpred(rbind(c(0,0),marginal.mean2$cumhaz),cc$time)[,2]

out <- list(based=dr,base1=base1,base2=base2,
	    base1.2=base1.2,base2.1=base2.1,
	    marginal.mean1=marginal.mean1,marginal.mean2=marginal.mean2,
	    prob1=sshed1,prob2=sshed2,
	    mean1risk=mean1risk,mean2risk=mean2risk)

### marginal sum_k int_0^t G(s) k P(N1(t-)==k|D>t) \lambda_{2,N1=k}(s) ds 
### strata og count skal passe sammen
  # {{{
  strat <- dr$strata[dr$jumps]
  Gt <- exp(-dr$cumhaz[,2])
  ###
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
###
  x <- base1.2
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  xstrata <- xx$strata
  cumhazDR <- cbind(xx$time,cumsumstrata(vals2[xstrata+1]*St*ssshed2*S0i,rep(0,lss),1))
  x <- base1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazIDR <- cbind(xx$time,cumsumstrata(St*mean2risk*S0i,rep(0,lss),1))
  mu1.i <- cumhazIDR[,2]
  mu1.2 <- cumhazDR[,2]
###
  x <- base2.1
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazDR <- cbind(xx$time,cumsumstrata(vals1[xx$strata+1]*St*ssshed1*S0i,rep(0,lss),1))
###
  x <- base2
  xx <- x$cox.prep
  lss <- length(xx$strata)
  S0i2 <- S0i <- rep(0,lss)
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazIDR <- cbind(xx$time,cumsumstrata(St*mean1risk*S0i,rep(0,lss),1))
  mu2.i <- cumhazIDR[,2]
  mu2.1 <- cumhazDR[,2]
# }}}

  out=c(out,list(EN1N2= mu1.2+mu2.1,mu2.1=mu2.1,mu1.2=mu1.2, 
		 mu2.i=mu2.i,mu1.i=mu1.i,
		 EIN1N2=mu2.i+mu1.i,EN1EN2=mu1*mu2,time=cc$time))

  class(out) <- "covariance.recurrent"

  return(out)
}# }}}

##' @export
plot.covariance.recurrent <- function(x,main="Covariance",...) 
{# {{{

legend <- NULL # to avoid R-check 

with(x, { plot(time,mu1.2,type="l",ylim=range(c(mu1.2,mu1.i)),...) 
          lines(time,mu1.i,col=2) })
legend("topleft",c(expression(integral(N[2](s)*dN[1](s),0,t)),"independence"),lty=1,col=1:2) 
title(main=main)
###
with(x, { plot(time,mu2.1,type="l",ylim=range(c(mu2.1,mu2.i)),...) 
          lines(time,mu2.i,col=2) 
	      })
legend("topleft",c(expression(integral(N[1](s)*dN[2](s),0,t)),"independence"),lty=1,col=1:2) 
title(main=main)
###
with(x,  { plot(time,EN1N2,type="l",lwd=2,ylim=range(c(EN1N2,EN1EN2,EIN1N2)),...) 
           lines(time,EN1EN2,col=2,lwd=2) 
           lines(time,EIN1N2,col=3,lwd=2) })
legend("topleft",c("E(N1N2)", "E(N1) E(N2) ", "E_I(N1 N2)-independence"),lty=1,col=1:3)
title(main=main)
} # }}}

meanRisk <- function(base1,base1.2)
{# {{{
cc <- base1.2$cox.prep
S0 <- rep(0,length(cc$strata))
mid <- max(cc$id)+1
risk2 <- revcumsumidstratasum(cc$sign,cc$id,mid,cc$strata,cc$nstrata)$sumidstrata

means <- .Call("meanriskR",cc$sign,cc$id,mid,cc$strata,cc$nstrata)
mean2risk <- means$meanrisk
mean2risk[is.na(mean2risk)] <- 0
risk <- means$risk
ssshed2 <- risk2/risk
ssshed2[is.na(ssshed2)] <- 0
vals2 <- unique(cc$id)

means2  <- list(cumhaz=cbind(cc$time,mean2risk),
	      strata=cc$strata,nstrata=cc$nstrata,
              jumps=1:length(cc$time),
	      strata.name="meansrisk",
	      strata.level=base1.2$strata.level)
sshed2  <- list(cumhaz=cbind(cc$time,ssshed2),
	      strata=cc$id,nstrata=mid,
              jumps=1:length(cc$time),
	      strata.name="prob",
	      strata.level=paste(vals2),real.strata=cc$strata)

return(list(meanrisk=means2,vals=vals2,probs=sshed2,jumps=cc$jumps+1))
}# }}}

intN2dN1 <- function(dr,base1,base1.2,pm)
{# {{{
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
###
  x <- base1.2
  xx <- x$cox.prep
  xstrata <- xx$id
  jumps <- xx$jumps+1
  mid <- max(xx$id)+1
  ## risk after both id~Count and strata 
  risk2 <- revcumsumidstratasum(xx$sign,xx$id,mid,xx$strata,xx$nstrata)$sumidstrata
  S0i <-  1/risk2[jumps]
  vals <- xx$id[jumps]
  St <- St[jumps]
  probs <- pm$probs$cumhaz[jumps,2]
  cumhazDR <- cbind(xx$time[jumps],cumsumstrata(St*vals*probs*S0i,xx$strata[jumps],xx$nstrata))
  x <- base1
  xx <- x$cox.prep
  S0i <-  c(1/x$S0)
  meanrisk <- pm$meanrisk$cumhaz[jumps,2]
  cumhazIDR <- cbind(xx$time[jumps],cumsumstrata(St*meanrisk*S0i,xx$strata[jumps],xx$nstrata))
  mu1.i <- cumhazIDR[,2]
  mu1.2 <- cumhazDR[,2]
  return(list(cumhaz=cumhazDR,cumhazI=cumhazIDR,mu1.i=mu1.i,mu1.2=mu1.2,
	      time=cumhazDR[,1],
	      strata=xx$strata[jumps],nstrata=xx$nstrata,jumps=1:length(mu1.i),
	      strata.name="intN2dN1",strata.level=x$strata.level))
}# }}}

recmarg2 <- function(recurrent,death,...)
{# {{{
  xr <- recurrent
  dr <- death 

  ### marginal expected events  int_0^t G(s) \lambda_r(s) ds 
  # {{{
  x <- dr
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <- 1/x$S0
  cumhazD <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  St      <- exp(-cumhazD[,2])
  ###
  x <- xr
  xx <- x$cox.prep
###  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i <-  1/x$S0
  jumps <- xx$jumps+1
  cumhazDR <- cbind(xx$time[jumps],cumsumstrata(St[jumps]*S0i,xx$strata[jumps],xx$nstrata))
  mu <- cumhazDR[,2]
# }}}

 varrs <- data.frame(mu=mu,time=cumhazDR[,1],strata=xr$strata[jumps],St=St[jumps])
 out <- list(mu=varrs$mu,times=varrs$time,St=varrs$St,cumhaz=cumhazDR,
     strata=varrs$strata,nstrata=xr$nstrata,jumps=1:nrow(varrs),
     strata.name=xr$strata.name)
 return(out)
}# }}}

##' @export
covarianceRecurrentS <- function(data,type1,type2,times=NULL,status="status",death="death",
                   	 start="start",stop="stop",id="id",names.count="Count",
			 strata="NULL",plot=0,output="matrix")
{# {{{


if (is.null(times)) times <- seq(0,max(data[,stop]),100)

### passing strata as id to be able to use for stratified calculations 
if (is.null(strata)) stop("must give strata, for example one strata\n"); 
## uses Counts1 as cluster to pass to risk set calculations 


formdr <- as.formula(paste("Surv(",start,",",stop,",",death,")~strata(",strata,")",sep=""))
form1 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~strata(",strata,")",sep=""))
form2 <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~strata(",strata,")",sep=""))
###
form1C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type1,")~cluster(",names.count,type2,")+strata(",strata,")",sep=""))
form2C <- as.formula(paste("Surv(",start,",",stop,",",status,"==",type2,")~cluster(",names.count,type1,")+strata(",strata,")",sep=""))


dr <- phreg(formdr,data=data)
###basehazplot.phreg(dr)
###
base1   <- phreg(form1,data=data)
base1.2 <- phreg(form1C,data=data)
###
base2   <- phreg(form2,data=data)
base2.1 <- phreg(form2C,data=data)

rm1 <- recmarg2(base1,dr)
rm2 <- recmarg2(base2,dr)
if (plot==1) {
basehazplot.phreg(rm1)
basehazplot.phreg(rm2)
}


pm1 <- meanRisk(base2,base2.1)
pm2 <- meanRisk(base1,base1.2)
if (plot==1) {
basehazplot.phreg(pm1$meanrisk)
basehazplot.phreg(pm2$meanrisk)
}


### marginal sum_k int_0^t G(s) k P(N1(t-)==k|D>t) \lambda_{2,N1=k}(s) ds 
###          sum_k int_0^t G(s) E(N1(t-)==k|D>t)   \lambda_{2}(s) ds 
iN2dN1 <-  intN2dN1(dr,base1,base1.2,pm2)
iN1dN2 <-  intN2dN1(dr,base2,base2.1,pm1)

###print("hej")
if (plot==1) {
par(mfrow=c(2,2))
basehazplot.phreg(iN2dN1)
plot(iN2dN1$time,iN2dN1$cumhazI[,2],type="l")
basehazplot.phreg(iN1dN2)
}

#### writing output in matrix form for each strata for the times
mu1g <- matrix(0,length(times),rm1$nstrata)
mu2g <- matrix(0,length(times),rm1$nstrata)
mu1.2 <- matrix(0,length(times),rm1$nstrata)
mu2.1 <- matrix(0,length(times),rm1$nstrata)
mu1.i <- matrix(0,length(times),rm1$nstrata)
mu2.i <- matrix(0,length(times),rm1$nstrata)
mu1 <- matrix(0,length(times),rm1$nstrata)
mu2 <- matrix(0,length(times),rm1$nstrata)

if (output=="matrix") {
all <- c()
i <- 1
### going through strata 
for (i in 1:rm1$nstrata) {
	j <- i-1
	mu1[,i]   <- timereg::Cpred(rm1$cumhaz[rm1$strata==j,],times)[,2]
	mu2[,i]   <- timereg::Cpred(rm2$cumhaz[rm2$strata==j,],times)[,2]
	mu1.2[,i] <- timereg::Cpred( iN2dN1$cumhaz[rm1$strata==j,],times)[,2]
	mu1.i[,i] <- timereg::Cpred(iN2dN1$cumhazI[rm1$strata==j,],times)[,2]
	mu2.1[,i] <- timereg::Cpred( iN1dN2$cumhaz[rm2$strata==j,],times)[,2]
	mu2.i[,i] <- timereg::Cpred(iN1dN2$cumhazI[rm2$strata==j,],times)[,2]
}
mu1mu2 <- mu1*mu2

out=c(list(EN1N2=mu1.2+mu2.1, EIN1N2=mu1.i+mu2.i, EN1EN2=mu1mu2,
	   mu1.2=mu1.2, mu1.i=mu1.i, mu2.1=mu2.1, mu2.i=mu2.i,
	   mu1=mu1,mu2=mu2, nstrata=rm1$nstrata, time=times))

} else out <- list(iN2dN1=iN2dN1,iN1dN2=iN1dN2,rm1=rm1,rm2=rm2,pm1=pm1,pm2=pm2)


  return(out)
}# }}} 

##' @export
BootcovariancerecurrenceS <- function(data,type1,type2,status="status",death="death",
	 start="start",stop="stop",id="id",names.count="Count",times=NULL,K=100)
{# {{{


  if (is.null(times)) times <- seq(0,max(data[,stop]),length=100)
  mu1.2 <- matrix(0,length(times),K)
  mu2.1 <- matrix(0,length(times),K)
  mu1.i <- matrix(0,length(times),K)
  mu2.i <- matrix(0,length(times),K)
  mupi <- matrix(0,length(times),K)
  mupg <- matrix(0,length(times),K)
  mu1mu2 <- matrix(0,length(times),K)
  n <- length(unique(data[,id]))

  formid <- as.formula(paste("~",id))
  rrb <- blocksample(data, size = n*K, formid)
  rrb$strata <- floor((rrb$id-0.01)/n)
  rrb$jump <- (rrb[,status] %in% c(type1,type2)) | (rrb[,death]==1)
  rrb <- tie.breaker(rrb,status="jump",start=start,stop=stop,id=id)

  mm <- covarianceRecurrentS(rrb,type1,type2,status=status,death=death,
 	               start=start,stop=stop,id=id,names.count=names.count,
		       strata="strata",times=times)

  mm <- c(mm,list(se.mui=apply(mm$EIN1N2,1,sd),se.mug=apply(mm$EN1N2,1,sd)))
  return(mm)

}# }}}

##' @export
Bootcovariancerecurrence <- function(data,type1,type2,status="status",death="death",
	 start="start",stop="stop",id="id",names.count="Count",times=NULL,K=100)
{# {{{

  strata <- NULL # to avoid R-check 

  if (is.null(times)) times <- seq(0,max(data[,stop]),length=100)
  mu1.2 <- matrix(0,length(times),K)
  mu2.1 <- matrix(0,length(times),K)
  mu1.i <- matrix(0,length(times),K)
  mu2.i <- matrix(0,length(times),K)
  mupi <- matrix(0,length(times),K)
  mupg <- matrix(0,length(times),K)
  mu1mu2 <- matrix(0,length(times),K)
  n <- length(unique(data[,id]))

  formid <- as.formula(paste("~",id))
  rrb <- blocksample(data, size = n*K, formid)
  rrb$strata <- floor((rrb$id-0.01)/n)
## rrb$jump <- (rrb[,status]!=0) | (rrb[,death]==1)
  rrb$jump <- (rrb[,status] %in% c(type1,type2)) | (rrb[,death]==1)
  rrb <- tie.breaker(rrb,status="jump",start=start,stop=stop,id=id)

  for (i in 1:K)
  {
     rrbs <- subset(rrb,strata==i-1)
     errb <- covarianceRecurrent(rrbs,type1,type2,status=status,death=death,
 	                          start=start,stop=stop,id=id,names.count=names.count)
     all <- timereg::Cpred(cbind(errb$time,errb$EIN1N2,errb$EN1N2,errb$EN1EN2,
			errb$mu1.2,errb$mu2.1,errb$mu1.i,errb$mu2.i),times)
     mupi[,i]   <- all[,2]
     mupg[,i]   <- all[,3]
     mu1mu2[,i] <- all[,4]
     mu1.2[,i]  <- all[,5]; 
     mu2.1[,i]  <- all[,6]; 
     mu1.i[,i]  <- all[,7]; 
     mu2.i[,i]  <- all[,8]
  }

return(list(mupi=mupi,mupg=mupg,mu1mu2=mu1mu2,time=times,
	    mu1.2=mu1.2,mu1.i=mu1.i,mu2.1=mu2.1,mu2.i=mu1.i,
	    mup=apply(mupi,1,mean),mug=apply(mupg,1,mean),
	    dmupg=apply(mupg-mupi,1,mean),mmu1mu2=apply(mu1mu2,1,mean), 
	    se.mui=apply(mupi,1,sd),se.mug=apply(mupg,1,sd)
###	    se.dmug=apply(mupg-mu1mu2,1,sd),se.dmui=apply(mupi-mu1mu2,1,sd),
###	    se.difmugmup=apply(mupg-mupi,1,sd), se.mu1mu2=apply(mu1mu2,1,sd)
	    ))
}# }}}

iidCovarianceRecurrent <-  function (rec1,death,xrS,xr,means)
{# {{{
    ### makes iid decompition for covariance under independence between events
    axr <- rec1
    adr <- death
    St <- exp(-adr$cum[, 2])
    timesr <- axr$cum[, 1]
    timesd <- adr$cum[, 1]
    times <- c(timesr[-1], timesd[-1])
    or <- order(times)
    times <- times[or]
    meano <- cbind(means$time,means$mean2risk)
    imeano <- sindex.prodlim(means$time, times, strict = FALSE)
    meano <- meano[imeano,2]
    keepr <- order(or)[1:length(timesr[-1])]
    rid <- sindex.prodlim(timesd, times, strict = FALSE)
    rir <- sindex.prodlim(timesr, times, strict = FALSE)
    Stt <- St[rid]
    ariid <- axr$cum[rir, 2]
    mu <- cumsum(meano * Stt * diff(c(0, ariid)))
    muS <- cumsum( Stt * diff(c(0, ariid)))
    nc <- length(axr$B.iid)
    muiid <- matrix(0, length(times), nc)

   cc <- xrS$cox.prep
   rrs <- .Call("riskstrataR",cc$sign*cc$strata,cc$id,max(cc$id)+1)$risk
   rr <- .Call("riskstrataR",cc$sign,cc$id,max(cc$id)+1)$risk
   ntot <- revcumsumstrata(cc$sign,rep(0,nrow(rr)),1)
   rr  <- apply(rr,2,revcumsumstrata,rep(0,nrow(rr)),1)
   rrs <- apply(rrs,2,revcumsumstrata,rep(0,nrow(rr)),1)

   xrid <- sindex.prodlim(cc$time, times, strict = FALSE)
   rr <- rr[xrid,]
   rrs <- rrs[xrid,]
   ntot <- ntot[xrid]
   rrs <- apply(Stt* diff(c(0,ariid))*rrs,2,cumsum)
###   dmu <- diff(c(0,mu))
   rrcum <- apply(Stt*meano*diff(c(0,ariid))*rr,2,cumsum)
   miid <- (rrs-rrcum)/ntot

    for (i in 1:nc) {
        mriid <- axr$B.iid[[i]]
        mdiid <- adr$B.iid[[i]]
        mriid <- mriid[rir]
        mdiid <- mdiid[rid]
        dmridd <- diff(c(0, mriid))
        dmdidd <- diff(c(0, mdiid))
        muiid[, i] <- cumsum(Stt * meano* dmridd) - mu * cumsum(dmdidd) +
            cumsum(mu * dmdidd) +  miid[,i]
    }
    var1 <- apply(muiid^2, 1, sum)
    se.mu <- var1[keepr]^0.5
    mu = mu[keepr]
    timeso <- times
    times <- times[keepr]

    out = list(iidtimes=timeso,muiid=muiid,times=times, 
        mu = mu, var.mu = var1[keepr], se.mu = se.mu, St = St, Stt = Stt[keepr], 
	cumhaz=cbind(times,mu),se.cumhaz=cbind(times,se.mu), 
	nstrata=1,strata=rep(0,length(mu)),jumps=1:length(mu)) 

}# }}} 

