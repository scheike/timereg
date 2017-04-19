#' Simulation of Piecewise constant hazard model (Cox).
#' 
#' Simulates data from piecwise constant baseline hazard that can also be of
#' Cox type. Censor data at highest value of the break points.
#' 
#' 
#' @param cumhazard cumulative hazard, or piece-constant rates for periods
#' defined by first column of input.
#' @param rr number of simulations or vector of relative risk for simuations.
#' @param cum.hazard specifies wheter input is cumulative hazard or rates.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' 
#' rates <-  c(0,0.01,0.052,0.01,0.04)
#' breaks <- c(0,10,   20,  30,   40)
#' haz <- cbind(breaks,rates)
#' n <- 1000
#' X <- rbinom(n,1,0.5)
#' beta <- 0.2
#' rrcox <- exp(X * beta)
#' cumhaz <- cumsum(c(0,diff(breaks)*rates[-1]))
#' cumhaz <- cbind(breaks,cumhaz)
#' 
#' pctime <- pc.hazard(haz,1000,cum.hazard=FALSE)
#' 
#' par(mfrow=c(1,2))
#' ss <- aalen(Surv(time,status)~+1,data=pctime,robust=0)
#' plot(ss)
#' lines(cumhaz,col=2,lwd=2)
#' 
#' pctimecox <- pc.hazard(cumhaz,rrcox)
#' pctime <- cbind(pctime,X)
#' ssx <- cox.aalen(Surv(time,status)~+prop(X),data=pctimecox,robust=0)
#' plot(ssx)
#' lines(cumhaz,col=2,lwd=2)
#' 
#' ### simulating data with hazard as real data 
#' data(TRACE)
#' 
#' par(mfrow=c(1,2))
#' ss <- cox.aalen(Surv(time,status==9)~+prop(vf),data=TRACE,robust=0)
#' par(mfrow=c(1,2))
#' plot(ss)
#' ###
#' pctime <- pc.hazard(ss$cum,1000)
#' ###
#' sss <- aalen(Surv(time,status)~+1,data=pctime,robust=0)
#' lines(sss$cum,col=2,lwd=2)
#' 
#' pctime <- pc.hazard(ss$cum,rrcox)
#' pctime <- cbind(pctime,X)
#' ###
#' sss <- cox.aalen(Surv(time,status)~+prop(X),data=pctime,robust=0)
#' summary(sss)
#' plot(ss)
#' lines(sss$cum,col=3,lwd=3)
#' 
#' @export
#' @aliases pchazard.sim cause.pchazard.sim 
pc.hazard <- function(cumhazard,rr,cum.hazard=TRUE)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE
  if (length(rr)==1) rr<-rep(1,rr)
  breaks <- cumhazard[,1]
  rates <- cumhazard[,2][-1]
  mm <- tail(breaks,1)
  if (cum.hazard==FALSE) {
        cumh <- cumsum(c(0,diff(breaks)*rates))
        cumhazard <- cbind(breaks,cumh)
  } else cumh <- cumhazard[,2] 
   n <- length(rr)
   ttt <- rexp(n)/rr
###
   ri <- sindex.prodlim(cumh,ttt)
   rr <- (ttt-cumh[ri])/(cumh[ri+1]-cumh[ri])
   rrx <- rr*(breaks[ri+1]-breaks[ri])+breaks[ri]
   status <- rep(1,n)
   status[is.na(rrx)] <- 0
   rrx[is.na(rrx)] <- mm
   dt <- data.frame(time=rrx,status=status)
   attr(dt,"cumhaz") <- cumhazard
   return(dt)
}# }}}

#' @export
pchazard.sim <- function(cumhazard,rr,cens=NULL,rrc=NULL,cens.cum.hazard=TRUE)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE

   ptt <- pc.hazard(cumhazard,rr,cum.hazard=TRUE)

   if (length(rr)==1) n<-rr else n <- length(rr)
   if (is.null(rrc)) {
	   if (length(rr)==1) rrc<-rr else rrc <- length(rr)
   }
   if (is.matrix(cens)) {
	   pct <- pc.hazard(cumhazard,rr,cum.hazard=cens.cum.hazard)
	   pct <- pct$time
   }
   else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
	   }
   }
   dt <- data.frame(time=pmin(ptt$time,pct), status=ifelse(ptt$time<pct,ptt$status,0))

   return(dt)
}# }}}

#' @export
cause.pchazard.sim <- function(cumhaz1,cumhaz2,rr1,rr2,cens=NULL,rrc=NULL,cens.cum.hazard=TRUE)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE

   if (length(rr1)==1) n<-rr1 else n <- length(rr1)
   if (missing(rr2)) rr2 <- n

   ptt1 <- pc.hazard(cumhaz1,rr1,cum.hazard=TRUE)
   ptt2 <- pc.hazard(cumhaz2,rr2,cum.hazard=TRUE)
   ptt <- data.frame(time=pmin(ptt1$time,ptt2$time),status=ifelse(ptt1$time<=ptt2$time,ptt1$status,ptt2$status*2))

   if (!is.null(cens)) {
      if (is.null(rrc)) rrc <- n 
	   if (is.matrix(cens)) {
		   pct <- pc.hazard(cens,rrc,cum.hazard=cens.cum.hazard)
		   pct <- pct$time
	   } else {
		   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
		      chaz <-sum(ptt1$status+ptt2$status)/sum(ptt1$time+ptt2$time)  ## hazard averate T haz 
		      pct<- rexp(n)/chaz 
		   }
	   }
	   ptt <- data.frame(time=pmin(ptt$time,pct), status=ifelse(ptt$time<pct,ptt$status,0))
   }

   return(ptt)
}# }}}


#' Simulation of output from Cox model.
#' 
#' Simulates data that looks like fit from Cox model. Censor data automatically
#' for highest value of the break points.
#' 
#' 
#' @param cox output form coxph or cox.aalen model fitting cox model.
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed
#' covariates).
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' given then takes average rate of in simulated data from cox model.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' data(TRACE)
#' 
#' cox <-  coxph(Surv(time,status==9)~vf+chf+wmi,data=TRACE)
#' sim1 <- sim.cox(cox,1000,data=TRACE)
#' cc <- coxph(Surv(time,status)~vf+chf+wmi,data=sim1)
#' cbind(cox$coef,cc$coef)
#' 
#' cor(sim1[,c("vf","chf","wmi")])
#' cor(TRACE[,c("vf","chf","wmi")])
#' ###library(mets)
#' ###dcor(sim1,~vf+chf+wmi)
#' ###dcor(TRACE,~vf+chf+wmi)
#' 
#' cox <-  cox.aalen(Surv(time, status==9) ~ prop(vf)+prop(chf)+prop(wmi),TRACE,robust=0)
#' sim2 <- sim.cox(cox,1000,data=TRACE)
#' cc <-  cox.aalen(Surv(time, status)~prop(vf)+prop(chf)+prop(wmi),data=sim2,robust=0)
#' ###
#' plot(cox)
#' lines(cc$cum,type="s",col=2)
#' cbind(cox$gamma,cc$gamma)
#' 
#' @export
sim.cox <- function(cox,n,data=NULL,cens=NULL,rrc=NULL)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE

if (class(cox)=="coxph")
{
   base <- basehaz(cox,centered=FALSE)[,c(2,1)] 
   mt <- model.frame(cox)
   Y <- model.extract(mt, "response")
   if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
   if (attr(Y, "type") == "right") {
    time <- Y[, "time"]; 
    status <- Y[, "status"]
   } else stop("Expected right-censored data.");
  Z <- na.omit(model.matrix(cox))
  nn <- colnames(Z)
  rownames(Z) <- NULL
  jtime <- sort(time[status==1])
  cumhazard <- Cpred(base,jtime)
  rr <- exp( Z %*% matrix(coef(cox),ncol=1))
  xid <- sample(1:nrow(Z),n,replace=TRUE)
  rr <- rr[xid]
  Z <- Z[xid,]
  colnames(Z) <- nn
}
if (class(cox)=="cox.aalen")
{
   formula <- attr(cox, "Formula")
###   Z1 <- na.omit(model.matrix(cox,data=data))
   p <- nrow(cox$gamma)
   Z <- na.omit(get_all_vars(formula,data=data))
   nz <- ncol(Z)
   Z <- Z[,seq(nz-p+1,nz)]
   lrr <- as.matrix(Z) %*% cox$gamma
   cumhazard <- cox$cum
   rr <- exp(lrr)
   xid <- sample(1:nrow(Z),n,replace=TRUE)
   rr <- rr[xid]
   Z <- Z[xid,]
}

   ptt <- pc.hazard(cumhazard,rr)

   if (length(rr)==1) n<-rr else n <- length(rr)
   if (is.null(rrc)) {
	   if (length(rr)==1) rrc<-rr else rrc <- length(rr)
   }
   if (is.matrix(cens)) {
	   pct <- pc.hazard(cens,rr,cum.hazard=TRUE)
	   pct <- pct$time
   }
   else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
	   }
   }
   dt <- cbind(data.frame(time=pmin(ptt$time,pct),
			  status=ifelse(ptt$time<pct,ptt$status,0)),Z)

   return(dt)
}# }}}


