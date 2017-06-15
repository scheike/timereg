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
   names(dt) <- c("time","status")
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

#' Simulation of Piecewise constant hazard models with two causes (Cox).
#' 
#' Simulates data from piecwise constant baseline hazard that can also be of
#' Cox type. Censor data at highest value of the break points.
#' 
#' @param cumhaz1 cumulative hazard of cause 1
#' @param cumhaz2 cumulative hazard of cause 1
#' @param rr1 number of simulations or vector of relative risk for simuations.
#' @param rr2 number of simulations or vector of relative risk for simuations.
#' @param cens 
#' @param rrc retlativ risk for censoring.
#' @param cens.cum.hazard possible cumulative hazard for censoring.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#
#' library(timereg)
#'
#' library(mets)
#' data(TRACE)
#' dtable(TRACE,~status)
#' dsummary(TRACE,time~status)
#'
#' cox1 <- cox.aalen(Surv(time,status==9)~prop(vf)+prop(chf)+prop(wmi),
#'             data=TRACE,robust=0)
#' cox2 <-  cox.aalen(Surv(time,status==0)~prop(vf)+prop(chf)+prop(wmi),
#'             data=TRACE,robust=0)
#'
#' X1 <- TRACE[,c("vf","chf","wmi")]
#' n <- 100000
#' xid <- sample(1:nrow(X1),n,replace=TRUE)
#' Z1 <- X1[xid,]
#' Z2 <- X1[xid,]
#' rr1 <- exp(as.matrix(Z1) %*% cox1$gamma)
#' rr2 <- exp(as.matrix(Z2) %*% cox2$gamma)
#'
#' cumhaz1 <- cox1$cum
#' cumhaz2 <- cox2$cum
#' d <-  cause.pchazard.sim(cox1$cum,cox2$cum,rr1,rr2)
#' d <- cbind(dd,Z0)
#' c1 <-   cox.aalen(Surv(time,status==1)~prop(vf)+prop(chf)+prop(wmi),
#'                   data=dd,robust=0)
#' cbind(sc1$gamma, cox1$gamma)
#' c2 <-  cox.aalen(Surv(time,status==2)~prop(vf)+prop(chf)+prop(wmi),
#'                   data=dd,robust=0)
#' cbind(sc2$gamma, cox2$gamma)
#' par(mfrow=c(1,2))
#' plot(cox1); lines(sc1$cum,col=2)
#' plot(cox2$cum,type="l");
#' lines(sc2$cum,col=2)
#' @export
cause.pchazard.sim <- function(cumhaz1,cumhaz2,rr1,rr2,cens=NULL,rrc=NULL,cens.cum.hazard=TRUE)
{#'# {{{
##'## cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
##'## cumh=cbind(breaks,cumhazard) if cumh=TRUE
 
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
{# {{{
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
}# }}}
if (class(cox)=="cox.aalen")
{# {{{
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
}# }}}
if (class(cox)=="phreg")
{# {{{
###   formula <- attr(cox, "Formula")
######   Z1 <- na.omit(model.matrix(cox,data=data))
   p <- length(cox$coef)
###   Z <- na.omit(get_all_vars(formula,data=data))
   Z  <-  cox$model.matrx[,-1]
   nz <- ncol(Z)
   lrr <- as.matrix(Z) %*% cox$coef
   cumhazard <- cox$cumhaz
   rr <- exp(lrr)
   xid <- sample(1:nrow(Z),n,replace=TRUE)
   rr <- rr[xid]
   Z <- Z[xid,]
}# }}}

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


#' Simulation of cause specific from Cox models.
#' 
#' Simulates data that looks like fit from cause specific Cox models. 
#' Censor data automatically. When censoring is given in the  list of causes this
#' will give censoring that looks like the data.  Covariates are drawn from data-set
#' with replacement. This gives covariates like the data. 
#' 
#' 
#' @param coxs list of cox models.
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed
#' covariates).
#' @param cens specifies censoring model, if NULL then only censoring for 
#'   	       each cause at end of last event of this type. 
#' 	       if "is.matrix" then uses cumulative. 
#'             hazard given, if "is.scalar" then uses rate for exponential, and if not
#'             given then takes average rate of in simulated data from cox model.
#'             But censoring can also be given as a cause.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' data(TRACE)
#' cox1 <- cox.aalen(Surv(time,status==9)~prop(vf)+prop(chf)+prop(wmi),
#'              data=TRACE,robust=0)
#'  cox2 <-  cox.aalen(Surv(time,status==0)~prop(vf)+prop(chf)+prop(wmi),
#'              data=TRACE,robust=0)
#'  coxs <- list(cox1,cox2)
#' 
#'  dd <- sim.cause.cox(coxs,1000,data=TRACE)
#'  scox1 <- cox.aalen(Surv(time,status==1)~prop(vf)+prop(chf)+prop(wmi),
#'              data=dd,robust=0)
#'  scox2 <-  cox.aalen(Surv(time,status==2)~prop(vf)+prop(chf)+prop(wmi),
#'              data=dd,robust=0)
#'  ### 
#'  cbind(cox1$gamma,scox1$gamma)
#'  cbind(cox2$gamma,scox2$gamma)
#'  par(mfrow=c(1,2))
#'  plot(cox1); lines(scox1$cum,col=2)
#'  plot(cox2$cum,type="l");
#'  lines(scox2$cum,col=2)
#' 
#' cox1 <- coxph(Surv(time,status==9)~vf+chf+wmi, data=TRACE)
#' cox2 <- coxph(Surv(time,status==0)~vf+chf+wmi, data=TRACE)
#'  
#' coxs <- list(cox1,cox2)
#' dd <- sim.cause.cox(coxs,100000,data=TRACE)
#' scox1 <- cox.aalen(Surv(time,status==1)~prop(vf)+prop(chf)+prop(wmi),
#'               data=dd,robust=0)
#' scox2 <-  cox.aalen(Surv(time,status==2)~prop(vf)+prop(chf)+prop(wmi),
#'               data=dd,robust=0)
#' cbind(coef(cox1),scox1$gamma)
#' cbind(coef(cox2),scox2$gamma)
#'  
#' @export
sim.cause.cox <- function(coxs,n,data=NULL,cens=NULL,rrc=NULL)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE

if (!is.list(coxs)) stop("Cox models in list form\n"); 

if (class(coxs[[1]])=="coxph")
{# {{{
   sims <- c()
   simscovs <- c()
   cox <- coxs[[1]]

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
  simcovs <- Z
  ptt <- pc.hazard(cumhazard,rr)

  i=2
   if (length(coxs)>=2)
   for (i in 2:length(coxs)) {
   cox <- coxs[[i]]

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
###  xid <- sample(1:nrow(Z),n,replace=TRUE)
  rr <- rr[xid]
  colnames(Z) <- nn

  Z1 <- as.matrix(Z[xid,])
  ptt1 <- pc.hazard(cumhazard,rr)
  ### removes covariates already included 
  rmn <- match(names(Z1),names(simcovs))
  Z1e <- Z1[,-rmn]
  simcovs <- cbind(simcovs,Z1e)

  ptt <- data.frame(time=pmin(ptt$time,ptt1$time),status=ifelse(ptt$time<=ptt1$time,ptt$status,ptt1$status*i))
   }

dt <- cbind(ptt,simcovs)

}# }}}

if (class(coxs[[1]])=="cox.aalen")
{# {{{

   simcovs <- c()
   cox <- coxs[[1]]
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
   ptt <- pc.hazard(cumhazard,rr)
   simcovs <- Z

   if (length(coxs)>=2)
   for (i in 2:length(coxs)) {# {{{
   cox <- coxs[[i]]
   formula <- attr(cox, "Formula")
###   Z1 <- na.omit(model.matrix(cox,data=data))
   p <- nrow(cox$gamma)
   Z <- na.omit(get_all_vars(formula,data=data))
   nz <- ncol(Z)
   Z <- Z[,seq(nz-p+1,nz)]
   lrr <- as.matrix(Z) %*% cox$gamma
   cumhazard <- cox$cum
   rr <- exp(lrr)
   rr <- rr[xid]
   Z1 <- Z[xid,]
   ptt1 <- pc.hazard(cumhazard,rr)
   ### removes covariates already included 
   rmn <- match(names(Z1),names(simcovs))
   Z1e <- Z1[,-rmn]
   simcovs <- cbind(simcovs,Z1e)

   ptt <- data.frame(time=pmin(ptt$time,ptt1$time),status=ifelse(ptt$time<=ptt1$time,ptt$status,ptt1$status*i))
   }# }}}

dt <- cbind(ptt,simcovs)

}# }}}


   if (length(rr)==1) n<-rr else n <- length(rr)
   if (is.null(rrc)) {
	   if (length(rr)==1) rrc<-rr else rrc <- length(rr)
   }
   if (!is.null(cens)) {
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
			  status=ifelse(ptt$time<pct,ptt$status,0)),simcovs)
   }


   return(dt)
}# }}}

