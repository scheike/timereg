#' Simulation of Piecewise constant hazard model (Cox).
#' 
#' Simulates data from piecwise constant baseline hazard that can also be of
#' Cox type. Censor data at highest value of the break points.
#' 
#' For a piecewise linear cumulative hazard the inverse is easy to compute with 
#' and delayed entry x we compute 
#' \deqn{\Lambda^{-1}(\Lambda(x) + E/RR)}, 
#' where RR are the relative risks and E is exponential with mean 1.
#' This quantity has survival function 
#' \deqn{P(T > t | T>x) = exp(-RR (\Lambda(t) - \Lambda(x)))}. 
#' 
#' @param cumhazard cumulative hazard, or piece-constant rates for periods defined by first column of input.
#' @param rr relative risk for simulations, alternatively when rr=1 specify n 
#' @param n number of simulation if rr not given 
#' @param entry delayed entry time for simuations.
#' @param cum.hazard specifies wheter input is cumulative hazard or rates.
#' @param cause name of cause 
#' @param extend  to extend piecewise constant with constant rate. Default is average rate over time from cumulative (when TRUE), if numeric then uses given rate.
#' @author Thomas Scheike
#' @keywords survival
#' @examples
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
#' pctime <- rchaz(haz,n=1000,cum.hazard=FALSE)
#' 
#' par(mfrow=c(1,2))
#' ss <- aalen(Surv(time,status)~+1,data=pctime,robust=0)
#' plot(ss)
#' lines(cumhaz,col=2,lwd=2)
#' 
#' pctimecox <- rchaz(cumhaz,rrcox)
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
#' pctime <- rchaz(ss$cum,n=1000)
#' ###
#' sss <- aalen(Surv(time,status)~+1,data=pctime,robust=0)
#' lines(sss$cum,col=2,lwd=2)
#' 
#' pctime <- rchaz(ss$cum,rrcox)
#' pctime <- cbind(pctime,X)
#' ###
#' sss <- cox.aalen(Surv(time,status)~+prop(X),data=pctime,robust=0)
#' summary(sss)
#' plot(ss)
#' lines(sss$cum,col=3,lwd=3)
#' 
#' @export 
#' @aliases pc.hazard simrchaz addCums
rchaz <- function(cumhazard,rr,n=NULL,entry=NULL,cum.hazard=TRUE,cause=1,extend=FALSE)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE
  if (!is.null(n)) rr <- rep(1,n)

  breaks <- cumhazard[,1]
  rates <- cumhazard[,2][-1]
  mm <- tail(breaks,1)
  if (cum.hazard==FALSE) {
        cumh <- cumsum(c(0,diff(breaks)*rates))
        cumhazard <- cbind(breaks,cumh)
  } else cumh <- cumhazard[,2] 
   n <- length(rr)
   ttt <- rexp(n)/rr
   if (cumhazard[1,2]>0)  { ## start cumulative hazard with a 0
	   warning("Safest to start with cumulative hazard 0 to avoid problems\n"); 
	   cumhazard <- rbind(c(0,0),cumhazard)
	   cumh <- c(0,cumh)
   }
   ###
   if (!is.null(entry)) {
	   if (length(entry)==1) entry <- rep(entry,n) else entry <- entry
	   cumentry <- lin.approx(entry,cumhazard,x=1)
   } else { entry <- cumentry <- rep(0,n) }
   ###
   ttte <- ttt+cumentry
   rrx <- lin.approx(ttte,cumhazard,x=-1)
   rrx <- ifelse(rrx>mm,mm,rrx)
   status <- rep(0,n)
   status <- ifelse(rrx<mm,cause,status)
   tcum <- tail(cumhazard,1)
   extend.rate <- NULL
   if (is.logical(extend))
	   if (extend) extend.rate <- tcum[2]/tcum[1]
   if (is.numeric(extend)) { extend.rate <- extend; extend <- TRUE;}
   if (extend) {
	   whoe <- which(status==0)
	   nn <- length(whoe)
	   if (nn>0) {
	     textend <- rexp(nn)/(rr[whoe]*extend.rate)
	     rrx[whoe] <- mm+textend
	     status[whoe] <- 1
	   }
   }
   dt <- data.frame(entry=entry,time=rrx,status=status,rr=rr)
   colnames(dt) <- c("entry","time","status","rr"); 
   attr(dt,"cumhaz") <- cumhazard
   attr(dt,"extend.rate") <- extend.rate
   return(dt)
}# }}}

pc.hazard <- function(cumhazard,rr,...)
{# {{{
rchaz(cumhazard,rr,...)
}# }}}

#' @export
lin.approx <- function(x2,xfx,x=1)
{# {{{
   ### x=1   gives  f(x2) 
   ### x=-1  gives  f^-1(x2) 
   breaks <- xfx[,x]
   fx     <- xfx[,-x]
   ri <- sindex.prodlim(breaks,x2)
   maxindex <- which(ri==length(breaks))
   rip1 <- ri+1
   rip1[maxindex] <- length(breaks)
   rrr <- (x2-breaks[ri])/(breaks[rip1]-breaks[ri])
   rrr[maxindex] <- 0
   res <- rrr*(fx[rip1]-fx[ri])+fx[ri]
   res[is.na(res)] <- tail(fx,1)
   return(res)
}# }}}

#' @export
addCums <- function(cumB,cumA,max=NULL)
{# {{{
 ## take times
 times <- sort(unique(c(cumB[,1],cumA[,1])))
 if (!is.null(max)) times <- times[times<max]
 cumBjx <- lin.approx(times,cumB,x=1)
 cumAjx <- lin.approx(times,cumA,x=1)
 cumBA <- cumBjx+cumAjx
 return(cbind(times,cumBA))
}# }}}

#' @export
simrchaz <- function(cumhazard,rr,n=NULL,cens=NULL,rrc=NULL,...)
{# {{{
###   adds censoring to to rchaz call
   if (!is.null(n)) rr <- rep(1,n)
   n <- length(rr)  

   ptt <- rchaz(cumhazard,rr,cum.hazard=TRUE,...)

   if (is.null(rrc)) {
	   rrc <- rep(1,n)
   }
   if (is.matrix(cens)) {
	   pct <- rchaz(cens,rrc,...)
	   pct <- pct$time
   } else {
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
#' Cox type. Censor data at highest value of the break points for either of the
#' cumulatives.
#' 
#' @param cumhaz1 cumulative hazard of cause 1
#' @param cumhaz2 cumulative hazard of cause 1
#' @param rr1 number of simulations or vector of relative risk for simuations.
#' @param rr2 number of simulations or vector of relative risk for simuations.
#' @param n number of simulation if rr not given 
#' @param cens to censor further , rate or cumumlative hazard
#' @param rrc retlativ risk for censoring.
#' @param ... arguments for rchaz 
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' data(TRACE)
#'
#' cox1 <- cox.aalen(Surv(time,status==9)~prop(vf)+prop(chf)+prop(wmi),
#'             data=TRACE,robust=0)
#' cox2 <-  cox.aalen(Surv(time,status==0)~prop(vf)+prop(chf)+prop(wmi),
#'             data=TRACE,robust=0)
#'
#' X1 <- TRACE[,c("vf","chf","wmi")]
#' n <- 1000
#' xid <- sample(1:nrow(X1),n,replace=TRUE)
#' Z1 <- X1[xid,]
#' Z2 <- X1[xid,]
#' rr1 <- exp(as.matrix(Z1) %*% cox1$gamma)
#' rr2 <- exp(as.matrix(Z2) %*% cox2$gamma)
#'
#' cumhaz1 <- cox1$cum
#' cumhaz2 <- cox2$cum
#' d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2)
#' dd <- cbind(d,Z1)
#' sc1 <-   cox.aalen(Surv(time,status==1)~prop(vf)+prop(chf)+prop(wmi),
#'                   data=dd,robust=0)
#' cbind(sc1$gamma, cox1$gamma)
#' sc2 <-  cox.aalen(Surv(time,status==2)~prop(vf)+prop(chf)+prop(wmi),
#'                   data=dd,robust=0)
#' cbind(sc2$gamma, cox2$gamma)
#' par(mfrow=c(1,2))
#' plot(cox1); lines(sc1$cum,col=2)
#' plot(cox2$cum,type="l");
#' lines(sc2$cum,col=2)
#' 
#' @export 
#' @aliases cause.pchazard.sim 
rcrisk <-function(cumhaz1,cumhaz2,rr1,rr2,n=NULL,cens=NULL,rrc=NULL,...)
{#'# {{{
 
if (!is.null(n)) { rr1 <- rep(1,n); rr2 <- rep(1,n) }
n <- length(rr1); 
if (missing(rr2)) rr2 <-rep(1,n)

ptt1 <- rchaz(cumhaz1,rr1,cum.hazard=TRUE,...)
ptt2 <- rchaz(cumhaz2,rr2,cum.hazard=TRUE,...)
ptt <- data.frame(time=pmin(ptt1$time,ptt2$time),status=ifelse(ptt1$time<=ptt2$time,ptt1$status,ptt2$status*2),entry=ptt1$entry)

if (!is.null(cens)) {
	if (is.null(rrc)) rrc <- rep(1,n)
	if (is.matrix(cens)) {
		pct <- rchaz(cens,rrc,...)
		pct <- pct$time
	} else {
		if (is.numeric(cens)) pct<- rexp(n)/cens  else {
			chaz <-sum(ptt1$status+ptt2$status)/sum(ptt1$time+ptt2$time)  ## hazard averate T haz 
			pct<- rexp(n)/chaz 
		}
	}
	ptt <- data.frame(time=pmin(ptt$time,pct),status=ifelse(ptt$time<pct,ptt$status,0),entry=ptt$entry)
}

return(ptt)
}# }}}

#' @export 
cause.pchazard.sim<-function(cumhaz1,cumhaz2,rr1,rr2,cens=NULL,rrc=NULL,...)
{
	rcrisk(cumhaz1,cumhaz2,rr1,rr2,cens=NULL,rrc=NULL,...)
}

#' Simulation of output from Cox model.
#' 
#' Simulates data that looks like fit from Cox model. Censor data automatically
#' for highest value of the event times by using cumulative hazard. 
#' 
#' @param cox output form coxph or cox.aalen model fitting cox model.
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed
#' covariates).
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' given then takes average rate of in simulated data from cox model.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param entry delayed entry variable for simulation.
#' @param ... arguments for rchaz, for example entry-time
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' data(TRACE)
#' cox <-  coxph(Surv(time,status==9)~vf+chf+wmi,data=TRACE)
#' sim1 <- sim.cox(cox,1000,data=TRACE)
#' cc <- coxph(Surv(time,status)~vf+chf+wmi,data=sim1)
#' cbind(cox$coef,cc$coef)
#' 
#' cor(sim1[,c("vf","chf","wmi")])
#' cor(TRACE[,c("vf","chf","wmi")])
#' 
#' \donttest{
#' ### do not test to avoid dependence on mets
#' library(mets)
#' cox <-  phreg(Surv(time, status==9)~vf+chf+wmi,data=TRACE)
#' sim3 <- sim.cox(cox,1000,data=TRACE)
#' cc <-  phreg(Surv(time, status)~vf+chf+wmi,data=sim3)
#' cbind(cox$coef,cc$coef)
#' basehazplot.phreg(cox,se=TRUE)
#' lines(cc$cumhaz,col=2)
#' 
#' cox <-  phreg(Surv(time,status==9)~strata(chf)+vf+wmi,data=TRACE)
#' sim3 <- sim.cox(cox,1000,data=TRACE)
#' cc <-   phreg(Surv(time, status)~strata(chf)+vf+wmi,data=sim3)
#' cbind(cox$coef,cc$coef)
#' basehazplot.phreg(cox)
#' basehazplot.phreg(cc,add=TRUE)
#' }
#' 
#' @export
#' @aliases sim.cox read.fit  sim.base
sim.cox <- function(cox,n,data=NULL,cens=NULL,rrc=NULL,entry=NULL,...)
{# {{{
	### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
	### cumh=cbind(breaks,cumhazard) if cumh=TRUE
	des <- 	read.fit(cox,n,data=data,...)
	Z <- des$Z; cumhazard <- des$cum; rr <- des$rr; 

	if (!inherits(cox,"phreg")) {
		if (cumhazard[1,2]>0) cumhazard <- rbind(c(0,0),cumhazard)
		ptt <- rchaz(cumhazard,rr,entry=entry) 
		ptt <- cbind(ptt,Z)
	} else {
		ptt <- data.frame()
		stratj <- cox$strata[cox$jumps]
		if (cox$nstrata>1) {
			for (j in 0:(cox$nstrata-1)) {
				cumhazardj <- rbind(c(0,0),cox$cumhaz[stratj==j,])
				if (!is.null(entry)) entry <- entry[des$strataid==j]
				pttj <- rchaz(cumhazardj,rr[des$strataid==j],entry=entry) 
				Zj <- Z[des$strataid==j,,drop=FALSE]
				pttj <- cbind(pttj,Zj)
				ptt  <-  rbind(ptt,pttj)
			}
		} else {
			if (cumhazard[1,2]>0) cumhazard <- rbind(c(0,0),cumhazard)
			ptt <- rchaz(cumhazard,rr,entry=entry) 
			ptt <- cbind(ptt,Z)
		}
	}

	if (!is.null(cens))  {
		if (is.null(rrc)) rrc <- rep(1,n)
		if (is.matrix(cens)) {
			pct <- rchaz(cens,rrc,entry=entry)
			pct <- pct$time
		}
		else {
			if (is.numeric(cens)) pct<- rexp(n)/cens  else {
				chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
				pct<- rexp(n)/chaz 
				if (!is.null(entry)) pct  <- entry + pct
			}
		}
		ptt$time <- pmin(ptt$time,pct)
		ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
	} 

	attr(ptt,"id") <- des$id
	return(ptt)
}# }}}

#' @export
sim.base <- function(cumhaz,n,data=NULL,strata=NULL,cens=NULL,entry=NULL,...)
{# {{{
	### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
	### cumh=cbind(breaks,cumhazard) if cumh=TRUE
	if (!is.null(strata)) {
		us <- unique(strata)
		nus <- length(us)
		if (length(n)!=nus) warning("n must be given for each strata as vector \n")
		if (length(n)!=nus) n <- rep(round(n/nus),nus)
	}

	nt <- 0
	if (is.null(strata))  {
		if (cumhaz[1,2]>0) cumhaz <- rbind(c(0,0),cumhaz)
		ptt <- rchaz(cumhaz,n=n,entry=entry) 
		nt <- n
	} else {
		ptt <- c()
		nstrata <- length(unique(strata)) 
			for (i in seq_along(unique(strata))) {
				j <- unique(strata)[i]
				cumhazardj <- rbind(c(0,0),cumhaz[strata==j,])
				if (!is.null(entry)) entry <- entry[strata==j]
				ns <- n[i]
				nt <- nt+ns
				pttj <- cbind(rchaz(cumhazardj,n=ns,entry=entry),j)
				colnames(pttj)[5] <- "strata"
				ptt  <-  rbind(ptt,pttj)
			}
	} 

	if (!is.null(cens))  {
		if (is.matrix(cens)) {
			pct <- rchaz(cens,n=nt,entry=entry)$time
		}
		else {
				pct<- rexp(nt)/(cens) 
	             }
		ptt$time <- pmin(ptt$time,pct)
		ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
	} 

	return(ptt)
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
#' @param ... arguments for rchaz, for example entry-time
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' nsim <- 100
#' 
#' data(bmt)
#' # coxph          
#' cox1 <- coxph(Surv(time,cause==1)~tcell+platelet,data=bmt)
#' cox2 <- coxph(Surv(time,cause==2)~tcell+platelet,data=bmt)
#' coxs <- list(cox1,cox2)
#' dd <- sim.cause.cox(coxs,nsim,data=bmt)
#' scox1 <- coxph(Surv(time,status==1)~tcell+platelet,data=dd)
#' scox2 <- coxph(Surv(time,status==2)~tcell+platelet,data=dd)
#' cbind(cox1$coef,scox1$coef)
#' cbind(cox2$coef,scox2$coef)
#' 
#' \donttest{
#' ### do not test to avoid dependence on mets
#' library(mets)          
#' data(bmt)
#' cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
#' cox2 <- phreg(Surv(time,cause==2)~tcell+platelet,data=bmt)
#' coxs <- list(cox1,cox2)
#' dd <- sim.cause.cox(coxs,nsim,data=bmt)
#' scox1 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
#' scox2 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
#' cbind(cox1$coef,scox1$coef)
#' cbind(cox2$coef,scox2$coef)
#' par(mfrow=c(1,2))
#' basehazplot.phreg(cox1); basehazplot.phreg(scox1,add=TRUE); 
#' basehazplot.phreg(cox2); basehazplot.phreg(scox2,add=TRUE); 
#'
#' cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=bmt)
#' cox2 <- phreg(Surv(time,cause==2)~strata(tcell)+platelet,data=bmt)
#' coxs <- list(cox1,cox2)
#' dd <- sim.cause.cox(coxs,nsim,data=bmt)
#' scox1 <- phreg(Surv(time,status==1)~strata(tcell)+platelet,data=dd)
#' scox2 <- phreg(Surv(time,status==2)~strata(tcell)+platelet,data=dd)
#' cbind(cox1$coef,scox1$coef)
#' cbind(cox2$coef,scox2$coef)
#' par(mfrow=c(1,2))
#' basehazplot.phreg(cox1); basehazplot.phreg(scox1,add=TRUE); 
#' basehazplot.phreg(cox2); basehazplot.phreg(scox2,add=TRUE); 
#' 
#' }  
#' @export
sim.cause.cox <- function(coxs,n,data=NULL,cens=NULL,rrc=NULL,...)
{# {{{
### cumh=cbind(breaks,rates), first rate is 0 if cumh=FALSE
### cumh=cbind(breaks,cumhazard) if cumh=TRUE

if (!is.list(coxs)) stop("Cox models in list form\n"); 

  ptt <- sim.cox(coxs[[1]],n,data=data)
  simcovs <- ptt[,(5:ncol(ptt))]

  i=2
  if (length(coxs)>=2)
  for (i in 2:length(coxs)) {
  cox <- coxs[[i]]
  ptt1 <- sim.cox(coxs[[i]],n,data=data,id=attr(ptt,"id"))
  Zi <- ptt1[,(5:ncol(ptt1))]
  pm <- match(names(Zi),names(simcovs))
  Zi <- Zi[,-pm]
  simcovs <- cbind(simcovs,Zi)
  dt <- data.frame(time=pmin(ptt$time,ptt1$time),
                    status=ifelse(ptt$time<=ptt1$time,ptt$status,ptt1$status*i))
  }

  dt <- cbind(dt,simcovs)

   if (!is.null(cens)) {
   if (is.matrix(cens)) {
	   pct <- rchaz(cens,rrc)
	   pct <- pct$time
   }
   else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
	   }
   }

   dt <- cbind(data.frame(time=pmin(ptt$time,pct),status=ifelse(ptt$time<pct,ptt$status,0)),simcovs)
   }

   return(dt)
}# }}}


#' Simulation from subdistribution function assuming piecwise linearity 
#' 
#' Simulation from subdistribution function assuming piecwise linearity for Fine-Gray or logistic link.
#' 
#' @param cumhazard matrix that specified times and values of some cumulative hazard.
#' @param rr "relative risk" terms 
#' @param n number of simulation if rr not given 
#' @param entry not implemented yet 
#' @param type either cloglog or logistic 
#' @param startcum c(0,0) to make cumulativ start in 0 with 0 cumhazard.
#' @param ... further arguments 
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' data(sTRACE)
#' 
#' cif <- comp.risk(Event(time,status)~const(vf),data=sTRACE,cause=9,model="logistic2")
#' cumhaz <- cif$cum
#' 
#' ## 1000 logistic without covariates with baseline from model fit  
#' sim1 <- simsubdist(cumhaz,n=1000,type="logistic")
#' ###
#' cifs <- comp.risk(Event(time,status)~+1,data=sim1,cause=1,model="logistic2")
#' ###
#' plot(cifs)
#' lines(cifs$cum,col=2)
#'
#' ## 1000 logistic with covariates with baseline from model fit  
#' x <- rbinom(1000,1,0.5)
#' rr <- exp(x*0.3)
#' sim1 <- simsubdist(cumhaz,rr,type="logistic")
#' sim1$x <- x
#'
#' cifs <- comp.risk(Event(time,status)~+const(x),data=sim1,cause=1,model="logistic2")
#' ###
#' cifs$gamma
#' plot(cifs)
#' lines(cumhaz,col=2)
#'
#' ##################################################################
#' ### simulation of cumulative incidence with specified functions
#' ##################################################################
#' F1logit<-function(t,lam0=0.2,beta=0.3,x=0) 
#' { pt <- t*lam0; rr <- exp(x*beta); return(pt*rr/(1+pt*rr)); } 
#' 
#' F1p<-function(t,lam0=0.4,beta=0.3,x=0) # proportional version 
#' { return( 1 - exp(-(t*lam0)*exp(x*beta))) }
#' 
#' n=1000
#' tt=seq(0,3,by=.01)
#' tt=seq(0,3,by=.01)
#' t1 <- invsubdist(cbind(tt,F1p(tt)),runif(n))
#' t2 <- invsubdist(cbind(tt,F1p(tt,lam0=0.1)),runif(n))
#' rt <- rbinom(n,1,(F1p(3)+F1p(3,lam0=0.1)))
#' rb <- rbinom(n,1,F1p(3)/(F1p(3)+F1p(3,lam0=0.1)))
#' cause=ifelse(rb==1,1,2)
#' time=ifelse(cause==1,t1$time,t2$time)
#' cause <- rt*cause
#' time[cause==0] <- 3
#' 
#' datC=data.frame(time=time,cause=cause)
#' p1=comp.risk(Event(time,cause)~+1,data=datC,cause=1)
#' p2=comp.risk(Event(time,cause)~+1,data=datC,cause=2)
#' pp1=predict(p1,X=1,se=0)
#' pp2=predict(p2,X=1,se=0)
#' par(mfrow=c(1,2))
#' plot(pp1)
#' lines(tt,F1p(tt),col=2)
#' plot(pp2)
#' lines(tt,F1p(tt,lam0=0.1),col=2)
#' 
#' #to avoid dependencies when checking 
#' #library(prodlim)
#' #pp=prodlim(Hist(time,cause)~+1)
#' #par(mfrow=c(1,2))
#' #plot(pp,cause="1")
#' #lines(tt,F1p(tt),col=2)
#' #plot(pp,cause="2")
#' #lines(tt,F1p(tt,lam0=0.1),col=2)
#' @aliases simsubdist
#' @export
simsubdist <- function(cumhazard,rr,n=NULL,entry=NULL,type="cloglog",startcum=c(0,0),...)
{# {{{
  ## Fine-Gray model cloglog F1= 1-exp(-cum(t)*rr)
  ## logistic                F1= cum(t)*rr/(1+cum(t)*rr)
  ## identity                F1= cum(t)*rr
  ## rr=exp(X^t beta) 
  if (!is.null(n)) rr <- rep(1,n)
  entry=NULL

  breaks <- cumhazard[,1]
  rates <- cumhazard[,2][-1]
  mm <- tail(breaks,1)
  cumh <- cumhazard[,2] 
  n <- length(rr)

  if (cumh[1]>0) {
	  cumh <-    c(startcum[2],cumh)
	  breaks <-  c(startcum[1],breaks)
  }
  if (type=="cloglog") {
      F1tau <- 1-exp(-tail(cumh,1)*rr)
      ttt <- -log(1-runif(n)*F1tau)/rr
  } else if (type=="logistic") {
     F1tau <- tail(cumh,1)*rr/(1+tail(cumh,1)*rr)
     v <- runif(n)*F1tau
     ttt <- exp(logit(v))/rr; 
  }  else if (type=="identity" | type=="cif") {
     F1tau <- tail(cumh,1)
     ttt <- runif(n)*F1tau
     ## rr only affects binomial draw 
     F1tau <- F1tau*rr 
  } else stop(" cloglog or logistic or give function (fun=) \n"); 
  ###
   entry <- cumentry <- rep(0,n)
   ttte <- ttt+cumentry
   rrx <- lin.approx(ttt,cbind(breaks,cumh),x=-1)
   timecause <- rrx
   ###
   rrx <- ifelse(rrx>mm,mm,rrx)
   status <- rbinom(n,1,F1tau) 
   rrx[status==0] <- mm
   dt <- data.frame(entry=entry,time=rrx,status=status,rr=rr,F1tau=F1tau,timecause=timecause)
   attr(dt,"cumhaz") <- cumhazard
   return(dt)
}# }}}

#' Finds inverse of piecwise linear sub-distribution 
#' 
#' Finds inverse of piecwise linear sub-distribution to be used for simulation of subdistributions
#' 
#' @param F1 matrix with x, and F1(x) 
#' @param u points for which to compute inverse
#' @param entry possible delayed entry points
#' @param cond 1 indcates that we draw given that this subdistribution is used, so draws from F1(t)/F1(tau) 
#' @param ptrunc possible trunction weigth for delayed entry, if NULL then uses ptrunc=1-F1(entry)
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' F1 <- cbind(c(0,5,8,10),c(0,0.1,0.3,0.9))
#' plot(F1,type="l")
#' u <- runif(100)
#' Fiu <- invsubdist(F1,u,cond=0)
#' points(Fiu$time,u,pch="x")
#' 
#' F1cond <- F1
#' F1cond[,2] <- F1cond[,2]/0.9
#' plot(F1cond,type="l")
#' u <- runif(100)
#' Ficond <- invsubdist(F1cond,u,cond=0)
#' points(Ficond$time,u,pch="-")
#' Fiu <- invsubdist(F1,u,cond=1)
#' points(Fiu$time,u,pch="x")
#' 
#' 
#' ## simulation of distribution with delayed entry starting at 3
#' par(mfrow=c(1,1))
#' F1 <- cbind(c(0,5,8,10),c(0,0.5,0.6,0.9))
#' F1
#' plot(F1,ylim=c(0,1),type="l")
#' 
#' n <- 100000
#' entry <- c(rep(3,10000),runif(n)*7+3)
#' ###entry <- rep(3,n)
#' u <- runif(n+10000)
#' ###
#' Fiu <- invsubdist(F1,u,cond=0,entry=entry)
#' ###
#' # library(prodlim)
#' # pp <- prodlim(Hist(time,status,entry=entry)~+1,data=Fiu)
#' # plot(pp,xlim=c(3,10))
#' ###
#' entry <- 3
#' ###
#' F1entry <- subdist(F1,entry)[,2]
#' ptrunc <- 1-F1entry
#' ###
#' F1entry5 <- F1
#' F1entry5[,1] <- F1entry5[,1]-entry
#' F1entry5[,2] <- (F1entry5[,2]-F1entry)/ptrunc
#' pos <- F1entry5[,1]>=0
#' F1entry5 <- rbind(c(0,0),F1entry5[pos,])
#' #
#' # lines(entry+F1entry5[,1],1-F1entry5[,2],col=2)
#' 
#' ##############################################################
#' ## Simulations of two cumulative incidence functions with truncation 
#' ##############################################################
#' par(mfrow=c(1,1))
#' F1 <- cbind(c(0,5,8,10),c(0,0.5,0.6,0.9)*0.3)
#' F2 <- cbind(c(0,5,8,10),c(0,0.5,0.6,0.9)*0.5)
#' plot(F1,ylim=c(0,1),type="l")
#' lines(F2,col=2)
#' entry1 <- 3
#' ###
#' F1entry <- subdist(F1,entry1)[,2]
#' F2entry <- subdist(F2,entry1)[,2]
#' ptrunc <- 1-F1entry-F2entry
#' ###
#' F1e <- F1
#' F1e[,1] <- F1e[,1]-entry1
#' F1e[,2] <- (F1e[,2]-F1entry)/ptrunc
#' pos <- F1e[,1]>=0
#' F1e <- rbind(c(0,0),F1e[pos,])
#' F2e <- F2
#' F2e[,1] <- F2e[,1]-entry1
#' F2e[,2] <- (F2e[,2]-F2entry)/ptrunc
#' pos <- F2e[,1]>=0
#' F2e <- rbind(c(0,0),F2e[pos,])
#' #
#' # truncated identifiable version 
#' lines(entry1+F1e[,1],F1e[,2],col=1)
#' lines(entry1+F2e[,1],F2e[,2],col=2)
#' 
#' n <- 10000
#' entry <- c(rep(entry1,10000),runif(n)*(10-entry1)+entry1)
#' u <- runif(n+10000)
#' ###
#' F1entry <- subdist(F1,entry)[,2]
#' F2entry <- subdist(F2,entry)[,2]
#' ptrunc <- 1-( F1entry+F2entry)
#' Fiu1 <- invsubdist(F1,u,cond=1,entry=entry,ptrunc=ptrunc)
#' Fiu2 <- invsubdist(F1,u,cond=1,entry=entry,ptrunc=ptrunc)
#' ###
#' ptot <- (tail(F1[,2],1)+tail(F2[,2],1)-F1entry-F2entry)/(ptrunc)
#' rt <- rbinom(n+10000,1,ptot)
#' p1 <- ((tail(F1[,2],1)-F1entry)/ptrunc)
#' p2 <- ((tail(F2[,2],1)-F2entry)/ptrunc)
#' rb <- rbinom(n+10000,1,p1/ptot)
#' cause=ifelse(rb==1,1,2)
#' time=ifelse(cause==1,Fiu1$time,Fiu2$time)
#' cause <- rt*cause
#' time[cause==0] <- 10
#' ### simulated data, now checking that things are working 
#' # pp <- prodlim(Hist(time,cause,entry=entry)~+1)
#' # plot(pp,xlim=c(entry1,10),cause=1)
#' # plot(pp,xlim=c(entry1,10),cause=2,add=TRUE)
#' ###
#' # lines(entry1+F1e[,1],F1e[,2],col=2)
#' # lines(entry1+F2e[,1],F2e[,2],col=2)
#' 
#' @export
invsubdist <- function(F1,u,entry=NULL,cond=1,ptrunc=NULL)
{# {{{
  ### computes inverse subdist F1^{-1}(u) 
  ### computes inverse subdist F1^{-1}(u)
  F1max <- tail(F1[,2],1)
  mmax <- tail(F1[,1],1)

  if (!is.null(entry))  F1entry <- subdist(F1,entry)[,2];  
  if (!is.null(entry)) if (is.null(ptrunc)) ptrunc <- 1- F1entry
  if (cond==0 & !is.null(entry))  u<-u*ptrunc+F1entry
  if (cond==1 & !is.null(entry))  u<-u*(F1max-F1entry)+F1entry 
  if (cond==1 & is.null(entry))   u<-u*F1max
  rrx  <- lin.approx(u,F1,x=-1)
  status <- rep(1,length(u))
  status[u>=F1max] <- 0
  dt <- data.frame(time=rrx,status=status,y=u)
  if (!is.null(entry)) dt <- cbind(dt,entry)
  attr(dt,"F1") <- F1
  attr(dt,"Fmax") <- F1max
  return(dt)
}# }}}

#' @export
subdist <- function(F1,times)
{# {{{
  rrx <-   lin.approx(times,F1,x=1)
  dt <- cbind(times,rrx)
  colnames(dt) <- c("time","subdist")
  return(dt)
}# }}}

#' Simulation of output from Cumulative incidence regression model 
#' 
#' Simulates data that looks like fit from fitted cumulative incidence model
#' 
#' @param cif output form prop.odds.subdist or ccr (cmprsk), can also call invsubdist with
#'  with cumulative and linear predictor 
#' @param n number of simulations.
#' @param data to extract covariates for simulations (draws from observed
#' covariates).
#' @param Z to use these covariates for simulation rather than drawing new ones. 
#' @param drawZ to random sample from Z or not 
#' @param cens specifies censoring model, if "is.matrix" then uses cumulative
#' hazard given, if "is.scalar" then uses rate for exponential, and if not
#' given then takes average rate of in simulated data from cox model.
#' @param rrc possible vector of relative risk for cox-type censoring.
#' @param cumstart to start cumulatives at time 0 in 0. 
#' @param ... arguments for invsubdist
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' data(TRACE)
#' 
#' ## Logit link for proportional odds model, using comp.risk to save time 
#' #' cif <-  prop.odds.subdist(Event(time,status)~vf+chf+wmi,data=TRACE,cause=9)
#' cif <-  comp.risk(Event(time,status)~const(vf)+const(chf)+const(wmi),
#'                   data=TRACE,cause=9,model="logistic2")
#' sim1 <- sim.cif(cif,100,data=TRACE)
#' #' cc <-  prop.odds.subdist(Event(time,status)~vf+chf+wmi,data=sim1,cause=1)
#' cc <-  comp.risk(Event(time,status)~const(vf)+const(chf)+const(wmi),
#'                   data=sim1,cause=1,model="logistic2")
#' cbind(cif$gamma,cc$gamma)
#' 
#' plot(cif) 
#' lines(cc$cum)
#' 
#' #################################################################
#' ## Fine-Gray model model, using comp.risk to avoid dependcies
#' #################################################################
#' cif <-  comp.risk(Event(time,status)~const(vf)+const(chf)+const(wmi),
#'                   data=TRACE,cause=9)
#' sim1 <- sim.cif(cif,100,data=TRACE)
#' #' cc <-  crr 
#' cc <-  comp.risk(Event(time,status)~const(vf)+const(chf)+const(wmi),
#'                   data=sim1,cause=1)
#' cbind(cif$gamma,cc$gamma)
#' plot(cif) 
#' lines(cc$cum)
#' 
#' # faster/better with  mets package 
#' # library(mets)
#' # scif <-  cifreg(Event(time,status)~vf+chf+wmi,data=sim1,cause=1,prop=NULL)
#' #   
#' # plot(scif$cum,type="l")
#' # lines(cif$cum,col=2)
#' # cbind(cif$gamma,scif$coef)
#' # 
#' ################################################################
#' #  simulating several causes with specific cumulatives 
#' ################################################################
#' data(bmt)
#' cif1 <-  comp.risk(Event(time,cause)~const(tcell)+const(age),
#'                   data=bmt,cause=1,model="logistic2")
#' cif2 <-  comp.risk(Event(time,cause)~const(tcell)+const(age),
#'                   data=bmt,cause=2,model="logistic2")
#' 
#' ## must look at same time-scale
#' cifs <- pre.cifs(list(cif1,cif2))
#' plot(cifs[[1]]$cum,type="l")
#' lines(cifs[[2]]$cum,col=2)
#' legend("topleft",c("cause1","cause2"),lty=1,col=1:2)
#' 
#'  n <- 100
#'  sim1 <- sim.cif(cifs[[1]],n,data=bmt)
#'  Z <- sim1[,c("tcell","age")]
#'  sim2 <- sim.cif(cifs[[2]],n,data=bmt,Z=Z,drawZ=FALSE)
#'  ###
#'  rt <- rbinom(n,1,(sim1$F1tau+sim2$F1tau))
#'  rb <- rbinom(n,1,sim1$F1tau/(sim1$F1tau+sim2$F1tau))
#'  cause=ifelse(rb==1,1,2)
#'  time=ifelse(cause==1,sim1$timecause,sim2$timecause)
#'  cause <- rt*cause
#'  time[cause==0] <- tail(cifs[[1]]$cum[,1],1)
#' 
#'  bt <- data.frame(time=time,cause=cause,tcell=sim1$tcell,age=sim1$age)
#'  scif1 <-  comp.risk(Event(time,cause)~const(tcell)+const(age),
#'                    data=bt,cause=1,model="logistic2")
#'  scif2 <-  comp.risk(Event(time,cause)~const(tcell)+const(age),
#'                    data=bt,cause=2,model="logistic2")
#' 
#'  plot(scif1$cum,type="l")
#'  lines(scif2$cum,col=1,lty=2)
#'  legend("topleft",c("cause1","cause2"),lty=1:2,col=1:1)
#'  lines(cifs[[1]]$cum,col=2)
#'  lines(cifs[[2]]$cum,col=2,lty=2)
#'
#' #  Everyhing wrapped in call assuming covariates work in the same way for two models
#' dd <- sim.cifs(list(cif1,cif2),100,data=bmt)
#' scif1 <-  comp.risk(Event(time,cause)~const(tcell)+const(age),
#'                   data=dd,cause=1,model="logistic2")
#' scif2 <-  comp.risk(Event(time,cause)~const(tcell)+const(age),
#'                   data=dd,cause=2,model="logistic2")
#'
#' plot(scif1$cum,type="l")
#' lines(scif2$cum,col=1,lty=2)
#' legend("topleft",c("cause1","cause2"),lty=1:2,col=1:1)
#' lines(cifs[[1]]$cum,col=2)
#' lines(cifs[[2]]$cum,col=2,lty=2)
#'
#' # Everyhing wrapped in call assuming covariates work in the same way for two models
#' # but now draws cif1 to be of correct model, but model 2 is adapted 
#' #(if needed) to make constraints satisfied F1+F2 <=1
#' # see doubleFG of mets package for paramtrization
#' # and drawns as "if not cause1" then distribute according to cause 2
#' # dd <- sim.cifsRestrict(list(cif1,cif2),200,data=bmt)
#'
#' # faster with mets package 
#' # dd <- sim.cifs(list(cif1,cif2),1000,data=bmt)
#' # scif1 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=1)
#' # scif2 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=2)
#' #   
#' # plot(scif1$cum,type="l")
#' # legend("topleft",c("cause1","cause2"),lty=1:2,col=1:1)
#' # lines(cifs[[1]]$cum,col=2)
#' # lines(cifs[[2]]$cum,col=2,lty=2)
#' # 
#' @export
#' @aliases sim.cif sim.cifs subdist pre.cifs sim.cifsRestrict setup.cif
sim.cif <- function(cif,n,data=NULL,Z=NULL,drawZ=TRUE,cens=NULL,rrc=NULL,cumstart=c(0,0),...)
{# {{{
## also extracts coefficients and baseline from coxph, cox.aalen, phreg
## and uses them assuming that the model is cloglog unless otherwise
if (!inherits(cif,"defined"))  {
des <- 	read.fit(cif,n,data=data,Z=Z,drawZ=drawZ,...)
Z <- des$Z; cumhaz <- des$cum; rr <- des$rr; 
id <- des$id
cumhaz <- rbind(cumstart,cumhaz)
znames <- names(Z); 
model <- des$model
} else {
cumhaz <- cif$cumhaz
if (is.null(Z)) stop("must give Z")
rr <- exp(as.matrix(Z) %*% cif$coef)
if (is.data.frame(Z)) znames <- names(Z) else znames <- colnames(Z)
model <- cif$model
id <- 1:nrow(Z)
}

if (!inherits(cif,"phreg")) {
if (model=="logistic2" | model=="logistic") ptt <- simsubdist(cumhaz,rr,type="logistic",...) else ptt <- simsubdist(cumhaz,rr,type="cloglog",...)
    ptt <- cbind(ptt,Z)
} else { ### phreg class# {{{
	ptt <- data.frame()
	if (cif$nstrata>1) {
		stratj <- cif$strata[cif$jumps]
		for (j in 0:(cif$nstrata-1)) {
                        cumhazardj <- rbind(c(0,0),cif$cumhaz[stratj==j,])
			if (model[3]) 
		          pttj <- simsubdist(cumhazardj,rr[des$strataid==j],type="cloglog") else 
	                  pttj <- simsubdist(cumhazardj,rr[des$strataid==j],type="logistic") 
			Zj <- Z[des$strataid==j,,drop=FALSE]
			pttj <- cbind(pttj,Zj)
			ptt  <-  rbind(ptt,pttj)
			ptt  <-  rbind(ptt,pttj)
		}
	} else {
	if (model[3]) ptt <- simsubdist(cumhaz,rr,type="cloglog") else 
	              ptt <- simsubdist(cumhaz,rr,type="logistic") 
		ptt <- cbind(ptt,Z)
	}
 }# }}}

  ### adds censoring 
   if (!is.null(cens))  {# {{{
      if (is.null(rrc)) rrc <- rep(1,n)
      if (is.matrix(cens)) {
	   pct <- pc.hazard(cens,rrc,cum.hazard=TRUE,...)
	   pct <- pct$time
      }
      else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
           }
      }
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
   } # }}}

   attr(ptt,"model") <- model
   attr(ptt,"id") <-  id
   attr(ptt,"znames") <- znames
   return(ptt)
}# }}}

#' @export
sim.cifs <- function(cifs,n,data=NULL,Z=NULL,cens=NULL,rrc=NULL,max.times=NULL,causes=c(1,2),...)
{# {{{

if (!is.list(cifs)) stop("Cif models in list form\n"); 

  ## must consider all models out to last observation times or max.times
  cifs <- pre.cifs(cifs,max.times=max.times)

  tau <- tail(cifs[[1]]$cum[,1],1)
  sim1 <- sim.cif(cifs[[1]],n,data=data,Z=Z)
  Z <- sim1[,attr(sim1,"znames")]
  if (!inherits(cifs[[2]],"defined"))  {
      sim2p <- read.fit(cifs[[2]],1,data=data,Z=NULL)
      Z2 <- data[attr(sim1,"id"),names(sim2p$Z)]
      sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2,drawZ=FALSE)
  } else { 
        Z2 <- Z[,attr(cifs[[2]],"znames")]
	sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2)
  }

  ptot <- sim1$F1tau+sim2$F1tau
  ###
  rt <- rbinom(n,1,pmin(ptot,1))
  rb <- rbinom(n,1,sim1$F1tau/ptot)
  cause=ifelse(rb==1,1,2)
  time=ifelse(cause==causes[1],sim1$timecause,sim2$timecause)
  cause <- rt*cause
  time[cause==0] <- tau

  ptt <- data.frame(time=time,status=cause,cause=cause,ptot=ptot)
  Zcovs <- cbind(Z,Z2)
  samecovs <-  match(names(Z2),names(Z)) 
  Ze <- Zcovs[,-samecovs]
  ptt <- cbind(ptt,Ze)

   if (!is.null(cens))  {# {{{
      if (is.null(rrc)) rrc <- rep(1,n)
      if (is.matrix(cens)) {
	   pct <- pc.hazard(cens,rrc,cum.hazard=TRUE,...)
	   pct <- pct$time
      }
      else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
           }
      }
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
   } # }}}

   return(ptt)
}# }}}

#' @export
sim.cifsRestrict <- function(cifs,n,data=NULL,Z=NULL,cens=NULL,rrc=NULL,max.times=NULL,causes=c(1,2),...)
{# {{{

if (!is.list(cifs)) stop("Cif models in list form\n"); 

  ## must consider all models out to last observation times or max.times
  cifs <- pre.cifs(cifs,max.times=max.times)

  tau <- tail(cifs[[1]]$cum[,1],1)
  sim1 <- sim.cif(cifs[[1]],n,data=data,Z=Z)
  Z <- sim1[,attr(sim1,"znames")]
  if (!inherits(cifs[[2]],"defined"))  {
      sim2p <- read.fit(cifs[[2]],1,data=data,Z=NULL)
      Z2 <- data[attr(sim1,"id"),names(sim2p$Z)]
      sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2,drawZ=FALSE)
  } else { 
        Z2 <- Z[,attr(cifs[[2]],"znames")]
	sim2 <- sim.cif(cifs[[2]],n,data=data,Z=Z2)
  }


  ptot <- sim1$F1tau+sim2$F1tau*(1-sim1$F1tau)
  ###
  rt <- rbinom(n,1,pmin(ptot,1))
  rb <- rbinom(n,1,sim1$F1tau/ptot)
  cause=ifelse(rb==1,1,2)
  time=ifelse(cause==causes[1],sim1$timecause,sim2$timecause)
  cause <- rt*cause
  time[cause==0] <- tau

  ptt <- data.frame(time=time,status=cause,cause=cause,ptot=ptot)
  ptt <- cbind(ptt,Z)
  samecovs <-  match(names(Z2),names(Z)) 
  Ze <- Z2[,-samecovs]
  ptt <- cbind(ptt,Ze)

   if (!is.null(cens))  {# {{{
      if (is.null(rrc)) rrc <- rep(1,n)
      if (is.matrix(cens)) {
           cum.hazard <- TRUE
	   pct <- pc.hazard(cens,rrc,cum.hazard=TRUE,...)
	   pct <- pct$time
      }
      else {
	   if (is.numeric(cens)) pct<- rexp(n)/cens  else {
	      chaz <-sum(ptt$status)/sum(ptt$time)  ## hazard averate T haz 
	      pct<- rexp(n)/chaz 
           }
      }
      ptt$time <- pmin(ptt$time,pct)
      ptt$status <- ifelse(ptt$time<pct,ptt$status,0)
   } # }}}

   return(ptt)
}# }}}

#' @export
pre.cifs <- function(cifs,n,data=NULL,pprint=FALSE,max.times=NULL,...)
{# {{{

if (!is.list(cifs)) stop("Cif models in list form\n"); 

## iff needed put cumulative baseline in argument cifs$cum
## extracts maximum time of all cumulatives
maxtimes <- rep(0,length(cifs))
  for (i in 1:length(cifs)) { 
     if (inherits(cifs[[i]],"coxph")) {
	     stop("Use phreg of mets instead\n"); 
     } else  if (inherits(cifs[[i]],"crr")) {
	cifs[[i]]$cum <- rbind(c(0,0),cbind(cifs[[i]]$uftime,cumsum(cifs[[i]]$bfitj)))
        maxtimes[i] <- max(cifs[[i]]$uftime)
     } else if (inherits(cifs[[i]],"phreg")) {
	cifs[[i]]$cumhaz <- cifs[[i]]$cumhaz
        maxtimes[i] <- max(cifs[[i]]$cumhaz[,1])
     }  else if (inherits(cifs[[i]],"defined")) {
	cifs[[i]]$cumhaz <- cifs[[i]]$cumhaz
        maxtimes[i] <- max(cifs[[i]]$cumhaz[,1])
     } else {
        maxtimes[i] <- max(cifs[[i]]$cum[,1])
     }
  }

  mmtimes <- min(maxtimes)
  if (is.null(max.times)) mtimes <- mmtimes else mtimes <- max.times
  if (mtimes > mmtimes)  {
	  warning("max.time is further out than max for some cumulatives\n"); 
	  cat(" Max times for cumulatives in cifs \n"); 
	  print(maxtimes) 
  }

  for (i in 1:length(cifs)) { 
	cums <- cifs[[i]]$cum
        keep <- cums[,1]<mtimes
        Fmm <- subdist(cums,mtimes)
	cums <- cums[keep,,drop=FALSE]
	cums <- rbind(cums,Fmm)
	cifs[[i]]$cum <- cums
	cifs[[i]]$cumhaz <- cums
	if (pprint) {
	print(head(cums))
	print(tail(cums))
	}
  }

  return(cifs)
}# }}}


## reads a coxph, cox.aalen, phreg, crr, comprisk, prop.odds.subdist object
## and returns cumulative hazard, linear predictor of a resample of size
## n of data, for coxph cox.aalen, comprisk the design Z is constructed
## using the data, but Z can also be given which is needed for crr
## if drawZ is true the covariates in Z are used but multplied coefficient 
## based on column names, if fixZ=TRUE the matrix Z is multiplied coefficients 
## as is Z %*% coef(x)  to form linear predictor
read.fit <- function(cox,n,data=NULL,Z=NULL,drawZ=TRUE,fixZ=FALSE,id=NULL)
{# {{{

if (inherits(cox,"coxph"))
{# {{{
   ## basehaz(cox,centered=FALSE)[,c(2,1)] 
   sfit <- survfit(cox, se.fit=FALSE)   
   zcoef <- ifelse(is.na(coef(cox)), 0, coef(cox))
   offset <- sum(cox$means * zcoef)
   chaz <- sfit$cumhaz * exp(-offset)
   base <- cbind(sfit$time,chaz)  
   mt <- model.frame(cox)
   Y <- model.extract(mt, "response")
   if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
   if (attr(Y, "type") == "right") {
    time <- Y[, "time"]; 
    status <- Y[, "status"]
   } else stop("Expected right-censored data.");
  if (is.null(Z)) Z <- na.omit(model.matrix(cox))
  nn <- colnames(Z)
  if (fixZ) Z <- Z else Z <- Z[,names(coef(cox)),drop=FALSE] 
  rownames(Z) <- NULL
  jtime <- sort(time[status==1])
  cumhazard <- rbind(c(0,0),Cpred(base,jtime))
  rr <- exp( Z %*% matrix(coef(cox),ncol=1))
  if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
  if (!is.null(id)) xid <- id
  rr <- rr[xid]
  Z <- Z[xid,,drop=FALSE]
  colnames(Z) <- nn
  model <- "fg"
}# }}}
if (inherits(cox,"cox.aalen"))
{# {{{
   formula <- attr(cox, "Formula")
###   Z1 <- na.omit(model.matrix(cox,data=data))
   p <- nrow(cox$gamma)
   if (is.null(Z)) {
   Z <- na.omit(get_all_vars(formula,data=data))
   nz <- ncol(Z)
   Z <- Z[,seq(nz-p+1,nz)]
   }
   lrr <- as.matrix(Z) %*% cox$gamma
   cumhazard <- cox$cum
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   model <- "fg"
   if (cox$prop.odds==TRUE) model <- "logistic"
}# }}}
if (inherits(cox,"phreg"))
{# {{{
   p <- length(cox$coef)
   if (is.null(Z)) {
      Z <- cox$model.frame[,-1,drop=FALSE]
      if (cox$nstrata>1) {
	   ms <- match(cox$strata.name,names(Z))
           stratname <-  substring(cox$strata.name,8,nchar(cox$strata.name)-1)
	   strata <- cox$strata
	   Z  <-  Z[,-ms,drop=FALSE]
      }
   }
   nz <- ncol(Z)
   if (fixZ) Z <- Z else Z <- Z[,names(cox$coef),drop=FALSE] 
   lrr <- as.matrix(Z) %*% cox$coef
   cumhazard <- rbind(c(0,0),cox$cumhaz)
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   if (cox$nstrata>1) stratid <- strata[xid] else stratid <- NULL
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   if (cox$nstrata>1) {
	   Z[,stratname] <- stratid
   }
   model <-c(class(cox),is.null(cox$propodds))
   ## if (cox$prop.odds==TRUE) cox$model <- "logistic"
}# }}}
if (inherits(cox,"crr"))
{# {{{
   p <- length(cox$coef)
   if (is.null(Z)) stop("must give covariates for crr simulations\n");
   nz <- ncol(Z)
   rownames(Z) <- NULL
   if (fixZ) Z <- Z else Z <- Z[,names(cox$coef),drop=FALSE] 
   lrr <- as.matrix(Z) %*% cox$coef
   cumhazard <- rbind(c(0,0),cbind(cox$uftime,cumsum(cox$bfitj)))
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   model <- "fg"
}# }}}
if (inherits(cox,"comprisk"))
{# {{{
   p <- length(cox$gamma)
   formula <- attr(cox, "Formula")
###   Z1 <- na.omit(model.matrix(cox,data=data))
   p <- nrow(cox$gamma)
   if (is.null(Z)) {
	   Z <- na.omit(get_all_vars(formula,data=data))
	   nz <- ncol(Z)
	   Z <- Z[,seq(nz-p+1,nz),drop=FALSE]
   }
   nz <- ncol(Z)
###   if (fixZ) Z <- Z else Z <- Z[,names(cox$gamma),drop=FALSE] 
   lrr <- as.matrix(Z) %*% cox$gamma
###   lrr <- as.matrix(Z) %*% cox$gamma
   cumhazard <- cox$cum
   rr <- exp(lrr)
   if (drawZ==TRUE) xid <- sample(1:nrow(Z),n,replace=TRUE) else xid <- 1:nrow(Z)
   if (!is.null(id)) xid <- id
   rr <- rr[xid]
   Z <- Z[xid,,drop=FALSE]
   if (!(cox$model %in% c("fg","logistic2"))) stop("musg be fg or logistic2"); 
   model <- cox$model
}# }}}

out <- list(Z=Z,cum=cumhazard,rr=rr,id=xid,model=model)
if (inherits(cox,"phreg"))
	if (cox$nstrata>1) out <- c(out,list(strataid=stratid,strataname=stratname))
return(out)

}# }}}

