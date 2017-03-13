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

sim.pc.hazard <- function(cumhazard,rr,cens=NULL,rrc=NULL,cens.cum.hazard=TRUE)
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

sim.cause.pc.hazard <- function(cumhaz1,cumhaz2,rr1,rr2,cens=NULL,rrc=NULL,cens.cum.hazard=TRUE)
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


