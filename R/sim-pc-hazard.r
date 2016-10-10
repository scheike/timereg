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

