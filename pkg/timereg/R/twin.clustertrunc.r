twin.clustertrunc <- function(survformula,data=data,theta.des=NULL,clusters=NULL,
			       entry="v",exit="time",status="status")
{ ## {{{ 
pweight <- rep(1,nrow(data))
if (!is.null(clusters)) clusters <- data[,clusters] else stop("must give clusters\n"); 
if (is.null(theta.des)) ptheta <- 0 else ptheta <- rep(0,ncol(theta.des))

for (i in 1:10)
{ ## {{{ 
  assign("pweight",pweight,env=environment(survformula))
  aout <- cox.aalen(survformula,data=data,weights=1/pweight,robust=0,n.sim=0)
  tout <- two.stage(aout,data=data,clusters=clusters,theta.des=theta.des)
  if (!is.null(theta.des)) theta <- c(theta.des %*% tout$theta) else theta <- tout$theta
  if (attr(tout, "var.link") == 1) theta <- exp(tout$theta) 
###
  data$thetades  <- theta
  dd <- fast.reshape(data,id=clusters)
  v1 <- dd[,paste(entry,"1",sep="")]; v2 <- dd[,paste(entry,"2",sep="")]
  time1 <- dd[,paste(exit,"1",sep="")]; time2 <- dd[,paste(exit,"2",sep="")]
  status1 <- dd[,paste(status,"1",sep="")]; status2 <- dd[,paste(status,"2",sep="")]
  ptruncv1t2 <- predict.two.stage(tout,times=v1,times2=time2,theta=theta)
   ptrunct1v2 <- predict.two.stage(tout,times=time1,times2=v2,theta=theta)
  nn <- nrow(dd)
  ppv1t2 <- .Call("claytonoakesR",dd$thetades1,
      rep(0,nn),dd$status2,ptruncv1t2$S1, ptruncv1t2$S2)$like
  ppt1v2 <- .Call("claytonoakesR",dd$thetades1,
      dd$status1,rep(0,nn),ptrunct1v2$S1, ptrunct1v2$S2)$like
  ppt1v2[dd$status1==0] <- ppt1v2[status1==0]/ptrunct1v2$S1[status1==0]
  ppv1t2 <- .Call("claytonoakesR",dd$thetades1,
      dd$status1,rep(0,nn),ptrunct1v2$S1, ptrunct1v2$S2)$like
  ppv1t2[dd$status2==0] <- ppv1t2[status2==0]/ptruncv1t2$S2[status2==0]
###
  dd$weight1 <- c(ppt1v2)
  dd$weight2 <- c(ppv1t2)
  dd2 <- fast.reshape(dd)
  pweight <- dd2$weight
  dtheta <- tout$theta-ptheta
  if (sum(abs(dtheta)) < 0.001) break; 
  ptheta <- tout$theta
} ## }}} 

res <- list(marg=aout,two=tout,marg.weights=pweight,dtheta=dtheta)

return(res)
} ## }}} 

