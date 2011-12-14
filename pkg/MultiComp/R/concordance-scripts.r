
conc2case <- function(conc,marg,test=1)
{ ## {{{
out <- conc
time1 <- conc$time
time2 <- marg$time

mintime <- max(time1[1],time2[1])
maxtime <- min(max(time1),max(time2))
timer <- seq(mintime, maxtime,length=100)
dtimer <- timer[2]-timer[1]

out$time <- timer
margtime <- Cpred(cbind(marg$time,c(marg$P1)),timer)[,2]
concP1 <- Cpred(cbind(conc$time,c(conc$P1)),timer)[,2]
out$P1 <- concP1/margtime

if (!is.null(out$se.P1)) out$se.P1 <-
	as.matrix(Cpred(cbind(conc$time,c(conc$se.P1)),timer)[,2]/margtime,ncol=100)
if (!is.null(out$P1.iid)) out$P1.iid <- Cpred(cbind(conc$time,conc$P1.iid[1,,]),timer)[,-1]/margtime
if (!is.null(out$P1.iid)) {
        margP1.iid  <- Cpred(cbind(marg$time,marg$P1.iid[1,,]),timer)[,-1]
	diff.iid <- out$P1.iid-margP1.iid
	sdiff.iid <- apply(diff.iid*dtimer,2,sum)
	sd.pepem <- sum(sdiff.iid^2)^.5
	diff.pepem <- sum((out$P1-margtime)*dtimer)
	z.pepem <- diff.pepem/sd.pepem
	pval.pepem <- 2*pnorm(-abs(z.pepem))
	cat("Pepe-Mori type test for H_0: casewise(t)= marginal-cif(t)\n")
		cat(paste("Time.range =",signif(mintime,3),"--",signif(maxtime,3),"\n\n"));
        outtest <- cbind(diff.pepem,sd.pepem,z.pepem,pval.pepem)
        colnames(outtest) <- c("int case1.1(t) - p_1(t) dt","sd","z","pval") 			 
	rownames(outtest) <- "pepe-mori"
	print(outtest,4)
###        prmatrix(outtest,3)			
} 

attr(out,"class") <- rev(attr(out,"class")) 
return(out)
} ## }}}

test.conc <- function(conc1,conc2,same.cluster=FALSE)
{ ## {{{
time  <- time1 <- conc1$time
time2 <- conc2$time
mintime <- max(time1[1],time2[1])
maxtime <- min(max(time1),max(time2))
timer <- seq(mintime, maxtime,length=100)
dtimer <- timer[2]-timer[1]

conc2timer <- Cpred(cbind(conc2$time,c(conc2$P1)),timer)[,2]
conc1timer <- Cpred(cbind(conc1$time,c(conc1$P1)),timer)[,2]

if (!is.null(conc1$P1.iid)) 
if (!is.null(conc2$P1.iid)) {
	### iid version af integraler
        conc2P1.iid  <- Cpred(cbind(conc2$time,conc2$P1.iid[1,,]),timer)[,-1]
        conc1P1.iid  <- Cpred(cbind(conc1$time,conc1$P1.iid[1,,]),timer)[,-1]
	if (same.cluster==TRUE) {
	diff.iid <- conc1P1.iid-conc2P1.iid
	sdiff.iid <- apply(diff.iid*dtimer,2,sum)
	sd.pepem <- sum(sdiff.iid^2)^.5
	} else {
	diff2.iid <- conc2P1.iid
	sdiff2.iid <- apply(diff2.iid*dtimer,2,sum)
	var2.pepem <- sum(sdiff2.iid^2)
	diff1.iid <- conc1P1.iid
	sdiff1.iid <- apply(diff1.iid*dtimer,2,sum)
	var1.pepem <- sum(sdiff1.iid^2)
	sd.pepem <- (var1.pepem+var2.pepem)^.5
	}
	diff.pepem <- sum(conc2timer-conc1timer)*dtimer
###	print(cbind(conc2timer,conc1timer))
	z.pepem <- diff.pepem/sd.pepem
	pval.pepem <- 2*pnorm(-abs(z.pepem))
	cat("\n")
	cat("Pepe-Mori type test for H_0: conc_1(t)= conc_2(t)\n")
	if (same.cluster==TRUE) cat("Assuming same clusters for the two functions\n") else 
		cat("Assuming independence for estimators\n");
		cat(paste("Time.range =",signif(mintime,3),"--",signif(maxtime,3),"\n\n"));
        outtest <- cbind(diff.pepem,sd.pepem,z.pepem,pval.pepem)
        colnames(outtest) <- c("int conc_1(t) - conc_2(t) dt","sd","z","pval") 			 
	rownames(outtest) <- "pepe-mori"
	print(outtest,4)
###        prmatrix(outtest,3)			
} 

###attr(out,"class") <- rev(attr(out,"class")) 
return(outtest)
} ## }}}

back2timereg <- function(obj)
{ ## {{{
	out <- obj
attr(out,"class") <- rev(attr(out,"class")) 
return(out)
} ## }}}

