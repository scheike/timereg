##' GOF for Cox PH regression
##'
##' Cumulative score process residuals for Cox PH regression
##' @param x is phreg object 
##' @param n.sim number of simulations for score processes
##' @param ... Additional arguments to lower level funtions
##' @author THomas Scheike and Klaus K. Holst
##' @export
##' @aliases gof.phreg 
##' @examples
##' data(TRACE)
##' 
##' m1 <- phreg(Surv(time,status==9)~vf+chf+diabetes,data=TRACE) 
##' 
##' gg <- gof.phreg(m1)
##' par(mfrow=c(1,3))
##' plot(gg)
##' 
##' m1 <- phreg(Surv(time,status==9)~strata(vf)+chf+diabetes,data=TRACE) 
##' gg <- gof.phreg(m1)
##' @export
gof.phreg  <- function(x,n.sim=1000)
{# {{{
### test for proportionality 
p <- length(x$coef)
nnames <- names(x$coef)
ii <- solve(x$hessian)
if (is.list(x$jumptimes)) jumptimes <- unlist(x$jumptimes) else jumptimes <- x$jumptimes
ord <- order(jumptimes)
Pt <- x$hessianttime[ord,]
U <-  x$U[ord,]
Pt <- apply(Pt,2,cumsum)
Ut <- apply(U,2,cumsum)

nd <- nrow(x$U)
Pt <-  .Call("CubeMat",Pt,ii,PACKAGE="mets")$XXX
hatt <- Pt
sup <- matrix(0,n.sim,nrow(ii))
hatti <- matrix(0,nd,nrow(ii))
obs <- apply(abs(Ut),2,max)

tt <- system.time(simcox1<-.Call("PropTestCox",U,Pt,10,obs,PACKAGE="mets"))
prt <- n.sim*tt[3]/(10*60)
if (prt>1) cat(paste("Predicted time minutes",signif(prt,2),"\n"))

simcox <-  .Call("PropTestCox",U,Pt,n.sim,obs,PACKAGE="mets")
sup <-  simcox$supUsim
res <- matrix(0,p,2)
res[,2] <- simcox$pvals
colnames(res) <- c("Sup|U(t)|","pval")

cat("Cumulative score process test for Proportionality:\n")
print(res)

out <- list(jumptimes=x$jumptimes,supUsim=sup,res=res,supU=obs,
	    pvals=pvals,score=Ut,simUt=simcox$simUt)
class(out) <- "gof.phreg"
return(out)
}# }}}

##' @export
plot.gof.phreg <-  function(x)
{# {{{
p <- ncol(x$score)
for (i in 1:p)
{
simU <- x$simUt[,(0:49)*p+i]
rsU <- max(abs(simU))
rsU <- max(rsU,abs(x$score[,i]))
plot(x$jumptimes,x$score[,i],type="l",ylim=c(-rsU,rsU),xlab="",ylab="")
title(main=rownames(x$res)[i])
matlines(x$jumptimes,simU,type="l",lwd=0.6,col=3)
lines(x$jumptimes,x$score[,i],lwd=1.5)
}

}# }}}

##' @export
summary.gof.phreg <-  function(x)
{# {{{
cat("Cumulative score process test for Proportionality:\n")
print(x$res)
} # }}}

