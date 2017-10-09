##' GOF for Cox PH regression
##'
##' Cumulative score process residuals for Cox PH regression
##' p-values based on Lin, Wei, Ying resampling.
##' @param x is phreg object 
##' @param n.sim number of simulations for score processes
##' @param silent to show timing estimate will be produced for longer jobs
##' @param ... Additional arguments to lower level funtions
##' @author THomas Scheike and Klaus K. Holst
##' @export
##' @aliases gof.phreg 
##' @examples
##' data(TRACE)
##' 
##' m1 <- phreg(Surv(time,status==9)~vf+chf+diabetes,data=TRACE) 
##' gg <- gof(m1)
##' par(mfrow=c(1,3))
##' plot(gg)
##' 
##' m1 <- phreg(Surv(time,status==9)~strata(vf)+chf+diabetes,data=TRACE) 
##' gg <- gof(m1)
##' @export
gof.phreg  <- function(x,n.sim=1000,silent=0,...)
{# {{{

### test for proportionality 
p <- length(x$coef)
nnames <- names(x$coef)
ii <- solve(x$hessian)
jumptimes <- x$jumptimes
Pt <- x$hessianttime
U <-  x$U
Pt <- apply(Pt,2,cumsum)
Ut <- apply(U,2,cumsum)

nd <- nrow(x$U)
Pt <-  .Call("CubeMat",Pt,ii,PACKAGE="mets")$XXX
sup <- matrix(0,n.sim,nrow(ii))
hatti <- matrix(0,nd,nrow(ii))
obs <- apply(abs(Ut),2,max)

tt <- system.time(simcox1<-.Call("PropTestCox",U,Pt,10,obs,PACKAGE="mets"))
prt <- n.sim*tt[3]/(10*60)
if (prt>1 & silent==0) cat(paste("Predicted time minutes",signif(prt,2),"\n"))

simcox <-  .Call("PropTestCox",U,Pt,n.sim,obs,PACKAGE="mets")
sup <-  simcox$supUsim
res <- cbind(obs,simcox$pval)
colnames(res) <- c("Sup|U(t)|","pval")
rownames(res) <- nnames 

if (silent==0) {
cat("Cumulative score process test for Proportionality:\n")
prmatrix(round(res,digits=2))
}

out <- list(jumptimes=x$jumptimes,supUsim=sup,res=res,supU=obs,
	    pvals=simcox$pval,score=Ut,simUt=simcox$simUt,type="prop")
class(out) <- "gof.phreg"
return(out)
}# }}}


##' GOF for Cox covariates in  PH regression
##'
##' Cumulative residuals after model matrix for Cox PH regression
##' p-values based on Lin, Wei, Ying resampling.
##' @param formula formula for cox regression 
##' @param data data for model
##' @param offset offset 
##' @param weights weights 
##' @param modelmatrix  matrix for cumulating residuals
##' @param n.sim number of simulations for score processes
##' @param silent to keep it absolutely silent, otherwise timing estimate will be prduced for longer jobs.
##' @param ... Additional arguments to lower level funtions
##' @author THomas Scheike and Klaus K. Holst
##' @export
##' @examples
##' data(TRACE)
##' 
##' dcut(TRACE)  <- ~. 
##' mm <- model.matrix(~-1+factor(wmicat.4),data=TRACE)
##' m1 <- gofM.phreg(Surv(time,status==9)~vf+chf+wmi,data=TRACE,modelmatrix=mm) 
##' summary(m1)
##' par(mfrow=c(2,2))
##' plot(m1)
##' 
##' m1 <- gofM.phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACE,modelmatrix=mm) 
##' summary(m1)
##' @export
gofM.phreg  <- function(formula,data,offset=NULL,weights=NULL,modelmatrix=NULL,
			n.sim=1000,silent=1,...)
{# {{{

if (is.null(modelmatrix)) stop(" must give matrix for cumulating residuals\n"); 

cox1 <- phreg(formula,data,offset=NULL,weights=NULL,Z=modelmatrix,cumhaz=FALSE,...) 
offsets <- as.matrix(cox1$model.frame[,names(cox1$coef)]) %*% cox1$coef
if (!is.null(offset)) offsets <- offsets*offset
if (!is.null(cox1$strata)) 
     coxM <- phreg(cox1$model.frame[,1]~modelmatrix+strata(cox1$strata),data,offset=offsets,weights=weights,no.opt=TRUE,cumhaz=FALSE,...)
else coxM <- phreg(cox1$model.frame[,1]~modelmatrix,data,offset=offsets,weights=weights,no.opt=TRUE,cumhaz=FALSE,...)
nnames <- colnames(modelmatrix)

Ut <- apply(coxM$U,2,cumsum)
jumptimes <- coxM$jumptimes
U <- coxM$U
Ubeta <- cox1$U
ii <- -solve(cox1$hessian)
EE <- .Call("vecMatMat",coxM$E,cox1$E)$vXZ; 
Pt <- cox1$ZX - EE
Pt <- apply(Pt,2,cumsum)
betaiid <- t(ii %*% t(Ubeta))
obs <- apply(abs(Ut),2,max)
simcox <-  .Call("ModelMatrixTestCox",U,Pt,betaiid,n.sim,obs,PACKAGE="mets")

sup <-  simcox$supUsim
res <- cbind(obs,simcox$pval)
colnames(res) <- c("Sup|U(t)|","pval")
rownames(res) <- nnames 

if (silent==0) {
   cat("Cumulative score process test for modelmatrix:\n")
   prmatrix(round(res,digits=2))
}

out <- list(jumptimes=jumptimes,supUsim=simcox$supUsim,res=res,supU=obs,
	    pvals=simcox$pval,score=Ut,simUt=simcox$simUt,type="modelmatrix")
class(out) <- "gof.phreg"

return(out)
}# }}}

##' Stratified baseline graphical GOF test for Cox covariates in PH regression
##'
##' Looks at stratified baseline in Cox model and plots all baselines versus each
##' other to see if lines are straight, with 50 resample versions under the 
##' assumptiosn that the stratified Cox is correct 
##'
##' @param x phreg object
##' @param sim to simulate som variation from cox model to put on graph
##' @param silent to keep it absolutely silent 
##' @param ... Additional arguments to lower level funtions
##' @author THomas Scheike and Klaus K. Holst
##' @export
##' @examples
##' data(TRACE)
##' 
##' m1 <- phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACE) 
##' m2 <- phreg(Surv(time,status==9)~vf+strata(chf)+wmi,data=TRACE) 
##' par(mfrow=c(2,2))
##' gofG.phreg(m1)
##' gofG.phreg(m2)
##' 
##' @export
gofG.phreg  <- function(x,sim=1,silent=1,...)
{# {{{

p <- length(x$coef)
nnames <- names(x$coef)
strata <- x$strata[x$jumps]
nstrata <- x$nstrata
jumptimes <- x$jumptimes
cumhaz <- x$cumhaz

ms <- match(x$strata.name,names(x$model.frame))
lstrata <- levels(x$model.frame[,ms])
stratn <-  substring(x$strata.name,8,nchar(x$strata.name)-1)
stratnames <- paste(stratn,lstrata,sep=":")

if (is.null(cumhaz)) stop("Must run phreg with cumhaz=TRUE (default)"); 
if (nstrata==1) stop("Stratified Cox to look at baselines");
###

for (i in 0:(nstrata-2))
for (j in (i+1):(nstrata-1)) { 
      iij <- which(strata %in% c(i,j))
      ii <- which(strata %in% i)
      ij <- which(strata %in% j)
      dijjumps  <- jumptimes[iij] 
      cumhazi <- Cpred(cumhaz[ii,],dijjumps,strict=FALSE)
      cumhazj <- Cpred(cumhaz[ij,],dijjumps,strict=FALSE)

      plot(cumhazj[,2],cumhazi[,2],type="s",lwd=2,
	   xlab=stratnames[j+1],ylab=stratnames[i+1])
      title(paste("Stratified baselines for",stratn))
      legend("topleft",c("Nonparametric","lm","Stratified-Cox-Sim"),lty=1,col=1:3)
      ab <- lm(cumhazi[,2]~-1+cumhazj[,2])
      if (sim==1) {
             Pt <- DLambeta.t <- apply(x$E/c(x$S0),2,cumsumstrata,strata,nstrata)
             II <- -solve(x$hessian)
             betaiid <- t(II %*% t(x$U))
	     simband <-  .Call("simBandCumHazCox",1/x$S0,Pt,betaiid,50,rep(1,nrow(Pt)),PACKAGE="mets")
	     simU <-simband$simUt
	     for (k in 1:50)
	     {
	      di <- Cpred(cbind(jumptimes[ii],simU[ii,k]),dijjumps,strict=FALSE)[,2]
	      dj <- Cpred(cbind(jumptimes[ij],simU[ij,k]),dijjumps,strict=FALSE)[,2]
	      lines(cumhazj[,2]+dj,cumhazi[,2]+di,type="s",lwd=0.1,col=3)
	     }
      }
      lines(cumhazj[,2],cumhazi[,2],type="s",lwd=2,col=1)
      abline(c(0,coef(ab)),col=2,lwd=2)
}

}# }}}

##' @export
plot.gof.phreg <-  function(x,col=3)
{# {{{
p <- ncol(x$score)
for (i in 1:p)
{
simU <- x$simUt[,(0:49)*p+i]
rsU <- max(abs(simU))
rsU <- max(rsU,abs(x$score[,i]))
plot(x$jumptimes,x$score[,i],type="s",ylim=c(-rsU,rsU),xlab="",ylab="")
title(main=rownames(x$res)[i])
matlines(x$jumptimes,simU,type="s",lwd=0.3,col=col)
lines(x$jumptimes,x$score[,i],lwd=1.5)
}

}# }}}

##' @export
summary.gof.phreg <-  function(x)
{# {{{
if (x$type=="prop")
     cat("Cumulative score process test for Proportionality:\n")
else cat("Cumulative residuals versus modelmatrix :\n")
print(x$res)
} # }}}

##' @export
print.gof.phreg <-  function(x)
{# {{{
if (x$type=="prop")
     cat("Cumulative score process test for Proportionality:\n")
else cat("Cumulative residuals versus modelmatrix :\n")
print(x$res)
} # }}}

