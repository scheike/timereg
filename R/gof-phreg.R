##' GOF for Cox PH regression
##'
##' Cumulative score process residuals for Cox PH regression
##' p-values based on Lin, Wei, Ying resampling.
##' @param object is phreg object 
##' @param n.sim number of simulations for score processes
##' @param silent to show timing estimate will be produced for longer jobs
##' @param robust to control wether robust dM_i(t) or dN_i  are used for simulations
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
##' ## to get Martingale ~ dN based simulations
##' gg <- gof(m1)
##' 
##' ## to get Martingale robust simulations, specify cluster in  call 
##' m1 <- phreg(Surv(time,status==9)~chf+diabetes+cluster(id),data=TRACE) 
##' gg <- gof(m1)
##' 
##' @export
gof.phreg  <- function(object,n.sim=1000,silent=1,robust=NULL,...)
{# {{{

### test for proportionality 
p <- length(object$coef)
nnames <- names(object$coef)
ii <- solve(object$hessian)
jumptimes <- object$jumptimes
Pt <- object$hessianttime
U <-  object$U
Pt <- apply(Pt,2,cumsum)
Ut <- apply(U,2,cumsum)

nd <- nrow(object$U)
Pt <-  .Call("CubeMat",Pt,ii,PACKAGE="mets")$XXX
sup <- matrix(0,n.sim,nrow(ii))
hatti <- matrix(0,nd,nrow(ii))
obs <- apply(abs(Ut),2,max)

if (is.null(robust)) 
     if (!is.null(object$id)) robust <- TRUE else robust <- FALSE

### cluster call or robust \hat M_i(t) based  
if (robust) {
        xx <- object$cox.prep
        S0i <- rep(0,length(xx$strata))
        S0i[xx$jumps+ 1] <- 1/object$S0
        Z <- xx$X
        ZdN <- U <- E <- matrix(0, nrow(xx$X), object$p)
        E[xx$jumps + 1, ] <- object$E
        U[xx$jumps + 1, ] <- object$U
        cumhaz <- c(cumsumstrata(S0i, xx$strata, xx$nstrata))
        EdLam0 <- apply(E * S0i, 2, cumsumstrata, xx$strata, xx$nstrata)
        rr <- c(xx$sign * exp(Z %*% coef(object) + xx$offset))
        MGt <- U[, drop = FALSE] - (Z * cumhaz - EdLam0) * rr * c(xx$weights)
        ### also weights 
        w <- c(xx$weights)
        nn <- nrow(Z)

	tt <- system.time(simcox1<- .Call("PropTestCoxClust",MGt,Pt,w*rr,Z,cumhaz,EdLam0,10,obs,nn,xx$id,xx$strata,xx$nstrata,xx$jumps))

	prt <- n.sim*tt[3]/(10*60)
	if (prt>1 & silent==0) cat(paste("Predicted time minutes",signif(prt,2),"\n"))

	simcox <- .Call("PropTestCoxClust",MGt,Pt,w*rr,Z,cumhaz,EdLam0,n.sim,obs,nn,xx$id,rep(0,nn),1,xx$jumps)
} else {
### no  or dN_i based  
	tt <- system.time(simcox1<-.Call("PropTestCox",U,Pt,10,obs,PACKAGE="mets"))
	prt <- n.sim*tt[3]/(10*60)
	if (prt>1 & silent==0) cat(paste("Predicted time minutes",signif(prt,2),"\n"))

	simcox <-  .Call("PropTestCox",U,Pt,n.sim,obs,PACKAGE="mets")
}
sup <-  simcox$supUsim
res <- cbind(obs,simcox$pval)
colnames(res) <- c("Sup|U(t)|","pval")
rownames(res) <- nnames 

if (silent==0) {
cat("Cumulative score process test for Proportionality:\n")
prmatrix(round(res,digits=2))
}

out <- list(jumptimes=object$jumptimes,supUsim=sup,res=res,supU=obs,
	    pvals=simcox$pval,score=Ut,simUt=simcox$simUt,type="prop",robust=robust)
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
##' library(mets)
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
##' 
##' ## cumulative sums in covariates, via design matrix mm 
##' mm <- cumContr(TRACE$wmi,breaks=10,equi=TRUE)
##' m1 <- gofM.phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACE,
##' 		  modelmatrix=mm,silent=0)
##' summary(m1)
##' 
##' ## cumulative sums in covariates, via design matrix mm 
##' m1 <- gofZ.phreg(Surv(time,status==9)~strata(vf)+chf+wmi+age,data=TRACE,
##'                  vars=c("wmi","age"))
##' summary(m1)
##' plot(m1,type="z")
##' 
##' @aliases cumContr gofZ.phreg
##' @export
gofM.phreg  <- function(formula,data,offset=NULL,weights=NULL,modelmatrix=NULL,
			n.sim=1000,silent=1,...)
{# {{{

if (is.null(modelmatrix)) stop(" must give matrix for cumulating residuals\n"); 

cox1 <- phreg(formula,data,offset=NULL,weights=NULL,Z=modelmatrix,cumhaz=FALSE,...) 
offsets <- cox1$X %*% cox1$coef
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
EE <- .Call("vecMatMat",coxM$E,cox1$E,PACKAGE="mets")$vXZ; 
###print(dim(EE))
###print(cox1$ZX)
Pt <- cox1$ZX - EE
Pt <- apply(Pt,2,cumsum)
betaiid <- t(ii %*% t(Ubeta))
obs <- apply(abs(Ut),2,max)
simcox <-  .Call("ModelMatrixTestCox",U,Pt,betaiid,n.sim,obs,PACKAGE="mets")

sup <-  simcox$supUsim
res <- cbind(obs,simcox$pval)
colnames(res) <- c("Sup_t |U(t)|","pval")
rownames(res) <- nnames 

if (silent==0) {
   cat("Cumulative score process test for modelmatrix:\n")
   prmatrix(round(res,digits=2))
}

 ## pvals efter z i model.matrix sup_z | M(z,tau) | 
 Utlast <- max(abs(tail(Ut,1)))
 maxlast <- apply(abs(simcox$last),1,max)
 pval.last <- mean(maxlast>=Utlast)
 res.last <- matrix(c(Utlast,pval.last),1,2)
 colnames(res.last) <- c("Sup_z |U(tau,z)|","pval")
 rownames(res.last) <- "matrixZ"

out <- list(jumptimes=jumptimes,supUsim=simcox$supUsim,res=res,supU=obs,
	    pvals=simcox$pval,score=Ut,simUt=simcox$simUt,
	    simUtlast=simcox$last,Utlast=Utlast,pval.last=pval.last,
	    res.last=res.last, type="modelmatrix")
class(out) <- "gof.phreg"

return(out)
}# }}}

##' @export
gofZ.phreg  <- function(formula,data,vars,offset=NULL,weights=NULL,breaks=20,equi=TRUE,
			n.sim=1000,silent=1,...)
{# {{{
 gres <- list()
 res <- matrix(0,length(vars),2)
 colnames(res) <- c("Sup_z |U(tau,z)|","pval")
 rownames(res) <- vars

i <- 1
for (vv in vars) {
 modelmatrix <- cumContr(data[,vv],breaks=breaks,equi=equi)
 lres <- gofM.phreg(formula,data,modelmatrix=modelmatrix) 

 res[i,] <- c(lres$Utlast,lres$pval.last)
 i <- i+1

 lres <- list(lres)
 names(lres) <- vv
 gres <- c(gres,lres)
}

out <- list(res=res,Zres=gres,type="Zmodelmatrix")
class(out) <- "gof.phreg"

return(out)
}# }}}

gofZ.phregGam  <- function(formula,data,vars,offset=NULL,weights=NULL,breaks=40,equi=TRUE,
			n.sim=1000,silent=1,...)
{# {{{

 res <- matrix(0,length(vars),2)
 colnames(res) <- c("Sup_z |U(tau,z)|","pval")
 rownames(res) <- vars

 i <- 1
for (vv in vars) {
 modelmatrix <- cumContr(data[,vv],breaks=breaks,equi=equi)

cox1 <- phreg(formula,data,offset=NULL,weights=NULL,Z=modelmatrix,cumhaz=FALSE,...) 
offsets <- cox1$X %*% cox1$coe
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
EE <- .Call("vecMatMat",coxM$E,cox1$E,PACKAGE="mets")$vXZ; 
Pt <- cox1$ZX - EE
Pt <- apply(Pt,2,cumsum)
betaiid <- t(ii %*% t(Ubeta))
obs <- apply(abs(Ut),2,max)
simcox <-  .Call("ModelMatrixTestCox",U,Pt,betaiid,n.sim,obs,PACKAGE="mets")

 ## pvals efter z i model.matrix sup_z | M(z,tau) | 
 Utlast <- max(abs(tail(Ut,1)))
 maxlast <- apply(abs(simcox$last),1,max)
 pval.last <- mean(maxlast>=Utlast)
 res[i,] <- c(Utlast,pval.last)
 i <- i+1
}

out <- list(res=res, type="modelmatrix")
class(out) <- "gof.phreg"

return(out)
}# }}}

##' @export
cumContr <- function(data, breaks = 4, probs = NULL,equi = TRUE,na.rm=TRUE,...)
 {# {{{
 if (is.vector(data)) {
        if (is.list(breaks))
            breaks <- unlist(breaks)
        if (length(breaks) == 1) {
            if (!is.null(probs)) {
                breaks <- quantile(data, probs, na.rm = na.rm,...)
	        breaks <- breaks[-1]
            }
            else {
                if (!equi) {
                  probs <- seq(0, 1, length.out = breaks + 1)
                  breaks <- quantile(data, probs, na.rm = na.rm, ...)
	          breaks <- breaks[-1]
                }
                if (equi) {
                  rr <- range(data, na.rm = na.rm)
                  breaks <- seq(rr[1], rr[2], length.out = breaks + 1)
	          breaks <- breaks[-1]
                }
            }
        }
        if (sum(duplicated(breaks)) == 0) {
		i <- 0; 
		gm <- matrix(0,length(data),length(breaks))
            for (bb in breaks)  {
		    i <- i+1 
		    gm[,i] <- (data <= bb)*1
	    }
	} else {
            wd <- which(duplicated(breaks))
            mb <- min(diff(breaks[-wd]))
            breaks[wd] <- breaks[wd] + (mb/2) * seq(length(wd))/length(wd)
            i <- 0; gm <- matrix(0,length(data),length(breaks))
            for (bb in breaks)  {
		    i <- i+1 
		    gm[,i] <- (data <= bb)*1
	    }
            warning(paste("breaks duplicated"))
        }
        colnames(gm) <- paste("<=",breaks,sep="")
	attr(gm,"breaks") <- breaks
        return(gm)
    }
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
##' @param lm addd line to plot, regressing the cumulatives on each other  
##' @param ... Additional arguments to lower level funtions
##' @author THomas Scheike and Klaus K. Holst
##' @export
##' @examples
##' data(TRACE)
##' 
##' m1 <- phreg(Surv(time,status==9)~strata(vf)+chf+wmi,data=TRACE) 
##' m2 <- phreg(Surv(time,status==9)~vf+strata(chf)+wmi,data=TRACE) 
##' par(mfrow=c(2,2))
##' 
##' gofG.phreg(m1)
##' gofG.phreg(m2)
##' 
##' bplot(m1,log="y")
##' bplot(m2,log="y")
##' @export
gofG.phreg  <- function(x,sim=0,silent=1,lm=TRUE,...)
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

if (is.null(x$opt) | is.null(x$coef)) fixbeta<- 1 else fixbeta <- 0

for (i in 0:(nstrata-2))
for (j in (i+1):(nstrata-1)) { 
      iij <- which(strata %in% c(i,j))
      ii <- which(strata %in% i)
      ij <- which(strata %in% j)
      dijjumps  <- jumptimes[iij] 
      cumhazi <- Cpred(cumhaz[ii,],dijjumps,strict=FALSE)
      cumhazj <- Cpred(cumhaz[ij,],dijjumps,strict=FALSE)

      plot(cumhazj[,2],cumhazi[,2],type="s",lwd=2,xlab=stratnames[j+1],ylab=stratnames[i+1])
      graphics::title(paste("Stratified baselines for",stratn))
      if ((fixbeta==0 | sim==0) & lm ) 
      graphics::legend("topleft",c("Nonparametric","lm"),lty=1,col=1:2)
###      graphics::legend("topleft",c("Nonparametric","Stratified-Cox-Sim"),lty=1,col=1:3)
      ab <- lm(cumhazi[,2]~-1+cumhazj[,2])
      if (sim==1 & fixbeta==0) {
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
      if (lm==TRUE) abline(c(0,coef(ab)),col=2,lwd=2)
}

}# }}}

##' @export
plot.gof.phreg <-  function(x,col=3,type="time",...)
{# {{{

if (type=="time") {
p <- ncol(x$score)
for (i in 1:p)
{
simU <- x$simUt[,(0:49)*p+i]
rsU <- max(abs(simU))
rsU <- max(rsU,abs(x$score[,i]))
plot(x$jumptimes,x$score[,i],type="s",ylim=c(-rsU,rsU),xlab="",ylab="")
title(main=rownames(x$res)[i])
matlines(x$jumptimes,simU,type="s",lwd=0.3,col=col)
lines(x$jumptimes,x$score[,i],type="s",lwd=1.5)
}
} else {
	if (x$type=="modelmatrix") {
		obsz <- c(tail(x$score,1))
		times <- 1:length(obsz)
		rsU <- max(max(abs(obsz)),max(abs(x$simUtlast[1:50,])))
		plot(times,obsz,type="l",ylim=c(-rsU,rsU),xlab="",ylab="")
		matlines(times,t(x$simUtlast[1:50,]),type="l",lwd=0.3,col=col)
	} else {
	   for (i in 1:length(x$Zres))
	   {
	    xr <- x$Zres[[i]]
	    obsz <- c(tail(xr$score,1))
	    times <- 1:length(obsz)
	    rsU <- max(max(abs(obsz)),max(abs(xr$simUtlast[1:50,])))
	    plot(times,obsz,type="l",ylim=c(-rsU,rsU),xlab="",ylab="")
	    matlines(times,t(xr$simUtlast[1:50,]),type="l",lwd=0.3,col=col)
           }
	}
}

}# }}}

##' @export
summary.gof.phreg <-  function(object,...)
{# {{{
if (object$type=="prop")
     cat("Cumulative score process test for Proportionality:\n")
else cat("Cumulative residuals versus modelmatrix :\n")
print(object$res)

if (!is.null(object$res.last)) {
   cat("\n")
   cat("Cumulative score process versus covariates (discrete z via model.matrix):\n")
   print(object$res.last)
}


} # }}}

##' @export
print.gof.phreg <-  function(x,...)
{# {{{
if (x$type=="prop")
     cat("Cumulative score process test for Proportionality:\n")
else cat("Cumulative residuals versus modelmatrix :\n")
print(x$res)

if (!is.null(x$res.last)) {
   cat("\n")
   cat("Cumulative score process versus covariates (discrete z via model.matrix):\n")
   print(x$res.last)
}

} # }}}

