#' Estimates restricted residual mean for Cox or Aalen model
#' 
#' The restricted means are the \deqn{ \int_0^\tau S(t) dt } the standard
#' errors are computed using the i.i.d. decompositions from the cox.aalen (that
#' must be called with the argument "max.timpoint.sim=NULL") or aalen function.
#' 
#' must have computed iid decomposition of survival models for standard errors
#' to be computed. Note that competing risks models can be fitted but then the
#' interpretation is not clear.
#' 
#' @param out an "cox.aalen" with a Cox model or an "aalen" model.
#' @param x matrix with covariates for Cox model or additive hazards model
#' (aalen).
#' @param tau restricted residual mean.
#' @param iid if iid=1 then uses iid decomposition for estimation of standard
#' errors.
#' @return Returns an object. With the following arguments:
#' \item{mean}{restricted mean for different covariates.}
#' \item{var.mean}{variance matrix.} \item{se}{standard errors.}
#' \item{S0tau}{estimated survival functions on time-range [0,tau].}
#' \item{timetau}{vector of time arguments for S0tau.}
#' @author Thomas Scheike
#' @references D. M. Zucker, Restricted mean life with covariates: Modification
#' and extension of a useful survival analysis method, J. Amer. Statist. Assoc.
#' vol. 93 pp. 702-709, 1998.
#' 
#' Martinussen and Scheike, Dynamic Regression Models for Survival Data,
#' Springer (2006).
#' @keywords survival
#' @examples
#' 
#' \donttest{
#' ### this example runs slowly and is therefore donttest
#' data(sTRACE)
#' sTRACE$cage <- scale(sTRACE$age)
#' # Fits Cox model  and aalen model 
#' out<-cox.aalen(Surv(time,status>=1)~prop(sex)+prop(diabetes)+prop(chf)+
#' 	       prop(vf),data=sTRACE,max.timepoint.sim=NULL,resample.iid=1)
#' outa<-aalen(Surv(time,status>=1)~sex+diabetes+chf+vf,
#' data=sTRACE,resample.iid=1)
#' 
#' coxrm <- restricted.residual.mean(out,tau=7,
#'    x=rbind(c(0,0,0,0),c(0,0,1,0),c(0,0,1,1),c(0,0,0,1)),iid=1)
#' plot(coxrm)
#' summary(coxrm)
#' 
#' ### aalen model not optimal here 
#' aalenrm <- restricted.residual.mean(outa,tau=7,
#'    x=rbind(c(1,0,0,0,0),c(1,0,0,1,0),c(1,0,0,1,1),c(1,0,0,0,1)),iid=1)
#' with(aalenrm,matlines(timetau,S0tau,type="s",ylim=c(0,1)))
#' legend("bottomleft",c("baseline","+chf","+chf+vf","+vf"),col=1:4,lty=1)
#' summary(aalenrm)
#' 
#' mm <-cbind(coxrm$mean,coxrm$se,aalenrm$mean,aalenrm$se)
#' colnames(mm)<-c("cox-res-mean","se","aalen-res-mean","se")
#' rownames(mm)<-c("baseline","+chf","+chf+vf","+vf")
#' mm
#' }
#' 
#' @export
restricted.residual.mean <- function(out,x=0,tau=10,iid=0)
{ ## {{{ 
  if ((!inherits(out,c('cox.aalen',"aalen","survfit")))) 
      stop ("Must be output from cox.aalen or aalen function\n") 

if (inherits(out,"survfit")) { ## {{{ 
   fit.table  <-  as.matrix(summary(out, rmean=tau)$table)
if (ncol(fit.table)==1) fit.table <- t(fit.table)
   ee  <-  fit.table[,"*rmean"]                     
   se  <- fit.table[,"*se(rmean)"]        
   variid <- diag(se^2)
   S0t <- NULL
   timetau <- NULL
} ## }}}

if (inherits(out,"cox.aalen")) { ## {{{ 
time <- out$cum[,1]
cumhaz <- out$cum[,2]
beta <- out$gamma

if (is.matrix(x)!=TRUE) x <- matrix(x,nrow=1)
timetau <- c(time[time<tau],tau)
RR <- c( exp( x %*% beta) ) 
cumhaz<- matrix(1,nrow=nrow(x)) %*% t(matrix(cumhaz,ncol=(ncol(out$cum)-1)))
S0 <- exp(-cumhaz*RR)
S0t <- Cpred(cbind(time,t(S0)),timetau,strict=FALSE)[,-1,drop=FALSE]
Lam0t <- Cpred(cbind(time,t(cumhaz)),timetau,strict=FALSE)[,-1,drop=FALSE]
ll <- length(timetau)
ee <- apply(diff(timetau)*S0t[-ll,,drop=FALSE],2,sum)
deet <- apply(Lam0t[-ll,]*diff(timetau)*S0t[-ll,,drop=FALSE],2,sum)
RRDeet <- RR* deet

nn <- nrow(out$gamma.iid)
etiid <- c()
if (iid==1) {
   for (j in 1:nn) { 
      gamiid <- x %*% out$gamma.iid[j,] 
      gamiid <- RRDeet *  gamiid 
      Biidj<-Cpred(cbind(time,out$B.iid[[j]]),timetau)[,2]
      baseiid <- RR* 
      apply(diff(timetau)*S0t[-ll,,drop=FALSE]*Biidj[-ll],2,sum)
	   etiid <- rbind(etiid,c(gamiid+ baseiid))
   }
   variid <- t(etiid) %*% etiid 
   se  <- diag(variid)^.5
} else { variid <- se <- NULL }
} ## }}}

if (inherits(out,"aalen"))  ## {{{ 
{
time <- out$cum[,1]
cumhaz <- out$cum[,-1]

timetau <- c(time[time<tau],tau)
cumhaz<- as.matrix(x) %*% t(cumhaz)
S0 <- exp(-cumhaz)
S0t <- Cpred(cbind(time,t(S0)),timetau,strict=FALSE)[,-1,drop=FALSE]
Lam0t <- Cpred(cbind(time,t(cumhaz)),timetau)[,-1,drop=FALSE]
ll <- length(timetau)
ee <- apply(diff(timetau)*S0t[-ll,,drop=FALSE],2,sum)
deet <- apply(Lam0t[-ll,]*diff(timetau)*S0t[-ll,,drop=FALSE],2,sum)
RRDeet <- deet

nn <- length(out$B.iid)
etiid <- c()
if (iid==1) {
   for (j in 1:nn) { 
###	   gamiid <- x %*% out$gamma.iid[j,] 
###	   gamiid <- RRDeet *  gamiid 
	   Biidj<-Cpred(cbind(time,out$B.iid[[j]]),timetau)[,-1]
	   Biidj<- t(x %*% t(Biidj))
	   baseiid <- apply(diff(timetau)*S0t[-ll,,drop=FALSE]*Biidj[-ll,],2,sum)
	   etiid <- rbind(etiid,c(baseiid))
   }
   variid <- t(etiid) %*% etiid 
   se  <- diag(variid)^.5
} else { variid <- se <- NULL }
} ## }}}

out <- list(mean=ee,var.mean=variid,se=se,S0tau=S0t,timetau=timetau)
class(out) <- "restricted.residual.mean"
return(out)
} ## }}} 

#' @export
plot.restricted.residual.mean <- function(x,...)
{ ## {{{ 
matplot(x$timetau,x$S0tau,type="s")
} ## }}} 

###print.restricted.residual.mean <- function(object,digits=3)
###{ ## {{{ 
###out <- cbind(object$mean,object$se)
###colnames(out) <- c("mean","se")
###
###prmatrix(signif(out,digits))
###cat("\n"); 
###} ## }}} 

#' @export
summary.restricted.residual.mean <- function(object, digits=3,...)
{ ## {{{ 
if (!is.null(object$se)) {
out <- cbind(object$mean,object$se)
colnames(out) <- c("mean","se") } else {
out <- matrix(object$mean,ncol=1)
colnames(out) <- "mean"} 

prmatrix(signif(out,digits))
cat("\n"); 
} ## }}} 

