##' Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data 
##' using marginals that are on Cox or addtive form. 
##' If clusters contain more than two times, the algoritm uses a compososite likelihood
##' based on the pairwise bivariate models.
##'
##' The reported standard errors are based on the estimated information from the 
##' likelihood assuming that the marginals are known. 
##'
##'
##' @export
##' @references
##' Clayton-Oakes and Plackett bivariate survival distributions,
##' @examples
##' data(diabetes)
##'
##' data(diabetes)
##' # Marginal Cox model  with treat as covariate
##' margph <- coxph(Surv(time,status)~treat,data=diabetes)
##' ### Clayton-Oakes, from timereg 
##' fitco1<-two.stage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id)
##' summary(fitco1)
##' ### Plackett model 
##' fitp<-twostage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,var.link=1)
##' summary(fitp)
##' ### Clayton-Oakes
##' fitco2<-twostage(margph,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,var.link=1,model="clayton.oakes")
##' summary(fitco2)
##' 
##' ### without covariates using Aalen for marginals
##' marg <- aalen(Surv(time,status)~+1,data=diabetes,n.sim=0,max.clust=NULL,robust=0)
##' fitpa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,score.method="optimize")
##' summary(fitpa)
##' 
##' fitcoa<-twostage(marg,data=diabetes,theta=1.0,detail=0,Nit=40,clusters=diabetes$id,var.link=1,model="clayton.oakes")
##' summary(fitcoa)
##' 
##' diab1 <- surv.boxarea(c(10,10),c(40,40),diabetes) ### data truncated at (10,10) and right-censored (40,40)
##' diab2 <- surv.boxarea(c(40,40),c(90,90),diabetes) ### data truncated at (40,40) and right-censored (90,90)
##' 
##' ### dependence for early events
##' marg1 <- aalen(Surv(left,time,status)~+1,data=diab1,n.sim=0,max.clust=NULL,robust=0)
##' fit1<-twostage(marg1,data=diab1,clusters=diab1$id,model="clayton.oakes")
##' summary(fit1)
##' 
##' ### dependence for later events
##' marg2 <- aalen(Surv(left,time,status)~+1,data=diab2,n.sim=0,max.clust=NULL,robust=0)
##' fit2<-twostage(marg2,data=diab2,clusters=diab2$id,model="clayton.oakes")
##' summary(fit2)
##' 
##' fit1p<-twostage(marg1,data=diab1,clusters=diab1$id,score.method="optimize")
##' fit2p<-twostage(marg2,data=diab2,clusters=diab2$id,score.method="optimize")
##' summary(fit1p)
##' 
##' outc <- piecewise.twostage(c(0,15,40,90),data=diabetes,score.method="optimize",iid=1)
##' summary(outc)
##' outp <- piecewise.twostage(c(0,15,40,90),data=diabetes,score.method="optimize",iid=1)
##' summary(outp)
##' @keywords survival
##' @author Thomas Scheike
##' @export
twostage <- function(margsurv,data=sys.parent(),score.method="nlminb",
Nit=60,detail=0,clusters=NULL,silent=1,weights=NULL,
control=list(),theta=NULL,theta.des=NULL,var.link=1,iid=0,
step=0.5,notaylor=0,model="plackett")
{ ## {{{
## {{{ seting up design and variables
rate.sim <- 1; sym=1; 
if (model=="clayton.oakes") dep.model <- 1
else if (model=="plackett") dep.model <- 2
else stop("Model must by either clayton.oakes or plackett \n"); 

if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen") { ## {{{
	 formula<-attr(margsurv,"Formula");
	 beta.fixed <- attr(margsurv,"beta.fixed")
	 if (is.null(beta.fixed)) beta.fixed <- 1; 
	 ldata<-aalen.des(formula,data=data,model="cox.aalen");
	 id <- attr(margsurv,"id"); 
	 mclusters <- attr(margsurv,"cluster.call")
	 X<-ldata$X; time<-ldata$time2; Z<-ldata$Z;  status<-ldata$status;
	 time2 <- attr(margsurv,"stop"); start <- attr(margsurv,"start")
	 antpers<-nrow(X);
         if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
         if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
	 px<-ncol(X);

         if (is.null(clusters) && is.null(mclusters)) stop("No cluster variabel specified in marginal or twostage call\n"); 
         if (is.null(clusters)) clusters <- mclusters else if (sum(abs(clusters-mclusters))>0) 
         cat("Warning: Clusters for marginal model different than those specified for two.stage\n"); 
         if (!is.null(attr(margsurv,"max.clust")))
         if (attr(margsurv,"max.clust")< attr(margsurv,"orig.max.clust")) 
		  cat("Warning: Probably want to estimate marginal model with max.clust=NULL\n"); 
	 if (nrow(X)!=length(clusters)) stop("Length of Marginal survival data not consistent with cluster length\n"); 
## }}}
} else { ### coxph ## {{{
	  notaylor <- 1
	  antpers <- margsurv$n
	  id <- 0:(antpers-1)
	  mt <- model.frame(margsurv)
	  Y <- model.extract(mt, "response")
	  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
	  if (attr(Y, "type") == "right") {
	      time2 <- Y[, "time"]; 
	      status <- Y[, "status"]
		start <- rep(0,antpers);
		} else {
		 start <- Y[, 1]; 
		 time2 <- Y[, 2];
		 status <- Y[, 3];
		}
###	   Z <- na.omit(model.matrix(margsurv)[,-1]) ## Discard intercept
	   Z <- matrix(1,antpers,length(coef(margsurv)));

	   if (is.null(clusters)) stop("must give clusters for coxph\n");
	   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
	   px <- 1; pz <- ncol(Z); 
	   if (sum(abs(start))>0) lefttrunk <- 1 else lefttrunk <- 0
	   start <- rep(0,antpers);
	   beta.fixed <- 0
	   semi <- 1
###	   start.time <- 0
} ## }}}

  if (sum(abs(start))>0) lefttrunk <- 1  else lefttrunk <- 0;  
  RR <-  rep(1,antpers); 
  ptrunc <- rep(1,antpers);

  if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen")  { ## {{{
         resi <- residualsTimereg(margsurv,data=data) 
         RR  <- resi$RR
	 if (lefttrunk==1) ptrunc <- exp(-resi$cumhazleft); 
	 psurvmarg <- exp(-resi$cumhaz); 
###     if ((attr(margsurv,"residuals")!=2) || (lefttrunk==1)) { 
###         resi <- residualsTimereg(margsurv,data=data) 
###         residuals <- resi$residuals; 
######	 cumhaz <- resi$cumhaz; 
######	 cumhazleft <- resi$cumhazleft; 
###	 RR <- resi$RR
###	 ptrunc <- exp(-cumhazleft); 
###	 psurvmarg <- exp(-resi$cumhaz); 
###     } else { 
###	    if (is.null(margsurv$residuals$dM))
###	     residuals <- margsurv$residuals$dM; 
###             cumhaz <- status-residuals; 
###	      psurvmarg <- exp(-cumhaz); 
###	      if (class(margsurv)=="cox.aalen") RR  <- exp( Z %*% margsurv$gamma)
###     }
  } ## }}}
  else if (class(margsurv)=="coxph") {  ## {{{
       notaylor <- 1
       residuals <- residuals(margsurv)
       cumhaz <- status-residuals
       psurvmarg <- exp(-cumhaz); 
       cumhazleft <- rep(0,antpers)
       RR<- exp(margsurv$linear.predictors-margsurv$means*coef(margsurv))
        if ((lefttrunk==1)) { 
         baseout <- basehaz(margsurv,centered=FALSE); 
         cum <- cbind(baseout$time,baseout$hazard)
	 cum <- Cpred(cum,start)[,2]
	 ptrunc <- exp(-cum * RR)
	}
  } ## }}}
  if (is.null(weights)==TRUE) weights <- rep(1,antpers); 

  out.clust <- cluster.index(clusters);  
  clusters <- out.clust$clusters
  maxclust <- out.clust$maxclust 
  antclust <- out.clust$antclust
  clusterindex <- out.clust$idclust
  clustsize <- out.clust$cluster.size

###  update <- 2;
###  if (update==1) 
###  if (class(margsurv)=="aalen" || class(margsurv)=="cox.aalen")  { ## {{{
###     if ((attr(margsurv,"residuals")!=2) || (lefttrunk==1)) { 
###         resi <- residualsTimereg(margsurv,data=data) 
###         residuals <- resi$residuals; 
###	 cumhaz <- resi$cumhaz; 
###	 cumhazleft <- resi$cumhazleft; 
###	 RR <- resi$RR
###     } else { residuals <- margsurv$residuals$dM; 
###              cumhaz <- status-residuals; 
###	      if (class(margsurv)=="cox.aalen") RR  <- exp( Z %*% margsurv$gamma)
###     }
###  }
###  else if (class(margsurv)=="coxph") {
###       notaylor <- 1
###       residuals <- residuals(margsurv)
###       cumhaz <- status-residuals
###       psurvmarg <- exp(-cumhaz); ### survival
###       trunkp <- rep(1,length(status)); 
###       cumhazleft <- rep(0,antpers)
######    RR<- exp(Z %*% coef(margsurv))
###       RR<- exp(margsurv$linear.predictors-margsurv$means*coef(margsurv))
###        if ((lefttrunk==1)) { 
###         baseout <- basehaz(margsurv,centered=FALSE); 
###         cum <- cbind(baseout$time,baseout$hazard)
###	 cum <- Cpred(cum,start)[,2]
###	 cumhazleft <- cum * RR 
###	}
###  } ## }}}

  ratesim<-rate.sim; pxz <- px + pz;

###  Biid<-c(); gamma.iid <- 0; 
###  if (is.null(margsurv$B.iid)) notaylor <- 1; 
###  if (notaylor==0) {
###    if (!is.null(margsurv$B.iid))
###    for (i in 1:antclust) Biid<-cbind(Biid,margsurv$B.iid[[i]]); 
###    if (!is.null(margsurv$gamma.iid)) gamma.iid<-margsurv$gamma.iid;

  if (is.null(theta.des)==TRUE) ptheta<-1; 
  if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) else
  theta.des<-as.matrix(theta.des); 
  ptheta<-ncol(theta.des); 
  if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");

  if (is.null(theta)==TRUE && var.link==1) theta<-rep(0,ptheta); 
  if (is.null(theta)==TRUE && var.link==0) theta<-rep(1,ptheta); 
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

  if (maxclust==1) stop("No clusters, maxclust size=1\n"); 
  ## }}}

  par.func <- NULL
  loglike <- function(par) 
  { ## {{{

    if (is.null(par.func)==TRUE) {
       Xtheta <- theta.des %*% matrix(c(par),ptheta,1); 
       DXtheta <- array(0,c(1,1,1));
    } else { ## {{{
	  Xtheta <- c(); DXtheta <- array(0,c(length(times),antpers,dimpar));
	  if (is.null(dpar.func))
		  stop("Must also provide derivative of function wrt to parameters")
	  s <- 0
	  for (t in times) {
		  s <- 1+s
		  Xttheta <- par.func(par,t,theta.des)
		  if (length(Xttheta)!=antpers) stop("parfunc(par,t,theta.des) must length n"); 
		  Xtheta <- cbind(Xtheta,Xttheta)
		  Dttheta <- dpar.func(par,t,theta.des); 
		  if (dim(Dttheta)[1]!=antpers || dim(Dttheta)[2]!=dimpar) 
			  stop("dparfunc must return matrix n x dimpar when called on dparfunc(par,t,theta.des)"); 
		  DXtheta[s,,] <- Dttheta
	  }
    } ## }}}

###      dyn.load("twostage.so")
      outl<-.Call("twostageloglike", ## {{{
      icause=status,ipmargsurv=psurvmarg, 
      itheta=c(par),iXtheta=Xtheta,iDXtheta=DXtheta,idimDX=dim(DXtheta),ithetades=theta.des,
      icluster=clusters,iclustsize=clustsize,iclusterindex=clusterindex,
      ivarlink=var.link, iiid=iid,iweights=weights, isilent=silent, idepmodel=dep.model,
      itrunkp=ptrunc, DUP=FALSE) 
      ## }}}

      if (detail==3) print(c(par,outl$loglike))

    attr(outl,"gradient") <-outl$score 
    if (oout==0) ret <- c(-1*outl$loglike) else if (oout==1) ret <- sum(outl$score^2) else ret <- outl
    return(ret)
  } ## }}}

  if (score.method=="optimize" && ptheta!=1) score.method <- "nlminb"

  theta.iid <- NULL
  p <- theta
  if (score.method=="fisher.scoring") { ## {{{
    oout <- 2;  ### output control for obj
    if (Nit>0) 
    for (i in 1:Nit)
    {
        out <- loglike(p)
	hess <- out$Dscore
	if (!is.na(sum(hess))) hessi <- lava::Inverse(out$Dscore) else hessi <- hess 
        if (detail==1) {## {{{
          cat(paste("Fisher-Scoring ===================: it=",i,"\n")); 
          cat("theta:");print(c(p))
          cat("loglike:");cat(c(out$loglike),"\n"); 
          cat("score:");cat(c(out$score),"\n"); 
	  cat("hess:\n"); cat(out$Dscore,"\n"); 
        }## }}}
        delta <- hessi %*% out$score *step 
        p <- p+delta* step
	if (is.nan(sum(out$score))) break; 
        if (sum(abs(out$score))<0.00001) break; 
        if (max(theta)>20) break; 
    }
    if (!is.nan(sum(p))) { 
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    out <- loglike(p) 
    logl <- out$loglike
    score <- out$score
    oout <- 0; 
    score1 <- jacobian(loglike,p)
    hess <- hessian(loglike,p)
    if (iid==1) theta.iid <- out$theta.iid
    }
    if (detail==1 & Nit==0) {## {{{
          cat(paste("Fisher-Scoring ===================: final","\n")); 
          cat("theta:");print(c(p))
          cat("loglike:");cat(c(out$loglike),"\n"); 
          cat("score:");cat(c(out$score),"\n"); 
	  cat("hess:\n"); cat(out$Dscore,"\n"); 
    }## }}}
    if (!is.na(sum(hess))) hessi <- lava::Inverse(hess) else hessi <- diag(nrow(hess))
    ## }}}
  } else if (score.method=="nlminb") { ## {{{ nlminb optimizer
    oout <- 0; 
    tryCatch(opt <- nlminb(theta,loglike,control=control),error=function(x) NA)
    if (detail==1) print(opt); 
    library(numDeriv)
    score <- jacobian(loglike,opt$par)
    hess <- hessian(loglike,opt$par)
    hessi <- lava::Inverse(hess); 
    theta <- opt$par
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$par)
    logl <- out$loglike
    score1 <- out$score
    if (iid==1) theta.iid <- out$theta.iid
  ## }}}
  } else if (score.method=="optimize" && ptheta==1) { ## {{{  optimizer
    oout <- 0; 
    if (var.link==1) {mino <- -20; maxo <- 10;} else {mino <- 0.001; maxo <- 100;}
    tryCatch(opt <- optimize(loglike,c(mino,maxo)));
    if (detail==1) print(opt); 
    library(numDeriv)
    opt$par <- opt$minimum
    score <- jacobian(loglike,opt$par)
    hess <- hessian(loglike,opt$par)
    hessi <- lava::Inverse(hess); 
    theta <- opt$par
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$par)
    logl <- out$loglike
    score1 <- out$score
    if (iid==1) theta.iid <- out$theta.iid
  ## }}}
  } else if (score.method=="nlm") { ## {{{ nlm optimizer
    iid <- 0; oout <- 0; 
    tryCatch(opt <- nlm(loglike,theta,hessian=TRUE,print.level=detail),error=function(x) NA)
    iid <- 1; 
###    library(numDeriv)
    hess <- opt$hessian
    score <- opt$gradient
    if (detail==1) print(opt); 
    hessi <- lava::Inverse(hess); 
    theta <- opt$estimate
    if (detail==1 && iid==1) cat("iid decomposition\n"); 
    oout <- 2
    out <- loglike(opt$estimate)
    logl <- out$loglike
    score1 <- out$score
    if (iid==1) theta.iid <- out$theta.iid
  ## }}}
  }  else stop("score.methods = optimize(dim=1) nlm nlminb fisher.scoring\n"); 


## {{{ handling output
  robvar.theta <- NULL
  if (iid==1) {
     theta.iid <- out$theta.iid %*% hessi
     robvar.theta  <- (t(theta.iid) %*% theta.iid) 
  }
  var.theta <- hessi
  if (!is.null(colnames(theta.des))) thetanames <- colnames(theta.des) else thetanames <- rep("intercept",ptheta)
  ud <- list(theta=theta,score=score,hess=hess,hessi=hessi,var.theta=var.theta,model=model,robvar.theta=robvar.theta,
             theta.iid=theta.iid,thetanames=thetanames,loglike=-logl,score1=score1,Dscore=out$Dscore,margsurv=psurvmarg); 
  class(ud)<-"twostage" 
  attr(ud, "Formula") <- formula
  attr(ud, "Clusters") <- clusters
  attr(ud,"sym")<-sym; 
  attr(ud,"var.link")<-var.link; 
  attr(ud,"antpers")<-antpers; 
  attr(ud,"antclust")<-antclust; 
  attr(ud, "Type") <- model
  return(ud);
  ## }}}

} ## }}}

##' @S3method summary twostage
summary.twostage <-function (object,digits = 3,...) { ## {{{
  if (!(inherits(object,"twostage"))) stop("Must be a Two-Stage object")
  
  var.link<-attr(object,"var.link");
  if (object$model=="plackett") cat("Dependence parameter for Plackett model \n"); 
  if (object$model=="clayton.oakes") cat("Dependence parameter for Clayton-Oakes model \n"); 

  if (sum(abs(object$score)>0.0001) ) {
	  cat("    Variance parameters did not converge, allow more iterations.\n"); 
	  cat(paste("    Score:",object$score,"  \n")); 
  }

  coefs <- coef.twostage(object,...);

  res <- list(estimates=coefs, type=attr(object,"Type"))
  class(res) <- "summary.twostage"
  res
} ## }}}

##' @S3method coef twostage
coef.twostage <- function(object,var.link=NULL,...)
{ ## {{{
  theta <- object$theta
  if (is.null(var.link))
     if (attr(object,"var.link")==1) vlink <- 1 else vlink <- 0
     else vlink <- var.link
  se<-diag(object$var.theta)^0.5
  res <- cbind(theta, se )
  wald <- theta/se
  waldp <- (1 - pnorm(abs(wald))) * 2
  library(numDeriv)
  if (object$model=="plackett") {
  spearman <- alpha2spear(theta,link=vlink)
  Dspear <- jacobian(alpha2spear,theta,link=vlink) 
  var.spearman <- Dspear %*% object$var.theta %*%  Dspear
  se.spearman <- diag(var.spearman)^.5
  res <- as.matrix(cbind(res, wald, waldp,spearman,se.spearman))
  if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Spearman Corr.","SE")
  else colnames(res) <- c("Coef.", "SE","z", "P-val","Spearman Corr.","SE")
  if (!is.null(object$thetanames)) rownames(res)<-object$thetanames
  }
  if (object$model=="clayton.oakes") {
  kendall <- alpha2kendall(theta,link=vlink)
  Dken <- jacobian(alpha2kendall,theta,link=vlink) 
  var.kendall<- Dken %*% object$var.theta %*%  Dken
  se.kendall <- diag(var.kendall)^.5
  res <- as.matrix(cbind(res, wald, waldp,kendall,se.kendall))
  if (vlink==1) colnames(res) <- c("log-Coef.", "SE","z", "P-val","Kendall tau","SE")
  else colnames(res) <- c("Coef.", "SE","z", "P-val","Kendall tau","SE")
  if (!is.null(object$thetanames)) rownames(res)<-object$thetanames
  }

  return(res)
} ## }}}

##' @export
alpha2spear <- function(theta,link=1) {
   if (link==1) theta <- exp(theta)
   if (theta!=1) return( (theta+1)/(theta-1) -2* theta* log(theta)/ (theta-1)^2)
   else return(0)
}

##' @export
alpha2kendall <- function(theta,link=0) { 
   if (link==1) theta <- exp(theta)
   return(1/(1+2/theta)) 
}

##' @S3method print twostage
print.twostage<-function(x,digits=3,...)
{ ## {{{
  print(x$call); 
  cat("\n")
  print(summary(x)); 
} ## }}}

##' @S3method plot twostage
plot.twostage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
		level=0.05, 
		start.time=0,stop.time=0,add.to.plot=FALSE,mains=TRUE,
                xlab="Time",ylab ="Cumulative regression function",...) 
{ ## {{{
  if (!(inherits(x, 'two.stage'))) stop("Must be a Two-Stage object")
  object <- x; rm(x);  
 
  B<-object$cum; V<-object$var.cum; p<-dim(B)[[2]]; 
  if (robust>=1) V<-object$robvar.cum; 

  if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1
  if (stop.time==0) stop.time<-max(B[,1]);

  med<-B[,1]<=stop.time & B[,1]>=start.time
  B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
  V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
  Vrob<-object$robvar.cum; 
  Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 

  c.alpha<- qnorm(1-level/2)
  for (v in comp) { 
    c.alpha<- qnorm(1-level/2)
    est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
    if (add.to.plot==FALSE) 
      {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
        if (mains==TRUE) title(main=colnames(B)[v]); }
    else lines(B[,1],est,type="s"); 
    if (pointwise.ci>=1) {
      lines(B[,1],ul,lty=pointwise.ci,type="s");
      lines(B[,1],nl,lty=pointwise.ci,type="s"); }
    if (robust>=1) {
      lines(B[,1],ul,lty=robust,type="s"); 
      lines(B[,1],nl,lty=robust,type="s"); }
    abline(h=0); 
  }
}  ## }}}

##' @S3method predict twostage
predict.twostage <- function(object,X=NULL,Z=NULL,times=NULL,times2=NULL,theta.des=NULL,diag=TRUE,...)
{ ## {{{
time.coef <- data.frame(object$cum)
if (!is.null(times)) {
cum <- Cpred(object$cum,times);
cum2 <- Cpred(object$cum,times);
} else { cum <- object$cum; cum2 <- object$cum }
if (!is.null(times2)) cum2 <- Cpred(object$cum,times2);

if (is.null(X)) X <- 1;
if (is.null(X) & (!is.null(Z))) { Z <- as.matrix(Z);  X <- matrix(1,nrow(Z),1)}
if (is.null(Z) & (!is.null(X)))  {X <- as.matrix(X);  Z <- matrix(0,nrow(X),1); gamma <- 0}

if (diag==FALSE) {
   time.part <-  X %*% t(cum[,-1]) 
   time.part2 <-  X %*% t(cum2[,-1]) 
   if (!is.null(object$gamma)) { RR <- exp( Z %*% gamma ); 
       cumhaz <- t( t(time.part) * RR ); cumhaz2 <- t( t(time.part2) * RR )}
	    else { cumhaz <- time.part;  cumhaz2 <- time.part2;   }
} else { 
	time.part <-  apply(as.matrix(X*cum[,-1]),1,sum) 
	time.part2 <-  apply(as.matrix(X*cum2[,-1]),1,sum) 
}

if (!is.null(object$gamma)) {
	RR<- exp(Z%*%gamma); 
	cumhaz <- t( t(time.part) * RR );  
	cumhaz2 <- t( t(time.part2) * RR )} else {
		cumhaz <- time.part;  cumhaz2 <- time.part2; 
} 
S1 <- exp(-cumhaz); S2 <- exp(-cumhaz2)

if (attr(object,"var.link")==1) theta  <- exp(object$theta) else theta <- object$theta
if (!is.null(theta.des)) theta <- c(theta.des %*% object$theta)

if (diag==FALSE) St1t2<- (outer(c(S1)^{-(1/theta)},c(S2)^{-(1/theta)},FUN="+") - 1)^(-(theta)) else 
St1t2<- ((S1^{-(1/theta)}+S2^{-(1/theta)})-1)^(-(theta))

out=list(St1t2=St1t2,S1=S1,S2=S2,times=times,times2=times2,theta=theta)
return(out)
} ## }}}

##' @export
surv.boxarea <- function(left.trunc,right.cens,data,timevar="time",status="status",id="id",covars=NULL,num=NULL)
{ ## {{{

  if (is.null(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,num] <- unlist(lapply(idtab,seq_len))
  }

  timevar2 <- paste(timevar,1:2,sep=".")
  status2 <- paste(status,1:2,sep=".")
  covars2 <- NULL; 
  if (length(covars)>0) covars2 <- paste(covars,1,sep=".")
 
  ww0 <- reshape(data[,c(timevar,status,covars,id,num)],direction="wide",idvar=id,timevar=num)[,c(timevar2,status2,covars2,id)]
  mleft <- (1*(ww0[,timevar2[1] ] > left.trunc[1])+ 1*(ww0[,timevar2[2] ] > left.trunc[2])) ==2
###  mleft <- apply(1*(ww0[,timevar2 ] > left.trunc),1,sum)==2
  ww0 <- ww0[!is.na(mleft),]
  mleft <- mleft[!is.na(mleft)]
  if (sum(mleft)==0) stop("No data selected\n"); 
  ww0 <- ww0[mleft,]
  right  <- (ww0[,timevar2] > right.cens)
  ww0[,timevar2[1]][right[,1]] <- right.cens[1]
  ww0[,timevar2[2]][right[,2]] <- right.cens[2]
  ww0[,status2[1]][right[,1]] <- 0
  ww0[,status2[2]][right[,2]] <- 0
  truncvar2 <- c("left.1","left.2")
  ww0[,truncvar2[1]] <- left.trunc[1]
  ww0[,truncvar2[2]] <- left.trunc[2]
###  t(ww0[,truncvar2]) <- left.trunc

  lr.data <- reshape(ww0,direction="long",
		     varying=list(c(timevar2),c(status2),c(truncvar2)),idvar="id",v.names=c(timevar,status,"left"))

return(lr.data)
} ## }}}

##' @export
piecewise.twostage <- function(cut1,cut2,data=sys.parent(),timevar="time",status="status",id="id",covars=NULL,num=NULL,
                          score.method="optimize",
                          Nit=60,detail=0,clusters=NULL,silent=1,weights=NULL,
                          control=list(),theta=NULL,theta.des=NULL,var.link=1,iid=0,
                          step=0.5,model="plackett",data.return=0)
{ ## {{{

ud <- list()
if (missing(cut2)) cut2 <- cut1; 
nc1 <- length(cut1); 
nc2 <- length(cut2)
names1 <- names2 <- c()
theta.mat <- matrix(0,nc1-1,nc2-1); 
se.theta.mat <- matrix(0,nc1-1,nc2-1); 
cor.mat <- matrix(0,nc1-1,nc2-1); 
se.cor.mat <- matrix(0,nc1-1,nc2-1); 
idi <- unique(data[,id]); 
if (iid==1) theta.iid <- matrix(0,length(idi),(nc1-1)*(nc2-1)) else theta.iid <- NULL

k <- 0; 
for (i1 in 2:nc1)
for (i2 in 2:nc2)
{
k <-(i1-2)*(nc2-1)+(i2-1)
 datalr <- surv.boxarea(c(cut1[i1-1],cut2[i2-1]),c(cut1[i1],cut2[i2]),data) 
 boxlr <- list(left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]))
### marg1 <- aalen(Surv(datalr$left,datalr[,timevar],datalr[,status])~+1,data=datalr,n.sim=0,max.clust=NULL,robust=0)
###ud=eval(parse(text=paste("Surv(datalr$left,datalr$",timevar,",datalr$",status,")",sep=""))) 
### marg1 <- aalen(Surv(left,datalr[,timevar],datalr[,status])~+1,data=datalr,n.sim=0,max.clust=NULL,robust=0)
### print(ud)
datalr$tstime <- datalr[,timevar]
datalr$tsstatus <- datalr[,status]
datalr$tsid <- datalr[,id]
marg1 <- aalen(Surv(left,tstime,tsstatus)~+1,data=datalr,n.sim=0,max.clust=NULL,robust=0)
fitlr<-twostage(marg1,data=datalr,clusters=datalr$tsid,model=model,score.method=score.method,
                 Nit=Nit,detail=detail,silent=silent,weights=weights,
                 control=control,theta=theta,theta.des=theta.des,var.link=var.link,iid=iid,step=step)
####
coef <- coef(fitlr)
theta.mat[i1-1,i2-1] <- fitlr$theta
se.theta.mat[i1-1,i2-1] <- fitlr$var.theta^.5
cor.mat[i1-1,i2-1] <- coef[1,5]
se.cor.mat[i1-1,i2-1] <- coef[1,6]
if (data.return==0) 
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr)
if (data.return==1) 
ud[[k]] <- list(index=c(i1,i2),left=c(cut1[i1-1],cut2[i2-1]),right=c(cut1[i1],cut2[i2]),fitlr=fitlr,data=datalr)
if (i2==2) names1 <- c(names1, paste(cut1[i1-1],"-",cut1[i1]))
if (i1==2) names2 <- c(names2, paste(cut2[i2-1],"-",cut2[i2]))
theta <- c(theta,fitlr$theta)

if (iid==1) theta.iid[idi %in% unique(datalr$tsid),k] <-  fitlr$theta.iid 
}

var.thetal <- NULL
if (iid==1)  var.thetal <- t(theta.iid) %*% theta.iid

colnames(cor.mat) <-  colnames(se.cor.mat)  <- colnames(se.theta.mat) <- colnames(theta.mat) <- names1; 
rownames(cor.mat) <-  rownames(se.cor.mat) <-  rownames(se.theta.mat) <- rownames(theta.mat) <- names2; 

ud <- list(model.fits=ud,theta=theta.mat,var.theta=se.theta.mat^2,
	   se.theta=se.theta.mat,thetal=theta,thetal.iid=theta.iid,var.thetal=var.thetal,model=model,
	   cor=cor.mat,se.cor=se.cor.mat); 
class(ud)<-"pc.twostage" 
attr(ud,"var.link")<-var.link; 
attr(ud, "Type") <- model
return(ud);
} ## }}}

##' @S3method summary pc.twostage
summary.pc.twostage <- function(object,var.link=NULL,...)
{ ## {{{
  if (!(inherits(object,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
  
  res <- list(estimates=object$theta,se=object$se.theta,cor=object$cor,se.cor=object$se.cor,
	      model=object$model)
  class(res) <- "summary.pc.twostage"
  attr(res,"var.link")<-attr(object,"var.link"); 
  attr(res, "Type") <- object$model
  res
} ## }}}

##' @S3method print pc.twostage
print.pc.twostage <- function(object,var.link=NULL,...)
{ ## {{{
   if (!(inherits(object,"pc.twostage"))) stop("Must be a Piecewise constant two-Stage object")
   print( summary(object,var.link=var.link,...))
} ## }}}

##' @S3method print summary.pc.twostage
print.summary.pc.twostage <- function(object,var.link=NULL, digits=3,...)
{ ## {{{
  
  if (is.null(var.link)) { if (attr(object,"var.link")==1) vlink <- 1 else vlink <- 0; } else vlink <- var.link
  print(vlink)

  if (object$model=="plackett") cat("Dependence parameter for Plackett model \n"); 
  if (object$model=="clayton.oakes") cat("Dependence parameter for Clayton-Oakes model \n"); 

  if (vlink==1) cat("log-coefficient for dependence parameter (SE) \n")  else cat("Dependence parameter (SE) \n");
  myest <- round(10^digits*(object$estimate))/10^digits;
  myest <- paste(ifelse(myest<0,""," "),myest,sep="")
  mysd <- round(10^digits*(object$se))/10^digits;  
  res <- matrix(paste(myest," (",mysd,")",sep=""),ncol=ncol(object$cor))
  dimnames(res) <- dimnames(object$cor)
  colnames(res) <- paste("",colnames(res))
  print(res,quote=FALSE)  
  cat("\n") 

  if (object$model=="plackett") {cat("Spearman Correlation (SE) \n");cor.type <- "Spearman Correlation"; }
  if (object$model=="clayton.oakes") {cat("Kendall's tau (SE) \n"); cor.type <- "Kendall's tau";}

  myest <- round(10^digits*(object$cor))/10^digits;
  myest <- paste(ifelse(myest<0,""," "),myest,sep="")
  mysd <- round(10^digits*(object$se.cor))/10^digits;  
  res <- matrix(paste(myest," (",mysd,")",sep=""),ncol=ncol(object$cor))
  dimnames(res) <- dimnames(object$cor);
  colnames(res) <- paste("",colnames(res))
  print(res,quote=FALSE)  
  cat("\n") 
} ## }}}

