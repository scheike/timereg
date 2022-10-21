
#' Fit Clayton-Oakes-Glidden Two-Stage model
#' 
#' Fit Clayton-Oakes-Glidden Two-Stage model with Cox-Aalen marginals and
#' regression on the variance parameters.
#' 
#' The model specifikatin allows a regression structure on the variance of the
#' random effects, such it is allowed to depend on covariates fixed within
#' clusters \deqn{ \theta_{k} = Q_{k}^T \nu }. This is particularly useful to
#' model jointly different groups and to compare their variances.
#' 
#' Fits an Cox-Aalen survival model.  Time dependent variables and counting
#' process data (multiple events per subject) are not possible !
#' 
#' The marginal baselines are on the Cox-Aalen form \deqn{ \lambda_{ki}(t) =
#' Y_{ki}(t) ( X_{ki}^T(t) \alpha(t) ) \exp(Z_{ki}^T \beta ) }
#' 
#' The model thus contains the Cox's regression model and the additive hazards
#' model as special cases. (see cox.aalen function for more on this).
#' 
#' The modelling formula uses the standard survival modelling given in the
#' \bold{survival} package. Only for right censored survival data.
#' 
#' The data for a subject is presented as multiple rows or 'observations', each
#' of which applies to an interval of observation (start, stop].  For counting
#' process data with the )start,stop] notation is used the 'id' variable is
#' needed to identify the records for each subject. Only one record per subject
#' is allowed in the current implementation for the estimation of theta.  The
#' program assumes that there are no ties, and if such are present random noise
#' is added to break the ties.
#' 
#' Left truncation is dealt with. Here the key assumption is that the maginals
#' are correctly estimated and that we have a common truncation time within
#' each cluster.
#' 
#' @param margsurv fit of marginal survival cox.aalen model with residuals=2,
#' and resample.iid=1 to get fully correct standard errors. See notaylor below.
#' @param data a data.frame with the variables.
#' @param start.time start of observation period where estimates are computed.
#' @param max.time end of observation period where estimates are computed.
#' Estimates thus computed from [start.time, max.time]. Default is max of data.
#' @param id For timevarying covariates the variable must associate each record
#' with the id of a subject.
#' @param clusters cluster variable for computation of robust variances.
#' @param robust if 0 then totally omits computation of standard errors.
#' @param Nit number of iterations for Newton-Raphson algorithm.
#' @param detail if 0 no details is printed during iterations, if 1 details are
#' given.
#' @param theta starting values for the frailty variance (default=0.1).
#' @param theta.des design for regression for variances. The defauls is NULL
#' that is equivalent to just one theta and the design with only a baseline.
#' @param var.link default "0" is that the regression design on the variances
#' is without a link, and "1" uses the link function exp.
#' @param step step size for Newton-Raphson.
#' @param notaylor if 1 then ignores variation due to survival model, this is
#' quicker and then resample.iid=0 and residuals=0 is ok for marginal survival
#' model that then is much quicker.
#' @param se.clusters cluster variable for sandwich estimator of variance.
#' @return returns an object of type "two.stage". With the following arguments:
#' \item{cum}{cumulative timevarying regression coefficient estimates are
#' computed within the estimation interval.} \item{var.cum}{the martingale
#' based pointwise variance estimates.} \item{robvar.cum}{robust pointwise
#' variances estimates.} \item{gamma}{estimate of parametric components of
#' model.} \item{var.gamma}{variance for gamma.} \item{robvar.gamma}{robust
#' variance for gamma.} \item{D2linv}{inverse of the derivative of the score
#' function from marginal model.} \item{score}{value of score for final
#' estimates.} \item{theta}{estimate of Gamma variance for frailty.}
#' \item{var.theta}{estimate of variance of theta.} \item{SthetaInv}{inverse of
#' derivative of score of theta.} \item{theta.score}{score for theta
#' parameters.}
#' @author Thomas Scheike
#' @references Glidden (2000), A Two-Stage estimator of the dependence
#' parameter for the Clayton Oakes model.
#' 
#' Martinussen and Scheike, Dynamic Regression Models for Survival Data,
#' Springer (2006).
#' @keywords survival
#' @examples
#' 
#' library(timereg)
#' data(diabetes)
#' # Marginal Cox model  with treat as covariate
#' marg <- cox.aalen(Surv(time,status)~prop(treat)+prop(adult)+
#' 	  cluster(id),data=diabetes,resample.iid=1)
#' fit<-two.stage(marg,data=diabetes,theta=1.0,Nit=40)
#' summary(fit)
#' 
#' # using coxph and giving clusters, but SE wittout cox uncetainty
#' margph <- coxph(Surv(time,status)~treat,data=diabetes)
#' fit<-two.stage(margph,data=diabetes,theta=1.0,Nit=40,clusters=diabetes$id)
#' 
#' 
#' # Stratification after adult 
#' theta.des<-model.matrix(~-1+factor(adult),diabetes);
#' des.t<-model.matrix(~-1+factor(treat),diabetes);
#' design.treat<-cbind(des.t[,-1]*(diabetes$adult==1),
#'                     des.t[,-1]*(diabetes$adult==2))
#' 
#' # test for common baselines included here 
#' marg1<-cox.aalen(Surv(time,status)~-1+factor(adult)+prop(design.treat)+cluster(id),
#'  data=diabetes,resample.iid=1,Nit=50)
#' 
#' fit.s<-two.stage(marg1,data=diabetes,Nit=40,theta=1,theta.des=theta.des)
#' summary(fit.s)
#' 
#' # with common baselines  and common treatment effect (although test reject this)
#' fit.s2<-two.stage(marg,data=diabetes,Nit=40,theta=1,theta.des=theta.des)
#' summary(fit.s2)
#' 
#' # test for same variance among the two strata
#' theta.des<-model.matrix(~factor(adult),diabetes);
#' fit.s3<-two.stage(marg,data=diabetes,Nit=40,theta=1,theta.des=theta.des)
#' summary(fit.s3)
#' 
#' # to fit model without covariates, use beta.fixed=1 and prop or aalen function
#' marg <- aalen(Surv(time,status)~+1+cluster(id),
#' 	 data=diabetes,resample.iid=1,n.sim=0)
#' fita<-two.stage(marg,data=diabetes,theta=0.95,detail=0)
#' summary(fita)
#' 
#' # same model but se's without variation from marginal model to speed up computations
#' marg <- aalen(Surv(time,status) ~+1+cluster(id),data=diabetes,
#' 	      resample.iid=0,n.sim=0)
#' fit<-two.stage(marg,data=diabetes,theta=0.95,detail=0)
#' summary(fit)
#' 
#' # same model but se's now with fewer time-points for approx of iid decomp of marginal 
#' # model to speed up computations
#' marg <- cox.aalen(Surv(time,status) ~+prop(treat)+cluster(id),data=diabetes,
#' 	      resample.iid=1,n.sim=0,max.timepoint.sim=5,beta.fixed=1,beta=0)
#' fit<-two.stage(marg,data=diabetes,theta=0.95,detail=0)
#' summary(fit)
#' 
##' @export
two.stage<-function(margsurv,data=parent.frame(),
Nit=60,detail=0,start.time=0,max.time=NULL,id=NULL,clusters=NULL,
robust=1,theta=NULL,theta.des=NULL,var.link=0,step=0.5,notaylor=0,se.clusters=NULL)
{ ## {{{
## {{{ seting up design and variables
rate.sim <- 1; secluster <- NULL

if (!inherits(margsurv,"coxph")) { ## {{{ 
 formula<-attr(margsurv,"Formula");
 beta.fixed <- attr(margsurv,"beta.fixed")
 if (is.null(beta.fixed)) beta.fixed <- 1; 
 ldata<-aalen.des(formula,data=data,model="cox.aalen");
 id <- attr(margsurv,"id"); 
 mclusters <- attr(margsurv,"cluster")
 mclustind <- attr(margsurv,"cluster")
 cluster.call <- attr(margsurv,"cluster.call")
 X<-ldata$X; time<-ldata$time2; Z<-ldata$Z;  status<-ldata$status;
 time2 <- attr(margsurv,"stop"); 
 start <- attr(margsurv,"start")
 antpers<-nrow(X);
 if (beta.fixed==1) Z <- NULL; 
 if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
 if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
 px<-ncol(X);

 if (is.null(clusters) && is.null(mclusters)) 
	 stop("No cluster variabel specified in marginal or twostage call \n"); 
 if (is.null(clusters)) { clusters <- mclusters; cluster.call <- cluster.call} else {cluster.call <- clusters;}

 if (is.null(se.clusters)) secluster <- clusters;
 antsecluster <- length(unique(secluster))
 if (is.numeric(secluster)) secluster <-  sindex.prodlim(unique(secluster),secluster)-1 else  {
      seclusters <- as.integer(factor(clusters, labels = 1:antsecluster))-1
 }

### print("two-stage"); print(head(cluster.call))
 if (is.null(cluster.call)) notaylor <- 1
 if (is.null(margsurv$gamma.iid)) notaylor <- 1

 ## }}} 
} else { ## coxph ## {{{ 
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
	 start <- Y[, 1]; time2 <- Y[, 2];status <- Y[, 3];
        }
   Z <- matrix(1,antpers,length(coef(margsurv)));

   if (is.null(clusters)) stop("must give clusters for coxph\n");
   cluster.call <- clusters 
   X <- matrix(1,antpers,1); ### Z <- matrix(0,antpers,1); ### no use for these
   px <- 1; pz <- ncol(Z); 
   beta.fixed <- 0
   semi <- 1
   start.time <- 0
}  ## }}} 


  if (any(is.na(clusters))) stop("Missing values in cluster varaibles\n"); 
  out.clust <- cluster.index.timereg(clusters);  
  clusters <- out.clust$clusters
  maxclust <- out.clust$maxclust 
  antclust <- out.clust$antclust
  idiclust <- out.clust$idclust
  cluster.size <- out.clust$cluster.size
###  if (anyNA(idiclust)) idiclust[is.na(idiclust)] <- 0

  ### setting secluster after cluster.index call to deal with characters 
  if (inherits(margsurv,"coxph")) {
  if (is.null(se.clusters) & is.null(secluster) ) secluster <- clusters;
  antsecluster <- length(unique(secluster))
  if (is.numeric(secluster)) secluster <-  sindex.prodlim(unique(secluster),secluster)-1 else  {
       clusters <- as.integer(factor(clusters, labels = 1:antsecluster))-1
  }
  }

 if (length(clusters)!=length(secluster)) stop("length of se.clusters not consistent with cluster length\n"); 
   
  if (sum(abs(start))>0) lefttrunk <- 1  else lefttrunk <- 0;  
  cumhazleft <- 0; 
  RR <-  rep(1,antpers); 

  update <- 1;
  if (update==0) { ## {{{
  if ((attr(margsurv,"residuals")!=2) || (lefttrunk==1)) { ### compute cum hazards in time point infty; 
	  nn <- nrow(margsurv$cum) 
	  cum <- Cpred(margsurv$cum,time2)[,-1]
	  if (npar==TRUE) cumhaz <- apply(cum*X,1,sum)
	  if (npar==FALSE) cumhaz <- apply(cum*X,1,sum)*exp( Z %*% margsurv$gamma)
	  if (lefttrunk==1) {
	     cum <- Cpred(margsurv$cum,start)[,-1]
	     cumhazleft <- apply(cum*X,1,sum)
	     if (npar==TRUE) cumhazleft <-  cumhazleft
	     if (npar==FALSE) cumhazleft <- cumhazleft * exp( Z %*% margsurv$gamma)
	  } 
  } else { residuals<-margsurv$residuals$dM; cumhaz<-status-residuals; }
  } ## }}}

  if (update==1) 
  if (inherits(margsurv,c("aalen","cox.aalen")))  { ## {{{
     if ((attr(margsurv,"residuals")!=2) || (lefttrunk==1)) { 
         resi <- residualsTimereg(margsurv,data=data) 
         residuals <- resi$residuals; 
	 cumhaz <- resi$cumhaz; 
	 cumhazleft <- resi$cumhazleft; 
	 RR <- resi$RR
     } else { residuals <- margsurv$residuals$dM; 
              cumhaz <- status-residuals; 
	      if (inherits(margsurv,"cox.aalen")) RR  <- exp( Z %*% margsurv$gamma)
     }
  }
  else if (inherits(margsurv,"coxph")) {
       notaylor <- 1
       residuals <- residuals(margsurv)
       cumhaz <- status-residuals
       cumhazleft <- rep(0,antpers)
       RR<- exp(margsurv$linear.predictors-sum(margsurv$means*coef(margsurv)))
        if ((lefttrunk==1)) { 
###           baseout <- basehaz(margsurv,centered=FALSE); 
             sfit <- survfit(margsurv, se.fit=FALSE)   
             zcoef <- ifelse(is.na(coef(margsurv)), 0, coef(margsurv))
             offset <- sum(margsurv$means * zcoef)
             chaz <- sfit$cumhaz * exp(-offset)
             cum <- cbind(sfit$time,chaz)  
	   cum <- Cpred(cum,start)[,2]
	   cumhazleft <- cum * RR 
	}
  } ## }}}

###  print(head(cbind(residuals,cumhaz,RR,time2,status)))

  ratesim<-rate.sim; 
  inverse<-var.link
  pxz <- px + pz;
  times<-c(start.time,time2[status==1]); 
  times<-sort(times);
  if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
  times<-times[times<maxtimes]
  Ntimes <- sum(status[time2<maxtimes])+1; 


  Biid<-c(); gamma.iid <- 0; 
  if (is.null(margsurv$B.iid)) notaylor <- 1; 
  if (notaylor==0) {
    nBiid <- length(margsurv$B.iid)
    if (nBiid!=antsecluster) stop("Number of clusters for marginal models must be consistent with se.cluster for standard error sandwich\n"); 
    for (i in 1:nBiid) Biid<-cbind(Biid,margsurv$B.iid[[i]]); 
    if (!is.null(margsurv$gamma.iid)) gamma.iid<-margsurv$gamma.iid;
    if (is.null(margsurv$time.sim.resolution)) { 
	   time.group <- (1:nrow(Biid))-1; 
           maxtimesim <- nrow(Biid); 
	   timereso <- margsurv$cum[,1] 
    }  
    else {
      timereso <- margsurv$time.sim.resolution
       qqc <- cut(times, breaks = margsurv$time.sim.resolution, include.lowest = TRUE)    
       time.group <- as.integer(factor(qqc, labels = 1:(nrow(Biid)-1)))
       maxtimesim <- nrow(Biid); 
    } 
    if (inherits(margsurv,"cox.aalen"))  { 
       times <- margsurv$time.sim.resolution 
       Ntimes <-length(times)
    }
  } else {time.group <- 1; maxtimesim <- 1; timereso <- 1}

  if (is.null(theta.des)==TRUE) ptheta<-1; 
  if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) else
  theta.des<-as.matrix(theta.des); 
  ptheta<-ncol(theta.des); 
  if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");

  if (is.null(theta)==TRUE) {
      if (var.link==1) theta<- rep(-1,ptheta); 
      if (var.link==0) theta<- rep(exp(-5),ptheta); 
  }
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
  theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

  if (maxclust==1) stop("No clusters !, maxclust size=1\n"); 
  theta.iid <- matrix(0,antsecluster,ptheta)
  ## }}}


###  print(ant.clust)
###  print(dim(idiclust))
###  print(sum(is.na(idiclust)))
###  print(table(clusters))
###  print(table(secluster))
  
  DUbeta <- matrix(0,pz,ptheta); 
  nparout <- .C("twostagereg", 
        as.double(times), as.integer(Ntimes), as.double(X),
   	as.integer(antpers), as.integer(px), as.double(Z), 
	as.integer(antpers), as.integer(pz), as.integer(antpers),         ## 9 
	as.double(start),as.double(time2), as.integer(Nit), 
	as.integer(detail), as.integer(id), as.integer(status),           ## 15
	as.integer(ratesim), as.integer(robust), as.integer(clusters),    
	as.integer(antclust), as.integer(beta.fixed), as.double(theta),
	as.double(var.theta), as.double(theta.score), as.integer(inverse), 
	as.integer(cluster.size),                                          ## 25
	as.double(theta.des), as.integer(ptheta), as.double(Stheta),
	as.double(step), as.integer(idiclust), as.integer(notaylor),
	as.double(gamma.iid),as.double(Biid),as.integer(semi), as.double(cumhaz) ,
	as.double(cumhazleft),as.integer(lefttrunk),as.double(RR),
	as.integer(maxtimesim),as.integer(time.group),as.integer(secluster),
	as.integer(antsecluster),as.double(theta.iid), as.double(timereso),as.double(DUbeta),PACKAGE = "timereg")

## {{{ handling output
   gamma <- margsurv$gamma
   Varbeta <- margsurv$var.gamma; RVarbeta <- margsurv$robvar.gamma;
   score <- margsurv$score; Iinv <- margsurv$D2linv;
   cumint <- margsurv$cum; vcum <- margsurv$var.cum; Rvcu <- margsurv$robvar.cum;

   theta<-matrix(nparout[[21]],ptheta,1);  
   var.theta<-matrix(nparout[[22]],ptheta,ptheta); 
   theta.score<-nparout[[23]]; 
   SthetaI<-matrix(nparout[[28]],ptheta,ptheta); 
   theta.iid  <- matrix(nparout[[43]],antsecluster,ptheta); 
   theta.iid <- theta.iid %*% SthetaI

###  DUbeta <- matrix(nparout[[45]],pz,ptheta);  
###  if (is.null(call.secluster) & is.null(max.clust)) rownames(theta.iid) <- unique(cluster.call) else rownames(theta.iid) <- unique(se.clusters)

   ud <- list(cum = cumint, var.cum = vcum, robvar.cum = Rvcu, 
       gamma = gamma, var.gamma = Varbeta, robvar.gamma = RVarbeta, 
       D2linv = Iinv, score = score,  theta=theta,var.theta=var.theta,
       SthetaInv=SthetaI,theta.score=theta.score,theta.iid=theta.iid)

  ptheta<-length(ud$theta); 
  if (ptheta>1) {
                rownames(ud$theta)<-colnames(theta.des);
                names(ud$theta.score)<-colnames(theta.des); } else { 
		names(ud$theta.score)<- rownames(ud$theta)<-"intercept" } 


  attr(ud,"Call")<- match.call(); 
  class(ud)<-"two.stage"
  attr(ud,"Formula")<-formula;
  attr(ud,"id")<-id;
  attr(ud,"cluster")<-clusters;
  attr(ud,"cluster.call")<-cluster.call;
  attr(ud,"secluster")<-secluster;
  attr(ud,"start")<-start; 
  attr(ud,"time2")<-time2; 
  attr(ud,"var.link")<-var.link
  attr(ud,"beta.fixed")<-beta.fixed
  attr(ud,"marg.model")<-class(margsurv)
###  attr(ud,"DUbeta")<-DUbeta

  return(ud) 
  ## }}} 
} ## }}} 
  
##' @export
summary.two.stage<-function (object,digits=3,...) { ## {{{ 

  if (!(inherits(object, 'two.stage') )) stop("Must be a Two-Stage object")
  prop<-TRUE; 
  if (is.null(object$prop.odds)==TRUE) p.o<-FALSE else p.o<-TRUE
    
  var.link<-attr(object,"var.link");
  cat("Dependence parameter for Clayton-Oakes-Glidden  model\n"); 

  if (sum(abs(object$theta.score)>0.000001) ) 
    cat("Variance parameters did not converge, allow more iterations\n\n"); 

  resdep <- coef.two.stage(object,...)

  prmatrix(resdep[,1:6,drop=FALSE]); cat("   \n");  

  prmatrix(resdep[,7:9,drop=FALSE]); cat("   \n");  

  if (attr(object,"marg.model")!="coxph")
  if (attr(object,"beta.fixed")==0) { ## {{{ 
###  cat("Marginal Cox-Aalen model fit\n\n"); 
  if (sum(abs(object$score)>0.000001) && sum(object$gamma)!=0) 
    cat("Marginal model did not converge, allow more iterations\n\n"); 
###  if (prop) {
###    if (p.o==FALSE) cat("Proportional Cox terms :  \n") else  cat("Covariate effects \n")
###
###    out=coef.two.stage(object,digits=digits);
###    out=signif(out,digits=digits)
###    print(out)
###
###  }

  } ## }}} 
###   cat("   \n");  cat("  Call: \n"); dput(attr(object, "Call")); 
  cat("\n");
} ## }}}

##' @export
print.two.stage <- function (x,digits = 3,...) { ## {{{
	summary.two.stage(x,digits=digits,...)
###  if (!(inherits(x, 'two.stage') )) stop("Must be a Two-Stage object")
###  cat(" Two-stage estimation for Clayton-Oakes-Glidden  model\n"); 
###  cat(" Marginals of Cox-Aalen form, dependence by variance of Gamma distribution\n\n");  
###  object <- x; rm(x);
###  
###  cat(" Nonparametric components : "); 
###  cat(colnames(object$cum)[-1]); cat("   \n");  
###  if (!is.null(object$gamma)) {
###    cat(" Parametric components :  "); cat(rownames(object$gamma)); 
###    cat("   \n");
###  } 
###  cat("   \n");  
###
###  cat(" Call: \n");
###  print(attr(object,'Call'))
} ## }}}

##' @export
vcov.two.stage <- function(object, ...) {
  rv <- object$robvar.gamma
  if (!identical(rv, matrix(0, nrow = 1L, ncol = 1L))) rv # else return NULL
}

##' @export
coef.two.stage<-function(object,digits=3,d2logl=1,alpha=0.05,...) { ## {{{ 

  if (!(inherits(object, 'two.stage') )) stop("Must be a Two-Stage object")
  var.link <- attr(object,"var.link")
  ptheta<-nrow(object$theta)
  sdtheta<-diag(object$var.theta)^.5

  if (var.link==0) {
      vari<-object$theta
      sdvar<-diag(object$var.theta)^.5
      upper <- object$theta-qnorm(alpha/2)*sdvar
      lower <- object$theta+qnorm(alpha/2)*sdvar
  }
  else {
      vari<-exp(object$theta)
      sdvar<-vari*diag(object$var.theta)^.5
      upper <- exp(object$theta-qnorm(alpha/2)*sdtheta)
      lower <- exp(object$theta+qnorm(alpha/2)*sdtheta)
  }
  dep<-cbind(object$theta[,1],sdtheta)
  walddep<-object$theta[,1]/sdtheta; 
  waldpdep<-(1-pnorm(abs(walddep)))*2

  kendall<-1/(1+2/vari) 
  kendall.ll<-1/(1+2/(object$theta+qnorm(alpha/2)*sdvar)) 
  kendall.ul<-1/(1+2/(object$theta-qnorm(alpha/2)*sdvar)) 
  if (var.link==0) resdep<-signif(as.matrix(cbind(dep,lower,upper,walddep,waldpdep,kendall,kendall.ll,kendall.ul)),digits)
  else resdep<-signif(as.matrix(cbind(dep,lower,upper,walddep,waldpdep,vari,sdvar,kendall,kendall.ll,kendall.ul)),digits);

  slower <- paste("lower",signif(100*alpha/2,2),"%",sep="")
  supper <- paste("upper",signif(100*(1-alpha/2),3),"%",sep="")
  if (var.link==0) colnames(resdep) <- c("Variance","SE",slower,supper,"z","P-val","Kendall's tau",slower,supper) 
  else colnames(resdep)<-c("log(Variance)","SE",slower,supper,"z","P-val","Variance","SE Var.","Kendall's tau",slower,supper)

###  prmatrix(resdep); cat("   \n");  
  return(resdep)
} ## }}} 

##' @export
plot.two.stage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
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

##' @export
predict.two.stage <- function(object,X=NULL,Z=NULL,times=NULL,times2=NULL,X2=NULL,Z2=NULL,
			      theta=NULL,theta.des=NULL,diag=TRUE,...)
{ ## {{{
time.coef <- data.frame(object$cum)
if (!is.null(times))   cum <- Cpred(object$cum,times)  else cum <- object$cum; 
if (!is.null(times2)) cum2 <- Cpred(object$cum,times2) else cum2 <- object$cum;

if (is.null(Z) & (!is.null(object$gamma))) Z <- matrix(0,1,nrow(object$gamma));
if (is.null(X) & (!is.null(Z))) { Z <- as.matrix(Z);  X <- matrix(1,nrow(Z),1)}
if (is.null(Z) & (!is.null(X)))  {X <- as.matrix(X);  Z <- matrix(0,nrow(X),1); gamma <- 0}
if (is.null(X)) X <- 1;

X2 <- X;
Z2 <- Z2;

if (diag==FALSE) {
   time.part <-  X %*% t(cum[,-1]) 
   time.part2 <-  X2 %*% t(cum2[,-1]) 
   if (!is.null(object$gamma)) { 
	   gamma <- object$gamma
	   RR <- exp( Z %*% gamma ); 
	   RR2 <- exp( Z2 %*% gamma ); 
       cumhaz <- t( t(time.part) * RR ); cumhaz2 <- t( t(time.part2) * RR2 )}
	    else { cumhaz <- time.part;  cumhaz2 <- time.part2;   }
} else { 
	time.part <-  apply(as.matrix(X*cum[,-1]),1,sum) 
	time.part2 <-  apply(as.matrix(X2*cum2[,-1]),1,sum) 
}

if (!is.null(object$gamma)) {
	RR<- exp(Z%*%object$gamma); 
	RR2 <- exp(Z2 %*%object$gamma); 
	cumhaz <- c((time.part) * RR) ;  
	cumhaz2 <- c((time.part2) * RR2); 
} else {
	cumhaz <- c(time.part);  cumhaz2 <- c(time.part2); 
} 
S1 <- pmin(1,exp(-cumhaz)); S2 <- pmin(1,exp(-cumhaz2))
###print(length(S1))
###print(length(S2))

if (is.null(theta))  theta <- object$theta
if (!is.null(theta.des)) theta <- c(theta.des %*% theta)
if (attr(object,"var.link")==1) theta  <- exp(theta) 

### theta is variance 
if (diag==FALSE) St1t2<- (outer(c(S1)^{-(theta)},c(S2)^{-(theta)},FUN="+") - 1)^(-(1/theta)) else 
St1t2<- ((S1^{-(theta)}+S2^{-(theta)})-1)^(-(1/theta))
###St1t2<- ((S1^{-(1/theta)}+S2^{-(1/theta)})-1)^(-(theta))

out=list(St1t2=St1t2,S1=S1,S2=S2,times=times,times2=times2,theta=theta)
return(out)
} ## }}}

