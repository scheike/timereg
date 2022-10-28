#' Identifies the multiplicative terms in Cox-Aalen model and proportional
#' excess risk model
#' 
#' Specifies which of the regressors that belong to the multiplicative part of
#' the Cox-Aalen model
#' 
#' \deqn{ \lambda_{i}(t) = Y_i(t) ( X_{i}^T(t) \alpha(t) ) \exp(Z_{i}^T(t)
#' \beta ) }{} for this model prop specified the covariates to be included in
#' \eqn{Z_{i}(t)}
#' @param x variable
#' 
#' @author Thomas Scheike
#' @keywords survival
prop<-function(x) x



#' Fit Cox-Aalen survival model
#' 
#' Fits an Cox-Aalen survival model.  Time dependent variables and counting
#' process data (multiple events per subject) are possible.
#' 
#' \deqn{ \lambda_{i}(t) = Y_i(t) ( X_{i}^T(t) \alpha(t) ) \exp(Z_{i}^T \beta ) }{}
#' 
#' The model thus contains the Cox's regression model as special case.
#' 
#' To fit a stratified Cox model it is important to parametrize the baseline
#' apppropriately (see example below).
#' 
#' Resampling is used for computing p-values for tests of time-varying effects.
#' Test for proportionality is considered by considering the score processes
#' for the proportional effects of model.
#' 
#' The modelling formula uses the standard survival modelling given in the
#' \bold{survival} package.
#' 
#' The data for a subject is presented as multiple rows or 'observations', each
#' of which applies to an interval of observation (start, stop].  For counting
#' process data with the )start,stop] notation is used, the 'id' variable is
#' needed to identify the records for each subject. The program assumes that
#' there are no ties, and if such are present random noise is added to break
#' the ties.
#' 
#' @param formula a formula object with the response on the left of a '~'
#' operator, and the independent terms on the right as regressors. The response
#' must be a survival object as returned by the `Surv' function. Terms with a
#' proportional effect are specified by the wrapper prop(), and cluster
#' variables (for computing robust variances) by the wrapper cluster().
#' @param data a data.frame with the variables.
#' @param start.time start of observation period where estimates are computed.
#' @param max.time end of observation period where estimates are computed.
#' Estimates thus computed from [start.time, max.time]. Default is max of data.
#' @param robust to compute robust variances and construct processes for
#' resampling. May be set to 0 to save memory and time, in particular for
#' rate.sim=1.
#' @param id For timevarying covariates the variable must associate each record
#' with the id of a subject.
#' @param clusters cluster variable for computation of robust variances.
#' @param n.sim number of simulations in resampling.
#' @param weighted.test to compute a variance weighted version of the
#' test-processes used for testing time-varying effects.
#' @param residuals to returns residuals that can be used for model validation
#' in the function cum.residuals. Estimated martingale increments (dM) and
#' corresponding time vector (time). When rate.sim=1 returns estimated
#' martingales, dM_i(t) and if rate.sim=0, returns a matrix of dN_i(t).
#' @param covariance to compute covariance estimates for nonparametric terms
#' rather than just the variances.
#' @param resample.iid to return i.i.d. representation for nonparametric and
#' parametric terms. based on counting process or martingale resduals
#' (rate.sim).
#' @param beta starting value for relative risk estimates.
#' @param Nit number of iterations for Newton-Raphson algorithm.
#' @param detail if 0 no details is printed during iterations, if 1 details are
#' given.
#' @param weights weights for weighted analysis.
#' @param rate.sim rate.sim=1 such that resampling of residuals is based on
#' estimated martingales and thus valid in rate case, rate.sim=0 means that
#' resampling is based on counting processes and thus only valid in intensity
#' case.
#' @param beta.fixed option for computing score process for fixed relative risk
#' parameter
#' @param max.clust sets the total number of i.i.d. terms in i.i.d.
#' decompostition. This can limit the amount of memory used by coarsening the
#' clusters. When NULL then all clusters are used.  Default is 1000 to save
#' memory and time.
#' @param exact.deriv if 1 then uses exact derivative in last iteration, if 2
#' then uses exact derivate for all iterations, and if 0 then uses
#' approximation for all computations and there may be a small bias in the
#' variance estimates. For Cox model always exact and all options give same
#' results.
#' @param silent if 1 then opppresses some output.
#' @param max.timepoint.sim considers only this resolution on the time scale
#' for simulations, see time.sim.resolution argument
#' @param basesim 1 to get simulations for cumulative baseline, including tests
#' for contant effects.
#' @param offsets offsets for analysis on log-scale. RR=exp(offsets+ x beta).
#' @param strata future option for making strata in a different day than
#' through X design in cox-aalen model (~-1+factor(strata)).
#' @param propodds if 1 will fit the proportional odds model. Slightly less
#' efficient than prop.odds() function but much quicker, for large data this
#' also works.
#' @param caseweight these weights have length equal to number of jump times,
#' and are multiplied all jump times dN.  Useful for getting the program to fit
#' for example the proportional odds model or frailty models.
#' @return returns an object of type "cox.aalen". With the following arguments:
#' \item{cum}{cumulative timevarying regression coefficient estimates are
#' computed within the estimation interval. } \item{var.cum}{the martingale
#' based pointwise variance estimates.  } \item{robvar.cum}{robust pointwise
#' variances estimates.  } \item{gamma}{estimate of parametric components of
#' model.  } \item{var.gamma}{variance for gamma sandwhich estimator based on
#' optional variation estimator of score and 2nd derivative.}
#' \item{robvar.gamma}{robust variance for gamma.  } \item{residuals}{list with
#' residuals.} \item{obs.testBeq0}{observed absolute value of supremum of
#' cumulative components scaled with the variance.}
#' \item{pval.testBeq0}{p-value for covariate effects based on supremum test.}
#' \item{sim.testBeq0}{resampled supremum values.} \item{obs.testBeqC}{observed
#' absolute value of supremum of difference between observed cumulative process
#' and estimate under null of constant effect.} \item{pval.testBeqC}{p-value
#' based on resampling.} \item{sim.testBeqC}{resampled supremum values.}
#' \item{obs.testBeqC.is}{observed integrated squared differences between
#' observed cumulative and estimate under null of constant effect.}
#' \item{pval.testBeqC.is}{p-value based on resampling.}
#' \item{sim.testBeqC.is}{resampled supremum values.}
#' \item{conf.band}{resampling based constant to construct robust 95\% uniform
#' confidence bands. } \item{test.procBeqC}{observed test-process of difference
#' between observed cumulative process and estimate under null of constant
#' effect over time.  } \item{sim.test.procBeqC}{list of 50 random realizations
#' of test-processes under null based on resampling.}
#' \item{covariance}{covariances for nonparametric terms of model.}
#' \item{B.iid}{Resample processes for nonparametric terms of model.}
#' \item{gamma.iid}{Resample processes for parametric terms of model.}
#' \item{loglike}{approximate log-likelihood for model, similar to Cox's
#' partial likelihood. Only computed when robust=1.} \item{D2linv}{inverse of
#' the derivative of the score function.} \item{score}{value of score for final
#' estimates.} \item{test.procProp}{observed score process for proportional
#' part of model.} \item{var.score}{variance of score process (optional
#' variation estimator for beta.fixed=1 and robust estimator otherwise).}
#' \item{pval.Prop}{p-value based on resampling.} \item{sim.supProp}{re-sampled
#' absolute supremum values.} \item{sim.test.procProp}{list of 50 random
#' realizations of test-processes for proportionality under the model based on
#' resampling.}
#' @author Thomas Scheike
#' @references Martinussen and Scheike, Dynamic Regression Models for Survival
#' Data, Springer (2006).
#' @keywords survival
#' @examples
#' 
#' library(timereg)
#' data(sTRACE)
#' # Fits Cox model 
#' out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
#' prop(vf)+prop(chf)+prop(diabetes),data=sTRACE)
#' 
#' # makes Lin, Wei, Ying test for proportionality
#' summary(out)
#' par(mfrow=c(2,3))
#' plot(out,score=1) 
#' 
#' # Fits stratified Cox model 
#' out<-cox.aalen(Surv(time,status==9)~-1+factor(vf)+ prop(age)+prop(sex)+
#' 	       prop(chf)+prop(diabetes),data=sTRACE,max.time=7,n.sim=100)
#' summary(out)
#' par(mfrow=c(1,2)); plot(out); 
#' # Same model, but needs to invert the entire marix for the aalen part: X(t) 
#' out<-cox.aalen(Surv(time,status==9)~factor(vf)+ prop(age)+prop(sex)+
#' 	       prop(chf)+prop(diabetes),data=sTRACE,max.time=7,n.sim=100)
#' summary(out)
#' par(mfrow=c(1,2)); plot(out); 
#' 
#' 
#' # Fits Cox-Aalen model 
#' out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
#'                vf+chf+prop(diabetes),data=sTRACE,max.time=7,n.sim=100)
#' summary(out)
#' par(mfrow=c(2,3))
#' plot(out)
#' 
cox.aalen<-function(formula=formula(data),data=parent.frame(),
beta=NULL,Nit=20,detail=0,start.time=0,max.time=NULL, id=NULL, 
clusters=NULL, n.sim=500, residuals=0,robust=1,
weighted.test=0,covariance=0,resample.iid=1,weights=NULL,
rate.sim=1,beta.fixed=0,max.clust=1000,exact.deriv=1,silent=1,
max.timepoint.sim=100,basesim=0,offsets=NULL,strata=NULL,propodds=0,caseweight=NULL)
{ ## {{{
# {{{ set up variables 
  if (n.sim == 0) sim <- 0 else sim <- 1
  if (resample.iid==1 & robust==0) {resample.iid <- 0;}
  if (covariance==1 & robust==0) {covariance<-0;cat("Covariance of baseline only for robust=1\n"); }
  if (robust==0 ) { n.sim <- 0; sim<-0;}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  if (beta.fixed==1) Nit<-1; 
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  m$robust<-m$start.time<- m$scaleLWY<-m$weighted.test<-m$beta<-m$Nit<-m$detail<-
	  m$max.time<-m$residuals<-m$n.sim<-m$id<-m$covariance<-m$resample.iid<-
	  m$clusters<-m$rate.sim<-m$beta.fixed<- m$max.clust <- m$exact.deriv <- 
	  m$silent <- m$max.timepoint.sim <- m$silent <- m$basesim  <- 
	  m$offsets <- m$strata <- m$propodds <- m$caseweight <- NULL

  special <- c("prop","cluster")
  Terms <- if(missing(data)) terms(formula, special)
  else  terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  m <- na.omit(m)
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms,model="cox.aalen")
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
  pxz <- px + pz;

###  if ( (nrow(Z)!=nrow(data)) && (!is.null(id))) stop("Missing values in design matrix not allowed with id\n"); 
###  if (nrow(Z)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 
  if (!is.null(id)) {
	  if (length(id)!=nrow(Z)) stop("id length and data not the same\n"); 
  }

  ### if clusters=null perhaps given through cluster() special that can be NULL also 
  if (is.null(clusters)) clusters <- des$clusters  
  cluster.call<-clusters; 

  survs<-read.surv(m,id,npar,clusters,start.time,max.time,model="cox.aalen",silent=silent)
  times<-survs$times;
  id<-id.call<-survs$id.cal;
  ### if no clusters then return "id" as cluster variable
  clusters<- gclusters <- survs$clusters; 

  start.call <- start <-  survs$start; 
  stop.call <- time2 <- survs$stop; 
  status<-survs$status;
  status.call <- status
  orig.max.clust <- survs$antclust
  nobs <- nrow(X); 
  if (is.null(weights)) weights <- rep(1,nrow(X));  
###  weights <- rep(1,nrow(X)); 
  if (length(weights)!=nrow(X)) stop("Lengths of weights and data do not match\n"); 
  if (is.null(offsets)) offsets <- rep(0,nrow(X));  
  offsets.call <- offsets; 
  weights.call <- weights; 
  if (length(offsets)!=nrow(X)) stop("Lengths of offsets and data do not match\n"); 
  if (!is.null(strata))  {
  if (length(strata)!=nrow(X)) stop("Lengths of strata and data do not match\n"); 
    iids <- unique(strata)
    antiid <- length(iids)
    if (is.numeric(strata)) strata <-  sindex.prodlim(iids,strata)-1
    else strata<- as.integer(factor(strata, labels = seq(antiid)))-1
  } 

  if (rate.sim==1 && robust==1) 
  if ((!is.null(max.clust))) if (max.clust<survs$antclust) {
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	gclusters <- clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
###	clusters <- as.integer(factor(qqc, labels = 1:max.clust)) -1
	survs$antclust <- max.clust    
  }                                                         

### if rate.sim==0 and no clusters then it suffices with jump processes and the rest is 0
  Ntimes <- sum(status)
  if (rate.sim==0 && is.null(cluster.call) && (robust==1)) {
     clusters <-rep(nobs+1,nobs)
     clusters[status==1] <- (1:nobs)[status==1]
     max.clust <- Ntimes+1
     gclusters <-  clusters <- as.integer(factor(clusters, labels = 1:max.clust)) -1
###   print(clusters) ; print(max.clust); print(c(max.clust,survs$antclust))
     survs$antclust <- max.clust
  }

  if ((length(beta)!=pz) && (is.null(beta)==FALSE)) beta <- rep(beta[1],pz); 
  if ((is.null(beta))) {
        if ( (attr(m[, 1], "type") == "right" ) ) 
        beta<-coxph(Surv(survs$stop,survs$status)~Z)$coef
        else { 
          if (survs$antpers< 10000) beta<-coxph(Surv(survs$start,survs$stop,survs$status)~Z)$coef 
	  else beta<-coxph(Surv(survs$stop,survs$status)~Z)$coef;  
	}
	if (detail>=2) {cat("starting values (coxph) \n");  print(beta);}
  }

if ( (attr(m[, 1], "type") == "right" ) ) {  ## {{{
   # order in time, status=0 first for ties
   ot<-order(-time2,status==1); 
   time2<-time2[ot]; 
   status<-status[ot]; 
   X<-as.matrix(X[ot,])
   if (npar==FALSE) Z<-as.matrix(Z[ot,])
   stop<-time2;
   clusters<-clusters[ot]
   id<-id[ot];
   weights <- weights[ot]
   offsets <- offsets[ot]
   entry=rep(-1,nobs); 
  if (!is.null(strata)) strata <- strata[ot]
  } else {
        eventtms <- c(survs$start,time2)
        status <- c(rep(0, nobs), status)
### order: strata , time, status=0 first for ties
###        if (!is.null(strata)) ix<-order(strata,-eventtms,status==1) else 
	ix <- order(-eventtms,status==1)
        etimes    <- eventtms[ix]  # Entry/exit times
	status <- status[ix]
        stop  <- etimes; 
        start <- rep(survs$start,2)[ix]; 
        tdiff    <- c(-diff(etimes),start.time) # Event time differences
        entry  <- c(rep(c(1, -1), each = nobs))[ix]
        weights <- rep(weights, 2)[ix]
        X  <- X[rep(1:nobs, 2)[ix],]
	if (npar==FALSE) Z <- Z[rep(1:nobs,2)[ix],]
	id <- rep(id,2)[ix]
	clusters <- rep(clusters,2)[ix]
        offsets <- rep(offsets,2)[ix]
        if (!is.null(strata)) strata <- rep(strata,2)[ix] 
} ## }}}

ldata<-list(start=start,stop=stop,antpers=survs$antpers,antclust=survs$antclust);
## }}}

###  if (npar==FALSE) covar<-data.matrix(cbind(X,Z)) else 
  if (npar==TRUE) stop("Both multiplicative and additive model needed");
  Ntimes <- sum(status); 

  if (px==0) stop("No nonparametric terms (needs one!)");

  ud<-cox.aalenBase(times,ldata,X,Z,
            status,id,clusters,Nit=Nit,detail=detail,beta=beta,weights=weights,
            sim=sim,antsim=n.sim,residuals=residuals,robust=robust,
            weighted.test=weighted.test,ratesim=rate.sim,
            covariance=covariance,resample.iid=resample.iid,namesX=covnamesX,
	    namesZ=covnamesZ,beta.fixed=beta.fixed,entry=entry,basesim=basesim,
	    offsets=offsets,exactderiv=exact.deriv,max.timepoint.sim=max.timepoint.sim,silent=silent,
	    strata=strata,propodds=propodds,caseweight=caseweight)

  ## {{{ output handling
  colnames(ud$test.procProp)<-c("time",covnamesZ)
  if (beta.fixed==1) colnames(ud$var.score)<-c("time",covnamesZ)
  if (robust==1 & beta.fixed==0) colnames(ud$var.score)<-c("time",covnamesZ)

  if (px>0) {
    colnames(ud$cum)<-colnames(ud$var.cum)<- c("time",covnamesX)
    if (robust==1) colnames(ud$robvar.cum)<- c("time",covnamesX)
    if (sim==1) {
      names(ud$pval.Prop)<- covnamesZ
    if (basesim[1]>0) {
      names(ud$conf.band)<- names(ud$pval.testBeq0)<- names(ud$pval.testBeqC)<-
      names(ud$obs.testBeq0)<- names(ud$obs.testBeqC)<- 
      colnames(ud$sim.testBeq0)<- covnamesX; 
    }
    } 
  }
  covariance<-ud$covariance

  rownames(ud$gamma)<-c(covnamesZ); colnames(ud$gamma)<-"estimate"; 
  rownames(ud$score)<-c(covnamesZ); colnames(ud$score)<-"score"; 
  namematrix(ud$var.gamma,covnamesZ); 
  namematrix(ud$robvar.gamma,covnamesZ); 
  if (beta.fixed==1) {ud$var.gamma<-matrix(0,pz,pz); 
                      ud$robvar.gamma<-matrix(0,pz,pz);
  }
  namematrix(ud$D2linv,covnamesZ); 
  ud$prop.odds <- propodds

  class(ud)<-"cox.aalen"
  attr(ud,"Call")<-call; 
  attr(ud,"stratum")<-ud$stratum; 
  attr(ud,"Formula")<-formula;
  attr(ud,"formula")<-formula;
  attr(ud,"rate.sim")<-rate.sim;
  attr(ud,"id.call")<-id.call;
  attr(ud,"id")<-id.call;
  attr(ud,"cluster.call")<-cluster.call;
  attr(ud,"cluster")<-gclusters; 
  attr(ud,"time2")<-time2; 
  attr(ud,"start.time")<-start.time; 
  attr(ud,"start")<-start.call; 
  attr(ud,"stop")<-stop.call; 
  attr(ud,"weights")<-weights.call; 
  attr(ud,"offsets")<-offsets.call; 
  attr(ud,"propodds")<-propodds 
  attr(ud,"type")<-"survival"
  attr(ud,"beta.fixed")<-beta.fixed
  attr(ud,"status")<-status.call; 
  attr(ud,"residuals")<-residuals; 
  attr(ud,"max.clust")<-max.clust; 
  attr(ud,"max.time")<-max.time; 
  attr(ud,"n")<-ldata$antpers; 
  attr(ud,"orig.max.clust")<- orig.max.clust 
  attr(ud,"max.timepoint.sim")<-max.timepoint.sim; 
  ud$call<-call
  ## }}}

  return(ud); 
} ## }}}

"plot.cox.aalen" <-  function (x,pointwise.ci=1, hw.ci=0,
sim.ci=0, robust.ci=0, col=NULL, specific.comps=FALSE,level=0.05, start.time = 0,
stop.time = 0, add.to.plot=FALSE,main=NULL,mains=TRUE,xlab="Time",score=FALSE,
ylab="Cumulative coefficients",...)
{ ## {{{
  object <- x; rm(x);  
  if (!inherits(object,'cox.aalen') ) stop ("Must be output from Cox-Aalen function")
  if (ylab=="Cumulative coefficients" && (1*score)>=1) ylab <- "Cumulative MG-residuals"

  if (score==FALSE) plot.cums(object, pointwise.ci=pointwise.ci,
        hw.ci=hw.ci, sim.ci=sim.ci, robust.ci=robust.ci, col=col, specific.comps=specific.comps,level=level,
        start.time = start.time, stop.time = stop.time, add.to.plot=add.to.plot,
	main=main, mains=mains, xlab=xlab,ylab=ylab,...)
  else plotScore(object, specific.comps=specific.comps, mains=mains,
		 main=main,xlab=xlab,ylab=ylab,...);
} ## }}}

"print.cox.aalen" <- function (x,...) 
{ ## {{{
summary.cox.aalen(x,...) 
###  cox.aalen.object <- x; rm(x);
###  if (!inherits(cox.aalen.object, 'cox.aalen')) 
###    stop ("Must be an aalen object")
###if (is.null(cox.aalen.object$prop.odds)==TRUE) p.o<-FALSE else p.o<-TRUE
###if (is.null(cox.aalen.object$gamma)==TRUE) prop<-FALSE else prop<-TRUE
###    
###  # We print information about object:  
###  cat("Cox-Aalen Model \n\n")
###  cat("Additive Aalen terms : "); 
###  cat(colnames(cox.aalen.object$cum)[-1]); cat("   \n");  
###  if (prop) {
###  cat("Proportional Cox terms :  "); 
###  cat(rownames(cox.aalen.object$gamma)); 
###  cat("   \n");  }
###  cat("   \n");  
###
###  cat("  Call: \n")
###  dput(attr(cox.aalen.object,"Call"))
###  cat("\n")
} ## }}}

"summary.cox.aalen" <- function (object,digits = 3,...) 
{ ## {{{
  cox.aalen.object <- object; rm(object);
  obj<-cox.aalen.object
  if (!inherits(cox.aalen.object, 'cox.aalen')) 
    stop ("Must be a Cox-Aalen object")
  
  prop<-TRUE; 
  if (is.null(cox.aalen.object$gamma)==TRUE) stop(" No regression terms"); 
  if (cox.aalen.object$prop.odds==0) p.o<-FALSE else p.o<-TRUE
    
  if (p.o==FALSE) cat("Cox-Aalen Model \n\n") else cat("Proportional Odds model \n\n")

  if (sum(abs(cox.aalen.object$score)>0.000001)) 
    cat("Did not converge, allow more iterations\n\n"); 

  if (p.o==FALSE) cat("Test for Aalen terms \n") else  cat("Test for baseline \n")

  if (is.null(obj$conf.band)==TRUE)  mtest<-FALSE else mtest<-TRUE; 
  if (mtest==FALSE) cat("Test not computed, sim=0 \n\n")
  if (mtest==TRUE) { 

    timetest(obj,digits=digits)
  }

  if (prop) {
    if (p.o==FALSE) cat("Proportional Cox terms :  \n") else  cat("Covariate effects \n")

    out=coef.cox.aalen(obj); 
    out=signif(out,digits=digits)
    print(out)

    if (p.o==FALSE) cat("Test of Proportionality \n") else  cat("Test of Goodness-of-fit \n")
    if (is.null(obj$pval.Prop)==TRUE)  ptest<-FALSE else ptest<-TRUE; 
    if (ptest==FALSE) cat("Test not computed, sim=0 \n\n")
    if (ptest==TRUE) { 
      testP<-cbind(
                   apply(abs(as.matrix(obj$test.procProp[,-1])),2,max),
                   obj$pval.Prop)
      testP<-as.matrix(testP); 
      colnames(testP) <- c("sup|  hat U(t) |","p-value H_0 ")
      prmatrix(signif(testP,digits)); cat("\n"); }  }

###  cat("   \n");  
###  cat("  Call: \n")
###  dput(attr(obj, "Call"))
###  cat("\n")
} ## }}}

coef.cox.aalen<-function(object,digits=3,d2logl=1,...) {
   coefBase(object,digits=digits, d2logl=d2logl,...)
}

vcov.cox.aalen <- function(object,robust=0, ...) {
  if (robust==0) rv <- object$var.gamma else  rv <- object$robvar.gamma
  if (!identical(rv, matrix(0, nrow = 1L, ncol = 1L))) rv # else return NULL
}


