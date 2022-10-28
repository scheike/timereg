#' Fit additive hazards model
#' 
#' Fits both the additive hazards model of Aalen and the semi-parametric
#' additive hazards model of McKeague and Sasieni.  Estimates are un-weighted.
#' Time dependent variables and counting process data (multiple events per
#' subject) are possible.
#' 
#' Resampling is used for computing p-values for tests of time-varying effects.
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
#' operator, and the independent terms on the right as regressors.The response
#' must be a survival object as returned by the `Surv' function. Time-
#' invariant regressors are specified by the wrapper const(), and cluster
#' variables (for computing robust variances) by the wrapper cluster().
#' @param data a data.frame with the variables.
#' @param start.time start of observation period where estimates are computed.
#' @param max.time end of observation period where estimates are computed.
#' Estimates thus computed from [start.time, max.time]. Default is max of data.
#' @param robust to compute robust variances and construct processes for
#' resampling. May be set to 0 to save memory.
#' @param id For timevarying covariates the variable must associate each record
#' with the id of a subject.
#' @param clusters cluster variable for computation of robust variances.
#' @param n.sim number of simulations in resampling.
#' @param weighted.test to compute a variance weighted version of the
#' test-processes used for testing time-varying effects.
#' @param residuals to returns residuals that can be used for model validation
#' in the function cum.residuals
#' @param covariance to compute covariance estimates for nonparametric terms
#' rather than just the variances.
#' @param resample.iid to return i.i.d. representation for nonparametric and
#' parametric terms.
#' @param deltaweight uses weights to estimate semiparametric model, under
#' construction, default=1 is standard least squares estimates
#' @param silent set to 0 to print warnings for non-inverible design-matrices
#' for different timepoints, default is 1.
#' @param weights weights for estimating equations.
#' @param max.clust sets the total number of i.i.d. terms in i.i.d.
#' decompostition. This can limit the amount of memory used by coarsening the
#' clusters. When NULL then all clusters are used.  Default is 1000 to save
#' memory and time.
#' @param gamma fixes gamme at this value for estimation.
#' @param offsets offsets for the additive model, to make excess risk
#' modelling.
#' @param caseweight caseweight: mutiplied onto dN for score equations.
#' @return returns an object of type "aalen". With the following arguments:
#' \item{cum}{cumulative timevarying regression coefficient estimates are
#' computed within the estimation interval. } \item{var.cum}{the martingale
#' based pointwise variance estimates for cumulatives.}
#' \item{robvar.cum}{robust pointwise variances estimates for cumulatives.}
#' \item{gamma}{estimate of parametric components of model.  }
#' \item{var.gamma}{variance for gamma.  } \item{robvar.gamma}{robust variance
#' for gamma.  } \item{residuals}{list with residuals. Estimated martingale
#' increments (dM) and corresponding time vector (time).}
#' \item{obs.testBeq0}{observed absolute value of supremum of cumulative
#' components scaled with the variance.} \item{pval.testBeq0}{p-value for
#' covariate effects based on supremum test.} \item{sim.testBeq0}{resampled
#' supremum values.} \item{obs.testBeqC}{observed absolute value of supremum of
#' difference between observed cumulative process and estimate under null of
#' constant effect.} \item{pval.testBeqC}{p-value based on resampling.}
#' \item{sim.testBeqC}{resampled supremum values.}
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
#' \item{deviance}{Least squares of increments.}
#' @author Thomas Scheike
#' @references Martinussen and Scheike, Dynamic Regression Models for Survival
#' Data, Springer (2006).
#' @keywords survival
#' @examples
#' 
#' data(sTRACE)
#' # Fits Aalen model 
#' out<-aalen(Surv(time,status==9)~age+sex+diabetes+chf+vf,
#' sTRACE,max.time=7,n.sim=100)
#' 
#' summary(out)
#' par(mfrow=c(2,3))
#' plot(out)
#' 
#' # Fits semi-parametric additive hazards model 
#' out<-aalen(Surv(time,status==9)~const(age)+const(sex)+const(diabetes)+chf+vf,
#' sTRACE,max.time=7,n.sim=100)
#' 
#' summary(out)
#' par(mfrow=c(2,3))
#' plot(out)
#' 
#' ## Excess risk additive modelling 
#' data(mela.pop)
#' dummy<-rnorm(nrow(mela.pop));
#' 
#' # Fits Aalen model  with offsets 
#' out<-aalen(Surv(start,stop,status==1)~age+sex+const(dummy),
#' mela.pop,max.time=7,n.sim=100,offsets=mela.pop$rate,id=mela.pop$id,
#' gamma=0)
#' summary(out)
#' par(mfrow=c(2,3))
#' plot(out,main="Additive excess riks model")
#' 
#' # Fits semi-parametric additive hazards model  with offsets 
#' out<-aalen(Surv(start,stop,status==1)~age+const(sex),
#' mela.pop,max.time=7,n.sim=100,offsets=mela.pop$rate,id=mela.pop$id)
#' summary(out)
#' plot(out,main="Additive excess riks model")
#' 
##' @export
aalen<-function (formula = formula(data),
     data = parent.frame(), start.time = 0, max.time = NULL, 
     robust=1, id=NULL, clusters=NULL, residuals = 0, n.sim = 1000,  
     weighted.test= 0,covariance=0,resample.iid=0,
     deltaweight=1,silent=1,weights=NULL,max.clust=1000,
     gamma=NULL,offsets=0,caseweight=NULL){ ## {{{
## {{{ setting up variables 
  if (n.sim == 0) sim <- 0 else sim <- 1
  if (resample.iid==1 & robust==0) { robust <- 1;}
  if (covariance==1 & robust==0) { covariance<-0;}
  if (sim==1 & robust==0) { n.sim<-0;sim <- 0}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
    m$start.time <- m$weighted.test <- m$max.time <- m$robust <- 
    m$weights <- m$residuals <- m$n.sim <- m$id <- m$covariance <- 
    m$resample.iid <- m$clusters <- m$deltaweight<-m$silent <- 
    m$max.clust <- m$gamma <- m$offsets <- m$caseweight <- NULL
  special <- c("const","cluster") 
  Terms <- if (missing(data)){
    terms(formula, special)
  } else {
    terms(formula, special, data = data)
  }
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
  if(is.null(clusters)) clusters <- des$clusters ##########
  pxz <- px + pz; 

  if ((nrow(X)!=nrow(data) && (!is.null(id)))) stop("Missing values in design matrix not allowed with id \n"); 
###  if (nrow(Z)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 
  if (!is.null(id)) {
	  if (length(id)!=nrow(X)) stop("id length and data not the same\n"); 
  }

  cluster.call<- clusters; 
  survs<-read.surv(m,id,npar,clusters,start.time,max.time,silent=silent)
  times<-survs$times; 
  id<-survs$id.cal; 
  id.call<-id; 
  clusters<-gclusters <- survs$clusters; 
  stop.call <- time2<-survs$stop
  start.call <- survs$start
  status<-survs$status; 
  orig.max.clust <- survs$antclust
  dtimes <- sort(survs$stop[survs$status==1])
  nobs <- nrow(X); 
  if (is.null(weights)) weights <- rep(1,nrow(X)); 
  weights.call <- weights; 

  if ((!is.null(max.clust))) if (max.clust<survs$antclust) {
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	clusters <- gclusters <-  as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
###	clusters <- as.integer(factor(qqc, labels = 1:max.clust)) -1
	survs$antclust <- max.clust    
  }                                                         

  if ( (attr(m[, 1], "type") == "right" ) ) {  ## {{{
   ot <- order(-time2,status==1); 
   time2<-time2[ot]; 
   status<-status[ot]; 
   X<-as.matrix(X[ot,])
   if (npar==FALSE) Z<-as.matrix(Z[ot,])
   survs$stop<-time2;
   clusters<-clusters[ot]
   id<-id[ot];
   entry=rep(-1,nobs); 
   weights <- weights[ot]
   if (sum(offsets)!=0) offsets <- offsets[ot]
  } else {
        eventtms <- c(survs$start,time2)
        status <- c(rep(0, nobs), status)
        ix <- order(-eventtms,status==1)
        etimes    <- eventtms[ix]  # Entry/exit times
	status <- status[ix]
        survs$stop  <- etimes; 
        survs$start <- c(survs$start,survs$start)[ix]; 
        tdiff    <- c(-diff(etimes),start.time) # Event time differences
        entry  <- c(rep(c(1, -1), each = nobs))[ix]
	weights <- rep(weights, 2)[ix]
	X        <- X[rep(1:nobs, 2)[ix],]
	if (npar==FALSE) Z <- Z[rep(1:nobs,2)[ix],]
	id <- rep(id,2)[ix]
	clusters <- rep(clusters,2)[ix]
	if (sum(offsets)!=0) offsets <- rep(offsets,2)[ix]
} ## }}}

ldata<-list(start=survs$start,stop=survs$stop,antpers=survs$antpers,antclust=survs$antclust);

## }}}

if (npar==FALSE & residuals==1)  ## warning for residuals 
   cat("Residuals only computed for non-parametric aalen model\n");

  if (npar== TRUE) { ## {{{ Aalen model 
   ud <- aalenBase(times, ldata, X, status, id, clusters, robust = robust, 
   sim = sim, retur = residuals, antsim = n.sim,
   weighted.test = weighted.test,covariance=covariance,
   resample.iid=resample.iid,namesX=covnamesX,
   silent=silent,weights=weights,entry=entry,offsets=offsets,caseweight=caseweight)

    colnames(ud$cum) <- colnames(ud$var.cum) <- c("time",covnamesX)
    if (robust == 1) colnames(ud$robvar.cum) <- c("time", covnamesX)
    if (sim >= 1) {
      colnames(ud$test.procBeqC) <- c("time", covnamesX)
      names(ud$conf.band) <- names(ud$pval.testBeq0) <- names(ud$pval.testBeqC) <- names(ud$pval.testBeqC.is) <- names(ud$obs.testBeq0) <- names(ud$obs.testBeqC) <- names(ud$obs.testBeqC.is) <- colnames(ud$sim.testBeq0) <- colnames(ud$sim.testBeqC) <- colnames(ud$sim.testBeqC.is) <- covnamesX
      ud$sim.testBeqC.is <- ud$sim.testBeqC <- FALSE
    }
  } ## }}}
  else { ## {{{ Semiparametric additive risk model 
    if (px == 0) stop("No nonparametric terms (needs one!)")
    ud<-semiaalen(times, ldata, X, Z, 
    status, id , clusters, robust = robust, sim = sim, antsim = n.sim, 
    weighted.test = weighted.test, retur =
	residuals,covariance=covariance,
	resample.iid=resample.iid,namesX=covnamesX,namesZ=covnamesZ,
	deltaweight=deltaweight,gamma=gamma,
	silent=silent,weights=weights,entry=entry,offsets=offsets,caseweight=caseweight)

    if (px > 0) {
      colnames(ud$cum) <- colnames(ud$var.cum) <- c("time", covnamesX)
      if (robust == 1) 
        colnames(ud$robvar.cum) <- c("time", covnamesX)
      if (sim >= 1) {
        colnames(ud$test.procBeqC) <- c("time", covnamesX)
        names(ud$conf.band) <- names(ud$pval.testBeq0) <- names(ud$pval.testBeqC) <- names(ud$pval.testBeqC.is) <- names(ud$obs.testBeqC.is) <- names(ud$obs.testBeq0) <- names(ud$obs.testBeqC) <- colnames(ud$sim.testBeq0) <- colnames(ud$sim.testBeqC.is) <- colnames(ud$sim.testBeqC) <- covnamesX
        ud$sim.testBeqC.is <- ud$sim.testBeqC <- FALSE
      }
    }
    ud$gamma<-as.matrix(ud$gamma);
    rownames(ud$gamma) <- c(covnamesZ)
    rownames(ud$intZHdN) <- c(covnamesZ)
    colnames(ud$gamma) <- "estimate"
    colnames(ud$var.gamma) <- c(covnamesZ)
    rownames(ud$var.gamma) <- c(covnamesZ)
    colnames(ud$robvar.gamma) <- c(covnamesZ)
    colnames(ud$intZHZ) <- c(covnamesZ)
    rownames(ud$var.gamma) <- c(covnamesZ)
  } ## }}}

  attr(ud,"stratum")<-ud$stratum; 
  attr(ud, "Call") <- call
  attr(ud, "Formula") <- formula
  attr(ud, "id") <- id.call
  attr(ud, "cluster.call") <- cluster.call
  attr(ud, "cluster") <- gclusters
  attr(ud, "start.time") <- start.time
  attr(ud, "stop") <- stop.call
  attr(ud, "start") <- start.call
  attr(ud, "status") <- survs$status
  attr(ud, "residuals") <- residuals
  attr(ud, "max.clust") <- max.clust; 
  attr(ud, "max.time") <- max.time; 
  attr(ud, "weights") <- weights.call; 
  attr(ud, "orig.max.clust") <- orig.max.clust 
  class(ud) <- "aalen"
  ud$call<-call
  return(ud)
} ## }}}



#' Plots estimates and test-processes
#' 
#' This function plots the non-parametric cumulative estimates for the additive
#' risk model or the test-processes for the hypothesis of time-varying effects
#' with re-sampled processes under the null.
#' 
#' 
#' @aliases plot.aalen plot.cox.aalen plot.timecox plot.prop.excess
#' @param x the output from the "aalen" function.
#' @param pointwise.ci if >1 pointwise confidence intervals are plotted with
#' lty=pointwise.ci
#' @param hw.ci if >1 Hall-Wellner confidence bands are plotted with lty=hw.ci.
#' Only 0.95 \% bands can be constructed.
#' @param sim.ci if >1 simulation based confidence bands are plotted with
#' lty=sim.ci. These confidence bands are robust to non-martingale behaviour.
#' @param robust.ci robust standard errors are used to estimate standard error
#' of estimate, otherwise martingale based standard errors are used.
#' @param col specifice colors of different components of plot, in order:
#' c(estimate,pointwise.ci,robust.ci,hw.ci,sim.ci) so for example, when we ask
#' to get pointwise.ci, hw.ci and sim.ci we would say c(1,2,3,4) to use colors
#' as specified.
#' @param specific.comps all components of the model is plotted by default, but
#' a list of components may be specified, for example first and third "c(1,3)".
#' @param level gives the significance level.
#' @param start.time start of observation period where estimates are plotted.
#' @param stop.time end of period where estimates are plotted. Estimates thus
#' plotted from [start.time, max.time].
#' @param add.to.plot to add to an already existing plot.
#' @param mains add names of covariates as titles to plots.
#' @param xlab label for x-axis.
#' @param ylab label for y-axis.
#' @param score to plot test processes for test of time-varying effects along
#' with 50 random realization under the null-hypothesis.
#' @param ... unused arguments - for S3 compatibility
#' @author Thomas Scheike
#' @references Martinussen and Scheike, Dynamic Regression models for Survival
#' Data, Springer (2006).
#' @keywords survival
#' @examples
#' 
#' # see help(aalen) 
#' data(sTRACE)
#' out<-aalen(Surv(time,status==9)~chf+vf,sTRACE,max.time=7,n.sim=100)
#' par(mfrow=c(2,2))
#' plot(out,pointwise.ci=1,hw.ci=1,sim.ci=1,col=c(1,2,3,4))
#' par(mfrow=c(2,2))
#' plot(out,pointwise.ci=0,robust.ci=1,hw.ci=1,sim.ci=1,col=c(1,2,3,4))
#' 
##' @export
plot.aalen <-  function (x, pointwise.ci=1, hw.ci=0,
sim.ci=0, robust.ci=0, col=NULL,
specific.comps=FALSE,level=0.05, start.time = 0, 
stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
ylab ="Cumulative coefficients",score=FALSE,...) 
{ ## {{{
  object <- x; rm(x);
  if (!inherits(object,'aalen') ) 
    stop ("Must be output from Aalen function") 
 
  if (score==FALSE) plot.cums(object, pointwise.ci=pointwise.ci, 
        hw.ci=hw.ci, sim.ci=sim.ci, robust.ci=robust.ci, col=col, 
	specific.comps=specific.comps,level=level, 
        start.time = start.time, stop.time = stop.time, add.to.plot=add.to.plot, 
        mains=mains, xlab=xlab, ylab =ylab,...) 
  else plotScore(object, specific.comps=specific.comps, mains=mains,
                  xlab=xlab,ylab =ylab,...); 
} ## }}}


##' @export
vcov.aalen <- function(object,robust=0, ...) {
  if (robust==0) rv <- object$var.gamma else  rv <- object$robvar.gamma
  if (!identical(rv, matrix(0, nrow = 1L, ncol = 1L))) rv # else return NULL
}


#' Prints call
#' 
#' Prints call for object. Lists nonparametric and parametric terms of model
#' 
#' 
#' @aliases print.aalen print.cox.aalen print.comprisk print.prop.excess
#' print.dynreg print.timecox print.cum.residuals
#' @param x an aalen object
#' @param ... unused arguments - for S3 compatibility
#' @author Thomas Scheike
#' @keywords survival
##' @export
"print.aalen" <- function (x,...) 
{ ## {{{
summary.aalen(x,...)
###  object <- x; rm(x);
###  if (!inherits(object, 'aalen')) stop ("Must be an aalen object")
###
###  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
###    
###                                        # We print information about object:
###  
###  cat("Additive Aalen Model \n\n")
###  cat(" Nonparametric terms : "); cat(colnames(object$cum)[-1]); cat("   \n");  
###  if (semi) {
###    cat(" Parametric terms :  "); cat(rownames(object$gamma)); 
###    cat("   \n");  } 
###  cat("   \n");  
###
###  cat("  Call: \n"); dput(attr(object, "Call")); cat("\n"); 
} ## }}}


#' Prints summary statistics
#' 
#' Computes p-values for test of significance for nonparametric terms of model,
#' p-values for test of constant effects based on both supremum and integrated
#' squared difference.
#' 
#' Returns parameter estimates and their standard errors.
#' 
#' 
#' @aliases summary.aalen summary.cox.aalen summary.prop.excess summary.timecox
#' summary.dynreg
#' @param object an aalen object.
#' @param digits number of digits in printouts.
#' @param ... unused arguments - for S3 compatibility
#' @author Thomas Scheike
#' @references Martinussen and Scheike,
#' @keywords survival
#' @examples
#' 
#' ### see help(aalen)
#' 
##' @export
"summary.aalen" <-
function (object,digits = 3,...) 
{ ## {{{
  aalen.object <- object; rm(object);
  
  obj<-aalen.object
  if (!inherits(aalen.object, 'aalen')) stop ("Must be an aalen object")
  
  if (is.null(aalen.object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Additive Aalen Model \n\n")
 #cat("Nonparametric terms : "); cat(colnames(aalen.object$cum)[-1]);
 #cat("   \n");  

  timetest(obj,digits=digits); 

  if (semi) {
    cat("Parametric terms :  ");  #cat(rownames(aalen.object$gamma)); 
  }
  cat("   \n");  

  if (semi) {
    out=coef.aalen(aalen.object,digits=digits); 
    out=signif(out,digits=digits)
    print(out)

    if (is.null(aalen.object$pstest.pval)==FALSE) {
      res<-cbind(aalen.object$sup.pscore,aalen.object$pstest.pval);
      colnames(res)<-c("sup of pseudo-score test","p-value H_0: B(t)=b t");  
      cat(" \n");  
      cat("Test for time invariant effects \n")
      prmatrix(signif(res,digits))
    }
  }
  cat("   \n");  

  cat("  Call: \n")
  dput(attr(aalen.object, "Call"))
  cat("\n")
} ## }}}

##' @export
coef.aalen <- function(object, digits=3,...) {
   coefBase(object,digits=digits,...)
}
