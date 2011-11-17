prop<-function(x) x

cox.aalen<-function(formula=formula(data),data=sys.parent(),
beta=NULL,Nit=10,detail=0,start.time=0,max.time=NULL, id=NULL, 
clusters=NULL, n.sim=500, residuals=0,robust=1,
weighted.test=0,covariance=0,resample.iid=0,weights=NULL,
rate.sim=1,beta.fixed=0,max.clust=1000,exact.deriv=1,silent=0,
max.timepoint.sim=100)
{ ## {{{
offsets=0; 
## {{{ set up variables 
  if (n.sim == 0) sim <- 0 else sim <- 1
  if (resample.iid==1 & robust==0) { robust <- 1;}
  if (covariance==1 & robust==0) { covariance<-0;}
  if (sim==1 & robust==0) { n.sim <- 0; sim<-0;}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  if (beta.fixed==1) Nit<-1; 
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  m$robust<-m$start.time<- m$scaleLWY<-m$weighted.test<-m$beta<-m$Nit<-m$detail<-m$max.time<-m$residuals<-m$n.sim<-m$id<-
	                   m$covariance<-m$resample.iid<-m$clusters<-m$rate.sim<-m$beta.fixed<-
                           m$max.clust <- m$exact.deriv <- m$silent <- m$max.timepoint.sim <- NULL
  special <- c("prop","cluster")
  Terms <- if(missing(data)) terms(formula, special)
  else          terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms,model="cox.aalen")
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
  if(is.null(clusters)) clusters <- des$clusters  
  pxz <- px + pz;

  survs<-read.surv(m,id,npar,clusters,start.time,max.time,model="cox.aalen",silent=silent)
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; 
  start.call <- start <-  survs$start; 
  stop.call <- time2 <- survs$stop; 
  status<-survs$status;
  orig.max.clust <- survs$antclust
  nobs <- nrow(X); 
  if (is.null(weights)) weights <- rep(1,nrow(X));  
  weights <- rep(1,nrow(X)); 
  if (sum(abs(offsets))!=0) stop("no offsets in this version \n"); 

  if ( (!is.null(max.clust)) )  {  
     if (max.clust < survs$antclust)   {
       qq <- quantile(clusters, probs = seq(0, 1, by = 1/max.clust))       
       qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
       clusters <- as.integer(factor(qqc, labels = 1:max.clust)) -1
       survs$antclust <- max.clust    
     }
  }                                                         
  cluster.call<-clusters; 

  if ((length(beta)!=pz) && (is.null(beta)==FALSE)) beta <- rep(beta[1],pz); 
  if ((is.null(beta))) beta<-coxph(Surv(survs$start,survs$stop,survs$status)~Z)$coef; 

if ( (attr(m[, 1], "type") == "right" ) ) {  ## {{{
   ot<-order(-time2,status==1); # order in time, status=0 first for ties
   time2<-time2[ot]; status<-status[ot]; 
   X<-as.matrix(X[ot,])
   if (npar==FALSE) Z<-as.matrix(Z[ot,])
   stop<-time2;
   clusters<-clusters[ot]
   id<-id[ot];
   weights <- weights[ot]
   if (sum(abs(offsets))!=0) offsets <- offsets[ot]
   entry=rep(-1,nobs); 
  } else {
        eventtms <- c(survs$start,time2)
        status <- c(rep(0, nobs), status)
        ix <- order(-eventtms,status==1)
        etimes    <- eventtms[ix]  # Entry/exit times
	status <- status[ix]
        stop  <- etimes; 
        start <- c(survs$start,survs$start)[ix]; 
        tdiff    <- c(-diff(etimes),start.time) # Event time differences
        entry  <- c(rep(c(1, -1), each = nobs))[ix]
        weights <- rep(weights, 2)[ix]
        X        <- X[rep(1:nobs, 2)[ix],]
	if (npar==FALSE) Z <- Z[rep(1:nobs,2)[ix],]
	id <- rep(id,2)[ix]
	clusters <- rep(clusters,2)[ix]
	if (sum(abs(offsets))!=0) offsets <- rep(offsets,2)[ix]
    } ## }}}

ldata<-list(start=start,stop=stop, antpers=survs$antpers,antclust=survs$antclust);
## }}}

  if (npar==FALSE) covar<-data.matrix(cbind(X,Z)) else 
  stop("Both multiplicative and additive model needed");
  Ntimes <- sum(status); 

  if (px==0) stop("No nonparametric terms (needs one!)");
  ud<-cox.aalenBase(times,ldata,X,Z,
            status,id,clusters,Nit=Nit,detail=detail,beta=beta,weights=weights,
            sim=sim,antsim=n.sim,residuals=residuals,robust=robust,
            weighted.test=weighted.test,ratesim=rate.sim,
            covariance=covariance,resample.iid=resample.iid,namesX=covnamesX,
	    namesZ=covnamesZ,beta.fixed=beta.fixed,entry=entry,
	    offsets=0,exactderiv=exact.deriv,max.timepoint.sim=max.timepoint.sim);

  colnames(ud$test.procProp)<-c("time",covnamesZ)
  if (beta.fixed==1) colnames(ud$var.score)<-c("time",covnamesZ)
  if (robust==1 & beta.fixed==0) colnames(ud$var.score)<-c("time",covnamesZ)

  if (px>0) {
    colnames(ud$cum)<-colnames(ud$var.cum)<- c("time",covnamesX)
    if (robust==1) colnames(ud$robvar.cum)<- c("time",covnamesX)
    if (sim==1) {
      names(ud$pval.Prop)<- covnamesZ
      names(ud$conf.band)<- names(ud$pval.testBeq0)<-
        names(ud$pval.testBeqC)<-
          names(ud$obs.testBeq0)<- names(ud$obs.testBeqC)<-
            colnames(ud$sim.testBeq0)<- covnamesX; 
    } }
  covariance<-ud$covariance

  rownames(ud$gamma)<-c(covnamesZ); colnames(ud$gamma)<-"estimate"; 
  rownames(ud$score)<-c(covnamesZ); colnames(ud$score)<-"score"; 
  namematrix(ud$var.gamma,covnamesZ); 
  namematrix(ud$robvar.gamma,covnamesZ); 
  if (beta.fixed==1) {ud$var.gamma<-matrix(0,pz,pz); 
                      ud$robvar.gamma<-matrix(0,pz,pz);
  }
  namematrix(ud$D2linv,covnamesZ); 

  attr(ud,"Call")<-sys.call(); 
  class(ud)<-"cox.aalen"
  attr(ud,"Formula")<-formula;
  attr(ud,"id")<-id.call;
  attr(ud,"cluster")<-cluster.call;
  attr(ud,"time2")<-time2; 
  attr(ud,"start.time")<-start.time; 
  attr(ud,"start")<-start.call; 
  attr(ud,"stop")<-stop.call; 
  attr(ud,"beta.fixed")<-beta.fixed
  attr(ud,"status")<-survs$status; 
  attr(ud,"residuals")<-residuals; 
  attr(ud,"max.clust")<-max.clust; 
  attr(ud,"orig.max.clust")<- orig.max.clust 
  attr(ud,"max.timepoint.sim")<-max.timepoint.sim; 
  ud$call<-call

  return(ud); 
} ## }}}

"plot.cox.aalen" <-  function (x,pointwise.ci=1, hw.ci=0,
sim.ci=0, robust=0, specific.comps=FALSE,level=0.05, start.time = 0,
stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
ylab ="Cumulative coefficients",score=FALSE,...)
{ ## {{{
  object <- x; rm(x);  
  if (!inherits(object,'cox.aalen') ) stop ("Must be output from Cox-Aalen function")

  if (score==FALSE) plot.cums(object, pointwise.ci=pointwise.ci,
        hw.ci=hw.ci,
        sim.ci=sim.ci, robust=robust, specific.comps=specific.comps,level=level,
        start.time = start.time, stop.time = stop.time, add.to.plot=add.to.plot,
        mains=mains, xlab=xlab, ylab =ylab)
  else plotScore(object, specific.comps=specific.comps, mains=mains,
                  xlab=xlab,ylab =ylab);
} ## }}}

"print.cox.aalen" <- function (x,...) 
{ ## {{{
  cox.aalen.object <- x; rm(x);
  if (!inherits(cox.aalen.object, 'cox.aalen')) 
    stop ("Must be an aalen object")

if (is.null(cox.aalen.object$prop.odds)==TRUE) p.o<-FALSE else p.o<-TRUE

if (is.null(cox.aalen.object$gamma)==TRUE) prop<-FALSE else prop<-TRUE
    
  # We print information about object:  
  cat("Cox-Aalen Model \n\n")
  cat("Additive Aalen terms : "); 
  cat(colnames(cox.aalen.object$cum)[-1]); cat("   \n");  
  if (prop) {
  cat("Proportional Cox terms :  "); 
  cat(rownames(cox.aalen.object$gamma)); 
  cat("   \n");  }
  cat("   \n");  

  cat("  Call: \n")
  dput(attr(cox.aalen.object,"Call"))
  cat("\n")
} ## }}}

"summary.cox.aalen" <- function (object,digits = 3,...) 
{ ## {{{
  cox.aalen.object <- object; rm(object);
  obj<-cox.aalen.object
  if (!inherits(cox.aalen.object, 'cox.aalen')) 
    stop ("Must be a Cox-Aalen object")
  
  prop<-TRUE; 
  if (is.null(cox.aalen.object$gamma)==TRUE) stop(" No regression terms"); 
  if (is.null(cox.aalen.object$prop.odds)==TRUE) p.o<-FALSE else p.o<-TRUE
    
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

    if (p.o==FALSE) cat("Test for Proportionality \n") else  
    cat("Test for Goodness-of-fit \n")
    if (is.null(obj$pval.Prop)==TRUE)  ptest<-FALSE else ptest<-TRUE; 
    if (ptest==FALSE) cat("Test not computed, sim=0 \n\n")
    if (ptest==TRUE) { 
      testP<-cbind(
                   apply(abs(as.matrix(obj$test.procProp[,-1])),2,max),
                   obj$pval.Prop)
      testP<-as.matrix(testP); 
      colnames(testP) <- c("sup|  hat U(t) |","p-value H_0 ")
      prmatrix(signif(testP,digits)); cat("\n"); }  }

  cat("   \n");  
  cat("  Call: \n")
  dput(attr(obj, "Call"))
  cat("\n")
} ## }}}

coef.cox.aalen<-function(object,digits=3,d2logl=1,...) {
   coefBase(object,digits=digits, d2logl=d2logl,...)
}
