.packageName <- "timereg"

.First.lib <- function(lib, pkg) {
  library.dynam("timereg", pkg, lib)
  cat("aalen.test    : semiparametric additive model with
improved test for constant effects and offsets \n")
  cat("aalen         : semiparametric additive model\n")
  cat("timecox       : semiparametric multiplicative model\n")
  cat("cox.aalen     : Cox-Aalen model\n")
  cat("prop.excess   : proportional excess model\n")
  cat("prop.odds     : semiparametric proportional odds-model\n")
  cat("Gprop.odds    : generalized semiparametric proportional odds-model\n")
  cat("pe.sasieni    : The proportional excess model by Peter Sasieni\n")
  cat("dynreg        : longitudinal data model\n")
  cat("cum.residuals : cumulative residuals for goodness-of-fit\n")
  cat("two.stage     : Clayton-Oakes-Glidden Two-Stage estimation\n")
  cat("comp.risk     : For flexible regression modellling for competing risks data\n"); 
  cat(" plus more \n"); 
}

.Last.lib <- function(lib){
  library.dynam.unload("timereg",lib)
}

aalen<-function (formula = formula(data),
                 data = sys.parent(), start.time = 0, max.time = NULL, 
                 robust=1, id=NULL, clusters=NULL, residuals = 0, n.sim = 1000, 
                 weighted.test= 0,covariance=0,resample.iid=0,deltaweight=1,
                 silent=0){
###deltaweight<-1; # always default
  if (n.sim == 0) sim <- 0 else sim <- 1
  if (resample.iid==1 & robust==0) {
    cat("When robust=0 no iid representaion computed\n"); resample.iid<-0;}
  if (covariance==1 & robust==0) {
    cat("When robust=0 no covariance computed \n"); 
    cat("Covariance based on robust iid representation\n")
    covariance<-0;}
  if (sim==1 & robust==0) {
    cat("When robust=0, No simulations \n"); cat("n.sim set to 0\n"); n.sim<-0;}
#if (residuals==1 & robust==0) {
#cat("When robust=0, no martingale residuals \n"); 
#residuals<-0;}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  call <- match.call()
  m <- match.call(expand = FALSE)
  m$start.time <- m$weighted.test <- m$max.time <- m$robust <- 
  m$sim <- m$residuals <- m$n.sim <- m$id <- m$covariance <- 
  m$resample.iid <- m$clusters <- m$deltaweight <- m$silent<- NULL
  special <- c("const","cluster") ###########
  Terms <- if (missing(data)){
    terms(formula, special)
  } else {
    terms(formula, special, data = data)
  }
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ

  scale<-0;
  if (scale==1) { X<-scale(as.matrix(X)); Z<-scale(as.matrix(Z)); 
     sXscale<- attr(X,"scaled:scale"); sXcenter<- attr(X,"scaled:center") 
     sZscale<- attr(Z,"scaled:scale"); sZcenter<- attr(Z,"scaled:center") 
     X[,attr(X,"scaled:scale")==0]<-1; Z[,attr(Z,"scaled:scale")==0]<-1; 
  }; 

  if(is.null(clusters)) clusters <- des$clusters ##########
  
  pxz <- px + pz; 

  survs<-read.surv(m,id,npar,clusters,start.time,max.time)
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; time2<-survs$stop
  status<-survs$status; 
  ldata<-list(start=survs$start,stop=survs$stop,
              antpers=survs$antpers,antclust=survs$antclust);

  if (npar== TRUE) { #cat("Nonparametric Additive Risk Model\n")
    ud <- aalenBase(times, ldata, X, status, id, clusters, robust = robust, 
                    sim = sim, retur = residuals, antsim = n.sim,
                    weighted.test = weighted.test,covariance=covariance,
                    resample.iid=resample.iid,namesX=covnamesX,silent=silent,scale=scale)
    colnames(ud$cum) <- colnames(ud$var.cum) <- c("time", covnamesX)
    if (robust == 1) colnames(ud$robvar.cum) <- c("time", covnamesX)
    if (sim >= 1) {
      colnames(ud$test.procBeqC) <- c("time", covnamesX)
      names(ud$conf.band) <- names(ud$pval.testBeq0) <- names(ud$pval.testBeqC) <- names(ud$pval.testBeqC.is) <- names(ud$obs.testBeq0) <- names(ud$obs.testBeqC) <- names(ud$obs.testBeqC.is) <- colnames(ud$sim.testBeq0) <- colnames(ud$sim.testBeqC) <- colnames(ud$sim.testBeqC.is) <- covnamesX
      ud$sim.testBeqC.is <- ud$sim.testBeqC <- FALSE
    }

  }
  else { #cat("Semiparametric Additive Risk Model\n")
    if (px == 0) 
      stop("No nonparametric terms (needs one!)")
    ud <- semiaalen(times, ldata, X, Z, 
                    status, id , clusters, robust = robust, sim = sim, antsim = n.sim, 
                    weighted.test = weighted.test, retur =
                    residuals,covariance=covariance,
                    resample.iid=resample.iid,namesX=covnamesX,namesZ=covnamesZ,
                    deltaweight=deltaweight,silent=silent,scale=scale)

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
  }
  attr(ud, "Call") <- sys.call()
  attr(ud, "Formula") <- formula
  attr(ud, "id") <- id.call
  attr(ud, "cluster") <- cluster.call
  attr(ud, "start") <- start.time
  attr(ud, "time2") <- time2
  class(ud) <- "aalen"
  ud$call<-call
  return(ud)
}

plot.aalen <-  function (x, pointwise.ci=1, hw.ci=0,
sim.ci=0, robust=0, specific.comps=FALSE,level=0.05, start.time = 0, 
stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
ylab ="Cumulative coefficients",score=FALSE,...) 
{
  object <- x; rm(x);
  if (!inherits(object,'aalen') ) 
    stop ("Must be output from Aalen function") 
 
  if (score==FALSE) plot.cums(object, pointwise.ci=pointwise.ci, 
        hw.ci=hw.ci,
        sim.ci=sim.ci, robust=robust, specific.comps=specific.comps,level=level, 
        start.time = start.time, stop.time = stop.time, add.to.plot=add.to.plot, 
        mains=mains, xlab=xlab, ylab =ylab) 
  else plot.score(object, specific.comps=specific.comps, mains=mains,
                  xlab=xlab,ylab =ylab); 
}

"print.aalen" <- function (x,...) 
{
  object <- x; rm(x);
  if (!inherits(object, 'aalen')) stop ("Must be an aalen object")

  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
                                        # We print information about object:
  
  cat("Additive Aalen Model \n\n")
  cat(" Nonparametric terms : "); cat(colnames(object$cum)[-1]); cat("   \n");  
  if (semi) {
    cat(" Parametric terms :  "); cat(rownames(object$gamma)); 
    cat("   \n");  } 
  cat("   \n");  

  cat("  Call: \n"); dput(attr(object, "Call")); cat("\n"); 
}


"summary.aalen" <-
function (object,digits = 3,...) 
{
  aalen.object <- object; rm(object);
  
  obj<-aalen.object
  if (!inherits(aalen.object, 'aalen')) stop ("Must be an aalen object")
  
  if (is.null(aalen.object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Additive Aalen Model \n\n")
  #cat(paste("loglikelihood : ",round(aalen.object$deviance,3),"\n"))
  #cat("Nonparametric terms : "); cat(colnames(aalen.object$cum)[-1]);
  #cat("   \n");  

  timetest(obj,digits=digits); 

  if (semi) {
    cat("Parametric terms :  ");  #cat(rownames(aalen.object$gamma)); 
  }
  cat("   \n");  

  if (semi) {
    coef.aalen(aalen.object,digits=digits); 

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
}

coef.aalen <- function(object, digits=3,...) {
  
   coefBase(object,digits=digits)
}
