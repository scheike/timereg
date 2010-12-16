prop<-function(x) x

cox.aalen<-function(formula=formula(data),data=sys.parent(),
beta=0,Nit=10,detail=0,start.time=0,max.time=NULL, id=NULL, 
clusters=NULL, n.sim=500, residuals=0,robust=1,
weighted.test=0,covariance=0,resample.iid=0,weights=NULL,
rate.sim=1,beta.fixed=0,max.clust=1000)
{
  ridge<-0; ratesim<-rate.sim; if (n.sim==0) sim<-0 else sim<-1; 
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  m$robust<-m$start.time<-m$scaleLWY<-m$weighted.test<-m$beta<-m$Nit<-m$detail<-m$max.time<-m$residuals<-m$n.sim<-m$id<-m$covariance<-m$resample.iid<-m$clusters<-m$rate.sim<-m$beta.fixed<-m$max.clust <- NULL

  if (resample.iid==1 & robust==0) {
    cat("When robust=0 no iid representaion computed\n"); 
    resample.iid<-0;}
  if (covariance==1 & robust==0) {
    cat("When robust=0 no covariance computed \n"); 
    cat("Covariance based on robust iid representation\n")
    covariance<-0;}
  if (sim==1 & robust==0) {
    cat("When robust=0, No simulations \n"); cat("n.sim set to 0\n"); n.sim <- 0; sim<-0;}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  if (beta.fixed==1) Nit<-1; 

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

  survs<-read.surv(m,id,npar,clusters,start.time,max.time,model="cox.aalen")
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; 
  time<-survs$start; time2<-survs$stop; status<-survs$status;
  ldata<-list(start=survs$start,stop=survs$stop,
              antpers=survs$antpers,antclust=survs$antclust);

      if ( (!is.null(max.clust)) )  {  
       if (max.clust < survs$antclust) {
       qq <- quantile(clusters, probs = seq(0, 1, by = 1/max.clust))       
       qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
       clusters <- as.integer(factor(qqc, labels = 1:max.clust)) -1
       survs$antclust <- max.clust    
       ldata$antclust <- max.clust    
       cluster.call <- clusters; 
       }
      }                                                         

  if (npar==FALSE) covar<-data.matrix(cbind(X,Z)) else 
  stop("Both multiplicative and additive model needed");
  Ntimes <- sum(status); 

  if ((sum(beta)==0) & (beta.fixed==0)) beta<-coxph(Surv(time,time2,status)~Z)$coef; 

  #cat("Cox-Aalen Survival Model"); cat("\n")
  if (px==0) stop("No nonparametric terms (needs one!)");
  ud<-cox.aalenBase(times,ldata,X,Z,
                    status,id,clusters,Nit=Nit,detail=detail,beta=beta,weights=weights,
                    sim=sim,antsim=n.sim,residuals=residuals,robust=robust,
                    weighted.test=weighted.test,ratesim=ratesim,
                    covariance=covariance,resample.iid=resample.iid,namesX=covnamesX,namesZ=covnamesZ,
                    beta.fixed=beta.fixed);

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
                      ud$robvar.gamma<-matrix(0,pz,pz);}
  namematrix(ud$D2linv,covnamesZ); 

  attr(ud,"Call")<-sys.call(); 
  class(ud)<-"cox.aalen"
  attr(ud,"Formula")<-formula;
  attr(ud,"id")<-id.call;
  attr(ud,"cluster")<-cluster.call;
  attr(ud,"start")<-start.time; 
  attr(ud,"time2")<-time2; 
  ud$call<-call

  return(ud); 
}

"plot.cox.aalen" <-  function (x,..., pointwise.ci=1, hw.ci=0,
sim.ci=0, robust=0, specific.comps=FALSE,level=0.05, start.time = 0,
stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
ylab ="Cumulative coefficients",score=FALSE)
{
  object <- x; rm(x);  
  if (!inherits(object,'cox.aalen') ) stop ("Must be output from Cox-Aalen function")

  if (score==FALSE) plot.cums(object, pointwise.ci=pointwise.ci,
        hw.ci=hw.ci,
        sim.ci=sim.ci, robust=robust, specific.comps=specific.comps,level=level,
        start.time = start.time, stop.time = stop.time, add.to.plot=add.to.plot,
        mains=mains, xlab=xlab, ylab =ylab)
  else plotScore(object, specific.comps=specific.comps, mains=mains,
                  xlab=xlab,ylab =ylab);
}


"print.cox.aalen" <-
function (x,...) 
{
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
}


"summary.cox.aalen" <-
function (object,digits = 3,...) 
{
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

    coef.cox.aalen(obj); 

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
}

coef.cox.aalen<-function(object,...,digits=3,d2logl=1) {
   coefBase(object,digits=digits, d2logl=d2logl)
}
