comp.risk<-function(formula,data=sys.parent(),cause,times=NULL,Nit=50,
clusters=NULL,gamma=0,n.sim=500,weighted=0,model="additive",
causeS=1,cens.code=0,detail=0,interval=0.01,resample.iid=1,
cens.model="KM",time.pow=NULL,time.pow.test=NULL,silent=1,conv=1e-6){
# trans=1 P_1=1-exp( - ( x' b(b)+ z' gam t) ), 
# trans=2 P_1=1-exp(-exp(x a(t)+ z` b )
# trans=6 P_1=1-exp(-x a(t) exp(z` b )) is not good numerically
# trans=3 logit(P_1)=(x a(t)+ z` b)
# trans=4 P_1=exp( ( x' b(b)+ z' gam t) ), 
# trans=5 P_1= (x' b(t)) exp( z' gam t), 
  if (model=="additive") trans<-1; 
  if (model=="prop")     trans<-2; 
  if (model=="logistic") trans<-3; 
  if (model=="rcif")     trans<-4; 
  if (model=="rcif2")    trans<-5; 
  if (model=="fg")       trans<-6; 
  line <- 0
  m<-match.call(expand = FALSE);
  m$gamma<-m$times<-m$cause<-m$Nit<-m$weighted<-m$n.sim<-
    m$model<-m$causeS<- m$detail<- m$cens.model<-m$time.pow<-m$silent<- 
    m$cens.code<-m$interval<- m$clusters<-m$resample.iid<-
    m$time.pow.test<-m$conv<-NULL
  special <- c("const","cluster")
  if (missing(data)) {
    Terms <- terms(formula, special)
  }  else {
    Terms <- terms(formula, special, data = data)
  }
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  if (attr(m[, 1], "type") == "right") {
    time2 <- m[, 1][, "time"]; time <- rep(0, length(time2))
    status <- m[, 1][, "status"]
  } else if (attr(m[, 1], "type") == "counting") {
    stop("only right censored data"); 
    time <- m[, 1][, 1];time2 <- m[, 1][, 2];status <- m[, 1][, 3];
  } else {
    stop("only right-censored or counting processes data")
  }

  if (n.sim==0) sim<-0 else sim<-1; antsim<-n.sim;
  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ;

  if(is.null(clusters)){ clusters <- des$clusters}

  if(is.null(clusters)){
    clusters <- 0:(nrow(X) - 1)
    antclust <- nrow(X)
  } else {
    clusters <- as.integer(factor(clusters))-1
    antclust <- length(unique(clusters))
  }
    
  pxz <-px+pz;

  if (is.null(times)) {times<-sort(unique(time2[cause==causeS])); 
                       ###times<-times[-c(1:5)];
  } else times <- sort(times); 

  n<-nrow(X); ntimes<-length(times);
  if (npar==TRUE) {Z<-matrix(0,n,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-pz;} 
  delta<-(cause!=cens.code)

  if (cens.model=="KM") {
    ud.cens<-survfit(Surv(time2,cause==cens.code)~+1); 
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time2)[,2];
    Gctimes<-Cpred(Gfit,times)[,2];
  } else if (cens.model=="cox") { 
    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
    ud.cens<-cox.aalen(Surv(time2,cause==cens.code)~prop(XZ),n.sim=0,robust=0);
    Gcx<-Cpred(ud.cens$cum,time2)[,2];
    RR<-exp(XZ %*% ud.cens$gamma)
    Gcx<-exp(-Gcx*RR)
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    Gctimes<-Cpred(Gfit,times)[,2];
    } else if (cens.model=="aalen") { 
    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
    ud.cens<-aalen(Surv(time2,cause==cens.code)~XZ,n.sim=0,robust=0);
    Gcx<-Cpred(ud.cens$cum,time2)[,-1];
    XZ<-cbind(1,XZ); 
    Gcx<-exp(-apply(Gcx*XZ,1,sum))
    Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    Gctimes<-Cpred(Gfit,times)[,2];
    } else { stop('Unknown censoring model') }

   times<-times[Gctimes>interval]; ntimes<-length(times); 

  if (resample.iid == 1) {
    biid <- double(ntimes* antclust * px);
    gamiid<- double(antclust *pg);
  } else {
    gamiid <- biid <- NULL;
  }

  ps<-px; betaS<-rep(0,ps); 

  if (model=="additive") est<-matrix(1/sum(cause==causeS),ntimes,ps+1) 
  else est<-matrix(0,ntimes,ps+1) 

  hess<-matrix(0,ps,ps); var<-score<-matrix(0,ntimes,ps+1); 
  if (sum(gamma)==0) gamma<-rep(0,pg); gamma2<-rep(0,ps); 
  test<-matrix(0,antsim,3*ps); testOBS<-rep(0,3*ps); unifCI<-c();
  testval<-c(); rani<--round(runif(1)*10000); 
  Ut<-matrix(0,ntimes,ps+1); simUt<-matrix(0,ntimes,50*ps);
  var.gamma<-matrix(0,pg,pg); 
  pred.covs.sem<-0

  if (is.null(time.pow)==TRUE & model=="prop" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="fg" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="additive")  time.pow<-rep(1,pg); 
  if (is.null(time.pow)==TRUE & model=="rcif" )     time.pow<-rep(1,pg); 
  if (is.null(time.pow)==TRUE & model=="rcif2" )     time.pow<-rep(1,pg); 
  if (is.null(time.pow)==TRUE & model=="logistic" ) time.pow<-rep(0,pg); 

  if (is.null(time.pow.test)==TRUE & model=="prop" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="fg" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="additive")  time.pow.test<-rep(1,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif" )    time.pow.test<-rep(1,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif2" )   time.pow.test<-rep(1,px); 
  if (is.null(time.pow.test)==TRUE & model=="logistic" ) time.pow.test<-rep(0,px); 

  silent <- c(silent,rep(0,ntimes-1));

  out<-.C("itfit",
          as.double(times),as.integer(ntimes),as.double(time2),
          as.integer(delta), as.integer(cause),as.double(Gcx),
          as.double(X),as.integer(n),as.integer(px),
          as.integer(Nit), as.double(betaS), as.double(score),
          as.double(hess), as.double(est), as.double(var),
          as.integer(sim),as.integer(antsim),as.integer(rani),
          as.double(test), as.double(testOBS), as.double(Ut),
          as.double(simUt),as.integer(weighted),as.double(gamma),
          as.double(var.gamma),as.integer(fixed),as.double(Z),
          as.integer(pg),as.integer(trans),as.double(gamma2),
          as.integer(causeS),as.integer(line),as.integer(detail),
          as.double(biid),as.double(gamiid),as.integer(resample.iid),
          as.double(time.pow),as.integer(clusters),as.integer(antclust),
          as.double(time.pow.test),as.integer(silent),
	  as.double(conv),PACKAGE="timereg")

  gamma<-matrix(out[[24]],pg,1); var.gamma<-matrix(out[[25]],pg,pg); 
  gamma2<-matrix(out[[30]],ps,1); 
  rownames(gamma2)<-covnamesX; 

  conv <- list(convp=out[[41]],convd=out[[42]]); 

  if (fixed==0) gamma<-NULL; 

  if (resample.iid==1)  {
    biid<-matrix(out[[34]],ntimes,antclust*px);
    if (fixed==1) gamiid<-matrix(out[[35]],antclust,pg) else gamiid<-NULL; 
    B.iid<-list();
    for (i in (0:(antclust-1))*px) {
    B.iid[[i/px+1]]<-matrix(biid[,i+(1:px)],ncol=px);
      colnames(B.iid[[i/px+1]])<-covnamesX; }
    if (fixed==1) colnames(gamiid)<-covnamesZ
  } else B.iid<-gamiid<-NULL;

  if (sim==1) {
    simUt<-matrix(out[[22]],ntimes,50*ps); UIt<-list();
    for (i in (0:49)*ps) UIt[[i/ps+1]]<-as.matrix(simUt[,i+(1:ps)]);
    Ut<-matrix(out[[21]],ntimes,ps+1);
    test<-matrix(out[[19]],antsim,3*ps); testOBS<-out[[20]];
    supUtOBS<-apply(abs(as.matrix(Ut[,-1])),2,max);
    p<-ps
    for (i in 1:(3*p)) testval<-c(testval,pval(test[,i],testOBS[i]))
    for (i in 1:p) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
    pval.testBeq0<-as.vector(testval[1:p]);
    pval.testBeqC<-as.vector(testval[(p+1):(2*p)]);
    pval.testBeqC.is<-as.vector(testval[(2*p+1):(3*p)]);
    obs.testBeq0<-as.vector(testOBS[1:p]);
    obs.testBeqC<-as.vector(testOBS[(p+1):(2*p)]);
    obs.testBeqC.is<-as.vector(testOBS[(2*p+1):(3*p)]);
    sim.testBeq0<-as.matrix(test[,1:p]);
    sim.testBeqC<-as.matrix(test[,(p+1):(2*p)]);
    sim.testBeqC.is<-as.matrix(test[,(2*p+1):(3*p)]);
  } else {test<-unifCI<-Ut<-UIt<-pval.testBeq0<-pval.testBeqC<-obs.testBeq0<-
          obs.testBeqC<- sim.testBeq0<-sim.testBeqC<-
          sim.testBeqC.is<- pval.testBeqC.is<-
          obs.testBeqC.is<-NULL;
        }

  est<-matrix(out[[14]],ntimes,ps+1); 
  score<-matrix(out[[12]],ntimes,ps+1); 
  var<-matrix(out[[15]],ntimes,ps+1); 
  colnames(var)<-colnames(est)<-c("time",covnamesX); 

  if (sim>=1) {
    colnames(Ut)<- c("time",covnamesX)
    names(unifCI)<-names(pval.testBeq0)<- names(pval.testBeqC)<- 
    names(pval.testBeqC.is)<- names(obs.testBeq0)<- names(obs.testBeqC)<- 
    names(obs.testBeqC.is)<- colnames(sim.testBeq0)<- colnames(sim.testBeqC)<- 
    colnames(sim.testBeqC.is)<- covnamesX;
  }

  if (fixed==1) { rownames(gamma)<-c(covnamesZ);
                  colnames(var.gamma)<- rownames(var.gamma)<-c(covnamesZ); }

  colnames(score)<-c("time",covnamesX);
  if (is.na(sum(score))==TRUE) score<-NA  else 
  if (sum(score[,-1])<0.00001) score<-sum(score[,-1]); 

  ud<-list(cum=est,var.cum=var,gamma=gamma,score=score,
           gamma2=gamma2,var.gamma=var.gamma,robvar.gamma=var.gamma,
           pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
           obs.testBeq0=obs.testBeq0,
           obs.testBeqC.is=obs.testBeqC.is,
           obs.testBeqC=obs.testBeqC,pval.testBeqC.is=pval.testBeqC.is,
           conf.band=unifCI,B.iid=B.iid,gamma.iid=gamiid,
           test.procBeqC=Ut,sim.test.procBeqC=UIt,conv=conv)

  ud$call<-call; ud$model<-model; ud$n<-n; 
  ud$formula<-formula; class(ud)<-"comprisk"; 
  attr(ud, "Call") <- sys.call()
  attr(ud, "Formula") <- formula
  attr(ud, "time.pow") <- time.pow
  attr(ud, "cause") <- cause
  attr(ud, "times") <- times
  return(ud); 
}

print.comprisk <- function (x,...) {
  object <- x; rm(x);
  if (!inherits(object, 'comprisk')) stop ("Must be an comprisk object")
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
   # We print information about object:  
  cat(paste("Competing risks model with",object$model,"\n\n"))
  cat(" Nonparametric terms : ");
  cat(colnames(object$cum)[-1]); cat("   \n");  
  if (semi) {
    cat(" Parametric terms :  ");
    cat(rownames(object$gamma)); 
    cat("   \n");
  } 

  if (object$conv$convd>=1) {
       cat("Warning problem with convergence for time points:\n")
       cat(object$cum[object$conv$convp>0,1])
       cat("\nReadjust analyses by removing points\n") }
  cat("   \n");  
}

coef.comprisk <- function(object, digits=3,...) {

   coefBase(object,digits=digits)
}

summary.comprisk <- function (object,digits = 3,...) {
  if (!inherits(object, 'comprisk')) stop ("Must be a comprisk object")
  
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Competing risks Model \n\n")
  
  modelType<-object$model
  #if (modelType=="additive" || modelType=="rcif") 
 
  if (sum(object$obs.testBeq0)==FALSE) cat("No test for non-parametric terms\n") else
  timetest(object,digits=digits); 

  if (semi) { cat("Parametric terms : \n"); coef(object); cat("   \n"); }

  if (object$conv$convd>=1) {
       cat("WARNING problem with convergence for time points:\n")
       cat(object$cum[object$conv$convp>0,1])
       cat("\nReadjust analyses by removing points\n\n") }

  cat("  Call: \n")
  dput(attr(object, "Call"))
  cat("\n")
}



