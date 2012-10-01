comp.risk<-function(formula,data=sys.parent(),cause,times=NULL,Nit=50,
clusters=NULL,est=NULL,fix.gamma=0,gamma=0,n.sim=500,weighted=0,model="fg",
causeS=1,cens.code=0,detail=0,interval=0.01,resample.iid=1,
cens.model="KM",cens.formula=NULL,time.pow=NULL,time.pow.test=NULL,silent=1,conv=1e-6,
weights=NULL,max.clust=1000,n.times=50,first.time.p=0.05,estimator=1,
trunc.p=NULL,entry.time=NULL,cens.weight=NULL,admin.cens=NULL,conservative=1) 
# {{{
{ 
## {{{
# trans=1 P_1=1-exp( - ( x' b(b)+ z' gam t) ), 
# trans=2 P_1=1-exp(-exp(x a(t)+ z` b )  Fine-Gray model, with baseline exp(x a(t)) 
# trans=3 P_1= exp(x a(t)+ z` b)/( exp(x a(t) + z' b) +1 );  logistic
# trans=4 P_1=exp( ( x' b(b)+ z' gam ) ), 
# trans=5 P_1= (x' b(t)) exp( z' gam ), 
# trans=6 P_1=1-exp(-(x a(t)) exp(z` b )) Fine-Gray model, with baseline x a(t) 
# trans=7 P_1= (x a(t)) exp( z` b)/( (x a(t) ) exp(z' b) +1 ); logistic2
  trans <- switch(model,additive=1,prop=2,logistic=3,rcif=4,rcif2=5,fg=6,logistic2=7)
###  if (model=="additive")  trans<-1; if (model=="prop")      trans<-2; if (model=="logistic")  trans<-3; 
###  if (model=="rcif")      trans<-4; if (model=="rcif2")     trans<-5; if (model=="fg")        trans<-6; 
###  if (model=="logistic2") trans<-7; 
  line <- 0
  m<-match.call(expand.dots=FALSE);
  m$gamma<-m$times<-m$n.times<-m$cause<-m$Nit<-m$weighted<-m$n.sim<-
    m$model<-m$causeS<- m$detail<- m$cens.model<-m$time.pow<-m$silent<- 
    m$cens.code<-m$cens.formula <- m$interval<- m$clusters<-m$resample.iid<-
    m$time.pow.test<-m$conv<- m$weights  <- m$max.clust <- m$first.time.p<- m$trunc.p <- 
    m$entry.time <- m$cens.weight <- m$admin.cens <- m$fix.gamma <- m$est  <- m$conservative <-     m$estimator <- NULL
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

  if (is.diag(t(X) %*% X)==TRUE) stratum <- 1 else stratum <- 0; 

  if(is.null(clusters)){ clusters <- des$clusters}
  if(is.null(clusters)){
    clusters <- 0:(nrow(X) - 1)
    antclust <- nrow(X)
  } else {
    antclust <- length(unique(clusters))
    clusters <- as.integer(factor(clusters,labels=1:antclust))-1
  }

  coarse.clust <- FALSE; 
  if ((!is.null(max.clust))) if (max.clust< antclust) {
        coarse.clust <- TRUE
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
	antclust <- max.clust    
  }                                                         
  cluster.call<-clusters; 

  pxz <-px+pz;

  if (is.null(times)) {
          timesc<-sort(unique(time2[cause==causeS])); 
	  if (!is.null(n.times)) {
	  if (length(timesc)> n.times) times <- quantile(timesc,prob=seq(first.time.p,1,length=n.times)) 
	  else times <- timesc
	  } else {times<-timesc; times<-times[times> quantile(timesc,prob=first.time.p)]; }
  } else times <- sort(times); 

  n<-nrow(X); ntimes<-length(times);
  if (npar==TRUE) {Z<-matrix(0,n,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-pz;} 
  delta<-(cause!=cens.code)

if (is.null(entry.time)) entry.call <- NULL else entry.call <- 0
if (is.null(entry.time)) entry <- rep(0,n) else entry <- entry.time
if (is.null(weights)==TRUE) weights <- rep(1,n); 
## }}}

## {{{ censoring and estimator 
if (!is.null(admin.cens)) estimator  <- 3;
Gcxe <- 1;  ordertime <- order(time2); 
###dcumhazcens <- rep(0,n); 

if (estimator==1 || estimator==4) {
if (is.null(cens.weight)) { ## {{{ censoring model stuff with possible truncation
  if (cens.model=="KM") { ## {{{
    if (is.null(entry.time)) { ud.cens<-survfit(Surv(time2,cause==cens.code)~+1);
###                               ud.censcum <- aalen(Surv(time2,cause==cens.code)~+1,n.sim=0,robust=0)
    }
    else {
	    ud.cens<-survfit(Surv(entry,time2,cause==cens.code)~+1); 
###            ud.censcum <- aalen(Surv(entry,time2,cause==cens.code)~+1,n.sim=0,robust=0)
    }
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time2)[,2];
###    cumhazcens<-Cpred(ud.censcum$cum,time2[ordertime])[,2];
###    dcumhazcens <- diff(c(0,cumhazcens));
    if (is.null(entry.time)==FALSE) Gcxe<-Cpred(Gfit,entry)[,2];
    Gcx <- Gcx/Gcxe; 
    Gctimes<-Cpred(Gfit,times)[,2]; ## }}}
  } else if (cens.model=="strat-KM") { ## {{{
    if (is.null(entry.time)) ud.cens<-survfit(Surv(time2,cause==cens.code)~XZ) 
    else ud.cens<-survfit(Surv(entry,time2,cause==cens.code)~+XZ); 
    stop("survfit based predictions strat-KM, under construction\n");
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time2)[,2];
    if (is.null(entry.time)==FALSE) Gcxe<-Cpred(Gfit,entry)[,2];
    Gcx <- Gcx/Gcxe; 
    Gctimes<-Cpred(Gfit,times)[,2]; ## }}}
  } else if (cens.model=="cox") { ## {{{
    if (!is.null(cens.formula)) { XZ <- model.matrix(cens.formula,data=data); 
    if (sum(XZ[,1])==nrow(XZ)) XZ <- as.matrix(XZ[,-1])
    } else {
    if (npar==TRUE) XZ<-X[,-1] else XZ <-cbind(X,Z)[,-1];
    }
    if (is.null(entry.time)) ud.cens<-coxph(Surv(time2,cause==cens.code)~XZ)
    else ud.cens<-coxph(Surv(entry,time2,cause==cens.code)~XZ);
    baseout <- basehaz(ud.cens,centered=FALSE); 
    baseout <- cbind(baseout$time,baseout$hazard)
    Gcx<-Cpred(baseout,time2)[,2];
    RR<-exp(XZ %*% coef(ud.cens))
    Gcx<-exp(-Gcx*RR)
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    if (is.null(entry.call)==FALSE) { Gcxe<-Cpred(Gfit,entry)[,2]; Gcxe<-exp(-Gcxe*RR); }
    Gcx <- Gcx/Gcxe; 
    Gctimes<-Cpred(Gfit,times)[,2]; 
    ## }}}
  } else if (cens.model=="aalen") {  ## {{{
    if (!is.null(cens.formula)) { XZ <- model.matrix(cens.formula,data=data); 
    } else {
    if (npar==TRUE) XZ <-X else XZ <-cbind(X,Z);
    }
    if (is.null(entry.time)) ud.cens<-aalen(Surv(time2,cause==cens.code)~-1+XZ+cluster(clusters),
					    n.sim=0,residuals=0,robust=0,silent=1)
    else ud.cens<-aalen(Surv(entry,time2,cause==cens.code)~-1+XZ+cluster(clusters),n.sim=0,silent=1,
			robust=0,residuals=0);
    Gcx <- Cpred(ud.cens$cum,time2)[,-1];
    Gcx<-exp(-apply(Gcx*XZ,1,sum))
    Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    if (is.null(entry.call)==FALSE) {
	    Gcxe<-Cpred(Gfit,entry)[,-1];
            Gcxe<-exp(-apply(Gcxe*XZ,1,sum))
            Gcxe[Gcxe>1]<-1; Gcxe[Gcxe<0]<-0
    }
    Gcx <- Gcx/Gcxe; 
    Gctimes<-Cpred(Gfit,times)[,2]; ## }}}
    } else  stop('Unknown censoring model') 
    cens.weight <- Gcx
    if ((min(Gcx[cause==causeS])< 0.00001) && (silent==0)) { 
	    cat("Censoring dist. zero for some points, summary cens:\n");
	    print(summary(Gcx)) 
    }
## }}}
} else { 
    if (length(cens.weight)!=n) stop("censoring weights must have length equal to nrow in data\n");  
    Gcx <- cens.weight
    Gctimes <- rep(1,length(times)); 
   }
} else { ## estimator==3 admin.cens 
    if (length(admin.cens)!=n) stop("censoring weights must have length equal to nrow in data\n");  
    Gcx <- admin.cens
    Gctimes <- rep(1,length(times)); 
}

   if (is.null(trunc.p)) trunc.p <- rep(1,n);  
   if (length(trunc.p)!=n) stop("truncation weights must have same length as data\n"); 
## }}}

## {{{ setting up more variables

  if (resample.iid == 1) {
    biid <- double(ntimes* antclust * px);
    gamiid<- double(antclust *pg);
  } else {
    gamiid <- biid <- NULL;
  }

  ps<-px; betaS<-rep(0,ps); 

  ## possible starting value for nonparametric components
  if (is.null(est)) { est<-matrix(0,ntimes,px+1); est[,1] <- times; }  
  if (nrow(est)!=length(times)) est <- Cpred(est,times); 

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
  if (is.null(time.pow)==TRUE & model=="rcif" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="rcif2" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="logistic" ) time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="logistic2" ) time.pow<-rep(0,pg); 
  if (length(time.pow)!=pg) time.pow <- rep(time.pow[1],pg); 

  if (is.null(time.pow.test)==TRUE & model=="prop" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="fg" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="additive")  time.pow.test<-rep(1,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif" )    time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif2" )   time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="logistic" ) time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="logistic2" ) time.pow.test<-rep(0,px); 
  if (length(time.pow.test)!=px) time.pow.test <- rep(time.pow.test[1],px); 

  if (ntimes>1) silent <- c(silent,rep(0,ntimes-1)) else silent <- c(silent,0) 
  ## }}}

###  dyn.load("comprisk.so")
  ssf <- 0; 
  out<-.C("itfit", ## {{{
          as.double(times),as.integer(ntimes),as.double(time2),
          as.integer(cens.code), as.integer(cause),as.double(Gcx),
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
	  as.double(conv),as.double(weights),as.double(entry),
	  as.double(trunc.p),as.integer(estimator),as.integer(fix.gamma),
	  as.integer(stratum),as.integer(ordertime-1),as.integer(conservative), 
	  as.double(ssf), PACKAGE="timereg") ## }}}

 ## {{{ handling output
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
           test.procBeqC=Ut,sim.test.procBeqC=UIt,conv=conv,cens.weight=cens.weight)

  ud$call<-call; 
  ud$model<-model; 
  ud$n<-n; 
  ud$formula<-formula; 
  class(ud)<-"comprisk"; 
  attr(ud, "Call") <- sys.call()
  attr(ud, "Formula") <- formula
  attr(ud, "time.pow") <- time.pow
  attr(ud, "cause") <- cause
  attr(ud, "clusters") <- clusters
  attr(ud, "coarse.clust") <- coarse.clust
  attr(ud, "max.clust") <- max.clust
  attr(ud, "causeS") <- causeS
  attr(ud, "cens.code") <- cens.code
  attr(ud, "times") <- times
  return(ud);  ## }}}
} ## }}}

print.comprisk <- function (x,...) { ## {{{
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
} ## }}}

coef.comprisk <- function(object, digits=3,...) { ## {{{
   coefBase(object,digits=digits)
} ## }}}

summary.comprisk <- function (object,digits = 3,...) {  ## {{{
  if (!inherits(object, 'comprisk')) stop ("Must be a comprisk object")
  
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Competing risks Model \n\n")
  
  modelType<-object$model
  #if (modelType=="additive" || modelType=="rcif") 
 
  if (sum(object$obs.testBeq0)==FALSE) cat("No test for non-parametric terms\n") else
  timetest(object,digits=digits); 

  if (semi) { cat("Parametric terms : \n"); 
              out=coef(object); print(signif(out,digits=digits)); cat("   \n"); 
  }

  if (object$conv$convd>=1) {
       cat("WARNING problem with convergence for time points:\n")
       cat(object$cum[object$conv$convp>0,1])
       cat("\nReadjust analyses by removing points\n\n") }

  cat("  Call: \n")
  dput(attr(object, "Call"))
  cat("\n")
} ## }}}

plot.comprisk <-  function (x, pointwise.ci=1, hw.ci=0,
                            sim.ci=0, specific.comps=FALSE,level=0.05, start.time = 0,
                            stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
                            ylab ="Coefficients",score=FALSE,...){
## {{{
  object <- x; rm(x);

  if (!inherits(object,'comprisk') ){
    stop ("Must be output from comp.risk function")
  }

  if (score==FALSE) {
    B<-object$cum;
    V<-object$var.cum;
    p<-dim(B)[[2]]; 

    if (sum(specific.comps)==FALSE){
      comp<-2:p
    } else {
      comp<-specific.comps+1
    }
    if (stop.time==0) {
      stop.time<-max(B[,1]);
    }

    med<-B[,1]<=stop.time & B[,1]>=start.time
    B<-B[med,];
    V<-V[med,]; 

    c.alpha<- qnorm(1-level/2)
    for (v in comp) { 
      c.alpha<- qnorm(1-level/2)
      est<-B[,v];
      ul<-B[,v]+c.alpha*V[,v]^.5;
      nl<-B[,v]-c.alpha*V[,v]^.5;
      if (add.to.plot==FALSE) {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab,...) 
        if (mains==TRUE) title(main=colnames(B)[v]);
      } else {
        lines(B[,1],est,type="s");
      }
      if (pointwise.ci>=1) {
        lines(B[,1],ul,lty=pointwise.ci,type="s");
        lines(B[,1],nl,lty=pointwise.ci,type="s");
      }
      if (hw.ci>=1) {
        if (level!=0.05){
          cat("Hall-Wellner bands only 95 % \n");
        }
        tau<-length(B[,1])
        nl<-B[,v]-1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="s"); 
        lines(B[,1],nl,lty=hw.ci,type="s");
      }
      if (sim.ci>=1) {
        if (is.null(object$conf.band)==TRUE){
          cat("Uniform simulation based bands only computed for n.sim> 0\n")
        }
        if (level!=0.05){
          c.alpha<-percen(object$sim.testBeq0[,v-1],1-level)
        } else {
          c.alpha<-object$conf.band[v-1];
        }
        nl<-B[,v]-c.alpha*V[,v]^.5;
        ul<-B[,v]+c.alpha*V[,v]^.5;
        lines(B[,1],ul,lty=sim.ci,type="s"); 
        lines(B[,1],nl,lty=sim.ci,type="s");
      }
      abline(h = 0)
    }
  } else {
    # plot score proces
    if (is.null(object$pval.testBeqC)==TRUE) {
      cat("Simulations not done \n"); 
      cat("To construct p-values and score processes under null n.sim>0 \n"); 
    } else {
      if (ylab=="Cumulative regression function"){ 
        ylab<-"Test process";
      }
      dim1<-ncol(object$test.procBeqC)
      if (sum(specific.comps)==FALSE){
        comp<-2:dim1
      } else {
        comp<-specific.comps+1
      }

      for (i in comp){
          ranyl<-range(object$test.procBeqC[,i]);
          for (j in 1:50){
            ranyl<-range(c(ranyl,(object$sim.test.procBeqC[[j]])[,i-1]));
          }
          mr<-max(abs(ranyl));

          plot(object$test.procBeqC[,1],
               object$test.procBeqC[,i],
               ylim=c(-mr,mr),lwd=2,xlab=xlab,ylab=ylab,type="s",...)
          if (mains==TRUE){
            title(main=colnames(object$test.procBeqC)[i]);
          }
          for (j in 1:50){
            lines(object$test.procBeqC[,1],
                  as.matrix(object$sim.test.procBeqC[[j]])[,i-1],col="grey",lwd=1,lty=1,type="s")
          }
          lines(object$test.procBeqC[,1],object$test.procBeqC[,i],lwd=2,type="s")
        }
    }
  }
} ## }}}

