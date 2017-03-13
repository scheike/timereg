haplo.cif<-function(formula,data=sys.parent(),cause,times,
designfuncX,designfuncZ,Nit=50,match=FALSE,
clusters=NULL,gamma=0,n.sim=500,weighted=0,model="additive",
causeS=1,cens.code=0,detail=0,interval=0.01,resample.iid=1,
cens.model="KM",time.pow=0,fix.haplofreq=0,haplo.freq=NULL,alpha.iid=NULL,
geno.setup=NULL,fit.haplofreq=NULL,design.test=0,
covnamesX=NULL,covnamesZ=NULL){
## {{{ setting up models 
# trans=1 P_1=1-exp(- ( x' b(b)+ z' gam t) ), 
# trans=2 P_1=1-exp(-exp(x a(t)+ z` b )
# trans=3 P_1= exp(x a(t)+ z` b)/( exp(x a(t) + z' b) +1 );  logistic
  trans <- switch(model,additive=1,prop=2,logistic=3) ### ,rcif=4,rcif2=5,fg=6,logistic2=7)
  if (trans==1) line<-1; if (trans==2) line<-0; if (trans==3) line<-0; 
# line=1 indicates that it is tested that "b(t) = gamma t".
# line=0 indicates that it is tested that "b(t) = gamma ".
  m<-match.call(expand = FALSE);
  m$match<-m$cause<-m$times<-m$designfuncX<-m$designfuncZ<-m$Nit<-
  m$clusters<-m$gamma<-m$n.sim<-m$weighted<-m$model<-
  m$causeS<-m$cens.code<-m$detail<-m$interval<-m$resample.iid<-
  m$cens.model<-m$time.pow<-m$fix.haplofreq<-m$haplo.freq<-m$alpha.iid<-
  m$geno.setup<-m$fit.haplofreq<-m$design.test<-NULL
  m$covnamesX<- m$covnamesZ<-NULL
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
  covnamesXin<-des$covnamesX; covnamesZin<-des$covnamesZ;
  antpers<-nrow(X); 

  if(is.null(clusters)){clusters<- des$clusters}

  if(is.null(clusters)){
    clusters <- 0:(nrow(X) - 1)
    antclust <- nrow(X)
  } else {
    clusters <- as.integer(factor(clusters))-1
    antclust <- length(unique(clusters))
  }
  if (npar==TRUE) {Z<-matrix(0,antpers,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-pz;} 
## }}}

if (match==FALSE) {
## {{{ haplo-designs for cif
if (fixed==0) {
sdesXcheck<- function(x,h) {
    out <- designfuncX(x,h)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
  }
xih<-sdesXcheck(X[1,],c(0,1))
} else  {
sdesXcheck<- function(x,z,h) {
    out <- designfuncX(x,z,h)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
  }
xih<-sdesXcheck(X[1,],Z[1,],c(0,1))
}
sdesZcheck<- function(x,z,h) {
    out <- designfuncZ(x,z,h)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
if (npar==FALSE) { zih<-sdesZcheck(X[1,],Z[1,],c(0,1));
dimzih<-length(zih); } else dimzih<-1; 
dimxih<-length(xih)
## }}}
} else {
## {{{ haplo-designs for matched cif
if (fixed==0) {
smdesXcheck<- function(x,hd,hp) {
    out <- designfuncX(x,hd,hp)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
  }
xih<-smdesXcheck(X[1,],c(0,1),c(0,1))
} else  {
smdesXcheck<- function(x,z,hd,hp) {
    out <- designfuncX(x,z,hd,hp)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
  }
xih<-smdesXcheck(X[1,],Z[1,],c(0,1),c(0,1))
}
smdesZcheck<- function(x,z,hd,hp) {
    out <- designfuncZ(x,z,hd,hp)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
if (npar==FALSE) { zih<-smdesZcheck(X[1,],Z[1,],c(0,1),c(0,1));
dimzih<-length(zih); } else dimzih<-1; 
dimxih<-length(xih)
#if (dimzih!=pz) {npar<-FALSE; fixed<-1; } 
#print(dim(Z)); print(fixed); print(dimzih); print(pg);
#print(pz); print(dimxih); print(px);
## }}}
}

 ## {{{ genostuff
if (is.null(geno.setup)) {
    stop("Must provide geno setup, see example in help file ");
} else setup<-geno.setup
    HPIordered <- setup$HPIordered
    oh <- unlist(HPIordered) -1 
    nph<-nHaps <-length(setup$uniqueHaploNames)
    nAllelesPerLocus <- setup$nAllelesPerLocus
    unorderedAlleles <- setup$unorderedAlleles
    nPeople <- setup$nPeople
    if (nPeople!=antpers)stop("Number of subjects not the same for genotype and survival data"); 
    nLoci <- setup$nLoci
    nPossHapPairsPerPerson <- sapply(HPIordered, length)
    if (is.null(fit.haplofreq)) {
    if (is.null(haplo.freq)) 
      stop("Must provide either haplo frequencies from fit or fixed haplo frequencies"); 
    } else
    { 
      if (is.null(haplo.freq)) {  
      haplo.freq<-fit.haplofreq$haplo.freq; 
      alpha.iid<-fit.haplofreq$alpha.iid
      } else haplo.freq<-haplo.freq; 
    }
     haplo.mass<-sum(haplo.freq);
     hapdim<-ncol(alpha.iid); 

   if (fix.haplofreq==1) { hapdim<-1; haplo.design<-diag(hapdim);} else 
   haplo.design<-fit.haplofreq$haplo.design
## }}}
  
## {{{ setting censoring weights + variables 

  pxz <-px+pz;
  ps<-dimxih; betaS<-rep(0,ps); 

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
    } else { stop('Unknown censoring model') }

  times<-times[Gctimes>interval]; ntimes<-length(times); 

  if (resample.iid == 1) {
    biid <- matrix(0, ntimes, antclust * dimxih);
    gamiid<- matrix(0, antclust ,dimzih);
  } else { gamiid <- biid <- NULL; }


  if (model=="additive") est<-matrix(1/sum(cause==causeS),ntimes,ps+1) 
  else est<-matrix(0,ntimes,ps+1) 

  hess<-matrix(0,ps,ps); var<-score<-matrix(0,ntimes,ps+1); 
  if (sum(gamma)==0) gamma<-rep(0,dimzih); gamma2<-rep(0,ps); 
  test<-matrix(0,antsim,3*ps); testOBS<-rep(0,3*ps); unifCI<-c();
  testval<-c(); rani<--round(runif(1)*10000); 
  Ut<-matrix(0,ntimes,ps+1); simUt<-matrix(0,ntimes,50*ps);
  var.gamma<-matrix(0,dimzih,dimzih); 

  if (sum(time.pow)==0 & model=="prop") time.pow<-rep(0,dimzih); 
  if (sum(time.pow)==0 & model=="additive") time.pow<-rep(1,dimzih); 
## }}}

#dyn.load("haplo.so")
if (match==FALSE) { 
## {{{ calling c haplocif 
  out<-.C("haplocif",
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
          as.integer(fix.haplofreq),as.double(haplo.freq),as.double(alpha.iid),
          as.integer(hapdim),
          as.integer(nph), as.integer(oh), as.integer(nPossHapPairsPerPerson),
          body(sdesXcheck), body(sdesZcheck),new.env(),as.integer(dimxih),
          as.integer(dimzih),as.double(haplo.design),as.integer(design.test) 
          ,PACKAGE="HaploSurvival")
## }}}
} else 
{
  out<-.C("cifhaplomatch", ## {{{
         as.double(times),as.integer(ntimes),as.double(time2), 
         as.integer(delta), as.integer(cause),as.double(Gcx),
         as.double(X),as.integer(n),as.integer(px),
         as.integer(Nit), as.double(betaS), as.double(score),
         as.double(hess), as.double(est), as.double(var),
         as.integer(sim),as.integer(antsim),as.integer(rani),  ### 5 
         as.double(test), as.double(testOBS), as.double(Ut),
         as.double(simUt),as.integer(weighted),as.double(gamma),
         as.double(var.gamma),as.integer(fixed),as.double(Z),
         as.integer(pg),as.integer(trans),as.double(gamma2),
         as.integer(causeS),as.integer(line),as.integer(detail), ### 10 
         as.double(biid),as.double(gamiid),as.integer(resample.iid),
         as.double(time.pow),as.integer(clusters),as.integer(antclust),
         as.integer(fix.haplofreq),as.double(haplo.freq),as.double(alpha.iid),
         as.integer(hapdim),as.integer(nph),as.integer(oh), 
         as.integer(nPossHapPairsPerPerson),body(smdesXcheck),body(smdesZcheck),##15
         new.env(),as.integer(dimxih), as.integer(dimzih),
         as.double(haplo.design),as.integer(design.test),PACKAGE="HaploSurvival")
## }}}
}

## {{{ output handling 
    gamma<-matrix(out[[24]],dimzih,1); 
    var.gamma<-matrix(out[[25]],dimzih,dimzih); 
    gamma2<-matrix(out[[30]],ps,1); 
    if (fixed==0) gamma<-NULL; 
    if (resample.iid==1)  {
        biid<-matrix(out[[34]],ntimes,antclust*dimxih);
        if (fixed==1) gamiid<-matrix(out[[35]],antclust,dimzih) else gamiid<-NULL; 
        B.iid<-list();
        for (i in (0:(antclust-1))*dimxih) {
            B.iid[[i/dimxih+1]]<-as.matrix(biid[,i+(1:dimxih)]);
        }
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
    } else {test<-FALSE; unifCI<-FALSE; Ut<-FALSE; UIt<-FALSE;
        pval.testBeq0<-FALSE;pval.testBeqC<-FALSE; obs.testBeq0<-FALSE;obs.testBeqC<-FALSE;
        sim.testBeq0<-FALSE;sim.testBeqC<-FALSE;
        sim.testBeqC.is<-FALSE; pval.testBeqC.is<-FALSE;
        obs.testBeqC.is<-FALSE;
    }
    est<-matrix(out[[14]],ntimes,ps+1); 
    score<-matrix(out[[12]],ntimes,ps+1); 
    var<-matrix(out[[15]],ntimes,ps+1); 
## }}}

## {{{ naming output
  ## best guess  for standard designs  
  haploX.names<-c()
  if (px!=dimxih) {
    haploef.x<-dimxih-px; 
    if (haploef.x>0) haploX.names<-rep("Haplo effect",haploef.x)
  }
  haploZ.names<-c()
  if (npar==FALSE) if (pz!=dimzih) {
    haploef.z<-dimzih-pz; 
    if (haploef.z>0) haploZ.names<-rep("Haplo effect",haploef.z)
    
  }

  if (is.null(covnamesX)==TRUE)  
     covnamesXuse<-c(covnamesXin,haploX.names)  
  else covnamesXuse<-c(covnamesX)

  if (is.null(covnamesZ)==TRUE)  
     covnamesZuse<-c(covnamesZin,haploZ.names)  
  else covnamesZuse<-c(covnamesZ)

  if (length(covnamesXuse)==ps) {
     colnames(var)<-colnames(est)<-c("time",covnamesXuse); 
     rownames(gamma2)<-c(covnamesXuse); 
     colnames(score)<-c("time",covnamesXuse);
  }

  if ((sim>=1) & (length(covnamesXuse)==ps)) {
    colnames(Ut)<- c("time",covnamesXuse)
    names(unifCI)<-names(pval.testBeq0)<-
    names(pval.testBeqC)<- names(pval.testBeqC.is)<-
    names(obs.testBeq0)<- names(obs.testBeqC)<- names(obs.testBeqC.is)<-
    colnames(sim.testBeq0)<-
    colnames(sim.testBeqC)<- colnames(sim.testBeqC.is)<-covnamesXuse;
  }

  if ((length(covnamesZuse)==dimzih) & (fixed==1))
  { rownames(gamma)<-colnames(var.gamma)<- rownames(var.gamma)<-covnamesZuse; 
  }

  if (is.na(sum(score))==TRUE) score<-NA  else 
  if (sum(score[,-1])<0.00001) score<-sum(score[,-1]); 

  ud<-list(cum=est,var.cum=var,gamma=gamma,score=score,
           gamma2=gamma2,var.gamma=var.gamma,robvar.gamma=var.gamma,
           pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
           obs.testBeq0=obs.testBeq0,
           obs.testBeqC.is=obs.testBeqC.is,
           obs.testBeqC=obs.testBeqC,pval.testBeqC.is=pval.testBeqC.is,
           conf.band=unifCI,B.iid=B.iid,gamma.iid=gamiid,
           test.procBeqC=Ut,sim.test.procBeqC=UIt,KMweights=Gcx)
## }}}

    ud$conv$convd <- 0 ### to make compatible with comp.risk output
    if (model=="prop")      time.pow<-rep(0,pg); 
    if (model=="additive")  time.pow<-rep(1,pg); 

    ud$call<-call; ud$model<-model; ud$n<-antpers; 
    ud$formula<-formula; class(ud)<-"comprisk"; 
    attr(ud, "Call") <- sys.call(); 
    attr(ud, "Formula") <- formula
    attr(ud, "time.pow") <- time.pow
    return(ud); 
}
