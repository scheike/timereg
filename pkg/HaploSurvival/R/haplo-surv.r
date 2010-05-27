haplo.surv<- function (formula = formula(data),data=sys.parent(),
    designfuncX,designfuncZ, beta = 0, match=FALSE, 
    Nit = 10, detail = 0, start.time = 0, max.time = NULL, id = NULL, n.sim = 500, 
    geno.type = NULL, geno.setup=NULL, haplo.freq = NULL,
    fix.beta = 0, fix.haplofreq = 0, two.stage = 0, weighted.test = 0,
    step = 1, lev.marq=1, min.lev.marq=0,
    haplo.design=NULL,haplo.baseline=NULL,alpha=NULL,resample.iid=1,
    covnamesX=NULL,covnamesZ=NULL)
{
  ## {{{ initialization  and setup
  sym <- 0 ; id.call <- id; residuals <- 0; robust <- 0
  call <- match.call()
  m <- match.call(expand = FALSE)
  m$geno.setup<-m$geno.type <- m$max.time <- m$start.time <- 
  m$weighted.test <- m$n.sim <- m$id <- m$Nit <- m$detail <-m$beta<- 
  m$haplo.freq <- m$fix.beta <- m$fix.haplofreq <-m$two.stage<-m$step <- 
  m$lev.marq<- m$min.lev.marq<- m$designfuncX<-m$designfuncZ<-
  m$haplo.baseline<-m$haplo.design<- m$alpha<-m$resample.iid<-m$match<-NULL
  m$covnamesX<- m$covnamesZ<-NULL
  if (n.sim == 0) sim <- 0 else sim <- 1
  antsim <- n.sim
  special <- c("prop")
  Terms <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  des <- read.design(m, Terms, model = "cox.aalen")
  X <- des$X; Z <- des$Z; 
  npar <- des$npar; px <- des$px; pz <- des$pz; pg <- pz; nx <- nrow(X); 
  covnamesXin <- des$covnamesX; covnamesZin <- des$covnamesZ; 
  clusters <- des$clusters; pxz <- px + pz; 
  survs <- read.surv(m, id, npar, clusters, start.time, max.time,model="cox.aalen")
  times <- survs$times; id <- id.call <- survs$id.cal
  clusters <- cluster.call <- survs$clusters
  time <- survs$start; time2 <- survs$stop; status <- survs$status; 
  antpers <- survs$antpers; antclust <- survs$antclust; 
  clustsize<-as.vector(table(clusters));
  maxclust<-max(clustsize); clusterindex<-matrix(0,antclust,maxclust);
  if (antclust!=antpers) {
    for (i in 1:antclust) { index<-(((1:antpers)[clusters==i])-1);
                            clusterindex[i,1:length(index)]<-(((1:antpers)[clusters==i])-1) }
  } else clusterindex<-(0:(antpers-1));

  if (npar == TRUE) {Z<-matrix(0,nrow(X),1); beta<-0; pg<-pz<-1; }
  if (npar ==TRUE & fix.haplofreq==1) Nit<-1; 
  Ntimes <- length(times)   
  ## }}}

  ## {{{ genostuff
  if (is.null(geno.setup)) {
    setup <- geno.setup(geno.type,
                        haplo.baseline=haplo.baseline)
  } else setup<-geno.setup
  HPIordered <- setup$HPIordered
  uniqueHaplos <- setup$uniqueHaploNames
  nph<-nHaps <-length(uniqueHaplos)
  nAllelesPerLocus <- setup$nAllelesPerLocus
  unorderedAlleles <- setup$unorderedAlleles
  nPeople <- setup$nPeople
  if (nPeople!=antpers)stop("Number of subjects not the same for genotype and survival data"); 
  nLoci <- setup$nLoci
  nPossHaps <- sapply(HPIordered, length)
  if (is.null(haplo.freq)) {
    haplo.freq <- as.vector(table(unlist(setup$HPIordered))/sum(table(unlist(setup$HPIordered))))
  }
  else {
    haplo.mass<-sum(haplo.freq);
    haplo.freq[length(uniqueHaplos)]<-1-sum(haplo.freq[1:(length(uniqueHaplos)-1)])
  }

  rho<-NULL; if (is.null(rho)) rho <- rep(0, nPeople); 

  if (sum(haplo.freq<=0)>0) cat("Warning, some haplofrequencies are 0\n");

  if (fix.haplofreq==0) {
    i0<-(haplo.freq<=0); n0<-sum(i0)
    il0<-(haplo.freq>0); nl0<-sum(il0)
    haplo.freq[i0]<-0.0001
    haplo.freq[il0]<-haplo.freq[il0]-0.0001* sum(i0)/sum(il0); 
    haplo.pars<-log(haplo.freq[1:(nHaps - 1)]/haplo.freq[nHaps])
  }
  else haplo.pars<-rep(0,nph-1); 
    

  if (is.null(haplo.design)==TRUE)  {
    haplo.design<-diag(nph-1) 
    dimhap<-nph-1; alpha<-haplo.pars;   # full haplo-frequency model 
  } else { 
    dimhap<- ncol(haplo.design);        # haplo-freq parameters
    if (fix.haplofreq==0) {
      init.alpha<-(t(haplo.design) %*% haplo.pars) / 
        (t(haplo.design) %*% rep(1,nrow(haplo.design)) )
    } else init.alpha<-rep(0,dimhap);  
    if (is.null(alpha)==TRUE) alpha<-init.alpha; 
  }
  #print(haplo.freq); print(haplo.pars); print(alpha)

  oh <- orderedHaplos <- unlist(HPIordered) -1 
  nph <- length(uniqueHaplos)
  Rho <- rho
  LogLike <- numeric(1)
  survscoregeno <- scoregeno <- Score <- numeric(nph - 1)
  d2lgeno <- D2L <- numeric((nph-1)^2)
  ## }}}

if (match==FALSE) {
#print(paste("match er ",match))
## {{{ haplo-designs
sdesXcheck<- function(x,z,h) {
    out <- designfuncX(x,z,h)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
sdesZcheck<- function(x,z,h) {
    out <- designfuncZ(x,z,h)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}

#if (dimzih!=pz) {npar<-FALSE; fixed<-1; } 
#print(dim(Z)); print(fixed); print(dimzih); print(pg);
#print(pz); print(dimxih); print(px);
xih<-sdesXcheck(X[1,],Z[1,],c(0,1))
zih<-sdesZcheck(X[1,],Z[1,],c(0,1))
dimxih<-length(xih)
dimzih<-length(zih)
## }}}
} else {
## {{{ setup haplo-design 
smdesXcheck<- function(x,z,hd,hp) {
    out <- designfuncX(x,z,hd,hp)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
smdesZcheck<- function(x,z,hd,hp) {
    out <- designfuncZ(x,z,hd,hp)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}

xih<-smdesXcheck(X[1,],Z[1,],c(0,1),c(0,1))
zih<-smdesZcheck(X[1,],Z[1,],c(0,1),c(0,1))
dimxih<-length(xih)
dimzih<-length(zih)
## }}}
}

  ## {{{ setup of furhter variables
  dimpar <- (fix.beta == 0) * dimzih + (fix.haplofreq == 0) * dimhap
  if (dimpar == 0) dimpar <- 1
  if (fix.beta == 1 && fix.haplofreq == 1) { Nit <- 1 }
  if (is.null(Z) == TRUE) { Z <- matrix(0, nx, 1); beta <- 0 }
  if (is.null(beta)) beta <- rep(0, dimzih)
  if (sum(abs(beta))==0) beta <- rep(0, dimzih)
  if (residuals == 1) { 
    cumAi <- matrix(0, Ntimes, antpers * 1); 
    cumAiiid <- matrix(0, Ntimes, antpers * 1)
  }
  else { cumAi <- 0; cumAiiid <- 0 }
  cumint <- matrix(0, Ntimes, dimxih + 1)
  vcum <- matrix(0, Ntimes, dimxih + 1)
  Rvcu <- matrix(0, Ntimes, dimxih + 1)
  RVarbeta <- Varbeta <- matrix(0, dimpar, dimpar)
  score <- rep(0, dimpar)
  Iinv <- matrix(0, dimpar, dimpar)
  varpar <- matrix(0, dimpar, dimpar)
  d2score <- matrix(0, dimpar, dimpar)
  Uit <- FALSE
  if (sim == 1 && fix.beta==0) Uit <- matrix(0, Ntimes, 50 * dimzih) 
  test <- matrix(0, antsim, 2 * dimxih)
  testOBS <- rep(0, 2 * dimxih)
  unifCI <- c()
  testval <- c()
  Ut <- matrix(0, Ntimes, dimpar + 1)
  simUt <- matrix(0, antsim, dimpar)
  loglike <- 0

  if (resample.iid == 1) {
    biid <- matrix(0, Ntimes, antclust * dimxih); 
    gamiid<- matrix(0,antclust,dimzih); }
  else {gamiid <- biid <- NULL; }
  pars <- rep(0, dimpar)
  ## }}}

if (attr(m[,1],"type")=="right") { ## {{{ sorting for right censored case
 ot<-order(time2,status==0); # order in time, status=1 first for ties
 time2<-time2[ot]; status<-status[ot]
 X<-X[ot,]; Z<-Z[ot,]
 clusters<-survs$clusters<-survs$clusters[ot]
 index.dtimes<-which(time2 %in% times)-1; 

 ## also sort genotype data similarly
 newHPIordered<-list(); 
 for (i in 1:antpers) newHPIordered[[i]]<-setup$HPIordered[[ ot[i] ]]
 nPossHaps<-setup$nPossHaps[ot]
 oh <- orderedHaplos <- unlist(newHPIordered) -1 
}
## }}}

###dyn.load("haplo.so")

if (match==FALSE) {

if (attr(m[,1],"type")=="right") {
  ## {{{ calling c routine for cox-aalen model 
nparout <- .C("simplehaplosurvdes", 
 as.double(times), as.integer(Ntimes), as.double(X), as.integer(nx), as.integer(px), 
 as.double(Z), as.integer(nx), as.integer(pg), as.integer(antpers), as.double(time), 
 as.double(time2), as.double(beta),  as.integer(Nit), as.double(cumint), as.double(vcum),
 as.double(loglike), as.double(Iinv), as.double(Varbeta), as.integer(detail), as.integer(sim), 
 as.integer(antsim), as.double(Rvcu), as.double(RVarbeta), as.double(test), as.double(testOBS),
 as.double(Ut), as.double(simUt), as.double(Uit), as.integer(weighted.test), 
as.integer(index.dtimes), 
as.integer(status), as.double(score), as.double(cumAi), as.double(cumAiiid), as.integer(residuals),
 as.integer(sym), as.integer(nph), as.integer(oh), as.integer(nPossHaps), as.double(haplo.pars), 
as.integer(fix.beta), as.integer(fix.haplofreq), as.integer(dimzih), as.integer(dimxih), as.double(rho),
as.double(scoregeno), as.double(d2lgeno), as.double(survscoregeno), as.integer(two.stage), as.double(varpar), 
 as.double(d2score), as.double(step), as.double(pars), as.double(lev.marq), as.double(min.lev.marq),
 body(sdesXcheck), body(sdesZcheck),new.env(), as.double(haplo.design),as.double(alpha),
 as.integer(dimhap), as.double(haplo.freq),as.double(biid),as.double(gamiid),
 as.integer(resample.iid),PACKAGE="HaploSurvival") 
  ## }}}
} else {
  ## {{{ calling c routine for haplo-surv cox-aalen model 
  nparout <- .C("haplosurvdes", 
                as.double(times), as.integer(Ntimes), as.double(X), as.integer(nx), as.integer(px), 
                as.double(Z), as.integer(nx), as.integer(pg), as.integer(antpers), as.double(time), 
                as.double(time2), as.double(beta),  as.integer(Nit), as.double(cumint), as.double(vcum),
                as.double(loglike), as.double(Iinv), as.double(Varbeta), as.integer(detail), as.integer(sim), 
                as.integer(antsim), as.double(Rvcu), as.double(RVarbeta), as.double(test), as.double(testOBS),
                as.double(Ut), as.double(simUt), as.double(Uit), as.integer(weighted.test), as.integer(id), 
                as.integer(status), as.double(score), as.double(cumAi), as.double(cumAiiid), as.integer(residuals),
                as.integer(sym), as.integer(nph), as.integer(oh), as.integer(nPossHaps), as.double(haplo.pars), 
                as.integer(fix.beta), as.integer(fix.haplofreq), as.integer(dimzih), as.integer(dimxih), as.double(rho),
                as.double(scoregeno), as.double(d2lgeno), as.double(survscoregeno), as.integer(two.stage), as.double(varpar), 
                as.double(d2score), as.double(step), as.double(pars), as.double(lev.marq), as.double(min.lev.marq),
                body(sdesXcheck), body(sdesZcheck),new.env(), as.double(haplo.design),as.double(alpha),
                as.integer(dimhap), as.double(haplo.freq),as.double(biid),as.double(gamiid),
                as.integer(resample.iid),PACKAGE="HaploSurvival") 
  ## }}}
}
} else {
if (attr(m[,1],"type")=="right") {
  ## {{{ calling c routine for haplo match case 
 nparout <- .C("simplehaplosurvmatch", 
as.double(times), as.integer(Ntimes), as.double(X), as.integer(nx), as.integer(px), 
as.double(Z), as.integer(nx), as.integer(pg), as.integer(antpers), as.double(time), 
as.double(time2), as.double(beta),  as.integer(Nit), as.double(cumint), as.double(vcum),
as.double(loglike), as.double(Iinv), as.double(Varbeta), as.integer(detail), as.integer(sim), 
as.integer(antsim), as.double(Rvcu), as.double(RVarbeta), as.double(test), as.double(testOBS),
as.double(Ut), as.double(simUt), as.double(Uit), as.integer(weighted.test), 
as.integer(index.dtimes), 
as.integer(status), as.double(score), as.double(cumAi), as.double(cumAiiid), as.integer(residuals),
as.integer(sym), as.integer(nph), as.integer(oh), as.integer(nPossHaps), as.double(haplo.pars), 
as.integer(fix.beta), as.integer(fix.haplofreq), as.integer(dimzih), as.integer(dimxih), as.double(rho),
as.double(scoregeno), as.double(d2lgeno), as.double(survscoregeno), as.integer(two.stage), as.double(varpar), 
as.double(d2score), as.double(step), as.double(pars), as.double(lev.marq), as.double(min.lev.marq),
body(smdesXcheck), body(smdesZcheck),new.env(), as.double(haplo.design),as.double(alpha),
as.integer(dimhap), as.double(haplo.freq),as.double(biid),as.double(gamiid),
as.integer(resample.iid),PACKAGE="HaploSurvival") 
## }}}
}
else {
## {{{ calling c routine for haplo match case 
nparout <- .C("haplosurvmatch", 
as.double(times), as.integer(Ntimes), as.double(X), as.integer(nx), as.integer(px), 
as.double(Z), as.integer(nx), as.integer(pg), as.integer(antpers), as.double(time), 
as.double(time2), as.double(beta),  as.integer(Nit), as.double(cumint), as.double(vcum),
as.double(loglike), as.double(Iinv), as.double(Varbeta), as.integer(detail), as.integer(sim), 
as.integer(antsim), as.double(Rvcu), as.double(RVarbeta), as.double(test), as.double(testOBS),
as.double(Ut), as.double(simUt), as.double(Uit), as.integer(weighted.test), 
as.integer(id), 
as.integer(status), as.double(score), as.double(cumAi), as.double(cumAiiid), as.integer(residuals),
as.integer(sym), as.integer(nph), as.integer(oh), as.integer(nPossHaps), as.double(haplo.pars), 
as.integer(fix.beta), as.integer(fix.haplofreq), as.integer(dimzih), as.integer(dimxih), as.double(rho),
as.double(scoregeno), as.double(d2lgeno), as.double(survscoregeno), as.integer(two.stage), as.double(varpar), 
as.double(d2score), as.double(step), as.double(pars), as.double(lev.marq), as.double(min.lev.marq),
body(smdesXcheck), body(smdesZcheck),new.env(), as.double(haplo.design),as.double(alpha),
as.integer(dimhap), as.double(haplo.freq),as.double(biid),as.double(gamiid),
as.integer(resample.iid),PACKAGE="HaploSurvival") 
## }}}
}
}

  ## {{{ output 
  gamma <- matrix(nparout[[53]], dimpar, 1)
  cumint <- matrix(nparout[[14]], Ntimes, dimxih + 1)
  vcum <- matrix(nparout[[15]], Ntimes, dimxih + 1)
  Iinv <- matrix(nparout[[17]], dimpar, dimpar)
  Rvcu <- matrix(nparout[[22]], Ntimes, dimxih + 1)
  Varbeta <- matrix(nparout[[18]], dimpar, dimpar)
  RVarbeta <- matrix(nparout[[23]], dimpar, dimpar)
  score <- c(nparout[[32]])
  Ut <- matrix(nparout[[26]], Ntimes, dimpar + 1)
  loglike <- nparout[[16]]
  if (residuals == 1) {
    cumAi <- matrix(nparout[[33]], Ntimes, antpers * 1)
    cumAiiid <- matrix(nparout[[34]], Ntimes, antpers * 1)
    cumAi <- list(time = times, dmg = cumAi, dmg.iid = cumAiiid)
  }
  else cumAi <- FALSE
  if (sim == 1) {
    if (fix.beta==0)  {
      Uit <- matrix(nparout[[28]],Ntimes,50*dimzih)
      UIt <- list()
      for (i in (0:49) * dimzih) UIt[[i/dimzih+ 1]]<-as.matrix(Uit[,i+(1:dimzih)])
      simUt <-as.matrix( matrix(nparout[[27]], antsim, dimpar)[,1:dimzih])
      supUtOBS <- apply(abs(as.matrix(Ut[, 2:(dimzih+1)])), 2, max)
      testUt <- c()
      for (i in 1:dimzih) testUt<-c(testUt, pval(simUt[, i], supUtOBS[i]))
      sim.supUt <- as.matrix(simUt)
    } else {
      UIt<-simUt<-supUtOBS<-testUt<-sim.supUt<- NULL; 
    }

    test <- matrix(nparout[[24]], antsim, 2 * dimxih)
    testOBS <- nparout[[25]]
    for (i in 1:(2 * dimxih)) testval <- c(testval, pval(test[,i], testOBS[i]))
    for (i in 1:dimxih) unifCI <- c(unifCI, percen(test[, i], 0.95))
    pval.testBeq0 <- as.vector(testval[1:dimxih])
    pval.testBeqC <- as.vector(testval[(dimxih + 1):(2 * dimxih)])
    obs.testBeq0 <- as.vector(testOBS[1:dimxih])
    obs.testBeqC <- as.vector(testOBS[(dimxih + 1):(2 * dimxih)])
    sim.testBeq0 <- as.matrix(test[, 1:dimxih])
    sim.testBeqC <- as.matrix(test[, (dimxih + 1):(2 * dimxih)])
  }
  if (sim != 1) {
    testUt <- test <- unifCI <- supUtOBS <- UIt <- testOBS <- testval <- 
      pval.testBeq0 <- pval.testBeqC <- obs.testBeq0 <- obs.testBeqC <- 
        sim.testBeq0 <- sim.testBeqC <- testUt <- sim.supUt <- NULL
  }
  if (fix.haplofreq== 0) {
    smed<-(fix.beta==0)*dimzih+1; 
    smed<-smed:dimpar
    alpha<-gamma[smed]
    var.hap.alpha<-Varbeta[smed,smed]
    robvar.hap.alpha<-RVarbeta[smed,smed]
    Iinv.hap.alpha<- Iinv[smed,smed]
    haplo.pars<-haplo.design %*% alpha; 
    haplo.freqs = exp(c(haplo.pars,0))/sum(exp(c(haplo.pars, 0)))
  }
  else {robvar.hap.alpha<-var.hap.alpha<-var.haplo.pars<- 
          Iinv.hap.alpha<- haplo.pars<- 
            robvar.haplo.pars<- Iinv.haplo.pars<-NULL
        }
  var.gamma<-robvar.gamma<- Iinv.gamma<-matrix(0,1,1); 
  if (fix.beta==0) {
    gamma<-matrix(gamma[1:dimzih],dimzih,1); 
    var.gamma<-as.matrix(Varbeta[1:dimzih,1:dimzih])
    robvar.gamma<-as.matrix(RVarbeta[1:dimzih,1:dimzih])
    Iinv.gamma<- as.matrix(Iinv[1:dimzih,1:dimzih])
    Ut.gamma<-as.matrix(Ut[,1:(dimzih+1)]) 
  } else { Ut.gamma<-0; gamma<-matrix(0,1,1); 
         }

  if (resample.iid==1)  {
    gamiid<-matrix(nparout[[64]],antclust,dimzih);
    biid<-matrix(nparout[[63]],Ntimes,antclust*dimxih);
    B.iid<-list();
    for (i in (0:(antclust-1))*dimxih) {
      B.iid[[i/dimxih+1]]<-as.matrix(biid[,i+(1:dimxih)]);
     # colnames(B.iid[[i/dimxih+1]])<- c(covnamesX,haploX.names); 
    }
  } else B.iid<-gamiid<-NULL; 

  ud <- list( cum = cumint, var.cum = vcum, robvar.cum = Rvcu,
             gamma = gamma, var.gamma = var.gamma, robvar.gamma = robvar.gamma,
             haplo.alpha=alpha, haplo.pars=haplo.pars,haplo.freqs=haplo.freqs, 
             var.haplo.alpha=var.hap.alpha, robvar.haplo.alpha=robvar.hap.alpha,
             var.all=Varbeta, robvar.all=RVarbeta, D2linv= Iinv.gamma, 
             D2linv.haplo.alpha= Iinv.hap.alpha, D2linv.all = Iinv, 
             resid.dMG = cumAi,
             score = score, loglike = loglike,
             pval.testBeq0 = pval.testBeq0, pval.testBeqC = pval.testBeqC,
             obs.testBeq0 = obs.testBeq0, obs.testBeqC = obs.testBeqC,
             sim.testBeq0 = sim.testBeq0, sim.testBeqC = sim.testBeqC,
             conf.band = unifCI, 
             test.procProp = Ut.gamma, sim.test.procProp = UIt,
             pval.Prop = testUt, sim.supProp = sim.supUt, 
             test.procProp.all = Ut, t=Terms, 
             haplo.design=haplo.design,
             B.iid=B.iid, gamma.iid=gamiid )
  ## }}}

  ## {{{ defining names with haplostuff
  ## best guess  for standard designs  
  haploX.names<-c()
  if (px!=dimxih) {
    haploef.x<-dimxih-px; 
    if (haploef.x>0) haploX.names<-rep("Haplo effect",haploef.x)
  }
  haploZ.names<-c()
  if (pz!=dimzih) {
    haploef.z<-dimzih-pz; 
    if (haploef.z>0) haploZ.names<-rep("Haplo effect",haploef.z)
  }
  ## }}}

  ## {{{ adding names 
  if (is.null(covnamesX)==TRUE)  covnamesXuse<-c(covnamesXin,haploX.names)  
  else covnamesXuse<-c(covnamesX)

  if (is.null(covnamesZ)==TRUE)  covnamesZuse<-c(covnamesZin,haploZ.names)  
  else covnamesZuse<-c(covnamesZ)
#print(covnamesZuse); print(covnamesXuse)

  if (length(covnamesXuse)==dimxih) {
     colnames(ud$cum)<-colnames(ud$var.cum) <-
     colnames(ud$robvar.cum) <-c("time",covnamesXuse); 
  }

  if ((length(covnamesZuse)==dimzih))  {
      names(ud$score)<-covnamesZuse;
      if (fix.beta==0) colnames(ud$test.procProp)<-c("time",covnamesZuse)
      if (resample.iid==1) colnames(gamiid)<-covnamesZuse; 
      if ((sim>=1) &  (fix.beta==0)) names(ud$pval.Prop)<-covnamesZuse
  }
  if ((sim>=1) & (length(covnamesXuse)==dimxih)) {
    names(ud$conf.band)<- names(ud$pval.testBeq0)<-
    names(ud$pval.testBeqC)<- names(ud$obs.testBeq0)<- 
    names(ud$obs.testBeqC)<- colnames(ud$sim.testBeq0)<-covnamesXuse 
  } 
  if (sim==0) {
    ud$pval.Prop<- ud$conf.band<- ud$pval.testBeq0<- ud$pval.testBeqC<-
      ud$obs.testBeq0<- ud$obs.testBeqC<- ud$sim.testBeq0<- NULL 
  }

  if ((length(covnamesZuse)==dimzih)) {
     rownames(ud$gamma)<-covnamesZuse;
     colnames(ud$gamma)<-"estimate"; 

     if (fix.beta==1) {
        namematrix(ud$var.gamma,covnamesZuse); 
        namematrix(ud$robvar.gamma,covnamesZuse); 
        namematrix(ud$D2linv,covnamesZuse); 
     }
  } 
  if (fix.beta==1) {  ud$var.gamma<-matrix(0,pz,pz); 
                      ud$robvar.gamma<-matrix(0,pz,pz);
                    }
  ## }}}

  attr(ud, "Call") <- sys.call()
  attr(ud, "Formula") <- formula
  attr(ud, "id") <- id.call
  class(ud) <- "cox.aalen"
  return(ud)
}

haplo.freqs<- function (geno.type, geno.setup=NULL,
                        Nit = 10, detail = 0, haplo.freq = NULL, step = 1,
                        lev.marq=1,min.lev.marq=0,
                        haplo.design=NULL,haplo.baseline=NULL,alpha=NULL)
{
  ## {{{ genostuff
  if (is.null(geno.setup)) {
    setup <- geno.setup(geno.type,haplo.baseline=haplo.baseline)
  } else setup<-geno.setup
  HPIordered <- setup$HPIordered
  uniqueHaplos <- setup$uniqueHaploNames
  nph<-nHaps <-length(uniqueHaplos)
  nAllelesPerLocus <- setup$nAllelesPerLocus
  unorderedAlleles <- setup$unorderedAlleles
  nPeople <- antpers<- setup$nPeople
  nLoci <- setup$nLoci
  nPossHapPairsPerPerson <- sapply(HPIordered, length)
  if (is.null(haplo.freq)) {
    haplo.freq <- as.vector(table(unlist(setup$HPIordered))/sum(table(unlist(setup$HPIordered))))
  }
  else {
    haplo.mass<-sum(haplo.freq);
    haplo.freq[length(uniqueHaplos)]<-1-sum(haplo.freq[1:(length(uniqueHaplos)-1)])
  }

  rho <- NULL
  if (is.null(rho)) {rho <- rep(0, nPeople)}

  print(nHaps); 
  print(haplo.freq[1:(nHaps - 1)]/haplo.freq[nHaps])

  haplo.pars<-log(haplo.freq[1:(nHaps - 1)]/haplo.freq[nHaps])
  print(haplo.pars)

  if (is.null(haplo.design)==TRUE)  {
    haplo.design<-diag(nph-1) 
    dimhap<-nph-1; alpha<-haplo.pars;   # full haplo-frequency model 
  } else { 
    dimhap<- ncol(haplo.design);        # haplo-freq parameters
    init.alpha<-(t(haplo.design) %*% haplo.pars) / 
      (t(haplo.design) %*% rep(1,nrow(haplo.design)) )
    if (is.null(alpha)==TRUE) alpha<-init.alpha; 
  }
  if (is.na(sum(alpha))) cat("starting value of alpha is NA"); 

  oh <- orderedHaplos <- unlist(HPIordered) -1 
  Rho <- rho
  LogLike <- numeric(1)
  ## }}}

  score <- numeric(dimhap); 
  d2lgeno<- Varhap <- Iinv <- matrix(0, dimhap, dimhap)
  score <- rep(0, dimhap)
  iid<-matrix(0,antpers,dimhap); 

  #dyn.load("genoMLE-des.so")

  out <- .C("genoMLEdes", 
            as.integer(antpers), as.integer(Nit), as.integer(detail), 
            as.integer(nph), as.integer(oh), as.integer(nPossHapPairsPerPerson), 
            as.double(haplo.pars), as.double(rho), as.double(score), 
            as.double(d2lgeno), as.double(Iinv), as.double(Varhap), 
            as.double(step), as.double(lev.marq), as.double(min.lev.marq),
            as.double(haplo.design),as.double(alpha),as.integer(dimhap),
            as.double(haplo.freq),as.double(iid),PACKAGE="HaploSurvival") 
            #gsi=as.integer(genoStringIndices),
            #nppgs=as.integer(numPeoplePerGenotypeString),
            #ngs=as.integer(length(numPeoplePerGenotypeString)),
            #as.double(X),as.integer(perms),as.integer(mp),as.double(alpha)) 

  Varhap <- matrix(out[[12]], dimhap, dimhap)
  score <- out[[9]]
  haplo.pars <- out[[7]]
  Iinv <- matrix(out[[11]], dimhap, dimhap)
  alpha<-out[[17]]
  iid<-matrix(out[[20]],antpers,dimhap)

  haplo.pars<-haplo.design %*% alpha; 
  haplo.freqs = exp(c(haplo.pars,0))/sum(exp(c(haplo.pars, 0)))

  out<-list(haplo.alpha=alpha,haplo.pars=haplo.pars,
            haplo.freq=haplo.freqs, 
            var.haplo.alpha=Varhap,D2linv=Iinv,score=score,
            haplo.design=haplo.design,alpha.iid=iid,
            freq.names=setup$uniqueHaploNames)

  attr(out, "Call") <- sys.call()
  class(out) <- "haplo.freqs"
  return(out)
}

print.haplo.freqs <- function(x,...){

  cat("A haplo.freqs object:",fill = TRUE)
  cat("  Data for ",nrow(x$alpha.iid)," people",sep = "",fill = TRUE)
  cat("  Data for ",length(unlist(strsplit(x$freq.names[1],',')))," loci",sep = "",fill = TRUE)
  cat("  ",length(x$haplo.freq)," possible haplotypes",sep = "",fill = TRUE)
  names(x$haplo.freq) <- x$freq.names
  x$haplo.freq <- round(x$haplo.freq,3)
  cat("  MLE haplotype freqs: ",sep = "",fill = TRUE)
  maxNchar <- max(c(nchar(head(x$freq.names)),5))
  numSprintfArgs <- min(length(x$haplo.freq),6)
  
  cat('    ')
  for(i in 1:numSprintfArgs){
    cat(sprintf(paste(c("%",maxNchar,"s "),collapse=''),x$freq.names[i]))
  }
  if(length(x$haplo.freq) > 6){
    cat('...')
  }
  cat('',fill = TRUE)

  cat('    ')
  for(i in 1:numSprintfArgs){
    cat(sprintf(paste(c("%",maxNchar,"g "),collapse=''),x$haplo.freq[i]))
  }
  if(length(x$haplo.freq) > 6){
    cat('...')
  }
  cat('',fill = TRUE)
  
}

summary.haplo.freqs <- function(object,...){
  print.haplo.freqs(object,...)
}
