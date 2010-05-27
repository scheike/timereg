
aalen.test<-function (formula = formula(data),
data = sys.parent(), start.time = 0, max.time = NULL, 
robust=1, id=NULL, clusters=NULL, residuals = 0, n.sim = 1000,  
weighted.test=0,covariance=0,resample.iid=0,weights=0,offsets=0,
fix.gam=0,pseudo.score=0,approx="dt",gamma=0,silent=0)
{
  deltaweight<-1; 
  if (n.sim == 0) sim <- 0 else sim <- 1
  if (resample.iid==1 & robust==0) {
    cat("When robust=0 no iid representaion computed\n"); 
    resample.iid<-0;}
  if (covariance==1 & robust==0) {
    cat("When robust=0 no covariance computed \n"); 
    cat("Covariance based on robust iid representation\n")
    covariance<-0;}
  if (sim==1 & robust==0) {
    cat("When robust=0, No simulations \n"); 
    cat("n.sim set to 0\n"); n.sim<-0;}
  if (resample.iid==0 & pseudo.score>=1) {
    cat("When pseudo.score>=1, resample.iid set to 1 \n"); 
    resample.iid<-1;}
  if (residuals==1 & robust==0) {
    cat("When robust=0, no martingale residuals \n"); 
    residuals<-0;}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  call <- match.call()
  m <- match.call(expand = FALSE)
  m$pseudo.score<-m$start.time <- m$weighted.test <- m$max.time <-
    m$robust <- m$sim <- m$residuals <- m$n.sim <- m$id <-
      m$fix.gam<- m$covariance <- m$resample.iid <- m$clusters <- 
        m$weights<-m$offsets<- m$gamma<-m$deltaweight<-m$approx<-m$silent<-NULL
  special <- c("const","cluster")
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

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ

  if(is.null(clusters)) clusters <- des$clusters  

  pxz <- px + pz;

  if (approx=="death-times") npar.call<-TRUE else npar.call<-FALSE
  survs<-read.surv(m,id,npar.call,clusters,start.time,max.time)
  times<-survs$times;id<-id.call<-survs$id.cal;
  if (approx=="death-times" & des$npar==FALSE) npar.call<-FALSE; 

  id<-id.call<-survs$id; 
  clusters<-cluster.call<-survs$clusters; 
  time2 <- survs$stop;  status <- survs$status
  ldata <- list(start = survs$start, stop = survs$stop,
                antpers = survs$antpers, antclust = survs$antclust)

  if (npar== TRUE) {
                                        #cat("Nonparametric Additive Risk Model\n")
    ud <- aalenBaseTest(times, ldata, X, status, id, clusters, robust = robust, 
                        sim = sim, retur = residuals, antsim = n.sim,
                        weighted.test = weighted.test,covariance=covariance,
                        resample.iid=resample.iid,namesX=covnamesX,weights=weights,
                        offsets=offsets,silent=silent)
    colnames(ud$cum) <- colnames(ud$var.cum) <- c("time", 
                                                  covnamesX)
    if (robust == 1) 
      colnames(ud$robvar.cum) <- c("time", covnamesX)
    if (sim >= 1) {
      colnames(ud$test.procBeqC) <- c("time", covnamesX)
      names(ud$conf.band) <- names(ud$pval.testBeq0) <- names(ud$pval.testBeqC) <- names(ud$pval.testBeqC.is) <- names(ud$obs.testBeq0) <- names(ud$obs.testBeqC) <- names(ud$obs.testBeqC.is) <- colnames(ud$sim.testBeq0) <- colnames(ud$sim.testBeqC) <- colnames(ud$sim.testBeqC.is) <- covnamesX
      ud$sim.testBeqC.is <- ud$sim.testBeqC <- FALSE
    }
  }
  else {
                                        #   cat("Semiparametric Additive Risk Model\n")
    if (px == 0) 
      stop("No nonparametric terms (needs one!)")
    ud <- semiaalenTest(times, ldata, X, Z, 
                        status, id , clusters, robust = robust, sim = sim, antsim = n.sim, 
                        weighted.test = weighted.test, retur =
                        residuals,covariance=covariance,
                        resample.iid=resample.iid,namesX=covnamesX,namesZ=covnamesZ,
                        deltaweight=deltaweight,weights=weights,offsets=offsets,
                        fixgam=fix.gam,gamma=gamma,pseudo.score=pseudo.score,
                        silent=silent)
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
    if (pseudo.score>=1) {
      names(ud$pstest.pval)<- names(ud$sup.pscore)<-c(covnamesZ); 
      colnames(ud$obs.pscore)<-c("time",covnamesZ); 
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
