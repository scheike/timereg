additive.plsR<-function (formula = formula(data), ## {{{
                         data = sys.parent(), start.time = 0, max.time = NULL, 
                         id=NULL, clusters=NULL,   
                         pls.dim=1, semi.pls=1,silent=0){
## {{{ setting up variables
  call <- match.call()
  m <- match.call(expand = FALSE)
  m$start.time <-  m$max.time <- 
  m$id <- m$clusters <- m$pls.dim<- m$semipls<- m$silent<- NULL
  special <- c("const","cluster")
  Terms <- if (missing(data)) terms(formula, special)
  else terms(formula, special, data = data)
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

  if(is.null(clusters)) clusters <- des$clusters  
  
  pxz <- px + pz; 
  npar.pls<-npar; 

  survs<-read.surv(m,id,FALSE,clusters,start.time,max.time)
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; time2<-survs$stop
  status<-survs$status; 
  ldata<-list(start=survs$start,stop=survs$stop,
              antpers=survs$antpers,antclust=survs$antclust);

  px<-px-1; 
  beta<-matrix(0,pls.dim,px); pls.comp<-matrix(0,nrow(X),pls.dim); 

  Yorig<-as.matrix(X[,1]); 
  pls.cov<-as.matrix(X[,-1]); Zorig<-Z; 
  Xorig<-as.matrix(X[,-1]); 
  ## }}}

  pvals<-var.gamma<-c()
  for (j in 1:pls.dim)
    {
      cat(paste("PLS component ",j,"\n")) 
      for (i in 1:px)
        {
          if (npar==TRUE ) {
            if (j==1) Z<-as.matrix(Xorig[,i]) else Z<-
              cbind(Xorig[,i],pls.comp[,1:(j-1)]); }
          if (npar==FALSE) {
            if (j==1) Z<-cbind(Xorig[,i],Zorig) else 
            Z<-cbind(Xorig[,i],pls.comp[,1:(j-1)],Zorig); }

        ud <- semiaalen(times, ldata, Yorig , Z, 
        status, id , clusters, 
	robust = 0, sim = 0, antsim = 0, weighted.test = 1, retur = 0,
        covariance=0, resample.iid=0,namesX=covnamesX,namesZ=covnamesZ,
                          deltaweight=1,silent=silent)

          pvals<-c(pvals,1-pchisq(ud$gamma[1]^2/ud$var.gamma,1))
          var.gamma<-c(var.gamma,ud$var.gamma)

          beta[j,i]<-c(ud$gamma)[1];
          pls.comp[,j]<-pls.comp[,j]+beta[j,i]*Xorig[,i]; 
        } 
    }
  cnames<-c()
  for (j in 1:pls.dim) cnames<-c(cnames,paste("PLS-",j)); 
  colnames(pls.comp)<-cnames

  if (npar==TRUE) Z<-pls.comp else Z<-cbind(pls.comp,Zorig)
  udf <- semiaalen(times, ldata, Yorig , Z, 
         status, id , clusters, 
	robust = 0, sim = 0, antsim = 0, weighted.test = 1, retur = 0,
        covariance=0, resample.iid=0,namesX=covnamesX,namesZ=covnamesZ,
                          deltaweight=1,silent=silent)

ud<-list(pls.comp=pls.comp,beta.pls=beta,
         beta=udf$gamma,baseline=udf$cum,pvals=pvals,var.gamma=var.gamma)
class(ud) <- "pls"
attr(ud, "Call") <- sys.call()
attr(ud, "Formula") <- formula
ud$call <- call
return(ud)
} ## }}}

additive.pls<-function (formula = formula(data), data = sys.parent(), 
start.time=0,max.time=NULL,id=NULL,pls.dim=1,silent=1) 
{  ## {{{
## {{{ setting up variables
semipls=1; deltaweight<-1; constant<-1; # always !  
weighted.pls=0; 
call <- match.call()
    m <- match.call(expand = FALSE)
    m$start.time <- m$max.time <-  m$id <- m$pls.dim<- m$weighted.pls<- m$silent<- NULL
    special <- c("const")
    Terms <- if (missing(data)) terms(formula, special)
    else terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

    des<-read.design(m,Terms)
    X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
    npar <- des$npar; 
    covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
    pxz <- px + pz; 

    clusters=NULL;
    if (pls.dim>=1) survs<-read.surv(m,id,FALSE,clusters,start.time,max.time)
    else survs<-read.surv(m,id,TRUE,clusters,start.time,max.time)
    times<-survs$times;id<-id.call<-survs$id.cal;
    clusters<-cluster.call<-survs$clusters; time2<-survs$stop
    status<-survs$status; 

  nobs <- nrow(X); 
  weights <- rep(1,nrow(X)); 

if ( (attr(m[, 1], "type") == "right" ) ) {  ## {{{
   ot<-order(-time2,status==1); # order in time, status=0 first for ties
   time2<-time2[ot]; status<-status[ot]; 
   stop <- time2; 
   start <- rep(0,nobs); 
   X<-as.matrix(X[ot,])
   if (npar==FALSE) Z<-as.matrix(Z[ot,])
   clusters<-clusters[ot]
   id<-id[ot];
   weights <- weights[ot]; 
   entry=rep(-1,nobs); 
  } else {
        stop("Left truncation not implemented, use additive.plsR"); 
        eventtms <- c(survs$start,time2)
        start <- c(survs$start,survs$start)
        status <- c(rep(0, nobs), status)
        ix <- order(-eventtms,status==1)
        etimes    <- eventtms[ix]  # Entry/exit times
        stop <- etimes 
	status <- status[ix]
	start <- start[ix]
        tdiff    <- c(-diff(etimes),start.time) # Event time differences
        entry  <- c(rep(c(1, -1), each = nobs))[ix]
        weights <- rep(weights, 2)[ix]
        X        <- X[rep(1:nobs, 2)[ix],]
	if (npar==FALSE) Z <- Z[rep(1:nobs,2)[ix],]
	id <- rep(id,2)[ix]
	clusters <- rep(clusters,2)[ix]
    } ## }}}

    px<-px-1; 
    beta<-matrix(0,pls.dim,px); 
    pls.comp<-matrix(0,nrow(X),pls.dim); 

    Yorig<-as.matrix(X[,1]); colnames(Yorig)<-"Intercept"
    pls.cov<-as.matrix(X[,-1]); 
    Zorig<-Z;  Xorig<-as.matrix(X[,-1]); 
 ## }}}

    if (pls.dim>0) { 
    for (j in 1:pls.dim)
    {## {{{

    Xmat<-as.matrix(Yorig); 
    covnamesXo<-colnames(Xmat); 
    if (npar==TRUE ) {
    if (j==1) Z<-as.matrix(pls.cov[,1]) else  {
    if (constant==1)  {Z<-cbind(pls.comp[,1:(j-1)],pls.cov[,1]);} else
                      {Z<-as.matrix(pls.cov[,1]); 
                       Xmat<-as.matrix(cbind(Yorig,pls.comp[,1:(j-1)])); 
                       } } }
    if (npar==FALSE) {
    Xmat<-cbind(Yorig,Zorig); 
    covnamesXo<-colnames(Xmat); 
    if (j==1) Z<-as.matrix(pls.cov[,1]) else  {

    if (constant==1)  {Z<-cbind(pls.comp[,1:(j-1)],pls.cov[,1]);} else
                      {Z<-as.matrix(cbind(Zorig,pls.cov[,1])); 
                       Xmat<-as.matrix(cbind(Yorig,pls.comp[,1:(j-1)]));}}}

    Z<-as.matrix(Z); Xmat<-as.matrix(Xmat); 
    pg<-ncol(as.matrix(Z)); px<-ncol(as.matrix(Xmat)); nx<-nrow(Yorig); 

    dimplscov<-ncol(pls.cov);
    betapls<-rep(0,dimplscov); plscomp<-Yorig[,1]; 

    Nalltimes <- length(times);
    Ntimes<-sum(survs$status[(survs$stop>times[1]) & (survs$stop<=times[Nalltimes])])+1;
    cum<-matrix(0,Ntimes,px);  

###    dyn.load("pls.so"); 

    semiout<-.C("plssemiaddN",
    as.double(times),as.integer(Nalltimes),as.integer(Ntimes),
    as.double(Xmat),as.integer(nx),as.integer(px),
    as.double(Z),as.integer(nx),as.integer(pg),
    as.integer(survs$antpers),as.double(start),as.double(stop),
    as.integer(id),as.double(cum),as.integer(status),
    as.integer(deltaweight),as.double(pls.cov),
    as.integer(dimplscov),as.double(betapls),as.double(plscomp),
    as.integer(semipls),as.integer(weighted.pls),as.integer(silent),as.double(weights),
    as.integer(entry), PACKAGE="timereg")

    plscomp<-semiout[[20]]; 
    betapls<-semiout[[19]];
    beta[j,]<-betapls; 
    pls.comp[,j]<-pls.comp[,j]+ Xorig  %*%  beta[j,]; 
    }
    cnames<-c()
    for (j in 1:pls.dim) cnames<-c(cnames,paste("PLS-",j)); 
    colnames(beta)<-covnamesX[-1]; 
    Z<-as.matrix(pls.comp) 

   if ( (attr(m[, 1], "type") == "right" ) )
   udf <-  aalen(Surv(stop,status)~-1+Xmat+const(Z),robust=0,max.time=max.time)
   else udf <-  aalen(Surv(start,stop,status)~-1+Xmat+const(Z),robust=0,max.time=max.time)

    betaf<-c(udf$gamma)
    tbeta.pls<-apply(beta*betaf[1:pls.dim],2,sum)
    colnames(udf$cum)<-c("times",covnamesXo)
} else {  ## {{{ pls.dim=0
    betaf<-rep(0,1); pls.comp<-matrix(0,nrow(Yorig),1); 
    names(betaf)<-"dummy"; 
    tbeta.pls<-rep(0,px); 
    names(tbeta.pls)<-covnamesX[-1]; 

if ( (attr(m[, 1], "type") == "right" ) )
 udf <-  aalen(Surv(stop,status)~-1+Yorig,robust=0,weights=weights,max.time=max.time)
else udf <-  aalen(Surv(start,stop,status)~-1+Yorig,robust=0,weights=weights,max.time=max.time)
colnames(udf$cum)<-c("times","Intercept")
} ## }}}


ud<-list(pls.comp=pls.comp,beta.pls=beta,beta=betaf,
baseline=udf$cum,tbeta.pls=tbeta.pls)
class(ud) <- "pls"
attr(ud, "Call") <- sys.call()
attr(ud, "Formula") <- formula
ud$call <- call

return(ud)
} ## }}}

predict.pls<-function(object,Z=NULL,X=NULL,times=NULL,monotone=TRUE, ...)
{ ## {{{
if(!(inherits(object,"pls")))stop("Must be output from additive.pls function")

  ## Allow for the possibility that X is a vector
  if (is.null(Z))  stop("Must specify Z, and baseline covariates X if needed \n");
  if(!is.matrix(Z) && !is.data.frame(Z)){
    Z<-matrix(Z,ncol = ncol(object$beta.pls)); }
  nobs<-nrow(Z);
  if (is.null(X)) X <- matrix(1,nobs,1);
  if(!is.matrix(X) && !is.data.frame(X)){
    X<-matrix(X,ncol = ncol(object$baseline)-1); } else { X<-X; }

  risk<- object$tbeta.pls %*% t(Z)
  if (is.null(times)) {times<-object$baseline[,1]; cumbaselinet<-object$baseline } else
  cumbaselinet<-Cpred(object$baseline,times)

  cumhaz <- as.matrix(X) %*% t(cumbaselinet[, -1])
  surv<-as.matrix(exp(-cumhaz-t(rbind(c(risk)) %x% cbind(c(times)))))

  if (monotone==TRUE) {
     survm<- -t(apply(-surv,1,pava));
     survm[survm>1]<-1; survm[survm<0]<-0; } else survm<-surv; 

return(list(times=times,surv=survm))
} ## }}}

summary.pls<-function(object, ...)
{ ## {{{
if (!(inherits(object,"pls")))stop("Must be output from additive.pls function")
} ## }}}

print.pls<-function(x, ...)
{ ## {{{
if(!(inherits(x,"pls")))stop("Must be output from additive.pls function")
} ## }}}
