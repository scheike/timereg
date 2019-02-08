##' CIF regression
##'
##' CIF logistic for propodds=1 default
##' CIF Fine-Gray (cloglog) regression for propodds=NULL
##'
##'
##' @param formula formula with 'Event' outcome 
##' @param data data frame
##' @param cause of interest 
##' @param cens.code code of censoring 
##' @param offset offsets for cox model
##' @param weights weights for Cox score equations
##' @param Gc censoring weights for time argument, default is to calculate these with a Kaplan-Meier estimator 
##' @param propodds 1 is logistic model, NULL is fine-gray model 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @examples
##' ## data with no ties
##' data(bmt,package="timereg")
##' bmt$time <- bmt$time+runif(nrow(bmt))*0.01
##' bmt$id <- 1:nrow(bmt)
##' 
##' ## logistic link  OR interpretation
##' ll=cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1)
##' bplot(ll)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pll <- predict(ll,nd)
##' plot(pll)
##' 
##' ## Fine-Gray model
##' llfg=cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,
##'                   propodds=NULL)
##' bplot(ll)
##' nd <- data.frame(tcell=c(1,0),platelet=0,age=0)
##' pll <- predict(ll,nd)
##' plot(pll)
##'
##' @export
cifreg <- function(formula,data=data,cause=1,cens.code=0,
			weights=NULL,offset=NULL,Gc=NULL,propodds=1,...)
{# {{{

  cl <- match.call()# {{{
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (class(Y)!="Event") stop("Expected a 'Event'-object")
  if (ncol(Y)==2) {
    exit <- Y[,1]
    entry <- NULL ## rep(0,nrow(Y))
    status <- Y[,2]
  } else {
    entry <- Y[,1]
    exit <- Y[,2]
    status <- Y[,3]
  }
  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    pos.cluster <- ts$terms
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  } else pos.cluster <- NULL
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms  <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  if (!is.null(id)) {
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(exit))-1; 
  ### id from call coded as numeric 1 -> 
  id.orig <- id; 

# }}}

 res <- c(cifreg01(data,X,exit,status,id,strata,offset,weights,strata.name,
		   cause=cause,cens.code=cens.code,Gc=Gc,propodds=propodds,...),
   list(call=cl,model.frame=m,
	formula=formula,
	strata.pos=pos.strata,cluster.pos=pos.cluster)
   )

  class(res) <- c("phreg","cif.reg")
  return(res)
}# }}}

cifreg01 <- function(data,X,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
             strata.name=NULL,beta,stderr=TRUE,
	     method="NR",no.opt=FALSE,propodds=1,profile=0,
	     case.weights=NULL,cause=1,cens.code=0,Gc=NULL,...) {# {{{
##  setting up weights, strata, beta and so forth before the action starts# {{{
 p <- ncol(X)
 if (missing(beta)) beta <- rep(0,p)
 if (p==0) X <- cbind(rep(0,length(exit)))

 cause.jumps <- which(status==cause)
 max.jump <- max(exit[cause.jumps])
 other <- which((!(status %in% c(cens.code,cause)) ) & (exit< max.jump))

 entry <- NULL
 n <- length(exit)
 trunc <- (!is.null(entry))
  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }
 if (!trunc) entry <- rep(0,length(exit))
 if (is.null(offset)) offset <- rep(0,length(exit)) 
 if (is.null(weights)) weights <- rep(1,length(exit)) 
 if (is.null(case.weights)) case.weights <- rep(1,length(exit)) 
 strata.call <- strata
 Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 

  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1; 
   ## orginal id coding into integers 
   id.orig <- id+1; 
# }}}

 ### censoring weights constructed
 whereC <- which(status==cens.code)
 time <- exit
 if (is.null(Gc)) {
 cens.model <- km(Surv(exit,status==cens.code)~+1,data=data)
 wpredS <- fast.approx(c(0,cens.model$time),exit,type="left")
### wpredS <- timereg:::sindex.prodlim(c(0,cens.model$time),exit)
 Stime <- c(1,cens.model$surv)[wpredS]
 } else { 
        if (length(whereC)>0) Ctimes <- sort(unique(exit[whereC]))
	else Ctimes <- 0
        Stime <- Gc
 }

 ## setting up all jumps of type "cause", need S0, S1, S2 at jumps of "cause"
 stat1 <- 1*(status==cause)
 xx2 <- .Call("FastCoxPrepStrata",entry,exit,stat1,X,id, 
	     trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets")
 xx2$nstrata <- nstrata
 jumps <- xx2$jumps+1
 jumptimes <- xx2$time[jumps]
 Xj <- xx2$X[jumps,,drop=FALSE]

 ## G(T_j-)
 if (length(whereC)>0) {
	 whereJ <- fast.approx(c(0,cens.model$time),jumptimes,type="left")
	 Gjumps <- c(1,cens.model$surv)[whereJ]
 } else {
         Gjumps <- rep(1,length(jumptimes))
 }

 ## computing terms for those experiencing another cause, need S0, S1, S2
 if ( length(other)>=1) {# {{{
	 trunc <- TRUE
	 weightso <- 1/Stime[other]
	 timeoo <- rep(max(exit)+1,length(other))
	 statuso <- rep(0,length(other))
	 Xo <- X[other,,drop=FALSE]
	 offseto <- offset[other]
	 entryo <- exit[other]
	 ido <- id[other]
	 stratao <- strata[other]
	 ###
	 xx <- .Call("FastCoxPrepStrata",entryo,timeoo,statuso,Xo,
	     ido,trunc,stratao,weightso,offseto,Zcall,case.weights[other],PACKAGE="mets")
	 xx$nstrata <- nstrata

	 timeo  <- xx$time
	 ## use right here because T_jump is larger than the T_(other) som står i listen
	 ## timeo
	 where <- fast.approx(c(0,timeo),jumptimes,type="right")
 }# }}}

obj <- function(pp,all=FALSE) {# {{{

if (length(other)>=1) {
	rr <- c(xx$sign*exp(xx$X %*% pp + xx$offset)*xx$weights)
	S0no <- revcumsumstrata(rr,xx$strata,xx$nstrata)
	S1no  <- apply(xx$X*rr,2,revcumsumstrata,xx$strata,xx$nstrata); 
	S2no  <- apply(xx$XX*rr,2,revcumsumstrata,xx$strata,xx$nstrata); 

	## look at jumptimes
	S0no <- c(0,S0no)[where]
	S1no <- rbind(0,S1no)[where,,drop=FALSE]
	S2no <- rbind(0,S2no)[where,,drop=FALSE]
} else { Gjumps <- S0no <- S1no <-  S2no <- 0} 

rr2 <- c(xx2$sign*exp(xx2$X %*% pp + xx2$offset)*xx2$weights)
rr2now <- c(xx2$sign*exp(xx2$X %*% pp + xx2$offset))
S0oo <- revcumsumstrata(rr2,xx2$strata,xx2$nstrata)
S1oo  <- apply(xx2$X*rr2,2,revcumsumstrata,xx2$strata,xx2$nstrata); 
S2oo  <- apply(xx2$XX*rr2,2,revcumsumstrata,xx2$strata,xx2$nstrata); 
S0oo <- S0oo[jumps,]
S1oo <- S1oo[jumps,,drop=FALSE]
S2oo <- S2oo[jumps,,drop=FALSE]

S0 <- c(S0oo+S0no*Gjumps)
E <-   (S1oo+S1no*Gjumps)/S0
weightsJ <- xx2$weights[jumps]
strataJ <- xx2$strata[jumps]
rr2now <- rr2now[jumps]
U <- (Xj-E)
ploglik <- (log(rr2now)-log(S0))*weightsJ; 

if (!is.null(propodds)) {
   strataJ <- xx2$strata[jumps]
   pow <- c(.Call("cumsumstrataPOR",weightsJ,S0,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$pow); 
   DLam <-.Call("DLambetaR",weightsJ,S0,E,Xj,strataJ,nstrata,propodds,rr2now,PACKAGE="mets")$res; 
   Dwbeta <- DLam*rr2now+(pow-1)*Xj
   DUadj  <- .Call("vecMatMat",Dwbeta,U,PACKAGE="mets")$vXZ 
}

Ut <- weightsJ*U
## E^2, as n x (pxp)
Et2 <-  .Call("vecMatMat",E,E,PACKAGE="mets")$vXZ
S2S0 <- (S2oo+S2no*Gjumps)/S0
DUt <-  -(S2S0-Et2)

if (!is.null(propodds)) {
	Ut  <- pow*Ut
	S0 <- S0/pow
	DUt <- pow*DUt
	DUt <- DUt+DUadj
	if (profile==1) {
		Ut <- Ut+c(ploglik)*Dwbeta
		## not implemented 
		DUt <- DUt 
	}
	ploglik <- pow*ploglik
}

U  <- apply(Ut,2,sum)
DUt <- weightsJ*DUt
DU <- -matrix(apply(DUt,2,sum),p,p)

ploglik <- sum(ploglik)

out <- list(ploglik=ploglik,gradient=U,hessian=-DU,cox.prep=xx2,
	    hessiantime=DUt,weightsJ=weightsJ,
	    time=jumptimes,S0=S0/weightsJ,S2S0=S2S0,E=E,U=Ut,X=Xj,Gjumps=Gjumps)

if (all) return(out) else with(out,structure(-ploglik, gradient=-gradient, hessian=-hessian))
}# }}}

 if (length(jumps)==0) no.opt <- TRUE

 opt <- NULL
  if (p>0) {# {{{
  if (no.opt==FALSE) {
      if (tolower(method)=="nr") {
          opt <- lava::NR(beta,obj,...)
          opt$estimate <- opt$par
      } else {
          opt <- nlm(obj,beta,...)
          opt$method <- "nlm"
      }
      cc <- opt$estimate;  names(cc) <- colnames(X)
      if (!stderr) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
      } else val <- c(list(coef=beta),obj(beta,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
  }# }}}

### opt <- lava::NR(beta,obj); beta.s <- opt$par 
 beta.s <- val$coef 
 ## getting final S's 
 opt <-  obj(beta.s,all=TRUE)

 ### iid version 
# {{{
 ###  iid.phreg
 ##iid robust phreg
 Gt <- S0i <- rep(0,length(xx2$strata))
 S0i[jumps] <- 1/opt$S0
 Z <- xx2$X
 U <- E <- matrix(0,nrow(Z),p)
 E[jumps,] <- opt$E
 U[jumps,] <- opt$U
 Gt[jumps] <- Gjumps
 ###
 cumhazA <- cumsumstratasum(S0i,xx2$strata,xx2$nstrata,type="all")
 cumhaz <- c(cumhazA$sum)
 rr <- c(xx2$sign*exp(Z %*% beta.s + xx2$offset))
 if (!is.null(propodds)) {
    cumhazm <- c(cumhazA$lagsum)
    S0star <- cumsumstrata(rr/(1+rr*cumhazm),xx2$strata,xx2$nstrata)
 }
 EdLam0 <- apply(E*S0i,2,cumsumstrata,xx2$strata,xx2$nstrata)

 ### Martingale  as a function of time and for all subjects to handle strata 
 MGt <- U[,drop=FALSE]-(Z*cumhaz-EdLam0)*rr*c(xx2$weights)
 mid <- max(xx2$id)
 UU <- apply(MGt,2,sumstrata,xx2$id,mid+1)

 if (length(other)>=1) {
	 ### T_j for jumps of other type 
         where <- fast.approx(c(0,xx2$time),entryo,type="right")
	 ###
	 rcumhazGt <- c(revcumsumstrata(S0i*Gt,xx2$strata,xx2$nstrata))
	 rEdLam0Gt <- apply(E*S0i*Gt,2,revcumsumstrata,xx2$strata,xx2$nstrata)

         ###
	 rcumhazGtx <- c(rcumhazGt[1],rcumhazGt)[where]
	 rEdLam0Gtx <- rbind(rEdLam0Gt[1],rEdLam0Gt)[where,]  
	 ###
	 rrx <- c(exp(Xo %*% beta.s + offseto)*weightso)

         if (!is.null(propodds)) {

	 }
	 ###	
	 MGtGtx <- -(Xo*rcumhazGtx-rEdLam0Gtx)*rrx

	 UU2 <- apply(MGtGtx,2,sumstrata,ido,mid+1)
	 UU  <-  UU+UU2
 }
# }}}


if ((length(other)>=1) & (length(whereC)>0)) {
 ### Censoring adjustment for jumps of other type {{{

 where <- fast.approx(xx2$time,entryo,type="right")
 rrrx <- rep(0,length(xx2$strata))
 rrrx[where] <- rrx
 Xos <- matrix(0,length(xx2$time),ncol(Xo));
 Xos[where,] <- Xo
 ###
 Xos <- apply(Xos,2,cumsum)
 rro <- cumsum(rrrx)
 ###
 q <- -(Xos*rcumhazGt-rEdLam0Gt*rro)

 cens.mgs = phreg(Surv(exit,status==0)~+cluster(id),data=data,no.opt=TRUE)
 cxx <- cens.mgs$cox.prep
 ###
 Gt <- S0i <-  S0i2 <- rep(0,length(cxx$strata))
 S0i[cxx$jumps+1] <- 1/cens.mgs$S0
 S0i2[cxx$jumps+1] <- 1/cens.mgs$S0^2
 qc <- matrix(0,nrow(q),ncol(q)) 
 ## sort q after censoring times
 qc[cxx$jumps+1,] <- q[cxx$jumps+1]
 ###
 EdLam0q <- apply(qc*S0i2,2,cumsumstrata,cxx$strata,cxx$nstrata)
 ### Martingale  as a function of time and for all subjects to handle strata 
 MGc <- qc[,drop=FALSE]*S0i-EdLam0q
 MGc <- apply(MGc,2,sumstrata,cxx$id,mid+1)
# }}}
} else MGc <- 0

 iH <- - tryCatch(solve(opt$hessian),error=
	 function(e) matrix(0,nrow(opt$hessian),ncol(opt$hessian)) )
 Uiid <-  (UU+MGc) %*% iH
 UUiid <- UU %*% iH

 var1 <-  crossprod(UUiid) 
 varm <-  crossprod(Uiid) 

 strata <- xx2$strata[jumps]
 cumhaz <- cbind(opt$time,cumsumstrata(1/opt$S0,strata,nstrata))
 colnames(cumhaz)    <- c("time","cumhaz")

out <- list(coef=beta.s,var=varm,se.coef=diag(varm)^.5,UUiid=UUiid,Uiid=Uiid,
	    ihessian=iH,hessian=opt$hessian,var1=var1,se1.coef=diag(var1)^.5,
	    ploglik=opt$ploglik,gradient=opt$gradient,
	    cumhaz=cumhaz,strata=xx2$strata,nstrata=nstrata,strata.name=strata.name,
	    strata.level=strata.level,propodds=propodds,
	    S0=opt$S0,E=opt$E,S2S0=opt$S2S0,time=opt$time,
            jumps=jumps,II=iH,exit=exit,p=p,opt=opt,n=nrow(X),nevent=length(jumps)
	    )

return(out)
}# }}}

