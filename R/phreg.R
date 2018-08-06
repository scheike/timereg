###{{{ phreg0 

phreg0 <- function(X,entry,exit,status,id=NULL,strata=NULL,beta,stderr=TRUE,method="NR",...) {# {{{
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(exit)))
  if (!is.null(strata)) {# {{{
    stratalev <- levels(strata)
    strataidx <- lapply(stratalev,function(x) which(strata==x))
    if (!all(unlist(lapply(strataidx,function(x) length(x)>0))))
      stop("Strata without any observation")
    dd <- lapply(strataidx, function(ii) {
        entryi <- entry[ii]
        trunc <- !is.null(entryi)
        if (!trunc) entryi <- rep(0,length(exit[ii]))
                 .Call("FastCoxPrep",
                       entryi,exit[ii],status[ii],
                       as.matrix(X)[ii,,drop=FALSE],
                       id[ii],
                       trunc,
                       PACKAGE="mets")
                 })
    if (!is.null(id))
      id <- unlist(lapply(dd,function(x) x$id[x$jumps+1]))
      obj <- function(pp,U=FALSE,all=FALSE) {
      val <- lapply(dd,function(d)
                    with(d,
                         .Call("FastCoxPL",pp,X,XX,sign,jumps,PACKAGE="mets")))
      ploglik <-     Reduce("+",lapply(val,function(x) x$ploglik))
      gradient <-    Reduce("+",lapply(val,function(x) x$gradient))
      hessian <-     Reduce("+",lapply(val,function(x) x$hessian))
      if (all) {
        U <- do.call("rbind",lapply(val,function(x) x$U))
        hessiantime <- do.call("rbind",lapply(val,function(x) x$hessianttime))
        time <- lapply(dd,function(x) x$time[x$ord+1])
        ord <- lapply(dd,function(x) x$ord+1)
        jumps <- lapply(dd,function(x) x$jumps+1)
        jumptimes <- lapply(dd,function(x) x$time[x$ord+1][x$jumps+1])
        S0 <- lapply(val,function(x) x$S0)
        nevent  <- unlist(lapply(S0,length))
        return(list(ploglik=ploglik,gradient=gradient,hessian=hessian,
                    U=U,S0=S0,nevent=nevent,hessianttime=hessiantime,
                    ord=ord,time=time,jumps=jumps,jumptimes=jumptimes))
      }
      structure(-ploglik,gradient=-gradient,hessian=-hessian)
    }# }}}
  } else {# {{{
      trunc <- !is.null(entry)
      if (!trunc) entry <- rep(0,length(exit))
      system.time(dd <- .Call("FastCoxPrep",
                              entry,exit,status,X,
                              as.integer(seq_along(entry)),
                              !is.null(entry),
                              PACKAGE="mets"))

      if (!is.null(id))
          id <- dd$id[dd$jumps+1]
      obj <- function(pp,U=FALSE,all=FALSE) {
          val <- with(dd,
                      .Call("FastCoxPL",pp,X,XX,sign,jumps,PACKAGE="mets"))
          if (all) {
              val$time <- dd$time[dd$ord+1]
              val$ord <- dd$ord+1
              val$jumps <- dd$jumps+1
              val$jumptimes <- val$time[val$jumps]
              val$nevent <- length(val$S0)
              return(val)
          }
          with(val, structure(-ploglik,gradient=-gradient,hessian=-hessian))
      }
  }# }}}
  opt <- NULL
  if (p>0) {
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
  } else {
      val <- obj(0,all=TRUE)
      val[c("ploglik","gradient","hessian","U")] <- NULL
  }

  ### computes Breslow estimator 
  cumhaz <- NULL

  res <- c(val,
           list(strata=strata,
                entry=entry,
                exit=exit,
                status=status,                
                p=p,
                X=X,
                id=id, opt=opt,cum=cumhaz))
  class(res) <- "phreg"
  res
} # }}}

###}}} phreg0

###{{{ phreg01

phreg01 <- function(X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
		    strata.name=NULL,cumhaz=TRUE,
  beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,...) {
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(exit)))
  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1; strata.level <- NULL; } else {
	  strata.level <- levels(strata)
	  ustrata <- sort(unique(strata))
	  nstrata <- length(ustrata)
	  strata.values <- ustrata
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }
  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 
  strata.call <- strata
  Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
  if (!is.null(Z)) Zcall <- Z

  trunc <- (!is.null(entry))
  if (!trunc) entry <- rep(0,length(exit))

  id.orig <- id; 
  if (!is.null(id)) {
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1; 

   system.time(dd <- .Call("FastCoxPrepStrata",
		     entry,exit,status,X, id, ### as.integer(seq_along(entry)),
		     trunc,strata,weights,offset,Zcall,PACKAGE="mets"))

	 print("prep her"); 

   dd$nstrata <- nstrata
	obj <- function(pp,U=FALSE,all=FALSE) {# {{{
		if (is.null(propodds) & is.null(AddGam)) 
	  val <- with(dd,
		   .Call("FastCoxPLstrata",pp,X,XX,sign,jumps,
		    strata,nstrata,weights,offset,ZX,PACKAGE="mets"))
         else if (is.null(AddGam)) 
		 val <- with(dd,
		   .Call("FastCoxPLstrataPO",pp,X,XX,sign,jumps,
		    strata,nstrata,weights,offset,ZX,propodds,PACKAGE="mets"))
	 else val <- with(dd,
		   .Call("FastCoxPLstrataAddGam",pp,X,XX,sign,jumps,
		    strata,nstrata,weights,offset,ZX,
		    AddGam$theta,AddGam$dimthetades,AddGam$thetades,AddGam$ags,AddGam$varlink,AddGam$dimjumprv,AddGam$jumprv,AddGam$JumpsCauses,PACKAGE="mets"))
	 

	  if (all) {
	      val$time <- dd$time
	      val$ord <- dd$ord+1
	      val$jumps <- dd$jumps+1
	      val$jumptimes <- val$time[val$jumps]
	      val$strata.jumps <- val$strata[val$jumps]
	      val$nevent <- length(val$S0)
	      val$nstrata <- dd$nstrata
	      val$strata <- dd$strata
	      return(val)
	  }
	  with(val, structure(-ploglik,gradient=-gradient,hessian=-hessian))
	}# }}}

  opt <- NULL
  if (p>0) {
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
###      val[c("ploglik","gradient","hessian","U")] <- NULL
  }

  se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL
  II <- NULL
  ### computes Breslow estimator 
  if (cumhaz==TRUE) { # {{{
	 if (no.opt==FALSE & p!=0) II <- -solve(val$hessian) else II <- matrix(0,p,p)
	 strata <- val$strata[val$jumps]
	 nstrata <- val$nstrata
	 jumptimes <- val$jumptimes

	 ## Brewslow estimator
	 cumhaz <- cbind(jumptimes,cumsumstrata(1/val$S0,strata,nstrata))
	 if (no.opt==FALSE & p!=0) { 
	     DLambeta.t <- apply(val$E/c(val$S0),2,cumsumstrata,strata,nstrata)
	     varbetat <-   rowSums((DLambeta.t %*% II)*DLambeta.t)
	 ### covv <-  apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
	 } else varbetat <- 0
	 var.cumhaz <- cumsumstrata(1/val$S0^2,strata,nstrata)+varbetat
	 se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)

	 colnames(cumhaz)    <- c("time","cumhaz")
	 colnames(se.cumhaz) <- c("time","se.cumhaz")
 } # }}} 
 else {cumhaz <- se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL}

  res <- c(val,
           list(cox.prep=dd,
		strata.call=strata.call, strata.level=strata.level,
                entry=entry,
                exit=exit,
                status=status,                
                p=p,
                X=X,
###             id.orig=id.orig, 
                id=id.orig, 
		opt=opt, 
		cumhaz=cumhaz, se.cumhaz=se.cumhaz,
		lcumhaz=lcumhaz, lse.cumhaz=lse.cumhaz,
		II=II,strata.name=strata.name,propodds=propodds))
  class(res) <- "phreg"
  res
}

###}}} phreg0

#####' @export
###iid.phreg <- function(x) {# {{{
###  if (!is.null(seed))
###      set.seed(seed)
###  m <- lvm()
###  regression(m,T~X1+X2) <- beta
###  distribution(m,~T+C) <- coxWeibull.lvm(scale=1/100)
###  distribution(m,~entry) <- coxWeibull.lvm(scale=1/10)
###  m <- eventTime(m,time~min(T,C=0),"status")
###  d <- sim(m,n);
###  if (!entry) d$entry <- 0
###  else d <- subset(d, time>entry,select=-c(T,C))
###  return(d)
###} # }}}


###{{{ simcox

##' @export
simCox <- function(n=1000, seed=1, beta=c(1,1), entry=TRUE) {
  if (!is.null(seed))
      set.seed(seed)
  m <- lvm()
  regression(m,T~X1+X2) <- beta
  distribution(m,~T+C) <- coxWeibull.lvm(scale=1/100)
  distribution(m,~entry) <- coxWeibull.lvm(scale=1/10)
  m <- eventTime(m,time~min(T,C=0),"status")
  d <- sim(m,n);
  if (!entry) d$entry <- 0
  else d <- subset(d, time>entry,select=-c(T,C))
  return(d)
}

###}}} simcox


###{{{ phreg

##' Fast Cox PH regression
##'
##' Fast Cox PH regression
##' Robust variance is default variance with the summary. 
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param offset offsets for cox model
##' @param weights weights for Cox score equations
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases phreg phreg.par robust.phreg KM 
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' ## tracesim <- timereg::sim.cox(out1,1000)
##' ## sout1 <- phreg(Surv(time,status==1)~vf+chf+strata(wmicat.4),data=tracesim)
##' ## robust standard errors default 
##' summary(out1)
##' 
##' par(mfrow=c(1,2))
##' bplot(out1)
##' ## bplot(sout1,se=TRUE)
##' 
##' ## computing robust variance for baseline
##' rob1 <- robust.phreg(out1)
##' bplot(rob1,se=TRUE,robust=TRUE)
##' 
##' ## making iid decomposition of regression parameters
##' betaiiid <- iid(out1)
##' 
##' @export
phreg <- function(formula,data,offset=NULL,weights=NULL,...) {
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster","offset")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
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
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else strata.name <- NULL
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms  <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  res <- c(phreg01(X,entry,exit,status,id,strata,offset,weights,strata.name,...),list(call=cl,model.frame=m,formula=formula))
  class(res) <- "phreg"
  
  res
}
###}}} phreg


###{{{ vcov

##' @export
vcov.phreg  <- function(object,...) {    
  res <- crossprod(ii <- iid(object,...))
  attributes(res)$ncluster <- attributes(ii)$ncluster
  attributes(res)$invhess <- attributes(ii)$invhess
  colnames(res) <- rownames(res) <- names(coef(object))
  res
}

###}}} vcov

###{{{ coef

##' @export
coef.phreg  <- function(object,...) {
  object$coef
}

###}}} coef

###{{{ iid & Robust variances 

##' @export
iid.phreg  <- function(x,type="robust",all=FALSE,...) {# {{{
  invhess <- -solve(x$hessian)
  if (type=="robust") {	
	  xx <- x$cox.prep
	  ii <- invhess 
	  S0i <- rep(0,length(xx$strata))
	  S0i[xx$jumps+1] <- 1/x$S0
	  Z <- xx$X
	  U <- E <- matrix(0,nrow(xx$X),x$p)
	  E[xx$jumps+1,] <- x$E
	  U[xx$jumps+1,] <- x$U
	  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
	  EdLam0 <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
	  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
	  ### Martingale  as a function of time and for all subjects to handle strata 
	  MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-EdLam0)*rr*c(xx$weights)
	  orig.order <- (1:nrow(xx$X))[xx$ord+1]
	  ooo <- order(orig.order)
	  ### back to order of data-set
	  MGt <- MGt[ooo,,drop=FALSE]
	  id <- xx$id[ooo]
  } else  { 
     MGt <- x$U; MG.base <- 1/x$S0; 
  }

  ncluster <- NULL
  if (type=="robust" & (!is.null(x$id) | any(x$entry>0))) {
    if (type=="martingale") id <- x$id[x$jumps]
    ###  ii <- mets::cluster.index(id)
    UU <- apply(MGt,2,sumstrata,id,max(id)+1)
###    for (i in seq(ii$uniqueclust)) {
###      UU[i,] <- colSums(MGt[ii$idclustmat[i,seq(ii$cluster.size[i])]+1,,drop=FALSE])
###    }
    ncluster <- nrow(UU)
  } else {
     UU <- MGt
  }
  
  structure(UU%*%invhess,invhess=invhess,ncluster=ncluster)
} # }}}

##' @export
robust.basehaz.phreg  <- function(x,type="robust",fixbeta=NULL,...) {# {{{

  ### sets fixbeta based on  wheter xr has been optimized in beta (so cox case)
  if (is.null(fixbeta)) 
  if (is.null(x$opt) | is.null(x$coef)) fixbeta<- 1 else fixbeta <- 0

  if (fixbeta==0) 
  invhess <- -solve(x$hessian)
  xx <- x$cox.prep
  S0i2 <- S0i <- rep(0,length(xx$strata))
  S0i[xx$jumps+1] <-  1/x$S0
  S0i2[xx$jumps+1] <- 1/x$S0^2
  if (fixbeta==0) {
	  Z <- xx$X
	  U <- E <- matrix(0,nrow(xx$X),x$p)
	  E[xx$jumps+1,] <- x$E
	  U[xx$jumps+1,] <- x$U
	  Ht <- apply(E*S0i,2,cumsumstrata,xx$strata,xx$nstrata)
  }
  ###    
  cumhaz <- cbind(xx$time,cumsumstrata(S0i,xx$strata,xx$nstrata))
  cumS0i2 <-    cumsumstrata(S0i2,xx$strata,xx$nstrata)
  if (fixbeta==0) 
  rr <- c(xx$sign*exp(Z %*% coef(x) + xx$offset))
  else rr <- c(xx$sign*exp( xx$offset))
  id <-   xx$id
  mid <- max(id)+1
  ### also weights 
  w <- c(xx$weights)
  xxx <- w*(S0i-rr*c(cumS0i2))

  ssf <- cumsumidstratasum(xxx,id,mid,xx$strata,xx$nstrata)$sumsquare
  ss <- c(revcumsumidstratasum(w*rr,id,mid,xx$strata,xx$nstrata)$lagsumsquare)*c(cumS0i2^2)
  covv <- covfridstrata(xxx,w*rr,id,mid,xx$strata,xx$nstrata)$covs*c(cumS0i2)
  varA <- c(ssf+ss+2*covv)

  if (fixbeta==0) {# {{{
      MGt <- U[,drop=FALSE]-(Z*cumhaz[,2]-Ht)*rr*c(xx$weights)
      UU <- apply(MGt,2,sumstrata,id,max(id)+1)
      betaiid <- UU %*% invhess
     vbeta <- crossprod(betaiid)
     varbetat <-   rowSums((Ht %*% vbeta)*Ht)
     ### writing each beta for all individuals 
     betakt <- betaiid[id+1,,drop=FALSE]
     ###
     covk1 <- apply(xxx*betakt,2,cumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="sum")
     covk2 <- apply(w*rr*betakt,2,revcumsumidstratasum,id,mid,xx$strata,xx$nstrata,type="lagsum")
     covk2 <- c(covk2)*c(cumS0i2)
     ###
     varA <- varA+varbetat-2*apply((covk1-covk2)*Ht,1,sum)
  }# }}}
  varA <- varA[x$jumps]

  strata <- xx$strata[x$jumps]
  cumhaz <- x$cumhaz
  se.cumhaz <- cbind(cumhaz[,1],varA^.5)
  colnames(se.cumhaz) <- c("time","se.cumhaz")
  
  return(list(cumhaz=cumhaz,se.cumhaz=se.cumhaz,strata=strata))
} # }}}


##' @export robust.phreg
robust.phreg  <- function(x,fixbeta=NULL,...) {

  if (is.null(fixbeta)) 
  if (is.null(x$opt) | is.null(x$coef)) fixbeta<- 1 else fixbeta <- 0

 if (fixbeta==0)  {
    gamma.iid <- iid.phreg(x) 
    robvar <- crossprod(gamma.iid)
 } else robvar <- gamma.iid <- NULL
 baseline <- robust.basehaz.phreg(x,fixbeta=fixbeta,...); 
 ## add arguments so that we can call basehazplot.phreg
 return(c(x,list(gamma.iid=gamma.iid,robvar=robvar,
		 robse.cumhaz=baseline$se.cumhaz)))
}

###}}}

###{{{ summary

##' @export
summary.phreg <- function(object,type=c("robust","martingale"),...) {
  cc <- ncluster <- V <- NULL
  if (length(object$p)>0 & object$p>0 & !is.null(object$opt)) {
    I <- -solve(object$hessian)
    V <- vcov(object,type=type[1])
    cc <- cbind(coef(object),diag(V)^0.5,diag(I)^0.5)
    cc  <- cbind(cc,2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
    colnames(cc) <- c("Estimate","S.E.","dU^-1/2","P-value")
    if (!is.null(ncluster <- attributes(V)$ncluster))
    rownames(cc) <- names(coef(object))
  } 
  Strata <- levels(object$strata)
  if (!is.null(Strata)) {
    n <- unlist(lapply(object$time,length))
  } else {
    n <- length(object$time)    
  }  
  res <- list(coef=cc,n=n,nevent=object$nevent,strata=Strata,ncluster=ncluster,var=V)
  class(res) <- "summary.phreg"
  res
}

###}}} summary

###{{{ print.summary

##' @export
print.summary.phreg  <- function(x,max.strata=5,...) {
  cat("\n")
  nn <- cbind(x$n, x$nevent)
  rownames(nn) <- levels(x$strata); colnames(nn) <- c("n","events")
  if (is.null(rownames(nn))) rownames(nn) <- rep("",NROW(nn))
  if (length(x$strata)>max.strata) {
      nn <- rbind(c(colSums(nn),length(x$strata)));
      colnames(nn) <- c("n","events","stratas")
      rownames(nn) <- ""
  } 
  print(nn,quote=FALSE)  
  if (!is.null(x$ncluster)) cat("\n ", x$ncluster, " clusters\n",sep="")
  if (!is.null(x$coef)) {
    cat("\n")
    printCoefmat(x$coef,...)
  }
  cat("\n")
}

###}}} print.summary

##' @export
sumstrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("sumstrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}


##' @export
cumsumstrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("cumsumstrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}

##' @export
revcumsumstrata <- function(x,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("revcumsumstrataR",x,strata,nstrata,PACKAGE="mets")$res
return(res)
}# }}}

##' @export
revcumsum <- function(x)
{# {{{
res <- .Call("revcumsumR",x,PACKAGE="mets")$res
return(res)
}# }}}


##' @export
revcumsumstratasum <- function(x,strata,nstrata,type="all")
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")    res <- .Call("revcumsumstratasumR",x,strata,nstrata)$sum
if (type=="lagsum") res <- .Call("revcumsumstratasumR",x,strata,nstrata)$lagsum
if (type=="all")    res <- .Call("revcumsumstratasumR",x,strata,nstrata)
return(res)
}# }}}

##' @export
cumsumstratasum <- function(x,strata,nstrata,type="all")
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")    res <- .Call("cumsumstratasumR",x,strata,nstrata)$sum
if (type=="lagsum") res <- .Call("cumsumstratasumR",x,strata,nstrata)$lagsum
if (type=="all")    res <- .Call("cumsumstratasumR",x,strata,nstrata)
return(res)
}# }}}

##' @export
matdoubleindex <- function(x,rows,cols)
{# {{{
ncols <- ncol(x)
nrows <- nrow(x)
if (any(rows>nrows) | any(cols>ncols)) stop("indeces out of matrix \n"); 
if (length(cols)!=length(rows)) stop("rows and cols different lengths\n"); 
res <- .Call("Matdoubleindex",x,rows-1,cols-1,length(cols))$mat
return(res)
}# }}}

##' @export
covfr  <- function(x,y,strata,nstrata)
{# {{{
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("covrfR",x,y,strata,nstrata)
return(res)
}# }}}

##' @export
revcumsumidstratasum <- function(x,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")    res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)$sum
if (type=="lagsum")    res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)$lagsum
if (type=="lagsumsquare") res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)$lagsumsquare
if (type=="all")    res <- .Call("revcumsumidstratasumR",x,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
revcumsumidstratasumCov <- function(x,y,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")    res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$sum
if (type=="lagsum")    res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$lagsum
if (type=="lagsumsquare") res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$lagsumsquare
if (type=="all")    res <- .Call("revcumsumidstratasumCovR",x,y,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
cumsumidstratasumCov <- function(x,y,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")   res <- .Call("cumsumidstratasumCovR",x,y,id,nid,strata,nstrata)$sum
else res <- .Call("cumsumidstratasumCovR",x,y,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
cumsumidstratasum <- function(x,id,nid,strata,nstrata,type="all")
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
if (type=="sum")   res <- .Call("cumsumidstratasumR",x,id,nid,strata,nstrata)$sum
else res <- .Call("cumsumidstratasumR",x,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
covfridstrata  <- function(x,y,id,nid,strata,nstrata)
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("covrfstrataR",x,y,id,nid,strata,nstrata)
return(res)
}# }}}

##' @export
covfridstrataCov  <- function(x,y,x1,y1,id,nid,strata,nstrata)
{# {{{
if (any(id<0) | any(id>nid-1)) stop("id index not ok\n"); 
if (any(strata<0) | any(strata>nstrata-1)) stop("strata index not ok\n"); 
res <- .Call("covrfstrataCovR",x,y,x1,y1,id,nid,strata,nstrata)
return(res)
}# }}}


##' Kaplan-Meier with robust standard errors 
##'
##' Kaplan-Meier with robust standard errors 
##' Robust variance is default variance with the summary. 
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param conf.type transformation 
##' @param conf.int level of confidence intervals 
##' @param robust for robust standard errors based on martingales 
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @aliases km 
##' @examples
##' data(TRACE)
##' TRACE$cluster <- sample(1:100,1878,replace=TRUE)
##' out1 <- km(Surv(time,status==9)~strata(vf,chf),data=TRACE)
##' out2 <- km(Surv(time,status==9)~strata(vf,chf)+cluster(cluster),data=TRACE)
##' 
##' par(mfrow=c(1,2))
##' bplot(out1,se=TRUE)
##' bplot(out2,se=TRUE)
##' @export
km <- function(formula,data=data,conf.type="log",conf.int=0.95,robust=TRUE)
{# {{{
 coxo <- phreg(formula,data=data)
 coxo <- robust.phreg(coxo)

 chaz <-     coxo$cumhaz[,2]
 time <-     coxo$cumhaz[,1]
 if (robust) std.err <-  coxo$robse.cumhaz[,2]
 else std.err <-  coxo$se.cumhaz[,2]
 strat <-    coxo$strata[coxo$jumps]

 S0i  <-  1/coxo$S0
 kmt <- exp(cumsumstrata(log(1-S0i),strat,coxo$nstrata))
 temp <- list(surv=kmt)

 zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)

 if (conf.type == "plain") {
    temp1 <- temp$surv + zval * std.err * temp$surv
    temp2 <- temp$surv - zval * std.err * temp$surv
    temp <- c(temp, list(upper = pmin(temp1, 1), lower = pmax(temp2,
	0), conf.type = "plain", conf.int = conf.int))
 }
 if (conf.type == "log") {
    xx <- ifelse(temp$surv == 0, 1, temp$surv)
    temp1 <- ifelse(temp$surv == 0, NA, exp(log(xx) + zval * std.err))
    temp2 <- ifelse(temp$surv == 0, NA, exp(log(xx) - zval * std.err))
    temp <- c(temp, list(upper = pmin(temp1, 1), lower = temp2,
	conf.type = "log", conf.int = conf.int))
 }
 if (conf.type == "log-log") {
    who <- (temp$surv == 0 | temp$surv == 1)
    temp3 <- ifelse(temp$surv == 0, NA, 1)
    xx <- ifelse(who, 0.1, temp$surv)
    temp1 <- exp(-exp(log(-log(xx)) + zval * std.err/log(xx)))
    temp1 <- ifelse(who, temp3, temp1)
    temp2 <- exp(-exp(log(-log(xx)) - zval * std.err/log(xx)))
    temp2 <- ifelse(who, temp3, temp2)
    temp <- c(temp, list(upper = temp1, lower = temp2,
	conf.type = "log-log", conf.int = conf.int))
 }

 ### to use basehazplot.phreg
 temp <- c(temp,
   list(cumhaz=cbind(time,kmt),se.cumhaz=cbind(time,kmt*std.err),
	time=time,
	strata=strat,nstrata=coxo$nstrata,
        jumps=1:length(kmt), strata.name=coxo$strata.name, strata.level=coxo$strata.level))
 class(temp) <- c("km","phreg")
 return(temp)
}# }}}

##' Cumulative incidence with robust standard errors 
##'
##' Cumulative incidence with robust standard errors 
##' Robust variance is default variance with the summary. 
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param cause NULL looks at all, otherwise specify which cause to consider
##' @param cens.code censoring code "0" is default
##' @param ... Additional arguments to lower level funtions
##' @author Thomas Scheike
##' @aliases cif  
##' @examples
##' data(TRACE)
##' TRACE$cluster <- sample(1:100,1878,replace=TRUE)
##' out1 <- cif(Event(time,status)~+1,data=TRACE,cause=9)
##' out2 <- cif(Event(time,status)~+1+cluster(cluster),data=TRACE,cause=9)
##' 
##' out1 <- cif(Event(time,status)~strata(vf,chf),data=TRACE,cause=9)
##' out2 <- cif(Event(time,status)~strata(vf,chf)+cluster(cluster),data=TRACE,cause=9)
##' 
##' par(mfrow=c(1,2))
##' bplot(out1,se=TRUE)
##' bplot(out2,se=TRUE)
##' @export
cif <- function(formula,data=data,cause=1,cens.code=0,...)
{# {{{

  cl <- match.call()
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
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else strata.name <- NULL
  if (!is.null(offsetpos <- attributes(Terms)$specials$offset)) {
    ts <- survival::untangle.specials(Terms, "offset")
    Terms  <- Terms[-ts$terms]
    offset <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  id.orig <- id; 
  if (!is.null(id)) {
	  ids <- sort(unique(id))
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(exit))-1; 


  statusE <- 1*(status==cause)
  statusD <- 1*(status!=cens.code)
  if (ncol(Y)==3) {
	  if (!is.null(strata)) {
  formE <- as.formula(paste("Surv(entry=entry,exit,statusE)~strata(strata)+cluster(id_æø_)",sep=""))
  formD <- as.formula(paste("Surv(entry=entry,exit,statusD)~strata(strata)+cluster(id_æø_)",sep=""))
	  } else {
  formE <- as.formula(paste("Surv(entry=entry,exit,statusE)~1+cluster(id_æø_)",sep=""))
  formD <- as.formula(paste("Surv(entry=entry,exit,statusD)~1+cluster(id_æø_)",sep=""))

	  }
  } else {
	  if (!is.null(strata)) {
  formE <- as.formula(paste("Surv(exit,statusE)~strata(strata)+cluster(id_æø_)",sep=""))
  formD <- as.formula(paste("Surv(exit,statusD)~strata(strata)+cluster(id_æø_)",sep=""))
	  } else {
  formE <- as.formula(paste("Surv(exit,statusE)~cluster(id_æø_)",sep=""))
  formD <- as.formula(paste("Surv(exit,statusD)~cluster(id_æø_)",sep=""))
	  } 
  }

  data$id_æø_ <- id

  if (sum(statusE)==0) warning("No events of type 1\n"); 

  coxE <- phreg(formE,data=data,...)
  coxS <- phreg(formD,data=data,...)

  ### cif 
  cifo <- recurrentMarginal(coxE,coxS)

  ### to use basehazplot.phreg
  class(cifo) <- c("cif","phreg")
  return(cifo)
}# }}}


###{{{ predict with se for baseline

predictPhreg <- function(x,jumptimes,S0,beta,time=NULL,X=NULL,surv=FALSE,band=FALSE,...) {
    strata <- x$strata[x$jumps]
    nstrata <- x$nstrata
    
    ## Brewslow estimator
    if (is.null(x$cumhaz)) {
        ##II <- x$II
        ##x$jumptimes
        II <- -solve(x$hessian)
        chaz <- cbind(jumptimes,cumsumstrata(1/S0,strata,nstrata))
        DLambeta.t <- apply(x$E/c(x$S0),2,cumsumstrata,strata,nstrata)
        varbetat <- apply((DLambeta.t %*%  II)*DLambeta.t,1,sum)
        se.chaz <- cbind(jumptimes,(cumsumstrata(1/S0^2,strata,nstrata)+varbetat)^.5)
    } else {
        chaz <- x$cumhaz
        se.chaz <- x$se.cumhaz
    }

    if (!is.null(time)) {
	### do within strata
        chaz <- Cpred(chaz,time)
        se.chaz <- Cpred(se.chaz,time)
    }
    colnames(chaz) <- c("time","chaz")
    colnames(se.chaz) <- c("time","se.chaz")

    if (band==TRUE) { ## on log-scale  for one strata# {{{
      ii <- -solve(x$hessian)
      Ubeta <- x$U
      betaiid <- t(ii %*% t(Ubeta))
      cumhaz <-  x$cumhaz[,1,drop=FALSE]
      se.chaz <- x$se.cumhaz[,1]
      ###
      rr <- c( exp(sum(c(X) %*% x$coef)))
      Pt      <- outer(cumhaz[,2],c(X))
      DLambeta.t <- apply(x$E/c(x$S0),2,cumsumstrata,strata,nstrata)
      Pt <- DLambeta.t  - Pt
      ### se of cumulaive hazard for this covariate , can use different versions of variance for beta
      varbetat <- rowSums((Pt %*% ii)*Pt)
      se.chazexb <- cbind(jumptimes,rr*(cumsumstrata(1/S0^2,strata,nstrata)+varbetat)^.5)
###      sig <- 0.95
###      n.sim <- 1000
###      simband <-  .Call("simBandCumHazCox",rr/x$S0,Pt,betaiid,n.sim,sig,se.chazexb[,2],PACKAGE="mets")
###      ## coefficients for uniform bands on log-scale 
###      uband <- apply(simband$supUsim,2,percen,per=1-sig);
    }# }}}

    if (!is.null(X)) {
      H <- exp(X%*%beta)
      if (nrow(chaz)==length(H)) {
        chaz[,2] <- chaz[,2]*H
      } else {
        chaz2 <- c()
        X <- rbind(X)
        for (i in seq(nrow(X)))
          chaz2 <- rbind(chaz2,
                         cbind(chaz[,1],chaz[,2]*H[i],
                               rep(1,nrow(chaz))%x%X[i,,drop=FALSE]))
        chaz <- chaz2;
        nn <- c("time","chaz",names(beta))
        colnames(chaz) <- nn
      }
    }
    if (surv) {    
      chaz[,2] <- exp(-chaz[,2])
      colnames(chaz)[2] <- "surv"
    }
    return(chaz)
}

##' Predictions from proportional hazards model
##'
##' @export
##' @param object phreg object
##' @param data data.frame
##' @param surv If TRUE predictions are provided on probability scale
##' @param time Time variable
##' @param X Design matrix
##' @param strata Strata variable
##' @param ... ADditional arguments to lower level functions
##' @aliases predict.phreg revcumsumstrata revcumsumstratasum cumsumstrata sumstrata covfr covfridstrata covfridstrataCov cumsumidstratasum cumsumidstratasumCov cumsumstratasum revcumsumidstratasum revcumsumidstratasumCov robust.basehaz.phreg matdoubleindex
predict.phreg  <- function(object,data,surv=FALSE,time=object$exit,X=object$X,strata=object$strata,...) {
    if (object$p==0) X <- NULL
    if (!is.null(object$strata)) {
        lev <- levels(object$strata)
        if (!is.null(object$strata) &&
            !(is.list(time) & !is.data.frame(time)) &&
            !(is.list(X) & !is.data.frame(X))) {
            X0 <- X
            time0 <- time
            X <- time <- c()
            for (i in seq(length(lev))) {
                idx <- which(strata==lev[i])
                X <- c(X,list(X0[idx,,drop=FALSE]))
                time <- c(time,list(time0[idx]))
            }
        }
        chaz <- c()
        for (i in seq(length(lev)))
            chaz <- c(chaz,list(predictPhreg(object$jumptimes[[i]],
                                             object$S0[[i]],
                                             coef(object),
                                             time[[i]],X[[i]],surv)))
        names(chaz) <- lev    
    } else {
        chaz <- predictPhreg(object$jumptimes,object$S0,coef(object),time,X,surv)
    }
    return(chaz)
}

###}}} predict

###{{{ plot

##' Plotting the baslines of stratified Cox 
##'
##' Plotting the baslines of stratified Cox 
##' @param x phreg object
##' @param se to include standard errors
##' @param time to plot for specific time variables
##' @param add to add to previous plot 
##' @param ylim to give ylim 
##' @param xlim to give xlim 
##' @param lty to specify lty of components
##' @param col to specify col of components
##' @param legend to specify col of components
##' @param ylab to specify ylab 
##' @param polygon to get standard error in shaded form
##' @param level of standard errors
##' @param stratas wich strata to plot 
##' @param robust to use robust standard errors if possible
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @export
##' @aliases basehazplot.phreg  bplot  basecumhaz
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' 
##' par(mfrow=c(2,2))
##' bplot(out1)
##' bplot(out1,stratas=c(0,3))
##' bplot(out1,stratas=c(0,3),col=2:3,lty=1:2,se=TRUE)
##' bplot(out1,stratas=c(0),col=2,lty=2,se=TRUE,polygon=FALSE)
##' bplot(out1,stratas=c(0),col=matrix(c(2,1,3),1,3),
##'             lty=matrix(c(1,2,3),1,3),se=TRUE,polygon=FALSE)
##' @export
basehazplot.phreg  <- function(x,se=FALSE,time=NULL,add=FALSE,ylim=NULL,xlim=NULL,
    lty=NULL,col=NULL,legend=TRUE,ylab=NULL,
    polygon=TRUE,level=0.95,stratas=NULL,robust=FALSE,...) {# {{{
	if (class(x)[1]=="phreg" & is.null(ylab)) ylab <- "Cumulative hazard"
	if (class(x)[1]=="km" & is.null(ylab)) ylab <- "Survival probability"
   level <- -qnorm((1-level)/2)
###   if (log==FALSE) 
   rr <- range(x$cumhaz[,-1]) ### else rr <- range(log(x$cumhaz[,-1]),na.rm=TRUE)
   strat <- x$strata[x$jumps]
   ylimo <- ylim
   if (is.null(ylim)) ylim <- rr
   if (is.null(xlim)) xlim <- range(x$cumhaz[,1])
   if (se==TRUE) {
	   if (is.null(x$se.cumhaz) & is.null(x$robse.cumhaz) ) 
		   stop("phreg must be with cumhazard=TRUE\n"); 
       rrse <- range(c(x$cumhaz[,-1]+level*x$se.cumhaz[,-1])) 
       if (class(x)[1]=="km") rrse <- c(min(x$lower),1)
###       else rrse <- range(log(c(x$cumhaz[,-1]+level*x$se.cumhaz[,-1])),na.rm=TRUE)
       if (is.null(ylimo)) ylim <- rrse
   }

   ## all strata
   if (is.null(stratas)) stratas <- 0:(x$nstrata-1) 

   ltys <- lty
   cols <- col

   if (length(stratas)>0 & x$nstrata>1) { ## with strata
   lstrata <- x$strata.level[(stratas+1)]
   stratn <-  substring(x$strata.name,8,nchar(x$strata.name)-1)
   stratnames <- paste(stratn,lstrata,sep=":")
      if (!is.matrix(lty)) {
         if (is.null(lty)) ltys <- 1:length(stratas) else if (length(lty)!=length(stratas)) ltys <- rep(lty[1],length(stratas))
      } else ltys <- lty
      if (!is.matrix(col)) {
         if (is.null(col)) cols <- 1:length(stratas) else 
		 if (length(col)!=length(stratas)) cols <- rep(col[1],length(stratas))
      } else cols <- col
   } else { 
     stratnames <- "Baseline" 
     if (is.matrix(col))  cols <- col
     if (is.null(col)) cols <- 1  else cols <- col[1]
     if (is.matrix(lty))  ltys <- lty
     if (is.null(lty)) ltys <- 1  else ltys <- lty[1]
   }

  if (!is.matrix(ltys))  ltys <- cbind(ltys,ltys,ltys)
  if (!is.matrix(cols))  cols <- cbind(cols,cols,cols)

   i <- 1
   j <- stratas[i]
        cumhazard <- x$cumhaz[strat==j,,drop=FALSE]

   if (!is.null(cumhazard)) {
   if (nrow(cumhazard)>1) {
    if (add) {
         lines(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1],...)
    } else {
         plot(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1],ylim=ylim,ylab=ylab,xlim=xlim,...)
    }
    if (se==TRUE) {
	    if (robust==TRUE) secumhazard  <- x$robse.cumhaz[strat==j,,drop=FALSE]
	    else secumhazard <- x$se.cumhaz[strat==j,,drop=FALSE]
      ul <-cbind(cumhazard[,1],cumhazard[,2]+level*secumhazard[,2])
      nl <-cbind(cumhazard[,1],cumhazard[,2]-level*secumhazard[,2])
      if (class(x)[1]=="km") { ul[,2] <- x$upper[x$strata==j]; 
                               nl[,2] <- x$lower[x$strata==j];
      }
      if (!polygon) {
      lines(nl,type="s",lty=ltys[i,2],col=cols[i,2])
      lines(ul,type="s",lty=ltys[i,3],col=cols[i,3])
      } else {
         tt <- c(nl[,1],rev(ul[,1]))
         yy <- c(nl[,2],rev(ul[,2]))
         col.alpha<-0.1
         col.ci<-cols[j+1]
         col.trans <- sapply(col.ci, FUN=function(x) 
                   do.call(grDevices::rgb,as.list(c(grDevices::col2rgb(x)/255,col.alpha))))
	 polygon(tt,yy,lty=ltys[i,2],col=col.trans)
      }
    }
   }
   }

    if (length(stratas)>1)  {
    for (i in 2:length(stratas)) {
	j <- stratas[i]
        cumhazard <- x$cumhaz[strat==j,,drop=FALSE]
        if (!is.null(cumhazard)) {
	if (nrow(cumhazard)>1) {
        lines(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1])   
        if (se==TRUE) {
	    if (robust==TRUE) secumhazard  <- x$robse.cumhaz[strat==j,,drop=FALSE]
	    else secumhazard <- x$se.cumhaz[strat==j,,drop=FALSE]
		 ul <-cbind(cumhazard[,1],cumhazard[,2]+level*secumhazard[,2])
		 nl <-cbind(cumhazard[,1],cumhazard[,2]-level*secumhazard[,2])
		      if (class(x)[1]=="km") { ul[,2] <- x$upper[x$strata==j]; 
					       nl[,2] <- x$lower[x$strata==j];
		      }
	      if (!polygon) {
	      lines(nl,type="s",lty=ltys[i,2],col=cols[i,2])
	      lines(ul,type="s",lty=ltys[i,3],col=cols[i,3])
	      } else {
		 tt <- c(nl[,1],rev(ul[,1]))
		 yy <- c(nl[,2],rev(ul[,2]))
		 col.alpha<-0.1
		 col.ci<-cols[j+1]
		 col.trans <- sapply(col.ci, FUN=function(x) 
			   do.call(grDevices::rgb,as.list(c(grDevices::col2rgb(x)/255,col.alpha))))
		 polygon(tt,yy,lty=ltys[i,2],col=col.trans)
	      }
        }
        }
        }
    }
    }

    where <- "topleft"; 
    if (class(x)[1]=="km") where <-  "topright"
    if (legend & (!add)) 
    graphics::legend(where,legend=stratnames,col=cols[,1],lty=ltys[,1])

}# }}}

##' @export
bplot <- function(x,...) basehazplot.phreg(x,...)

##' @export
basecumhaz <- function(x,type="matrix",robust=FALSE,...) {# {{{
   ## all strata
   strat <- x$strata[x$jumps]
   stratas <- 0:(x$nstrata-1) 

###   se.cum <- cum <- x$cumhaz
   se.cum <- cum <- c()
   strata <- rep(0,nrow(x$cumhaz))
   if (type=="matrix") { se.cum <- cum <- x$cumhazard }
   if (robust==TRUE) secum <- x$robse.cumhaz else secum  <- x$se.cumhaz
   if (is.null(secum)) nose <- TRUE else nose <- FALSE

   start <- 1
   for (i in stratas) {
	   cumhazard <- x$cumhaz[strat==i,,drop=FALSE]
	   if (!is.null(cumhazard)) {
		   nr <- nrow(cumhazard)
		   if (nr>=1) {
		   slut <- start-1+nr
	###	   cum[start:slut,] <- cumhazard
		   cum <- rbind(cum,cumhazard)
	###	   if (!nose) se.cum[start:slut,] <- secum[strat==i,]
		   if (!nose) se.cum <- rbind(se.cum,secum[strat==i,])
		   strata[start:slut] <- i
		   start <- slut+1
	      }
	   }
   }

   list(cumhaz=cum,se.cumhaz=se.cum,strata=strata)
}# }}}


##' @export
lines.phreg <- function(x,...,add=TRUE) plot(x,...,add=add)

###}}} plot

###{{{ plot

##' @export
plot.phreg  <- function(x,surv=TRUE,X=NULL,time=NULL,add=FALSE,...) {
    if (!is.null(X) && nrow(X)>1) {
        P <- lapply(split(X,seq(nrow(X))),function(xx) predict(x,X=xx,time=time,surv=surv))
    } else {
        P <- predict(x,X=X,time=time,surv=surv)
    }
    if (!is.list(P)) {
        if (add) {
            lines(P,type="s",...)
        } else {
            plot(P,type="s",...)
        }
        return(invisible(P))
    }

    if (add) {
        lines(P[[1]][,1:2],type="s",lty=1,col=1,...)
    } else {
        plot((P[[1]])[,1:2],type="s",lty=1,col=1,...)
    }
    for (i in seq_len(length(P)-1)+1) {
        lines(P[[i]][,1:2],type="s",lty=i,col=i,...)   
    }
    return(invisible(P))
}

##' @export
lines.phreg <- function(x,...,add=TRUE) plot(x,...,add=add)

###}}} plot

###{{{ print
##' @export
print.phreg  <- function(x,...) {
  cat("Call:\n")
  dput(x$call)
  print(summary(x),...)
}
###}}} print

