###{{{ phreg0 

phreg0 <- function(X,entry,exit,status,id=NULL,strata=NULL,beta,stderr=TRUE,method="NR",...) {
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
}

###}}} phreg0

###{{{ phreg01

phreg01 <- function(X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
		    strata.name=NULL,cumhaz=TRUE,
  beta,stderr=TRUE,method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,...) {
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(exit)))
  if (is.null(strata)) { strata <- rep(0,length(exit)); nstrata <- 1} else {
	  ustrata <- unique(strata)
	  nstrata <- length(ustrata)
	  strata.values <- unique(strata)
      if (is.numeric(strata)) strata <-  fast.approx(ustrata,strata)-1 else  {
      strata <- as.integer(factor(strata,labels=seq(nstrata)))-1
    }
  }
  if (is.null(offset)) offset <- rep(0,length(exit)) 
  if (is.null(weights)) weights <- rep(1,length(exit)) 

   Zcall <- matrix(1,1,1) ## to not use for ZX products when Z is not given 
   if (!is.null(Z)) Zcall <- Z

   trunc <- (!is.null(entry))
   if (!trunc) entry <- rep(0,length(exit))

   system.time(dd <- .Call("FastCoxPrepStrata",
			      entry,exit,status,X, as.integer(seq_along(entry)),
			      trunc,strata,weights,offset,Zcall,PACKAGE="mets"))
   dd$nstrata <- nstrata

	if (!is.null(id))
	  id <- dd$id[dd$jumps+1]
	obj <- function(pp,U=FALSE,all=FALSE) {
		if (is.null(propodds)) 
	  val <- with(dd,
		   .Call("FastCoxPLstrata",pp,X,XX,sign,jumps,
		    strata,nstrata,weights,offset,ZX,PACKAGE="mets"))
         else val <- with(dd,
		   .Call("FastCoxPLstrataPO",pp,X,XX,sign,jumps,
		    strata,nstrata,weights,offset,ZX,propodds,PACKAGE="mets"))
	  if (all) {
	      val$time <- dd$time[dd$ord+1]
	      val$ord <- dd$ord+1
	      val$jumps <- dd$jumps+1
	      val$jumptimes <- val$time[val$jumps]
	      val$nevent <- length(val$S0)
	      val$nstrata <- dd$nstrata
	      val$strata <- dd$strata
	      return(val)
	  }
	  with(val, structure(-ploglik,gradient=-gradient,hessian=-hessian))
	}
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

	 II <- -solve(val$hessian)
	 strata <- val$strata[val$jumps]
	 nstrata <- val$nstrata
	 jumptimes <- val$jumptimes

	 ## Brewslow estimator
	 cumhaz <- cbind(jumptimes,cumsumstrata(1/val$S0,strata,nstrata))
	 DLambeta.t <- apply(val$E/c(val$S0),2,cumsumstrata,strata,nstrata)
	 varbetat <-   rowSums((DLambeta.t %*% II)*DLambeta.t)
	 ### covv <-  apply(covv*DLambeta.t,1,sum) Covaraince is "0" by construction
	 var.cumhaz <- cumsumstrata(1/val$S0^2,strata,nstrata)+varbetat
	 se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)

	 colnames(cumhaz)    <- c("time","cumhaz")
	 colnames(se.cumhaz) <- c("time","se.cumhaz")

###    if (nstrata>1) {
###	 lcumhaz <- lse.cumhaz <- list()
### 	 cumhaz    <- cbind(cumhaz,val$strata[val$jumps])
###         se.cumhaz <- cbind(se.cumhaz,val$strata[val$jumps])
###	 for (i in 0:(nstrata-1)) {
###		 ii <- (val$strata[val$jumps]==i)
###		 if (length(ii)>1) {
###			 lcumhaz[[i+1]]  <-  cumhaz[ii,1:2]
###			 lse.cumhaz[[i+1]]  <- se.cumhaz[ii,1:2]
###		 }
###	 }
###    }
 } # }}} 
 else {cumhaz <- se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL}

  res <- c(val,
           list(strata=strata,
                entry=entry,
                exit=exit,
                status=status,                
                p=p,
                X=X,
                id=id, 
		opt=opt, 
		cumhaz=cumhaz, se.cumhaz=se.cumhaz,
		lcumhaz=lcumhaz, lse.cumhaz=lse.cumhaz,
		II=II,strata.name=strata.name,propodds=propodds))
  class(res) <- "phreg"
  res
}

###}}} phreg0

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
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst
##' @export
##' @aliases phreg phreg.par
##' @examples
##' simcox <- function(n=1000, seed=1, beta=c(1,1), entry=TRUE) {
##'   if (!is.null(seed))
##'     set.seed(seed)
##'   library(lava)
##'   m <- lvm()
##'   regression(m,T~X1+X2) <- beta
##'   distribution(m,~T+C) <- coxWeibull.lvm(scale=1/100)
##'   distribution(m,~entry) <- coxWeibull.lvm(scale=1/10)
##'   m <- eventTime(m,time~min(T,C=0),"status")
##'   d <- sim(m,n);
##'   if (!entry) d$entry <- 0
##'   else d <- subset(d, time>entry,select=-c(T,C))
##'   return(d)
##' }
##' \dontrun{
##' n <- 10;
##' d <- mets:::simCox(n); d$id <- seq(nrow(d)); d$group <- factor(rbinom(nrow(d),1,0.5))
##' 
##' (m1 <- phreg(Surv(entry,time,status)~X1+X2,data=d))
##' (m2 <- coxph(Surv(entry,time,status)~X1+X2+cluster(id),data=d))
##' (coef(m3 <-cox.aalen(Surv(entry,time,status)~prop(X1)+prop(X2),data=d)))
##' 
##' 
##' (m1b <- phreg(Surv(entry,time,status)~X1+X2+strata(group),data=d))
##' (m2b <- coxph(Surv(entry,time,status)~X1+X2+cluster(id)+strata(group),data=d))
##' (coef(m3b <-cox.aalen(Surv(entry,time,status)~-1+group+prop(X1)+prop(X2),data=d)))
##' 
##' m <- phreg(Surv(entry,time,status)~X1*X2+strata(group)+cluster(id),data=d)
##' m
##' plot(m,ylim=c(0,1))
##' }
phreg <- function(formula,data,...) {
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster")
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
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  res <- c(phreg0(X,entry,exit,status,id,strata,...),list(call=cl,model.frame=m))
  class(res) <- "phreg"
  
  res
}
###}}} phreg

###{{{ phreg1
##' Fast Cox PH regression
##'
##' Fast Cox PH regression
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param weights weights for Cox score equations
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @export
##' @aliases phreg phreg.par
##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg1(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' tracesim <- sim.cox(out1,1000)
##' sout1 <- phreg1(Surv(time,status==1)~vf+chf+strata(wmicat.4),data=tracesim)
##' 
##' par(mfrow=c(1,2))
##' baseplot.phreg(out1)
##' baseplot.phreg(sout1)
##' 
##' @export
phreg1 <- function(formula,data,offset=NULL,weights=NULL,...) {
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

  res <- c(phreg01(X,entry,exit,status,id,strata,offset,weights,strata.name,...),list(call=cl,model.frame=m))
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

###{{{ iid

##' @export
iid.phreg  <- function(x,...) {
    invhess <- solve(x$hessian)
  ncluster <- NULL
  if (!is.null(x$id)) {
    ii <- mets::cluster.index(x$id)
    UU <- matrix(nrow=ii$uniqueclust,ncol=ncol(invhess))
    for (i in seq(ii$uniqueclust)) {
      UU[i,] <- colSums(x$U[ii$idclustmat[i,seq(ii$cluster.size[i])]+1,,drop=FALSE])
    }
    ncluster <- nrow(UU)
  } else {
      UU <- x$U
  }
  structure(UU%*%invhess,invhess=invhess,ncluster=ncluster)
}

###}}}

###{{{ summary

##' @export
summary.phreg <- function(object,se="robust",...) {
  cc <- ncluster <- NULL
  if (length(object$p)>0 && object$p>0) {
    I <- -solve(object$hessian)
    V <- vcov(object)
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
  res <- list(coef=cc,n=n,nevent=object$nevent,
              strata=Strata,ncluster=ncluster)
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

cumsumstrata <- function(x,strata,nstrata)
{# {{{
res <- .Call("cumsumstrataR",x,strata,nstrata)$res
return(res)
}# }}}


###{{{ predict with se for baseline

predictPhreg <- function(x,jumptimes,S0,beta,time=NULL,X=NULL,surv=FALSE,...) {
###    x <- x$jumptimes
###    S0 <- x$S0
###    II <- x$II ###  -solve(x$hessian)
###    strata <- x$strata[x$jumps]
###    nstrata <- x$nstrata

    ## Brewslow estimator
    if (is.null(x$cumhaz)) {
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
      rr <- c( exp(sum(c(Z) %*% x$coef)))
      Pt      <- outer(cumhaz[,2],c(Z))
      DLambeta.t <- apply(x$E/c(x$S0),2,cumsumstrata,strata,nstrata)
      Pt <- DLambeta.t  - Pt
      ### se of cumulaive hazard for this covariate , can use different versions of variance for beta
      varbetat <- rowSums((Pt %*% ii)*Pt)
###   varbetat <- rowSums((Pt %*% vcov(x))*Pt)
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

##' @export
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
##' @param lty to specify lty of components
##' @param col to specify col of components
##' @param legend to specify col of components
##' @param ylab to specify ylab 
##' @param polygon to get standard error in shaded form
##' @param level of standard errors
##' @param stratas wich strata to plot 
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @export
##' @aliases phreg phreg.par

##' @examples
##' data(TRACE)
##' dcut(TRACE) <- ~.
##' out1 <- phreg1(Surv(time,status==9)~vf+chf+strata(wmicat.4),data=TRACE)
##' 
##' par(mfrow=c(2,2))
##' baseplot.phreg(out1)
##' baseplot.phreg(out1,stratas=c(0,3))
##' baseplot.phreg(out1,stratas=c(0,3),col=2:3,lty=1:2,se=TRUE)
##' baseplot.phreg(out1,stratas=c(0),col=2,lty=2,se=TRUE,polygon=FALSE)
##' baseplot.phreg(out1,stratas=c(0),col=matrix(c(2,1,3),1,3),
##'                                  lty=matrix(c(1,2,3),1,3),se=TRUE,polygon=FALSE)
##' @export
baseplot.phreg  <- function(x,se=FALSE,time=NULL,add=FALSE,ylim=NULL,
			    lty=NULL,col=NULL,legend=TRUE,
			    ylab="Cumulative hazard",polygon=TRUE,level=0.95,stratas=NULL,...) {# {{{

   level <- -qnorm((1-level)/2)
   rr <- range(x$cumhaz[,-1])
   strat <- x$strata[x$jumps]
   if (is.null(ylim)) ylim <- rr
   if (se==TRUE) {
	   if (is.null(x$se.cumhaz)) stop("phreg must be with cumhazard=TRUE\n"); 
       rrse <- range(c(x$cumhaz[,-1]+level*x$se.cumhaz[,-1]))
       ylim <- rrse
   }

   ## all strata
   if (is.null(stratas)) stratas <- 0:(x$nstrata-1) 

   ltys <- lty
   cols <- col
   if (length(stratas)>0 & x$nstrata>1) { ## with strata
      ms <- match(x$strata.name,names(x$model.frame))
      lstrata <- levels(x$model.frame[,ms])[(stratas+1)]
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
   cumhazard <- x$cumhaz[strat==j,]
    if (add) {
        lines(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1],...)
    } else {
         plot(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1],ylim=ylim,ylab=ylab,...)
    }
    if (se==TRUE) {
      secumhazard <- x$se.cumhaz[strat==j,]
      ul <-cbind(cumhazard[,1],cumhazard[,2]+level*secumhazard[,2])
      nl <-cbind(cumhazard[,1],cumhazard[,2]-level*secumhazard[,2])
      if (!polygon) {
      lines(nl,type="s",lty=ltys[i,2],col=cols[i,2])
      lines(ul,type="s",lty=ltys[i,3],col=cols[i,3])
      } else {
         tt <- c(nl[,1],rev(ul[,1]))
         yy <- c(nl[,2],rev(ul[,2]))
         col.alpha<-0.1
         col.ci<-cols[j+1]
         col.trans <- sapply(col.ci, FUN=function(x) 
                   do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
	 polygon(tt,yy,lty=ltys[i,2],col=col.trans)
      }
    }

    if (length(stratas)>1)  {
    for (i in 2:length(stratas)) {
	j <- stratas[i]
        cumhazard <- x$cumhaz[strat==j,]
        lines(cumhazard,type="s",lty=ltys[i,1],col=cols[i,1])   
        if (se==TRUE) {
         secumhazard <- x$se.cumhaz[strat==j,]
         ul <-cbind(cumhazard[,1],cumhazard[,2]+level*secumhazard[,2])
         nl <-cbind(cumhazard[,1],cumhazard[,2]-level*secumhazard[,2])
      if (!polygon) {
      lines(nl,type="s",lty=ltys[i,2],col=cols[i,2])
      lines(ul,type="s",lty=ltys[i,3],col=cols[i,3])
      } else {
         tt <- c(nl[,1],rev(ul[,1]))
         yy <- c(nl[,2],rev(ul[,2]))
         col.alpha<-0.1
         col.ci<-cols[j+1]
         col.trans <- sapply(col.ci, FUN=function(x) 
                   do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
	 polygon(tt,yy,lty=ltys[1,2],col=col.trans)
      }
    }
    }
    }

    if (legend)
    legend("topleft",legend=stratnames,col=cols[,1],lty=ltys[,1])

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

