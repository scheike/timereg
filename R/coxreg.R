###{{{ coxreg0 
coxreg0 <- function(X,entry,exit,status,id=NULL,strata=NULL,beta,stderr=TRUE,...) {
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(exit)))
  if (!is.null(strata)) {
    stratalev <- levels(strata)
    strataidx <- lapply(stratalev,function(x) which(strata==x))
    if (!all(unlist(lapply(strataidx,function(x) length(x)>0))))
      stop("Strata without any observation")
    dd <- lapply(strataidx, function(ii)                      
                 .Call("FastCoxPrep",
                       entry[ii],exit[ii],status[ii],
                       as.matrix(X)[ii,,drop=FALSE],
                       id[ii],
                       package="mets"))
    if (!is.null(id))
      id <- unlist(lapply(dd,function(x) x$id[x$jumps+1]))
    obj <- function(pp,U=FALSE,all=FALSE) {
      val <- lapply(dd,function(d)
                    with(d,
                         .Call("FastCoxPL",pp,X,XX,sign,jumps,package="mets")))
      ploglik <- do.call("+",lapply(val,function(x) x$ploglik))
      gradient <- do.call("+",lapply(val,function(x) x$gradient))
      hessian <- do.call("+",lapply(val,function(x) x$hessian))
      if (all) {
        U <- do.call("rbind",lapply(val,function(x) x$U))
        time <- lapply(dd,function(x) x$time[x$ord+1])
        ord <- lapply(dd,function(x) x$ord+1)
        jumps <- lapply(dd,function(x) x$jumps+1)
        jumptimes <- lapply(dd,function(x) x$time[x$ord+1][x$jumps+1])
        S0 <- lapply(val,function(x) x$S0)
        nevent  <- unlist(lapply(S0,length))
        return(list(ploglik=ploglik,gradient=gradient,hessian=hessian,
                    U=U,S0=S0,nevent=nevent,
                    ord=ord,time=time,jumps=jumps,jumptimes=jumptimes))
      }
      structure(-ploglik,gradient=-gradient,hessian=-hessian)
    }
  } else {
    system.time(dd <- .Call("FastCoxPrep",entry,exit,status,X,id,package="mets"))
    if (!is.null(id))
      id <- dd$id[dd$jumps+1]
    obj <- function(pp,U=FALSE,all=FALSE) {
      val <- with(dd,
                  .Call("FastCoxPL",pp,X,XX,sign,jumps,package="mets"))
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
  }
  if (p>0) {
    opt <- nlm(obj,beta)
    cc <- opt$estimate; names(cc) <- colnames(X)
    if (!stderr) return(cc)    
    val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
  } else {
    val <- obj(0,all=TRUE)
    val[c("ploglik","gradient","hessian","U")] <- NULL
  }
  res <- c(val,
           list(strata=strata,
                entry=entry,
                exit=exit,
                status=status,
                p=p,
                X=X,
                id=id))
  class(res) <- "coxreg"
  res
}
###}}} coxreg0

###{{{ coxreg
##' Fast Cox PH regression
##'
##' Fast Cox PH regression
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst
##' @export
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
##' 
##' n <- 1e3;
##' d <- simcox(n); d$id <- seq(nrow(d)); d$group <- factor(rbinom(nrow(d),1,0.5))
##' 
##' (m1 <- coxreg(Surv(entry,time,status)~X1+X2,data=d))
##' (m2 <- coxph(Surv(entry,time,status)~X1+X2+cluster(id),data=d))
##' (coef(m3 <-cox.aalen(Surv(entry,time,status)~prop(X1)+prop(X2),data=d)))
##' 
##' \dontrun{
##' (m1b <- coxreg(Surv(entry,time,status)~X1+X2+strata(group),data=d))
##' (m2b <- coxph(Surv(entry,time,status)~X1+X2+cluster(id)+strata(group),data=d))
##' (coef(m3b <-cox.aalen(Surv(entry,time,status)~-1+group+prop(X1)+prop(X2),data=d)))
##' }
##' 
##' m <- coxreg(Surv(entry,time,status)~X1*X2+strata(group)+cluster(id),data=d)
##' m
##' plot(m,ylim=c(0,5))
coxreg <- function(formula,data,...) {
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
    ts <- untangle.specials(Terms, "cluster")
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- untangle.specials(Terms, "strata")
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)
  res <- c(coxreg0(X,entry,exit,status,id,strata,...),list(call=cl))
  class(res) <- "coxreg"
  res
}
###}}} coxreg

###{{{ vcov
##' @S3method vcov coxreg
vcov.coxreg  <- function(object,...) {
  I <- -solve(object$hessian)
  ncluster <- NULL
  if (!is.null(object$id)) {
    ii <- mets::cluster.index(object$id)
    UU <- matrix(nrow=ii$uniqueclust,ncol=ncol(I))
    for (i in seq(ii$uniqueclust)) {
      UU[i,] <- colSums(object$U[ii$idclustmat[i,seq(ii$cluster.size[i])]+1,,drop=FALSE])
    }
    ncluster <- nrow(UU)
    J <- crossprod(UU)
  } else {
    J <- crossprod(object$U)
  }
  res <- I%*%J%*%t(I)
  attributes(res)$ncluster <- ncluster
  colnames(res) <- rownames(res) <- names(coef(object))
  res
}
###}}} vcov

###{{{ coef
##' @S3method coef coxreg
coef.coxreg  <- function(object,...) {
  object$coef
}
###}}} coef

###{{{ summary
##' @S3method summary coxreg
summary.coxreg <- function(object,se="robust",...) {
  cc <- NULL
  if (object$p>0) {
    I <- -solve(object$hessian)
    V <- vcov(object)
    cc <- cbind(coef(object),diag(V)^0.5,diag(I)^0.5)
    cc  <- cbind(cc,2*(1-pnorm(abs(cc[,1]/cc[,2]))))
    colnames(cc) <- c("Estimate","S.E.","dU^-1/2","P-value")
    if (!is.null(attributes(V)$ncluster))
    rownames(cc) <- names(coef(object))
  }
  Strata <- levels(object$strata)
  if (!is.null(Strata)) {
    n <- unlist(lapply(object$time,length))
  } else {
    n <- length(object$time)    
  }  
  res <- list(coef=cc,n=n,nevent=object$nevent,
              strata=Strata,ncluster=attributes(V)$ncluster)
  class(res) <- "summary.coxreg"
  res
}
###}}} summary

###{{{ print.summary
##' @S3method print summary.coxreg
print.summary.coxreg  <- function(x,...) {
  cat("\n")
  nn <- cbind(x$n, x$nevent)
  rownames(nn) <- x$strata; colnames(nn) <- c("n","events")
  if (is.null(rownames(nn))) rownames(nn) <- rep("",NROW(nn))
  print(nn,quote=FALSE)
  if (!is.null(x$ncluster)) cat("\n ", x$ncluster, " clusters\n",sep="")
  if (!is.null(x$coef)) {
    cat("\n")
    printCoefmat(x$coef,...)
  }
  cat("\n")
}
###}}} print.summary

###{{{ predict
##' @S3method predict coxreg
predict.coxreg  <- function(object,surv=FALSE,...) {
  if (!is.null(object$strata)) {
    lev <- levels(object$strata)
    chaz <- c()
    for (i in seq(length(lev))) {
      ## Brewslow estimator
      chaz0 <- cbind(object$jumptimes[[i]],cumsum(1/object$S0[[i]]))
      colnames(chaz0) <- c("time","chaz")
      chaz <- c(chaz,list(chaz0))
    }
    names(chaz) <- lev
  } else {
    chaz <- cbind(object$jumptimes,cumsum(1/object$S0))
    colnames(chaz) <- c("time","chaz")
  }
  return(chaz)
}
###}}} predict

###{{{ plot
##' @S3method plot coxreg
plot.coxreg  <- function(x,surv=FALSE,add=FALSE,...) {
  P <- predict(x)
  if (!is.list(P)) {
    if (add) {
      lines(P,type="s",...)
    } else {
      plot(P,type="s",...)
    }
    return(invisible(P))
  }
  plot(P[[1]],type="s",lty=1,col=1,...)
  for (i in seq_len(length(P)-1)+1) {
    lines(P[[i]],type="s",lty=i,col=i,...)   
  }  
}
###}}} plot

###{{{ print
##' @S3method print coxreg
print.coxreg  <- function(x,...) {
  cat("Call:\n")
  dput(x$call)
  print(summary(x),...)
}
###}}} print


