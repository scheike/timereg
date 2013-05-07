fast.cox0 <- function(X,entry,exit,status,id=NULL,strata=NULL,beta,...) {
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
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
    obj <- function(pp,U=FALSE) {
      val <- lapply(dd,function(d)
                    with(d,
                         .Call("FastCoxPL",pp,X,XX,sign,jumps,package="mets")))
      ploglik <- do.call("+",lapply(val,function(x) x$ploglik))
      gradient <- do.call("+",lapply(val,function(x) x$gradient))
      hessian <- do.call("+",lapply(val,function(x) x$hessian))
      if (U) U <- do.call("rbind",lapply(val,function(x) x$U))
      structure(-ploglik,gradient=-gradient,hessian=-hessian,U=U)
    }
  } else {  
    dd <- .Call("FastCoxPrep",entry,exit,status,as.matrix(X),id,package="mets")
    if (!is.null(id))
      id <- dd$id[dd$jumps+1]
    obj <- function(pp,U=FALSE) {
      val <- with(dd,
                  .Call("FastCoxPL",pp,X,XX,sign,jumps,package="mets"))
      if (U) U <- val$U
      with(val, structure(-ploglik,gradient=-gradient,hessian=-hessian,U=U))
    }
  }
  opt <- nlm(obj,beta)
  val <- obj(opt$estimate,TRUE)
  cc <- opt$estimate; names(cc) <- colnames(X)

  res <- list(coef=cc,
              score=-1*attributes(val)$U,
              hessian=-1*attributes(val)$hessian,
              logLik=-1*attributes(val)$ploglik,
              strata=strata,
              id=id)
  class(res) <- "fast.cox"
  res
}

##' Fast Cox PH regression
##'
##' Fast Cox PH regression
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param ... Additional arguments to lower level funtions
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
##' d <- simcox(n); d$id <- seq(nrow(d)); d$group <- rbinom(nrow(d),1,0.5)
##' fast.cox(Surv(entry,time,status)~X1*X2+strata(group)+cluster(id),data=d)
fast.cox <- function(formula,data,...) {
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  special <- c("strata", "cluster")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  Y <- model.extract(m, "response")
  if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
  if (ncol(Y)==2) {
    exit <- Y[,1]
    entry <- rep(0,nrow(Y))
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
  if (ncol(X)==0) return(NULL)
  fast.cox0(X,entry,exit,status,id,strata,...)
}


##' @S3method vcov fast.cox
vcov.fast.cox  <- function(object,...) {
  I <- -solve(object$hessian)
  if (!is.null(object$id)) {
    ii <- mets::cluster.index(object$id)
    UU <- matrix(nrow=ii$uniqueclust,ncol=ncol(I))
    for (i in seq(ii$uniqueclust)) {
      UU[i,] <- colSums(object$score[ii$idclustmat[i,seq(ii$cluster.size[i])]+1,,drop=FALSE])
    }
    J <- crossprod(UU)
  } else {
    J <- crossprod(object$score)
  }
  res <- I%*%J%*%t(I)
  colnames(res) <- rownames(res) <- names(coef(object))
  res
}

##' @S3method coef fast.cox
coef.fast.cox  <- function(object,...) {
  object$coef
}

##' @S3method summary fast.cox
summary.fast.cox <- function(object,...) {
  I <- -solve(object$hessian)
  V <- vcov(object)
  cc <- cbind(coef(object),diag(I)^0.5,diag(V)^0.5)
  colnames(cc) <- c("Estimate","Naive S.E.","Robust S.E.")
  rownames(cc) <- names(coef(object))
  res <- list(coef=cc)
  class(res) <- "summary.fast.cox"
  res
}

##' @S3method print summary.fast.cox
print.summary.fast.cox  <- function(x,...) {
  cat("\n")
  printCoefmat(x$coef,...)
  cat("\n")
}

##' @S3method print fast.cox
print.fast.cox  <- function(x,...) {
  print(summary(x),...)
}



