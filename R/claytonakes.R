##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title Clayton-Oakes model with piece-wise constant hazards
##' @param formula formula specifying the marginal proportional (piecewise constant) hazard structure with the right-hand-side being a survival object (Surv) specifying the entry time (optional), the follow-up time, and event/censoring status at follow-up. The clustering can be specified using the special function \code{id} (see example below).
##' @param data Data frame
##' @param id Variable defining the clustering (if not given in the formula)
##' @param var.formula Formula specifying the variance component structure (if not given via the id special function in the formula) using a linear model with log-link.
##' @param cuts Cut points defining the piecewise constant hazard
##' @param type Type of estimation (Clayton-Oakes or conditional frailty model)
##' @param start Optional starting values
##' @param control Control parameters to the optimization routine
##' @param var.link Link function for variance structure
##' @param ... Additional arguments
##' @author Klaus K. Holst
##' @examples
##' set.seed(1)
##' d <- simClaytonOakes(2000,4,2,1,stoptime=2,left=0.5)
##' e <- ClaytonOakes(Surv(lefttime,time,status)~x1+id(~1,cluster),cuts=c(0,0.5,1,2),data=subset(d,!truncated))
##' e
##' plot(e,add=FALSE)
##' @export
ClaytonOakes <- function(formula,data=parent.frame(),id,var.formula=~1,cuts=NULL,type="co",start,control=list(),var.link="log",...) {
  
  mycall <- match.call()
  formulaId <- Specials(formula,"id") 
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- "~."
  if (!is.null(formulaId)) {
    var.formulaId <- ~1
    if (length(formulaId)>1) {
      var.formula <- as.formula(formulaId[[1]])
      formulaId <- formulaId[[2]]
    }
    id <- formulaId
    mycall$id <- id
    formulaSt <- paste(formulaSt,paste("-id(",paste(var.formula,collapse=""),
                                       ",",formulaId,")"))
  }
  formulaSt <- paste(formulaSt,paste("-strata(",paste(formulaStrata,collapse="+"),")"))
  formula <- update(formula,formulaSt)

  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  if (missing(id)) stop("Missing 'id' variable")
  
  timevar <- terms(formula)[[2]]
  if (is.call(timevar)) {
    delayedentry <- (length(timevar)==4)*1
    entry <- NULL
    if (delayedentry==1)
      entry <- as.character(timevar[[2]])
    causes <- timevar[[3+delayedentry]]
    timevar <- timevar[[2+delayedentry]]
  }  
  timevar <- as.character(timevar)
  causes <- as.character(causes)
  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  X <- NULL
  nbeta <- 0  
  if (length(covars)>0) {
##    X <- model.matrix(as.formula(paste("~-1+",paste(covars,collapse="+"))),data)
    X <- model.matrix(update(formula,.~.+1),data)[,-1,drop=FALSE]
    nbeta <- ncol(X)
  }
  ngamma <- 0
  Z <- model.matrix(var.formula,data)
  ngamma <- ncol(Z)
  
  if (is.data.frame(data)) {
    mydata <- data.frame(T=data[,timevar],status=data[,causes],cluster=data[,id],entry=0)
    if (!is.null(entry)) {
      mydata$entry <- data[,entry]
    }
  } else {
    mydata <- data.frame(T=get(timevar,envir=data),status=get(causes,envir=data),cluster=get(id,envir=data),entry=0)
    if (!is.null(entry))      
      mydata$entry <- get(entry,envir=data)
  }
  if (is.null(cuts)) {
    cuts <- c(0,max(mydata$T))
  }
  if (max(mydata$T)>tail(cuts,1)) stop("Interval does not embed time observations")
  if (any(with(mydata, T<entry))) stop("Entry time occuring after event")

  ucluster <- unique(mydata$cluster)
  
  npar <- length(cuts)-1
  if (!is.null(X)) npar <- npar+ncol(X)
  npar <- npar+ncol(Z)
  p0 <- rep(0,npar)
  if (!missing(start)) p0 <- c(start,rep(0,max(0,length(npar)-length(start))))
  
  obj <- function(p) {
    varpar <- p[seq(ngamma)]
    p <- p[-seq(ngamma)]
    ##    theta0 <- rep(exp(varpar),length(ucluster));
    theta0 <- exp(Z%*%varpar)
    multhaz <- rep(1,nrow(mydata))
    if (!is.null(X)) {
      nbeta <- ncol(X)
      beta <- p[seq(nbeta)]
      p <- p[-seq(nbeta)]
      multhaz <- exp(X%*%beta)
    }
    ##browser()
    res <- .Call("claytonoakes",
           ds=mydata$status,ts=mydata$T,es=mydata$entry,
           allcs=mydata$cluster,cs=ucluster,
           cuts=cuts,hs=exp(p),mulths=multhaz,
           var=theta0,DUP=FALSE)$logLik
    return(-res)
  }
  opt <- tryCatch(nlminb(p0,obj,control=control),error=function(x) NULL)
  if (is.null(opt)) stop("Critical optmization problem")
  if (any(is.na(opt$par)) | any(!is.finite(opt$par)) | any(is.nan(opt$par)) ) {
    V <- matrix(NA,length(p0),length(p0))
  } else {    
    I <- hessian(obj,opt$par)
    ee <- tryCatch(eigen(I),error=function(x) NULL); 
    if (!is.null(ee)) {
      threshold <- 1e-12
      idx <- ee$values>threshold
      ee$values[idx] <- 1/ee$values[idx];
      if (!all(idx))
        ee$values[!idx] <- 0
      V <- with(ee, vectors%*%diag(values)%*%t(vectors))
    } else {
      V <- matrix(NA,length(p0),length(p0))
    }
  }
  res <- list(coef=opt$par,vcov=V,cuts=cuts,nbeta=nbeta,ngamma=ngamma,betanames=colnames(X),gammanames=colnames(Z),opt=opt)
  class(res) <- "claytonoakes"
  return(res)
}

##################################################

##' @S3method print claytonoakes
print.claytonoakes <- function(x,...) {
  print(summary(x))
}

##' @S3method print summary.claytonoakes
print.summary.claytonoakes <- function(x,...) {
  print(x$coef[,c(1,3,4)])
}

##' @S3method summary claytonoakes
summary.claytonoakes <- function(object,...) {
  mycoef <- matrix(nrow=length(object$coef),ncol=4)
  mycoef[,1:2] <- cbind(object$coef,sqrt(diag(object$vcov)))
  mycoef[,3:4] <- cbind(mycoef[,1]-qnorm(0.975)*mycoef[,2],mycoef[,1]+qnorm(0.975)*mycoef[,2])
  colnames(mycoef) <- c("Estimate","Std.Err","2.5%","97.5%")
  if (length(object$cuts))
  cutnames <- levels(cut(0,breaks=object$cuts))
  rownames(mycoef) <- c(paste("log-Var:",object$gammanames,sep=""),object$betanames,cutnames)
  mycoef[-seq(object$ngamma),] <- exp(mycoef[-seq(object$ngamma),])
  res <- list(coef=mycoef)
  class(res) <- "summary.claytonoakes"
  res
}

##' @S3method plot claytonoakes
plot.claytonoakes <- function(x,chaz=TRUE,add=!is.null(dev.list()),col="darkblue",...) {
  haz <- summary(x)$coef[-seq(x$nbeta+x$ngamma),,drop=FALSE]
  t <- x$cuts
  L <- approxfun(t,f=1,cumsum(c(0,haz[,1]*diff(t))),method="linear")
  if (add) {
    lines(t,L(t),col=col,...)
  } else {
    plot(t,L(t),type="l",col=col,...)
  }
  invisible(x)  
}
