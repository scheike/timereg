##' Wild bootstrap for Cox PH regression
##'
##' wild bootstrap for uniform bands for Cox models
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param offset offsets for cox model
##' @param weights weights for Cox score equations
##' @param B bootstraps 
##' @param type distribution for multiplier
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst, Thomas Scheike
##' @aliases phreg phreg.par robust.phreg readPhreg 
##' @examples
##' 
##'  n <- 100
##'  x <- 4*rnorm(n)
##'  time1 <- 2*rexp(n)/exp(x*0.3)
##'  time2 <- 2*rexp(n)/exp(x*(-0.3))
##'  status <- ifelse(time1<time2,1,2)
##'  time <- pmin(time1,time2)
##'  rbin <- rbinom(n,1,0.5)
##'  cc <-rexp(n)*(rbin==1)+(rbin==0)*rep(3,n)
##'  status <- ifelse(time < cc,status,0)
##'  time  <- ifelse(time < cc,time,cc)
##'  data <- data.frame(time=time,status=status,x=x)
##'  out <- pred.cif.boot(b1,b2,c1,c2,gplot=0)
##' 
##'  b1 <- Bootphreg(Surv(time,status==1)~x,data,B=0)
##'  b2 <- Bootphreg(Surv(time,status==2)~x,data,B=1000)
##'  c1 <- phreg(Surv(time,status==1)~x,data)
##'  c2 <- phreg(Surv(time,status==2)~x,data)
##' 
##'  ### exp to make all bootstraps positive
##'  out <- pred.cif.boot(b1,b2,c1,c2,gplot=0)
##' 
##'  cif.true <- (1-exp(-out$time))*.5
##'  with(out,plotl(time,cif,ylim=c(0,1)))
##'  lines(out$time,cif.true)
##'  with(out,mets:::plot.conf.region(time,band.EE,col=1))
##'  with(out,mets:::plot.conf.region(time,band.EE.log,col=3))
##'  with(out,mets:::plot.conf.region(time,band.EE.log.o,col=2))
##' 
##' @references
##' 
##' Wild bootstrap based confidence intervals for multiplicative hazards models, 
##' Dobler, Pauly, and Scheike (2018), 
##' 
##' @aliases pred.cif.boot 
##' @export
Bootphreg <- function(formula,data,offset=NULL,weights=NULL,B=1000,type=c("exp","poisson","normal"),...) {# {{{
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
    pos.cluster <- ts$terms
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  } else pos.cluster <- NULL
  if (!is.null(attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    pos.strata <- ts$terms
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
    strata.name <- ts$vars
  }  else { strata.name <- NULL; pos.strata <- NULL}
###  if (!is.null(attributes(Terms)$specials$offset)) {
###    ts <- survival::untangle.specials(Terms, "offset")
###    pos.offset <- ts$terms
###    Terms  <- Terms[-ts$terms]
###    offset <- m[[ts$vars]]
###  }  else pos.offset <- NULL
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  res <- Bootphreg01(X,entry,exit,status,id,strata,offset,weights,strata.name,B=B,
		     type=type,...)
###	   list(call=cl,model.frame=m,formula=formula,
###           strata.pos=pos.strata,cluster.pos=pos.cluster))
  class(res) <- "boot-phreg"
  
  res
}# }}}

###{{{ Bootpreg01

Bootphreg01 <- function(X,entry,exit,status,id=NULL,strata=NULL,offset=NULL,weights=NULL,
             strata.name=NULL,B=1000,type=c("normal","poisson","exp"),
	     cumhaz=TRUE,
             beta,stderr=TRUE,
	     method="NR",no.opt=FALSE,Z=NULL,propodds=NULL,AddGam=NULL,
	     case.weights=NULL,...) {
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

  ## possible casewights to use for bootstrapping and other things
  if (is.null(case.weights)) case.weights <- rep(1,length(exit)) 

  trunc <- (!is.null(entry))
  if (!trunc) entry <- rep(0,length(exit))

  id.orig <- id; 
  if (!is.null(id)) {
	  ids <- unique(id)
	  nid <- length(ids)
      if (is.numeric(id)) id <-  fast.approx(ids,id)-1 else  {
      id <- as.integer(factor(id,labels=seq(nid)))-1
     }
   } else id <- as.integer(seq_along(entry))-1; 

   system.time(dd <- .Call("FastCoxPrepStrata",
		     entry,exit,status,X, id, ### as.integer(seq_along(entry)),
		     trunc,strata,weights,offset,Zcall,case.weights,PACKAGE="mets"))

   dd$nstrata <- nstrata
   nj <- length(exit)

   res <- list()
   for (i in 1:B) {

      g <- switch(type[1],
                    "normal" =  {rnorm(nj)+1},
                    "exp"     = {rexp(nj)},
                    "poisson" = {rpois(nj,1)}
      )
      ## follows id that may be given in call with cluster(id)
      dd$caseweights[dd$id+1] <- g

   obj <- function(pp,U=FALSE,all=FALSE) {# {{{
		if (is.null(propodds) & is.null(AddGam)) 
	  val <- with(dd,
		   .Call("FastCoxPLstrata",pp,X,XX,sign,jumps,
		    strata,nstrata,weights,offset,ZX,caseweights,PACKAGE="mets"))
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
      val <- obj(beta,all=TRUE)
###      val[c("ploglik","gradient","hessian","U")] <- NULL
  }

  se.cumhaz <- lcumhaz <- lse.cumhaz <- NULL
  II <- NULL
  ### computes Breslow estimator 
### if (no.opt==FALSE & p!=0) II <- -solve(val$hessian) else II <- matrix(0,p,p)
 strata <- val$strata[val$jumps]
 nstrata <- val$nstrata
 jumptimes <- val$jumptimes

 ## Brewslow estimator
 cumhaz <- cbind(jumptimes,cumsumstrata(1/val$S0,strata,nstrata))
### varbetat <- 0
### if (no.opt==FALSE & p!=0) { 
###     DLambeta.t <- apply(val$E/c(val$S0),2,cumsumstrata,strata,nstrata)
###     varbetat <-   rowSums((DLambeta.t %*% II)*DLambeta.t)
### }
 ### covv <-  apply(covv*DLambeta.t,1,sum) Covariance is "0" by construction
### var.cumhaz <- cumsumstrata(1/val$S0^2,strata,nstrata)+varbetat
### se.cumhaz <- cbind(jumptimes,(var.cumhaz)^.5)

 colnames(cumhaz)    <- c("time","cumhaz")
### colnames(se.cumhaz) <- c("time","se.cumhaz")

 res[[i]] <- list(coef=cc,cumhaz=cumhaz[,2])
 }
 names(res) <- 1:B

 return(res)
}

###}}} phreg0

##' @export
pred.cif.boot <- function(b1,b2,c1,c2,gplot=1)
{# {{{
B <- length(b1)

times1 <- c1$cumhaz[,1]
times2 <- c2$cumhaz[,1]
coef1 <- do.call("rbind",lapply(b1,function(x) x$coef))
coef2 <- do.call("rbind",lapply(b2,function(x) x$coef))
###
###cums1 <- do.call("cbind",lapply(b1,function(x) Cpred(rbind(c(0,0),x$cumhaz),xx)[,2]))
###cums2 <- do.call("cbind",lapply(b2,function(x) Cpred(rbind(c(0,0),x$cumhaz),xx)[,2]))
bcums1 <- do.call("cbind",lapply(b1,function(x) x$cumhaz))
bcums2 <- do.call("cbind",lapply(b2,function(x) x$cumhaz))

where2 <- sindex.prodlim(c(0,times2),times1,strict=TRUE)
cums2 <- c(0,c2$cumhaz[,2])
cums2 <- cums2[where2]
cums1 <- c1$cumhaz[,2]

bcums2 <- rbind(0,bcums2)[where2,]

n <- length(times1)
cif1 <- cumsum( exp(-c(0,cums1[-n])-cums2)*diff(c(0,cums1)))

ccoef1 <- c1$coef
ccoef2 <- c2$coef
###

bcifs <- apply(exp(-rbind(0,bcums1[-n,])-bcums2)*apply(rbind(0,bcums1),2,diff),2,cumsum)

if (gplot==1) {
   matplot(times1,bcifs,type="s",lwd=0.2)
   lines(times1,cif1,type="s",lwd=2)
}

    cumx <- cif1
    ccums <- bcifs-cif1
    sdcumb <- apply(ccums,1,sd)
    zcums   <- ccums/sdcumb

    ### log-scale 
    lcum <- log(cif1)
    lccums <- log(bcifs)-lcum
    sdlogcumb <- apply(lccums,1,sd)
    zlogcums   <- lccums/sdlogcumb
  
    cumx.inv <- 1/cumx
    # in order to not divide by 0
    cumx.inv[cumx.inv == 0] <- 1
    
    # In fact: here we use the same quantiles independent of log or not log. 
    # Therefore: A division by cumx is required in the definition of the log bands.
    pcumsdb.EE <- percen(c(apply(abs(zcums),2,max, na.rm=TRUE)),0.95)
    pcumsdb.EE.log <- percen(apply(abs(zcums),2,max, na.rm=TRUE),0.95)
    pcumsdb.EE.log.o <- percen(apply(abs(zlogcums),2,max, na.rm=TRUE),0.95)

    band.EE <-     cbind( cumx - sdcumb * pcumsdb.EE , cumx + sdcumb * pcumsdb.EE)
    band.EE.log <- cbind( cumx*exp(- pcumsdb.EE.log * sdcumb * cumx.inv)  ,cumx*exp( pcumsdb.EE.log * sdcumb * cumx.inv ))
    band.EE.log.o <- cbind(exp(lcum - sdlogcumb* pcumsdb.EE.log.o), 
			   exp(lcum + sdlogcumb* pcumsdb.EE.log.o))

 return(list(time=times1,cif=cif1,
	     sdcif=sdcumb,sdlogcif=sdlogcumb,bcifs=bcifs,
	     band.EE=band.EE,band.EE.log=band.EE.log,
	     band.EE.log.o=band.EE.log.o))

}# }}}

