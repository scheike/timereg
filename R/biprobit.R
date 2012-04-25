##' Bivariate Probit model
##'
##' .. content for \details{} ..
##' @title Bivariate Probit model
##' @param formula 
##' @param data 
##' @param id 
##' @param time 
##' @param strata 
##' @param eqmarg 
##' @param indep 
##' @param weight 
##' @param biweight 
##' @param samecens 
##' @param randomeffect 
##' @param vcov 
##' @param pairsonly 
##' @param allmarg 
##' @param control 
##' @param bound 
##' @param messages 
##' @param p 
##' @param ... 
##' @return \code{biprobit}
##' @author Klaus K. Holst
##' @keywords regression
##' @export
biprobit <- function(formula, data, id, time, strata=NULL, eqmarg=TRUE,
                     indep=FALSE, weight=NULL,
                     biweight=function(x) 1/min(x),
                     samecens=TRUE, randomeffect=FALSE, vcov="robust",
                     pairsonly=FALSE,
                     allmarg=samecens&!is.null(weight),
                     control=list(trace=0,method="qausi"),
                     bound=FALSE,
                     messages=1,
                     p,...) {

  mycall <- match.call()
  formulaId <- Specials(formula,"cluster")
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- paste("~.-cluster(",formulaId,")-strata(",paste(formulaStrata,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) strata <- formulaStrata
  mycall$formula <- formula

  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    fit <- lapply(seq(length(dd)),function(i) {
      if (messages>0) message("Strata '",names(dd)[i],"'")
      mycall$data <- dd[[i]]
      eval(mycall)
    })
    res <- list(model=fit)
    res$strata <- names(res$model) <- names(dd)
    class(res) <- c("biprobit.strata","biprobit")
    res$coef <- unlist(lapply(res$model,coef))
    res$vcov <- blockdiag(lapply(res$model,vcov.biprobit))
    res$N <- length(dd)
    res$idx <- seq(length(coef(res$model[[1]])))
    rownames(res$vcov) <- colnames(res$vcov) <- names(res$coef)
    return(res)
  }
  
  if (missing(id)) {    
    if (!is.null(weight)) {
      weights <- data[,weight]
      return(glm(formula,data=data,family=binomial(probit),weights=weights,...))
    }
    return(glm(formula,data=data,family=binomial(probit),...))    
  }

  mycall <- match.call()
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  if (pairsonly) {
    data <- data[which(as.character(data[,id])%in%names(idtab)[idtab==2]),]
    idtab <- table(data[,id])
  }
  
  if (missing(time)) {
    time <- "time"
    while (time%in%names(data)) time <- paste(time,"_",sep="")
    data[,time] <- unlist(lapply(idtab,seq))
  }
  ff <- paste(as.character(formula)[3],"+",time,"+",id)
  yvar <- paste(deparse(formula[[2]]),collapse="")
  if (!is.null(weight))
    ff <- paste(weight,"+",ff)
  ff <- paste("~",yvar,"+",ff)

  if (is.logical(data[,yvar])) data[,yvar] <- data[,yvar]*1
  if (is.factor(data[,yvar])) data[,yvar] <- as.numeric(data[,yvar])-1  
##  Y <- cbind(as.numeric(data[,yvar]))-(!is.numeric(data[,yvar]))
  
  formula0 <- as.formula(ff)
  Data <- model.matrix(formula0,data,na.action=na.pass)
  rnames1 <- setdiff(colnames(Data),c(yvar,time,id,weight))
  X0 <- as.matrix(Data[,rnames1])
  
  nx <- length(rnames1)
  if (nx==0) stop("Zero design not allowed")
  
  ##  data0 <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  ##  data0 <- data
  ##  data0 <- data0[order(data0[,id]),]
  midx1 <- seq(nx)
  midx2 <- midx1+nx
  midx <- seq(2*nx)
  plen <- ifelse(eqmarg,nx+1,2*nx+1)

  Wide <- reshape(as.data.frame(Data),idvar=id,timevar=time,direction="wide")
  W0 <- NULL
  yidx <- paste(yvar,1:2,sep=".")
  rmidx <- c(id,yidx)
  if (!is.null(weight)) {
    W <- cbind(data[,weight])
    widx <- paste(weight,1:2,sep=".")
    W0 <- as.matrix(Wide[,widx])
    rmidx <- c(rmidx,widx)
  }
  Y0 <- as.matrix(Wide[,yidx])
  XX0 <- as.matrix(Wide[,setdiff(colnames(Wide),rmidx)])
  XX0[is.na(XX0)] <- 0

  datanh <- function(r) 1/(1-r^2)
  dtanh <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
  vartr <- tanh
  dvartr <- dtanh; varitr <- atanh
  trname <- "tanh"; itrname <- "atanh"    
  Sigma1 <- diag(2)  
  Sigma2 <- matrix(c(0,1,1,0),2,2)
  dS0 <- rbind(c(0,1,1,0))
  varcompname <- "Tetrachoric correlation"
  msg <- "Variance of latent residual term = 1 (standard probit link)"
  if (randomeffect) {
    dS0 <- rbind(rep(1,4))
    vartr <- dvartr <- exp; inv <- log
    trname <- "exp"; itrname <- "log"
    Sigma2 <- 1
    varcompname <- NULL
  }
  model <- list(tr=vartr,name=trname,inv=itrname,invname=itrname,deriv=dvartr,varcompname=varcompname,dS=dS0,eqmarg=eqmarg)

  MyData <- ExMarg(Y0,XX0,W0,dS0,midx1,midx2,eqmarg=eqmarg,allmarg=allmarg)
  if (samecens & !is.null(weight)) {
    MyData$W0 <- cbind(apply(MyData$W0,1,biweight))
    if (!is.null(MyData$Y0_marg)) {
      MyData$W0_marg <- cbind(apply(MyData$W0_marg,1,biweight))
    }
  }
  
  SigmaFun <- function(p,...) {
    Sigma <- Sigma1+Sigma2*vartr(p[plen])
    if (indep) Sigma <- diag(2)
    return(Sigma)
  }

  U <- function(p,indiv=FALSE) {
    if (bound) p[plen] <- min(p[plen],20)
    Sigma <- SigmaFun(p)
    lambda <- eigen(Sigma)$values
    if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
    Mu_marg <- NULL
    if (eqmarg) {
      B <- cbind(p[midx1])
      Mu <- with(MyData,
                 cbind(XX0[,midx1,drop=FALSE]%*%B,XX0[,midx2,drop=FALSE]%*%B))     
##      Mu <- with(MyData, matrix(X0%*%B,ncol=2,byrow=TRUE))
      if (!is.null(MyData$Y0_marg)) 
        Mu_marg <- with(MyData, XX0_marg%*%B)
    } else {
      B1 <- cbind(p[midx1])
      B2 <- cbind(p[midx2])
      Mu <- with(MyData,
                 cbind(XX0[,midx1,drop=FALSE]%*%B1,XX0[,midx2,drop=FALSE]%*%B2))
      if (!is.null(MyData$Y0_marg))
        Mu_marg <- with(MyData, rbind(X0_marg1%*%B1,X0_marg2%*%B2))
    }

    U <- with(MyData, .Call("biprobit2",
                             Mu,XX0,
                             Sigma,dS0*dvartr(p[plen]),Y0,W0,
                             !is.null(W0),TRUE,eqmarg))
    
    if (!is.null(MyData$Y0_marg)) {
      ## U_marg <- uniprobit(Mu_marg[,1],XX0_marg,
      ##                     Sigma[1,1],dS0_marg*dvartr(p[plen]),Y0_marg,
      ##                     W0_marg,indiv=TRUE)
      ## U$score <- rbind(U$score,U_marg)
      ## U$loglik <- c(U$loglik,attributes(U_marg)$logLik)
      ##      W0_marg <- rep(1,nrow(XX0_marg))
      ##      browser()
      U_marg <- with(MyData, .Call("uniprobit",
                                   Mu_marg,XX0_marg,
                                   Sigma[1,1],dS0_marg*dvartr(p[plen]),Y0_marg,
                                   W0_marg,!is.null(W0_marg),TRUE))
      U$score <- rbind(U$score,U_marg$score)
      U$loglik <- c(U$loglik,U_marg$loglik)
    }

    if (indiv) {
      val <- U$score[MyData$id,,drop=FALSE]
      N <- length(MyData$id)
      idxs <- seq_len(N)
      for (i in seq_len(N)) {
        idx <- which((MyData$idmarg)==(MyData$id[i]))+N
        idxs <- c(idxs,idx)
        val[i,] <- val[i,]+colSums(U$score[idx,,drop=FALSE])
      }
      val <- rbind(val, U$score[-idxs,,drop=FALSE])
      attributes(val)$logLik <- U$loglik
      return(val)
    }
    val <- colSums(U$score)
    attributes(val)$logLik <- sum(U$loglik)
    return(val)
  }
  
  p0 <- rep(0,plen)
  if (!is.null(control$start)) {
    p0 <- control$start
    control$start <- NULL
  }
##  browser()

  if (!missing(p)) return(U(p,indiv=FALSE))

  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))
  h0 <- function(p) crossprod(U(p,indiv=TRUE))

  if (is.null(control$method)) {
    ##    control$method <- ifelse(samecens & !is.null(weight), "bhhh","quasi")
    control$method <- "quasi"
  }
  control$method <- tolower(control$method)
  if (control$method=="score") {
    control$method <- NULL
    op <- nlminb(p0,f,control=control,...)
  } else if (control$method=="quasi") {
    control$method <- NULL
    op <- nlminb(p0,f0,gradient=g0,control=control,...)
  } else if (control$method=="bhhh") {
    controlnr <- list(stabil=FALSE,
                      gamma=0.1,
                      gamma2=1,
                      ngamma=5,
                      iter.max=200,
                      epsilon=1e-12,
                      tol=1e-9,
                      trace=1,
                      stabil=FALSE)
    controlnr[names(control)] <- control
    op <- lava:::NR(start=p0,NULL,g0, h0,control=controlnr)
  } else {
    control$method <- NULL
    op <- nlminb(p0,f0,control=control,...)
  }
  UU <- U(op$par,indiv=TRUE)
  J <- crossprod(UU)
  ##  iJ <- Inverse(J)
  iI <- Inverse(-numDeriv::jacobian(U,op$par))
  V <- switch(vcov,
              robust=,
              sandwich=iI%*%J%*%iI,##iJ%*%I%*%iJ,
              score=,
              outer=Inverse(J),
              hessian=iI              
              )
  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  rho <- ifelse(itrname=="log","U","rho")
  if (!eqmarg)
    rownames(cc) <- c(paste(rnames1,rep(c(1,2),each=length(rnames1)),sep="."),
                      paste(itrname,"(",rho,")",sep=""))
  else
    rownames(cc) <- c(rnames1,paste(itrname,"(",rho,")",sep=""))
  rownames(V) <- colnames(V) <- rownames(cc)
  npar <- list(intercept=attributes(terms(formula))$intercept,
              pred=nrow(attributes(terms(formula))$factor)-1)
  if (!eqmarg) npar <- lapply(npar,function(x) x*2)
  npar$var <- 1##nrow(cc)-sum(unlist(npar))
  N <- with(MyData, c(n=nrow(XX0)*2+length(margidx), pairs=nrow(XX0)))
  val <- list(coef=cc,N=N,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall, model=model,msg=msg,npar=npar,
              SigmaFun=SigmaFun)
  class(val) <- "biprobit"
  return(val)
}

