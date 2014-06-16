##' @export
biprobit <- function(x,...) UseMethod("biprobit")

##' @export
biprobit.table <- function(x,eqmarg=TRUE,...) {
    x
}

##' @export
biprobit.default <- function(x,id,
                             weight=NULL,biweight=function(x) { u <- min(x); ifelse(u==0,0,1/min(u,1e-2)) },
                             eqmarg=TRUE,ref,...) {
    
    x <- as.factor(x)
    lev <- levels(x)
    if (missing(ref)) ref <- lev[2]
    if (!is.null(weight)) {
        dd <- na.omit(fast.reshape(data.frame(x,weight,id),id=id))
    } else {
        dd <- na.omit(fast.reshape(data.frame(x,id),id=id))
    }
    M <- table(dd[,c("x1","x2")])
    ## if (is.null(weight)) return(biprobit(M))
    
    Y0 <- matrix(c(0,0, 0,1, 1,0, 1,1),ncol=2,byrow=TRUE)
    XX0 <- matrix(1,ncol=2,nrow=4)
    MyData <- ExMarg(Y0,XX0,W0=NULL,NULL,
                     midx1=1,midx2=2,eqmarg=eqmarg,allmarg=FALSE)

    datanh <- function(r) 1/(1-r^2)
    dtanh <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
    vartr <- tanh
    dvartr <- dtanh; varitr <- atanh
    trname <- "tanh"; itrname <- "atanh"    
    dS <- rbind(c(0,1,1,0))
    varcompname <- "Tetrachoric correlation"    
    ##msg <- "Variance of latent residual sterm = 1 (standard probit link)"
    msg <- NULL
    model <- list(tr=vartr,name=trname,inv=itrname,invname=itrname,deriv=dvartr,varcompname=varcompname,dS=dS,eqmarg=eqmarg,...)
    SigmaFun <- function(p,...) {
        Sigma1 <- diag(2)  
        Sigma2 <- matrix(c(0,1,1,0),2,2)
        Sigma <- Sigma1+Sigma2*vartr(p)
        ##if (indep) Sigma <- diag(2)
        attributes(Sigma)$dvartr <- dvartr
        return(Sigma)
    }

    if (length(weight)<2) {
        w0 <- as.vector(t(M))
    } else {
        w <- apply(dd[,c("weight1","weight2")],1,biweight)
        pos <- (dd[,"x1"]==ref)*1+(dd[,"x2"]==ref)*2+1    
        w0 <- as.vector(by(w,pos,sum))
    }
    p0 <- c(rep(M[2,1]/(M[2,1]+M[1,1]),1+!eqmarg),0.5)
    
    U <- function(p,w0) {
        val <- Ubiprobit(p,SigmaFun,dS,eqmarg,1,MyData,indiv=TRUE)
        logl <- w0*attributes(val)$logLik
        score <- apply(val,2,function(x) w0*x)
        return(structure(score,logLik=logl))
    }
    f0 <- function(p) -sum(attributes(U(p,w0))$logLik)
    g0 <- function(p) -as.numeric(colSums(U(p,w0)))
    op <- nlminb(p0,f0,gradient=g0,...)

    iI <- Inverse(numDeriv::jacobian(g0,op$par))
    V <- iI
    UU <- U(op$par,w0)
    logLik <- sum(attributes(UU)$logLik)
    if (length(weight)>1) {
        UU <- U(op$par,1)
        UU <- apply(UU[pos,],2,function(x) w*x)
        meat <- crossprod(UU)
        V <- iI%*%meat%*%iI
    }    
    
    mycall <- match.call()    
    cc <- cbind(op$par,sqrt(diag(V)))
    cc <- cbind(cc,cc[,1]/cc[,2],2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
    colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
    nn <- c("Intercept")
    if (!eqmarg) nn <- paste(nn,1:2,sep=".")
    rhonam <- ifelse(itrname=="log","U","rho")
    nn <- c(nn,paste(itrname,"(",rhonam,")",sep=""))
    rownames(cc) <- nn   
    
    val <- list(coef=cc,
                N=c(n=2*sum(M),pairs=sum(M)),
                vcov=V,
                bread=iI,
                score=rbind(UU),
                logLik=logLik,
                opt=op,
                call=mycall,
                model=model,
                msg=msg,
                table=M,
                npar=list(intercept=1+!eqmarg,pred=NULL,var=1),
                SigmaFun=SigmaFun)
    class(val) <- "biprobit"
    return(val)
}


##' @export
biprobit.formula <- function(x, data, id, rho=~1, num=NULL, strata=NULL, eqmarg=TRUE,
                             indep=FALSE, weight=NULL,
                             biweight=function(x) 1/min(x),
                             samecens=TRUE, randomeffect=FALSE, vcov="robust",
                             pairsonly=FALSE,                             
                             allmarg=samecens&!is.null(weight),
                             control=list(trace=0,method="quasi"),
                             constrain,
                             bound=FALSE,
                             messages=1,
                             table=TRUE,
                             p,...) {

  mycall <- match.call()
  formulaId <- unlist(Specials(x,"cluster"))
  if (is.null(formulaId)) {
    formulaId <- unlist(Specials(x,"id"))
  }
  formulaStrata <- unlist(Specials(x,"strata"))
  formulaSt <- paste("~.-cluster(",formulaId,
                     ")-id(",formulaId,
                     ")-strata(",paste(formulaStrata,collapse="+"),")")
  formula <- update(x,formulaSt)
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
    class(res) <- c("twinlm.strata","biprobit")
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

  yx <- getoutcome(formula)
  
  if (length(attributes(yx)$x)==0 && pairsonly && table) {
      if (!is.null(weight)) weight <- data[,weight]      
      return(biprobit(data[,yx],data[,id],weight,biweight=biweight,eqmarg=eqmarg,...))
  }
  
  mycall <- match.call()
  DD <- procdatabiprobit(formula,data,id,num=num,weight=weight,pairsonly=pairsonly,rho,...)
  rnames1 <- DD$rnames1

  nx <- length(rnames1)
  if (nx==0) stop("Zero design not allowed")  
  midx1 <- seq(nx)
  midx2 <- midx1+nx
  midx <- seq(2*nx)
  blen <- ifelse(eqmarg,nx,2*nx)
  zlen <- ncol(DD$Z0)      
  plen <- blen+zlen
  
  datanh <- function(r) 1/(1-r^2)
  dtanh <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
  vartr <- tanh
  dvartr <- dtanh; varitr <- atanh
  trname <- "tanh"; itrname <- "atanh"    
  Sigma1 <- diag(2)  
  Sigma2 <- matrix(c(0,1,1,0),2,2)
  dS0 <- rbind(c(0,1,1,0))
  
  varcompname <- "Tetrachoric correlation"
  msg <- NULL
  if (randomeffect) msg <- "Variance of latent residual term = 1 (standard probit link)"
  if (randomeffect) {
      dS0 <- rbind(rep(1,4))
      vartr <- dvartr <- exp; inv <- log
      trname <- "exp"; itrname <- "log"
      Sigma2 <- 1
      varcompname <- NULL
  }
  model <- list(tr=vartr,name=trname,inv=itrname,invname=itrname,deriv=dvartr,varcompname=varcompname,dS=dS0,eqmarg=eqmarg,randomeffect=randomeffect,blen=blen,zlen=zlen)

  MyData <- with(DD,ExMarg(Y0,XX0,W0,dS0,midx1,midx2,eqmarg=eqmarg,allmarg=allmarg,Z0))
  if (samecens & !is.null(weight)) {
      MyData$W0 <- cbind(apply(MyData$W0,1,biweight))
      if (!is.null(MyData$Y0_marg)) {
          MyData$W0_marg <- cbind(apply(MyData$W0_marg,1,biweight))
      }
  }
  
  SigmaFun <- function(p,Z=MyData$Z0,cor=!randomeffect,...) {
      if (!cor) {
          r <- vartr(p[1])
          Sigma <- Sigma1+Sigma2*vartr(p)
          if (indep) Sigma <- diag(2)
          attributes(Sigma)$dvartr <- dvartr
          return(Sigma)
      }
      val <- Z%*%p
      dr <- apply(Z,2,function(x) x*dvartr(val))
      structure(list(rho=vartr(val),lp=val,drho=dr),dvartr=dvartr,vartr=vartr)
  }

  U <- function(p,indiv=FALSE) {
      gamma <- p[seq(zlen)+blen]
      if (bound) gamma <- min(gamma,20)
      Sigma <- SigmaFun(gamma)
      if (randomeffect) {
          lambda <- eigen(Sigma)$values
          if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
      }
      Mu_marg <- NULL
      if (eqmarg) {
          B <- cbind(p[midx1])
          Mu <- with(MyData,
                     cbind(XX0[,midx1,drop=FALSE]%*%B,XX0[,midx2,drop=FALSE]%*%B))     
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

      if (randomeffect) {
          U <- with(MyData, .Call("biprobit2",
                                  Mu,XX0,
                                  Sigma,dS0*attributes(Sigma)$dvartr(p[plen]),Y0,W0,
                                  !is.null(W0),TRUE,eqmarg,FALSE))
      } else {
          U <- with(MyData, .Call("biprobit2",
                                  Mu,XX0,
                                  Sigma$rho,Sigma$drho,
                                  Y0,W0,
                                  !is.null(W0),TRUE,eqmarg,TRUE))
      }
      
      if (!is.null(MyData$Y0_marg)) {
          if (randomeffect) {
              U_marg <- with(MyData, .Call("uniprobit",
                                           Mu_marg,XX0_marg,
                                           Sigma[1,1],dS0_marg*attributes(Sigma)$dvartr(p[plen]),Y0_marg,
                                           W0_marg,!is.null(W0_marg),TRUE))
          } else {
              U_marg0 <- matrix(0,length(MyData$Y0_marg),ncol=plen)
              U_marg <- with(MyData, .Call("uniprobit",
                                           Mu_marg,XX0_marg,
                                           1,matrix(ncol=0,nrow=0),Y0_marg,
                                           W0_marg,!is.null(W0_marg),TRUE))
              U_marg0[,seq(blen)] <-  U_marg[[1]]
              U_marg[[1]] <- U_marg0
          }
          U$score <- rbind(U$score,U_marg$score)
          U$loglik <- c(U$loglik,U_marg$loglik)
      }
      
      if (indiv) {
          ## val <- with(MyData, cluster.index(c(id,idmarg),mat=U$score))
          ## attributes(val)$logLik <- U$loglik
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
  } else {
      g <- suppressWarnings(glm(formula,data,family=binomial(probit)))
      p0[midx1] <- coef(g)
      if (!eqmarg) p0[midx2] <- coef(g)    
  }

  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))
  h0 <- function(p) crossprod(U(p,indiv=TRUE))

  if (is.null(control$method)) {
    ##    control$method <- ifelse(samecens & !is.null(weight), "bhhh","quasi")
    control$method <- "quasi"
  }
  control$method <- tolower(control$method)
  
  if (missing(p)) {
      if (control$method=="score") {
          control$method <- NULL
          op <- nlminb(p0,f,control=control,...)
      } else if (control$method=="quasi") {
          control$method <- NULL
          op <- nlminb(p0,f0,gradient=g0,control=control,...)
          ## }
  ## else if (control$method=="bhhh") {
  ##   controlnr <- list(stabil=FALSE,
  ##                     gamma=0.1,
  ##                     gamma2=1,
  ##                     ngamma=5,
  ##                     iter.max=200,
  ##                     epsilon=1e-12,
  ##                     tol=1e-9,
  ##                     trace=1,
  ##                     stabil=FALSE)
  ##   controlnr[names(control)] <- control
  ##   op <- lava:::NR(start=p0,NULL,g0, h0,control=controlnr)
      } else {
          control$method <- NULL
          op <- nlminb(p0,f0,control=control,...)
      }
  } else op <- list(par=p)
  
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
  cc <- cbind(cc,cc[,1]/cc[,2],2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")

  p1 <- "("; p2 <- ")"
  if (itrname=="log") rhonam <- "U" else {
      rhonam <- DD$znames
      p1 <- p2 <- ""; itrname <- "r:" 
  }
  if (!eqmarg)
    rownames(cc) <- c(paste(rnames1,rep(c(1,2),each=length(rnames1)),sep="."),
                      paste(itrname,p1,rhonam,p2,sep=""))
  else
    rownames(cc) <- c(rnames1,paste(itrname,p1,rhonam,p2,sep=""))
  rownames(V) <- colnames(V) <- rownames(cc)
  npar <- list(intercept=attributes(terms(formula))$intercept,
              pred=nrow(attributes(terms(formula))$factor)-1)
  if (!eqmarg) npar <- lapply(npar,function(x) x*2)
  npar$var <- 1##nrow(cc)-sum(unlist(npar))
  N <- with(MyData, c(n=nrow(XX0)*2+length(margidx), pairs=nrow(XX0)))
  val <- list(coef=cc,N=N,vcov=V,bread=iI,score=UU,logLik=attributes(UU)$logLik,opt=op, call=mycall, model=model,msg=msg,npar=npar,
              SigmaFun=SigmaFun)
  class(val) <- "biprobit"
  return(val)
}



procdatabiprobit <- function(formula,data,id,num=NULL,weight=NULL,pairsonly=FALSE,rho=~1,...) {

    data <- data[order(data[,id]),]
    idtab <- table(data[,id])
    if (pairsonly) {
        data <- data[which(as.character(data[,id])%in%names(idtab)[idtab==2]),]
        idtab <- table(data[,id])
    }
  
    ff <- paste(as.character(formula)[3],"+",
                paste(c(id,num),collapse="+"))
    yvar <- paste(deparse(formula[[2]]),collapse="")
    if (!is.null(weight))
        ff <- paste(weight,"+",ff)
    ff <- paste("~",yvar,"+",ff)

    if (is.logical(data[,yvar])) data[,yvar] <- data[,yvar]*1
    if (is.factor(data[,yvar])) data[,yvar] <- as.numeric(data[,yvar])-1  
                                        
    formula0 <- as.formula(ff)  
    opt <- options(na.action="na.pass")
    Data <- model.matrix(formula0,data)
    options(opt)
    rnames1 <- setdiff(colnames(Data),c(yvar,num,id,weight))
    X0 <- as.matrix(Data[,rnames1])

    ex <- 1+!is.null(num)
    rho <- update(rho,paste("~.+",
                           paste(c(id,num),collapse="+")))
    Z0 <- model.matrix(rho,data);
    znames <- setdiff(colnames(Z0),c(id,num))
    Z0 <- as.matrix(subset(fast.reshape(Z0,id=id),select=-c(id,num))[,seq(ncol(Z0)-ex),drop=FALSE])
    Wide <- fast.reshape(as.data.frame(Data),id=id,num=num,sep=".",labelnum=TRUE)
    
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

    list(Y0=Y0,XX0=XX0,W0=W0,Z0=Z0,znames=znames,rnames1=rnames1)
}



Ubiprobit <- function(p,S,dS,eqmarg,nx,MyData,indiv=FALSE) {
    midx1 <- seq(nx)
    midx2 <- midx1+nx
    midx <- seq(2*nx)
    plen <- ifelse(eqmarg,nx+1,2*nx+1)
    Sigma <- S(p[plen])
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
                            Sigma,dS*attributes(Sigma)$dvartr(p[plen]),Y0,W0,
                            !is.null(W0),TRUE,eqmarg,FALSE))
    
    if (!is.null(MyData$Y0_marg)) {
        U_marg <- with(MyData, .Call("uniprobit",
                                     Mu_marg,XX0_marg,
                                     Sigma[1,1],dS0_marg*attributes(Sigma)$dvartr(p[plen]),Y0_marg,
                                     W0_marg,!is.null(W0_marg),TRUE))
        U$score <- rbind(U$score,U_marg$score)
        U$loglik <- c(U$loglik,U_marg$loglik)
    }

    if (indiv) {
        val <- with(MyData, cluster.index(c(id,idmarg),mat=U$score))
        attributes(val)$logLik <- U$loglik
        return(val)
    }
    val <- colSums(U$score)
    attributes(val)$logLik <- sum(U$loglik)
    return(val)
}



##' @export
biprobit.time <- function(formula,data,id,...,
                          breaks=Inf,pairsonly=TRUE,
                          cens.formula,cens.model="aalen",weight="w") {

    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

    censtime <- model.extract(M, "response")
    status <- censtime[,2]
    time <- censtime[,1]
    outcome <- as.character(terms(formula)[[2]])    
    if (is.null(breaks)) breaks <-  quantile(time,c(0.25,0.5,0.75,1))

    outcome0 <- paste(outcome,"_dummy")
    res <- list()
    for (tau in breaks) {
        if (length(breaks)>1) message(tau)
        data0 <- data
        time0 <- time
        cond0 <- time0>tau
        status0 <- status
        status0[cond0 & status==1] <- 3 ## Censored
        data0[cond0,outcome] <- FALSE
        time0[cond0] <- tau
        data0$S <- Surv(time0,status0==1)        
        dataw <- ipw(update(cens.formula,S~.), data=data0, cens.model=cens.model,
                     cluster=id,weightname=weight,obsonly=TRUE)
        suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weight=weight, pairsonly=pairsonly,...))
        res <- c(res,list(summary(b)))
    }
    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,coef=lapply(res,function(x) x$all),summary=res,call=m,type="time")
    class(res) <- "multitwinlm"
    return(res)    
}

biprobit.time2 <- function(formula,data,id,...,
                          breaks=Inf,pairsonly=TRUE,
                          cens.formula,cens.model="aalen",weight="w") {

 
    m <- match.call(expand.dots = TRUE)[1:3]
    Terms <- terms(cens.formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m)
    censtime <- model.extract(M, "response")
    status <- censtime[,2]
    time <- censtime[,1]
    outcome <- as.character(terms(formula)[[2]])    
    if (is.null(breaks)) breaks <-  quantile(time,c(0.25,0.5,0.75,1))

    outcome0 <- paste(outcome,"_dummy")
    res <- list()
    for (tau in breaks) {
        if (length(breaks)>1) message(tau)
        data0 <- data
        time0 <- time
        cond0 <- time0>tau
        status0 <- status
        status0[cond0 & status==1] <- 3 ## Not censored anymore (status=1,censored in original sample) 
        data0[cond0,outcome] <- FALSE 
        ##        time0[cond0] <- tau
        time0[cond0 & status==1] <- tau
        data0$S <- Surv(time0,status0==1)        
        dataw <- ipw(update(cens.formula,S~.), data=data0, cens.model=cens.model,
                     cluster=id,weightname=weight,obsonly=TRUE)
        suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weight=weight, pairsonly=pairsonly,...))
        res <- c(res,list(summary(b)))
    }
    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,coef=lapply(res,function(x) x$all),summary=res,call=m,type="time")
    class(res) <- "multitwinlm"
    return(res)    
}
