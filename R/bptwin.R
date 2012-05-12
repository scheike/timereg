##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title Liability model for twin data
##' @param formula 
##' @param data 
##' @param id 
##' @param zyg 
##' @param DZ 
##' @param OS 
##' @param weight 
##' @param biweight 
##' @param strata 
##' @param messages 
##' @param control 
##' @param type 
##' @param eqmean 
##' @param param 
##' @param pairsonly 
##' @param stderr 
##' @param robustvar 
##' @param p 
##' @param indiv 
##' @param constrain 
##' @param samecens 
##' @param allmarg 
##' @param bound 
##' @param debug 
##' @param ... 
##' @author Klaus K. Holst
##' @export
bptwin <- function(formula, data, id, zyg, DZ, OS,
                   weight=NULL,
                   biweight=function(x) 1/min(x),
                   strata=NULL,
                   messages=1,
                   control=list(trace=0),
                   type="ace",
                   eqmean=TRUE,
                   param=0,
                   pairsonly=FALSE,
                   stderr=TRUE,                  
                   robustvar=TRUE,                   
                   p, indiv=FALSE,
                   constrain,
                   samecens=TRUE,
                   allmarg=samecens&!is.null(weight),
                   bound=FALSE,
                   debug=FALSE,...) {

###{{{ setup

  mycall <- match.call()
  formulaId <- unlist(Specials(formula,"cluster"))
  formulaStrata <- unlist(Specials(formula,"strata"))
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
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
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
  }

##################################################
### No strata
    if (is.null(control$method)) {
    control$method <- "gradient"
    if (!samecens & !is.null(weight)) control$method <- "bhhh"
  }
  if (length(grep("flex",tolower(type)))>0) { type <- "u"; eqmean <- FALSE }

  yvar <- paste(deparse(formula[[2]]),collapse="")
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  if (pairsonly)
    data <- data[as.character(data[,id])%in%names(idtab)[idtab==2],]
  if (is.logical(data[,yvar])) data[,yvar] <- data[,yvar]*1
  if (is.factor(data[,yvar])) data[,yvar] <- as.numeric(data[,yvar])-1  

  ## Y <- suppressWarnings(model.frame(formula,data,na.action="na.pass")[,1])    
  ## Blen <- 0
  ## Bord <- c()
  ## if (degree>0)
  ##   if (is.Surv(Y) | !missing(time) | !is.null(B)) {
  ##     if (is.Surv(Y)) {
  ##       if (ncol(Y)==3) {
  ##         entry <- Y[,1]
  ##         time <- Y[,2]
  ##         Y <- Y[,3]*1 # event
  ##       } else {
  ##         time <- Y[,1]
  ##         Y <- Y[,2]*1
  ##       }
  ##     } else {
  ##       if (is.character(time))
  ##         time <- data[,time]
  ##       if (is.character(entry))
  ##         entry <- data[,entry]
  ##     }
  ##     if (is.null(B)) {
  ##       B <- bs(time,degree=degree)
  ##       Bord <- cbind(time,B)[order(time),,drop=FALSE]
  ##       colnames(B) <- paste("B","_",1:ncol(B),sep="")
  ##       Blen <- ncol(B)   
  ##     }
  ##   }

  idtab <- table(data[,id])
  ##  ii0 <- which(as.character(data[,id])%in%names(idtab)[idtab==2])              
  ##  data0 <- data[ii0,]  
  idx2 <- NULL
  if (!missing(DZ)) {
    if (!missing(OS))
      idx2 <- data[,zyg]==OS
    idx1 <- data[,zyg]==DZ
    idx0 <- data[,zyg]!=DZ
    data[,zyg] <- (data[,zyg]!=DZ)*1
  } else {
    if (!missing(OS))
      idx2 <- (as.factor(data[,zyg])==levels(as.factor(data[,zyg]))[3]) # DZos
    DZlev <- levels(as.factor(data[,zyg]))[1]
    idx1 <- (as.factor(data[,zyg])==DZlev) # DZ
    idx0 <- (as.factor(data[,zyg])==levels(as.factor(data[,zyg]))[2]) # MZ
    message("Using '",DZlev,"' as DZ",sep="")
    data[,zyg] <- (data[,zyg]!=levels(as.factor(data[,zyg]))[1])
  }
  
  
  time <- "time"
  while (time%in%names(data)) time <- paste(time,"_",sep="")
  data[,time] <- unlist(lapply(idtab,seq))
  
  ##  ff <- as.formula(paste("~",paste(attributes(terms(formula))$term.labels,
  ##                                   collapse="+"),"+1",sep=""))
  ff <- paste(as.character(formula)[3],"+",time,"+",id,"+",zyg)
  if (!is.null(weight))
    ff <- paste(weight,"+",ff)
  ff <- paste("~",yvar,"+",ff)
  formula0 <- as.formula(ff)
  Data <- model.matrix(formula0,data,na.action=na.pass)
  rnames1 <- setdiff(colnames(Data),c(yvar,time,id,weight,zyg))
  ##  X0 <- as.matrix(Data[,rnames1])
  nx <-length(rnames1) ##ncol(Data)-4 + !is.null(weight)
  if (nx==0) stop("Zero design not allowed")
  
  ##  X <- cbind(model.matrix(ff,data0)
  ##  nx <- ncol(X)
  ##  if (nx==0) stop("Zero design not allowed")
##################################################


  bidx0 <- seq(nx)
  midx0 <- bidx0; midx1 <- midx0+nx  
  dS0 <- rbind(rep(1,4),rep(1,4),rep(1,4)) ## MZ
  dS1 <- rbind(c(1,.5,.5,1),rep(1,4),c(1,.25,.25,1)) ## DZ

  ##mytr <- function(x) x; dmytr <- function(x) 1
  ##mytr <- function(x) x^2; dmytr <- function(x) 2*x
  ##mytr <- function(z) 1/(1+exp(-z)); dmytr <- function(z) exp(-z)/(1+exp(-z))^2
  mytr <- exp; dmytr <- exp; myinvtr <- log
  trname <- "exp"; invtrname <- "log"      
  ACDU <- sapply(c("a","c","d","e","u"),function(x) length(grep(x,tolower(type)))>0)
  if (ACDU["u"]) {
    ##      datanh <- function(r) 1/(1-r^2)
    dmytr <- function(z) 4*exp(2*z)/(exp(2*z)+1)^2
    mytr <- tanh;  myinvtr <- atanh
    trname <- "tanh"; invtrname <- "atanh"    
    dS0 <- rbind(c(0,1,1,0))
    vidx0 <- 1
    vidx1 <- 2
    dS1 <- dS0
    nvar <- length(vidx0)+length(vidx1)
  } else {
    nvar <- sum(ACDU[1:3])
    vidx0 <- vidx1 <- 1:nvar
    dS0 <- dS0[ACDU[1:3],,drop=FALSE]
    dS1 <- dS1[ACDU[1:3],,drop=FALSE]
  }  
  if (eqmean) {
    bidx1 <- bidx0
##    midx0 <- midx1 <- midx    
  } else {
    bidx1 <- bidx0+nx
##    midx0 <- 1:nx; midx1 <- midx0+nx
    nx <- 2*nx;
  }
  
  vidx0 <- vidx0+nx; vidx1 <- vidx1+nx
  vidx <- nx+seq_len(nvar)
  midx <- seq_len(nx)
  plen <- nx+nvar

  Am <- matrix(c(1,.5,.5,1),ncol=2)
  Dm <- matrix(c(1,.25,.25,1),ncol=2)
  Em <- diag(2)

##################################################

  Wide <- reshape(as.data.frame(Data),idvar=c(id,zyg),timevar=time,direction="wide")
  yidx <- paste(yvar,1:2,sep=".")
  rmidx <- c(id,yidx,zyg)
  ## if (!is.null(weight)) {
  ##   W <- cbind(data[,weight])
  ##   widx <- paste(weight,1:2,sep=".")
  ##   WW <- as.matrix(Wide[,widx])
  ##   
  ## }
  W0 <- W1 <- W2 <- NULL
  if (!is.null(weight)) {
    widx <- paste(weight,1:2,sep=".")
    rmidx <- c(rmidx,widx)
    W0 <- as.matrix(Wide[Wide[,zyg]==1,widx,drop=FALSE])
    W1 <- as.matrix(Wide[Wide[,zyg]==0,widx,drop=FALSE])
  }
  XX <- as.matrix(Wide[,setdiff(colnames(Wide),rmidx)])
  XX[is.na(XX)] <- 0
  Y0 <- as.matrix(Wide[Wide[,zyg]==1,yidx,drop=FALSE])
  Y1 <- as.matrix(Wide[Wide[,zyg]==0,yidx,drop=FALSE])
  XX0 <- XX[Wide[,zyg]==1,,drop=FALSE]
  XX1 <- XX[Wide[,zyg]==0,,drop=FALSE]
  
##################################################

###}}} setup

###{{{ Mean/Var function

  ##Marginals etc.
  MyData0 <- ExMarg(Y0,XX0,W0,dS0,eqmarg=TRUE,allmarg=allmarg)
  MyData1 <- ExMarg(Y1,XX1,W1,dS1,eqmarg=TRUE,allmarg=allmarg) 
  N <- cbind(sum(idx0),sum(idx1),sum(idx2)); 
  if (missing(OS)) N <- N[,-3,drop=FALSE]
  N <- cbind(N,
             2*nrow(MyData0$Y0)+nrow(MyData0$Y0_marg),
             2*nrow(MyData1$Y0)+nrow(MyData1$Y0_marg),
             nrow(MyData0$Y0),nrow(MyData1$Y0))

  if (samecens & !is.null(weight)) {
    MyData0$W0 <- cbind(apply(MyData0$W0,1,biweight))
    if (!is.null(MyData0$Y0_marg))
      MyData0$W0_marg <- cbind(apply(MyData0$W0_marg,1,biweight))
  }
  if (samecens & !is.null(weight)) {
    MyData1$W0 <- cbind(apply(MyData1$W0,1,biweight))
    if (!is.null(MyData1$Y0_marg))
      MyData1$W0_marg <- cbind(apply(MyData1$W0_marg,1,biweight))
  }
  rm(Y0,XX0,W0,Y1,XX1,W1)
  ##  suppressMessages(browser())
  ##  N <- rbind(##c("","","Complete","","Complete pairs",""),
  ##             rep(c("MZ","DZ"),ncol(N)/2), N)
  colnames(N) <- c("Total.MZ","Total.DZ","Complete.MZ","Complete.DZ","Complete pairs.MZ","Complete pairs.DZ")[seq(ncol(N))]
  rownames(N) <- rep("",nrow(N))
##  print(N,quote=FALSE)
  
  ## Mu <- function(p0) {
  ##   b0 <- cbind(p0[midx0])
  ##   b1 <- cbind(p0[midx1])
  ##   b00 <- b0; b11 <- b1
  ##   if (Bconstrain) {
  ##     b00 <- trMean(b0,Blen); b11 <- trMean(b1,Blen)
  ##   }
  ##     mu0 <- with(MyData0, X0%*%b00)
  ##     mu1 <- with(MyData1, X0%*%b11)
  ##   return(list(mu0=mu0,mu1=mu1))
  ## }
  Sigma <- function(p0) {    
    p0[vidx] <- mytr(p0[vidx])    
    if (ACDU["u"]) {     
      ##      Sigma0 <- Em+p0[plen-1]; Sigma1 <- Em+p0[plen]
      Sigma0 <- diag(2) + p0[plen-1]*matrix(c(0,1,1,0),2,2)
      Sigma1 <- diag(2) + p0[plen]*matrix(c(0,1,1,0),2,2)
    } else {    
      pv <- ACDU*1;  pv[which(ACDU[1:3])] <- p0[vidx]
      Sigma0 <- Em*pv["e"] + pv["a"] + pv["c"] + pv["d"]
      Sigma1 <- Em*pv["e"] + pv["a"]*Am + pv["c"] + pv["d"]*Dm
    }
    return(list(Sigma0=Sigma0,Sigma1=Sigma1))
  }

  ###}}} Mean/Var function
  
###{{{ U  

  U <- function(p,indiv=FALSE) {
    b0 <- cbind(p[bidx0])
    b1 <- cbind(p[bidx1])
    b00 <- b0; b11 <- b1
    ## if (Bconstrain) {
    ##   b00 <- trMean(b0,Blen); b11 <- trMean(b1,Blen)
    ## }
    if (bound) p[vidx] <- min(p[vidx],20)
    S <- Sigma(p)
    lambda <- eigen(S$Sigma0)$values
    if (any(lambda<1e-12 | lambda>1e9)) stop("Variance matrix out of bounds")
    ##    browser()
    ##mu0 <- with(MyData0, X0%*%b00)
    ##    Mu0 <- matrix(mu0,ncol=2,byrow=TRUE)
    ##        browser()
    
    Mu0 <- with(MyData0, cbind(XX0[,midx0,drop=FALSE]%*%b00,
                               XX0[,midx1,drop=FALSE]%*%b00))
    U0 <- with(MyData0, .Call("biprobit0",
                             Mu0,
                             S$Sigma0,dS0,Y0,XX0,W0,!is.null(W0),samecens))

    if (!is.null(MyData0$Y0_marg)) {
      mum <- with(MyData0, XX0_marg%*%b00)
       U_marg <- with(MyData0, .Call("uniprobit",
                                   mum,XX0_marg,
                                   S$Sigma0[1,1],t(dS0_marg),Y0_marg,
                                   W0_marg,!is.null(W0_marg),TRUE))
      U0$score <- rbind(U0$score,U_marg$score)
      U0$loglik <- c(U0$loglik,U_marg$loglik)
    }

##    mu1 <- with(MyData1, X0%*%b11) 
##    Mu1 <- matrix(mu1,ncol=2,byrow=TRUE)
    Mu1 <- with(MyData1, cbind(XX0[,midx0,drop=FALSE]%*%b11,
                               XX0[,midx1,drop=FALSE]%*%b11))

    U1 <- with(MyData1, .Call("biprobit0",
                             Mu1,
                             S$Sigma1,dS1,Y0,XX0,W0,!is.null(W0),samecens))
    if (!is.null(MyData1$Y0_marg)) {
      mum <- with(MyData1, XX0_marg%*%b11)
      U_marg <- with(MyData1, .Call("uniprobit",
                                    mum,XX0_marg,
                                    S$Sigma1[1,1],t(dS0_marg),Y0_marg,
                                    W0_marg,!is.null(W0_marg),TRUE))
      U1$score <- rbind(U1$score,U_marg$score)
      U1$loglik <- c(U1$loglik,U_marg$loglik)
    }

    if (indiv) {
      ll0 <- U0$loglik
      ll1 <- U1$loglik
      val0 <- U0$score[MyData0$id,,drop=FALSE]
      val1 <- U1$score[MyData1$id,,drop=FALSE]
      N0 <- length(MyData0$id)
      idxs0 <- seq_len(N0)
      if (length(MyData0$margidx)>0) {
        for (i in seq_len(N0)) {
          idx0 <- which((MyData0$idmarg)==(MyData0$id[i]))+N0
          idxs0 <- c(idxs0,idx0)
          val0[i,] <- val0[i,]+colSums(U0$score[idx0,,drop=FALSE])
        }
        val0 <- rbind(val0, U0$score[-idxs0,,drop=FALSE])
        ll0 <- c(ll0,ll0[-idxs0])
      }
      N1 <- length(MyData1$id)
      idxs1 <- seq_len(N1)
      if (length(MyData1$margidx)>0) {
        for (i in seq_len(N1)) {
          idx1 <- which((MyData1$idmarg)==(MyData1$id[i]))+N1
          idxs1 <- c(idxs1,idx1)
          val1[i,] <- val1[i,]+colSums(U1$score[idx1,,drop=FALSE])
        }
        val1 <- rbind(val1, U1$score[-idxs1,,drop=FALSE])
        ll1 <- c(ll1,ll1[-idxs1])        
      }
      val <- matrix(0,ncol=plen,nrow=nrow(val0)+nrow(val1))
      val[seq_len(nrow(val0)),c(bidx0,vidx0)] <- val0
      val[nrow(val0)+seq_len(nrow(val1)),c(bidx1,vidx1)] <- val1
      for (ii in vidx) {
        val[,ii] <- val[,ii]*dmytr(p[ii])
      }
      attributes(val)$logLik <- c(U0$loglik,U1$loglik)
      return(val)
      
#########      
      ## val <- matrix(0,ncol=plen,nrow=nrow(U0$score)+nrow(U1$score))
      ## val[seq_len(nrow(U0$score)),c(bidx0,vidx0)] <- U0$score
      ## val[nrow(U0$score)+seq_len(nrow(U1$score)),c(bidx1,vidx1)] <- U1$score
      ## for (ii in vidx) {
      ##   val[,ii] <- val[,ii]*dmytr(p[ii])
      ## }
      ## ## if (Bconstrain & Blen>0) {
      ## ##   Bidx <- attributes(b00)$idx
      ## ##   val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(b00)$D[Bidx,Bidx,drop=FALSE])
      ## ##   if (!eqmean) {
      ## ##     Bidx <- bidx1[attributes(b11)$idx]
      ## ##     val[,Bidx] <- as.numeric((val[,Bidx,drop=FALSE])%*%attributes(b11)$D[attributes(b11)$idx,attributes(b11)$idx])
      ## ##   }
      ## ## }      
      ## attributes(val)$logLik <- c(U0$loglik,U1$loglik)
      ## return(val)
      
    }
    val <- numeric(plen)
    val[c(bidx0,vidx0)] <- colSums(U0$score)
    val[c(bidx1,vidx1)] <- val[c(bidx1,vidx1)]+colSums(U1$score)
    for (ii in vidx)
      val[ii] <- val[ii]*dmytr(p[ii])
    ## if (Bconstrain & Blen>0) {
    ##   Bidx <- attributes(b00)$idx
    ##   val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b00)$D[Bidx,Bidx])
    ##   if (!eqmean) {
    ##     Bidx <- bidx1[attributes(b11)$idx]
    ##     val[Bidx] <- as.numeric(val[Bidx]%*%attributes(b11)$D[attributes(b11)$idx,attributes(b11)$idx])
    ##   }
    ## }    
    attributes(val)$logLik <- sum(U0$loglik)+sum(U1$loglik)
    return(val)
  }

###}}} U

###{{{ optim

  p0 <- rep(-1,plen); ##p0[vidx] <- 0
  if (type=="u")
    p0[vidx] <- 0.3
  if (!is.null(control$start)) {
    p0 <- control$start
    control$start <- NULL
  } else {
    X <- rbind(MyData0$XX0[,midx0,drop=FALSE],MyData0$XX0[,midx1,drop=FALSE])
    Y <- rbind(MyData0$Y0[,1,drop=FALSE],MyData0$Y0[,2,drop=FALSE])
    g <- suppressWarnings(glm(Y~-1+X,family=binomial(probit)))
    p0[midx] <- coef(g)
    ## if (Blen>0) {
    ##   pB <- p0[tail(midx,Blen)]
    ##   pB[1] <- ifelse(pB[1]<0,-2,log(pB[1]))
    ##   if (Blen>1) {
    ##     pB[seq_len(Blen-1)+1] <- -2
    ##   }
    ##   p0[tail(midx,Blen)] <- pB
    ## }
  }
 
  if (!missing(p)) return(U(p,indiv=indiv))


  f <- function(p) crossprod(U(p))[1]
  f0 <- function(p) -sum(attributes(U(p))$logLik)
  g0 <- function(p) -as.numeric(U(p))
  h0 <- function(p) crossprod(U(p,indiv=TRUE))

  
  if (!missing(constrain)) {
    freeidx <- is.na(constrain)
    f <- function(p) {      
      p1 <- constrain; p1[freeidx] <- p
      res <- U(p1)[freeidx]
      crossprod(res)[1]
    }
    f0 <- function(p) {
      p1 <- constrain; p1[freeidx] <- p
      -sum(attributes(U(p1))$logLik)
    }
    g0 <- function(p) {
      p1 <- constrain; p1[freeidx] <- p
      -as.numeric(U(p1)[freeidx])
    }
    p0 <- p0[is.na(constrain)]    
  }


  controlstd <- list(hessian=0)
  controlstd[names(control)] <- control
  control <- controlstd
  
  nlminbopt <- intersect(names(control),c("eval.max","iter.max","trace","abs.tol","rel.tol","x.tol","step.min"))
  ucminfopt <- intersect(names(control),c("trace","grtol","xtol","stepmax","maxeval","grad","gradstep","invhessian.lt"))
  optimopt <- names(control) 

  if (debug) browser()
  op <- switch(tolower(control$method),
               nlminb=nlminb(p0,f0,gradient=g0,control=control[nlminbopt]),
               optim=optim(p0,fn=f0,gr=g0,control=control[ucminfopt]),
               ucminf=,
               quasi=,
               gradient=ucminf(p0,fn=f0,gr=g0,control=control[ucminfopt],hessian=0),
               ## ,
               ## bhhh={
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
               ##   lava:::NR(start=p0,NULL,g0, h0,control=controlnr)
##               },
               ##                 op <- switch(mycontrol$method,
               ##                              ucminf=ucminf(p0,f,control=mycontrol[ucminfopt],hessian=F),
               ##                optim=optim(p0,f,control=mycontrol[ucminfopt],...),
                 nlminb(p0,f,control=control[nlminbopt]))
  ##  op <- nlm(f,p0,print.level=2)
  ##  op <- spg(p0,f,control=control)
  

  if (stderr) {
    
    ## WW <- rbind(W0,W1)
    ## ff <- function(p) {
    ##   UU <- U(p,indiv=TRUE)
    ##   s <- apply(UU,2,function(x) x) ##*WW[,1,drop=FALSE])
    ##   l <- attributes(UU)$logLik##*WW[,1]
    ##   res <- structure(sum(l), grad=colSums(s))
    ##   res
    ## }
    
    UU <- U(op$par,indiv=TRUE)
    I <- -numDeriv::jacobian(U,op$par)
    tol <- 1e-15
    V <- Inverse(I,tol)
    sqrteig <- attributes(V)$sqrteig
    J <- NULL
    if (robustvar) {
      J <- crossprod(UU)
      V <- V%*%J%*%V
    }
    if (any(sqrteig<tol)) warning("Near-singular covariance matrix (pseudo-inverse used)")
  } else {
    UU <- matrix(NA,ncol=length(op$par),nrow=1)
    I <- J <- V <- matrix(NA,ncol=length(op$par),nrow=length(op$par))
  }

  ###}}} optim

###{{{ return

  cc <- cbind(op$par,sqrt(diag(V)))
  cc <- cbind(cc,cc[,1]/cc[,2],2*(1-pnorm(abs(cc[,1]/cc[,2]))))
  colnames(cc) <- c("Estimate","Std.Err","Z","p-value")
  vnames1 <- NULL
  trnam <- " "
  if (debug) browser()
  if (!eqmean) {
    rnames1 <- c(paste(rnames1,"MZ",sep=trnam),paste(rnames1,"DZ",sep=trnam))
  }

  if (ACDU["u"]) {
##    rnames <- c(rnames1,paste(c("log(var(U))","log(var(U))"),c("MZ","DZ"),sep=trnam))
    rnames <- c(rnames1,paste(c("atanh(rho)","atanh(rho)"),c("MZ","DZ"),sep=trnam))
  } else {
    rnames <- c(rnames1,c("log(var(A))","log(var(C))","log(var(D))")[ACDU[1:3]])
  }
  if (!missing(constrain)) rnames <- rnames[freeidx]
  rownames(cc) <- rnames
  rownames(V) <- colnames(V) <- rnames
  S <- Sigma(op$par)

  npar <- list(intercept=attributes(terms(formula))$intercept,
               pred=nrow(attributes(terms(formula))$factor)-1,
               var=sum(ACDU[-4]),
               ACDU=ACDU[-4]*1)
  
  npar[unlist(lapply(npar,length))==0] <- 0
##  npar$var <- nrow(cc)-sum(unlist(npar))
  
  val <- list(coef=cc,vcov=V,score=UU,logLik=attributes(UU)$logLik,opt=op, Sigma0=S$Sigma0, Sigma1=S$Sigma1, dS0=dS0, dS1=dS1, N=N, midx0=midx0, midx1=midx1, vidx0=vidx0, vidx1=vidx1, eqmean=eqmean, I=I,J=J, robustvar=robustvar,
              transform=list(tr=mytr, invtr=myinvtr, dtr=dmytr,
                name=trname, invname=invtrname),
              SigmaFun=Sigma, ##MuFun=Mu,
              npar=npar
              )
  class(val) <- c("bptwin","biprobit")
  return(val)
}

###}}} return


