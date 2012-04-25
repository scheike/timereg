##' @S3method summary bptwin
summary.bptwin <- function(object,level=0.05,...) {
  logit <- function(p) log(p/(1-p))
  tigol <- function(z) 1/(1+exp(-z))
  dlogit <- function(p) 1/(p*(1-p))
  trnam <- " "
  vcoef1 <- paste("log(var(",c("A","C","D"),"))",sep="")
  vcoef2 <- paste("atanh(",
                  c(paste("rho)","MZ",sep=trnam),
                    paste("rho)","DZ",sep=trnam)),sep="")
  idx1 <- na.omit(match(vcoef1,names(coef(object))))
  idx2 <- na.omit(match(vcoef2,names(coef(object))))
  CIs <- c()
  alpha <- level/2
  CIlab <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  V <- c()
  if (length(idx2)>0) {    
    idx <- idx2
    V <- vcov(object)[idx,idx]
    arho <- coef(object)[idx2[1:2]]
    mz <- multinomlogit(coef(object)[idx2[1]]); names(mz) <- c("U","E")
    dz <- multinomlogit(coef(object)[idx2[2]]); names(dz) <- c("U","E")
    cc <- tanh(arho)
  ##    cc <- c(mz[1],dz[1]) ##,mz[2],dz[2])
    names(cc) <- c("Correlation MZ","Correlation DZ")
##    corMZ <- mz[1]; corDZ <- dz[1]
    corMZ <- cc[1]; corDZ <- cc[2]
##    D <- (cbind(c(attributes(mz)$gradient[1],0),c(0,attributes(dz)$gradient[1])))
    D <- diag(object$tr$dtr(arho))
    h <- function(x) 2*(x[1]-x[2])
    dh <- function(x) c(2,-2)
    i1 <- 1:2
    corr <- NULL
  }
  if (length(idx1)>0) {
    idx <- idx1
    V <- vcov(object)[idx,idx]
    ACD <- match(names(coef(object))[idx1],vcoef1)
    nn <- c(c("A","C","D")[ACD],"E")
    dzsc <- c(1/2,1,1/4)[ACD]
    pp <- coef(object)[idx1]
    cc <- multinomlogit(pp); names(cc) <- nn
    D <- attributes(cc)$gradient  
    ##    p <- coef(object)[2:3,1]
    ##    F <- function(p) {
    ##      logit(multinomlogit(p))
    ##    }    
    cc2 <- logit(cc)
    D2 <- diag(dlogit(cc))
    DD <- D2%*%D
    Vc2 <- DD%*%V%*%t(DD)
    CIs <- tigol(cc2%x%cbind(1,1)+diag(Vc2)^0.5%x%cbind(-1,1)*qnorm(1-alpha))
    K <- length(ACD)
    Ki <- seq_len(K)
    corMZ <- sum(cc[Ki]); corDZ <- sum(cc[Ki]*dzsc)
    i1 <- seq_len(length(dzsc))
    h <- function(x) 2*(sum(x[i1])-sum(x[i1]*dzsc))
    dh <- function(x) 2*(1-dzsc)
    
  }
  Vc <- D%*%V%*%t(D)
  datanh <- function(r) 1/(1-r^2)
  if (length(idx1)>0) {
    pp <- coef(object)[idx]
    b <- cbind(rep(1,K))
    corMZ.sd <- (t(b)%*%Vc[Ki,Ki]%*%b)[1]^0.5
    corDZ.sd <- (t(dzsc)%*%Vc[Ki,Ki]%*%dzsc)[1]^0.5    
    corr <- rbind(c(corMZ,corMZ.sd),c(corDZ,corDZ.sd))
    zrho <- atanh(corr[,1])
    zrho.var <- datanh(corr[,1])^2*corr[,2]^2
    corr <- cbind(corr, tanh(zrho%x%cbind(1,1)+zrho.var^0.5%x%cbind(-1,1)*qnorm(1-alpha)))
    rownames(corr) <- c("Correlation MZ","Correlation DZ")    
  } else {   
    zrho <- atanh(cc)
    zrho.var <- datanh(cc)^2*diag(Vc)
    CIs <- tanh(zrho%x%cbind(1,1)+zrho.var^0.5%x%cbind(-1,1)*qnorm(1-alpha))
  }
  newcoef <- rbind(cbind(cc,diag(Vc)^0.5,CIs),corr);
  ##  CIs <- rbind(CIs,c(NA,NA),c(NA,NA))
  ##  newcoef <- cbind(newcoef,CIs)
  colnames(newcoef) <- c("Estimate","Std.Err",CIlab)
  logith <- function(x) logit(h(x))
  dlogith <- function(x) dlogit(h(x))*dh(x)
  Dlh <- dlogith(cc[i1])
  sdlh <- (t(Dlh)%*%Vc[i1,i1]%*%(Dlh))[1]^0.5
  H <- h(cc[i1])
  hstd <- t(dh(cc[i1]))%*%Vc[i1,i1]%*%dh(cc[i1])
  ci <- tigol(logith(cc[i1]) + qnorm(1-alpha)*c(-1,1)*sdlh)  
  
  concordance <-  conditional <- marg <- c()

  probs <- function(p,idx=1) {
    if (idx==0) {
      m <- p[1]
      ##else m <- p[length(object$midx0)+1]
      S <- object$SigmaFun(p)
      conc1 <- pmvnorm(upper=c(m,m),sigma=S[[1]],algorithm="Miwa")
      conc2 <- pmvnorm(upper=c(m,m),sigma=S[[2]],algorithm="Miwa")
      marg <- pnorm(m,sd=S[[1]][1,1]^0.5)      
      return(logit((conc1-conc2)/(marg*(1-marg))))
    }
    S <- (object$SigmaFun(p))[[idx]]
    m <- 0
    if((object$npar$intercept==1 & idx==1) | object$eqmean) m <- p[1]
    else m <- p[length(object$midx0)+1]
    mu.cond <- function(x) m+S[1,2]/S[2,2]*(x-m)
    var.cond <- S[1,1]-S[1,2]^2/S[2,2]    
    conc <- pmvnorm(upper=c(m,m),sigma=S,algorithm="Miwa")
    marg <- pnorm(m,sd=S[1,1]^0.5)
    cond <- conc/marg
    logit(c(conc,cond,marg))
  }
   
  mycoef <- coef(object)
  formals(probs) <- alist(p=,idx=0)
  hp <- probs(mycoef)
  Dhp <- grad(probs,mycoef)
  shp <- diag(t(Dhp)%*%vcov(object)%*%(Dhp))^0.5
  
  formals(probs) <- alist(p=,idx=1)
  probMZ <- probs(mycoef)
  
  Dp0 <- jacobian(probs,mycoef)
  formals(probs) <- alist(p=,idx=2)
  probDZ <- probs(mycoef)
  Dp1 <- jacobian(probs,mycoef)
  sprobMZ <- diag((Dp0)%*%vcov(object)%*%t(Dp0))^0.5
  sprobDZ <- diag((Dp1)%*%vcov(object)%*%t(Dp1))^0.5
  probMZ <- tigol(cbind(probMZ,probMZ-qnorm(1-alpha)*sprobMZ,probMZ+qnorm(1-alpha)*sprobMZ))
  probDZ <- tigol(cbind(probDZ,probDZ-qnorm(1-alpha)*sprobDZ,probDZ+qnorm(1-alpha)*sprobDZ))
  rownames(probMZ) <- rownames(probDZ) <- c("Concordance","Conditional","Marginal")
  colnames(probMZ) <- colnames(probDZ) <- c("Estimate",CIlab)
 
  ## mu <- coef(object)[c(object$bidx0[1],object$bidx1[1])]
  ## Sigma <- list(object$Sigma0,object$Sigma1)
  ## for (i in 1:2) {
  ##   conc <- function()
  ##   mu.cond <- function(x) mu+Sigma[[i]][1,2]/Sigma[[i]][2,2]*(x-mu[i])
  ##   var.cond <- Sigma[[i]][1,1]-Sigma[[i]][1,2]^2/Sigma[[i]][2,2]    
  ##   cc0 <- pmvnorm(upper=c(mu[i],mu[i]),sigma=Sigma[[i]])
  ##   px <- pnorm(mu[i],sd=Sigma[[i]][2,2]^0.5)
  ##   concordance <- c(concordance,cc0)
  ##   marg <- c(marg,px)
  ##   conditional <- c(conditional,cc0/px)
  ## }
  ## names(concordance) <- names(conditional) <- c("MZ","DZ")

  hval <- rbind(c(H,hstd^0.5,ci)); colnames(hval) <- c("Estimate","Std.Err",CIlab); 

  hval <- rbind(hval, tigol(c(hp,NA,hp-qnorm(1-alpha)*shp,hp+qnorm(1-alpha)*shp)))
  rownames(hval) <- c("Heritability","h^2 pr")
  
  Nstr <- object$N
  nN <- ncol(object$N)
  npos <- seq(nN/2)
  Nstr <- rbind(paste(Nstr[npos*2-1],Nstr[npos*2],sep="/"))
  rownames(Nstr) <- ""
  colnames(Nstr) <- unlist(lapply(strsplit(colnames(object$N)[npos*2-1],".",fixed=TRUE),
                                  function(x) paste(x[1], "MZ/DZ")))
  res <- list(object=object, h=hval,
              probMZ=probMZ, probDZ=probDZ, Nstr=Nstr,
              coef=newcoef) ##, concordance=concordance, conditional=conditional)
  class(res) <- "summary.bptwin"
  res
}

print.summary.bptwin <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  print(x$object,digits=digits,...)
  cat("\n")
  x$Nstr <- x$Nstr[,which((colnames(x$Nstr)!="Complete MZ/DZ")),drop=FALSE]
  print(x$Nstr,quote=FALSE)
  cat("\n")
  print(RoundMat(x$coef[,-2,drop=FALSE],digits=digits),quote=FALSE)
  cat("\nMZ:\n");
  print(RoundMat(x$probMZ,digits=digits),quote=FALSE)
  cat("DZ:\n")
  print(RoundMat(x$probDZ,digits=digits),quote=FALSE)
##  cat("\nConcordance (MZ; DZ):\t\t", x$concordance,"\n")
##  cat("Case-wise concordance (MZ; DZ):\t", x$conditional,"\n\n")
  cat("\n")
  print(RoundMat(x$h[,-2,drop=FALSE],digits=digits),quote=FALSE)
  cat("\n")
}
