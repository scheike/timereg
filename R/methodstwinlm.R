###{{{ print.twinlm

##' @S3method print twinlm
print.twinlm <- function(x,...) {
  print(summary(x,...))
  invisible(x)
}

###}}} print.twinlm

###{{{ summary.twinlm

##' @S3method summary twinlm
summary.twinlm <- function(object,...) {
  e <- object$estimate
  zygtab <- with(object, table(data[,zyg]))
  theta <- pars(e)
  theta.sd <- sqrt(diag(e$vcov))
  myest <- cbind(theta,theta.sd,(Z <- theta/theta.sd),2*(1-pnorm(abs(Z))))
  colnames(myest) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")
  
  if (length(grep("u",object$type))>0) {

    cc <- coef(e)
    pi <- seq_len(nrow(myest))
    pp <- modelPar(e$model,pi)$p
    nn <- rep(NA,nrow(myest))
    nn[pp[[1]]] <- coef(e$model$lvm[[1]])
    nn.na <- which(is.na(nn))
    i2 <- 0
    for (i in 1:length(pp)) {
      if (nn.na[1]%in%pp[[i]]) {
        i2 <- i
        break;
      }      
    }
    pp.i2 <- which(pp[[i2]]%in%nn.na)
    parnum <- pp[[i2]][pp.i2]
    nn[nn.na] <- coef(e$model$lvm[[i2]])[pp.i2]
    nn <- gsub(".1","",nn,fixed=TRUE)
    nn <- gsub(".2","",nn,fixed=TRUE)
    u.idx <- c(grep("<-u",nn))
    e.idx <- c(grep("<-e",nn))

    lastpos <- all(parnum>=u.idx) ## i2-coef positioined after coef. of model 1
    isMZ <- "sdu1"%in%parlabels(e$model$lvm[[i2]])
    if ((lastpos & isMZ)|(!lastpos & !isMZ)) {
      u.idx <- rev(u.idx); e.idx <- rev(e.idx)
    }
    nn <- c(nn); nn[u.idx] <- c("MZ:sd(U)","DZ:sd(U)")
    if (length(e.idx)==1) {
      nn[e.idx] <- "MZ:sd(E)"
    } else {
      nn[e.idx] <- c("MZ:sd(E)","DZ:sd(E)")
    }    
    rownames(myest) <- nn
    neword <- c(setdiff(seq_len(nrow(myest)),c(e.idx,u.idx)),e.idx,u.idx)

    MZcc <- DZcc <- NULL
##     if (object$binary) {
##       MZest <- myest[-match("DZ:sd(U)",rownames(myest)),1]
##       DZest <- myest[-match(c("MZ:sd(E)","MZ:sd(U)"),rownames(myest)),1]
##       M1 <- moments(object$model.mz,MZest,cond=FALSE)
##       M2 <- moments(object$model.dz,DZest,cond=FALSE)
##       myidx <- index(object$model.mz)$endo.obsidx
##       MZcc <- with(M1, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
##       myidx <- index(object$model.dz)$endo.obsidx
##       DZcc <- with(M2, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
## ##      browser()
##     }
    
    K <- 0 ##ifelse (object$binary,object$probitscale,0)
    ##L <- binomial("logit")
    logit <- function(p) log(p/(1-p))
    tigol <- function(z) 1/(1+exp(-z))
    if (length(e.idx)==1) {
      corDZ <- function(x) (x[2]^2)/(x[2]^2+K)
      corMZ <- function(x) (x[1]^2)/(x[1]^2+x[3]^2)
    } else {
      corDZ <- function(x) (x[2]^2)/(x[2]^2+x[4]^2)
      corMZ <- function(x) (x[1]^2)/(x[1]^2+x[3]^2)
    }
    h <- function(x) 2*(corMZ(x)-corDZ(x))
    
    dh <- function(x) {
      x2 <- x^2;
      s2 <- ifelse(length(e.idx)==1,x2[2]+K^2,x2[2]+x2[4])
      s1 <- ifelse(length(e.idx)==1,x2[1]+x2[3],x2[1]+x2[3])
      if (length(e.idx)==1) {
        De <- -4*x[3]*x2[1]/s1^2
      } else {
        De <- c(-4*x[3]*x2[1]/s1^2,4*x[4]*x2[2]/s2^2)
      }      
      c(2*c(1/s1^2,-1/s2^2)*(2*x[1:2]*c(s1,s2)-x2[1:2]*2*x[1:2]),De)      
    }

    
    dlogit <- function(p) 1/(p*(1-p))
    logith <- function(x) logit(h(x))
    dlogith <- function(x) dlogit(h(x))*dh(x)

    V <- e$vcov[c(u.idx,e.idx),c(u.idx,e.idx)]
    b <- pars(e)[c(u.idx,e.idx)]
    Db <- dlogith(b)
    sd. <- (t(Db)%*%V%*%Db)^0.5
    hval <- c(h(b),(t(dh(b))%*%V%*%(dh(b)))^0.5)
    hci <- tigol(logith(b)+qnorm(0.975)*c(-1,1)*sd.)
    names(hci) <- c("2.5%","97.5%")
    res <- list(estimate=myest[neword,], zyg=zygtab,
                varEst=NULL, varSigma=NULL, heritability=hval, hci=hci,
                corMZ=corMZ(b), corDZ=corDZ(b),
                concMZ=MZcc, concDZ=DZcc)                
    class(res) <- "summary.twinlm"
    return(res)
  }

  rownames(myest) <- gsub(".1","",coef(Model(Model(e))[[1]],
                                       mean=e$meanstructure, silent=TRUE),
                          fixed=TRUE)
  rownames(myest) <- gsub(".2","",rownames(myest),fixed=TRUE)
##  rownames(myest) <- coef(Model(Model(e))[[1]],
##                                      mean=e$meanstructure, silent=TRUE)
                     
  lambda.idx <- sapply(c("<-a1","<-c1","<-d1","<-e1"),function(x) grep(x,rownames(myest)))
  lambda.w <- which(sapply(lambda.idx, function(x) length(x)>0))
  rownames(myest)[unlist(lambda.idx)] <- paste("sd(",c("A)","C)","D)","E)"),sep="")[lambda.w]


  varEst <- rep(0,4)
  varEst[lambda.w] <- myest[unlist(lambda.idx),1]
  varSigma <- matrix(0,4,4);
  varSigma[lambda.w,lambda.w] <- e$vcov[unlist(lambda.idx),unlist(lambda.idx)]

  L <- binomial("logit")
  varcomp <- c()
  genpos <- c()
  pos <- 0
  ## if ("e1"%in%latent(e) & object$binary) {
  ##   if (length(object$probitscale)>0)
  ##     varEst[4] <- object$probitscale
  ##   if ("a1"%in%latent(e)) { varcomp <- "lambda[a]"; pos <- pos+1; genpos <- c(genpos,pos) }
  ##   if ("c1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[c]"); pos <- pos+1 }
  ##   if ("d1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[d]"); pos <- pos+1;
  ##                            genpos <- c(genpos,pos) }
  ##   f <- paste("h2~",paste(varcomp,collapse="+"))
  ##   constrain(e, as.formula(f)) <- function(x) L$linkfun(sum(x[genpos]^2)/sum(c(x^2,varEst[4])))
    
  ## } else
  {
    varcomp <- c()
        if ("a1"%in%latent(e)) { varcomp <- "lambda[a]"; pos <- pos+1;
                                 genpos <- c(genpos,pos) }
    if ("c1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[c]"); pos <- pos+1 }
    if ("d1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[d]"); pos <- pos+1;
                             genpos <- c(genpos,pos) }
    if ("e1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[e]"); pos <- pos+1 }
    f <- paste("h2~",paste(varcomp,collapse="+"))

    constrain(e, as.formula(f)) <- function(x) L$linkfun(sum(x[genpos]^2)/sum(x^2))
  }
  ci.logit <- L$linkinv(constraints(e)["h2",5:6])
  
  h <- function(x) (x[1]^2)/(sum(x^2))
  dh <- function(x) {
    y <- x^2
    cbind(1/sum(y)^2*c(2*x[1]*(sum(y)-y[1]), -2*x[2:4]*y[1]))
  }
  ##  f(x1,x2,x3,x4) = x1/(x1+x2+x3+x4) = x1/s
  ## Quotient rule: (u/v)' = (u'v - uv')/v^2
  ##  f1(x1,x2,x3,x4) = (s - x1)/s^2 = (x2+x3+x4)/s^2
  ##  f2(x1,x2,x3,x4) = -1/(s)^2
  h2 <- function(x) (x[1]^2+x[3]^2)/sum(x^2)
  dh2 <- function(x) {
    y <- x^2
    cbind(1/sum(y)^2*c(
                       2*x[1]*(sum(y)-(y[1]+y[3])),
                       -2*x[2]*(y[1]+y[3]),
                       2*x[3]*(sum(y)-(y[1]+y[3])),
                       -2*x[4]*(y[1]+y[3])
                       ))
  }
  hval <- cbind(h(varEst), (t(dh(varEst))%*%varSigma%*%(dh(varEst)))^0.5)
  colnames(hval) <- c("Estimate", "Std.Err"); rownames(hval) <- "h squared"
  h2val <- cbind(h2(varEst), (t(dh2(varEst))%*%varSigma%*%(dh2(varEst)))^0.5)
  colnames(h2val) <- c("Estimate", "Std.Err"); rownames(h2val) <- "h squared"

  corMZ <- sum(varEst[1:3]^2)/sum(varEst^2)
  corDZ <- sum(varEst[1:3]^2*c(0.5,1,0.25))/sum(varEst^2)

  MZcc <- DZcc <- NULL
##   if (object$binary) {    
##     M1 <- moments(object$model.mz,coef(object),cond=FALSE)
##     M2 <- moments(object$model.dz,coef(object),cond=FALSE)
##     myidx <- index(object$model.mz)$endo.obsidx
##     MZcc <- with(M1, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
##     DZcc <- with(M2, pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=xi[myidx,1],sigma=C[myidx,myidx]))
## ##    browser()
##   }
  
  res <- list(estimate=myest, zyg=zygtab, varEst=varEst, varSigma=varSigma, hval=hval, heritability=h2val, hci=ci.logit, corMZ=corMZ, corDZ=corDZ, acde=acde.twinlm(object),
              logLik=logLik(e), AIC=AIC(e), BIC=BIC(e),
              concMZ=MZcc, concDZ=DZcc)                              
  class(res) <- "summary.twinlm"
  return(res)
}

###}}} summary.twinlm

###{{{ print.summary.twinlm

##' @S3method print summary.twinlm
print.summary.twinlm <- function(x,signif.stars=FALSE,...) {
  printCoefmat(x$estimate,signif.stars=signif.stars,...)
  cat("\n")
  myzyg <- with(x,zyg)
  names(myzyg) <- c("Group 1 (DZ)", "Group 2 (MZ)")
  print(myzyg)
##   cat("\n")
##   mynames <- c("sigmaA","sigmaC","sigmaD","sigmaE")
##   for (i in 1:4) {
##     cat(mynames[i], "=", x$varEst[i], "\n")
##   }  
  cat("\nVariance decomposition:\n")
  print(x$acde)
  cat("\n")
  ##  cat("hn2 = ", hn2, "\t hb2 = ", hb2, "\n\n")
##  cat("Narrow-sense heritability (additive genetic factors):\n")
##  print(x$hval)
##  cat("\n")  
  cat("Broad-sense heritability (total genetic factors):\n")
  h <- with(x, c(heritability[1],hci));
  names(h) <- c("Estimate",names(x$hci))
  h <- na.omit(h)
  print(h)  
  cat("\n")  
  cat("Correlation within MZ:", x$corMZ, "\n")
  cat("Correlation within DZ:", x$corDZ, "\n")
  cat("\n")
  if (!is.null(x$concMZ)) {
    cat("Concordance MZ:", x$concMZ, "\n")
    cat("Concordance DZ:", x$concDZ, "\n")
    cat("\n")
  }
  print(x$logLik)
  cat("AIC:", x$AIC, "\n")
  cat("BIC:", x$BIC, "\n")
}

###}}} print.summary.twinlm

###{{{ compare.twinlm

##' @S3method compare twinlm
compare.twinlm <- function(object,...) {
  objects <- list(object,...)
  if (length(objects)<2)
    return(summary(objects))
  res <- list()
  for (i in 1:(length(objects)-1)) {
    res <- c(res, list(compare(objects[[i]]$estimate,objects[[i+1]]$estimate)))
  }
  if (length(res)==1)
    return(res[[1]])
  return(res)
}
###}}} compare.twinlm

###{{{ plot.twinlm

##' @S3method plot twinlm
plot.twinlm <- function(x,diag=TRUE,labels=TRUE,...) {
  op <- par(mfrow=c(2,1))
  plot(x$model,...)
  par(op)
}
###}}}

###{{{ vcov.twinlm

##' @S3method vcov twinlm
vcov.twinlm <- function(object,...) {
  return(object$vcov)
}
###}}} vcov.twinlm

###{{{ logLik.twinlm

##' @S3method logLik twinlm
logLik.twinlm <- function(object,...) logLik(object$estimate,...)
###}}} logLik.twinlm

###{{{ model.frame.twinlm

##' @S3method model.frame twinlm
model.frame.twinlm <- function(formula,...) {
  return(formula$estimate$model$data)
}
###}}} model.frame.twinlm

###{{{ acde

##"acde" <- function(x,...) UseMethod("acde")
acde.twinlm <- function(x,...) {
  m <- x$estimate$model$lvm[[1]]
  lambdas <- c("lambda[a]","lambda[c]","lambda[d]","lambda[e]")
  ACDE <- lambdas%in%as.vector(m$par)
  lcur <- lambdas[ACDE]
  for (l in lcur) {
    pos <- which(lcur%in%l)
    par <- substr(strsplit(l,"[",fixed=TRUE)[[1]][2],1,1)
    f <- as.formula(paste(par,"~",paste(lcur,collapse="+")))
    myfun <- eval(parse(text=paste("function(x) qnorm(x[",pos,"]^2/sum(x^2))")))
    constrain(x$estimate,f) <- myfun ##function(x) x[get("pos")]^2/sum(x^2)
  }
  M <- pnorm(constraints(x$estimate)[,c(1,5,6),drop=FALSE])
  rownames(M) <- toupper(rownames(M))
  M
}

###}}} acde
