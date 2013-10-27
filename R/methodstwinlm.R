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
    browser()
  e <- object$estimate    
  ## counts <- function(dd) {
  ##   dd0 <- apply(dd,2,function(x) !is.na(x))
  ##   pairs <- sum(dd0[,1]*dd0[,2])
  ##   singletons <- sum((!dd0[,1])*dd0[,2] + (!dd0[,2])*dd0[,1])
  ##   return(c(pairs,singletons))
  ## }
  ## mz  <- counts(object$data.mz[,object$outcomes])
  ## dz  <- counts(object$data.dz[,object$outcomes])
  ## if (!object$estimate$model$missing) {
  ##   zygtab <- c("MZ-pairs"=mz[1],"DZ-pairs"=dz[1])
  ## } else {
  ##   zygtab <- c(paste(mz,collapse="/"),paste(dz,collapse="/"))
  ##   names(zygtab) <- c("MZ-pairs/singletons","DZ-pairs/singletons")
  ## }
    zygtab <- NULL
  ## if (object$OS) {
  ##   os  <- counts(object$data.os[,object$outcomes])
  ##   if (!object$estimate$model$missing) {
  ##     zygtab <- c(zygtab,"OS-pairs"=os[1])
  ##   } else {
  ##     zygtab <- c(zygtab,paste(os,collapse="/"))
  ##     names(zygtab)[3] <- "OS-pairs/singletons"
  ##   }
  ## } 


  ## zygtab <- with(object, table(data[,zyg]))
  ## names(zygtab) <- paste(names(zygtab),"pairs",sep="-")

  theta <- pars(e)
  theta.sd <- sqrt(diag(e$vcov))
  myest <- cbind(theta,theta.sd,(Z <- theta/theta.sd),2*(pnorm(abs(Z),lower.tail=FALSE)))
  colnames(myest) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)")

  if (object$type%in%c("u","flex","sat")) {
    corMZ <- corDZ <- NULL
    if (object$constrain) {
      i1 <- lava:::parpos.multigroup(object$estimate$model,p="atanh(rhoMZ)")[1]
      i2 <- lava:::parpos.multigroup(object$estimate$model,p="atanh(rhoDZ)")[1]
      ## i3 <- lava:::parpos.multigroup(object$estimate$model,p="atanh(rhoOS)")[1]
      if (length(i1)>0) {
        ## corest <- coef(object$estimate,level=0)[c(i1,i2,i3)]
        ## sdest <- vcov(object$estimate)[cbind(c(i1,i2,i3),c(i1,i2,i3))]^0.5
        corest <- coef(object$estimate,level=0)[c(i1,i2)]
        sdest <- vcov(object$estimate)[cbind(c(i1,i2),c(i1,i2))]^0.5

        ciest <- tanh(cbind(corest,corest)+qnorm(0.975)*cbind(-sdest,sdest))
        corest <- tanh(corest)
        corMZ <- c(corest[1],ciest[1,])
        corDZ <- c(corest[2],ciest[2,])
        corOS <- NULL
        ## if (object$OS) {
        ##   corOS <- c(corest[3],ciest[3,])
        ## }
      }
    }
    aa <- capture.output(e)
    res <- list(estimate=aa, zyg=zygtab,
                varEst=NULL, varSigma=NULL, heritability=NULL,
                corMZ=corMZ, corDZ=corDZ, 
                logLik=logLik(e), AIC=AIC(e), BIC=BIC(e), type=object$type
                )                
    class(res) <- "summary.twinlm"
    return(res)
  }


  ## kinshipOS  <- NULL
  ## if(object$OS) {
  ##   KinshipOS <- constraints(e,k=3)[,c(1,5,6),drop=FALSE]
  ##   rownames(KinshipOS)[1] <- "OS Kinship"
  ##   if (nrow(KinshipOS)>1) {
  ##     rownames(KinshipOS) <- c("OS Kinship(A)","OS Kinship(D)")
  ##   }
  ##   myest <- myest[-lava:::parpos.multigroup(object$estimate$model,p=c("ra","rd")),,drop=FALSE]
  ## }

  r1 <- gsub(".1","",coef(Model(Model(e))[[1]],
                          mean=e$meanstructure, silent=TRUE),
                          fixed=TRUE)
  rownames(myest) <- seq(nrow(myest))
  idx <- seq(nrow(myest)); 
  rownames(myest)[idx] <- r1

  coefpost <- paste(lava.options()$symbol[1],c("a1","c1","d1","e1"),sep="")
  coefpost2 <- paste(lava.options()$symbol[1],c("a2","c2","d2","e2"),sep="")
  myest.varpos <- unlist(sapply(coefpost,function(x) grep(x,rownames(myest))))  
  lambda.idx <- sapply(coefpost,function(x) grep(x,names(coef(e))))
  lambda.idx2 <- sapply(coefpost2,function(x) grep(x,names(coef(e))))
  for (k in seq_len(length(lambda.idx)))
    if (length(lambda.idx[[k]])==0) lambda.idx[[k]] <- lambda.idx2[[k]] 

  lambda.w <- which(sapply(lambda.idx, function(x) length(x)>0))

  rownames(myest)[myest.varpos] <- paste("sd(",c("A)","C)","D)","E)"),sep="")[lambda.w]

  varEst <- rep(0,4)
  varEst[lambda.w] <- myest[myest.varpos,1]
  varSigma <- matrix(0,4,4);
  varSigma[lambda.w,lambda.w] <- e$vcov[unlist(lambda.idx),unlist(lambda.idx)]

  L <- binomial("logit")
  varcomp <- c()
  genpos <- c()
  pos <- 0L
  {
    varcomp <- c()
    if ("a1"%in%latent(e)) { varcomp <- "lambda[a]"; pos <- pos+1
                             genpos <- c(genpos,pos) }
    if ("c1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[c]"); pos <- pos+1 }
    if ("d1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[d]"); pos <- pos+1;
                             genpos <- c(genpos,pos) }
    if ("e1"%in%latent(e)) { varcomp <- c(varcomp,"lambda[e]"); pos <- pos+1 }
    f <- paste("h2~",paste(varcomp,collapse="+"))

    constrain(e, as.formula(f)) <- function(x) {
      L$linkfun(sum(x[genpos]^2)/sum(x^2))
    }

  }
  ci.logit <- L$linkinv(constraints(e,k=1)["h2",5:6])
  
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
 
  atanhcorMZf <- function(x) atanh(sum(x[1:3]^2)/sum(x^2))
  atanhcorDZf <- function(x) atanh(sum(x[1:3]^2*c(0.5,1,0.25))/sum(x^2))
 
  e1 <- atanhcorMZf(varEst)
  D1 <- numDeriv::grad(atanhcorMZf,varEst)
  s1 <- (t(D1)%*%varSigma%*%(D1))^0.5
  ci1 <- e1+qnorm(0.975)*c(-1,1)*s1
  e2 <- atanhcorDZf(varEst)
  D2 <- numDeriv::grad(atanhcorDZf,varEst)
  s2 <- (t(D2)%*%varSigma%*%(D2))^0.5
  ci2 <- e2+qnorm(0.975)*c(-1,1)*s2
  
  corMZ <- c(tanh(c(e1,ci1)))
  corDZ <- c(tanh(c(e2,ci2)))  

  acde <- acde.twinlm(object)
  coef <- rbind(acde)  

  ## corOS <- NULL
  ## if (object$OS)
  ##   acde <- rbind(acde,KinshipOS)

  
  hrow <- rbind(c(h2val,ci.logit));
  rownames(hrow) <- "Broad-sense heritability"
  colnames(hrow)[1:2] <- c("Estimate","Std.Err")
  
  all <- rbind(hrow[,c(1,3,4),drop=FALSE],coef,corMZ,corDZ) ##,corOS)
  
  res <- list(estimate=myest, zyg=zygtab, varEst=varEst,
              varSigma=varSigma, heritability=hrow, corMZ=corMZ, corDZ=corDZ,
              acde=acde, logLik=logLik(e), AIC=AIC(e), BIC=BIC(e),
              type=object$type, coef=coef, all=all);
  class(res) <- "summary.twinlm"
  return(res)
}

###}}} summary.twinlm

###{{{ print.summary.twinlm

##' @S3method print summary.twinlm
print.summary.twinlm <- function(x,signif.stars=FALSE,...) {
  
  if (x$type%in%c("flex","u","sat")) {
    cat(x$estimate,sep="\n")
  } else {

    printCoefmat(x$estimate,signif.stars=signif.stars,...)
    cat("\n")
    print(x$zyg,quote=FALSE)

    if (!is.null(x$acde)) {
      cat("\nVariance decomposition:\n")
      printCoefmat(x$acde,...)
    }
    cat("\n\n")
    ##    cat("Broad-sense heritability (total genetic factors):\n")
    h <- with(x, heritability[,c(1,3,4),drop=FALSE]);
    h <- na.omit(h)
    printCoefmat(h,...)  
    cat("\n")
  }
  if (!is.null(x$corMZ)) {
    cc <- with(x, rbind(corMZ,corDZ)) ##,corOS))
    rownames(cc)[1:2] <- c("Correlation within MZ:","Correlation within DZ:")
    ## if (!is.null(x$corOS)) rownames(cc)[3] <- "Correlation within OS:"
    colnames(cc) <- c("Estimate","2.5%","97.5%")
    printCoefmat(cc,signif.stars=FALSE,...)
  }
  
  cat("\n")
  print(x$logLik)
  cat("AIC:", x$AIC, "\n")
  cat("BIC:", x$BIC, "\n")
  invisible(x)
}

###}}} print.summary.twinlm

###{{{ compare.twinlm

##' @S3method compare twinlm
compare.twinlm <- function(object,...) {
  if (length(list(...))==0) return(compare(object$estimate))
  lava:::compare.default(object,...)
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

###{{{ score.twinlm
##' @S3method score twinlm
score.twinlm <- function(x,...) score(x$estimate,...)
###}}} score.twinlm

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
  M <- pnorm(constraints(x$estimate,k=1)[,c(1,5,6),drop=FALSE])
##  if (x$OS) M <- M[-nrow(M),,drop=FALSE]
  rownames(M) <- toupper(rownames(M))
  M
}

###}}} acde

coef.twinlm <- function(object,...) coef(object$estimate,...)
