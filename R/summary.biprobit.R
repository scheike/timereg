##' @export
summary.biprobit <- function(object,level=0.05,transform,...) {
  alpha <- level/2
  varcomp <- object$coef[length(coef(object)),1:2]
  varcomp <- rbind(object$model$tr(c(varcomp[1],varcomp[1]%x%cbind(1,1) + qnorm(1-alpha)*varcomp[2]%x%cbind(-1,1))))
  colnames(varcomp)[2:3] <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  rownames(varcomp) <- ifelse(is.null(object$model$varcompname),"Variance component",object$model$varcompname)


  h <- function(p) log(p/(1-p)) ## logit
  ih <- function(z) 1/(1+exp(-z)) ## expit
  ##dlogit <- function(p) 1/(p*(1-p))
  if (!missing(transform)) {
      h <- asin; ih <- sin      
      if (is.null(transform)) {
          h <- ih <- identity
      }
      if (is.list(transform)) {
          h <- transform[[1]]; ih <- transform[[2]]
      }     
  }
  probs <- function(p) {
    ##    S <- diag(2); S[1,2] <- S[2,1] <- exp(tail(p,1))
    S <- object$SigmaFun(p[length(p)])
    m <- c(0,0)
    if (object$npar$intercept==1) m[1:2] <- p[1]
    if (object$npar$intercept==2) {
        idx <- 1:2
        if (length(object$npar$pred)>0) idx <- c(1,1+object$npar$pred/2+1)
        m[1:2] <- p[idx]
    }
    mu.cond <- function(x) m[1]+S[1,2]/S[2,2]*(x-m[2])
    var.cond <- S[1,1]-S[1,2]^2/S[2,2]
    conc <- pmvn(upper=m,sigma=S)
    disconc <- pmvn(lower=c(-Inf,m[1]),upper=c(m[2],Inf),sigma=S)
    marg <- pnorm(m[1],sd=S[1,1]^0.5)
    cond <- conc/marg
    lambda <- cond/marg
    discond <- disconc/(1-marg)
    logOR <- log(cond)-log(1-cond)-log(discond)+log(1-discond)
    c(h(c(conc,cond,marg)),lambda,logOR)
  }
  alpha <- level/2
  CIlab <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  mycoef <- coef(object)
  prob <- probs(mycoef)
  Dprob <- numDeriv::jacobian(probs,mycoef)
  sprob <- diag((Dprob)%*%vcov(object)%*%t(Dprob))^0.5
  pp <- cbind(prob,prob-qnorm(1-alpha)*sprob,prob+qnorm(1-alpha)*sprob)
  pp[1:3,] <- ih(pp[1:3,])
  nn <- c("Concordance","Casewise Concordance","Marginal","Rel.Recur.Risk","log(OR)")
  if (nrow(pp)-length(nn)>0) nn <- c(nn,rep("",nrow(pp)-length(nn)))
  rownames(pp) <- nn
  colnames(pp) <- c("Estimate",CIlab)
  
  res <- list(all=rbind(varcomp,pp),varcomp=varcomp,prob=pp,coef=object$coef,score=colSums(object$score),logLik=object$logLik,msg=object$msg,N=object$N)
  class(res) <- "summary.biprobit"
  res
}
