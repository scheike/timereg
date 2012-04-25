simcif <- function(n,theta=list(c(-10,1,0.1),c(-10,0.1,.2)),
                   Sigma, pr,
                   range=c(0,100),
                   a=function(t,theta) theta[1]+theta[2]*t+theta[3]*exp(t),
                   ia, cens,
                   X,
                   bX,
                   Z,
                   ZtIdx
                ) {
  if (!is.list(theta)) theta <- list(theta)
  ncauses <- length(theta)
  if (missing(Sigma)) {
    Sigma <- diag(2*ncauses)+1
    for (i in seq_len(ncauses)) {
      Sigma[1:2+(i-1)*2,1:2+(i-1)*2] <- Sigma[1:2+(i-1)*2,1:2+(i-1)*2]+i
    }    
  }
  np <- (ncauses+1)*(ncauses)/2
  if (missing(pr)) {
    pr <- 1/2^seq_len(np-1)
  }
  pr <- c(pr,1-sum(pr))
  P <- matrix(ncol=ncauses,nrow=ncauses)
  if (ncauses==1) {
    P[1,1] <- 1
  } else {    
    P[upper.tri(P,diag=TRUE)] <- pr
    P[lower.tri(P)] <- P[upper.tri(P)] <- P[upper.tri(P)]/2
  }

  invexpl <- !missing(ia)
  t. <- seq(range[1],range[2],length.out=1000)
  pa <- matrix(ncol=1+ncauses,nrow=length(t.))
  pa[,1] <- t.
  aa <- pa

  simcauses <- rmultinom(1,n,as.vector(P))
  T <- Y <- matrix(ncol=4,nrow=n)
  N <- 0
  pos <- 0

  for (j in seq_len(ncauses)) {
    for (i in seq_len(ncauses)) {
      pos <- pos+1
      if (simcauses[pos]>0) {
        Nprev <- N+1
        N <- N + simcauses[pos]
        pseq <- seq(Nprev,N)
        T[pseq,3] <- Y[pseq,3] <- i
        T[pseq,4] <- Y[pseq,4] <- j
        idx <- c(1+(i-1)*2,2+(j-1)*2)
        S <- Sigma[idx,idx,drop=FALSE]
        val <- rnorm(simcauses[pos],sd=S[1,1]^0.5)
        myMean <- c(0,0)
        if (!missing(X)) {
          if (missing(bX)) bX <- rep(1,ncol(X))
          myMean <- (X%*%bX) %x% cbind(1,1)
          val <- val + myMean[,1]
        }
        cm <- CondMom(myMean,S,2,X=cbind(val))
        Y[pseq,1] <- val
        Y[pseq,2] <- with(cm, rnorm(simcauses[pos],mean=as.vector(mean),sd=var^0.5))
        
        theta0 <- theta[[i]];
        a. <- a(t.,theta0)
        if (!invexpl) {
          ia <- function(x,...) {
            fastapprox(a.,x,t.)$t
          }
        }
        T[pseq,1] <- ia(Y[pseq,1],theta0)
        theta0 <- theta[[j]]; # subject two
        a. <- a(t.,theta0)
        if (!invexpl) {
          ia <- function(x,...) {
            fastapprox(a.,x,t.)$t
          }
        }      
        T[pseq,2] <- ia(Y[pseq,2],theta0)        
      }
    }
  }
  
  for (i in seq_len(ncauses)) { ## Marginals
    theta0 <- theta[[i]];
    a. <- a(t.,theta0)
    aa[,i+1] <- a.
    pa[,i+1] <- pnorm(a.,sd=Sigma[i*2,i*2]^0.5)
  }
  
  if (!missing(cens)) {    
    if (is.function(cens)) {
      cens <- cens(T[,1:2])
    }    
    T[T[,1]>=cens[,1],3] <- 0
    T[T[,2]>=cens[,2],4] <- 0
    T[T[,3]==0,1] <- cens[T[,3]==0,1]
    T[T[,4]==0,2] <- cens[T[,4]==0,2]
  }
  colnames(T) <- c("t1","t2","cause1","cause2")
  res <- list(data=T,prob=pa,var=Sigma,P=P,a=aa)
  class(res) <- "simt"
  return(res)
}

print.simt <- function(x,...) {
  ##  cat("(t1,t2,cause1,cause2):\n")
  with(x, print(data))
  cat("Variance:\n")  
  with(x, print(var))
  cat("Prob. causes:\n")  
  with(x, print(P))
  invisible(x)
}
