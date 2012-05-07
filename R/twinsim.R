##' Simulate twin data from a linear normal ACE/ADE/AE model.
##'
##' @title Simulate twin data
# ##' @return A \code{list} with the following elements
# ##'     \item{data}{Data in long format (one row for each individual.}
# ##'     \item{model}{The multigroup structural equation model.}
# ##'     \item{wide}{A list of data.frames (MZ data and DZ data) in wide format
# ##'     (one row for each pair) eady for a structural equation model.}
##' @author Klaus K. Holst
##' @export
##' @seealso \code{\link{twinlm}}
##' @keywords models
##' @keywords regression
##' @param nMZ Number of monozygotic twin pairs
##' @param nDZ Number of dizygotic twin pairs
##' @param b1 Effect of covariates (labelled x1,x2,...) of type 1. One
##'     distinct covariate value for each twin/individual.
##' @param b2 Effect of covariates (labelled g1,g2,...) of type 2. One
##'     covariate value for each twin pair.
##' @param mu Intercept parameter.
##' @param acde Variance of random effects (in the order A,C,D,E)
##' @param randomslope Logical indicating wether to include random slopes of
##'     the variance components w.r.t. x1,x2,...
##' @param threshold Treshold used to define binary outcome y0
##' @param cens Logical variable indicating whether to censor outcome
##' @param wide Logical indicating if wide data format should be returned
##' @param ... Additional arguments parsed on to lower-level functions
twinsim <- function(nMZ=100,nDZ=nMZ,b1=c(),b2=c(),mu=0,acde=c(1,1,0,1),randomslope=NULL,threshold=0,cens=FALSE,wide=FALSE,...) {
  n <- nMZ+nDZ
  sA <- acde[1]^0.5; sC <- acde[2]^0.5; sD <- acde[3]^0.5; sE <- acde[4]^0.5; 
  A.MZ <- rnorm(nMZ,sd=sA); A.MZ <- cbind(A.MZ,A.MZ)
  D.MZ <- rnorm(nMZ,sd=sD); D.MZ <- cbind(D.MZ,D.MZ)  
  S2 <- matrix(c(0,1,1,0),2)
  A.DZ <- sA*rmvnorm(nDZ,sigma=diag(2)+S2*0.5)
  D.DZ <- sD*rmvnorm(nDZ,sigma=diag(2)+S2*0.25)
  C.MZ <- rnorm(nMZ,sd=sC)
  C.DZ <- rnorm(nDZ,sd=sC)
  yMZ <- mu + A.MZ + cbind(C.MZ,C.MZ) + D.MZ + cbind(rnorm(nMZ,sd=sE),rnorm(nMZ,sd=sE))
  yDZ <- mu + A.DZ + cbind(C.DZ,C.DZ) + D.DZ + cbind(rnorm(nDZ,sd=sE),rnorm(nDZ,sd=sE))
  y <- rbind(yMZ,yDZ)
  if (length(b1)>0) {
    x1 <- rmvnorm(n,rep(0,length(b1)),diag(length(b1)))
    x2 <- rmvnorm(n,rep(0,length(b1)),diag(length(b1)))
    y <- y+cbind(x1%*%b1,x2%*%b1)
  }
  if (length(b2)>0) {
    g <- rmvnorm(n,rep(0,length(b2)),diag(length(b2)))
    ge <- g%*%b2
    y <- y+cbind(ge,ge)
  }
  Cens <- ifelse(cens,rep(Inf,nMZ+nDZ),rnorm(nMZ+nDZ,threshold+1))
  d <- data.frame(id=seq(n),y=y,zyg=c(rep("MZ",nMZ),rep("DZ",nDZ)),
                  cens=Cens)
  vary <- list(c("y1","y2"))
  colnames(d)[2:3] <- vary[[1]]
  if (length(b1)>0) {
    d <- cbind(d,x1=x1,x2=x2);
    vary <- c(vary,lapply(seq(length(b1)),FUN=function(x) paste(c("x1","x2"),x,sep="")))
  }
  if (length(b2)>0) { d <- cbind(d,g=g) }
  colnames(d) <- sub(".","",colnames(d),fixed=TRUE)
  if (wide) return(d)
  dd <- reshape(d,direction="long",varying=vary)
  dd <- transform(dd,y=(y1>threshold & y1<cens)*1,y0=(y1>threshold),status=y1<cens) 
  ## S.MZ <- diag(2)*vE+vC+vA
  ## S.DZ <- diag(2)*(vE+vA) + rho*S2*vA + vC
  ## Mu <- (threshold-mu)
  ## probs <- c(marginal=pnorm(mu,threshold,S.MZ[1]^0.5),concordance=pmvnorm(lower=c(threshold,threshold),mean=c(mu,mu),sigma=S.MZ))
  ## probs <- c(probs,casewise=probs[2]/probs[1])
  ## probs
  return(dd)
}

