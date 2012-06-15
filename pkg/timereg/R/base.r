coefBase<- function(object, digits=3, d2logl=0) { ## {{{
  res <- cbind(object$gamma,
               diag(object$var.gamma)^0.5,
               diag(object$robvar.gamma)^0.5)
  if (d2logl==1) res<-cbind(res,diag(object$D2linv)^.5)
  wald <- object$gamma/diag(object$robvar.gamma)^0.5
  waldp <- (1 - pnorm(abs(wald))) * 2
  res <- as.matrix(cbind(res, wald, waldp))
  if (d2logl==1) colnames(res) <- c("Coef.", "SE", "Robust SE","D2log(L)^-1","z","P-val") else colnames(res) <- c("Coef.", "SE", "Robust SE", "z", "P-val")
###  prmatrix(signif(res, digits))
###  cat("\n")
  
  return(res)
} ## }}}

wald.test <- function(object,contrast,coef.null=NULL,Sigma=NULL,null=NULL)
{ ## {{{
  if (is.null(Sigma)) {
     if (class(object)=="cor") Sigma <- object$var.theta else Sigma <- object$var.gamma;
  }
  coefs <- coefBase(object)[,1]
  nl <- length(coefs)
  if (missing(contrast)) {
      contrast <- rep(1,length(coefs))
      contrast <- diag(1,nl);
  }
  if (!is.null(coef.null)) {
	  contrast <- c()
	  for (i in coef.null) 
      contrast <- rbind(contrast,c((1:nl)==i)*1)
  }
  if (missing(null)) null <- 0

  ### Wald test
  B <- contrast
  p <- coefs
  if (is.vector(B)) { B <- rbind(B); colnames(B) <- names(contrast) }

### if (ncol(B)<length(p)) {
###    nn <- colnames(B)
###    myidx <- parpos(Model(object),p=nn)
###    B0 <- matrix(0,nrow=nrow(B),ncol=length(coef(object)))
###    B0[,myidx] <- B[,attributes(myidx)$ord]
###    B <- B0
### }
 Q <- t(B%*%p-null)%*%solve(B%*%Sigma%*%t(B))%*%(B%*%p-null)
 df <- qr(B)$rank; names(df) <- "df"
 attributes(Q) <- NULL; names(Q) <- "chisq";
 pQ <- ifelse(df==0,NA,1-pchisq(Q,df))
 method = "Wald test";
 ##    hypothesis <-
 res <- list(##data.name=hypothesis,
  statistic = Q, parameter = df, p.value=pQ, method = method)
  class(res) <- "htest"
  attributes(res)$B <- B
return(res)
} ## }}}

timetest<-function(object,digits=3,hyp.label="p-value H_0:constant effect",out=0)
{  ## {{{
  cat("Test for nonparametric terms \n")
  if (is.null(object$conf.band)==TRUE)  mtest<-FALSE else mtest<-TRUE;
  if (mtest==FALSE) cat("Test not computed, sim=0 \n\n")
  if (mtest==TRUE) {
  test0<-cbind(object$obs.testBeq0,object$pval.testBeq0)
  testC<-cbind(object$obs.testBeqC,object$pval.testBeqC)
  colnames(test0)<-c("Supremum-test of significance","p-value H_0: B(t)=0")
  colnames(testC)<-c("      Kolmogorov-Smirnov test",hyp.label)
  if (is.null(object$obs.testBeqC.is)!=TRUE)  {
  testCis<-cbind(object$obs.testBeqC.is,object$pval.testBeqC.is)
  colnames(testCis) <-
                   c("        Cramer von Mises test",hyp.label)
  }
  cat("\n")
  cat("Test for non-significant effects \n")
  prmatrix(signif(test0,digits))
  cat("\n")
  cat("Test for time invariant effects \n")
  prmatrix(signif(testC,digits))
  if (is.null(object$obs.testBeqC.is)!=TRUE)  prmatrix(signif(testCis,digits))
  cat("\n")
  if (out==1) return(cbind(test0,testC)); 
}
} ## }}}

is.diag <- function(m)
{
p <- nrow(m)
adiag <- min(diag(m)*1)
if (adiag==0) ud <- FALSE else ud <- TRUE
dm <- diag(p); diag(dm) <- diag(m); 
ndiag <- sum(abs(c(m - dm)))
if (ndiag>0.0000001) ud <- FALSE;
return(ud)
}

cluster.index <- function(clusters,index.type=FALSE)
{ ## {{{
 antpers <- length(clusters)
if (index.type==FALSE)  {
	max.clust <- length(unique(clusters))
	clusters <- as.integer(factor(clusters, labels = 1:max.clust)) -1
}
 nclust <- .C("nclusters",
	as.integer(antpers), as.integer(clusters), as.integer(rep(0,antpers)), 
	as.integer(0), as.integer(0), package="timereg")
  maxclust <- nclust[[5]]
  antclust <- nclust[[4]]
  cluster.size <- nclust[[3]][1:antclust]
  clustud <- .C("clusterindex",as.integer(clusters),
		as.integer(antclust),as.integer(antpers),
                as.integer(rep(0,antclust*maxclust)),as.integer(rep(0,antclust)),
	  package="timereg")
idclust <- matrix(clustud[[4]],antclust,maxclust)

out <- list(clusters=clusters,maxclust=maxclust,antclust=antclust,idclust=idclust,cluster.size=cluster.size)
} ## }}}

residualsTimereg <- function(object,data=data)
{ ## {{{

if (class(object)!="cox.aalen" & class(object)!="aalen") stop("Computes residuals for Aalen or Cox.aalen object") 
else {
 formula<-attr(object,"Formula");
 beta.fixed <- attr(object,"beta.fixed")
 if (is.null(beta.fixed)) beta.fixed <- 1; 
 model <- class(object); 
 ldata<-aalen.des(formula,data=data,model=model);
 id <- attr(object,"id"); 
 mclusters <- attr(object,"cluster")
 X<-ldata$X; 
 time<-ldata$time2; 
 Z<-ldata$Z;  
 status<-ldata$status;
 time2 <- attr(object,"stop"); 
 start <- attr(object,"start");
 status <- attr(object,"status");
 if (!is.null(attr(object,"max.time"))) status <- status*(time2< attr(object,"max.time")); 
 antpers<-nrow(X);
 if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else { Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
 if (npar==TRUE) {Z<-matrix(0,antpers,1); pz<-1; fixed<-0;} else {fixed<-1;pz<-ncol(Z);}
 px<-ncol(X);

 if (sum(abs(start))>0) lefttrunk <- 1  else lefttrunk <- 0;  cumhazleft <- 0; 
 nn <- nrow(object$cum) 

 cum <- Cpred(object$cum,time2)[,-1]
 cumhaz0 <- apply(cum*X,1,sum)
 cumhazleft <- rep(0,antpers)
 RR <- rep(1,antpers); 

if (class(object)=="cox.aalen")
{ ## {{{
  RR <- exp(Z %*% object$gamma); 
  cumhaz <- cumhaz0 * RR;
  if (lefttrunk==1) {
      cum <- Cpred(object$cum,start)[,-1]
      cumhazleft <- apply(cum*X,1,sum)
      cumhazleft <- cumhazleft * RR;
  }
} ## }}}

if (class(object)=="aalen")
{#{{{
  if (npar==FALSE) { ## semi-parametric risk model
      ex.haz <- (Z %*% object$gamma) ; 
      cumhaz <- cumhaz0+ex.haz*time2
     if (lefttrunk==1) {
	 cum <- Cpred(object$cum,start)[,-1]
	 cumhazleft <- apply(cum*X,1,sum)
	 cumhazleft  <-  cumhazleft+ex.haz*start
     }

  } else {  ## Aalen model
	  cumhaz <- cumhaz0
          if (lefttrunk==1) {
	     cum <- Cpred(object$cum,start)[,-1]
	     cumhazleft <- apply(cum*X,1,sum)
	     if (npar==TRUE) cumhazleft <-  cumhazleft
	     ex.haz <- (Z %*% object$gamma) * time2 ; 
          }

  }
} #}}}

} 

residuals <- status- cumhaz
out <- list(residuals=c(residuals),status=c(status),cumhaz=c(cumhaz),cumhazleft=c(cumhazleft),RR=RR)
} ## }}}
