Grandom.cif<-function(cif,data=sys.parent(),cause,
causeS=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL,theta=NULL,theta.des=NULL,step=1,design.rv=NULL,link=0,
notaylor=1,same.cens=FALSE,entry=NULL,trunkp=1)
{
if (is.null(clusters)) stop("Clusters must be specified to estimate correlation\n");
### link=1 uses
## {{{ setting up variables
 multi<-0; inverse<-link; sdscore=1; 
 # trans=1 P_1=1-exp(- ( x' b(b)+ z' gam diag(t^time.pow))), 
 formula<-attr(cif,"Formula"); 
 ldata<-aalen.des(formula,data,model="aalen");
 X<-ldata$X; time<-ldata$time2; Z<-ldata$Z;  status<-ldata$status;
 #print(dim(X)); print(dim(Z)); 
  antpers<-nrow(X); 
  if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else {Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
  if (npar==TRUE) {Z<-matrix(0,antpers,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-ncol(Z);} 
  ng<-antpers;px<-ncol(X);  
 #print(dim(X)); print(dim(Z)); 
  times<-cif$cum[,1]; delta<-(cause!=cens.code)
 if (semi==1) gamma<-c(cif$gamma)  else gamma<-0; 

if (is.null(entry)) entry.call <- NULL else entry.call <- 0
if (is.null(entry)) { entry <- rep(0,antpers);  cif1lin <- entry;} else {
   cum1<-Cpred(rbind(rep(0,px+1),cif$cum),entry)[,-1];
   cif1lin  <-  (Z %*% gamma )*entry + apply(X*cum1,1,sum) 
}
Gcxe <- 1; 
if (length(trunkp)==1) trunkp <- rep(1,antpers)
## }}}

if (is.null(design.rv)==TRUE) design.rv <- matrix(1,antpers,1); 
dim.rv <- ncol(design.rv); 
###print(head(design.rv)); print(head(theta.des)) 

## {{{ censoring weights
  if (cens.model=="KM") {
    ud.cens<-survfit(Surv(time,cause==cens.code)~+1); 
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time)[,2];
    if (is.null(entry.call)==FALSE) Gcxe<-Cpred(Gfit,entry)[,2];
    if (min(Gcx)< 0.00001) { 
	    cat("Censoring dist. zero for some points, summary cens:\n");
	    print(summary(Gcx)) 
    }
    Gcx <- Gcx/Gcxe; 
    Gctimes<-Cpred(Gfit,times)[,2];
  }
  if (cens.model=="cox") {
    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
    ud.cens<-cox.aalen(Surv(time,cause==cens.code)~prop(XZ),n.sim=0,robust=0);
    Gcx<-Cpred(ud.cens$cum,time)[,2];
    RR<-exp(XZ %*% ud.cens$gamma)
    Gcx<-exp(-Gcx*RR)
    Gfit<-rbind(c(0,1),cbind(time,Gcx)); 
    Gctimes<-Cpred(Gfit,times)[,2];
  }
  if (cens.model == "aalen") { 
	  if (npar == TRUE) XZ <- X[, -1] else XZ <- cbind(X, Z)[, -1]
     ud.cens <- aalen(Surv(time, cause == cens.code) ~ XZ,n.sim = 0, robust = 0)
     Gcx <- Cpred(ud.cens$cum, time)[, -1]
     XZ <- cbind(1, XZ) 
     Gcx <- exp(-apply(Gcx * XZ, 1, sum))
     Gcx[Gcx < 0] <- 0
     Gfit <- rbind(c(0, 1), cbind(time, Gcx))
     Gctimes <- Cpred(Gfit, times)[, 2] 
  }
  ntimes<-length(times); 
## }}}
###
##### {{{ censoring weights
###  if (cens.model=="KM") {
###    ud.cens<-survfit(Surv(time,cause==cens.code)~+1); 
###    Gfit<-cbind(ud.cens$time,ud.cens$surv)
###    Gfit<-rbind(c(0,1),Gfit); 
###    Gcx<-Cpred(Gfit,time)[,2];
###    Gctimes<-Cpred(Gfit,times)[,2];
###  }
###  if (cens.model=="cox") {
###    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
###    ud.cens<-cox.aalen(Surv(time,cause==cens.code)~prop(XZ),n.sim=0,robust=0);
###    Gcx<-Cpred(ud.cens$cum,time)[,2];
###    RR<-exp(XZ %*% ud.cens$gamma)
###    Gcx<-exp(-Gcx*RR)
###    Gfit<-rbind(c(0,1),cbind(time,Gcx)); 
###    Gctimes<-Cpred(Gfit,times)[,2];
###  }
###  ntimes<-length(times); 
##### }}}
###

## {{{ cluster + definining variables
  if (is.null(clusters)== TRUE) {clusters<-0:(antpers-1); antclust<-antpers;} else {
     clus<-unique(clusters); antclust<-length(clus);
     clusters <- as.integer(factor(clusters, labels = 1:(antclust)))-1; 
  }

  out <- cluster.index(clusters); 
  clustsize <- out$cluster.size
  maxclust <- out$maxclust
  clusterindex <- out$idclust
  if (maxclust==1) stop("No clusters given \n"); 

###  clustsize<-as.vector(table(clusters));
###  maxclust<-max(clustsize); clusterindex<-matrix(0,antclust,maxclust);
###   cs<- rep(1,antclust)
###   for (i in 1:antpers) {
###      clusterindex[clusters[i],cs[clusters[i]]]<-i-1;
###      cs[clusters[i]]<- cs[clusters[i]]+1;
###   }

 if (is.null(theta.des)==TRUE) ptheta<-dim.rv^2 else ptheta <- ncol(theta.des); ;
 if (is.null(theta.des)==TRUE) theta.des<-matrix(diag(dim.rv),antpers,dim.rv^2,byrow=TRUE) else theta.des<-as.matrix(theta.des,antpers,ptheta);
  ptheta <- ptheta/dim.rv; 

  if (is.null(theta)==TRUE) theta<-rep(0.5,ptheta);
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta);

  ##print(est); print(gamma); print(semi); 
  est<-cif$cum; if (semi==1) gamma<-cif$gamma  else gamma<-0; 
  Biid<-c()
  if (notaylor==0) { for (i in 1:antclust) Biid<-cbind(Biid,cif$B.iid[[i]])} else Biid <- 0; 
  if (npar==TRUE) gamma.iid<-0 else gamma.iid<-cif$gamma.iid;
  var.theta<-hess<-matrix(0,ptheta,ptheta); score<-rep(0,ptheta); 
  time.pow<-attr(cif,"time.pow"); 
  if (sum(time.pow)==0) time.pow<-rep(1,pg); 
## }}}

  out<-.C("rcifdes", ## {{{
  as.double(times),as.integer(ntimes),as.double(time),
  as.integer(delta), as.integer(cause), as.integer(causeS),
  as.double(Gcx), as.double(X),as.integer(antpers),
  as.integer(px), as.integer(Nit), as.double(score),
  as.double(hess), as.double(est), as.double(gamma), 
  as.integer(semi),as.double(Z), as.integer(pg),
  as.integer(detail), as.double(Biid),as.double(gamma.iid),
  as.double(time.pow), as.double(theta), as.double(var.theta), 
  as.double(theta.des), as.integer(ptheta), as.integer(antclust), 
  as.integer(clusters), as.integer(clustsize), as.integer(clusterindex),
  as.integer(maxclust), as.double(step) , 
  as.integer(inverse), as.integer(sdscore),
  as.double(design.rv),as.integer(dim.rv),as.integer(notaylor),
  as.integer(same.cens),
  as.double(trunkp), as.double(entry),as.double(cif1lin), 
  PACKAGE="MultiComp") ## }}}

## {{{ output 
theta<-matrix(out[[23]],ptheta,1);var.theta<-matrix(out[[24]],ptheta,ptheta); 
score<-matrix(out[[12]],ptheta,1); hess<-matrix(out[[13]],ptheta,ptheta); 

if (ptheta>1) names<-colnames(theta.des) else names<-"intercept"

if (is.null(names)==FALSE && length(names)==ptheta) {
   rownames(score)<-rownames(theta)<-rownames(var.theta)<-colnames(var.theta)<-
   rownames(hess)<-colnames(hess)<- colnames(hess)<-names;
   colnames(theta)<-"coefficient"; colnames(score)<-"score";
}


ud<-list(score=score,hess,theta=theta,var.theta=var.theta)
ud$call<-call; ud$antpers<-antpers; ud$antclust<-antclust; 
class(ud)<-"randomcifrv"; 
attr(ud, "Call") <- sys.call()
attr(ud, "Formula") <- formula
attr(ud, "Clusters") <- clusters
attr(ud,"inverse")<-inverse
attr(ud,"cause1")<-causeS
attr(ud,"cause2")<-causeS
return(ud); 
## }}}
}

summary.randomcifrv<-function (object, ...) 
{ ## {{{
    if (!inherits(object, "randomcifrv")) 
        stop("Must be a random.cifrv  object")
    cat("Random effect parameters for additive gamma random effects \n\n")
    cat("Cause", attr(object, "cause1"), "and cause", attr(object, 
        "cause2"), fill = TRUE)
    cat("\n")
    if (sum(abs(object$score)) > 1e-06) 
        cat("WARNING: check score for convergence")
    cat("\n")
    coef.randomcifrv(object, ...)
} ## }}}

coef.randomcifrv<- function (object, digits = 3, ...) 
{ ## {{{
    if (attr(object,"inverse")==1) elog <- 1 else elog  <- 0; 
    if (elog==1) theta <- exp(object$theta) else theta <- object$theta 
    se <- diag(object$var.theta)^0.5
    res <- cbind(object$theta, se)
    wald <- object$theta/se
    waldp <- (1 - pnorm(abs(wald))) * 2
    res <- as.matrix(cbind(res, wald, waldp))
    if (elog==0)  colnames(res) <- c("Coef.", "SE", "z", "P-val")
    if (elog==1) res <- cbind(res,exp(object$theta), exp(object$theta)^2*se)  
    if (elog==1) colnames(res) <- c("log-parameter","SE","z","P-val","exp(theta)","SE")

    if (is.null((rownames(res))) == TRUE) rownames(res) <- rep(" ", nrow(res))
    prmatrix(signif(res, digits))

    cat("\n\n Random effect variances for gamma random effects \n\n")
    varpar <- theta/sum(theta)^2 
    res <- as.matrix(varpar); 
    if (elog==0)  { var.theta <-   object$var.theta; 
                    df <- 0*var.theta; 
                    for (i in 1:nrow(var.theta))
                    df[i,] <- -theta[i]*2*theta; 
		    diag(df) <- diag(df)+sum(theta)^2
		    df <- df/sum(theta)^4
		    var.varpar <- df %*% var.theta %*% df
                  }
    if (elog==1)  { 
	            var.theta <-   object$var.theta; 
                    var.varpar <- var.theta
                  }
    res <- cbind(res,diag(var.varpar)^.5)  
    colnames(res) <- c("variance","SE")
    if (is.null((rownames(res))) == TRUE) rownames(res) <- rep(" ", nrow(res))
    prmatrix(signif(res, digits))

} ## }}}

print.randomcifrv<- function (x , digits = 3, ...) 
{ ## {{{
 summary(x, ...)
} ## }}}
