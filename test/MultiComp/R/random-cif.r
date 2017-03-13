random.cif<-function(cif,data=sys.parent(),cause,
causeS=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL,theta=NULL,theta.des=NULL,step=1)
{
## {{{ setting up variables
 multi<-0; inverse<-0; sdscore=1; 
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
## }}}

## {{{ censoring weights
  if (cens.model=="KM") {
    ud.cens<-survfit(Surv(time,cause==cens.code)~+1); 
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time)[,2];
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

## {{{ cluster + definining variables
  if (is.null(clusters)== TRUE) {clusters<-0:(antpers-1); antclust<-antpers;} else {
    clus<-unique(clusters); antclust<-length(clus);
    clusters <- as.integer(factor(clusters, labels = 1:(antclust)));
  }
  clustsize<-as.vector(table(clusters));
  maxclust<-max(clustsize); clusterindex<-matrix(0,antclust,maxclust); 

 cs<- rep(1,antclust)
 for (i in 1:antpers) { 
     clusterindex[clusters[i],cs[clusters[i]]]<-i-1;
     cs[clusters[i]]<- cs[clusters[i]]+1; 
  } 

##################################################################
 if (is.null(theta.des)==TRUE) ptheta<-1;
 if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) else 
 theta.des<-as.matrix(theta.des);
  ptheta<-ncol(theta.des);
  if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim"); 
  if (is.null(theta)==TRUE) theta<-rep(0.5,ptheta);
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta);
  est<-cif$cum; if (semi==1) gamma<-cif$gamma  else gamma<-0; 
  ##print(est); print(gamma); print(semi); 
  Biid<-c()
  for (i in 1:antclust) Biid<-cbind(Biid,cif$B.iid[[i]]); 
  if (npar==TRUE) gamma.iid<-0 else gamma.iid<-cif$gamma.iid;
  var.theta<-hess<-matrix(0,ptheta,ptheta); score<-rep(0,ptheta); 
  time.pow<-attr(cif,"time.pow"); 
  if (sum(time.pow)==0) time.pow<-rep(1,pg); 
## }}}

#dyn.load("random-cif.so");

  out<-.C("randomcif", ## {{{
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
  as.integer(inverse), as.integer(sdscore), PACKAGE="MultiComp") ## }}}

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
class(ud)<-"randomcif"; 
attr(ud, "Call") <- sys.call()
attr(ud, "Formula") <- formula
attr(ud, "Clusters") <- clusters
attr(ud,"inverse")<-inverse
attr(ud,"cause1")<-causeS
attr(ud,"cause2")<-causeS
return(ud); 
## }}}
}

summary.randomcif<-function (object, ...) 
{
    if (!inherits(object, "randomcif")) 
        stop("Must be a random.cif  object")
    cat("Random effect variance for variation due to clusters\n\n")
    cat("Cause", attr(object, "cause1"), "and cause", attr(object, 
        "cause2"), fill = TRUE)
    cat("\n")
    if (sum(abs(object$score)) > 1e-06) 
        cat("WARNING: check score for convergence")
    cat("\n")
    coef.randomcif(object, ...)
}

coef.randomcif<- function (object, digits = 3, ...) 
{
    res <- cbind(object$theta, diag(object$var.theta)^0.5)
    se <- diag(object$var.theta)^0.5
    wald <- object$theta/se
    waldp <- (1 - pnorm(abs(wald))) * 2
    cor <- object$theta + 1
    res <- as.matrix(cbind(res, wald, waldp, cor, se))
    colnames(res) <- c("Coef.", "SE", "z", "P-val", "Cross odds ratio", 
        "SE")
    if (is.null((rownames(res))) == TRUE) 
        rownames(res) <- rep(" ", nrow(res))
    prmatrix(signif(res, digits))
}

print.randomcif<- function (x , digits = 3, ...) 
{
}
