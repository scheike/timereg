random.cif2cause<-function(cif,data=sys.parent(),cause,
cause1=1,cause2=1,cens.code=0,cif2=NULL, 
cens.model="KM",Nit=40,detail=0,clusters=NULL,theta=NULL,
theta.des=NULL,step=1,squarepar=0,same.cens=FALSE)
{
## {{{ setting things up
if (cause1==cause2) stop("use random.cif when considering two equivalent causes"); 
if (is.null(cif2)==TRUE) stop("must give both marginal estimates"); 
 multi=0;inverse=0; sdscore=1; sym=1; 
 # trans=1 P_1=1-exp(- ( x' b(b)+ z' gam t) ), 
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
  est<-cif$cum; if (semi==1) gamma<-cif$gamma  else gamma<-0; 
if (cause1!=cause2) {
if (is.null(cif2)==TRUE) stop("Must provide marginal model for both causes"); 
 formula2<-attr(cif2,"Formula"); 
 ldata2<-aalen.des(formula2,data,model="aalen");
 X2<-ldata2$X; timec2<-ldata2$time2; Z2<-ldata$Z;  status2<-ldata2$status;
 if (is.null(Z2)==TRUE) {npar2<-TRUE; semi2<-0;}  else {Z2<-as.matrix(Z2); npar2<-FALSE; semi2<-1;}
 if (npar2==TRUE) {Z2<-matrix(0,antpers,1); pg2<-1; fixed2<-0;} else {fixed2<-1;pg2<-ncol(Z2);} 
 px2<-ncol(X2);  
 ntimes<-length(times); 
 est2<-cif2$cum; if (semi2==1) gamma2<-cif2$gamma  else gamma2<-0; 
 est<-Cpred(est,times); 
 est2<-Cpred(est2,times);
####      print(est); print(est2)
if (nrow(cif$cum)!=nrow(cif2$cum)) stop("Marginal models must be computed in same
time points \n"); 
}
## }}}

## {{{ censoring weights
  if (cens.model=="KM") {
    ud.cens<-survfit(Surv(time,cause==cens.code)~+1); 
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time)[,2];
    Gctimes<-Cpred(Gfit,times)[,2];
  } else if (cens.model=="cox") {
    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
    ud.cens<-cox.aalen(Surv(time,cause==cens.code)~prop(XZ),n.sim=0,robust=0);
    Gcx<-Cpred(ud.cens$cum,time)[,2];
    RR<-exp(XZ %*% ud.cens$gamma)
    Gcx<-exp(-Gcx*RR)
    Gfit<-rbind(c(0,1),cbind(time,Gcx)); 
    Gctimes<-Cpred(Gfit,times)[,2];
  } else stop("Censoring model not specified for this function \n"); 
## }}}

## {{{ clusters and definition of variables
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
###
### cs<- rep(1,antclust)
### for (i in 1:antpers) { 
###     clusterindex[clusters[i],cs[clusters[i]]]<-i-1;
###     cs[clusters[i]]<- cs[clusters[i]]+1; 
###  } 

 if (is.null(theta.des)==TRUE) ptheta<-1;
 if (is.null(theta.des)==TRUE) theta.des<-matrix(1,antpers,ptheta) else 
 theta.des<-as.matrix(theta.des);
  ptheta<-ncol(theta.des);
 if (nrow(theta.des)!=antpers) stop("Theta design does not have correct dim");
  if (is.null(theta)==TRUE) theta<-rep(0.5,ptheta);
  if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta);
  Biid<-c()
  for (i in 1:antclust) Biid<-cbind(Biid,cif$B.iid[[i]]); 
  B2iid<-c()
  for (i in 1:antclust) B2iid<-cbind(B2iid,cif2$B.iid[[i]]); 
  if (npar==TRUE) gamma.iid<-0 else gamma.iid<-cif$gamma.iid;
  if (npar2==TRUE) gamma2.iid<-0 else gamma2.iid<-cif2$gamma.iid;
  var.theta<-hess<-matrix(0,ptheta,ptheta); score<-rep(0,ptheta); 
  time.pow<-attr(cif,"time.pow"); 
###  if (sum(time.pow)==0) time.pow<-rep(1,pg); 
## }}}

##dyn.load("random-cif.so");

  out<-.C("randomcif2cause", ## {{{
  as.double(times),as.integer(ntimes),as.double(time), as.integer(delta), as.integer(cause), as.integer(cause1),
  as.double(Gcx), as.double(X),as.integer(antpers), as.integer(px), as.integer(Nit), as.double(score),
  as.double(hess), as.double(est), as.double(gamma), as.integer(semi),as.double(Z), as.integer(pg),
  as.integer(detail), as.double(Biid),as.double(gamma.iid), as.double(time.pow), as.double(theta), as.double(var.theta), 
  as.double(theta.des), as.integer(ptheta), as.integer(antclust), as.integer(clusters), as.integer(clustsize), as.integer(clusterindex),
  as.integer(maxclust), as.double(step) , as.integer(inverse),as.integer(cause2) , as.double(X2),as.integer(px2),
  as.integer(semi2),as.double(Z2), as.integer(pg2), as.double(est2), as.double(gamma2),as.double(B2iid),
as.double(gamma2.iid),as.integer(sdscore),as.integer(squarepar),
as.integer(1-sym),as.integer(same.cens),PACKAGE="MultiComp") 
## }}}

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
class(ud)<-"randomcif"; 
attr(ud, "Call") <- sys.call()
attr(ud, "Formula") <- formula
attr(ud, "Clusters") <- clusters
attr(ud,"inverse")<-inverse
attr(ud,"cause1")<-cause1
attr(ud,"cause2")<-cause2
attr(ud,"antclust")<-antclust
attr(ud,"antpers")<-antpers
return(ud); 
## }}}

}


