cor.cif<-function(cif,data=sys.parent(),cause,times=NULL,parfunc=NULL,dparfunc=NULL,
cause1=1,cause2=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL,theta=NULL,theta.des=NULL,step=1,sym=1,colnames=NULL,dimpar=NULL)
{ ## {{{
## {{{ set up data and design
multi=0;inverse=0; cif2<-NULL
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
 if (is.null(times)) times<-cif$cum[,1]; 
 delta<-(cause!=cens.code)
 est<-cif$cum; if (semi==1) gamma<-cif$gamma  else gamma<-0; 
 ntimes<-length(times); 

if (cause1!=cause2 && 3==2) {
if (is.null(cif2)==TRUE) stop("Must provide marginal model for both causes"); 
 formula2<-attr(cif2,"Formula"); 
 ldata2<-aalen.des(formula2,data,model="aalen");
 X2<-ldata2$X; timec2<-ldata2$time2; Z2<-ldata$Z;  status2<-ldata2$status;
 if (is.null(Z2)==TRUE) {npar2<-TRUE; semi2<-0;}  else {Z2<-as.matrix(Z2); npar2<-FALSE; semi2<-1;}
 if (npar2==TRUE) {Z2<-matrix(0,antpers,1); pg2<-1; fixed2<-0;} else {fixed2<-1;pg2<-ncol(Z2);} 
 px2<-ncol(X2);  
 est2<-cif2$cum; if (semi2==1) gamma2<-cif2$gamma  else gamma2<-0; 
 est<-Cpred(est,times); 
 est2<-Cpred(est2,times);
####      print(est); print(est2)
if (ncol(est)!=ncol(est2)) stop("Marginal models must be computed in same
time points \n"); 
} else { 
X2<-0; Z2<-0; pg2<-1; px2<-1;  semi2<-0; est2<-0; gamma2<-0; 
npar2<-FALSE
}
## }}}

## {{{ censoring model stuff
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
## }}}

## {{{ set up cluster + theta design + define iid variables 
  if (is.null(clusters)== TRUE) {clusters<-0:(antpers-1); antclust<-antpers;} else {
    clus<-unique(clusters); antclust<-length(clus);
    clusters <- as.integer(factor(clusters, labels = 0:(antclust-1)));}
  clustsize<-as.vector(table(clusters));
  maxclust<-max(clustsize); clusterindex<-matrix(0,antclust,maxclust); 
  if (antclust!=antpers) {
  for (i in 1:antclust) { index<-(((1:antpers)[clusters==i])-1); 
  clusterindex[i,1:length(index)]<-(((1:antpers)[clusters==i])-1) }
  } else clusterindex<-(0:(antpers-1)); 

 if (is.null(theta.des)==TRUE) { ptheta<-1; theta.des<-matrix(1,antpers,ptheta);} else 
  theta.des<-as.matrix(theta.des); ptheta<-ncol(theta.des);
if ( (!is.null(parfunc)) && is.null(dimpar) ) 
   stop("Must specify dimension of score when specifying R-functions\n")
  if (is.null(dimpar) && is.null(parfunc)) dimpar<-ptheta; 
  ### if (!is.null(dparfunc)) dimpar<-length(dparfunc(

  if (is.null(theta)==TRUE) theta<-rep(0.0,dimpar);
  if (length(theta)!=dimpar) theta<-rep(theta[1],dimpar);
  ##print(est); print(gamma); print(semi); 
Biid<-c()
for (i in 1:antclust) Biid<-cbind(Biid,cif$B.iid[[i]]); 
if (npar==TRUE) gamma.iid<-0 else gamma.iid<-cif$gamma.iid;
if (cause1!=cause2)  {
B2iid<-c()
for (i in 1:antclust) B2iid<-cbind(B2iid,cif2$B.iid[[i]]); 
if (npar2==TRUE) gamma2.iid<-0 else  gamma2.iid<-cif2$gamma.iid;
} else { B2iid<-Biid; gamma2.iid<-gamma.iid; }
  var.theta<-hess<-matrix(0,dimpar,dimpar); score<-rep(0,dimpar); 
  time.pow<-attr(cif,"time.pow"); 
  #if (sum(time.pow)==0) time.pow<-rep(1,pg); 
  theta.iid<-matrix(0,antclust,dimpar); 
## }}}

## {{{ function and derivative
if (is.null(parfunc)==FALSE)
{
flex.func<-1; # use flexible design
htheta<- function(theta,t,x) {
    out <- parfunc(theta,t,x)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
dhtheta<- function(theta,t,x) {
    out <- dparfunc(theta,t,x)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
}
else { 
   htheta<-function() {return(1)}; dhtheta<-function() {return(1)}; 
   dimpar<-ptheta; flex.func<-0; 
} 
## }}}

###dyn.load("cor.so");

  out<-.C("cor", ## {{{
  as.double(times),as.integer(ntimes),as.double(time), as.integer(delta), as.integer(cause), as.integer(cause1),
  as.double(Gcx), as.double(X),as.integer(antpers), as.integer(px), as.integer(Nit), as.double(score),
  as.double(hess), as.double(est), as.double(gamma), as.integer(semi),as.double(Z), as.integer(pg),
  as.integer(detail), as.double(Biid),as.double(gamma.iid), as.double(time.pow), as.double(theta), as.double(var.theta), 
  as.double(theta.des), as.integer(ptheta), as.integer(antclust), as.integer(clusters), as.integer(clustsize), as.integer(clusterindex),
  as.integer(maxclust), as.double(step) , as.integer(inverse),as.integer(cause2) , as.double(X2),as.integer(px2),
as.integer(semi2),as.double(Z2), as.integer(pg2), as.double(est2), as.double(gamma2),as.double(B2iid),
as.double(gamma2.iid), body(htheta),body(dhtheta),new.env(),as.integer(dimpar),as.integer(flex.func),
as.double(theta.iid),as.integer(sym) , PACKAGE="MultiComp")
## }}}

## {{{ output
theta<-matrix(out[[23]],dimpar,1);var.theta<-matrix(out[[24]],dimpar,dimpar); 
score<-matrix(out[[12]],dimpar,1); hess<-matrix(out[[13]],dimpar,dimpar); 
theta.iid<-matrix(out[[49]],antclust,dimpar); 

if (is.null(colnames)==FALSE) names<-colnames else {
   if (is.null(parfunc)) names<-colnames(theta.des)
   else  {
   if (dimpar>1) names<-paste("par",1:dimpar,sep="") else names<-"intercept"
   }
}

if (is.null(names)==FALSE && length(names)==dimpar) {
   rownames(score)<-rownames(theta)<-rownames(var.theta)<-colnames(var.theta)<-
   rownames(hess)<-colnames(hess)<- colnames(hess)<-names; 
   colnames(theta)<-"coefficient"; colnames(score)<-"score"; 
}

ud<-list(score=score,hess=hess,theta=theta,var.theta=var.theta,theta.iid=theta.iid)

class(ud)<-"cor"; 
attr(ud, "Call") <- sys.call()
attr(ud, "Formula") <- formula
attr(ud, "Clusters") <- clusters
attr(ud, "Type") <- "cor"
attr(ud,"cause1")<-cause1; attr(ud,"cause2")<-cause2
attr(ud,"sym")<-sym; 
attr(ud,"antpers")<-antpers; 
attr(ud,"antclust")<-antclust; 
return(ud); 
## }}}
} ## }}}

rr.cif<-function(cif,data=sys.parent(),cause,times=NULL,
cause1=1,cause2=1,cens.code=0,cif2=NULL,cens.model="KM",Nit=40,detail=0,
clusters=NULL,theta=NULL,theta.des=NULL,parfunc=NULL,dparfunc=NULL,
step=1,sym=1,colnames=NULL,dimpar=NULL)
{ ## {{{
## {{{ set up data and design
multi=0;inverse=0; 
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
 if (is.null(times)) times<-cif$cum[,1]; 
 delta<-(cause!=cens.code)
 est<-cif$cum; if (semi==1) gamma<-cif$gamma  else gamma<-0; 
 ntimes<-length(times); 

if (cause1!=cause2) {
if (is.null(cif2)==TRUE) stop("Must provide marginal model for both causes"); 
 formula2<-attr(cif2,"Formula"); 
 ldata2<-aalen.des(formula2,data,model="aalen");
 X2<-ldata2$X; timec2<-ldata2$time2; Z2<-ldata$Z;  status2<-ldata2$status;
 if (is.null(Z2)==TRUE) {npar2<-TRUE; semi2<-0;}  else {Z2<-as.matrix(Z2); npar2<-FALSE; semi2<-1;}
 if (npar2==TRUE) {Z2<-matrix(0,antpers,1); pg2<-1; fixed2<-0;} else {fixed2<-1;pg2<-ncol(Z2);} 
 px2<-ncol(X2);  
 est2<-cif2$cum; if (semi2==1) gamma2<-cif2$gamma  else gamma2<-0; 
 est<-Cpred(est,times); 
 est2<-Cpred(est2,times);
####      print(est); print(est2)
if (ncol(est)!=ncol(est2)) stop("Marginal models must be computed in same
time points \n"); 
} else { 
X2<-0; Z2<-0; pg2<-1; px2<-1;  semi2<-0; est2<-0; gamma2<-0; 
npar2<-FALSE
}
## }}}

## {{{ censoring model stuff
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
## }}}

## {{{ set up cluster + theta design + define iid variables 
  if (is.null(clusters)== TRUE) {clusters<-0:(antpers-1); antclust<-antpers;} else {
    clus<-unique(clusters); antclust<-length(clus);
    clusters <- as.integer(factor(clusters, labels = 0:(antclust-1)));}
  clustsize<-as.vector(table(clusters));
  maxclust<-max(clustsize); clusterindex<-matrix(0,antclust,maxclust); 
  if (antclust!=antpers) {
  for (i in 1:antclust) { index<-(((1:antpers)[clusters==i])-1); 
  clusterindex[i,1:length(index)]<-(((1:antpers)[clusters==i])-1) }
  } else clusterindex<-(0:(antpers-1)); 

 if (is.null(theta.des)==TRUE) { ptheta<-1; theta.des<-matrix(1,antpers,ptheta);} else 
  theta.des<-as.matrix(theta.des); ptheta<-ncol(theta.des);
if ( (!is.null(parfunc)) && is.null(dimpar) ) 
   stop("Must specify dimension of score when specifying R-functions\n")
  if (is.null(dimpar) && is.null(parfunc)) dimpar<-ptheta; 
  ### if (!is.null(dparfunc)) dimpar<-length(dparfunc(

  if (is.null(theta)==TRUE) theta<-rep(0.0,dimpar);
  if (length(theta)!=dimpar) theta<-rep(theta[1],dimpar);
  ##print(est); print(gamma); print(semi); 
Biid<-c()
for (i in 1:antclust) Biid<-cbind(Biid,cif$B.iid[[i]]); 
if (npar==TRUE) gamma.iid<-0 else gamma.iid<-cif$gamma.iid;
if (cause1!=cause2)  {
B2iid<-c()
for (i in 1:antclust) B2iid<-cbind(B2iid,cif2$B.iid[[i]]); 
if (npar2==TRUE) gamma2.iid<-0 else  gamma2.iid<-cif2$gamma.iid;
} else { B2iid<-Biid; gamma2.iid<-gamma.iid; }
  var.theta<-hess<-matrix(0,dimpar,dimpar); score<-rep(0,dimpar); 
  time.pow<-attr(cif,"time.pow"); 
  #if (sum(time.pow)==0) time.pow<-rep(1,pg); 
  theta.iid<-matrix(0,antclust,dimpar); 
## }}}

## {{{ function and derivative
if (is.null(parfunc)==FALSE)
{
flex.func<-1; # use flexible design
htheta<- function(theta,t,x) {
    out <- parfunc(theta,t,x)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
dhtheta<- function(theta,t,x) {
    out <- dparfunc(theta,t,x)
    if(!is.numeric(out)) stop("Need a numeric result")
    as.double(out)
}
}
else { 
   htheta<-function() {return(1)}; dhtheta<-function() {return(1)}; 
   dimpar<-ptheta; flex.func<-0; 
} 
## }}}

###dyn.load("cor.so");

  out<-.C("mcifrr", ## {{{
  as.double(times),as.integer(ntimes),as.double(time), as.integer(delta), as.integer(cause), as.integer(cause1),
  as.double(Gcx), as.double(X),as.integer(antpers), as.integer(px), as.integer(Nit), as.double(score),
  as.double(hess), as.double(est), as.double(gamma), as.integer(semi),as.double(Z), as.integer(pg),
  as.integer(detail), as.double(Biid),as.double(gamma.iid), as.double(time.pow), as.double(theta), as.double(var.theta), 
  as.double(theta.des), as.integer(ptheta), as.integer(antclust), as.integer(clusters), as.integer(clustsize), as.integer(clusterindex),
  as.integer(maxclust), as.double(step) , as.integer(inverse),as.integer(cause2) , as.double(X2),as.integer(px2),
as.integer(semi2),as.double(Z2), as.integer(pg2), as.double(est2), as.double(gamma2),as.double(B2iid),
as.double(gamma2.iid), body(htheta),body(dhtheta),new.env(),as.integer(dimpar),as.integer(flex.func),
as.double(theta.iid),as.integer(sym) , PACKAGE="MultiComp")
## }}}

## {{{ output
theta<-matrix(out[[23]],dimpar,1);var.theta<-matrix(out[[24]],dimpar,dimpar); 
score<-matrix(out[[12]],dimpar,1); hess<-matrix(out[[13]],dimpar,dimpar); 
theta.iid<-matrix(out[[49]],antclust,dimpar); 

if (is.null(colnames)==FALSE) names<-colnames else {
   if (is.null(parfunc)) names<-colnames(theta.des)
   else  {
   if (dimpar>1) names<-paste("par",1:dimpar,sep="") else names<-"intercept"
   }
}

if (is.null(names)==FALSE && length(names)==dimpar) {
   rownames(score)<-rownames(theta)<-rownames(var.theta)<-colnames(var.theta)<-
   rownames(hess)<-colnames(hess)<- colnames(hess)<-names; 
   colnames(theta)<-"coefficient"; colnames(score)<-"score"; 
}

ud<-list(score=score,hess=hess,theta=theta,var.theta=var.theta,theta.iid=theta.iid)

class(ud)<-"cor"; 
attr(ud, "Call") <- sys.call()
attr(ud, "Type") <- "RR"
attr(ud, "Formula") <- formula
attr(ud, "Clusters") <- clusters
attr(ud,"cause1")<-cause1; attr(ud,"cause2")<-cause2
attr(ud,"sym")<-sym; 
attr(ud,"antpers")<-antpers; 
attr(ud,"antclust")<-antclust; 
return(ud); 
## }}}
} ## }}}

summary.cor<-function(object,digits=3,...)
{ ## {{{
if (!inherits(object, "cor")) stop("Must be a cor.cif  object")

if (attr(object,"Type")=="cor") {
cat("Cross odds ratio dependence for competing risks\n\n")
cat("Effect of cause2=",attr(object,"cause2")," on cause1=",attr(object,"cause1"),
" under symmetry=",attr(object,"sym"),fill=TRUE,sep="")
} else {
cat("Ratio of joint and product of marginals for competing risks\n\n")
cat("Ratio of cumulative incidence for cause1=",attr(object,"cause1"),
		                   " and cause2=",attr(object,"cause2"),sep="")
}
cat("\n")
if (sum(abs(object$score))>0.000001) cat("WARNING: check score for convergence")
cat("\n")
prmatrix(signif(coef.cor(object,...),digits))
} ## }}}

coef.cor<-function(object,...)
{ ## {{{
 res <- cbind(object$theta, diag(object$var.theta)^0.5)
 se<-diag(object$var.theta)^0.5
 wald <- object$theta/se
 waldp <- (1 - pnorm(abs(wald))) * 2
 cor<-exp(object$theta)
 res <- as.matrix(cbind(res, wald, waldp,cor,se*cor))
if (attr(object,"Type")=="cor") 
 colnames(res) <- c("log-Coef.", "SE", "z", "P-val","Cross odds ratio","SE")
else colnames(res) <- c("log-ratio Coef.", "SE", "z", "P-val","Ratio","SE")
 if (is.null((rownames(res)))==TRUE) rownames(res)<-rep(" ",nrow(res))

return(res)
} ## }}}

print.cor<-function(x,digits=3,...)
{ ## {{{
print(attr(x,"Call")); 
cat("\n\n")
summary(x); 
} ## }}}

