dep.cif<-function(cif,data,cause,model="OR",cif2=NULL,times=NULL,
cause1=1,cause2=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL, theta=NULL,theta.des=NULL,parfunc=NULL,dparfunc=NULL,
step=1,sym=0,colnames=NULL,dimpar=NULL,weights=NULL,notaylor=1,
same.cens=FALSE,censoring.probs=NULL,silent=1,
entry=NULL,estimator=1,trunkp=1,admin.cens=NULL,...)
{ ## {{{
## {{{ set up data and design
 multi=0; dscore=1; stab.cens <- FALSE; entry.call <- entry
 formula<-attr(cif,"Formula"); 
 ldata<-aalen.des(formula,data,model="aalen");
 X<-ldata$X; time<-ldata$time2; Z<-ldata$Z;  status<-ldata$status;
 #print(dim(X)); print(dim(Z)); 
  antpers<-nrow(X); 
  if (is.null(Z)==TRUE) {npar<-TRUE; semi<-0;}  else {Z<-as.matrix(Z); npar<-FALSE; semi<-1;}
  if (npar==TRUE) {Z<-matrix(0,antpers,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-ncol(Z);} 
  ng<-antpers;px<-ncol(X);  
 if (is.null(times)) times<-cif$cum[,1]; 
 est<-cif$cum; if (semi==1) gamma<-cif$gamma  else gamma<-0; 
 ntimes<-length(times); 
 if ((cif$model!="additive") && (cif$model!="fg")) 
 stop("Marginal Cumulative incidence model, must be either additive or extended Fine-Gray model\n")
 cif.model <- switch(cif$model,additive=1,fg=2)
 if (missing(cause)) cause <- attr(cif,"cause"); 
 delta<-(cause!=cens.code)

 if (cause1!=attr(cif,"causeS")) cat("Cause for marginal model and correlation not the same\n"); 
if ((model!="COR") && (cause1[1]!=cause2[1])) {
 if (is.null(cif2)==TRUE) stop("Must provide marginal model for both causes"); 
 formula2<-attr(cif2,"Formula"); 
 ldata2<-aalen.des(formula2,data,model="aalen");
 X2<-ldata2$X; timec2<-ldata2$time2; Z2<-ldata$Z;  status2<-ldata2$status;
 if (is.null(Z2)==TRUE) {npar2<-TRUE; semi2<-0;}  else {Z2<-as.matrix(Z2); npar2<-FALSE; semi2<-1;}
 if (npar2==TRUE) {Z2<-matrix(0,antpers,1); pg2<-1; fixed2<-0;} else {fixed2<-1;pg2<-ncol(Z2);} 
 px2<-ncol(X2);  
 est2<-cif2$cum; if (semi2==1) gamma2<-cif2$gamma  else gamma2<-0; 
 est2<-Cpred(est2,times);
} else { 
X2<-0; Z2<-0; pg2<-1; px2<-1;  semi2<-0; est2<-0; gamma2<-0; 
npar2<-FALSE; est2 <- 0; 
}


### For truncation 
if (is.null(entry)) entry.call <- NULL else entry.call <- 0
if (is.null(entry)) entry <- rep(0,antpers); 
cum1<-Cpred(rbind(rep(0,px+1),cif$cum),entry)[,-1];
if (cif.model==1) cif1entry  <-  1-exp(-apply(X*cum1,1,sum)- (Z %*% gamma )*entry)
else if (cif.model==2) cif1entry  <-  1-exp(-apply(X*cum1,1,sum)*exp(Z %*% gamma ))
if ((model!="COR") && (cause1[1]!=cause2[1])) {
cum2<-Cpred(rbind(rep(0,px2+1),cif2$cum),entry)[,-1];
if (cif.model==1) cif2entry  <-  1-exp(-apply(X2*cum2,1,sum)- (Z2 %*% gamma2 )*entry)
else if (cif.model==2) cif2entry  <-  1-exp(-apply(X2*cum2,1,sum)*exp(Z2 %*% gamma2 ))
} else {cum2 <- cum1; cif2entry <- cif1entry;}
## }}}

## {{{ censoring model stuff
cens.weight <- cif$cens.weight ### censoring weights from cif function
Gcxe <- 1; 
if (cens.model!="user.weights") {
   if (is.null(cens.weight) ) {
   if (cens.model=="KM") { ## {{{
    ud.cens<-survfit(Surv(time,cause==cens.code)~+1); 
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    Gcx<-Cpred(Gfit,time)[,2];
    if (is.null(entry.call)==FALSE) Gcxe<-Cpred(Gfit,entry)[,2];
     Gcx <- Gcx/Gcxe; 
###  Gctimes<-Cpred(Gfit,times)[,2];
     Gctimes<- Gcx ## }}}
  } else if (cens.model=="cox") { ## {{{
    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
    ud.cens<-cox.aalen(Surv(time,cause==cens.code)~prop(XZ),n.sim=0,robust=0);
    Gcx<-Cpred(ud.cens$cum,time)[,2];
    RR<-exp(XZ %*% ud.cens$gamma)
    Gcx<-exp(-Gcx*RR)
    Gfit<-rbind(c(0,1),cbind(time,Gcx)); 
    if (is.null(entry.call)==FALSE) {
	    Gcxe<-Cpred(Gfit,entry)[,2];
            Gcxe<-exp(-Gcxe*RR)
    }
    Gcx <- Gcx/Gcxe; 
###  Gctimes<-Cpred(Gfit,times)[,2];
     Gctimes<- Gcx  ## }}}
  } else if (cens.model=="aalen") {  ## {{{
    if (npar==TRUE) XZ<-X[,-1] else XZ<-cbind(X,Z)[,-1];
    ud.cens<-aalen(Surv(time,cause==cens.code)~XZ,n.sim=0,robust=0);
    Gcx<-Cpred(ud.cens$cum,time)[,-1];
    XZ<-cbind(1,XZ); 
    Gcx<-exp(-apply(Gcx*XZ,1,sum))
    Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
    Gfit<-rbind(c(0,1),cbind(time,Gcx)); 
    if (is.null(entry.call)==FALSE) {
	    Gcxe<-Cpred(Gfit,entry)[,-1];
            Gcxe<-exp(-apply(Gcxe*XZ,1,sum))
            Gcxe[Gcxe>1]<-1; Gcxe[Gcxe<0]<-0
    }
    Gcx <- Gcx/Gcxe; 
###  Gctimes<-Cpred(Gfit,times)[,2];
     Gctimes<- Gcx  ## }}}
   }
   } else { Gcx <- cens.weight; Gctimes <- cens.weight;} 
} else {
    Gcx <- censoring.probs
    Gctimes <- cens.weight
} 
   ntimes<-length(times); 
## }}}

## {{{ set up cluster + theta design + define iid variables 
  if (is.null(clusters)== TRUE) {
	  ## take clusters from cif model 
	  clusters <- attr(cif,"clusters"); 
	  antclust<- length(unique(clusters)); 
  } else {
    clus<-unique(clusters); antclust<-length(clus);
    clusters <- as.integer(factor(clusters, labels = 1:(antclust)))-1;
  }

  outc <- cluster.index(clusters,index.type=TRUE); 
  clustsize <- outc$cluster.size
  maxclust <- outc$maxclust
  clusterindex <- outc$idclust
  if (maxclust==1) stop("No clusters given \n"); 

 if (is.null(theta.des)==TRUE) { ptheta<-1; theta.des<-matrix(1,antpers,ptheta);} else 
  theta.des<-as.matrix(theta.des); ptheta<-ncol(theta.des);
if ( (!is.null(parfunc)) && is.null(dimpar) ) 
   stop("Must specify dimension of score when specifying R-functions\n")
  if (is.null(dimpar) && is.null(parfunc)) dimpar<-ptheta; 
  ### if (!is.null(dparfunc)) dimpar<-length(dparfunc(

  if (is.null(theta)==TRUE) theta<-rep(0.0,dimpar);
  if (length(theta)!=dimpar) theta<-rep(theta[1],dimpar);
  ##print(est); print(gamma); print(semi); 
Biid<-c(); gamma.iid <- 0; B2iid<-c(); gamma2.iid <- 0; 
if (notaylor==0) {
for (i in 1:antclust) Biid<-cbind(Biid,cif$B.iid[[i]]); 
if (npar==TRUE) gamma.iid<-0 else gamma.iid<-cif$gamma.iid;
if ((model!="COR") && (cause1!=cause2))  {
B2iid<-c()
for (i in 1:antclust) B2iid<-cbind(B2iid,cif2$B.iid[[i]]); 
if (npar2==TRUE) gamma2.iid<-0 else  gamma2.iid<-cif2$gamma.iid;
} else { B2iid<-Biid; gamma2.iid<-gamma.iid; }
}
  var.theta<-hess<-matrix(0,dimpar,dimpar); score<-rep(0,dimpar); 
  time.pow<-attr(cif,"time.pow"); 
  #if (sum(time.pow)==0) time.pow<-rep(1,pg); 
  theta.iid<-matrix(0,antclust,dimpar); 

  if (nrow(theta.des)!=nrow(X)) stop("Dependence design not consistent with data\n"); 
 if (length(trunkp)!=antpers) trunkp <- rep(1,antpers)
if (is.null(weights)==TRUE) weights <- rep(1,antpers); 
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

print(clustsize)
print(clusters)
###dyn.load("cor.so");
if (model=="COR")
out<-.C("cor", ## {{{
  as.double(times),as.integer(ntimes),as.double(time), as.integer(delta), as.integer(cause), as.integer(cause1),
  as.double(Gcx), as.double(X),as.integer(antpers), as.integer(px), as.integer(Nit), as.double(score),
  as.double(hess), as.double(est), as.double(gamma), as.integer(semi),as.double(Z), as.integer(pg),
  as.integer(detail), as.double(Biid),as.double(gamma.iid), as.double(time.pow), as.double(theta), as.double(var.theta), 
  as.double(theta.des), as.integer(ptheta), as.integer(antclust), as.integer(clusters), as.integer(clustsize), as.integer(clusterindex),
  as.integer(maxclust), as.double(step) , as.integer(dscore),as.integer(cause2) , as.double(X2),as.integer(px2),
 as.integer(semi2),as.double(Z2), as.integer(pg2), as.double(est2), as.double(gamma2),as.double(B2iid),
 as.double(gamma2.iid),body(htheta),body(dhtheta),new.env(),as.integer(dimpar),as.integer(flex.func),
 as.double(theta.iid),as.integer(sym) , as.double(weights), as.integer(notaylor), 
 as.integer(same.cens),as.integer(stab.cens), as.double(Gctimes),
 as.integer(silent),as.integer(cif.model),PACKAGE="MultiComp") 
###as.integer(entryage),as.double(cifentry),as.double(trunkp),
## }}}
else if (model=="RR")
 out<-.C("mcifrr", ## {{{
  as.double(times),as.integer(ntimes),as.double(time), as.integer(delta), as.integer(cause), as.integer(cause1),
  as.double(Gcx), as.double(X),as.integer(antpers), as.integer(px), as.integer(Nit), as.double(score),
  as.double(hess), as.double(est), as.double(gamma), as.integer(semi),as.double(Z), as.integer(pg),
  as.integer(detail), as.double(Biid),as.double(gamma.iid), as.double(time.pow), as.double(theta), as.double(var.theta), 
  as.double(theta.des), as.integer(ptheta), as.integer(antclust), as.integer(clusters), as.integer(clustsize), as.integer(clusterindex),
  as.integer(maxclust), as.double(step) , as.integer(dscore),as.integer(cause2) , as.double(X2),as.integer(px2),
 as.integer(semi2),as.double(Z2), as.integer(pg2), as.double(est2), as.double(gamma2),as.double(B2iid),
 as.double(gamma2.iid), body(htheta),body(dhtheta),new.env(),as.integer(dimpar),as.integer(flex.func),
 as.double(theta.iid),as.integer(sym),as.double(weights), as.integer(notaylor),
 as.integer(same.cens),as.integer(estimator), 
 as.double(entry),as.double(cif1entry),as.double(cif2entry),as.double(trunkp),
 as.integer(silent),as.integer(cif.model),PACKAGE="MultiComp") 
## }}}
else if (model=="OR")
out<-.C("plackor", ## {{{
  as.double(times),as.integer(ntimes),as.double(time), as.integer(delta), as.integer(cause), as.integer(cause1),
  as.double(Gcx), as.double(X),as.integer(antpers), as.integer(px), as.integer(Nit), as.double(score),
  as.double(hess), as.double(est), as.double(gamma), as.integer(semi),as.double(Z), as.integer(pg),
  as.integer(detail), as.double(Biid),as.double(gamma.iid), as.double(time.pow), as.double(theta), as.double(var.theta), 
  as.double(theta.des), as.integer(ptheta), as.integer(antclust), as.integer(clusters), as.integer(clustsize), as.integer(clusterindex),
  as.integer(maxclust), as.double(step) , as.integer(dscore),as.integer(cause2) , as.double(X2),as.integer(px2),
  as.integer(semi2),as.double(Z2), as.integer(pg2), as.double(est2), as.double(gamma2),as.double(B2iid),
  as.double(gamma2.iid), body(htheta),body(dhtheta),new.env(),as.integer(dimpar),as.integer(flex.func),
  as.double(theta.iid),as.integer(sym) , as.double(weights), as.integer(notaylor),
  as.integer(same.cens),as.integer(estimator), 
  as.double(entry),as.double(cif1entry),as.double(cif2entry),as.double(trunkp),
  as.integer(silent),as.integer(cif.model),PACKAGE="MultiComp") 
## }}}

## {{{ output
theta<-matrix(out[[23]],dimpar,1);var.theta<-matrix(out[[24]],dimpar,dimpar); 
score<-matrix(out[[12]],dimpar,1); hess<-matrix(out[[13]],dimpar,dimpar); 
theta.iid<-matrix(out[[49]],antclust,dimpar); 

if (is.null(colnames)==FALSE) names<-colnames else { ## {{{
   if (is.null(parfunc)) names<-colnames(theta.des)
   else  {
   if (dimpar>1) names<-paste("par",1:dimpar,sep="") else names<-"intercept"
   }

  if (is.null(names)==FALSE && length(names)==dimpar) {
     rownames(score)<-rownames(theta)<-rownames(var.theta)<-
     colnames(var.theta)<-
     rownames(hess)<-colnames(hess)<- colnames(hess)<-names; 
     colnames(theta)<-"coefficient"; colnames(score)<-"score"; 
  }
} ## }}}

ud<-list(score=score,hess=hess,theta=theta,var.theta=var.theta,theta.iid=theta.iid)

class(ud)<-"cor"; 
attr(ud, "Call") <- sys.call()
###attr(ud, "Type") <- "cor"
attr(ud, "Formula") <- formula
attr(ud, "Clusters") <- clusters
attr(ud,"cause1")<-cause1; attr(ud,"cause2")<-cause2
attr(ud,"sym")<-sym; 
attr(ud,"antpers")<-antpers; 
attr(ud,"antclust")<-antclust; 
## }}}
if (model=="COR") attr(ud, "Type") <- "cor"
if (model=="RR") attr(ud, "Type") <- "RR"
if (model=="OR") attr(ud, "Type") <- "OR-cif"
return(ud); 
} ## }}}

cor.cif<-function(cif,data,cause,times=NULL,
cause1=1,cause2=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL, theta=NULL,theta.des=NULL,parfunc=NULL,dparfunc=NULL,
step=1,sym=0,colnames=NULL,dimpar=NULL,weights=NULL,notaylor=1,
same.cens=FALSE,censoring.probs=NULL,silent=1,...)
{ ## {{{
fit <- dep.cif(cif=cif,data=data,cause=cause,model="COR",cif2=NULL,times=times,
         cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
         clusters=clusters,theta=theta,theta.des=theta.des,parfunc=NULL,dparfunc=dparfunc,
         step=step,sym=sym,colnames=colnames,dimpar=dimpar,weights=weights,notaylor=notaylor,
         same.cens=same.cens,censoring.probs=censoring.probs,silent=silent,...)
    fit$call <- match.call()
    fit
} ## }}}

rr.cif<-function(cif,data,cause,cif2=NULL,times=NULL,
cause1=1,cause2=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL, theta=NULL,theta.des=NULL,parfunc=NULL,dparfunc=NULL,
step=1,sym=0,colnames=NULL,dimpar=NULL,weights=NULL,notaylor=1,
same.cens=FALSE,censoring.probs=NULL,silent=1,
entry=NULL,estimator=1,trunkp=1,admin.cens=NULL,...)
{ ## {{{
fit <- dep.cif(cif=cif,data=data,cause=cause,model="RR",cif2=cif2,times=times,
         cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
         clusters=clusters,theta=theta,theta.des=theta.des,parfunc=NULL,dparfunc=dparfunc,
         step=step,sym=sym,colnames=colnames,dimpar=dimpar,weights=weights,notaylor=notaylor,
         same.cens=same.cens,censoring.probs=censoring.probs,silent=silent,
	 entry=entry,estimator=estimator,trunkp=trunkp,admin.cens=admin.cens,...)
    fit$call <- match.call()
    fit
} ## }}}

or.cif<-function(cif,data,cause,cif2=NULL,times=NULL,
cause1=1,cause2=1,cens.code=0,cens.model="KM",Nit=40,detail=0,
clusters=NULL, theta=NULL,theta.des=NULL,parfunc=NULL,dparfunc=NULL,
step=1,sym=0,colnames=NULL,dimpar=NULL,weights=NULL,notaylor=1,
same.cens=FALSE,censoring.probs=NULL,silent=1,
entry=NULL,estimator=1,trunkp=1,admin.cens=NULL,...)
{ ## {{{
fit <- dep.cif(cif=cif,data=data,cause=cause,model="OR",cif2=cif2,times=times,
         cause1=cause1,cause2=cause2,cens.code=cens.code,cens.model=cens.model,Nit=Nit,detail=detail,
         clusters=clusters,theta=theta,theta.des=theta.des,parfunc=NULL,dparfunc=dparfunc,
         step=step,sym=sym,colnames=colnames,dimpar=dimpar,weights=weights,notaylor=notaylor,
         same.cens=same.cens,censoring.probs=censoring.probs,silent=silent,
	 entry=entry,estimator=estimator,trunkp=trunkp,admin.cens=admin.cens,...)
    fit$call <- match.call()
    fit
} ## }}}

print.summary.cor <- function(x,digits=3,...)
{ ## {{{
  if (x$type=="cor") {
    cat("Cross odds ratio dependence for competing risks\n\n")
    cat("Effect of cause1=",x$cause1," on cause2=",x$cause2,
        " under symmetry=",x$sym,fill=TRUE,sep="")
  } else if (x$type=="RR") {
    cat("Ratio of joint and product of marginals for competing risks\n\n")
    cat("Ratio of cumulative incidence for cause1=",x$cause1," and cause2=",x$cause2,sep=" ")
  } else if (x$type=="OR-cif") {
    cat("OR for dependence for competing risks\n\n")
    cat("OR of cumulative incidence for cause1=",x$cause1," and cause2=",x$cause2,sep=" ")
  }
  cat("\n")
  prmatrix(signif(x$estimates,digits))
  cat("\n")
  if (!is.null(x$marg)) {
    cat(paste("Marginal cumulative incidencen",signif(x$marg,digits),"\n"))
    prmatrix(signif(x$casewise,digits))
    prmatrix(signif(x$concordance,digits))
    cat("\n")
  }
  invisible(x)
} ## }}}

summary.cor<-function(object,digits=3,marg.cif=NULL,...)
{ ## {{{
  if (!inherits(object, "cor")) stop("Must be a cor.cif  object")
  if (sum(abs(object$score))>0.000001) warning("WARNING: check score for convergence\n")
  coefs <- coef.cor(object,...);

  outcase <- outconc <- NULL
  if (is.null(marg.cif)==FALSE) {
    marg.cif <- max(marg.cif)
    ## {{{
    if (attr(object,"Type")=="cor") {
      concordance <- exp(coefs[,1])*marg.cif^2/((1-marg.cif)+exp(coefs[,1])*marg.cif)
      conclower  <- exp(coefs[,1]-1.96*coefs[,2])*marg.cif^2/((1-marg.cif)+exp(coefs[,1]-1.96*coefs[,2])*marg.cif)
      concup  <- exp(coefs[,1]+1.96*coefs[,2])*marg.cif^2/((1-marg.cif)+exp(coefs[,1]+1.96*coefs[,2])*marg.cif)
      casewise <- concordance/marg.cif
      caselower  <- conclower/marg.cif
      caseup     <- concup/marg.cif
    } else if (attr(object,"Type")=="RR") {
      casewise<- exp(coefs[,1])*c(marg.cif)
      concordance <- exp(coefs[,1])*marg.cif^2
      caselower  <- marg.cif*exp(coefs[,1]-1.96*coefs[,2])
      caseup     <- marg.cif*exp(coefs[,1]+1.96*coefs[,2])
      conclower  <- marg.cif^2* exp(coefs[,1]-1.96*coefs[,2])
      concup     <- marg.cif^2*exp(coefs[,1]+1.96*coefs[,2])
    } else if (attr(object,"Type")=="OR-cif") {
      thetal <-  coefs[,1]-1.96*coefs[,2]
      thetau <-  coefs[,1]+1.96*coefs[,2]
      casewise<- plack.cif2(marg.cif,marg.cif,c(coefs[,1]))/marg.cif
      concordance <- plack.cif2(marg.cif,marg.cif,c(coefs[,1]))
      caselower  <- plack.cif2(marg.cif,marg.cif,thetal)/marg.cif
      caseup  <- plack.cif2(marg.cif,marg.cif,thetau)/marg.cif
      conclower  <- plack.cif2(marg.cif,marg.cif,thetal)
      concup     <- plack.cif2(marg.cif,marg.cif,thetau)
    }
    outcase <- cbind(casewise,caselower,caseup)
    outconc <- cbind(concordance,conclower,concup)
    rownames(outcase) <- rownames(outconc)  <-  rownames(coefs)
    colnames(outcase) <- c("casewise concordance","2.5 %","97.5%")
    colnames(outconc) <- c("concordance","2.5 %","97.5%")
  }
  ## }}}
  res <- list(casewise=outcase,concordance=outconc,estimates=coefs,marg=marg.cif,type=attr(object,"Type"),sym=attr(object,"sym"),cause1=attr(object,"cause1"),cause2=attr(object,"cause2"))
  class(res) <- "summary.cor"
  res
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

summary.cor<-function(object,marg.cif=NULL,marg.cif2=NULL,digits=3,...)
{ ## {{{
  if (!inherits(object, "cor")) stop("Must be a cor.cif  object")
  if (sum(abs(object$score))>0.000001) warning("WARNING: check score for convergence\n")
  coefs <- coef.cor(object,...);

  outcase <- outconc <- NULL
  if (is.null(marg.cif)==FALSE) {
      marg.cif <- max(marg.cif)
      if (attr(object,"cause2")==attr(object,"cause1")) marg.cif2=marg.cif 
      marg.cif2 <- max(marg.cif2)
      pmarg.cif <- marg.cif*marg.cif2
    ## {{{
    if (attr(object,"Type")=="cor") {
      concordance <- exp(coefs[,1])*pmarg.cif/((1-marg.cif)+exp(coefs[,1])*marg.cif)
      conclower  <- exp(coefs[,1]-1.96*coefs[,2])*pmarg.cif/((1-marg.cif)+exp(coefs[,1]-1.96*coefs[,2])*marg.cif)
      concup  <- exp(coefs[,1]+1.96*coefs[,2])*pmarg.cif/((1-marg.cif)+exp(coefs[,1]+1.96*coefs[,2])*marg.cif)
      casewise <- concordance/marg.cif
      caselower  <- conclower/marg.cif
      caseup     <- concup/marg.cif
    } else if (attr(object,"Type")=="RR") {
      casewise<- exp(coefs[,1])*c(marg.cif2)
      concordance <- exp(coefs[,1])*pmarg.cif
      caselower  <- marg.cif2*exp(coefs[,1]-1.96*coefs[,2])
      caseup     <- marg.cif2*exp(coefs[,1]+1.96*coefs[,2])
      conclower  <- pmarg.cif* exp(coefs[,1]-1.96*coefs[,2])
      concup     <- pmarg.cif*exp(coefs[,1]+1.96*coefs[,2])
    } else if (attr(object,"Type")=="OR-cif") {
      thetal <-  coefs[,1]-1.96*coefs[,2]
      thetau <-  coefs[,1]+1.96*coefs[,2]
      casewise<- plack.cif2(marg.cif,marg.cif,c(coefs[,1]))/marg.cif
      concordance <- plack.cif2(marg.cif,marg.cif,c(coefs[,1]))
      caselower  <- plack.cif2(marg.cif,marg.cif,thetal)/marg.cif
      caseup  <- plack.cif2(marg.cif,marg.cif,thetau)/marg.cif
      conclower  <- plack.cif2(marg.cif,marg.cif,thetal)
      concup     <- plack.cif2(marg.cif,marg.cif,thetau)
    }
    outcase <- cbind(casewise,caselower,caseup)
    outconc <- cbind(concordance,conclower,concup)
    rownames(outcase) <- rownames(outconc)  <-  rownames(coefs)
    colnames(outcase) <- c("casewise concordance","2.5 %","97.5%")
    colnames(outconc) <- c("concordance","2.5 %","97.5%")
  }
  ## }}}
  res <- list(casewise=outcase,concordance=outconc,estimates=coefs,
	      marg.cif=marg.cif, marg.cif2=marg.cif2,type=attr(object,"Type"),
	      sym=attr(object,"sym"),cause1=attr(object,"cause1"),cause2=attr(object,"cause2"))
  class(res) <- "summary.cor"
  res
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

concordance <- function(object,cif1,cif2=NULL,messages=TRUE,model=NULL,coefs=NULL,...)
{ ## {{{

if (is.null(model)) { 
if (!inherits(object, "cor")) stop("Must be a rr.cif, cor.cif or or.cif object")
model <- attr(object,"Type")
} 
if (is.null(coefs)) coefs <- coef(object)
if (is.null(coefs)) stop("Must give dependence parameters\n"); 

if (!is.null(object)) { 
	cause1 <-  attr(object,"cause1");
	cause2 <-  attr(object,"cause2"); 
} else cause1 <- cause2 <- 1

if (is.null(cif2)==TRUE) cif2 <- cif1; 

if (messages) {
if (model=="cor") { ## {{{
message("Cross odds ratio dependence for competing risks\n\n")
message("Odds of cause1=",cause1," given cause2=",cause2," relative to Odds of cause1=",cause1,"\n",fill=TRUE,sep="")
} else if (model=="RR") {
message("Ratio of joint and product of marginals for competing risks\n\n")
message("Ratio of cumulative incidence for cause1=",cause1," and cause2=",cause2,sep=" ")
} else if (model=="OR-cif") {
message("OR for dependence for competing risks\n\n")
message("OR of cumulative incidence for cause1=",cause1," and cause2=",cause2,sep=" ")
} ## }}}
}

out <- list()
for (k in 1:nrow(coefs)) {
## {{{
if (model=="cor") {
concordance <- exp(coefs[k,1])*cif1*cif2/((1-cif1)+exp(coefs[k,1])*cif1)
conclower  <- exp(coefs[k,1]-1.96*coefs[k,2])*cif1*cif2/((1-cif1)+exp(coefs[k,1]-1.96*coefs[k,2])*cif1)
concup  <- exp(coefs[k,1]+1.96*coefs[k,2])*cif1*cif2/((1-cif1)+exp(coefs[k,1]+1.96*coefs[k,2])*cif1)
casewise <- concordance/cif1
caselower  <- conclower/cif1
caseup     <- concup/cif1
} else if (model=="RR") {
casewise<- exp(coefs[k,1])*c(cif2)
concordance <- exp(coefs[k,1])*cif1*cif2
caselower  <- cif2*exp(coefs[k,1]-1.96*coefs[k,2])
caseup     <- cif2*exp(coefs[k,1]+1.96*coefs[k,2])
conclower  <- cif1*cif2* exp(coefs[k,1]-1.96*coefs[k,2])
concup     <- cif1*cif2*exp(coefs[k,1]+1.96*coefs[k,2])
} else if (model=="OR-cif") {
thetal <-  coefs[k,1]-1.96*coefs[k,2]
thetau <-  coefs[k,1]+1.96*coefs[k,2]
casewise<- plack.cif2(cif1,cif2,c(coefs[k,1]))/cif1
concordance <- plack.cif2(cif1,cif2,c(coefs[k,1]))
caselower  <- plack.cif2(cif1,cif2,thetal)/cif1
caseup  <- plack.cif2(cif1,cif2,thetau)/cif1
conclower  <- plack.cif2(cif1,cif2,thetal)
concup     <- plack.cif2(cif1,cif2,thetau)
}

outcase <- cbind(c(casewise),c(caselower),c(caseup))
outconc <- cbind(c(concordance),c(conclower),c(concup))
colnames(outcase) <- c("casewise concordance","2.5 %","97.5%")
colnames(outconc) <- c("concordance","2.5 %","97.5%")
## }}}

out[[k]] <- list(concordance=outconc,casewise.concordance=outcase)
names(out)[k] <- rownames(coefs)[k]
###k <- k+1
}

return(out)
} ## }}}

plack.cif <- function(cif1,cif2,object,X=1) 
{ ## {{{
coefs <- coef(object)
theta <- exp(object$theta); 
cif1 <- c(cif1); cif2 <- c(cif2)
cifs=cif1+cif2; 

valn=2*(theta-1); 
val1=(1+(theta-1)*(cifs))-( ((1+(theta-1)*cifs))^2-4*cif1*cif2*theta*(theta-1))^0.5; 
vali=cif1*cif2;
valr <- vali;
valr[valn!=0] <- val1/valn; 

valr <- matrix(valr,length(c(theta)),1)
rownames(valr)=colnames(coefs)
return(valr); 
} ## }}}

plack.cif2 <- function(cif1,cif2,theta,X=1) 
{ ## {{{
theta <- exp(c(theta))
cif1 <- c(cif1); cif2 <- c(cif2)
cifs=cif1+cif2; 

valn=2*(theta-1); 
val1=(1+(theta-1)*(cifs))-( ((1+(theta-1)*cifs))^2-4*cif1*cif2*theta*(theta-1))^0.5; 
vali=cif1*cif2;

valr <- vali;
valr[valn!=0] <- val1/valn; 

return(valr); 
} ## }}}
