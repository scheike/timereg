prop<-function(x) x

two.stage<-function(formula=formula(data),data=sys.parent(),
beta=0,Nit=60,detail=0,start.time=0,max.time=NULL,id=NULL, 
clusters=NULL, robust=1,
rate.sim=1,beta.fixed=0,theta=NULL,theta.des=NULL,var.link=0,step=1)
{
  ratesim<-rate.sim; inverse<-var.link
  call <- match.call()
  m <- match.call(expand=FALSE)
  m$robust<-m$start.time<-m$beta<-m$Nit<-m$detail<-m$max.time<-m$clusters<-m$rate.sim<-m$beta.fixed<-m$theta<-m$theta.des<-m$var.link<-m$step<-NULL

  if (robust==0) cat("When robust=0 no variance estimate\n"); 

  special <- c("prop","cluster")
  Terms <- if(missing(data)) terms(formula, special)
  else          terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms,model="cox.aalen")
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ

  if(is.null(clusters)) clusters <- des$clusters  
  
  pxz <- px + pz;

  survs<-read.surv(m,id,npar,clusters,start.time,max.time)
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; 
  time<-survs$start; time2<-survs$stop; status<-survs$status;
  ldata<-list(start=survs$start,stop=survs$stop,
              antpers=survs$antpers,antclust=survs$antclust);

  if (npar==FALSE) covar<-data.matrix(cbind(X,Z)) else 
  stop("Both multiplicative and additive model needed");

  Ntimes <- sum(status); 
  times<-c(start.time,time2[status==1]); times<-sort(times);
  if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
  times<-times[times<maxtimes]

  if ((sum(beta)==0) & (beta.fixed==0)) beta<-coxph(Surv(time,time2,status)~Z)$coef; 

  if (px==0) stop("No nonparametric terms (needs one!)");
  ud<-two.stageBase.reg(times,ldata,X,Z,
                        status,id,clusters,Nit=Nit,detail=detail,beta=beta,
                        robust=robust,ratesim=ratesim,namesX=covnamesX,
   namesZ=covnamesZ,beta.fixed=beta.fixed,theta=theta,theta.des=theta.des,
   inverse=var.link,step=step);

  if (px>0) {
    colnames(ud$cum)<-colnames(ud$var.cum)<- c("time",covnamesX)
    if (robust==1) colnames(ud$robvar.cum)<- c("time",covnamesX) }

  rownames(ud$gamma)<-c(covnamesZ); colnames(ud$gamma)<-"estimate"; 
  rownames(ud$score)<-c(covnamesZ); colnames(ud$score)<-"score"; 

  ptheta<-length(ud$theta); 
  if (ptheta>1) {
                rownames(ud$theta)<-colnames(theta.des);
                names(ud$theta.score)<-colnames(theta.des); } 
 else { names(ud$theta.score)<- rownames(ud$theta)<-"intercept" } 

 if (beta.fixed==1) {
    ud$gamma <- ud$var.gamma <- ud$robvar.gamma <- NULL
 }

  attr(ud,"Call")<-sys.call(); 
  class(ud)<-"two.stage"
  attr(ud,"Formula")<-formula;
  attr(ud,"id")<-id.call;
  attr(ud,"cluster")<-cluster.call;
  attr(ud,"start")<-start.time; 
  attr(ud,"time2")<-time2; 
  attr(ud,"var.link")<-var.link
  attr(ud,"beta.fixed")<-beta.fixed

  return(ud); 
}

two.stageBase.reg<-function (times, fdata, designX, designG, status,
id, clusters, Nit = 5, beta = 0, detail = 0, robust = 1, 
ratesim = 1, namesZ=NULL,namesX=NULL,beta.fixed=0,theta=NULL,
theta.des=NULL,inverse=0,step=1) 
{
    additive.resamp <-0; ridge <- 0; XligZ <- 0;
    Ntimes <- length(times)
    designX <- as.matrix(designX); designG <- as.matrix(designG)
    if (is.matrix(designX) == TRUE) px <- as.integer(dim(designX)[2])
    if (is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])
    if (is.matrix(designG) == TRUE) pg <- as.integer(dim(designG)[2])
    if (is.matrix(designG) == TRUE) ng <- as.integer(dim(designG)[1])
    if (nx != ng) print(" A design og B designs er ikke ens\n")

    cumint <- matrix(0, Ntimes, px + 1)
    vcum <- matrix(0, Ntimes, px + 1)
    Rvcu <- matrix(0, Ntimes, px + 1)
    if (sum(abs(beta)) == 0) betaS <- rep(0, pg) else betaS <- beta
    score <- betaS
    Varbeta <- matrix(0, pg, pg)
    Iinv <- matrix(0, pg, pg)
    RVarbeta <- matrix(0, pg, pg)
    if (is.null(theta.des)==TRUE) ptheta<-1; 
    if (is.null(theta.des)==TRUE) theta.des<-matrix(1,ng,ptheta) else
    theta.des<-as.matrix(theta.des); 
    ptheta<-ncol(theta.des); 
    if (is.null(theta)==TRUE) theta<-rep(0.1,ptheta); 
    if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
    theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

    cluster.size<-as.vector(table(clusters));
    maxclust<-max(cluster.size)
    idiclust<-matrix(0,fdata$antclust,maxclust); 
    for (i in 1:fdata$antclust) { idiclust[i,]<- which( clusters %in% (i-1))-1 }

    #dyn.load("two-stage-reg.so"); 

    nparout <- .C("twostagereg", 
        as.double(times), as.integer(Ntimes), 
        as.double(designX), as.integer(nx), as.integer(px), 
	as.double(designG), as.integer(ng), as.integer(pg), 
	as.integer(fdata$antpers),as.double(fdata$start),as.double(fdata$stop),
	as.double(betaS), as.integer(Nit), as.double(cumint), 
	as.double(vcum),  as.double(Iinv), as.double(Varbeta), 
	as.integer(detail), as.double(Rvcu), as.double(RVarbeta), 
         as.integer(id), as.integer(status), as.integer(ratesim), 
	as.double(score), as.integer(robust), as.integer(clusters),
        as.integer(fdata$antclust), as.integer(beta.fixed),
        as.double(theta),as.double(var.theta),as.double(theta.score),
        as.integer(inverse), as.integer(cluster.size), as.double(theta.des),
        as.integer(ptheta), as.double(Stheta),as.double(step),
        as.integer(idiclust),PACKAGE = "timereg")

    gamma <- matrix(nparout[[12]], pg, 1)
    cumint <- matrix(nparout[[14]], Ntimes, px + 1)
    vcum <- matrix(nparout[[15]], Ntimes, px + 1)
    Iinv <- matrix(nparout[[16]], pg, pg)
    Varbeta <- -matrix(nparout[[17]], pg, pg)
    Rvcu <- matrix(nparout[[19]], Ntimes, px + 1)
    RVarbeta <- -matrix(nparout[[20]], pg, pg)
    score <- matrix(nparout[[24]], pg, 1)


   theta<-matrix(nparout[[29]],ptheta,1);  
   var.theta<-matrix(nparout[[30]],ptheta,ptheta); 
   theta.score<-nparout[[31]]; 
   Stheta<-matrix(nparout[[32]],ptheta,ptheta); 

   ud <- list(cum = cumint, var.cum = vcum, robvar.cum = Rvcu, 
       gamma = gamma, var.gamma = Varbeta, robvar.gamma = RVarbeta, 
       D2linv = Iinv, score = score,  theta=theta,var.theta=var.theta,
       S.theta=Stheta,theta.score=theta.score)
   return(ud)
}

summary.two.stage<-function (object,digits = 3,...) {
  if (!(inherits(object, 'two.stage') )) stop("Must be a Two-Stage object")
  
  prop<-TRUE; 
  if (is.null(object$prop.odds)==TRUE) p.o<-FALSE else p.o<-TRUE
    
  var.link<-attr(object,"var.link");
  cat("Dependence parameter for Clayton-Oakes-Glidden  model\n"); 

  if (sum(abs(object$theta.score)>0.000001) ) 
    cat("Variance parameters did not converge, allow more iterations\n\n"); 

  ptheta<-nrow(object$theta)
  sdtheta<-diag(object$var.theta)^.5
  if (var.link==0) {
      vari<-object$theta
      sdvar<-diag(object$var.theta)^.5
  }
  else {
      vari<-exp(object$theta)
      sdvar<-vari*diag(object$var.theta)^.5
  }
  dep<-cbind(object$theta[,1],sdtheta)
  walddep<-object$theta[,1]/sdtheta; 
  waldpdep<-(1-pnorm(abs(walddep)))*2

  kendall<-1/(1+2/vari) 
  kendall.ll<-1/(1+2/(object$theta+1.96*sdvar)) 
  kendall.ul<-1/(1+2/(object$theta-1.96*sdvar)) 
  if (var.link==0) resdep<-signif(as.matrix(cbind(dep,walddep,waldpdep,kendall)),digits)
  else resdep<-signif(as.matrix(cbind(dep,walddep,waldpdep,vari,sdvar,kendall)),digits);

  if (var.link==0) colnames(resdep) <- c("Variance","SE","z","P-val","Kendall's tau") 
  else colnames(resdep)<-c("log(Variance)","SE","z","P-val","Variance","SE Var.",
                           "Kendall's tau")
  prmatrix(resdep); cat("   \n");  

  if (attr(object,"beta.fixed")==0) {
  cat("Marginal Cox-Aalen model fit\n\n"); 
  if (sum(abs(object$score)>0.000001) && sum(object$gamma)!=0) 
    cat("Marginal model did not converge, allow more iterations\n\n"); 

  if (prop) {
    if (p.o==FALSE) cat("Proportional Cox terms :  \n") else  cat("Covariate effects \n")

    coef.two.stage(object,digits=digits);
  }
  }
   cat("   \n");  cat("  Call: \n"); dput(attr(object, "Call")); cat("\n");
}

print.two.stage <- function (x,digits = 3,...) {
  if (!(inherits(x, 'two.stage') )) stop("Must be a Two-Stage object")
  cat(" Two-stage estimation for Clayton-Oakes-Glidden  model\n"); 
  cat(" Marginals of Cox-Aalen form, dependence by variance of Gamma distribution\n\n");  
  object <- x; rm(x);
  
  cat(" Nonparametric components : "); 
  cat(colnames(object$cum)[-1]); cat("   \n");  
  if (!is.null(object$gamma)) {
    cat(" Parametric components :  "); cat(rownames(object$gamma)); 
    cat("   \n");
  } 
  cat("   \n");  

  cat(" Call: \n");
  print(attr(object,'Call'))
}

coef.two.stage<-function(object,digits=3,d2logl=1,...) {
   coefBase(object,digits=digits,d2logl=d2logl,...)
}

plot.two.stage<-function(x,pointwise.ci=1,robust=0,specific.comps=FALSE,
		level=0.05, 
		start.time=0,stop.time=0,add.to.plot=FALSE,mains=TRUE,
                xlab="Time",ylab ="Cumulative regression function",...) {
  if (!(inherits(x, 'two.stage'))) stop("Must be a Two-Stage object")
  object <- x; rm(x);  
 
  B<-object$cum; V<-object$var.cum; p<-dim(B)[[2]]; 
  if (robust>=1) V<-object$robvar.cum; 

  if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1
  if (stop.time==0) stop.time<-max(B[,1]);

  med<-B[,1]<=stop.time & B[,1]>=start.time
  B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
  V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
  Vrob<-object$robvar.cum; 
  Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 

  c.alpha<- qnorm(1-level/2)
  for (v in comp) { 
    c.alpha<- qnorm(1-level/2)
    est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
    if (add.to.plot==FALSE) 
      {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
        if (mains==TRUE) title(main=colnames(B)[v]); }
    else lines(B[,1],est,type="s"); 
    if (pointwise.ci>=1) {
      lines(B[,1],ul,lty=pointwise.ci,type="s");
      lines(B[,1],nl,lty=pointwise.ci,type="s"); }
    if (robust>=1) {
      lines(B[,1],ul,lty=robust,type="s"); 
      lines(B[,1],nl,lty=robust,type="s"); }
    abline(h=0); 
}
}   


