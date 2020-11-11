#' Identifies proportional excess terms of model
#' 
#' Specifies which of the regressors that lead to proportional excess hazard
#' 
#' @param x variable
#' 
#' @author Thomas Scheike
#' @keywords survival
cox<-function(x) x



#' Fits Proportional excess hazards model
#' 
#' Fits proportional excess hazards model.
#' 
#' The models are written using the survival modelling given in the survival
#' package.
#' 
#' The program assumes that there are no ties, and if such are present random
#' noise is added to break the ties.
#' 
#' @param formula a formula object, with the response on the left of a `~'
#' operator, and the terms on the right.  The response must be a survival
#' object as returned by the `Surv' function.
#' @param data a data.frame with the variables.
#' @param excess specifies for which of the subjects the excess term is
#' present. Default is that the term is present for all subjects.
#' @param tol tolerance for numerical procedure.
#' @param max.time stopping considered time-period if different from 0.
#' Estimates thus computed from [0,max.time] if max.time>0. Default is max of
#' data.
#' @param n.sim number of simulations in re-sampling.
#' @param alpha tuning paramter in Newton-Raphson procedure. Value smaller than
#' one may give more stable convergence.
#' @param frac number between 0 and 1.  Is used in supremum test where observed
#' jump times t1, ..., tk is replaced by t1, ..., tl with l=round(frac*k).
#' @return Returns an object of type "prop.excess". With the following
#' arguments: \item{cum}{estimated cumulative regression functions. First
#' column contains the jump times, then follows the estimated components of
#' additive part of model and finally the excess cumulative baseline. }
#' \item{var.cum}{robust pointwise variance estimates for estimated
#' cumulatives. } \item{gamma}{estimate of parametric components of model. }
#' \item{var.gamma}{robust variance estimate for gamma. } \item{pval}{p-value
#' of Kolmogorov-Smirnov test (variance weighted) for excess baseline and Aalen
#' terms, H: B(t)=0. } \item{pval.HW}{p-value of supremum test (corresponding
#' to Hall-Wellner band) for excess baseline and Aalen terms, H: B(t)=0.
#' Reported in summary. } \item{pval.CM}{p-value of Cramer von Mises test for
#' excess baseline and Aalen terms, H: B(t)=0. } \item{quant}{95 percent
#' quantile in distribution of resampled Kolmogorov-Smirnov test statistics for
#' excess baseline and Aalen terms. Used to construct 95 percent simulation
#' band.  } \item{quant95HW}{95 percent quantile in distribution of resampled
#' supremum test statistics corresponding to Hall-Wellner band for excess
#' baseline and Aalen terms. Used to construct 95 percent Hall-Wellner band.  }
#' \item{simScoreProp}{observed scoreprocess and 50 resampled scoreprocesses
#' (under model). List with 51 elements. }
#' @author Torben Martinussen
#' @references Martinussen and Scheike, Dynamic Regression Models for Survival
#' Data, Springer Verlag (2006).
#' @keywords survival
#' @examples
#' 
#' ###working on memory leak issue, 3/3-2015
#' ###data(melanoma)
#' ###lt<-log(melanoma$thick)          # log-thickness 
#' ###excess<-(melanoma$thick>=210)    # excess risk for thick tumors
#' ###
#' #### Fits Proportional Excess hazards model 
#' ###fit<-prop.excess(Surv(days/365,status==1)~sex+ulc+cox(sex)+
#' ###                 cox(ulc)+cox(lt),melanoma,excess=excess,n.sim=100)
#' ###summary(fit)
#' ###par(mfrow=c(2,3))
#' ###plot(fit)
#' 
#' @export
prop.excess<-function(formula=formula(data),data=parent.frame(),excess=1,
tol=0.0001,max.time=NULL,n.sim=1000,alpha=1,frac=1)
  {
    id<-NULL;detail<-0
    robust<-0;  
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
    m$detail<-m$excess<-m$alpha<-m$frac<-m$id<-m$tol<-m$max.time<-
      m$n.sim<-NULL

    special <- c("cox")
    Terms <- if(missing(data)) terms(formula, special)
    else              terms(formula, special, data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    mt <- attr(m, "terms")
    intercept<-attr(mt, "intercept")
    Y <- model.extract(m, "response")

    des<-read.design(m,Terms,model="prop.excess")
    X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
    covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
    pxz <- px + pz;

    start.time<-0; clusters<-NULL; 
    survs<-read.surv(m,id,npar,clusters,start.time,max.time,model="prop-exs")
    times<-survs$times;id<-id.call<-survs$id.cal;
    clusters<-cluster.call<-survs$clusters; time2<-survs$stop
    status<-survs$status;
    ldata<-list(start=survs$start,stop=survs$stop,
                antpers=survs$antpers,antclust=survs$antclust);

    if (npar==FALSE) {covar<-cbind(X,Z);}
    else {stop("Both multiplicative and additive model needed");}

    Ntimes <- sum(status); 

    times<-c(0,time2[status==1]); times<-sort(times);
    if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
    times<-times[times<maxtimes] 
    status[time2>maxtimes]<-0

    #cat(" Proportional Excess Survival Model "); cat("\n")
    if (px==0) stop("No Aalen terms (need one!)");

    ud<-prop.excessBase(time2,status,X,Z,excess,alpha=1,frac=1,no.sim=n.sim)

    colnames(ud$cum)<-colnames(ud$var.cum)<- 
      c("time",covnamesX,"Excess baseline")
    rownames(ud$gamma)<-c(covnamesZ); 
    colnames(ud$gamma)<-"estimate"; 
    namematrix(ud$var.gamma,covnamesZ); 
                                        #namematrix(ud$robvar.gamma,covnamesZ); 
                                        #namematrix(ud$D2linv,covnamesZ); 

    attr(ud,"Call")<-call; 
    ud$call<-call

    class(ud)<-"prop.excess"
    return(ud); 
  }


#' @export
"plot.prop.excess" <- function(x, 
pointwise.ci=1, hw.ci=0,sim.ci=0,specific.comps=FALSE,level=0.95,start.time = 0, stop.time = 0, add.to.plot=FALSE, mains=TRUE, 
xlab="Time", ylab ="Cumulative regression function",score=FALSE,...) 
{
  prop.excess.object <- x; rm(x);
                                        # pointwise.ci   Pointwise confidence intervals 
                                        # hw.ci          Hall-Wellner confidence bands 95 \% 
                                        # robust         Robust variance used #not included here 
                                        # specific.comps Plots specific components c(2,3,4), e.g.
                                        # signi          Significance level
  signi<-1-level; 

  if (!inherits(prop.excess.object,'prop.excess') ) 
    stop ("Must be output from Proportional Excess Survival Model function") 
                                        #if (score==TRUE) {cat("Do not plot score processes"); score<-FALSE;}
 
  if (score==FALSE) {
    B<-prop.excess.object$cum; V<-prop.excess.object$var.cum; 
    p<-dim(B)[[2]]; 
                                        # if (robust==1) V<-prop.excess.object$robvar.cum; 

    if (sum(specific.comps)==FALSE) comp<-2:p else comp<-specific.comps+1

    if (stop.time==0) stop.time<-max(B[,1]);

    med<-B[,1]<=stop.time & B[,1]>=start.time
    B<-B[med,]; V<-V[med,]
    Bs<-B[1,];  Vs<-V[1,]
    B<-t(t(B)-Bs); V<-t( t(V)-Vs); 
    B[,1]<-B[,1]+Bs[1]

    c.alpha<- qnorm(1-signi/2)

    for (v in comp) { 
      ul<-B[,v]+c.alpha*V[,v]^.5; nl<-B[,v]-c.alpha*V[,v]^.5
      est<-B[,v]; 
      if (add.to.plot==FALSE) 
        {
          plot(B[,1],est,ylim=range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
          if (mains==TRUE) title(main=colnames(B)[v]); }
      else lines(B[,1],B[,v],type="s"); 
      if (pointwise.ci>=1) {
        lines(B[,1],ul,lty=pointwise.ci,type="s"); 
        lines(B[,1],nl,lty=pointwise.ci,type="s"); }
                                        #if (robust>=1) {
                                        #lines(B[,1],ul,lty=robust,type="s"); 
                                        #lines(B[,1],nl,lty=robust,type="s"); }
      if (hw.ci>=1) {
        if (signi!=0.05) cat("Hall-Wellner band only 95% \n"); 
        tau<-length(B[,1])
        nl<-B[,v]-prop.excess.object$quant95HW[v-1]*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+prop.excess.object$quant95HW[v-1]*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="s"); lines(B[,1],nl,lty=hw.ci,type="s"); }
      if (sim.ci>=1) {
        if (signi!=0.05) cat("Simulation based band only 95% \n");                                        
        V<-prop.excess.object$var.cum;                             
        nl<-B[,v]-prop.excess.object$quant95[v-1]*V[,v]^.5; 
        ul<-B[,v]+prop.excess.object$quant95[v-1]*V[,v]^.5;        
        lines(B[,1],ul,lty=sim.ci); lines(B[,1],nl,lty=sim.ci,type="s"); }

      abline(h=0)
    }
  } else {
                                        # plot score process
    dim1<-ncol(prop.excess.object$simScoreProp[[1]])
    if (sum(specific.comps)==FALSE) comp<-1:dim1 else comp<-specific.comps

    for (i in comp)
      {
        ranyl<-range(as.matrix(prop.excess.object$simScoreProp[[1]])[,i]);
        for (j in 2:51)
          ranyl<-range(c(ranyl,
                         as.matrix(prop.excess.object$simScoreProp[[j]])[,i]));
        mr<-max(abs(ranyl));

        plot(c(0,prop.excess.object$cum[,1]),
             c(0,as.matrix(prop.excess.object$simScoreProp[[1]])[,i]),type="s",
             ylim=c(-mr,mr),lwd=2,xlab=xlab,ylab="Scoreprocess")
        if (mains==TRUE) title(main=rownames(prop.excess.object$gamma)[i]); 
        for (j in 2:51)
          lines(c(0,prop.excess.object$cum[,1]),
                c(0,as.matrix(prop.excess.object$simScoreProp[[j]])[,i]),
                col="grey",lwd=1,lty=1,type="s")
        lines(c(0,prop.excess.object$cum[,1]),
              c(0,as.matrix(prop.excess.object$simScoreProp[[1]])[,i]),lwd=2,type="s")
      } }
}

#' @export
"print.prop.excess" <- function (x,...) {
  prop.excess.object <- x; rm(x);
  if (!inherits(prop.excess.object, 'prop.excess')) 
    stop ("Must be a Proportional Excess Survival Model object")

  if (is.null(prop.excess.object$gamma)==TRUE) cox<-FALSE else cox<-TRUE
    
                                        # We print information about object:  
  cat("Proportional Excess Survival Model \n\n")
  cat("Additive Aalen terms: "); 
  cat(colnames(prop.excess.object$cum)[-1]);
  cat("   \n");  
  if (cox) {
    cat("Proportional terms:  "); 
    cat(rownames(prop.excess.object$gamma)); 
    cat("   \n");  }
  cat("   \n");  

  cat("  Call: ")
  dput(attr(prop.excess.object, "Call"))
  cat("\n")
}


#' @export
"summary.prop.excess" <-
function (object,digits=3,...) 
{
    prop.excess.object <- object; rm(object);
  obj<-prop.excess.object
  if (!inherits(prop.excess.object, 'prop.excess')) 
    stop ("Must be a Proportional Excess Survival Model object")
  
  cox<-TRUE; 
  if (is.null(prop.excess.object$gamma)==TRUE) stop(" No proportional  terms"); 
    
                                        # We print information about object:  
  cat("Proportional Excess Survival Model \n\n")

                                        #if (sum(obj$conf.band)==FALSE)  mtest<-FALSE else mtest<-TRUE; 
                                        #if (mtest==FALSE) cat("Test not computed, sim=0 \n\n")
                                        #if (mtest==TRUE) { 
  test0<-cbind(obj$pval.HW,obj$pval.CM)
                                        # testC<-cbind(obj$obs.testBeqC,obj$pval.testBeqC) 
  colnames(test0)<- c("KS-test pval",
                      "CM-test pval")
  rownames(test0)<-colnames(prop.excess.object$cum)[-1]
                                        #rownames(test0)<- c("","p-value")
                                        #colnames(testC)<- c("sup| B(t) - (t/tau)B(tau)|","p-value H_0: B(t)=b t")  
  cat("Test for non-significant effects \n");cat("   \n");
  cat("Test for Aalen terms, H_0: B(t)=0  \n")
  prmatrix(signif(test0,digits));cat("   \n");  
                                        #cat("Test for time invariant effects \n")
                                        #prmatrix(signif(testC,digits))
                                        #cat("\n")
                                        #}

  if (cox) {
    cat("Proportional terms:  \n"); 

    res <- cbind(obj$gamma,
                 diag(obj$var.gamma)^.5) #,diag(obj$robvar.gamma)^.5,diag(obj$D2linv)^.5)
    z<-c((res[,1]/res[,2]))
    pval<-1-pchisq(z^2,1)
    res<-as.matrix(cbind(res,z,pval)); 
    colnames(res) <- c("coef", "se(coef)","z","p")
                                        #,#Robust Std.  Error","D2log(L)^-(1/2)")  
    prmatrix(signif(res, digits)); cat("   \n");  
  }

  cat("  Call: ")
  dput(attr(obj, "Call"))
  cat("\n")
}
