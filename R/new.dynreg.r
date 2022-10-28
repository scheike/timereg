#' Fit time-varying regression model
#' 
#' Fits time-varying regression model with partly parametric components.
#' Time-dependent variables for longitudinal data.  The model assumes that the
#' mean of the observed responses given covariates is a linear time-varying
#' regression model :
#' 
#' \deqn{ E( Z_{ij} | X_{ij}(t) ) = \beta^T(t) X_{ij}^1(t) + \gamma^T X_{ij}^2(t) } 
#' where \eqn{Z_{ij}} is the j'th measurement at time t for the
#' i'th subject with covariates \eqn{X_{ij}^1} and \eqn{X_{ij}^2}. Resampling
#' is used for computing p-values for tests of timevarying effects.
#
#' The data for a subject is presented as multiple rows or 'observations', each
#' of which applies to an interval of observation (start, stop].  For counting
#' process data with the )start,stop] notation is used the 'id' variable is
#' needed to identify the records for each subject. The program assumes that
#' there are no ties, and if such are present random noise is added to break
#' the ties.
#' @param formula a formula object with the response on the left of a '~'
#' operator, and the independent terms on the right as regressors.
#' @param data a data.frame with the variables.
#' @param start.time start of observation period where estimates are computed.
#' @param max.time end of observation period where estimates are computed.
#' Estimates thus computed from [start.time, max.time]. Default is max of data.
#' @param id For timevarying covariates the variable must associate each record
#' with the id of a subject.
#' @param n.sim number of simulations in resampling.
#' @param weighted.test to compute a variance weighted version of the
#' test-processes used for testing time-varying effects.
#' @param aalenmod Aalen model for measurement times. Specified as a survival
#' model (see aalen function).
#' @param bandwidth bandwidth for local iterations. Default is 50\% of the
#' range of the considered observation period.
#' @param bhat initial value for estimates. If NULL local linear estimate is
#' computed.
#' @param meansub if '1' then the mean of the responses is subtracted before
#' the estimation is carried out.
#' @param resample returns resample processes.
#' @return returns an object of type "dynreg". With the following arguments:
#' \item{cum}{the cumulative regression coefficients. This is the efficient
#' estimator based on an initial smoother obtained by local linear regression :
#' \deqn{ \hat B(t) = \int_0^t \tilde \beta(s) ds+ \hspace{4 cm}} 
#'  \deqn{\int_0^t X^{-} (Diag(z) -Diag( X^T(s) \tilde \beta(s)) ) dp(ds \times dz), } 
#' where \eqn{\tilde \beta(t)} is an initial estimate either
#' provided or computed by local linear regression. To plot this estimate use
#' type="eff.smooth" in the plot() command. } 
#' \item{var.cum}{the martingale based pointwise variance estimates.} 
#' \item{robvar.cum}{robust pointwise variances estimates.} 
#' \item{gamma}{estimate of semi-parametric components of model.} 
#' \item{var.gamma}{variance for gamma.} 
#' \item{robvar.gamma}{robust variance for gamma.} 
#' \item{cum0}{simple estimate of cumulative regression
#' coefficients that does not use use an initial smoothing based estimate
#' \deqn{ \hat B_0(t) = \int_0^t X^{-} Diag(z) dp(ds \times dz). } To plot this
#' estimate use type="0.mpp" in the plot() command. } 
#' \item{var.cum0}{the martingale based pointwise variance estimates of cum0.}
#' \item{cum.ms}{estimate of cumulative regression coefficients based on
#' initial smoother (but robust to this estimator).  \deqn{ \hat B_{ms}(t) =
#' \int_0^t X^{-} (Diag(z)-f(s)) dp(ds \times dz), } where \eqn{f} is chosen as
#' the matrix \deqn{ f(s) = Diag( X^T(s) \tilde \beta(s)) ( I - X_\alpha(s) X_\alpha^-(s) ), }
#' where \eqn{X_{\alpha}} is the design for the sampling
#' intensities.
#' This is also an efficient estimator when the initial estimator is consistent
#' for \eqn{\beta(t)} and then asymptotically equivalent to cum, but small
#' sample properties appear inferior. Its variance is estimated by var.cum.
#' To plot this estimate use type="ms.mpp" in the plot() command. }
#' \item{cum.ly}{estimator where local averages are subtracted. Special case of
#' cum.ms. To plot this estimate use type="ly.mpp" in plot.  }
#' \item{var.cum.ly}{the martingale based pointwise variance estimates. }
#' \item{gamma0}{estimate of parametric component of model.  }
#' \item{var.gamma0}{estimate of variance of parametric component of model. }
#' \item{gamma.ly}{estimate of parametric components of model. }
#' \item{var.gamma.ly}{estimate of variance of parametric component of model. }
#' \item{gamma.ms}{estimate of variance of parametric component of model. }
#' \item{var.gamma.ms}{estimate of variance of parametric component of model.}
#' \item{obs.testBeq0}{observed absolute value of supremum of cumulative
#' components scaled with the variance.} \item{pval.testBeq0}{p-value for
#' covariate effects based on supremum test.} \item{sim.testBeq0}{resampled
#' supremum values.} \item{obs.testBeqC}{observed absolute value of supremum of
#' difference between observed cumulative process and estimate under null of
#' constant effect.} \item{pval.testBeqC}{p-value based on resampling.}
#' \item{sim.testBeqC}{resampled supremum values.}
#' \item{obs.testBeqC.is}{observed integrated squared differences between
#' observed cumulative and estimate under null of constant effect.}
#' \item{pval.testBeqC.is}{p-value based on resampling.}
#' \item{sim.testBeqC.is}{resampled supremum values.}
#' \item{conf.band}{resampling based constant to construct robust 95\% uniform
#' confidence bands.} 
#' \item{test.procBeqC}{observed test-process of difference
#' between observed cumulative process and estimate under null of constant effect.} 
#' \item{sim.test.procBeqC}{list of 50 random realizations of
#' test-processes under null based on resampling.}
#' \item{covariance}{covariances for nonparametric terms of model.}
#' @author Thomas Scheike
#' @references Martinussen and Scheike, Dynamic Regression Models for Survival
#' Data, Springer (2006).
#' @keywords survival
#' @examples
#' 
#' \donttest{ 
#' ## this runs slowly and is therfore donttest
#' data(csl)
#' indi.m<-rep(1,length(csl$lt)) 
#' 
#' # Fits time-varying regression model 
#' out<-dynreg(prot~treat+prot.prev+sex+age,data=csl,
#' Surv(lt,rt,indi.m)~+1,start.time=0,max.time=2,id=csl$id,
#' n.sim=100,bandwidth=0.7,meansub=0)
#' summary(out)
#' par(mfrow=c(2,3))
#' plot(out)
#' 
#' # Fits time-varying semi-parametric regression model.
#' outS<-dynreg(prot~treat+const(prot.prev)+const(sex)+const(age),data=csl,
#' Surv(lt,rt,indi.m)~+1,start.time=0,max.time=2,id=csl$id,
#' n.sim=100,bandwidth=0.7,meansub=0)
#' summary(outS)
#' }
#' 
##' @export
dynreg<-function(formula=formula(data),data=parent.frame(),aalenmod,
bandwidth=0.5,id=NULL,bhat=NULL,start.time=0,
max.time=NULL,n.sim=500,meansub=1,weighted.test=0,resample=0)
{
  if (n.sim==0) sim<-0 else sim<-1; smoothXX<-0
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");} 

  b<-bandwidth
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  m$weighted.test<-m$meansub<-m$bandwidth<-m$aalenmod<-m$start.time<-m$max.time<-
  m$return.mg<-m$n.sim<-m$bhat<-m$id<-m$clusters<-m$resample<-NULL
  special <- c("const")
  Terms <- if(missing(data)) terms(formula, special)
  else              terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
  pxz <- px + pz;
  XZ<-cbind(X,Z); 

  clusters <- des$clusters ##########
####################################################################
### Aalen design is interpreted  
  udaal<-aalen.des(aalenmod,data=data); 
  time<-udaal$time; time2<-udaal$time2; 
  covarA<-data.matrix(udaal$X); status<-udaal$status; 
  pa<-ncol(covarA); 


####################################################################
  meanY<-mean(Y); 
  if (meansub==1) { Y<-Y-meanY;}
  covar<-data.matrix(cbind(Y,XZ,covarA)); 
####################################################################
  Ntimes <- sum(status); 

   # adds random noise to make survival times unique
  if (sum(duplicated(time2[status==1]))>0) {
    # cat("Non unique survival times: break ties ! \n")
    # cat("Break ties yourself\n");
    ties<-TRUE; 
    dtimes<-time2[status==1]
    index<-(1:length(time2))[status==1]
    ties<-duplicated(dtimes); nties<-sum(ties); index<-index[ties]
    dt<-diff(sort(time2)); dt<-min(dt[dt>0]);
    time2[index]<-time2[index]+runif(nties,0,min(0.001,dt/2));
  } else ties<-FALSE; 

  if (is.null(id)==TRUE) stop("Must specify subject id variable \n")

  Ntimes <- Ntimes+1; 
  times<-c(0,time2[status==1]); 
  if (is.null(max.time)==TRUE) max.time<-max(times)+0.1 else max.time<-min(max(times),max.time);   
  times<-times[times<max.time & times>start.time]; 
  times<-sort(times); Ntimes<-length(times); 
  clusters<-cluster.call<-des$clusters; 

  if (is.null(id)==TRUE) {antpers<-length(time); id<-0:(antpers-1); }
  else { pers<-unique(id); antpers<-length(pers); 
         id<-as.integer(factor(id,labels=1:(antpers)))-1; }

  ldata<-list(start=time,stop=time2,
              antpers=antpers,antclust=des$antclust);
###  X<-as.matrix(covar); 

  if (npar==TRUE) {
  #cat("Nonparametric Additive Model "); cat("\n")
  # local linear regression for preliminary estimates 
####################################################
    bandwidth<-(max.time-start.time)*bandwidth;

    if (is.null(bhat)==TRUE) {
      xval<-seq(times[2],times[Ntimes],length=30); 
      bhat<-localTimeReg(time2,Y[status==1],X[status==1,],xval,b,lin=1)[,1:(px+1)]
    }

    ud<-dynregBase(times,status,Y,ldata,X,covarA,id,clusters,
                   sim=sim,resample=resample,antsim=n.sim,b=b,bhat=bhat,
                   smoothXX=smoothXX,weighted.test=weighted.test);

    colnames(ud$cum.ly)<- colnames(ud$var.cum.ly)<-  
      colnames(ud$cum)<-colnames(ud$var.cum)<- colnames(ud$cum0)<- 
        colnames(ud$cum.ms)<-colnames(ud$robvar.cum)<-c("time",covnamesX)

    if (sim==1) {
      colnames(ud$test.procBeqC)<- c("time",covnamesX)
      names(ud$conf.band)<- names(ud$pval.testBeq0)<- names(ud$pval.testBeqC)<- 
        names(ud$obs.testBeq0)<- names(ud$obs.testBeqC)<- 
          names(ud$obs.testBeqC.is)<- names(ud$pval.testBeqC.is)<- 
            colnames(ud$sim.testBeq0)<- colnames(ud$sim.testBeqC) <- 
              colnames(ud$sim.testBeqC.is) <- covnamesX; }
  }
  else {
   #cat(" Semiparametric Additive Model"); cat("\n")

    if (is.null(bhat)==TRUE) {
      cat(" Computes initial estimates based on local regression\n")
      cat(" for efficient estimates, you may provide these\n");
      xval<-seq(times[2],times[Ntimes],length=30); 
      bhat<-localTimeReg(time2,Y[status==1],XZ[status==1,],xval,b,lin=1); 
      bhat<-bhat[,1:(pxz+1)]; 
      gamma<-apply(as.matrix(bhat[,((px+2):(pxz+1))]),2,mean); 
      bhat<-bhat[,1:(px+1)]; 
    } else {bhat<-cbind(xval,matrix(0,30,px)); gamma<-rep(0,pz);}

#if (is.null(bhat)==TRUE) {
#ud<-dynregBase(times,status,Y,ldata,
#X,covarA,id,bhat=bhat, 
#sim=0,retur=0,antsim=0,b=b,smoothXX=smoothXX,
#weighted.test=weighted.test);
##pcregci(ud$cum,ud$var.cum,0,3); 
#xval<-seq(times[1],times[Ntimes],length=30); 
#bhat<-CsmoothB(ud$cum,xval,b); 
#gamma<-ud$cum[signif(Ntimes*3/4),(px+2):(px+pz+1)]/
#       ud$cum[signif(Ntimes*3/4),1];
#};
#print(apply(cbind(X[,1:(pxz+1)],X[,(pxz+2):(pxz+pa+1)]),2,mean))
#print(ud$cum[200,]); 

ud<-semiregBase(times,status,Y,ldata,X,Z,covarA,id,clusters,
 bhat=bhat,sim=sim,antsim=n.sim,b=b,gamma=gamma,weighted.test=weighted.test,
 resample=resample);

    if (px>0) {
      colnames(ud$cum)<- colnames(ud$var.cum)<- colnames(ud$cum0)<- 
        colnames(ud$cum.ms)<- colnames(ud$robvar.cum)<-c("time",covnamesX)

      if (sim==1) {
        colnames(ud$test.procBeqC)<- c("time",covnamesX)
        names(ud$conf.band)<- names(ud$pval.testBeq0)<- 
          names(ud$pval.testBeqC)<- 
            names(ud$pval.testBeqC.is)<- names(ud$obs.testBeqC.is)<- 
              names(ud$obs.testBeq0)<- names(ud$obs.testBeqC)<- 
                colnames(ud$sim.testBeq0)<- 

                  colnames(ud$sim.testBeqC.is)<- 
                    colnames(ud$sim.testBeqC)<- covnamesX; }
    }

    ud$gamma<-nameestimate(ud$gamma,covnamesZ); 
    ud$gamma.ms<-nameestimate(ud$gamma.ms,covnamesZ); 
    ud$gamma0<-nameestimate(ud$gamma0,covnamesZ); 
    ud$gamma.ly<-nameestimate(ud$gamma.ly,covnamesZ); 
#ud$gamma.ef<-nameestimate(ud$gamma.ef,covnamesZ); 
#ud$gamma.efms<-nameestimate(ud$gamma.efms,covnamesZ); 

    ud$var.gamma<-namematrix(ud$var.gamma,covnamesZ); 
    ud$robvar.gamma<-namematrix(ud$robvar.gamma,covnamesZ); 
    ud$var.gamma.ms<-namematrix(ud$var.gamma.ms,covnamesZ); 
    ud$var.gamma.ly<-namematrix(ud$var.gamma.ly,covnamesZ); 
#ud$var.gamma.ef<-namematrix(ud$var.gamma.ef,covnamesZ); 
#ud$robvar.gamma.ef<-namematrix(ud$robvar.gamma.ef,covnamesZ); 
ud$mean.response<-meanY;
  }
  attr(ud,"Call")<-call; 
  class(ud)<-"dynreg"
  return(ud); 
}

namematrix<-function(mat,names)
{ colnames(mat)<-names; rownames(mat)<-names; return(mat); }
nameestimate<-function(mat,names)
{ colnames(mat)<-"estimate"; rownames(mat)<-names; return(mat); }



#' Plots estimates and test-processes
#' 
#' This function plots the non-parametric cumulative estimates for the additive
#' risk model or the test-processes for the hypothesis of constant effects with
#' re-sampled processes under the null.
#' 
#' @param x the output from the "dynreg" function.
#' @param type the estimator plotted. Choices "eff.smooth", "ms.mpp", "0.mpp"
#' and "ly.mpp". See the dynreg function for more on this.
#' @param pointwise.ci if >1 pointwise confidence intervals are plotted with
#' lty=pointwise.ci
#' @param hw.ci if >1 Hall-Wellner confidence bands are plotted with lty=hw.ci.
#' Only 0.95 \% bands can be constructed.
#' @param sim.ci if >1 simulation based confidence bands are plotted with
#' lty=sim.ci. These confidence bands are robust to non-martingale behaviour.
#' @param robust robust standard errors are used to estimate standard error of
#' estimate, otherwise martingale based estimate are used.
#' @param specific.comps all components of the model is plotted by default, but
#' a list of components may be specified, for example first and third "c(1,3)".
#' @param level gives the significance level.
#' @param start.time start of observation period where estimates are plotted.
#' @param stop.time end of period where estimates are plotted. Estimates thus
#' plotted from [start.time, max.time].
#' @param add.to.plot to add to an already existing plot.
#' @param mains add names of covariates as titles to plots.
#' @param xlab label for x-axis.
#' @param ylab label for y-axis.
#' @param score to plot test processes for test of time-varying effects along
#' with 50 random realization under the null-hypothesis.
#' @param ... unused arguments - for S3 compatibility
#' @author Thomas Scheike
#' @references Martinussen and Scheike, Dynamic Regression Models for Survival
#' Data, Springer (2006).
#' @keywords survival
#' @examples
#' 
#' \donttest{
#' ### runs slowly and therefore donttest 
#' data(csl)
#' indi.m<-rep(1,length(csl$lt))
#' 
#' # Fits time-varying regression model
#' out<-dynreg(prot~treat+prot.prev+sex+age,csl,
#' Surv(lt,rt,indi.m)~+1,start.time=0,max.time=3,id=csl$id,
#' n.sim=100,bandwidth=0.7,meansub=0)
#' 
#' par(mfrow=c(2,3))
#' # plots estimates 
#' plot(out)
#' # plots tests-processes for time-varying effects 
#' plot(out,score=TRUE)
#' }
#' 
##' @export
plot.dynreg<-function(x,type="eff.smooth",pointwise.ci=1,hw.ci=0,
sim.ci=0,robust=0,specific.comps=FALSE,level=0.05,start.time=0,stop.time=0,
add.to.plot=FALSE,mains=TRUE,xlab="Time",ylab ="Cumulative coefficients",score=FALSE,...)
{
  object <- x; rm(x);  
  if (!inherits(object, 'dynreg')) 
    stop ("Must be output from dynreg() function") 
 
  if (score==FALSE) {
    if (type=="eff.smooth") { B<-object$cum;V<-object$var.cum;}
    else if (type=="ms.mpp") {B<-object$cum.ms;V<-object$var.cum;}
    else if (type=="0.mpp") { B<-object$cum0; 
                              if (is.numeric(object$gamma)==FALSE) V<-object$var.cum0 else V<-B*0; }
    else if (type=="ly.mpp") {
      if (is.numeric(object$gamma)==FALSE) {
        B<-object$cum.ly; V<-object$var.cum.ly;} else 
      stop("Non-par estimates not computed for LY correction\n"); }
    else stop("not valid type"); 

    p<-dim(B)[[2]]; 
    if (is.null(V)==TRUE) robust<-1; 
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
      if (hw.ci>=1) {
        if (level!=0.05) cat("Hall-Wellner bands only 95 % \n");
        tau<-length(B[,1])
        nl<-B[,v]-1.27*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.27*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="s"); 
        lines(B[,1],nl,lty=hw.ci,type="s"); }
      if (sim.ci>=1) {
        if (sum(object$conf.band)==FALSE)
          cat("Uniform simulation based bands only computed for n.sim> 0\n")
        if (level!=0.05) c.alpha<-percen(object$sim.testBeq0[,v-1],1-level)
        else c.alpha<-object$conf.band[v-1];
        nl<-B[,v]-c.alpha*Vrob[,v]^.5; ul<-B[,v]+c.alpha*Vrob[,v]^.5;
        lines(B[,1],ul,lty=sim.ci,type="s"); 
        lines(B[,1],nl,lty=sim.ci,type="s"); }
      abline(h=0)
    }
  } else {
                                        # plot score proces
    dim1<-ncol(object$test.procBeqC)
    if (sum(specific.comps)==FALSE) comp<-2:dim1 else comp<-specific.comps+1

    for (i in comp)
      {
        ul<-2*max(abs(object$test.procBeqC[,i]));
        plot(object$test.procBeqC[,1],
             object$test.procBeqC[,i],type="l",ylim=c(-ul,ul),
             lwd=2,xlab=xlab,ylab=ylab)
        if (mains==TRUE) title(main=colnames(object$test.procBeqC)[i]); 
        for (j in 1:50)
          lines(object$test.procBeqC[,1],
                as.matrix(object$sim.test.procBeqC[[j]])[,i-1],lwd=1,col="grey",lty=1)
        lines(object$test.procBeqC[,1],object$test.procBeqC[,i],lwd=2)
      } }

}

##' @export
"print.dynreg" <-
function (x,...) 
{
  dynreg.object <- x; rm(x);
  if (!inherits(dynreg.object, 'dynreg')) 
    stop ("Must be a dynreg object")

if (is.null(dynreg.object$gamma0)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Dynamic Additive Regression Model \n\n")
  cat(" Nonparametric terms : "); cat(colnames(dynreg.object$cum)[-1]);
  cat("   \n");  
  if (semi) {
  cat(" Parametric terms :  "); cat(rownames(dynreg.object$gamma0)); 
  cat("   \n");  }
  cat("   \n");  

  cat("  Call: \n")
  dput(attr(dynreg.object, "Call"))
  cat("\n")
}


##' @export
"summary.dynreg" <- function(object,digits = 3,...) 
{
  dynreg.object <- object; rm(object);
  obj<-dynreg.object
  if (!inherits(dynreg.object, 'dynreg')) stop ("Must be an dynreg object")
  if (is.null(dynreg.object$gamma.ms)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Dynamic Additive Regression Model \n\n")
  cat(" Nonparametric terms : "); cat(colnames(dynreg.object$cum)[-1]);
  cat("   \n");  

  timetest(obj,digits=digits);

  if (semi) {
    cat(" Parametric terms :  "); cat(rownames(dynreg.object$gamma0)); 
    cat("   \n");  
    out=coef.dynreg(dynreg.object); 
    out=signif(out,digits=digits)
    print(out)
    cat("   \n");  
  }

###  cat("  Call: \n")
###  dput(attr(dynreg.object, "Call"))
  cat("\n")
}


##' @export
coef.dynreg<- function(object,...,digits=3) {
   coefBase(object,digits=digits)
}

##' @export
aalen.des<-function(formula=formula(data),data=parent.frame(),model="aalen")
{
  call <- match.call(); 
  m <- match.call(expand.dots=FALSE); 
  m$model<-NULL
  special <- c("cluster","prop","const")
  Terms <- if(missing(data)) terms(formula,special) else terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m,"response")

  des<-read.design(m,Terms,model=model)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px;
  pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ;
  clusters<-des$clusters;

  if (attr(m[, 1], "type") == "right") {
    type<-"right"; 
    status <- m[, 1][, "status"];
    time2  <- m[, 1][, "time"]; time   <- rep(0,length(time2)); 
  } else if (attr(m[, 1], "type") == "counting") {
    type<-"counting"; 
    time   <- m[, 1][,1]; time2  <- m[, 1][,2]; status <- m[, 1][,3];
  } else { stop("only right-censored or counting processes data") } 
return(list(type=type,time=time,time2=time2,status=status,
 X=X,Z=Z,px=px,pz=pz,npar=npar,
 covnamesX=covnamesX,covnamesZ=covnamesZ,clusters=clusters))
}

