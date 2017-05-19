#' Competings Risks Regression
#' 
#' Fits a semiparametric model for the cause-specific quantities : \deqn{ P(T <
#' t, cause=1 | x,z) = P_1(t,x,z) = h( g(t,x,z) ) } for a known link-function
#' \eqn{h()} and known prediction-function \eqn{g(t,x,z)} for the probability
#' of dying from cause 1 in a situation with competing causes of death.
#' 
#' We consider the following models : 1) the additive model where
#' \eqn{h(x)=1-\exp(-x)} and \deqn{ g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta
#' } 2) the proportional setting that includes the Fine & Gray (FG) "prop"
#' model and some extensions where \eqn{h(x)=1-\exp(-\exp(x))} and \deqn{
#' g(t,x,z) = (x^T A(t) + (diag(t^p) z)^T \beta) } The FG model is obtained
#' when \eqn{x=1}, but the baseline is parametrized as \eqn{\exp(A(t))}.
#' 
#' The "fg" model is a different parametrization that contains the FG model,
#' where \eqn{h(x)=1-\exp(-x)} and \deqn{ g(t,x,z) = (x^T A(t)) \exp((diag(t^p)
#' z)^T \beta) } The FG model is obtained when \eqn{x=1}.
#' 
#' 3) a "logistic" model where \eqn{h(x)=\exp(x)/( 1+\exp(x))} and \deqn{
#' g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta }
#' 
#' The "logistic2" is \deqn{ P_1(t,x,z) = x^T A(t) exp((diag(t^p) z)^T \beta)/
#' (1+ x^T A(t) exp((diag(t^p) z)^T \beta)) } The simple logistic model with
#' just a baseline can also be fitted by an alternative procedure that has
#' better small sample properties see prop.odds.subist().
#' 
#' 4) the relative cumulative incidence function "rcif" model where
#' \eqn{h(x)=\exp(x)} and \deqn{ g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta }
#' 
#' The "rcif2" \deqn{ P_1(t,x,z) = (x^T A(t)) \exp((diag(t^p) z)^T \beta) }
#' 
#' Where p by default is 1 for the additive model and 0 for the other models.
#' In general p may be powers of the same length as z.
#' 
#' Since timereg version 1.8.4. the response must be specified with the
#' \code{\link{Event}} function instead of the \code{\link{Surv}} function and
#' the arguments. For example, if the old code was
#' 
#' comp.risk(Surv(time,cause>0)~x1+x2,data=mydata,cause=mydata$cause,causeS=1)
#' 
#' the new code is
#' 
#' comp.risk(Event(time,cause)~x1+x2,data=mydata,cause=1)
#' 
#' Also the argument cens.code is now obsolete since cens.code is an argument
#' of \code{\link{Event}}.
#' 
#' @param formula a formula object, with the response on the left of a '~'
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the `Event' function. The status indicator is not important
#' here. Time-invariant regressors are specified by the wrapper const(), and
#' cluster variables (for computing robust variances) by the wrapper cluster().
#' @param data a data.frame with the variables.
#' @param cause For competing risk models specificies which cause we consider.
#' @param times specifies the times at which the estimator is considered.
#' Defaults to all the times where an event of interest occurs, with the first
#' 10 percent or max 20 jump points removed for numerical stability in
#' simulations.
#' @param Nit number of iterations for Newton-Raphson algorithm.
#' @param clusters specifies cluster structure, for backwards compability.
#' @param est possible starting value for nonparametric component of model.
#' @param fix.gamma to keep gamma fixed, possibly at 0.
#' @param gamma starting value for constant effects.
#' @param n.sim number of simulations in resampling.
#' @param weighted Not implemented. To compute a variance weighted version of
#' the test-processes used for testing time-varying effects.
#' @param model "additive", "prop"ortional, "rcif", or "logistic".
#' @param detail if 0 no details are printed during iterations, if 1 details
#' are given.
#' @param interval specifies that we only consider timepoints where the
#' Kaplan-Meier of the censoring distribution is larger than this value.
#' @param resample.iid to return the iid decomposition, that can be used to
#' construct confidence bands for predictions
#' @param cens.model specified which model to use for the ICPW, KM is
#' Kaplan-Meier alternatively it may be "cox"
#' @param cens.formula specifies the regression terms used for the regression
#' model for chosen regression model. When cens.model is specified, the default
#' is to use the same design as specified for the competing risks model.
#' @param time.pow specifies that the power at which the time-arguments is
#' transformed, for each of the arguments of the const() terms, default is 1
#' for the additive model and 0 for the proportional model.
#' @param time.pow.test specifies that the power the time-arguments is
#' transformed for each of the arguments of the non-const() terms. This is
#' relevant for testing if a coefficient function is consistent with the
#' specified form A_l(t)=beta_l t^time.pow.test(l). Default is 1 for the
#' additive model and 0 for the proportional model.
#' @param silent if 0 information on convergence problems due to non-invertible
#' derviates of scores are printed.
#' @param conv gives convergence criterie in terms of sum of absolute change of
#' parameters of model
#' @param weights weights for estimating equations.
#' @param max.clust sets the total number of i.i.d. terms in i.i.d.
#' decompostition. This can limit the amount of memory used by coarsening the
#' clusters. When NULL then all clusters are used.  Default is 1000 to save
#' memory and time.
#' @param first.time.p first point for estimation is pth percentile of cause
#' jump times.
#' @param n.times only uses 50 points for estimation, if NULL then uses all
#' points, subject to p.start condition.
#' @param estimator default estimator is 1.
#' @param trunc.p truncation weight for delayed entry, P(T > entry.time | Z_i),
#' typically Cox model.
#' @param cens.weights censoring weights can be given here rather than
#' calculated using the KM, cox or aalen models.
#' @param admin.cens censoring times for the administrative censoring
#' @param conservative set to 0 to compute correct variances based on censoring
#' weights, default is conservative estimates that are much quicker.
#' @param monotone monotone=0, uses estimating equations \deqn{ (D_\beta P_1)
#' w(t) ( Y(t)/G_c(t) - P_1(t,X)) and } montone 1 uses \deqn{ w(t) (
#' Y(t)/G_c(t) - P_1(t,X)) and }
#' @param step step size for Fisher-Scoring algorithm.
#' @return returns an object of type 'comprisk'. With the following arguments:
#' \item{cum}{cumulative timevarying regression coefficient estimates are
#' computed within the estimation interval.} \item{var.cum}{pointwise variances
#' estimates.  } \item{gamma}{estimate of proportional odds parameters of
#' model.} \item{var.gamma}{variance for gamma.  } \item{score}{sum of absolute
#' value of scores.} \item{gamma2}{estimate of constant effects based on the
#' non-parametric estimate. Used for testing of constant effects.}
#' \item{obs.testBeq0}{observed absolute value of supremum of cumulative
#' components scaled with the variance.} \item{pval.testBeq0}{p-value for
#' covariate effects based on supremum test.} \item{obs.testBeqC}{observed
#' absolute value of supremum of difference between observed cumulative process
#' and estimate under null of constant effect.} \item{pval.testBeqC}{p-value
#' based on resampling.} \item{obs.testBeqC.is}{observed integrated squared
#' differences between observed cumulative and estimate under null of constant
#' effect.} \item{pval.testBeqC.is}{p-value based on resampling.}
#' \item{conf.band}{resampling based constant to construct 95\% uniform
#' confidence bands.} \item{B.iid}{list of iid decomposition of non-parametric
#' effects.} \item{gamma.iid}{matrix of iid decomposition of parametric
#' effects.} \item{test.procBeqC}{observed test process for testing of
#' time-varying effects} \item{sim.test.procBeqC}{50 resample processes for for
#' testing of time-varying effects} \item{conv}{information on convergence for
#' time points used for estimation.}
#' @author Thomas Scheike
#' @references Scheike, Zhang and Gerds (2008), Predicting cumulative incidence
#' probability by direct binomial regression,Biometrika, 95, 205-220.
#' 
#' Scheike and Zhang (2007), Flexible competing risks regression modelling and
#' goodness of fit, LIDA, 14, 464-483.
#' 
#' Martinussen and Scheike (2006), Dynamic regression models for survival data,
#' Springer.
#' @keywords survival
#' @examples
#' 
#' data(bmt); 
#' 
#' clust <- rep(1:204,each=2)
#' addclust<-comp.risk(Event(time,cause)~platelet+age+tcell+cluster(clust),data=bmt,
#' cause=1,resample.iid=1,n.sim=100,model="additive")
#' ###
#' 
#' addclust<-comp.risk(Event(time,cause)~+1+cluster(clust),data=bmt,cause=1,
#' 		    resample.iid=1,n.sim=100,model="additive")
#' pad <- predict(addclust,X=1)
#' plot(pad)
#' 
#' add<-comp.risk(Event(time,cause)~platelet+age+tcell,data=bmt,
#' cause=1,resample.iid=1,n.sim=100,model="additive")
#' summary(add)
#' 
#' par(mfrow=c(2,4))
#' plot(add); 
#' ### plot(add,score=1) ### to plot score functions for test
#' 
#' ndata<-data.frame(platelet=c(1,0,0),age=c(0,1,0),tcell=c(0,0,1))
#' par(mfrow=c(2,3))
#' out<-predict(add,ndata,uniform=1,n.sim=100)
#' par(mfrow=c(2,2))
#' plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=1)
#' 
#' add<-comp.risk(Event(time,cause)~platelet+age+tcell,data=bmt,
#' 	       cause=1,resample.iid=0,n.sim=0,cens.model="cox",
#' 	       cens.formula=~factor(platelet),model="additive")
#' 
#' out<-predict(add,ndata,se=0,uniform=0)
#' par(mfrow=c(2,2))
#' plot(out,multiple=0,se=0,uniform=0,col=1:3,lty=1)
#' 
#' ## fits additive model with some constant effects 
#' add.sem<-comp.risk(Event(time,cause)~
#' const(platelet)+const(age)+const(tcell),data=bmt,
#' cause=1,resample.iid=1,n.sim=100,model="additive")
#' summary(add.sem)
#' 
#' out<-predict(add.sem,ndata,uniform=1,n.sim=100)
#' par(mfrow=c(2,2))
#' plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=0)
#' 
#' ## Fine & Gray model 
#' fg<-comp.risk(Event(time,cause)~
#' const(platelet)+const(age)+const(tcell),data=bmt,
#' cause=1,resample.iid=1,model="fg",n.sim=100)
#' summary(fg)
#' 
#' out<-predict(fg,ndata,uniform=1,n.sim=100)
#' 
#' par(mfrow=c(2,2))
#' plot(out,multiple=1,uniform=0,col=1:3,lty=1,se=0)
#' 
#' ## extended model with time-varying effects
#' fg.npar<-comp.risk(Event(time,cause)~platelet+age+const(tcell),
#' data=bmt,cause=1,resample.iid=1,model="prop",n.sim=100)
#' summary(fg.npar); 
#' 
#' out<-predict(fg.npar,ndata,uniform=1,n.sim=100)
#' head(out$P1[,1:5]); head(out$se.P1[,1:5])
#' 
#' par(mfrow=c(2,2))
#' plot(out,multiple=1,uniform=0,col=1:3,lty=1,se=0)
#' 
#' ## Fine & Gray model with alternative parametrization for baseline
#' fg2<-comp.risk(Event(time,cause)~const(platelet)+const(age)+const(tcell),data=bmt,
#' cause=1,resample.iid=1,model="prop",n.sim=100)
#' summary(fg2)
#' 
#' #################################################################
#' ## Delayed entry models, 
#' #################################################################
#' nn <- nrow(bmt)
#' entrytime <- rbinom(nn,1,0.5)*(bmt$time*runif(nn))
#' bmt$entrytime <- entrytime
#' times <- seq(5,70,by=1)
#' 
#' bmtw <- prep.comp.risk(bmt,times=times,time="time",entrytime="entrytime",cause="cause")
#' 
#' ## non-parametric model 
#' outnp <- comp.risk(Event(time,cause)~tcell+platelet+const(age),
#' 		   data=bmtw,cause=1,fix.gamma=1,gamma=0,
#'  cens.weights=bmtw$cw,weights=bmtw$weights,times=times,n.sim=0)
#' par(mfrow=c(2,2))
#' plot(outnp)
#' 
#' outnp <- comp.risk(Event(time,cause)~tcell+platelet,
#' 		   data=bmtw,cause=1,
#'  cens.weights=bmtw$cw,weights=bmtw$weights,times=times,n.sim=0)
#' par(mfrow=c(2,2))
#' plot(outnp)
#' 
#' 
#' ## semiparametric model 
#' out <- comp.risk(Event(time,cause)~const(tcell)+const(platelet),data=bmtw,cause=1,
#'  cens.weights=bmtw$cw,weights=bmtw$weights,times=times,n.sim=0)
#' summary(out)
#' 
#' 
##' @export
comp.risk<-function(formula,data=sys.parent(),cause,times=NULL,Nit=50,clusters=NULL,est=NULL,
		    fix.gamma=0,gamma=0,n.sim=0,weighted=0,model="fg",detail=0,interval=0.01,resample.iid=1,
                    cens.model="KM",cens.formula=NULL,time.pow=NULL,time.pow.test=NULL,silent=1,conv=1e-6,
                    weights=NULL,max.clust=1000,n.times=50,first.time.p=0.05,estimator=1,
		    trunc.p=NULL,cens.weights=NULL,admin.cens=NULL,conservative=1,monotone=0,step=NULL) 
    # {{{
{
    if (!missing(cause)){
       if (length(cause)!=1) stop("Argument cause has new meaning since 
   timereg version 1.8.4., it now specifies the cause of interest, see help(comp.risk) for details.")
   } 

    ## {{{
    # trans=1 P_1=1-exp( - ( x' b(b)+ z' gam t) ), 
    # trans=2 P_1=1-exp(-exp(x a(t)+ z` b )  Fine-Gray model, with baseline exp(x a(t)) 
    # trans=3 P_1= exp(x a(t)+ z` b)/( exp(x a(t) + z' b) +1 );  logistic
    # trans=4 P_1=exp( ( x' b(b)+ z' gam ) ), 
    # trans=5 P_1= (x' b(t)) exp( z' gam ), 
    # trans=6 P_1=1-exp(-(x a(t)) exp(z` b )) Fine-Gray model, with baseline x a(t) 
    # trans=7 P_1= (x a(t)) exp( z` b)/( (x a(t) ) exp(z' b) +1 ); logistic2
    trans <- switch(model,additive=1,prop=2,logistic=3,rcif=4,rcif2=5,fg=6,logistic2=7)
    ###  if (model=="additive")  trans<-1; if (model=="prop")      trans<-2; if (model=="logistic")  trans<-3; 
    ###  if (model=="rcif")      trans<-4; if (model=="rcif2")     trans<-5; if (model=="fg")        trans<-6; 
    ###  if (model=="logistic2") trans<-7; 
    line <- 0
    cause.call <- causeS <- cause
    m <- match.call(expand.dots=FALSE);
    m$gamma<-m$times<-m$n.times<-m$cause<-m$Nit<-m$weighted<-m$n.sim<-
             m$model<-m$detail<- m$cens.model<-m$time.pow<-m$silent<- m$step <- 
             m$cens.formula <- m$interval<- m$clusters<-m$resample.iid<- m$monotone <- 
             m$time.pow.test<-m$conv<- m$weights  <- m$max.clust <- m$first.time.p<- m$trunc.p <- 
             m$cens.weights <- m$admin.cens <- m$fix.gamma <- m$est  <- m$conservative <- m$estimator <- NULL
  
    if ((trans==2 || trans==3 || trans==7) && is.null(step)) step <- 0.5
    if (is.null(step)) step <- 1

    special <- c("const","cluster")
    if (missing(data)) {
        Terms <- terms(formula, special)
    }  else {
        Terms <- terms(formula, special, data = data)
    }
    m$formula <- Terms

    if (substr(as.character(m$formula)[2],1,4)=="Hist") {
       stop("Since timereg version 1.8.6.: The left hand side of the formula must be specified as 
       Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
    }

    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    if (NROW(m) == 0) stop("No (non-missing) observations")
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    event.history <- model.extract(m, "response")

  if (class(event.history)!="Event"){
       stop("Since timereg version 1.8.6.: The left hand side of the formula must be specified as 
       Event(time, event) or with non default censoring codes Event(time, event, cens.code=0).")
  }

   model.type <- "competing.risks"

   ## {{{ Event stuff
    cens.code <- attr(event.history,"cens.code")
    if (ncol(event.history)==2) {
	time2 <- eventtime <- event.history[,1]
	status <- delta  <- event.history[,2]
	entrytime <- rep(0,length(time2))
	left <- 0
    } else {
	time2 <- eventtime <- event.history[,2]
	status <- delta  <- event.history[,3]
	entrytime <- event.history[,1]
	left <- 1
	if (max(entrytime)==0) left <- 0
    }
    event <- (abs(status)==cause)
    if (sum(event)==0) stop("No events of interest in data\n"); 

    ## }}} 

  if (n.sim==0) sim<-0 else sim<-1; antsim<-n.sim;
  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ;

  if (nrow(X)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 

  if (is.diag(t(X) %*% X)==TRUE) stratum <- 1 else stratum <- 0; 

  ## {{{ cluster set up
  if(is.null(clusters)){ clusters <- des$clusters}
  if(is.null(clusters)){
    cluster.call<-clusters; 
    clusters <- 0:(nrow(X) - 1)
    antclust <- nrow(X)
  } else {
    cluster.call<-clusters; 
    antclust <- length(unique(clusters))
    clusters <- as.integer(factor(clusters,labels=1:antclust))-1
  }

    coarse.clust <- FALSE; 
    if ((!is.null(max.clust))) if (max.clust< antclust) {
        coarse.clust <- TRUE
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
	antclust <- max.clust    
    }                                                         
    ## }}} 

    pxz <-px+pz;
    
    if (is.null(times)) {
        timesc<-sort(unique(eventtime[event==1])); 
        if (!is.null(n.times)) {
            if (length(timesc)> n.times) times <- quantile(timesc,prob=seq(first.time.p,1,length=n.times)) 
            else times <- timesc
        } else {times<-timesc; times<-times[times> quantile(timesc,prob=first.time.p)]; }
    } else times <- sort(times); 

  n<-nrow(X); ntimes<-length(times);
  if (npar==TRUE) {Z<-matrix(0,n,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-pz;} 
  if (is.null(weights)==TRUE) weights <- rep(1,n); 
  ## }}}

  ## {{{ censoring and estimator 
  if (!is.null(admin.cens)) estimator  <- 3;
  Gcxe <- 1;  
  ordertime <- order(eventtime); 
  ###dcumhazcens <- rep(0,n); 

    if (estimator==1 || estimator==2) {
        if (is.null(cens.weights)) { ## {{{ censoring model stuff with possible truncation
            if (cens.model=="KM") { ## {{{
	        if (left==1) ud.cens<-survfit(Surv(entrytime,eventtime,delta==cens.code)~+1) else 
		ud.cens<-survfit(Surv(eventtime,delta==cens.code)~+1)
                Gfit<-cbind(ud.cens$time,ud.cens$surv)
                Gfit<-rbind(c(0,1),Gfit); 
                Gcx<-Cpred(Gfit,eventtime,strict=TRUE)[,2];
                Gcxe<-Cpred(Gfit,entrytime,strict=TRUE)[,2];
		### strictly before, but starts in 1. 
		Gcxe[Gcxe==0] <- 1
		### only conditional on L if trunc given 
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
                Gctimes<-Cpred(Gfit,times,strict=TRUE)[,2]; ## }}}
            } else if (cens.model=="stratKM") { ## {{{
	        XZ <- model.matrix(cens.formula,data=data); 
	        strata <- as.factor(XZ)
		Gcx <- pred.stratKM(data,time=eventtime,cause=delta,strata=strata)
		### only conditional on L if trunc given 
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
                Gctimes<-Cpred(Gfit,times)[,2]; ## }}}
            } else if (cens.model=="cox") { ## {{{
                if (!is.null(cens.formula)) { 
		      XZ <- model.matrix(cens.formula,data=data); 
                      if (sum(XZ[,1])==nrow(XZ)) XZ <- as.matrix(XZ[,-1])
                } else {
                       if (npar==TRUE) XZ<-X[,-1] else XZ <-cbind(X,Z)[,-1];
                }
		if (left==1) ud.cens<-coxph(Surv(entrytime,eventtime,delta==cens.code)~XZ)                
		else ud.cens<-coxph(Surv(eventtime,delta==cens.code)~XZ)                
		baseout <- basehaz(ud.cens,centered=FALSE); 
		baseout <- cbind(baseout$time,baseout$hazard)
		Gcx<-Cpred(baseout,eventtime,strict=TRUE)[,2];
		Gcxe<-Cpred(baseout,entrytime,strict=TRUE)[,2];
		Gcxe[Gcxe==0] <- 1
		RR<-exp(as.matrix(XZ) %*% coef(ud.cens))
		Gcx<-exp(-Gcx*RR)
		Gcxe<-exp(-Gcxe*RR)
		Gfit<-rbind(c(0,1),cbind(eventtime,Gcx)); 
		### only conditional on L if trunc given 
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
		Gctimes<-Cpred(Gfit,times,strict=TRUE)[,2]; 
                ## }}}
            } else if (cens.model=="aalen") {  ## {{{
                if (!is.null(cens.formula)) { 
			XZ <- model.matrix(cens.formula,data=data); 
                 } else {
                      if (npar==TRUE) XZ <-X else XZ <-cbind(X,Z);
                }
	        if (left==1) ud.cens<-aalen(Surv(entrytime,eventtime,delta==cens.code)~-1+XZ+cluster(clusters),
			       n.sim=0,residuals=0,robust=0,silent=1)
	        else ud.cens<-aalen(Surv(eventtime,delta==cens.code)~-1+XZ+cluster(clusters),
			       n.sim=0,residuals=0,robust=0,silent=1); 
                Gcx <- Cpred(ud.cens$cum,eventtime,strict=TRUE)[,-1];
                Gcx<-exp(-apply(Gcx*XZ,1,sum))
                Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-1
                Gcxe <- Cpred(ud.cens$cum,entrytime,strict=TRUE)[,2];
		Gcxe[Gcxe==0] <- 1
		if (!is.null(trunc.p)) Gcx <- Gcx/Gcxe; 
                Gfit<-rbind(c(0,1),cbind(eventtime,Gcx)); 
                Gctimes<-Cpred(Gfit,times,strict=TRUE)[,2]; ## }}}
            } else  stop('Unknown censoring model') 
            cens.weights <- Gcx
            if ((min(Gcx[event==1])< 0.00001) && (silent==0)) { 
                cat("Censoring dist. approx zero for some points, summary cens:\n");
                print(summary(Gcx)) 
            }
            ## }}}
        } else { 
            if (length(cens.weights)!=n) stop("censoring weights must have length equal to nrow in data\n");  
            Gcx <- cens.weights
	    ### for left truncation specification
            ord2 <- order(time2)
            Gctimes <- Cpred(cbind(time2[ord2],weights[ord2]),times)
        }
    } else { ## estimator==3 admin.cens 
        if (length(admin.cens)!=n) stop("censoring weights must have length equal to nrow in data\n");  
        Gcx <- admin.cens
        Gctimes <- rep(1,length(times)); 
    }

   if (left==1 & is.null(trunc.p) & is.null(cens.weights))  {  ## {{{ 
	 ### geskus weights: from mstate crprep 
	 stop("For left-truncated data call prep.comp.risk\n call with weights and cens.weights\n"); 
         n=length(time2)
         prec.factor <- 100
         prec <- .Machine$double.eps * prec.factor
         surv.trunc <- survfit(Surv(-time2,-entrytime+prec,rep(1,n)) ~ 1)
         trunc.dist <- summary(surv.trunc)
         trunc.dist$time <- rev(-trunc.dist$time)
         trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
         Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),time2)
         Lw <- Lfit[,2]
###	 weights <- 1/Lw
	 weights <- 1/((Lw)*Gcx); 
	 weights[delta==cens.code] <- 0
	 Gcx <- rep(1,n)
   } ## }}} 
   if (is.null(trunc.p)) trunc.p <- rep(1,n);  
   if (length(trunc.p)!=n) stop("truncation weights must have same length as data\n"); 
## }}}


## {{{ setting up more variables

  if (resample.iid == 1) {
    biid <- double(ntimes* antclust * px);
    gamiid<- double(antclust *pg);
  } else {
    gamiid <- biid <- NULL;
  }

  ps<-px; 
  betaS<-rep(0,ps); 

  ## possible starting value for nonparametric components
  if (is.null(est)) { 
	  est<-matrix(0.0+0.1,ntimes,px+1); 
          est[,1] <- times; 
  }  else {
      est <- as.matrix(est); 
  }
  if (nrow(est)!=length(times)) est <- Cpred(est,times); 

  hess<-matrix(0,ps,ps); var<-score<-matrix(0,ntimes,ps+1); 
  if (sum(gamma)==0) gamma<-rep(0,pg); gamma2<-rep(0,ps); 
  test<-matrix(0,antsim,3*ps); testOBS<-rep(0,3*ps); unifCI<-c();
  testval<-c(); rani<--round(runif(1)*10000); 
  Ut<-matrix(0,ntimes,ps+1); simUt<-matrix(0,ntimes,50*ps);
  var.gamma<-matrix(0,pg,pg); 
  pred.covs.sem<-0

  if (is.null(time.pow)==TRUE & model=="prop" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="fg" )       time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="additive")  time.pow<-rep(1,pg); 
  if (is.null(time.pow)==TRUE & model=="rcif" )     time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="rcif2" )    time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="logistic" ) time.pow<-rep(0,pg); 
  if (is.null(time.pow)==TRUE & model=="logistic2" )time.pow<-rep(0,pg); 
  if (length(time.pow)!=pg) time.pow <- rep(time.pow[1],pg); 

  if (is.null(time.pow.test)==TRUE & model=="prop" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="fg" )     time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="additive")  time.pow.test<-rep(1,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif" )    time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="rcif2" )   time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="logistic" ) time.pow.test<-rep(0,px); 
  if (is.null(time.pow.test)==TRUE & model=="logistic2" ) time.pow.test<-rep(0,px); 
  if (length(time.pow.test)!=px) time.pow.test <- rep(time.pow.test[1],px); 

  if (ntimes>1) silent <- c(silent,rep(0,ntimes-1))
  ## }}}

###    print(Gctimes); 

  ###  dyn.load("comprisk.so"0
  ssf <- step;  ## takes step size over 
  out<-.C("itfit", ## {{{
          as.double(times),as.integer(ntimes),as.double(eventtime),
          as.integer(cens.code), as.integer(status),as.double(Gcx),
          as.double(X),as.integer(n),as.integer(px),
          as.integer(Nit), as.double(betaS), as.double(score),
          as.double(hess), as.double(est), as.double(var),
          as.integer(sim),as.integer(antsim),as.integer(rani),
          as.double(test), as.double(testOBS), as.double(Ut),
          as.double(simUt),as.integer(weighted),as.double(gamma),
          as.double(var.gamma),as.integer(fixed),as.double(Z),
          as.integer(pg),as.integer(trans),as.double(gamma2),
          as.integer(cause),as.integer(line),as.integer(detail),
          as.double(biid),as.double(gamiid),as.integer(resample.iid),
          as.double(time.pow),as.integer(clusters),as.integer(antclust),
          as.double(time.pow.test),as.integer(silent), as.double(conv),
	  as.double(weights),as.double(entrytime),as.double(trunc.p),
	  as.integer(estimator),as.integer(fix.gamma), as.integer(stratum),
	  as.integer(ordertime-1),as.integer(conservative), as.double(ssf), 
	  as.double(Gctimes),as.double(rep(0,pg)),as.double(matrix(0,pg,pg)),
	  as.integer(monotone),PACKAGE="timereg") ## }}}
 
 ## {{{ handling output
  ssf <- out[[51]]; 
  gamma<-matrix(out[[24]],pg,1); 
  var.gamma<-matrix(out[[25]],pg,pg); 
  Dscore.gamma<-matrix(out[[54]],pg,pg); 
  gamma2<-matrix(out[[30]],ps,1); 
  rownames(gamma2)<-covnamesX; 
  
  conv <- list(convp=out[[41]],convd=out[[42]]);
    
  if (fixed==0) gamma<-NULL; 

  if (resample.iid==1)  {
    biid<-matrix(out[[34]],ntimes,antclust*px);
    if (fixed==1) gamiid<-matrix(out[[35]],antclust,pg) else gamiid<-NULL; 
    B.iid<-list();
    for (i in (0:(antclust-1))*px) {
    B.iid[[i/px+1]]<-matrix(biid[,i+(1:px)],ncol=px);
      colnames(B.iid[[i/px+1]])<-covnamesX; }
    if (fixed==1) colnames(gamiid)<-covnamesZ
  } else B.iid<-gamiid<-NULL;

  if (sim==1) {
    simUt<-matrix(out[[22]],ntimes,50*ps); UIt<-list();
    for (i in (0:49)*ps) UIt[[i/ps+1]]<-as.matrix(simUt[,i+(1:ps)]);
    Ut<-matrix(out[[21]],ntimes,ps+1);
    test<-matrix(out[[19]],antsim,3*ps); testOBS<-out[[20]];
    supUtOBS<-apply(abs(as.matrix(Ut[,-1])),2,max);
    p<-ps
    for (i in 1:(3*p)) testval<-c(testval,pval(test[,i],testOBS[i]))
    for (i in 1:p) unifCI<-as.vector(c(unifCI,percen(test[,i],0.95)));
    pval.testBeq0<-as.vector(testval[1:p]);
    pval.testBeqC<-as.vector(testval[(p+1):(2*p)]);
    pval.testBeqC.is<-as.vector(testval[(2*p+1):(3*p)]);
    obs.testBeq0<-as.vector(testOBS[1:p]);
    obs.testBeqC<-as.vector(testOBS[(p+1):(2*p)]);
    obs.testBeqC.is<-as.vector(testOBS[(2*p+1):(3*p)]);
    sim.testBeq0<-as.matrix(test[,1:p]);
    sim.testBeqC<-as.matrix(test[,(p+1):(2*p)]);
    sim.testBeqC.is<-as.matrix(test[,(2*p+1):(3*p)]);
  } else {test<-unifCI<-Ut<-UIt<-pval.testBeq0<-pval.testBeqC<-obs.testBeq0<-
          obs.testBeqC<- sim.testBeq0<-sim.testBeqC<-
          sim.testBeqC.is<- pval.testBeqC.is<-
          obs.testBeqC.is<-NULL;
  }

    est<-matrix(out[[14]],ntimes,ps+1);
    ## in case of no convergence set estimates to NA 
    est[conv$convp>0,-1] <- NA
    score<-matrix(out[[12]],ntimes,ps+1); 
    gamscore <- matrix(out[[53]],pg,1)
    scores <- list(score=score,gamscore=gamscore)
    var<-matrix(out[[15]],ntimes,ps+1);
    ## in case of no convergence set var to NA 
    var[conv$convp>0,-1] <- NA
    colnames(var)<-colnames(est)<-c("time",covnamesX); 

  if (sim>=1) {
    colnames(Ut)<- c("time",covnamesX)
    names(unifCI)<-names(pval.testBeq0)<- names(pval.testBeqC)<- 
    names(pval.testBeqC.is)<- names(obs.testBeq0)<- names(obs.testBeqC)<- 
    names(obs.testBeqC.is)<- colnames(sim.testBeq0)<- colnames(sim.testBeqC)<- 
    colnames(sim.testBeqC.is)<- covnamesX;
  }

    if (fixed==1) { rownames(gamma)<-c(covnamesZ);
                    colnames(var.gamma)<- rownames(var.gamma)<-c(covnamesZ); }

    colnames(score)<-c("time",covnamesX);
    if (is.na(sum(score))==TRUE) score<-NA  else 
    if (sum(score[,-1])<0.00001) score<-sum(score[,-1]); 
  ud<-list(cum=est,var.cum=var,gamma=gamma,score=score,
           gamma2=gamma2,var.gamma=var.gamma,robvar.gamma=var.gamma,
           pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
           obs.testBeq0=obs.testBeq0,
           obs.testBeqC.is=obs.testBeqC.is,
           obs.testBeqC=obs.testBeqC,pval.testBeqC.is=pval.testBeqC.is,
           conf.band=unifCI,B.iid=B.iid,gamma.iid=gamiid,ss=ssf,
           test.procBeqC=Ut,sim.test.procBeqC=UIt,conv=conv,
	   weights=weights,cens.weights=cens.weights,scores=scores,Dscore.gamma=Dscore.gamma,step=step)

    ud$call<- match.call()
    ud$model<-model; 
    ud$n<-n; 
    ud$clusters <- clusters
    ud$formula<-formula;
    ud$response <- event.history
    ud$cause <- status
    class(ud)<-"comprisk"; 
    attr(ud, "Call") <- match.call()
    attr(ud, "Formula") <- formula
    attr(ud, "time.pow") <- time.pow
    attr(ud, "causeS") <- causeS
    attr(ud, "cause") <- status
    attr(ud, "cluster.call") <- cluster.call
    attr(ud, "coarse.clust") <- coarse.clust
    attr(ud, "max.clust") <- max.clust
    attr(ud, "clusters") <- clusters
    attr(ud, "cens.code") <- cens.code
    attr(ud, "times") <- times
    return(ud);  ## }}}
} ## }}}

##' @export
print.comprisk <- function (x,...) { ## {{{
  object <- x; rm(x);
  if (!inherits(object, 'comprisk')) stop ("Must be an comprisk object")
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:
  causeS <- attr(object,"causeS")
  print(causeS)
  cat(paste("\nAnalysed cause:",causeS,"\n"))      
  cat(paste("\nLink _function:",object$model,"\n\n"))
  cat(" Nonparametric terms : ");
  cat(colnames(object$cum)[-1]); cat("   \n");  
  if (semi) {
      cat(" Parametric terms :  ");
      cat(rownames(object$gamma)); 
      cat("   \n");
  } 

  if (object$conv$convd>=1) {
      if (all(object$conv$convp==1)){
          if (NROW(object$cum)>1){
              cat("\nWarning: problem with convergence at all time points\n")
          } else{
              cat("\nWarning: problem with convergence at the evaluation time.\n")
          }
      }else{
          cat("Warning: problem with convergence at the following time points:\n")
          cat(object$cum[object$conv$convp>0,1])
          cat("\nYou may try to readjust analyses by removing these time points\n")
      }
  }
  cat("   \n");  
} ## }}}

##' @export
coef.comprisk <- function(object, digits=3,...) { ## {{{
   coefBase(object,digits=digits)
} ## }}}

##' @export
summary.comprisk <- function (object,digits = 3,...) {  ## {{{
  if (!inherits(object, 'comprisk')) stop ("Must be a comprisk object")
  
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Competing risks Model \n\n")
  
  modelType<-object$model
  #if (modelType=="additive" || modelType=="rcif") 
 
  if (sum(object$obs.testBeq0)==FALSE) cat("No test for non-parametric terms\n") else
  timetest(object,digits=digits); 

  if (semi) { 
         if (sum(abs(object$score)>0.000001)) 
         cat("Did not converge, allow more iterations\n\n"); 

	 cat("Parametric terms : \n"); 
	 prmatrix(coef(object,digits=digits))
	 cat("   \n"); 
  }

  if (object$conv$convd>=1) {
       cat("WARNING problem with convergence for time points:\n")
       cat(object$cum[object$conv$convp>0,1])
       cat("\nReadjust analyses by removing points\n\n") }

} ## }}}

##' @export
vcov.comp.risk <- function(object, ...) {
  rv <- object$robvar.gamma
  if (!identical(rv, matrix(0, nrow = 1L, ncol = 1L))) rv # else return NULL
}

##' @export
plot.comprisk <-  function (x, pointwise.ci=1, hw.ci=0,
                            sim.ci=0, specific.comps=FALSE,level=0.05, start.time = 0,
                            stop.time = 0, add.to.plot=FALSE, mains=TRUE, xlab="Time",
                            ylab ="Coefficients",score=FALSE,...){
## {{{
  object <- x; rm(x);

  if (!inherits(object,'comprisk') ){
    stop ("Must be output from comp.risk function")
  }

  if (score==FALSE) {
    B<-object$cum;
    V<-object$var.cum;
    p<-dim(B)[[2]]; 

    if (sum(specific.comps)==FALSE){
      comp<-2:p
    } else {
      comp<-specific.comps+1
    }
    if (stop.time==0) {
      stop.time<-max(B[,1]);
    }

    med<-B[,1]<=stop.time & B[,1]>=start.time
    B<-B[med,];
    V<-V[med,]; 

    c.alpha<- qnorm(1-level/2)
    for (v in comp) { 
      c.alpha<- qnorm(1-level/2)
      est<-B[,v];
      ul<-B[,v]+c.alpha*V[,v]^.5;
      nl<-B[,v]-c.alpha*V[,v]^.5;
      if (add.to.plot==FALSE) {
        plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab,...) 
        if (mains==TRUE) title(main=colnames(B)[v]);
      } else {
        lines(B[,1],est,type="s");
      }
      if (pointwise.ci>=1) {
        lines(B[,1],ul,lty=pointwise.ci,type="s");
        lines(B[,1],nl,lty=pointwise.ci,type="s");
      }
      if (hw.ci>=1) {
        if (level!=0.05){
          cat("Hall-Wellner bands only 95 % \n");
        }
        tau<-length(B[,1])
        nl<-B[,v]-1.27*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.27*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="s"); 
        lines(B[,1],nl,lty=hw.ci,type="s");
      }
      if (sim.ci>=1) {
        if (is.null(object$conf.band)==TRUE){
          cat("Uniform simulation based bands only computed for n.sim> 0\n")
        }
        if (level!=0.05){
          c.alpha<-percen(object$sim.testBeq0[,v-1],1-level)
        } else {
          c.alpha<-object$conf.band[v-1];
        }
        nl<-B[,v]-c.alpha*V[,v]^.5;
        ul<-B[,v]+c.alpha*V[,v]^.5;
        lines(B[,1],ul,lty=sim.ci,type="s"); 
        lines(B[,1],nl,lty=sim.ci,type="s");
      }
      abline(h = 0)
    }
  } else {
    # plot score proces
    if (is.null(object$pval.testBeqC)==TRUE) {
      cat("Simulations not done \n"); 
      cat("To construct p-values and score processes under null n.sim>0 \n"); 
    } else {
      if (ylab=="Cumulative regression function"){ 
        ylab<-"Test process";
      }
      dim1<-ncol(object$test.procBeqC)
      if (sum(specific.comps)==FALSE){
        comp<-2:dim1
      } else {
        comp<-specific.comps+1
      }

      for (i in comp){
          ranyl<-range(object$test.procBeqC[,i]);
          for (j in 1:50){
            ranyl<-range(c(ranyl,(object$sim.test.procBeqC[[j]])[,i-1]));
          }
          mr<-max(abs(ranyl));

          plot(object$test.procBeqC[,1],
               object$test.procBeqC[,i],
               ylim=c(-mr,mr),lwd=2,xlab=xlab,ylab=ylab,type="s",...)
          if (mains==TRUE){
            title(main=colnames(object$test.procBeqC)[i]);
          }
          for (j in 1:50){
            lines(object$test.procBeqC[,1],
                  as.matrix(object$sim.test.procBeqC[[j]])[,i-1],col="grey",lwd=1,lty=1,type="s")
          }
          lines(object$test.procBeqC[,1],object$test.procBeqC[,i],lwd=2,type="s")
        }
    }
  }
} ## }}}



#' Set up weights for delayed-entry competing risks data for comp.risk function
#' 
#' Computes the weights of Geskus (2011) modified to the setting of the
#' comp.risk function. The returned weights are
#' \eqn{1/(H(T_i)*G_c(min(T_i,tau)))} and tau is the max of the times argument,
#' here \eqn{H} is the estimator of the truncation distribution and \eqn{G_c}
#' is the right censoring distribution.
#' 
#' 
#' @param data data frame for comp.risk.
#' @param times times for estimating equations.
#' @param entrytime name of delayed entry variable, if not given computes
#' right-censoring case.
#' @param time name of survival time variable.
#' @param cause name of cause indicator
#' @param cname name of censoring weight.
#' @param tname name of truncation weight.
#' @param strata strata variable to obtain stratified weights.
#' @param nocens.out returns only uncensored part of data-frame
#' @param cens.formula censoring model formula for Cox models for the
#' truncation and censoring model.
#' @param cens.code code for censoring among causes.
#' @param prec.factor precision factor, for ties between censoring/even times,
#' truncation times/event times
#' @param trunc.mintau specicies wether the truncation distribution is
#' evaluated in death times or death times minimum max(times), FALSE makes the
#' estimator equivalent to Kaplan-Meier (in the no covariate case).
#' @return Returns an object. With the following arguments: \item{dataw}{a
#' data.frame with weights.}
#' 
#' The function wants to make two new variables "weights" and "cw" so if these
#' already are in the data frame it tries to add an "_" in the names.
#' @author Thomas Scheike
#' @references Geskus (2011), Cause-Specific Cumulative Incidence Estimation
#' and the Fine and Gray Model Under Both Left Truncation and Right Censoring,
#' Biometrics (2011), pp 39-49.
#' 
#' Shen (2011), Proportional subdistribution hazards regression for
#' left-truncated competing risks data, Journal of Nonparametric Statistics
#' (2011), 23, 885-895
#' @keywords survival
#' @examples
#' 
#' data(bmt)
#' nn <- nrow(bmt)
#' entrytime <- rbinom(nn,1,0.5)*(bmt$time*runif(nn))
#' bmt$entrytime <- entrytime
#' times <- seq(5,70,by=1)
#' 
#' ### adds weights to uncensored observations
#' bmtw <- prep.comp.risk(bmt,times=times,time="time",
#' 		       entrytime="entrytime",cause="cause")
#' 
#' #########################################
#' ### nonparametric estimates
#' #########################################
#' ## {{{ 
#' ### nonparametric estimates, right-censoring only 
#' out <- comp.risk(Event(time,cause)~+1,data=bmt,
#' 		 cause=1,model="rcif2",
#' 		 times=c(5,30,70),n.sim=0)
#' out$cum
#' ### same as 
#' ###out <- prodlim(Hist(time,cause)~+1,data=bmt)
#' ###summary(out,cause="1",times=c(5,30,70))
#' 
#' ### with truncation 
#' out <- comp.risk(Event(time,cause)~+1,data=bmtw,cause=1,
#'   model="rcif2",
#'   cens.weight=bmtw$cw,weights=bmtw$weights,times=c(5,30,70),
#'   n.sim=0)
#' out$cum
#' ### same as
#' ###out <- prodlim(Hist(entry=entrytime,time,cause)~+1,data=bmt)
#' ###summary(out,cause="1",times=c(5,30,70))
#' ## }}} 
#' 
#' #########################################
#' ### Regression 
#' #########################################
#' ## {{{ 
#' ### with truncation correction
#' out <- comp.risk(Event(time,cause)~const(tcell)+const(platelet),data=bmtw,
#'  cause=1,cens.weight=bmtw$cw,
#'  weights=bmtw$weights,times=times,n.sim=0)
#' summary(out)
#' 
#' ### with only righ-censoring, standard call
#' outn <- comp.risk(Event(time,cause)~const(tcell)+const(platelet),data=bmt,
#' 	  cause=1,times=times,n.sim=0)
#' summary(outn)
#' ## }}} 
#' 
#' 
##' @export
prep.comp.risk <- function(data,times=NULL,entrytime=NULL,
			   time="time",cause="cause",cname="cweight",tname="tweight",
			   strata=NULL,nocens.out=TRUE,cens.formula=NULL,cens.code=0,
			   prec.factor=100,trunc.mintau=FALSE)
{ ## {{{ 
## {{{  geskus weights, up to min(T_i,max(times))
   if (is.null(times)) times <- max(data[,time])
   if (is.null(entrytime)) entrytime <- rep(0,nrow(data)) else entrytime <- data[,entrytime]
   mtt <- max(times)
   prec.factor <- 100
   prec <- .Machine$double.eps * prec.factor
   trunc.model <- cens.model <- NULL ## output of Cox models for entry cens

   if (is.null(cens.formula)) { 
   if (is.null(strata)) { ## {{{ 
	   if (!is.null(entrytime)) {
	   surv.trunc <- 
	   survfit(Surv(-data[,time],-entrytime+prec,rep(1,nrow(data))) ~ 1) 
	   trunc.dist <- summary(surv.trunc)
	   trunc.dist$time <- rev(-trunc.dist$time)
	   trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
	   if (trunc.mintau==TRUE) Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),pmin(mtt,data[,time])) else 
	   Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),data[,time])
	   Lw <- Lfit[,2]
	   } else Lw <- 1
	   ud.cens<- survfit(Surv(entrytime,data[,time],data[,cause]==0)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,pmin(mtt,data[,time]),strict=TRUE)[,2];
           weights <- 1/(Lw*Gcx); 
	   cweights <-  Gcx; 
	   tweights <-  Lw; 
   ### ## }}} 
   } else { ## {{{ 
	   ### compute for each strata and combine 
	  vstrata <- as.numeric(data[,strata])
          weights <- rep(1,nrow(data))
          cweights <- rep(1,nrow(data))
          tweights <- rep(1,nrow(data))
	  for (i in unique(vstrata)) { ## {{{ for each strata
	       who <- (vstrata == i)
	       if (sum(who) <= 1) stop(paste("strata",i,"less than 1 observation\n")); 
	   datas <- subset(data,who)
	   if (!is.null(entrytime)) {
		   entrytimes <- entrytime[who]
		   surv.trunc <- 
		   survfit(Surv(-datas[,time],-entrytimes+prec,rep(1,nrow(datas))) ~ +1) 
		   trunc.dist <- summary(surv.trunc)
		   trunc.dist$time <- rev(-trunc.dist$time)
		   trunc.dist$surv <- c(rev(trunc.dist$surv)[-1], 1)
          	   if (trunc.mintau==TRUE) Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),pmin(mtt,datas[,time])) else 
	           Lfit <-Cpred(cbind(trunc.dist$time,trunc.dist$surv),datas[,time])
		   Lw <- Lfit[,2]
	   } else Lw <- 1
	   ud.cens<- survfit(Surv(entrytimes,datas[,time],datas[,cause]==0)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,pmin(mtt,datas[,time]),strict=TRUE)[,2];
	   weights[who]<-  1/(Lw*Gcx); 
	   cweights[who]<-  Gcx; 
	   tweights[who]<-  Lw; 
          } ## }}} 
   } ## }}} 
   } else { ### cens.formula Cox models  ## {{{
        X <- model.matrix(cens.formula,data=data)[,-1,drop=FALSE]; 

	if (!is.null(entrytime)) {
		trunc.model <- coxph(Surv(-data[,time],-entrytime+prec,rep(1,nrow(data))) ~ X) 
		baseout <- basehaz(trunc.model,centered=FALSE); 
		baseout <- cbind(rev(-baseout$time),rev(baseout$hazard))
	###

	   if (trunc.mintau==TRUE) Lfit <-Cpred(baseout,pmin(mtt,data[,time]))[,-1] else 
		   Lfit <-Cpred(baseout,data[,time])[,-1]
		RR<-exp(as.matrix(X) %*% coef(trunc.model))
		Lfit<-exp(-Lfit*RR)
		Lw <- Lfit
	   } else Lw <- 1
###
	cens.model <- coxph(Surv(entrytime,data[,time],data[,cause]==0)~+X) 
        baseout <- basehaz(cens.model,centered=FALSE); 
	baseout <- cbind(baseout$time,baseout$hazard)
	Gfit<-Cpred(baseout,pmin(mtt,data[,time]),strict=TRUE)[,2];
	RR<-exp(as.matrix(X) %*% coef(cens.model))
	Gfit<-exp(-Gfit*RR)
        weights <- 1/(Lw*Gfit); 
        cweights <- Gfit
        tweights <- Lw
   } ## }}} 
   data[,cname] <- cweights
   data[,tname] <- tweights

   if (!is.null(entrytime)) {
   mint <- min(tweights); maxt <- min(tweights) 
   if (mint<0 | mint>1) warning("min(truncation weights) strange, maybe prec.factor should be different\n")
   if (maxt<0 | maxt>1) warning("max(truncation weights) strange, maybe prec.factor should be different\n")
   }

   if ("weights" %in% names(data)) {
       warning("Weights in variable 'weights_' \n")
       wname<- "weights_"
       data[,wname] <- weights
   } else data[,"weights"] <- weights
###
   if ("cw" %in% names(data)) {
     warning("cw weights in variable 'cw_' \n")
     cwname<- "cw_"
     data[,cwname] <- 1
   } else data[,"cw"] <- 1
###
   if (nocens.out) {
     med <- ((data[,time]>mtt & data[,cause]==cens.code)) | (data[,cause]!=cens.code)
     data <- data[med,]
   } 

   attr(data,"trunc.model") <- trunc.model
   attr(data,"cens.model") <- cens.model 
## }}} 
   return(data)
} ## }}} 

##' @export
pred.stratKM <- function(data,entrytime=NULL,time="time",cause="cause",strata="strata",event.code=0)
{ ## {{{ 

     if (is.numeric(time)) time <- time else {
	     if (!is.null(data)) time <- data[,time] else stop("time not given\n"); 
     }
     if (is.numeric(cause)) cause <- cause else {
	     if (!is.null(data)) cause <- data[,cause] else stop("cause not given\n"); 
     }
     if (is.numeric(strata)) strata <- strata else {
	     if (!is.null(data)) strata <- data[,strata] else stop("strata not given\n"); 
     }
  if (is.null(entrytime)) entrytime <- rep(0,nrow(data)) else {
     if (is.numeric(entrytime)) entrytime <- entrytime else {
	     if (!is.null(data)) entrytime <- data[,entrytime] else stop("entrytime not given\n"); 
     }
  }
  vstrata <- as.numeric(strata)
  weights <- rep(1,length((data)))
  for (i in unique(vstrata)) { ## {{{ for each strata
	   who <- (vstrata == i)
	   if (sum(who) <= 1) stop(paste("strata",i,"less than 1 observation\n")); 
	   times <- time[who]
	   causes <- cause[who]
	   entrytimes <- entrytime[who]
	   ud.cens<- survfit(Surv(entrytimes,times,causes==event.code)~+1) 
	   Gfit<-cbind(ud.cens$time,ud.cens$surv)
	   Gfit<-rbind(c(0,1),Gfit); 
	   Gcx<-Cpred(Gfit,times,strict=TRUE)[,2];
	   weights[who]<-  Gcx; 
   } ## }}} 
   return(weights); 
} ## }}} 


