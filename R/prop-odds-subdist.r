#' Fit Semiparametric Proportional 0dds Model for the competing risks
#' subdistribution
#' 
#' Fits a semiparametric proportional odds model: \deqn{ logit(F_1(t;X,Z)) =
#' log( A(t)) + \beta^T Z }{} where A(t) is increasing but otherwise unspecified.
#' Model is fitted by maximising the modified partial likelihood.  A
#' goodness-of-fit test by considering the score functions is also computed by
#' resampling methods.
#' 
#' An alternative way of writing the model : \deqn{ F_1(t;X,Z) = \frac{ \exp(
#' \beta^T Z )}{ (A(t)) + \exp( \beta^T Z) } }{} such that \eqn{\beta} is the
#' log-odds-ratio of cause 1 before time t, and \eqn{A(t)} is the odds-ratio.
#' 
#' The modelling formula uses the standard survival modelling given in the
#' \bold{survival} package.
#' 
#' The data for a subject is presented as multiple rows or "observations", each
#' of which applies to an interval of observation (start, stop]. The program
#' essentially assumes no ties, and if such are present a little random noise
#' is added to break the ties.
#' 
#' @param formula a formula object, with the response on the left of a '~'
#' operator, and the terms on the right.  The response must be an object as
#' returned by the `Event' function.
#' @param data a data.frame with the variables.
#' @param cause cause indicator for competing risks.
#' @param beta starting value for relative risk estimates
#' @param Nit number of iterations for Newton-Raphson algorithm.
#' @param detail if 0 no details is printed during iterations, if 1 details are
#' given.
#' @param start.time start of observation period where estimates are computed.
#' @param max.time end of observation period where estimates are computed.
#' Estimates thus computed from [start.time, max.time].  This is very useful to
#' obtain stable estimates, especially for the baseline. Default is max of
#' data.
#' @param id For timevarying covariates the variable must associate each record
#' with the id of a subject.
#' @param n.sim number of simulations in resampling.
#' @param weighted.test to compute a variance weighted version of the
#' test-processes used for testing time-varying effects.
#' @param profile use profile version of score equations.
#' @param sym to use symmetrized second derivative in the case of the
#' estimating equation approach (profile=0).  This may improve the numerical
#' performance.
#' @param cens.model specifies censoring model. So far only Kaplan-Meier "KM".
#' @param cens.formula possible formula for censoring distribution covariates.
#' Default all !
#' @param clusters to compute cluster based standard errors.
#' @param max.clust number of maximum clusters to be used, to save time in iid
#' decomposition.
#' @param baselinevar set to 0 to save time on computations.
#' @param weights additional weights.
#' @param cens.weights specify censoring weights related to the observations.
#' @return returns an object of type 'cox.aalen'. With the following arguments:
#' 
#' \item{cum}{cumulative timevarying regression coefficient estimates are
#' computed within the estimation interval.} \item{var.cum}{the martingale
#' based pointwise variance estimates.  } \item{robvar.cum}{robust pointwise
#' variances estimates.  } \item{gamma}{estimate of proportional odds
#' parameters of model.} \item{var.gamma}{variance for gamma.  }
#' \item{robvar.gamma}{robust variance for gamma.  } \item{residuals}{list with
#' residuals. Estimated martingale increments (dM) and corresponding time
#' vector (time).} \item{obs.testBeq0}{observed absolute value of supremum of
#' cumulative components scaled with the variance.}
#' \item{pval.testBeq0}{p-value for covariate effects based on supremum test.}
#' \item{sim.testBeq0}{resampled supremum values.} \item{obs.testBeqC}{observed
#' absolute value of supremum of difference between observed cumulative process
#' and estimate under null of constant effect.} \item{pval.testBeqC}{p-value
#' based on resampling.} \item{sim.testBeqC}{resampled supremum values.}
#' \item{obs.testBeqC.is}{observed integrated squared differences between
#' observed cumulative and estimate under null of constant effect.}
#' \item{pval.testBeqC.is}{p-value based on resampling.}
#' \item{sim.testBeqC.is}{resampled supremum values.}
#' \item{conf.band}{resampling based constant to construct robust 95\% uniform
#' confidence bands. } \item{test.procBeqC}{observed test-process of difference
#' between observed cumulative process and estimate under null of constant
#' effect over time.} \item{loglike}{modified partial likelihood, pseudo
#' profile likelihood for regression parameters.} \item{D2linv}{inverse of the
#' derivative of the score function.} \item{score}{value of score for final
#' estimates.} \item{test.procProp}{observed score process for proportional
#' odds regression effects.} \item{pval.Prop}{p-value based on resampling.}
#' \item{sim.supProp}{re-sampled supremum values.}
#' \item{sim.test.procProp}{list of 50 random realizations of test-processes
#' for constant proportional odds under the model based on resampling.}
#' @author Thomas Scheike
#' @references Eriksson, Li, Zhang and Scheike (2014), The proportional odds
#' cumulative incidence model for competing risks, Biometrics, to appear.
#' 
#' Scheike, A flexible semiparametric transformation model for survival data,
#' Lifetime Data Anal. (2007).
#' 
#' Martinussen and Scheike, Dynamic Regression Models for Survival Data,
#' Springer (2006).
#' @keywords survival
#' @examples
#' 
#' library(timereg)
#' data(bmt)
#' # Fits Proportional odds model 
#' out <- prop.odds.subdist(Event(time,cause)~platelet+age+tcell,data=bmt,
#'  cause=1,cens.model="KM",detail=0,n.sim=1000)
#' summary(out) 
#' par(mfrow=c(2,3))
#' plot(out,sim.ci=2); 
#' plot(out,score=1) 
#' 
#' # simple predict function without confidence calculations 
#' pout <- predictpropodds(out,X=model.matrix(~platelet+age+tcell,data=bmt)[,-1])
#' matplot(pout$time,pout$pred,type="l")
#' 
#' # predict function with confidence intervals
#' pout2 <- predict(out,Z=c(1,0,1))
#' plot(pout2,col=2)
#' pout1 <- predictpropodds(out,X=c(1,0,1))
#' lines(pout1$time,pout1$pred,type="l")
#' 
#' # Fits Proportional odds model with stratified baseline, does not work yet!
#' ###out <- Gprop.odds.subdist(Surv(time,cause==1)~-1+factor(platelet)+
#' ###prop(age)+prop(tcell),data=bmt,cause=bmt$cause,
#' ###cens.code=0,cens.model="KM",causeS=1,detail=0,n.sim=1000)
#' ###summary(out) 
#' ###par(mfrow=c(2,3))
#' ###plot(out,sim.ci=2); 
#' ###plot(out,score=1) 
#' 
#' @export
prop.odds.subdist<-function(formula,data=parent.frame(),cause=1,beta=NULL,
Nit=10,detail=0,start.time=0,max.time=NULL,id=NULL,n.sim=500,weighted.test=0,
profile=1,sym=0,cens.model="KM",cens.formula=NULL,
clusters=NULL,max.clust=1000,baselinevar=1,weights=NULL,
cens.weights=NULL)
{
## {{{ 
## {{{ 
 if (!missing(cause)){
    if (length(cause)!=1) stop("Argument cause specifies the cause of interest, see help(prop.odds.subdist) for details.")
 } 

    ## {{{ reading formula
    rate.sim <- 1
    cause.call <- causeS <- cause
    m<-match.call(expand.dots=FALSE);
    if (n.sim==0) sim<-0 else sim<-1; 
    antsim<-n.sim; id.call<-id; 
    residuals<-0;  robust<-1; resample.iid <- 1 
    m$cens.model <- m$cause <- m$sym<-m$profile <- m$max.time<- m$start.time<- m$weighted.test<- m$n.sim<-
    m$id<-m$Nit<-m$detail<-m$beta <- m$baselinevar <- m$clusters <- m$max.clust <- m$weights <-  NULL
    m$cens.weights <- m$cens.formula  <- NULL

    special <- c("cluster")
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
    if (substr(as.character(m$formula)[2],1,4)=="Surv") stop("Must call with Event(time,cause) \n"); 

    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    if (NROW(m) == 0) stop("No (non-missing) observations")
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    event.history <- model.extract(m, "response")
    ## }}} 

    ## {{{  Hist to Event  stuff
    if (match("Hist",class(event.history),nomatch=0)==1){
        stop("Since timereg version 1.8.6., the right hand side of the formula must be specified as Event(time, event) or Event(time, event, cens.code=0).")
    }

    ## }}} 

   ## {{{ Event stuff
    if (match("Surv",class(event.history),nomatch=0)==1){

    } else {
    cens.code <- attr(event.history,"cens.code")
    time2 <- eventtime <- event.history[,1]
    status <- delta  <- event.history[,2]
    event <- (status==cause)
    entrytime <- rep(0,length(time2))
    if (sum(event)==0) stop("No events of interest in data\n"); 
    }

    ## }}} 

    ## {{{ 

    if (n.sim==0) sim<-0 else sim<-1; antsim<-n.sim;
    desX <- model.matrix(Terms, m)[,-1,drop=FALSE]; 
    covnamesX<-dimnames(desX)[[2]]; 
    ###desX<-as.matrix(X);
    if(is.matrix(desX) == TRUE) pg <- as.integer(dim(desX)[2])
    if(is.matrix(desX) == TRUE) nx <- as.integer(dim(desX)[1])
    px<-1; 
    if ( (nx!=nrow(data)) & (!is.null(id))) stop("Missing values in design matrix not allowed with id \n"); 

    # adds random noise to make survival times unique
    time2 <- eventtime
    jtimes <- time2[event==1]
    if (sum(duplicated(jtimes))>0) {
    ties<-TRUE
    index<-(1:length(time2))[event==1]
    ties<-duplicated(jtimes); 
    nties<-sum(ties); 
    index<-index[ties]
    dt<-diff(sort(jtimes)); 
    dt<-min(dt[dt>0]);
    time2[index]<-time2[index]+runif(nties,0,min(0.001,dt/2));
    } else ties<-FALSE; 

   times<-time2[event==1]; 
   index<-(1:length(time2))[event==1];
   index <- index[order(times)]; 
   times<-sort(times);
   times <- c(start.time,times)
   index <- c(0,index)

   start <- entrytime
   stop <- time2 

if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
times<-times[times<=maxtimes]
Ntimes <- length(times); 
## }}} 

## {{{ cluster and id set up
if (is.null(id)==TRUE) {antpers<-length(time2); id<-0:(antpers-1); }
else { pers<-unique(id); antpers<-length(pers); 
       id<-as.integer(factor(id,labels=1:(antpers)))-1; 
}

cluster.call<-clusters; 
if (is.null(clusters)== TRUE) {clusters<-id; antclust<-antpers;} else {
       clus<-unique(clusters); antclust<-length(clus); 
       clusters <- as.integer(factor(clusters, labels = 1:(antclust))) - 1;
}

if ((!is.null(max.clust))) if (max.clust<antclust) {
	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
	clusters <- as.integer(qqc)-1
	max.clusters <- length(unique(clusters))
###	clusters <- as.integer(factor(qqc, labels = 1:max.clust)) -1
	antclust <- max.clust    
  }                
## }}} 

## {{{ setting up more variables

  if (resample.iid == 1) {
    biid <- double(Ntimes* antclust );
    gamiid<- double(antclust *pg);
  } else {
    gamiid <- biid <- NULL;
  }


if ((is.null(beta)==FALSE)) {
	if (length(c(beta))!=pg) beta <- rep(beta[1],pg); 
} else {
      beta<-coxph(Surv(eventtime,event)~desX)$coef
     if (max(abs(beta))>5) {
	     cat("Warning, starting values from Cox model large, may set beta=\n")
     }
}

if (residuals==1) {
   cumAi<-matrix(0,Ntimes,antpers*1);
   cumAiiid<-matrix(0,Ntimes,antpers*1); 
} else { cumAi<-0; cumAiiid<-0; }

cumint<-matrix(0,Ntimes,px+1); 
vcum<-matrix(0,Ntimes,px+1);
Rvcu<-matrix(0,Ntimes,px+1);
score<-beta;
Varbeta<-matrix(0,pg,pg); Iinv<-matrix(0,pg,pg);
RVarbeta<-matrix(0,pg,pg);
if (sim==1) Uit<-matrix(0,Ntimes,50*pg) else Uit<-NULL;

test<-matrix(0,antsim,2*px); testOBS<-rep(0,2*px); unifCI<-c();
testval<-c();
rani<--round(runif(1)*10000); 
Ut<-matrix(0,Ntimes,pg+1); simUt<-matrix(0,antsim,pg);
loglike<-0; 
## }}}

## {{{ censoring and estimator

if (is.null(cens.weights)) { ## {{{ censoring model stuff with possible truncation
if (cens.model=="KM") { ## {{{
    ud.cens<-survfit(Surv(time2,delta==cens.code)~+1);
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    KMti<-Cpred(Gfit,time2)[,2];
    KMtimes<-Cpred(Gfit,times)[,2]; ## }}}
  } else if (cens.model=="cox") { ## {{{

     if (!is.null(cens.formula)) desXc <- model.matrix(cens.formula,data=data)[,-1] else desXc <- desX; 

    ud.cens<-cox.aalen(Surv(time2,delta==cens.code)~prop(desXc),n.sim=0,robust=0)
###  baseout <- basehaz(ud.cens,centered=FALSE); 
###  baseout <- cbind(baseout$time,baseout$hazard)
    Gcx<-Cpred(ud.cens$cum,time2)[,2];
    RR<-exp(desXc %*% ud.cens$gamma)
    KMti<-exp(-Gcx*RR)
    KMtimes<-Cpred(cbind(time2,KMti),times)[,2]; 
    ## }}}
  } else if (cens.model=="aalen") {  ## {{{

     if (!is.null(cens.formula)) desXc <- model.matrix(cens.formula,data=data) else desXc <- desX; 

    ud.cens<-aalen(Surv(time2,delta==cens.code)~desXc+
    cluster(clusters),n.sim=0,residuals=0,robust=0,silent=1)
    KMti <- Cpred(ud.cens$cum,time2)[,-1];
    Gcx<-exp(-apply(Gcx*desXc,1,sum))
    Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    KMti <- Gcx
    KMtimes<-Cpred(Gfit,times)[,2]; ## }}}
    } 
    else if (cens.model=="test") {    
    KMti <- 1-time2/6; KMtimes <- 1-times/6;  
    }
    else if (cens.model=="po") {    
    KMti <- rep(1,length(time2)); KMtimes <- rep(1,length(times));
    }
    else  { stop('Unknown censoring model') }
    }  else {
       if (length(cens.weights)!=nx) stop("censoring weights must have length equal to nrow in data\n");  
       KMti <- cens.weights
       Gctimes <- rep(1,length(times)); 
       ord2 <- order(time2) 
       KMtimes<-Cpred(cbind(time2[ord2],cens.weights[ord2]),times)[,2]; 
    }
## }}}

    if (is.null(weights)) weights <- rep(1,nx); 
    if (length(weights)!=nx) stop("weights must have same length as data\n"); 

    tmp <- status!=0 & status!=causeS
    if (any(tmp)==TRUE) if (min(KMti[tmp])==0) stop("Other causes have censoring weight equal to 0\n");  
 ## }}}

###print(table(status))
###print(head(stop))
###print(head(times))
###print(Ntimes)
###print(sum(KMtimes))
###print(sum(KMti))
###print(cens.code)
###print(table(status))
###print(causeS)

nparout<- .C("posubdist2",
	as.double(times),as.integer(Ntimes),as.double(desX),
	as.integer(nx),as.integer(pg),as.integer(antpers),
	as.double(start),as.double(stop), as.double(beta),
	as.integer(Nit), as.double(cumint), as.double(vcum),
	as.double(Iinv),as.double(Varbeta),as.integer(detail),
	as.integer(sim),as.integer(antsim),as.integer(rani),
	as.double(Rvcu),as.double(RVarbeta),as.double(test),
	as.double(testOBS),as.double(Ut),as.double(simUt),
	as.double(Uit),as.integer(id),as.integer(status),
	as.integer(weighted.test),as.integer(rate.sim),as.double(score), 
	as.double(cumAi),as.double(cumAiiid),as.integer(residuals), 
	as.double(loglike),as.integer(profile),as.integer(sym),
	as.double(KMtimes),as.double(KMti),as.double(time2),
	as.integer(causeS), as.integer(index-1), as.integer(baselinevar),
	as.integer(clusters), as.integer(antclust), as.integer(cens.code), 
        as.double(biid),as.double(gamiid),as.double(weights),PACKAGE="timereg");

## {{{ output handling

gamma<-matrix(nparout[[9]],pg,1);
cumint<-matrix(nparout[[11]],Ntimes,px+1);
vcum<-matrix(nparout[[12]],Ntimes,px+1);
Iinv<-matrix(nparout[[13]],pg,pg);
Varbeta<--matrix(nparout[[14]],pg,pg);
Rvcu<-matrix(nparout[[19]],Ntimes,px+1);
RVarbeta<--matrix(nparout[[20]],pg,pg);
score<-matrix(nparout[[30]],pg,1);
Ut<-matrix(nparout[[23]],Ntimes,pg+1);
loglike<-nparout[[34]]

if (residuals==1) {
cumAi<-matrix(nparout[[31]],Ntimes,antpers*1);
cumAiiid<-matrix(nparout[[32]],Ntimes,antpers*1);
cumAi<-list(time=times,dmg=cumAi,dmg.iid=cumAiiid);} else cumAi<-NULL;

 if (resample.iid==1)  {
    biid<-matrix(nparout[[46]],Ntimes,antclust);
    gamiid<-matrix(nparout[[47]],antclust,pg) 
    gamiid  <-  t(Iinv %*% t(gamiid))
    B.iid<-list();
    for (i in (1:antclust)) {
    B.iid[[i]]<-matrix(biid[,i],ncol=1);
    colnames(B.iid[[i]])<-"Baselineiid"; 
    }
      colnames(gamiid)<-covnamesX
  } else B.iid<-gamiid<-NULL;


if (sim==1) {
Uit<-matrix(nparout[[25]],Ntimes,50*pg); UIt<-list();
for (i in (0:49)*pg) UIt[[i/pg+1]]<-as.matrix(Uit[,i+(1:pg)]);
simUt<-matrix(nparout[[24]],antsim,pg);
test<-matrix(nparout[[21]],antsim,2*px); testOBS<-nparout[[22]];
supUtOBS<-apply(abs(as.matrix(Ut[,-1])),2,max);
for (i in 1:(2*px)) testval<-c(testval,pval(test[,i],testOBS[i]));
for (i in 1:px) unifCI<-c(unifCI,percen(test[,i],0.95));
testUt<-c();
for (i in 1:pg) testUt<-c(testUt,pval(simUt[,i],supUtOBS[i]));

pval.testBeq0<-as.vector(testval[1:px]);
pval.testBeqC<-as.vector(testval[(px+1):(2*px)]);
obs.testBeq0<-as.vector(testOBS[1:px]);
obs.testBeqC<-as.vector(testOBS[(px+1):(2*px)]);
sim.testBeq0<-as.matrix(test[,1:px]);
sim.testBeqC<-as.matrix(test[,(px+1):(2*px)]);
sim.supUt<-as.matrix(simUt);
}

if (sim!=1) {
testUt<-NULL;test<-NULL;unifCI<-NULL;supUtOBS<-NULL;UIt<-NULL;testOBS<-NULL;testval<-NULL;
pval.testBeq0<- pval.testBeqC<- obs.testBeq0<- obs.testBeqC<- sim.testBeq0<-
sim.testBeqC<-NULL; testUt<-NULL; sim.supUt<-NULL;
}

ud<-list(cum=cumint,var.cum=vcum,robvar.cum=Rvcu,
gamma=gamma,var.gamma=Varbeta,robvar.gamma=RVarbeta,
resid.dMG=cumAi,D2linv=Iinv,score=score,loglike=loglike,
pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
obs.testBeq0=obs.testBeq0,obs.testBeqC=obs.testBeqC,
sim.testBeq0= sim.testBeq0,sim.testBeqC=sim.testBeqC,
conf.band=unifCI,
test.procProp=Ut,sim.test.procProp=UIt,pval.Prop=testUt,
sim.supProp=sim.supUt,prop.odds=TRUE,B.iid=B.iid,gamma.iid=gamiid,cens.weights=KMti)

colnames(ud$cum)<-colnames(ud$var.cum)<- c("time","Baseline")
colnames(ud$robvar.cum)<- c("time","Baseline"); 

if (px>0) {
if (sim==1) {
colnames(ud$test.procProp)<-c("time",covnamesX)
names(ud$pval.Prop)<-covnamesX
names(ud$conf.band)<-names(ud$pval.testBeq0)<-
names(ud$pval.testBeqC)<-names(ud$obs.testBeq0)<- 
names(ud$obs.testBeqC)<-colnames(ud$sim.testBeq0)<-"Baseline";
} }

rownames(ud$gamma)<-c(covnamesX); colnames(ud$gamma)<-"estimate";
rownames(ud$score)<-c(covnamesX); colnames(ud$score)<-"score";
namematrix(ud$var.gamma,covnamesX);
namematrix(ud$robvar.gamma,covnamesX);
namematrix(ud$D2linv,covnamesX);
## }}} 

attr(ud,"Call")<-call; 
attr(ud,"Formula")<-formula; 
attr(ud,"id")<-id.call; 
attr(ud,"baselinevar") <- 1
if (cens.model!="po") attr(ud,"type") <- "comprisk" else attr(ud,"type") <- "survival"
class(ud)<-"cox.aalen"
return(ud); 
} ## }}} 
## }}} 

#' @export
predictpropodds <- function(out,X=NULL,times=NULL)
{  ## {{{ 
beta     <- out$gamma
baseline <- out$cum[,2]
btimes   <- out$cum[,1]

if (!is.null(times)) btimes <- times; 

pcum <- Cpred(out$cum,btimes)
RR <- matrix(X,ncol=length(beta),byrow=TRUE) %*% beta
HRR  <-  outer(pcum[,2],exp(RR),"*")[,,1]
pred <- HRR/(1+HRR)

return(list(pred=pred,time=btimes))
}  ## }}}

#' @export
prop.odds.subdist.ipw <- function(compriskformula,glmformula,data=parent.frame(),cause=1,
			 max.clust=NULL,ipw.se=FALSE,...)
{ ## {{{ 
  ggl <- glm(glmformula,family='binomial',data=data)
  mat <-  model.matrix(glmformula,data=data);
  glmcovs <- attr(ggl$terms,"term.labels")
  data$ppp <- predict(ggl,type='response')

  dcc <- data[ggl$y==1,]
  ppp <- dcc$ppp
  udca <- prop.odds.subdist(compriskformula,data=dcc,cause=cause,weights=1/ppp,n.sim=0,
		    max.clust=max.clust,...)  
  ### iid of beta for comprisk model 
  compriskiid <- udca$gamma.iid

if (ipw.se==TRUE)  { ## {{{ 
###requireNamespace("lava"); 
###requireNamespace("NumDeriv"); 
	glmiid <-   lava::iid(ggl)
	mat <- mat[ggl$y==1,]
	par <- coef(ggl)

	compriskalpha <- function(par)
	{ ## {{{ 
	  rr <- mat %*% par
	  pw <- c(exp(rr)/(1+exp(rr)))
	  assign("pw",pw,envir=environment(compriskformula))
	  ud <- prop.odds.subdist(compriskformula,data=dcc,
			  cause=cause,
			  weights=1/pw,baselinevar=0,beta=udca$gamma,
			  Nit=1,n.sim=0,...)  
	  ud$score
	} ## }}} 

	DU <-  numDeriv::jacobian(compriskalpha,par)
	IDU <-  udca$D2linv %*% DU 
	alphaiid <-t( IDU %*% t(glmiid))
	###
	iidfull <- alphaiid
	###
	iidfull[ggl$y==1,] <- compriskiid + alphaiid[ggl$y==1,]
	###
	var2 <- t(iidfull) %*% iidfull
	se <- cbind(diag(var2)^.5); colnames(se) <- "se"
} else { iidfull <- NULL; var2 <- NULL; se <- NULL} ## }}} 

var.naive=udca$robvar.gamma
se.naive=matrix(diag(var.naive)^.5,nrow(var.naive),1); 
colnames(se.naive) <- "se.naive"

res <- list(iid=iidfull,coef=udca$gamma,var.naive=var.naive,
	    se.naive=se.naive,var=var2,se=se,
	    comprisk.ipw=udca)
class(res) <- "comprisk.ipw"
return(res)
} ## }}} 





