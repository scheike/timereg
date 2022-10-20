#' Fit Generalized Semiparametric Proportional 0dds Model
#' 
#' Fits a semiparametric proportional odds model: \deqn{ logit(1-S_{X,Z}(t)) =
#' log(X^T A(t)) + \beta^T Z } where A(t) is increasing but otherwise
#' unspecified. Model is fitted by maximising the modified partial likelihood.
#' A goodness-of-fit test by considering the score functions is also computed
#' by resampling methods.
#' 
#' An alternative way of writing the model : \deqn{ S_{X,Z}(t)) = \frac{ \exp(
#' - \beta^T Z )}{ (X^T A(t)) + \exp( - \beta^T Z) } } such that \eqn{\beta} is
#' the log-odds-ratio of dying before time t, and \eqn{A(t)} is the odds-ratio.
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
#' operator, and the terms on the right.  The response must be a survival
#' object as returned by the `Surv' function.
#' @param data a data.frame with the variables.
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
#' @param beta starting value for relative risk estimates
#' @param Nit number of iterations for Newton-Raphson algorithm.
#' @param detail if 0 no details is printed during iterations, if 1 details are
#' given.
#' @param sym to use symmetrized second derivative in the case of the
#' estimating equation approach (profile=0).  This may improve the numerical
#' performance.
#' @param mle.start starting values for relative risk parameters.
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
#' @references Scheike, A flexible semiparametric transformation model for
#' survival data, Lifetime Data Anal. (to appear).
#' 
#' Martinussen and Scheike, Dynamic Regression Models for Survival Data,
#' Springer (2006).
#' @keywords survival
#' @examples
#' 
#' data(sTRACE)
#' \donttest{
#' ### runs slowly and is therefore donttest 
#' data(sTRACE)
#' # Fits Proportional odds model with stratified baseline
#' age.c<-scale(sTRACE$age,scale=FALSE); 
#' ## out<-Gprop.odds(Surv(time,status==9)~-1+factor(diabetes)+prop(age.c)+prop(chf)+
#' ##                 prop(sex)+prop(vf),data=sTRACE,max.time=7,n.sim=50)
#' ## summary(out) 
#' ## par(mfrow=c(2,3))
#' ## plot(out,sim.ci=2); plot(out,score=1) 
#' 
#' }
#' 
#' @export
Gprop.odds<-function(formula=formula(data),data=parent.frame(),beta=0,Nit=50,detail=0,start.time=0,max.time=NULL,id=NULL,n.sim=500,weighted.test=0,sym=0,mle.start=0)
{
id.call<-id; call<-match.call(); residuals<-0;  
robust<-0; ratesim<-0; profile<-0; exppar<-0; 
m<-match.call(expand.dots = FALSE); 
m$max.time<-m$start.time<-m$weighted.test<-m$n.sim<-m$id<-m$Nit<-m$detail<-m$beta<-m$sym<-m$mle.start<-NULL
if (n.sim==0) sim<-0 else sim<-1; 
antsim<-n.sim; 

  special <- c("prop")
  Terms <- if(missing(data)) terms(formula, special)
           else              terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")

   XZ<-model.matrix(Terms,m)[, drop = FALSE]
   cols<-attributes(XZ)$assign
   l.cols<-length(cols)

   semicov <- attr(Terms, "specials")$prop
   ZtermsXZ<-semicov-1
#   print("semiterms"); print(semicov-1); 

   if (length(semicov)) { renaalen<-FALSE; Zterms<-c();
   for (i in ZtermsXZ) Zterms<-c(Zterms,(1:l.cols)[cols==i]); 
    } else {renaalen<-TRUE;}

    if (length(semicov)) {
       X<-as.matrix(XZ[,-Zterms]);
       covnamesZ <- dimnames(XZ)[[2]][-Zterms];
       dimnames(X)[[2]]<-covnamesZ; 
       Z<-as.matrix(XZ[,Zterms]);
       covnamesX <- dimnames(XZ)[[2]][Zterms]; 
       dimnames(Z)[[2]]<-covnamesX; 
       }
    else {X<-as.matrix(XZ); covnamesZ<-dimnames(XZ)[[2]]}

 if ((nrow(X)!=nrow(data)) && (!is.null(id))) stop("Missing values in design matrix not allowed with id\n"); 

#    print("X"); print(X[1,]); print(covnamesX); 
#    if (renaalen==FALSE) {print("Z"); print(Z[1,]);
#    print(covnamesZ); 
#    }

    px <- ncol(X)
    if (renaalen == FALSE) pz <- ncol(Z) else pz <- 0
    pxz <- px + pz

  desX<-as.matrix(Z);
  if(is.matrix(desX) == TRUE) pg <- as.integer(dim(desX)[2])
  if(is.matrix(desX) == TRUE) nx <- as.integer(dim(desX)[1])

  desZ<-as.matrix(X); px<-ncol(desZ); 

  if (is.diag(  t(desZ) %*% desZ  )==TRUE) stratum <- 1 else stratum <- 0

if (!inherits(Y, "Surv")) stop("Response must be a survival object")

if (attr(m[, 1], "type") == "right") {
  time2  <- m[, 1][, "time"]; time   <- rep(0,length(time2));
  status <- m[, 1][, "status"]    } else 
if (attr(m[, 1], "type") == "counting") {
  time   <- m[, 1][,1]; time2  <- m[, 1][,2]; status <- m[, 1][,3]; } else {
  stop("only right-censored or counting processes data") } 

Ntimes <- sum(status); 

# adds random noise to make survival times unique
if (sum(duplicated(time2[status==1]))>0) {
#cat("Non unique survival times: break ties ! \n")
#cat("Break ties yourself\n");
ties<-TRUE
dtimes<-time2[status==1]; index<-(1:length(time2))[status==1];
ties<-duplicated(dtimes); nties<-sum(ties); index<-index[ties]
dt<-diff(sort(time2)); dt<-min(dt[dt>0]);
time2[index]<-time2[index]+runif(nties,0,min(0.001,dt/2));
} else ties<-FALSE; 

start<-time; stop<-time2; 
times<-c(start.time,time2[status==1]); times<-sort(times);
if (is.null(max.time)) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
times<-times[times<maxtimes]
Ntimes <- length(times); 

########################################################################
if (is.null(id)==TRUE) {antpers<-length(time); id<-0:(antpers-1); }
else { pers<-unique(id); antpers<-length(pers); 
       id<-as.integer(factor(id,labels=1:(antpers)))-1; }

if (is.null(beta)) beta<-rep(0,pg);
if (length(beta)!=pg) beta<-rep(0,pg);

if (residuals==1) {
cumAi<-matrix(0,Ntimes,antpers*1);
cumAiiid<-matrix(0,Ntimes,antpers*1); }
else { cumAi<-0; cumAiiid<-0; }

cumint<-matrix(0,Ntimes,px+1); vcum<-matrix(0,Ntimes,px+1);
Rvcu<-matrix(0,Ntimes,px+1);
score<-beta;
Varbeta<-matrix(0,pg,pg); Iinv<-matrix(0,pg,pg);
RVarbeta<-matrix(0,pg,pg);
if (sim==1) Uit<-matrix(0,Ntimes,50*pg) else Uit<-FALSE;

test<-matrix(0,antsim,2*px); testOBS<-rep(0,2*px); unifCI<-c();
testval<-c();
rani<--round(runif(1)*10000); 
Ut<-matrix(0,Ntimes,pg+1); simUt<-matrix(0,antsim,pg);
loglike<-0; 
########################################################################

#cat("Generalised Proportional odds model \n"); 
#dyn.load("Gprop-odds.so")

nparout<- .C("Gtranssurv",
as.double(times),as.integer(Ntimes),as.double(desZ),
as.integer(nx),as.integer(px),as.double(desX),
as.integer(nx),as.integer(pg),as.integer(antpers),
as.double(start),as.double(stop),as.double(beta),
as.integer(Nit),as.double(cumint),as.double(vcum),
as.double(loglike),as.double(Iinv),as.double(Varbeta),
as.integer(detail),as.integer(sim),as.integer(antsim),
as.integer(rani),as.double(Rvcu), as.double(RVarbeta),
as.double(test),as.double(testOBS), as.double(Ut),
as.double(simUt),as.double(Uit),as.integer(id),
as.integer(status),as.integer(weighted.test),as.double(score),
as.double(cumAi),as.double(cumAiiid),as.integer(residuals),
as.integer(exppar),as.integer(sym),as.integer(mle.start),as.integer(stratum),PACKAGE="timereg");

gamma<-matrix(nparout[[12]],pg,1);
cumint<-matrix(nparout[[14]],Ntimes,px+1);
vcum<-matrix(nparout[[15]],Ntimes,px+1);
Iinv<-matrix(nparout[[17]],pg,pg);
Varbeta<--matrix(nparout[[18]],pg,pg);
Rvcu<-matrix(nparout[[23]],Ntimes,px+1);
RVarbeta<--matrix(nparout[[24]],pg,pg);
score<-matrix(nparout[[33]],pg,1);
Ut<-matrix(nparout[[27]],Ntimes,pg+1);
loglike<-nparout[[16]]

if (residuals==1) {
cumAi<-matrix(nparout[[34]],Ntimes,antpers*1);
cumAiiid<-matrix(nparout[[35]],Ntimes,antpers*1);
cumAi<-list(time=times,dmg=cumAi,dmg.iid=cumAiiid);} else cumAi<-FALSE;

if (sim==1) {
Uit<-matrix(nparout[[29]],Ntimes,50*pg); UIt<-list();
for (i in (0:49)*pg) UIt[[i/pg+1]]<-as.matrix(Uit[,i+(1:pg)]);
simUt<-matrix(nparout[[28]],antsim,pg);
test<-matrix(nparout[[25]],antsim,2*px); testOBS<-nparout[[26]];
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
testUt<-FALSE;test<-FALSE;unifCI<-FALSE;supUtOBS<-FALSE;UIt<-FALSE;testOBS<-FALSE;testval<-FALSE;
pval.testBeq0<- pval.testBeqC<- obs.testBeq0<- obs.testBeqC<- sim.testBeq0<-
sim.testBeqC<-FALSE; testUt<-FALSE; sim.supUt<-FALSE;
}

ud<-list(t=Terms,cum=cumint,var.cum=vcum,robvar.cum=Rvcu,
gamma=gamma,var.gamma=Varbeta,robvar.gamma=RVarbeta,
resid.dMG=cumAi,D2linv=Iinv,score=score,loglike=loglike,
pval.testBeq0=pval.testBeq0,pval.testBeqC=pval.testBeqC,
obs.testBeq0=obs.testBeq0,obs.testBeqC=obs.testBeqC,
sim.testBeq0= sim.testBeq0,sim.testBeqC=sim.testBeqC,
conf.band=unifCI,
test.procProp=Ut,sim.test.procProp=UIt,pval.Prop=testUt,
sim.supProp=sim.supUt,prop.odds=TRUE)

colnames(ud$cum)<-colnames(ud$var.cum)<- c("time",covnamesZ)
if (robust==1) colnames(ud$robvar.cum)<- c("time",covnamesZ); 


if (px>0) {
if (sim==1) {
colnames(ud$test.procProp)<-c("time",covnamesX)
names(ud$pval.Prop)<-covnamesX
names(ud$conf.band)<-names(ud$pval.testBeq0)<-
names(ud$pval.testBeqC)<-names(ud$obs.testBeq0)<- 
names(ud$obs.testBeqC)<-colnames(ud$sim.testBeq0)<-covnamesZ;
} }

rownames(ud$gamma)<-c(covnamesX); colnames(ud$gamma)<-"estimate";
rownames(ud$score)<-c(covnamesX); colnames(ud$score)<-"score";
namematrix(ud$var.gamma,covnamesX);
namematrix(ud$robvar.gamma,covnamesX);
namematrix(ud$D2linv,covnamesX);


attr(ud,"Call")<-call; 
attr(ud,"Formula")<-formula; 
attr(ud,"id")<-id.call; 
class(ud)<-"cox.aalen"
return(ud); 
}
