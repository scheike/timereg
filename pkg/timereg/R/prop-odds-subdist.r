prop.odds.subdist<-function(formula,data=sys.parent(),cause=1,beta=NULL,
Nit=10,detail=0,start.time=0,max.time=NULL,id=NULL,n.sim=500,weighted.test=0,
profile=1,sym=0,cens.model="KM",clusters=NULL,max.clust=100,baselinevar=1)
{
## {{{ 

 if (!missing(cause)){
    if (length(cause)!=1) stop("Argument cause specifies the cause of interest, see help(prop.odds.subdist) for details.")
 } 

    ## {{{ reading formula
    m<-match.call(expand.dots=FALSE);
    id.call<-id; 
    residuals<-0;  robust<-0; ratesim<-0; 
    m$cens.model <- m$cause <- 
    m$sym<-m$profile <- m$max.time<- m$start.time<- m$weighted.test<- m$n.sim<-
    m$id<-m$Nit<-m$detail<-m$beta <- m$baselinevar <- m$clusters <- m$max.clust <- NULL
    if (n.sim==0) sim<-0 else sim<-1; 
    antsim<-n.sim; 
    causeS <- cause
  
    special <- c("cluster")
    if (missing(data)) {
        Terms <- terms(formula, special)
    }  else {
        Terms <- terms(formula, special, data = data)
    }
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    if (NROW(m) == 0) stop("No (non-missing) observations")
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    event.history <- model.extract(m, "response")

    if (match("Hist",class(event.history),nomatch=0)==0){
        stop("Since timereg version 1.8.4., the right hand side of the formula must be specified as Hist(time, event) or Hist(time, event, cens.code).")
    }
    ## if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    cens.type <- attr(event.history,"cens.type")
    cens.code <- attr(event.history,"cens.code")
    delayed <- !(is.null(attr(event.history,"entry.type"))) && !(attr(event.history,"entry.type")=="")
    model.type <- attr(event.history,"model")
    if (model.type %in% c("competing.risks","survival")){
        if (cens.type %in% c("rightCensored","uncensored")){
            delta <- event.history[,"status"]
            if (model.type=="competing.risks"){
                states <- attr(event.history,"states")
                ## if (missing(cause)) stop("Argument cause is needed for competing risks models.")
                if (missing(cause)) {
                    cause <- states[1]
                    warning(paste("Argument cause is missing, analysing cause 1: ",cause,". Other causes are:",paste(states[-1],collapse=","),sep=""))
                }else {
                    if (!(cause %in% states)) stop(paste("Cause",cause," is not among the causes in data; these are:",paste(states,collapse=",")))
                    cause <- match(cause,states,nomatch=0)
                }
                ## event is 1 if the event of interest occured and 0 otherwise
                event <-  event.history[,"event"] == cause
                if (sum(event)==0) stop(paste("No events of type:", cause, "in data."))
            }
            else
                event <- delta
            if (delayed){
                stop("Delayed entry is not (not yet) supported.")
                entrytime <- event.history[,"entry"]
                eventtime <- event.history[,"time"]
            }else{
                eventtime <- event.history[,"time"]
                entrytime <- rep(0,length(eventtime))
            }
        }else{
            stop("Works only for right-censored data")
        }
    }else{stop("Response is neither competing risks nor survival.")}

    start <- entrytime
    stop <- time2 <- eventtime
    states <- as.integer(states)
    cens.code <- as.integer(cens.code)
    ccc <- status <- event.history[,"event"]
    ### back to orignal coding 
    ccc[status==(length(states)+1)] <- cens.code
    for (i in 1:length(states)) ccc[status==i] <- states[i]
    status <- event <- ccc

    if (n.sim==0) sim<-0 else sim<-1; antsim<-n.sim;
###    des<-read.design(m,Terms)
###    X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
###    print(head(Z))
###    covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ;

###  if (nrow(X)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 

###  if(is.null(clusters)){ clusters <- des$clusters}
###  if(is.null(clusters)){
###    cluster.call<-clusters; 
###    clusters <- 0:(nrow(X) - 1)
###    antclust <- nrow(X)
###  } else {
###    cluster.call<-clusters; 
###    antclust <- length(unique(clusters))
###    clusters <- as.integer(factor(clusters,labels=1:antclust))-1
###  }
###
###    coarse.clust <- FALSE; 
###    if ((!is.null(max.clust))) if (max.clust< antclust) {
###        coarse.clust <- TRUE
###	qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
###	qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
###	clusters <- as.integer(qqc)-1
###	max.clusters <- length(unique(clusters))
###	antclust <- max.clust    
###    }                                                         
    
###    pxz <-px+pz;
    
###  n<-nrow(X); ntimes<-length(times);
###  if (npar==TRUE) {Z<-matrix(0,n,1); pg<-1; fixed<-0;} else {fixed<-1;pg<-pz;} 
###  if (is.null(weights)==TRUE) weights <- rep(1,n); 

  ## }}}
###
###  ## {{{ 
###call<-match.call(expand.dots=FALSE); 
###id.call<-id; 
###residuals<-0;  robust<-0; ratesim<-0; # profile<-0; 
###m<-match.call(expand.dots = FALSE);
###m$causeS <- m$cens.code <- m$cens.model <- m$cause <- 
###m$sym<-m$profile<- m$max.time<- m$start.time<- m$weighted.test<- m$n.sim<-
###m$id<-m$Nit<-m$detail<-m$beta <- m$baselinevar<-m$clusters <- m$max.clust <- NULL
###if (n.sim==0) sim<-0 else sim<-1; 
###antsim<-n.sim; 
###
###Terms <- if(missing(data)) terms(formula)
###           else              terms(formula, data=data)
###m$formula <- Terms
###m[[1]] <- as.name("model.frame")
###m <- eval(m, sys.parent())
###mt <- attr(m, "terms")
###intercept<-attr(mt, "intercept")
###Y <- model.extract(m, "response")
###if (!inherits(Y, "Surv")) stop("Response must be a survival object")
###
###if (attr(m[, 1], "type") == "right") {
###  time2  <- m[, 1][, "time"]; time   <- rep(0,length(time2));
###  status <- m[, 1][, "status"]    } else 
###if (attr(m[, 1], "type") == "counting") {
###  time   <- m[, 1][,1]; time2  <- m[, 1][,2]; status <- m[, 1][,3]; } else {
###  stop("only right-censored or counting processes data") } 
##### }}} 
###

desX <- model.matrix(Terms, m)[,-1,drop=FALSE]; 
covnamesX<-dimnames(desX)[[2]]; 
###desX<-as.matrix(X);
if(is.matrix(desX) == TRUE) pg <- as.integer(dim(desX)[2])
if(is.matrix(desX) == TRUE) nx <- as.integer(dim(desX)[1])
px<-1; 

if ( (nx!=nrow(data)) & (!is.null(id))) stop("Missing values in design matrix not allowed with id \n"); 
Ntimes <- sum(event==cause); 

# adds random noise to make survival times unique
if (sum(duplicated(time2[event==cause]))>0) {
ties<-TRUE
index<-(1:length(time2))[event==cause]
dtimes<-time2[event==cause]; ties<-duplicated(dtimes); 
ties<-duplicated(dtimes); nties<-sum(ties); index<-index[ties]
dt<-diff(sort(time2)); 
dt<-min(dt[dt>0]);
time2[index]<-time2[index]+runif(nties,0,min(0.001,dt/2));
} else ties<-FALSE; 

stop<-time2; 
dtimes<-time2[event==cause]; 
times<-time2[event==cause]; 
index<-(1:length(time2))[event==cause];
index <- index[order(times)]; 
times<-sort(times);
times <- c(start.time,times)
index <- c(0,index)

if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
times<-times[times<=maxtimes]
Ntimes <- length(times); 

########################################################################
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

if ((is.null(beta)==FALSE)) {
	if (length(beta)!=pg) beta <- rep(beta[1],pg); 
} else {
      beta<-coxph(Surv(time2,status==cause)~desX)$coef
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

if (cens.model=="KM") { ## {{{
    ud.cens<-survfit(Surv(time2,status==cens.code)~+1);
    Gfit<-cbind(ud.cens$time,ud.cens$surv)
    Gfit<-rbind(c(0,1),Gfit); 
    KMti<-Cpred(Gfit,time2)[,2];
    KMtimes<-Cpred(Gfit,times)[,2]; ## }}}
  } else if (cens.model=="cox") { ## {{{
    ud.cens<-coxph(Surv(time2,status==cens.code)~desX)
    aseout <- basehaz(ud.cens,centered=FALSE); 
    baseout <- cbind(baseout$time,baseout$hazard)
    Gcx<-Cpred(baseout,time2)[,2];
    RR<-exp(desX %*% coef(ud.cens))
    KMti<-exp(-Gcx*RR)
    KMtimes<-Cpred(Gfit,times)[,2]; 
    ## }}}
  } else if (cens.model=="aalen") {  ## {{{
    ud.cens<-aalen(Surv(time2,status==cens.code)~desX+cluster(clusters),n.sim=0,residuals=0,robust=0,silent=1)
    KMti <- Cpred(ud.cens$cum,time2)[,-1];
    Gcx<-exp(-apply(Gcx*desX,1,sum))
    Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
    Gfit<-rbind(c(0,1),cbind(time2,Gcx)); 
    KMti <- Gcx
    KMtimes<-Cpred(Gfit,times)[,2]; ## }}}
    } else  stop('Unknown censoring model') 
## }}}

###cat("Proportional odds model \n"); 
###dyn.load("Gprop-odds.so")

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
	as.integer(weighted.test),as.integer(ratesim),as.double(score),
	as.double(cumAi),as.double(cumAiiid),as.integer(residuals),
	as.double(loglike),as.integer(profile),as.integer(sym),
	as.double(KMtimes),as.double(KMti),as.double(time2),as.integer(causeS),
	as.integer(index-1),
	as.integer(baselinevar),as.integer(clusters),
       	as.integer(antclust), as.integer(cens.code), PACKAGE="timereg");

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
sim.supProp=sim.supUt,prop.odds=TRUE)

colnames(ud$cum)<-colnames(ud$var.cum)<- c("time","Baseline")
if (robust==1) colnames(ud$robvar.cum)<- c("time","Baseline"); 

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

attr(ud,"Call")<-sys.call(); 
attr(ud,"Formula")<-formula; 
attr(ud,"id")<-id.call; 
class(ud)<-"cox.aalen"
return(ud); 
}

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

