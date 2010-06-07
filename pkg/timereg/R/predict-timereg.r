slaaop<-function(z,time,cum)
{x<-time; y<-cum; index<-sum(x<=z);if (index==0) index<-1;
retur<-y[index]; return(retur); }

pred.cum<-function(x,time,cum) {ud<-sapply(x,slaaop,time,cum); 
 return(ud)}

"pred.des"<-function(formula,data=sys.parent())
{ ## {{{
  call <- match.call();
  m <- match.call(expand=FALSE);
  special <- c("const")
  Terms <- if(missing(data)) terms(formula, special)
           else          terms(formula, special,data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, data)
  mt <- attr(m, "terms")

  XZ<-model.matrix(Terms,m)[,drop = FALSE]
  cols<-attributes(XZ)$assign
  l.cols<-length(cols)
  semicov <- attr(Terms, "specials")$const
  ZtermsXZ<-semicov-1

  if (length(semicov)) {renaalen<-FALSE; Zterms<-c();
  for (i in ZtermsXZ) Zterms<-c(Zterms,(1:l.cols)[cols==i]);
   } else {renaalen<-TRUE;}

 if (length(semicov)) {
  X<-as.matrix(XZ[,-Zterms]);
  #covnamesX <- dimnames(XZ)[[2]][-Zterms]; dimnames(X)[[2]]<-covnamesX;
  Z<-as.matrix(XZ[,Zterms]);
  #covnamesZ <- dimnames(XZ)[[2]][Zterms];dimnames(Z)[[2]]<-covnamesZ;          
 }
  else {X<-as.matrix(XZ); #covnamesX <- dimnames(XZ)[[2]];
        Z<-FALSE; #dimnames(X)[[2]]<-covnamesX; 
  }
  X <- data.matrix(X); 
  return(list(covarX=X,covarZ=Z))
} ## }}}

predict.aalen <-  function(object,...){

  predict.comprisk(object,...)

}

predict.cox.aalen <-  function(object,...){

  predict.comprisk(object,...)

}

aalen.des2 <-  function(formula,data=sys.parent(),model=NULL,...){
  call <- match.call()
  m <- match.call(expand=FALSE)
  m$model <- NULL
###  special <- c("const","prop")
  Terms <- if(missing(data)) terms(formula )
  else              terms(formula, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
 intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")

  if (model=="cox.aalen") modela <- "cox.aalen" else modela <- "aalen"
  des<-read.design(m,Terms,model=modela)
  return(des) 
}

###pred2<-function(object,newdata=NULL,X=NULL,
predict.comprisk<-function(object,newdata=NULL,X=NULL,
                           Z=NULL,n.sim=500, uniform=TRUE,
                           se=TRUE,alpha=0.05,...){
## {{{
  if (!(inherits(object,'comprisk') || inherits(object,'aalen')
        || inherits(object,'cox.aalen')))
    stop ("Must be output from comp.risk function")

  if(inherits(object,'aalen')) { modelType <- 'aalen';
  } else if(inherits(object,'comprisk')) { modelType <- object$model;
  } else if(inherits(object,'cox.aalen')) { modelType <- 'cox.aalen'; }
  n <- length(object$B.iid) ## Number of clusters (or number of individuals
                            ## if no cluster structure is specified)

  if (is.null(object$B.iid)==TRUE & (se==TRUE | uniform==TRUE)) {
    stop("resample processes necessary for these computations, set resample.iid=1");
  }
  if (is.null(object$gamma)==TRUE) { semi<-FALSE } else { semi<-TRUE }

  ## {{{ extracts design based on the different specifications
  ## cox.aalen uses prop(...), while aalen and comp.risk use const(...)
  ## this accounts for the different number of characters
  if(inherits(object,'cox.aalen')){ indexOfFirstChar <- 6; } else { indexOfFirstChar <- 7; }
  
  constant.covs <- NULL
  if (!is.null(newdata)) {
    ##  The time-constant effects first
    formulao <- attr(object,"Formula")
    des <- aalen.des2(formula(delete.response(terms(formulao))),
		    data=newdata,model=modelType)
    time.vars <- des$X 

    if (semi==TRUE) {
      constant.covs <- des$Z 
      const <- c(object$gamma)
      names(const) <-substr(dimnames(object$gamma)[[1]],indexOfFirstChar,
                            nchar(dimnames(object$gamma)[[1]])-1)
    }

    ## Then extract the time-varying effects
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
    fittime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    nobs <- nrow(newdata)
  } else if ((is.null(Z)==FALSE) || (is.null(X)==FALSE)){

  if(is.null(Z) || is.matrix(Z) || is.data.frame(Z)){
    constant.covs<-Z;
    } else { constant.covs<-matrix(Z,ncol = nrow(object$gamma)); }

  if ((is.null(Z)==FALSE) && is.null(X))
  {
   if (ncol(object$cum)==2) X<-time.vars<-matrix(1,nrow(Z),1) else
   stop("When X is not specified we assume that it an intercept terms\n");
  }
      
    ## Allow for the possibility that X is a vector 
    if(!is.matrix(X) && !is.data.frame(X)){
      time.vars<-matrix(X,ncol = ncol(object$cum)-1);
    } else {
      time.vars<-X;
    }
    nobs<-nrow(time.vars);
    ## extract the time-varying effects
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
  } else {
    stop("Must specify either newdata or X, Z\n");
  }

  ## }}}

  ## {{{ predictions for competing risks and survival data

  cumhaz<-as.matrix(time.vars) %*% t(time.coef[,-1])
  time<-time.coef[,1]; if (semi==TRUE) pg <- ncol(constant.covs); 
  nt<-length(time);

  ### set up articial time.pow for aalen  and cox.aalen to make unified code for 
  ###  comp.risk and survival
  if(inherits(object,'aalen') & semi==TRUE) timepow <- rep(1,pg)
  if(inherits(object,'cox.aalen')) timepow <- rep(0,pg)
  if(inherits(object,'comprisk')) timepow <- attr(object,"time.pow")

  if (semi==TRUE)
  constant.part <- constant.covs %*% (matrix( rep(c(time),pg)^timepow,pg,nt) *c(object$gamma) )

  if (inherits(object,'comprisk')) { ## {{{ competing models
    if (modelType == "additive") {
      if (semi==FALSE){
        P1=1-exp(-cumhaz);
      } else {
        P1=1-exp(-cumhaz-constant.part )
      }
      RR<-NULL; 
    } else if (modelType == 'rcif') { # P1=exp(x^T b(t) + z^t t^p gamma) 
        if (semi==FALSE){
           P1=exp(cumhaz);
         } else {
         P1<-exp(cumhaz+constant.part);
       }
       RR<-1;
    } else if (modelType == 'prop') {# model proportional 
        if (semi==FALSE){
        RR<-exp(cumhaz);
      } else {
        RR<-exp(cumhaz+constant.part);
    }
    P1<-1-exp(-RR);
    } else if (modelType == 'logistic') { #model logistic
      if (semi==FALSE){ RR<-exp(cumhaz); }   else { RR<-exp(cumhaz+constant.part); }
      P1<-RR/(1+RR);
    } ## }}}
    } else {  # survival model  ## {{{
       if (modelType == "aalen") {    #Aalen model
          if (semi==FALSE){ S0=exp(-cumhaz); } else { S0=exp(-cumhaz-constant.part) }
       RR<-NULL; 
       } else if(modelType == 'cox.aalen'){  #Cox-Aalen model
       if(semi == FALSE){ RR <- NULL; S0 <- exp(-cumhaz); } else {
         RR <- exp(constant.part);
         S0 <- exp(-cumhaz * RR);
       }
     }
    } ## }}}

    ## }}}

  se.P1 <- NULL
  se.S0 <- NULL
  ## i.i.d decomposition for computation of standard errors  ## {{{
  if (se==1) {
    pg<-length(object$gamma); 
    delta<-c();
    for (i in 1:n) {
       tmp<- as.matrix(time.vars) %*% t(object$B.iid[[i]]) 

       if (semi==TRUE) {
###       if (inherits(object,'comprisk')) {
             gammai <- matrix(object$gamma.iid[i,],pg,1); 
             tmp.const<-constant.covs %*% (matrix( rep(time,pg)^timepow,pg,nt)*c(gammai) )
###          if (semi==TRUE) { tmp.const <- constant.covs %*% gammai; }
       } 
###       else { }

      if (i==0) {
        print(tmp.const);
        if (modelType=="additive" || modelType == 'aalen'){ 
          print(tmp.const %*% matrix(time,1,nt))
        } else if (modelType=="prop"){
          print(tmp.const %*% matrix(1,1,nt));
        } else if (modelType=="cox.aalen") {
          tmp <- RR * tmp + RR * cumhaz * matrix(tmp.const,nobs,nt);
        }
      }

      if (semi==TRUE){
        if(modelType=="additive" || modelType == "aalen") {
           # || modelType=="rcif") {
          tmp<-tmp+ tmp.const ## %*% matrix(time,1,nt)
        } else if (modelType=="prop" || modelType=="rcif") {
###          tmp<-RR*tmp+RR*matrix(tmp.const,nobs,nt);
          tmp<-RR*tmp+RR*tmp.const;
	  ## modification of jeremy's code RR
        } else if (modelType=="cox.aalen") {
###          tmp <- RR * tmp + RR * cumhaz * matrix(tmp.const,nobs,nt);
          tmp <- RR * tmp + RR * cumhaz * tmp.const
        }
      }
      delta<-cbind(delta,c(tmp)); 
    }
    se<-apply(delta^2,1,sum)^.5
    if(modelType == 'additive' || modelType == 'prop'){
      se.P1<-matrix(se,nobs,nt)*(1-P1) 
    } 
    else if(modelType == 'rcif'){
      se.P1<-matrix(se,nobs,nt)*(P1) 
    } 
    else if (modelType == 'logistic'){
      se.P1<-matrix(se,nobs,nt)*P1/(1+RR)
    } 
    else if (modelType == 'aalen' || modelType == 'cox.aalen'){
      se.S0<-matrix(se,nobs,nt)*S0
    }
    ## }}}

    ### uniform confidence bands, based on resampling  ## {{{
    if (uniform==1) {
      mpt <- .C('confBandBasePredict',
                delta = as.double(delta),
                nObs = as.integer(nobs),
                nt = as.integer(nt),
                n = as.integer(n),
                se = as.double(se),
                mpt = double(n.sim*nobs),
                nSims = as.integer(n.sim),
		PACKAGE="timereg")$mpt;
  
      mpt <- matrix(mpt,n.sim,nobs,byrow = TRUE);
      uband <- apply(mpt,2,percen,per=1-alpha);
    } else uband<-NULL; 
  } else {
    uband<-NULL;
  }

  if(modelType == 'additive' || modelType == 'prop' || modelType=="logistic"
     || modelType=='rcif'){
    P1<-matrix(P1,nrow=nobs);
  } else if (modelType == 'aalen' || modelType == 'cox.aalen'){
    S0<-matrix(S0,nrow=nobs);
  }

  out<-list(time=time,unif.band=uband,model=modelType,alpha=alpha,
            newdata=list(X = time.vars, Z = constant.covs),RR=RR,
            call=sys.calls()[[1]], initial.call = attr(object,'Call'));
  if(modelType == 'additive' || modelType == 'prop' || modelType=="logistic"
     || modelType=='rcif'){
    out$P1 <- P1;
    out$se.P1 <- se.P1;    
  } else if (modelType == 'aalen' || modelType == 'cox.aalen'){
    out$S0 <- S0;
    out$se.S0 <- se.S0;    
  }
   # e.g. for an compound risk model, className = predictComprisk
  className <- switch(class(object),aalen='predictAalen',cox.aalen='predictCoxAalen',comprisk='predictComprisk')
  class(out) <- className

  return(out)
} ## }}}


plot.predictAalen <-  function(x,...){

  plot.predictComprisk(x,...)

}

plot.predictCoxAalen <-  function(x,...){

  plot.predictComprisk(x,...)

}

pava = function(x, w=rep(1,length(x)))  # R interface to the compiled code
{ ## {{{
  n = length(x)
  if (n != length(w)) return (0)    # error
  result  = .C("pava",
        y = as.double(x),
        as.double(w),
        as.integer(n) )
  result[["y"]]
} ## }}}


plot.predictComprisk<-function(x,uniform=1,new=1,se=1,col=1,lty=1,lwd=2,multiple=0,specific.comps=0,
xlab="Time",ylab="Probability",transparency=FALSE,monotone=TRUE,...){
## {{{
  object <- x; rm(x);
  modelType <- object$model;
  time<-object$time;
  uband<-object$unif.band;
  nobs<-nrow(object$newdata$X);
  RR<-object$RR;
  alpha <- object$alpha;
  ### Here we use mainLine as the central line (between confidence
  ### intervals or bands), so that we don't have to distinguish
  ### between the case when we want to plot a predicted survival function
  ### and the case when we want to plot a predicted risk funcion
  
  if (modelType == 'aalen' || modelType == 'cox.aalen'){
    type<-"surv"
    mainLine <- object$S0;
    if (monotone==TRUE) { mainLine<--t(apply(as.matrix(-mainLine),1,pava)); 
    mainLine[mainLine<0]<-0; 
    mainLine[mainLine>1]<-1; 
    }
    mainLine.se <- object$se.S0;    
  } else if(modelType == 'additive' || modelType == 'prop' || modelType=="logistic" || modelType=="rcif"){
    type<-"cif"
    mainLine <- object$P1;
    if (monotone==TRUE) { mainLine<-t(apply(as.matrix(mainLine),1,pava)); 
                           mainLine[mainLine<0]<-0; 
                           mainLine[mainLine>1]<-1; 
    }
    mainLine.se <- object$se.P1;    
  }
  
  if (length(col)!=nobs){ col<-rep(col[1],nobs); }
  if (length(lty)!=nobs){ lty<-rep(lty[1],nobs); }
  if (length(lwd)!=nobs){ lwd<-rep(lwd[1],nobs); }
  if (sum(specific.comps)==0){
    comps<-1:nobs
  } else {
    comps<-specific.comps
  }

  for (i in comps) {
    if (new==1 & (multiple!=1 | i==comps[1])) {
      plot(time,mainLine[i,],type="s",ylim=c(0,1),xlab=xlab,ylab=ylab,col=col[i],lty=lty[i],lwd=lwd[i],...)
    } else {
      lines(time,mainLine[i,],type="s",col=col[i],lty=lty[i],lwd=lwd[i])
    }

    if (se==1 & is.null(mainLine.se)==FALSE ) {
      lower<-mainLine[i,]-qnorm(1-alpha/2)*mainLine.se[i,]
      upper<-mainLine[i,]+qnorm(1-alpha/2)*mainLine.se[i,]
       if (monotone==TRUE) { 
       if (type=="cif") { lower<- pava(lower); upper<- pava(upper); }
       if (type=="surv") { lower<- -pava(-lower); upper<- -pava(-upper); }
        lower[lower<0]<-0; lower[lower>1]<-1; 
        upper[upper<0]<-0; upper[upper>1]<-1; 
       }

      lines(time,lower,type="s",col=col[i],lty=3,lwd=lwd[i]/2);
      lines(time,upper,type="s",col=col[i],lty=3,lwd=lwd[i]/2);
    }

    if (uniform==1 & is.null(uband)==FALSE ) {
      #if (level!=0.05) c.alpha<-percen(object$sim.test[,i],1-level)
      #else c.alpha<-object$conf.band.cumz[i];
      c.alpha=uband[i]; 
      upper<-mainLine[i,]-uband[i]*mainLine.se[i,];
      lower<-mainLine[i,]+uband[i]*mainLine.se[i,];
       if (monotone==TRUE) { 
          if (type=="cif") { lower<- pava(lower); upper<- pava(upper); }
          if (type=="surv") { lower<- -pava(-lower); upper<- -pava(-upper); }
          lower[lower<0]<-0; lower[lower>1]<-1; 
          upper[upper<0]<-0; upper[upper>1]<-1; 
       }
      if (transparency==0 || transparency==2) {
      lines(time,upper,type="s",col=col[i],lty=2,lwd=lwd[i]/2);
      lines(time,lower,type="s",col=col[i],lty=2,lwd=lwd[i]/2);
      }

    ## Prediction polygons bandds ## {{{
    if (transparency>=1) {
     col.alpha<-0.2
     col.ci<-"darkblue"
     col.ci<-col[i]; 
     lty.ci<-2
      if (col.alpha==0) col.trans <- col.ci
      else
      col.trans <- sapply(col.ci, FUN=function(x) do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))

      #print(t); print(ci)
      n<-length(time)
      tt<-seq(time[1],time[n],length=n*10); 
      ud<-Cpred(cbind(time,upper,lower),tt)[,2:3]
      tt <- c(tt, rev(tt))
      yy <- c(upper, rev(lower))
#      tt <- c(time, rev(time))
#      yy <- c(upper, rev(lower))
     yy <- c(ud[,1], rev(ud[,2]))
      polygon(tt,yy, col=col.trans, lty=0)      
  } ## }}}

    }
  }
} ## }}}

print.predictAalen <- function(x,...){
  print.predictComprisk(x,...)
}

print.predictCoxAalen <- function(x,...){
  print.predictComprisk(x,...)
}

print.predictComprisk <- function(x,...){

  object <- x; rm(x);
  
  if(!(inherits(object,'predictAalen') ||
       inherits(object,'predictCoxAalen') ||
       inherits(object,'predictComprisk'))){
    stop('Wrong class of object');
  }

  if (is.null(object$newdata$Z)==TRUE) semi<-FALSE else semi<-TRUE
  
  modelType <- object$model;
  modelAnnouncement <- ' Predicted survival for'
  addTo <- switch(modelType,
                  cox.aalen = 'a Cox-Aalen',
                  aalen = 'an Aalen',
                  prop = 'a proportional competing risks',
                  additive = 'an additive competing risks')
  modelAnnouncement <- paste(modelAnnouncement,addTo,'model',sep = ' ')
  cat(modelAnnouncement,fill=TRUE)

  cat(" Nonparametric terms : "); cat(colnames(object$newdata$X)[-1]); cat("   \n");  
  if (semi == TRUE) {
    cat(" Parametric terms :  "); cat(rownames(object$newdata$Z)); 
    cat("   \n");  } 
  cat("   \n");  
  
  call <- object$call;
  cat('Call to predict:',fill=TRUE);
  print(call)
  call <- object$initial.call;
  cat('Initial call:',fill=TRUE);
  print(call)
    
}

summary.predictAalen <- function(object,...){
  summary.predictComprisk(object,...)
}

summary.predictCoxAalen <- function(object,...){
  summary.predictComprisk(object,...)
}

summary.predictComprisk <- function(object,...){
  if(!(inherits(object,'predictAalen') ||
       inherits(object,'predictCoxAalen') ||
       inherits(object,'predictComprisk'))){
    stop('Wrong class of object');
  }

  modelClass <- class(object)
  modelType <- object$model;
  time<-object$time;
  uband<-object$unif.band;
  nobs<-nrow(object$newdata$X);
  RR<-object$RR;
  alpha <- object$alpha;
  call <- object$call;
  if (modelType == 'aalen' || modelType == 'cox.aalen'){
    se <- object$se.S0;    
  } else if(modelType == 'additive' || modelType == 'prop'){
    se <- object$se.P1;    
  }
  
  modelAnnouncement <- 'Predicted survival for'
    addTo <- switch(modelType,
                    cox.aalen = 'a Cox-Aalen',
                    aalen = 'an Aalen',
                    prop = 'a proportional competing risks',
                    additive = 'an additive competing risks')
  modelAnnouncement <- paste(modelAnnouncement,addTo,'model',sep = ' ')
  cat(modelAnnouncement,fill=TRUE)
  timeStatement <- paste('At',length(time),'times:',paste(c(head(time),''),collapse = ', '),'...,',time[length(time)])
  cat(timeStatement,fill=TRUE)
  obsStatement <- paste('Given covariates for',nobs,'new observations');
  cat(obsStatement,fill=TRUE)
  if(is.null(se)){
    addTo <- " - not yet done";
  } else if(is.null(uband)){
    addTo <- " - only pointwise calculations have been done"
  } else {
    addTo <- " - pointwise CI and uniform confidence band available"
  }
  cat('Standard error calculations:',fill=TRUE);
  cat(addTo,fill=TRUE);
  cat('Call:',fill=TRUE);
  print(call)
  
}

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
        nl<-B[,v]-1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
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


"summary.comprisk" <- function (object,digits = 3,...) {
  if (!inherits(object, 'comprisk')) stop ("Must be a comprisk object")
  
  if (is.null(object$gamma)==TRUE) semi<-FALSE else semi<-TRUE
    
  # We print information about object:  
  cat("Competing risks Model \n\n")
  
  modelType<-object$model
  #if (modelType=="additive" || modelType=="rcif") 
 
  if (sum(object$obs.testBeq0)==FALSE) cat("No test for non-parametric terms\n") else
  timetest(object,digits=digits); 

  if (semi) { cat("Parametric terms : \n"); coef(object); cat("   \n"); }

  cat("  Call: \n")
  dput(attr(object, "Call"))
  cat("\n")
}

