##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title Inverse Probability of Censoring Weights
##' @param formula 
##' @param data 
##' @param cluster 
##' @param samecens 
##' @param obsonly 
##' @param weightname 
##' @param pairs 
##' @param response 
##' @author Klaus K. Holst
##' @export
ipw <- function(formula,data,cluster,samecens=FALSE,obsonly=TRUE,weightname="w",
                pairs=FALSE,response) {

  timevar <- as.character(terms(formula)[[2]][[2]])
  otimes <- data[,timevar]
  utimes <- sort(unique(otimes))
  delta <- min(diff(c(0,utimes)))/2
  ## We want prediction just before event
  ##  if (length(attributes(terms(formula))$term.labels)) {    
  ##    fit <- cph(formula,data=data,surv=TRUE,x=TRUE,y=TRUE)
  ##    pr <- survest(fit,what="parallel",newdata=data,
  ##                    times=otimes-delta)
  ##  } else { ## cph does not work without covariates.. Kaplan-Meier:  
  ##  }
  ##  fit <- survfit(formula,data=data);
  ##  sfit <- summary(fit)
  ##  stratas <- fit$strata
  ##  if (is.null(stratas)) {    
  ##  Gfit <- cbind(fit$time,fit$surv)
  ##  pr <- fastapprox(Gfit[,1],otimes-delta,Gfit[,2])[[1]]
  ##    Gfit2<-rbind(c(0,1),Gfit); 
  ##    pr<-Cpred(Gfit2,otimes)[,2];
  ## } else {
  ##   for (s in stratas) {      
  ##   }
  ## }
  
  XZ <- model.matrix(formula,data)
  ud.cens<- aalen(formula,n.sim=0,robust=0,data=data);
  Gcx<-Cpred(ud.cens$cum,otimes)[,-1];
  Gcx<-exp(-apply(Gcx*XZ,1,sum))
  Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
  pr <- Gcx  
  
  noncens <- with(data,!eval(terms(formula)[[2]][[3]]))
  data[,weightname] <- pr  
  if (samecens & !missing(cluster)) {
    message("Minimum weights...")
    myord <- order(data[,cluster])
    data <- data[myord,,drop=FALSE]
    id <-  table(data[,cluster])
    if (pairs) {
      gem <- data[,cluster]%in%(names(id)[id==2])
      id <- id[id==2]
      data <- data[gem,]
    }
    d0 <- subset(data,select=c(cluster,weightname))
    noncens <- with(data,!eval(terms(formula)[[2]][[3]]))
    d0[,"observed."] <- noncens
    timevar <- paste("_",cluster,weightname,sep="")
    d0[,timevar] <- unlist(lapply(id,seq))
    Wide <- reshape(d0,direction="wide",timevar=timevar,idvar=cluster)
    W <- apply(Wide[,paste(weightname,1:2,sep=".")],1,
               function(x) min(x,na.rm=TRUE))
    Wmarg <- d0[,weightname]
    data[,weightname] <- 1/Wmarg
    Wmin <- rep(W,id)
    ##    d0[,weightname] <- 0

##################################################
############################################################
##################################################

    ##      Wcomb <- (Wmin-Wmarg)/(Wmarg*Wmin)      
    obs1only <- rep(with(Wide, observed..1 & (is.na(observed..2) | !observed..2)),id)
    obs2only <- rep(with(Wide, observed..2 & (is.na(observed..1) | !observed..1)),id)
    obsOne <- which(na.omit(obs1only|obs2only))
    obsBoth <- rep(with(Wide, !is.na(observed..1) & !is.na(observed..2) & observed..2 & observed..1),id)

    data[obsBoth,weightname] <-
      ifelse(noncens[obsBoth],1/Wmin[obsBoth],0)    
    data[obsOne,weightname] <-
      ifelse(noncens[obsOne],1/Wmarg[obsOne],0)
  } 
  if (obsonly)
    data <- data[noncens,,drop=FALSE]
  return(data)    
}
