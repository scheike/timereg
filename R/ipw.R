##' Calculates Inverse Probability of Censoring Weights (IPCW) and
##' adds them to a data.frame
##'
##' @title Inverse Probability of Censoring Weights
##' @param formula Formula specifying the censoring model 
##' @param data data frame
##' @param cluster clustering variable
##' @param samecens For clustered data, should same censoring be assumed (bivariate probability calculated as mininum of the marginal probabilities)
##' @param obsonly Return data with uncensored observations only
##' @param weightname Name of weight variable in the new data.frame
##' @param cens.model Censoring model (default Aalens additive model)
##' @param pairs For paired data (e.g. twins) only the complete pairs are returned (With pairs=TRUE)
##' @param ... Additional arguments to censoring model 
##' @author Klaus K. Holst
##' @examples
##' data(prt)
##' prtw <- ipw(Surv(time,status==0)~country, data=prt[sample(nrow(prt),5000),],
##'             cluster="id",weightname="w")
##' plot(0,type="n",xlim=range(prtw$time),ylim=c(0,1),xlab="Age",ylab="Probability")
##' count <- 0
##' for (l in unique(prtw$country)) {
##'     count <- count+1
##'     prtw <- prtw[order(prtw$time),]
##'     with(subset(prtw,country==l),
##'          lines(time,w,col=count,lwd=2))
##' }
##' legend("topright",legend=unique(prtw$country),col=1:4,pch=-1,lty=1)
##' @export
ipw <- function(formula,data,cluster,
                 samecens=FALSE,obsonly=TRUE,weightname="w",
                 cens.model="aalen", pairs=FALSE, ...) {
                 ##iid=TRUE,

    XZ <- model.matrix(formula,data)
    cens.args <- c(list(formula,n.sim=0,robust=0,data=eval(data)),list(...))
    if (tolower(cens.model)%in%c("weibull","phreg.par","phreg.weibull")) {
        ud.cens <- do.call(mets::phreg.par,cens.args)
        pr <- predict(ud.cens)
        noncens <- which(!ud.cens$status)        
    } else {
        m <- match.call(expand.dots = TRUE)[1:3]
        Terms <- terms(formula, data = cens.args$data, env=parent.frame())
        m$formula <- Terms
        m$data <- cens.args$data
        m[[1]] <- as.name("model.frame")
        M <- eval(m)
        censtime <- model.extract(M, "response")
        status <- censtime[,2]
        otimes <- censtime[,1]
        noncens <- !status
        ud.cens <- do.call(cens.model,cens.args)
        Gcx<-Cpred(ud.cens$cum,otimes)[,-1];
        Gcx<-exp(-apply(Gcx*XZ,1,sum))
        Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
        pr <- Gcx
    }
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
