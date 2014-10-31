##' Internal function.
##' Calculates Inverse Probability of Censoring
##' Weights (IPCW) and adds them to a data.frame
##'
##' @title Inverse Probability of Censoring Weights
##' @param formula Formula specifying the censoring model 
##' @param data data frame
##' @param cluster clustering variable
##' @param same.cens For clustered data, should same censoring be assumed (bivariate probability calculated as mininum of the marginal probabilities)
##' @param obs.only Return data with uncensored observations only
##' @param weight.name Name of weight variable in the new data.frame
##' @param cens.model Censoring model (default Aalens additive model)
##' @param pairs For paired data (e.g. twins) only the complete pairs are returned (With pairs=TRUE)
##' @param ... Additional arguments to censoring model 
##' @author Klaus K. Holst
##' @examples
##' data(prt)
##' prtw <- ipw(Surv(time,status==0)~country, data=prt[sample(nrow(prt),5000),],
##'             cluster="id",weight.name="w")
##' plot(0,type="n",xlim=range(prtw$time),ylim=c(0,1),xlab="Age",ylab="Probability")
##' count <- 0
##' for (l in unique(prtw$country)) {
##'     count <- count+1
##'     prtw <- prtw[order(prtw$time),]
##'     with(subset(prtw,country==l),
##'          lines(time,w,col=count,lwd=2))
##' }
##' legend("topright",legend=unique(prtw$country),col=1:4,pch=-1,lty=1)
ipw <- function(formula,data,cluster,
                same.cens=FALSE,obs.only=TRUE,
                weight.name="w",
                trunc.prob=FALSE,weight.name2="wt",
                cens.model="aalen", pairs=FALSE,
                theta.formula=~1, ...) {
                 ##iid=TRUE,

    ##cens.args <- c(list(formula,n.sim=0,robust=0,data=eval(data)),list(...))
    if (tolower(cens.model)%in%c("weibull","phreg.par","phreg.weibull")) {
	ud.cens <- phreg.par(formula,data,...)
        pr <- predict(ud.cens)
        noncens <- which(!ud.cens$status)        
    } else {
        m <- match.call(expand.dots = FALSE)
        m <- m[match(c("", "formula", "data", "subset", "na.action"), 
                     names(m), nomatch = 0)]
        special <- c("strata", "factor", "NN", "cluster", "dummy")
        if (missing(data)) 
            Terms <- terms(formula)
        else Terms <- terms(formula, data = data)
        m$formula <- Terms
        m[[1]] <- as.name("model.frame")
        M <- eval(m, parent.frame())
        censtime <- model.extract(M, "response")
        if (ncol(censtime)==3) {
            status <- censtime[,3]
            otimes <- censtime[,2]
            ltimes <- censtime[,1]
        } else {
            status <- censtime[,2]
            otimes <- censtime[,1]
        }
        noncens <- !status
        if (is.null(attr(terms(formula,"prop"),"specials")$prop)) {
            ud.cens <- aalen(formula,n.sim=0,robust=0,data=data,...)
            XZ <- model.matrix(formula,data)
            Gcx <- ud.cens$cum[prodlim::sindex(ud.cens$cum[,1],otimes),-1]
            ##Gcx<-Cpred(ud.cens$cum,otimes)[,-1];            
            Gcx<-exp(-apply(Gcx*XZ,1,sum))            
        } else {
            ud.cens <- cox.aalen(formula,n.sim=0,robust=0,data=data,...)
            XZ <- model.matrix(formula,data)
            Gcx<-Cpred(ud.cens$cum,otimes)[,-1];
            Gcx<-exp(-apply(Gcx*XZ,1,sum))            
        }
        ##ud.cens <- do.call(cens.model,cens.args)
        Gcx[Gcx>1]<-1; Gcx[Gcx<0]<-0
        pr <- Gcx
    }
    if (trunc.prob & ncol(censtime)==3) { ## truncation
        browser()
        data$truncsurv <- Surv(ltimes,otimes,noncens)
        trunc.formula <- update(formula,truncsurv~.)        
        ud.trunc <- aalen(trunc.formula,data=data,robust=0,n.sim=0,residuals=0,silent=1,max.clust=NULL,
                          clusters=data[,cluster], ...)
        dependX0 <- model.matrix(theta.formula,data)        
        twostage.fit <-two.stage(ud.trunc,
                                 data=data,robust=0,detail=0,
                                 theta.des=dependX0)#,Nit=20,step=1.0,notaylor=1)
        X <- model.matrix(trunc.formula,data)
        Xnam <- colnames(X)
        ww <- fast.reshape(cbind(X,".num"=seq(nrow(X)),".lefttime"=ltimes),varying=c(".num",".lefttime"),id=data[,cluster])
        dependX <- model.matrix(theta.formula,ww)                
        Prob <- predict.two.stage(twostage.fit,X=ww[,Xnam],
                                  times=ww[,".lefttime1"],times2=ww[,".lefttime2"],
                                  theta.des=dependX)
        P0 <- numeric(nrow(X))
        P0[ww[,".num1"]] <- Prob$St1t2
        P0[ww[,".num2"]] <- Prob$St1t2
        data[,weight.name2] <- P0
    }
    data[,weight.name] <- pr
        
        
    if (same.cens & !missing(cluster)) {        
        message("Minimum weights...")
        myord <- order(data[,cluster])
        data <- data[myord,,drop=FALSE]
        id <-  table(data[,cluster])
        if (pairs) {
            gem <- data[,cluster]%in%(names(id)[id==2])
            id <- id[id==2]
            data <- data[gem,]
        }
        d0 <- subset(data,select=c(cluster,weight.name))
        noncens <- with(data,!eval(terms(formula)[[2]][[3]]))
        d0[,"observed."] <- noncens
        timevar <- paste("_",cluster,weight.name,sep="")
        d0[,timevar] <- unlist(lapply(id,seq))
        Wide <- reshape(d0,direction="wide",timevar=timevar,idvar=cluster)
        W <- apply(Wide[,paste(weight.name,1:2,sep=".")],1,
                   function(x) min(x,na.rm=TRUE))
        Wmarg <- d0[,weight.name]
        data[,weight.name] <- 1/Wmarg
        Wmin <- rep(W,id)
        
        obs1only <- rep(with(Wide, observed..1 & (is.na(observed..2) | !observed..2)),id)
        obs2only <- rep(with(Wide, observed..2 & (is.na(observed..1) | !observed..1)),id)
        obsOne <- which(na.omit(obs1only|obs2only))
        obsBoth <- rep(with(Wide, !is.na(observed..1) & !is.na(observed..2) & observed..2 & observed..1),id)
        
        data[obsBoth,weight.name] <-
            ifelse(noncens[obsBoth],1/Wmin[obsBoth],0)    
        data[obsOne,weight.name] <-
            ifelse(noncens[obsOne],1/Wmarg[obsOne],0)
    }

    if (obs.only)
        data <- data[noncens,,drop=FALSE]
    
    return(data)    
}
