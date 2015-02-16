##' @export
`lifetable` <- function(x,...) UseMethod("lifetable")

##' Create simple life table 
##' 
##' @title Life table
##' @param x time formula (Surv) or matrix/data.frame with columns time,status or entry,exit,status
##' @param strata Strata
##' @param data data.frame
##' @param breaks Time intervals
##' @param confint If TRUE 95\% confidence limits are calculated
##' @param ... Additional arguments to lower level functions
##' @author Klaus K. Holst
##' @aliases lifetable lifetable.matrix lifetable.formula
##' @usage
##'  \method{lifetable}{matrix}(x, strata = list(), breaks = c(),
##'    confint = FALSE, ...)
##'
##'  \method{lifetable}{formula}(x, data=parent.frame(), breaks = c(),
##'    confint = FALSE, ...)
##' @examples
##' library(timereg)
##' data(TRACE)
##' 
##' d <- with(TRACE,lifetable(Surv(time,status==9)~sex+vf,breaks=c(0,0.2,0.5,8.5)))
##' summary(glm(events ~ offset(log(atrisk))+factor(int.end)*vf + sex*vf,
##'             data=d,poisson))
##' @export
lifetable.matrix <- function(x,strata=list(),breaks=c(),confint=FALSE,...) {
    if (ncol(x)==3) {
        status <- x[,3]
        entry <- x[,1]
        time <- x[,2]
    } else {
        status <- x[,2]
        time <- x[,1]
        entry <- rep(0,length(time))
    }
    LifeTable(time,status,entry,strata,breaks,confint,...)
}

##' @export
lifetable.formula <- function(x,data=parent.frame(),breaks=c(),confint=FALSE,...) {
    cl <- match.call()
    mf <- model.frame(x,data)
    Y <- model.extract(mf, "response")
    Terms <- terms(x, data = data)
    if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- NULL ## rep(0,nrow(Y))
        status <- Y[,2]
    } else {
        entry <- Y[,1]
        exit <- Y[,2]
        status <- Y[,3]
    }
    strata <- list()
    X <- model.matrix(Terms, data)
    if (!is.null(intpos  <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)>0) {
        strata <- as.list(model.frame(Terms,data)[,-1,drop=FALSE])
    }
    LifeTable(exit,status,entry,strata,breaks,confint,...)       
}


LifeTable <- function(time,status,entry=NULL,strata=list(),breaks=c(),confint=FALSE,interval=TRUE,mesg=FALSE) {    
    if (is.null(entry)) entry <- rep(0,NROW(time))
    if (mesg) message(dim(time))
    if ((is.matrix(time) || is.data.frame(time)) && ncol(time)>1) {
        if (ncol(time)==3) {
            status <- time[,3]
            entry <- time[,1]
            time <- time[,2]
        } else {
            status <- time[,2]
            time <- time[,1]
            entry <- rep(0,length(time))
        }
    }
    if (length(strata)>0) {
        a <- by(cbind(entry,time,status), strata,
                FUN=LifeTable, breaks=breaks, confint=confint)
        cl <- lapply(strata,class)
        nulls <- which(unlist(lapply(a,is.null)))
        nonnulls <- setdiff(seq_along(a),nulls)
        nn <- do.call("expand.grid",attributes(a)$dimnames)
        if (length(nulls)>0) nn <- nn[-nulls,,drop=FALSE]
        nam <- nn[rep(seq(NROW(nn)),each=NROW(a[[nonnulls[1]]])),,drop=FALSE]
        xx <- list()
        for (i in seq(ncol(nam))) {
            if (cl[i]%in%c("numeric","integer"))
                xx <- c(xx,list(as.numeric(as.character(nam[,i]))))
            else                    
                xx <- c(xx, list(do.call(paste("as.",as.character(cl[i]),sep=""),list(nam[,i]))))
        }
        xx <- as.data.frame(xx); colnames(xx) <- colnames(nam)
        res <- Reduce("rbind",a)
        res <- cbind(res,xx)
        return(res)
    }
    if (length(breaks)==0) breaks <- c(0,max(time,na.rm=TRUE))
    if (length(breaks)==1) breaks <- c(0,breaks)
    breaks <- sort(unique(breaks))
    ## en <- matrix(unlist(lapply(c(-Inf, breaks),function(x) pmax(x,entry))),
    ##               ncol=length(breaks)+1)
    ## ex <- matrix(unlist(lapply(c(breaks, Inf),function(x) pmin(x,time))),
    ##               ncol=length(breaks)+1)
    en <- matrix(unlist(lapply(breaks[-length(breaks)],function(x) pmax(x,entry))),
                  ncol=length(breaks)-1)
    ex <- matrix(unlist(lapply(breaks[-1],function(x) pmin(x,time))),
                  ncol=length(breaks)-1)
    dur <- ex-en
    endur <- rbind(breaks[-1])%x%cbind(rep(1,nrow(en)))-en
    ##endur <- rbind(c(breaks,Inf))%x%cbind(rep(1,nrow(en)))-en
    dur[dur<0] <- NA
    enter <- colSums(!is.na(dur))
    atrisk <- colSums(dur,na.rm=TRUE)
    eventcens <- dur<endur
    eventcens <- rbind(apply(dur<endur,2,function(x) x*(status+1)))
    lost <- colSums(eventcens==1,na.rm=TRUE)
    events <- colSums(eventcens==2,na.rm=TRUE)
    res <- subset(data.frame(enter=enter,
                             atrisk=atrisk,
                             lost=lost,
                             events=events,
                             ## int.start=c(-Inf,breaks),
                             ## int.end=c(breaks,Inf),
                             int.start=breaks[-length(breaks)],
                             int.end=breaks[-1],
                             surv=0,
                             rate=events/atrisk))
    if (interval) res$interval <- factor(paste0("[",res$int.start,";",res$int.end,")"))
    cumsum.na <- function(x,...) { x[is.na(x)] <- 0; cumsum(x) }
    res$surv <- with(res, exp(-cumsum.na(rate*(int.end-int.start))))    
    if (confint) {
        ff <- events ~ offset(log(atrisk))
        if (length(breaks)>2) ff <- update(ff,.~.+factor(int.end)-1)
        g <- glm(ff,data=res,poisson)
        suppressMessages(ci <- rbind(exp(stats::confint(g))))
        res[,"2.5%"] <- ci[,1]
        res[,"97.5%"] <- ci[,2]
    }
    res
}



##' Summary for survival analyses via the 'lifetable' function
##'
##' Summary for survival analyses via the 'lifetable' function 
##' @title Extract survival estimates from lifetable analysis
##' @param object glm object (poisson regression)
##' @param ... Contrast arguments
##' @param timevar Name of time variable
##' @param time Time points (optional)
##' @param int.len Time interval length (optional)
##' @param confint If TRUE confidence limits are supplied
##' @param level Level of confidence limits
##' @param individual Individual predictions
##' @param length.out Length of time vector
##' @export
##' @author Klaus K. Holst
survpois <- function(object,...,timevar="int.end",time,int.len,confint=FALSE,level=0.95,individual=FALSE,length.out=25) {
    nn <- names(coef(object))
    timevar_re0 <- gsub("\\$|\\^","",glob2rx(timevar))
    timevar_re <- paste0(timevar_re0,"[0-9]+\\.*[0-9]*")
    idx <- regexpr(timevar_re,nn)
    dots <- list(...)
    if (all(idx<0)) {
        timevar_re0 <- gsub("\\$|\\^","",glob2rx(paste0("factor(",timevar,")")))
        timevar_re <- paste0(timevar_re0,"[0-9]+\\.*[0-9]*")
        idx <- regexpr(timevar_re,nn)
    }
    tvar <- unique(regmatches(nn,idx))
    if (missing(time)) {        
        i0 <- intersect(seq(nrow(object$data)),which(object$data[,"rate"]>0))
        for (i in seq_along(dots)) {
             i0 <- intersect(i0,which(object$data[,names(dots)[i]]==dots[[i]]))
        }
        time <- sort(unique(as.numeric(gsub(timevar_re0,"",tvar))))        
        rg <- range(object$data[i0,timevar],na.rm=TRUE)
        if (length(time)==0) {
            time <- seq(rg[1],rg[2],length.out=length.out)
        } else {
            time <- seq(1,rg[2])
        }
    }
    message(paste(time,collapse=","))
    if (missing(int.len)) {
        int.len <- diff(c(0,time))
    } else if (int.len==1) int.len <- rep(int.len,nrow(res))
    tt <- terms(object)
    ##if (attr(tt,"offset")) tt <- drop.terms(tt,dropx=attr(tt,"offset")-1)
    vv <- all.vars(formula(tt))
    dots[[timevar]] <- time
    if (individual) {
        newdata <- as.data.frame(dots)
        for (v in vv) {
            if (v%ni%names(dots)) newdata[,v] <- object$data[1,v]
        }        
    } else {
        for (v in vv) {
            if (v%ni%names(dots)) dots[[v]] <- object$data[1,v]
        }        
        args <- c(list(`_data`=model.frame(object)),dots)
        newdata <- do.call(Expand, args)
    }
    Terms <- delete.response(tt)    
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    ## NA-robust matrix product
    beta0 <- coef(object); beta0[is.na(beta0)] <- -.Machine$double.xmax    
    V <- vcov(object)    
    V0 <- structure(matrix(0,length(beta0),length(beta0)),dimnames=list(names(beta0),names(beta0)))
    idx <- which(rownames(V0)%in%rownames(V))
    V0[idx,idx] <- V    
    coefs <- X%*%beta0
    res <- cbind(time,exp(coefs))
    if (!confint) {
        res <- cbind(res,cumsum(res[,2]*int.len))
        res <- cbind(res,exp(-res[,3]))
    } else {
        S0 <- X%*%V0%*%t(X)
        system.time(e0 <- estimate(NULL,coef=coefs,vcov=S0,function(x) cumsum(exp(x)*int.len))) ## Cumulative hazard
        system.time(e1 <- estimate(e0, function(x) exp(-x),level=level)) ## Survival
        res <- cbind(res,e0$coefmat[,1],e1$coefmat[,c(1,3:4)])
    }
    colnames(res)[1:4] <- c("time","rate","chaz","surv") 
    structure(as.data.frame(cbind(res,newdata)),class=c("survpois","data.frame"))
}

##' @export
plot.survpois <- function(x,confint=TRUE,cont=TRUE,length.out=200,surv=TRUE,type=ifelse(surv,"l","s"),lty=c(1,2),add=FALSE,xlab="Time",ylab="Survival probability",xlim,ylim,...) {
    if (!add) {
        if (missing(xlim)) xlim <- c(0,max(x[,1]))
        if (missing(ylim)) {            
            ylim <- c(0,1)
            if (!surv) ylim[2] <- max(x$chaz)
        }
        plot(0,type="n",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,...)
    }
    if (cont) {
        tt <- seq(0,max(x$time),length.out=length.out)
        ff <- with(x, approxfun(c(0,time),c(0,chaz),method="linear"))
        if (surv) {
            lines(tt,exp(-ff(tt)),lty=lty[1],type=type,...)
        } else {
            lines(tt,ff(tt),lty=lty[1],type=type,...)
        }
        return(invisible(NULL))
    }
    if (surv) lines(x[,c(1,4)],lty=lty[1],type=type,...)
    else lines(x[,c(1,3)],lty=lty[1],type=type,...)
    if (ncol(x)>4 && confint && surv) {
        lines(x[,c(1,5)],lty=lty[2],type=type,...)
        lines(x[,c(1,6)],lty=lty[2],type=type,...)
    }
}


