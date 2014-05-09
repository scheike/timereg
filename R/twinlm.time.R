## d <- twinsim(1000,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
## b <- twinlm.strata(y~x1,data=d,id="id",zyg="zyg",DZ="DZ",var="x2",quantiles=c(0,0.25,0.5,0.75,1))
## plot(b,which=5:6,ylim=c(0,1),type="l")
## plot(b,which=5:6,ylim=c(0,1),col=c("darkred","darkblue"),legend=c("MZ","DZ"),lty=1:2)

##' @export
twinlm.strata <- function(formula,data,var,breaks,quantiles,...) {
    varname <- ""
    if (is.character(var) && length(var)==1) {
        varname <- var
        var <- data[,var]
    }
    ## 
    if (!missing(quantiles)) breaks <- quantile(var,quantiles,include.lowest=TRUE)
    if (!missing(breaks)) {
        var <- cut(var,c(breaks))
    }
    if (!inherits(var,c("character","factor","logical"))) stop("Expected factor, character, or logical vector")
    lev <- na.omit(unique(var))
    res <- list()
    for (i in seq_along(lev)) {
        tau <- lev[i]
        if (length(breaks)>1) message(tau)
        data0 <- data[var==tau,,drop=FALSE]
        suppressWarnings(b <- twinlm(formula,data=data0,...))
        res <- c(res,list(summary(b)))
    }
    if (length(breaks)==1) return(b)
    coef <- c(lapply(res,function(x) x$all),list(res[[length(res)]]$all))
    res <- list(varname=varname,var=breaks,coef=coef,summary=res,type="strata")
    class(res) <- "multitwinlm"
    return(res)
}


## data(prt)
## bb <- twinlm.time(cancer~country,data=prt,id="id",zyg="zyg",DZ="DZ",cens.formula=Surv(time,status==0)~zyg,breaks=seq(70,90,by=4))
## plot(bb,which=c(7,11),ylim=c(0,28),legendpos="topright",col=c("darkred","darkblue"),lty=c(1,2),legend=c("MZ","DZ"),ylab="Relative recurrence risk ratio")
## plot(bb,which=c(7,11),ylim=c(0,28),legendpos="topright",col=c("darkred","darkblue"),lty=c(1,2),legend=c("MZ","DZ"),ylab="Relative recurrence risk ratio",type="l")

##' @export
twinlm.time <- function(formula,data,id,type="u",...,
                        breaks=Inf,
                        cens.formula,cens.model="aalen",weight="w") {

    m <- match.call(expand.dots = TRUE)[1:3]
    Terms <- terms(cens.formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m)
    censtime <- model.extract(M, "response")
    status <- censtime[,2]
    time <- censtime[,1]
    outcome <- as.character(terms(formula)[[2]])    
    if (is.null(breaks)) breaks <-  quantile(time,c(0.25,0.5,0.75,1))

    outcome0 <- paste(outcome,"_dummy")
    res <- list()
    for (tau in breaks) {
        if (length(breaks)>1) message(tau)
        data0 <- data
        time0 <- time
        cond0 <- time0>tau
        status0 <- status
        status0[cond0 & status==1] <- 3
        data0[cond0,outcome] <- FALSE
        time0[cond0] <- tau
        data0$S <- Surv(time0,status0==1)        
        dataw <- ipw(S~zyg, data=data0,
                     cluster=id,weightname=weight,pairs=TRUE)
        b <- bptwin(formula, data=dataw, id=id, weight=weight,type=type,...)
        res <- c(res,list(summary(b)))
    }
    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,coef=lapply(res,function(x) x$all),summary=res,call=m,type="time")
    class(res) <- "multitwinlm"
    return(res)    
}

##' @S3method summary multitwinlm
summary.multitwinlm <- function(object,which=seq(nrow(object$coef[[1]])),...) {
    res <- list()
    for (i in which) {    
        rr <- matrix(unlist(lapply(object$coef,function(z) z[i,])),ncol=3,byrow=TRUE)
        colnames(rr) <- colnames(object$coef[[1]])
        rr <- cbind(object$var,rr)
        colnames(rr)[1] <- object$varname
        res <- c(res,list(rr))
    }
    names(res) <- rownames(object$coef[[1]])[which]
    return(res)
}


##' @S3method print multitwinlm
print.multitwinlm <- function(x,...) {
    res <- summary(x,...)
    for (i in seq_along(res)) {
        cat(i, ": ", names(res)[i], "\n",sep="")
        printCoefmat(res[[i]],...)
    }
    invisible(x)
}

##' @S3method plot multitwinlm
plot.multitwinlm <- function(x,...,which=1,
                            type="s",
                            lwd=2,lty=1,col,fillcol,alpha=0.2,
                            xlab=x$varname,
                            ylab="",idx=seq_along(x$var),
                            lasttick=TRUE,
                            legend=TRUE,legendpos="topleft") {
    ss <- summary(x,which)
    add <- FALSE
    if (missing(col)) col <- seq_along(which)
    if (length(col)==1) col <- rep(col,length(which))
    if (length(lwd)==1) lwd <- rep(lwd,length(which))
    if (length(lty)==1) lty <- rep(lty,length(which))
    if (alpha>0 & missing(fillcol)) fillcol <- Col(col,alpha)
    count <- 0
    for (tt in seq_along(which)) {
        count <- count+1
        zz <- ss[[tt]][idx,,drop=FALSE]
        if (!add) {    
            plot(zz[,1:2,drop=FALSE],type=type,lwd=lwd[count],lty=lty[count],
                 ylab=ylab,xlab=xlab,col=col[count],...)
            dev <- devcoords()
        }
        if (lasttick && type=="s") {
            zz2 <- rbind(zz,tail(zz,1))
            zz2[nrow(zz2),1] <- dev$fig.x2
            confband(zz2[,1],zz2[,3],zz2[,4],polygon=TRUE,step=(type=="s"),col=fillcol[count],border=0)
        } else {
            confband(zz[,1],zz[,3],zz[,4],polygon=TRUE,step=(type=="s"),col=fillcol[count],border=0)
        }
        add <- TRUE
        ## if (lasttick && type=="s") {
        ##     axis(4,at=zz[nrow(zz),3],labels=FALSE,
        ##          lwd=lwd[count],col=fillcol[count],tcl=.5,...)
        ##     axis(4,at=zz[nrow(zz),4],labels=FALSE,
        ##          lwd=lwd[count],col=fillcol[count],tcl=.5,...)
        ## }
        lines(zz[,1:2,drop=FALSE],lwd=lwd[count],lty=lty[count],col=col[count],type=type,...)        
    }
    if (!is.null(legend) || !legend[1]) {
        if (is.logical(legend) || length(legend)==1) legend <- rownames(x$coef[[1]])[which]
        graphics::legend(legendpos,legend=legend,col=col,lwd=lwd,lty=lty)
    }
    invisible(x)    
}
