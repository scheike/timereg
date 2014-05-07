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
    res <- list(time=breaks,coef=lapply(res,function(x) x$all),summary=res,call=m)
    class(res) <- "timetwinlm"
    return(res)    
}

##' @S3method summary timetwinlm
summary.timetwinlm <- function(object,which=seq(nrow(object$coef[[1]])),...) {
    res <- list()
    for (i in which) {    
        rr <- matrix(unlist(lapply(object$coef,function(z) z[i,])),ncol=3,byrow=TRUE)
        colnames(rr) <- colnames(object$coef[[1]])
        rr <- cbind(time=object$time,rr)
        res <- c(res,list(rr))
    }
    names(res) <- rownames(object$coef[[1]])[which]
    return(res)
}


##' @S3method print timetwinlm
print.timetwinlm <- function(x,...) {
    res <- summary(x,...)
    for (i in seq_along(res)) {
        cat(i, ": ", names(res)[i], "\n",sep="")
        printCoefmat(res[[i]],...)
    }
    invisible(x)
}

##' @S3method plot timetwinlm
plot.timetwinlm <- function(x,...,which=1,
                            type="s",
                            lwd=2,lty=1,col,fillcol,alpha=0.2,
                            xlab=x$timevar,
                            ylab="",idx=seq_along(x$time),
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
        }
        add <- TRUE
        confband(zz[,1],zz[,3],zz[,4],polygon=TRUE,step=(type=="s"),col=fillcol[count],border=0)
        lines(zz[,1:2,drop=FALSE],lwd=lwd[count],lty=lty[count],col=col[count],type=type,...)
    }
    if (!is.null(legend) || !legend[1]) {
        if (is.logical(legend) || length(legend)==1) legend <- rownames(x$coef[[1]])[which]
        graphics::legend(legendpos,legend=legend,col=col,lwd=lwd,lty=lty)
    }
    invisible(x)    
}
