##' @export
biprobit.time <- function(formula,data,id,...,
                          breaks=NULL,n.times=20,pairs.only=TRUE,fix.cens.weights=FALSE,
                          cens.formula,cens.model="aalen",weights="w",messages=FALSE,
                          return.data=FALSE,
                          estimator="biprobit", summary.function) {

    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

    censtime <- model.extract(M, "response")
    status <- censtime[,2]
    time <- censtime[,1]
    outcome <- as.character(terms(formula)[[2]])
    jj <- jumptimes(time,data[,outcome],data[,id],sample=n.times)
    lastjump <- tail(jj,1)
    if (is.null(breaks)) {
        breaks <- jj
    }
    if (any(breaks>lastjump)) {
        breaks <- unique(pmin(breaks,tail(jj,1)))        
        if (messages) message("Looking at event before time ",max(breaks))
    }    

    outcome0 <- paste(outcome,"_dummy")
    res <- list(); k <- 0
    breaks <- rev(breaks)
    for (tau in breaks) {
        if (length(breaks)>1 && messages) message(tau)
        ## construct min(T_i,tau) or T_i and related censoring variable, 
        ## thus G_c(min(T_i,tau)) or G_c(T_i) as weights
        if ((fix.cens.weights==1 & k==0) | (fix.cens.weights==0)) {
            data0 <- data
            time0 <- time
            status0 <- status
        }
        cond0 <- time0>tau
        if (!fix.cens.weights) {
            status0[cond0 & status==1] <- 3 ## Not-censored if T>tau
        }
        data0[,outcome] <- data[outcome]
        data0[cond0,outcome] <- FALSE ## Non-case if T>tau
        if (!fix.cens.weights) time0[cond0] <- tau
        if ((fix.cens.weights & k==0) | (!fix.cens.weights)) {
            data0$S <- Surv(time0,status0==1)
            ## data0$status0 <- status0
            ## data0$time0 <- time0
            ## data0$y0 <- data0[,outcome]
            dataw <- ipw(update(cens.formula,S~.), data=data0, cens.model=cens.model,
                         weight.name=weights,obs.only=TRUE)            
        }
        if (fix.cens.weights) {
            timevar <- dataw$S[,1]
            dataw[,outcome] <- (dataw[,outcome])*(timevar<tau)
        }
        k <- k+1
        if (return.data) return(dataw)
        args <- c(list(x=formula,data=dataw,id=id,weights=weights, pairs.only=pairs.only), list(...))
        suppressWarnings(b <- do.call(estimator, args))
        ## suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weights=weights, pairs.only=pairs.only,...))
        res <- c(res,list(summary(b,...)))
    }
    if (length(breaks)==1) {        
        return(structure(b,time=breaks))
    }
    if (missing(summary.function)) {
        summary.function <- function(x,...) x$all
    } 
    mycoef <- lapply(rev(res),function(x) summary.function(x))
    res <- list(varname="Time",var=rev(breaks),coef=mycoef,summary=rev(res),call=m,type="time")
    class(res) <- "timemets"
    return(res)    
}


biprobit.time2 <- function(formula,data,id,...,
                          breaks=Inf,pairs.only=TRUE,
                          cens.formula,cens.model="aalen",weights="w") {

    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

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
        status0[cond0 & status==1] <- 3 ## Censored
        data0[cond0,outcome] <- FALSE
        time0[cond0] <- tau
        data0$S <- Surv(time0,status0==1)        
        dataw <- ipw(update(cens.formula,S~.), data=data0, cens.model=cens.model,
                     cluster=id,weight.name=weights,obs.only=TRUE)
        message("control")
        suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weights=weights, pairs.only=pairs.only,...))
        res <- c(res,list(summary(b)))
    }
    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,coef=lapply(res,function(x) x$all),summary=res,call=m,type="time")
    class(res) <- "timemets"
    return(res)    
}


