##' @export
biprobit.time <- function(formula,data,id,...,
                          breaks=NULL,n.times=20,pairs.only=TRUE,fix.censweights=FALSE,
                          cens.formula,cens.model="aalen",weights="w",messages=FALSE,
                          return.data=FALSE) {

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
    if (is.null(breaks)) {
        breaks <- jumptimes(time,data[,outcome],data[,id],sample=n.times)
    }

    outcome0 <- paste(outcome,"_dummy")
    res <- list(); k <- 0
    breaks <- rev(breaks)
    for (tau in breaks) {
        if (length(breaks)>1 && messages) message(tau)
        ## construct min(T_i,tau) or T_i and related censoring variable, 
        ## thus G_c(min(T_i,tau)) or G_c(T_i) as weights
        if ((fix.censweights==1 & k==0) | (fix.censweights==0)) {
            data0 <- data
            time0 <- time
            status0 <- status
        }
        cond0 <- time0>tau
        if (!fix.censweights) status0[cond0 & status==1] <- 3 ## Censored
        data0[,outcome] <- data[outcome]
        data0[cond0,outcome] <- FALSE
        if (!fix.censweights) time0[cond0] <- tau
        if ((fix.censweights & k==0) | (!fix.censweights)) {
            data0$S <- Surv(time0,status0==1)
            dataw <- ipw(update(cens.formula,S~.), data=data0, cens.model=cens.model,
                         cluster=id,weightname=weights,obsonly=TRUE)
        }
        if (fix.censweights) {
            timevar <- dataw$S[,1]
            dataw[,outcome] <- (dataw[,outcome])*(timevar<tau)
        }
        k <- k+1
        if (return.data) return(dataw)
        suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weights=weights, pairs.only=pairs.only,...))
        res <- c(res,list(summary(b,...)))
    }
    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=rev(breaks),coef=lapply(rev(res),function(x) x$all),summary=rev(res),call=m,type="time")
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
                     cluster=id,weightname=weights,obsonly=TRUE)
        message("control")
        suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weights=weights, pairs.only=pairs.only,...))
        res <- c(res,list(summary(b)))
    }
    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,coef=lapply(res,function(x) x$all),summary=res,call=m,type="time")
    class(res) <- "timemets"
    return(res)    
}


