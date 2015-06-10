##' @export
folds<- function (n, folds = 10) 
{ ## {{{ 
   split(sample(1:n), rep(1:folds, length = n))
} ## }}} 

##' Split a data set and run function 
##'
##' @title Split a data set and run function 
##' @param func called function
##' @param data data-frame
##' @param size size of splits
##' @param splits number of splits (ignored if size is given)
##' @param id optional cluster variable
##' @param ... Additional arguments to lower level functions
##' @author Thomas Scheike, Klaus K. Holst
##' @export
##' @examples
##' library(timereg)
##' data(TRACE)
##' res <- divide.conquer(prop.odds,TRACE,
##' 	     formula=Event(time,status==9)~chf+vf+age,n.sim=0,size=200)
divide.conquer <- function(func=NULL,data,size,splits,id=NULL,...)
{ ## {{{
    nn <- nrow(data)
    if (!is.null(id)) {
        if (is.character(id) && length(id)==1) id <- data[,id]
        if (length(id)!=nn) stop("Wrong length of id variable")
        cc <- cluster.index(d$id)
        if (!missing(size)) splits <- round(cc$uniqueclust/size)
        splits <- min(splits,cc$uniqueclust)
        all.folds <- folds(cc$uniqueclust,splits)
        res <- lapply(all.folds, function(i) {
            idx <- na.omit(as.vector(cc$idclustmat[i,]+1))
            do.call(func, c(list(data=data[idx,]),...))
        })        
        return(res)
    }    
    splits <- round(nn/size)
    all.folds <- folds(nn,splits)
    res <- lapply(all.folds, function(i) 
	do.call(func, c(list(data=data[i,]),...)))
    res
} ## }}} 

##' Split a data set and run function of cox-aalen type and aggregate results 
##'
##' @title Split a data set and run function from timereg and aggregate
##' @param func called function
##' @param data data-frame
##' @param size size of splits
##' @param ... Additional arguments to lower level functions
##' @author Thomas Scheike, Klaus K. Holst
##' @export
##' @examples
##' library(timereg)
##' data(TRACE)
##' a <- divide.conquer.timereg(prop.odds,TRACE,
##'                             formula=Event(time,status==9)~chf+vf+age,n.sim=0,size=200)
##' coef(a)
##' a2 <- divide.conquer.timereg(prop.odds,TRACE,
##'                              formula=Event(time,status==9)~chf+vf+age,n.sim=0,size=500)
##' coef(a2)
##' 
##' if (interactive()) {
##' par(mfrow=c(1,1))
##' plot(a,xlim=c(0,8),ylim=c(0,0.01))
##' par(new=TRUE)
##' plot(a2,xlim=c(0,8),ylim=c(0,0.01))
##' }
divide.conquer.timereg <- function(func=NULL,data,size,...)
{ ## {{{ 

res <- divide.conquer(func=func,data,size,...)
K <- length(res)

gamma <- Reduce("+", lapply(res,function(x) x$gamma))/K
gammas <- do.call("cbind", lapply(res,function(x) x$gamma))
var.gamma <-    Reduce("+", lapply(res,function(x) x$var.gamma))/K^2
robvar.gamma <-    Reduce("+", lapply(res,function(x) x$robvar.gamma))/K^2 

times <- c(0,sort(unlist(lapply(res,function(x) x$cum[-1,1]))))
cum <- Reduce("+",lapply(res,function(x) Cpred(x$cum,times)))/K 
var.cum <-    Reduce("+", lapply(res,function(x)  Cpred(x$var.cum,times)))/K^2
robvar.cum <-    Reduce("+", lapply(res,function(x) Cpred(x$robvar.cum,times)))/K^2
robvar.cum[,1] <- var.cum[,1] <- times

res <- list(gamma=gamma, var.gamma=var.gamma, robvar.gamma=robvar.gamma,
	    gammas=gammas, sdgammas=apply(gammas,1,sd),
	    cum=cum,var.cum=var.cum,robvar.cum=robvar.cum,
            res=res,prop.odds=TRUE,score=0,D2linv=diag(nrow(gamma)))
class(res) <- "cox.aalen"

return(res)
} ## }}} 


