##' Calculate summary statistics grouped by variable
##'
##' Calculate summary statistics grouped by 
##' @title Calculate summary statistics grouped by 
##' @param data Data.frame
##' @param input.variable Input variables (character or formula)
##' @param ... functions
##' @param id.variable id
##' @param order.variable order
##' @param sort.order sort order (id+order variable)
##' @param combine.data If TRUE result is appended to data
##' @param no.check No sorting or check for missing data
##' @param args Optional list of arguments to functions (...)
##' @param names Optional vector of column names
##' @param fast if FALSE fallback to slower (safer) function evaluation
##' @export 
##' @author Klaus K. Holst and Thomas Scheike
##' @examples
##' n <- 4
##' k <- c(3,rbinom(n-1,3,0.5)+1)
##' N <- sum(k)
##' d <- data.frame(y=rnorm(N),x=rnorm(N),id=rep(seq(n),k),num=unlist(sapply(k,seq)))
##' d2 <- d[sample(nrow(d)),]
##'
##' dby(d, y~id, mean)
##' dby(d, y~id|num, cumsum)
##'
##' dby(d,y~id|num, dlag, args=list(k=1:2))
##'
##' dby(d, y~id|num, dlag)
##' dby(d, y~id|num, y1=dlag, args=list(k=1:2), names=c("y1","y2"))
##' dby(d, y~id|num, mean=mean, csum=cumsum, n=length)
##' dby(d2,y~id|num, a=cumsum, b=mean, N=length, l1=function(x) c(NA,x)[-length(x)])
##' 
##' dby(d, y~id|num, nn=seq, n=length)
##' 
##' f <- function(x) apply(x,1,min)
##' dby(d, y+x~id, min=f)
##' 
##' dby(d,y+x~id|num, function(x) x)
##' 
##' f <- function(x) { cbind(cumsum(x[,1]),cumsum(x[,2]))/sum(x)}
##' dby(d, y+x~id, f)
dby <- function(data,input.variable,...,id.variable=NULL,order.variable=NULL,sort.order=0,combine.data=TRUE,no.check=FALSE,args=NULL,names,fast=TRUE) {
    if (inherits(input.variable,"formula")) {
        input.variable <- procformdata(input.variable,sep="\\|",data=data,na.action=na.pass,do.filter=FALSE)
        if (is.null(id.variable)) id.variable <- input.variable$predictor
        if (is.null(order.variable)) {
            if (length(input.variable$group)==0) order.variable <- NULL
            else order.variable <- input.variable$group[[1]]
        }
        input.variable <- input.variable$response
    }
    funs <- list(...)
    if (length(funs)==0) stop("Please specify function")
    if (inherits(id.variable,"formula")) id.variable <- model.frame(id.variable,data=data,na.action=na.pass)
    if (inherits(order.variable,"formula")) order.variable <- model.frame(order.variable,data=data,na.action=na.pass)
    if (length(id.variable)==0) stop("id variable needed")
    if (NCOL(id.variable)>1) id.variable <- interaction(id.variable)
    if (is.data.frame(id.variable)) {
        id.variable <- as.integer(id.variable[,1,drop=TRUE])
    } else {
        id.variable <- as.integer(id.variable)
    }
    if (is.data.frame(order.variable)) {
        order.variable <- order.variable[,1,drop=TRUE]
    }
    if (!no.check) {
        if (length(order.variable)>0) {
            ord <- order(id.variable,order.variable,decreasing=sort.order,method="radix")
        } else {
            ord <- order(id.variable,decreasing=sort.order,method="radix")
        }
        input.variable <- as.matrix(input.variable[ord,,drop=FALSE])
        id.variable <- id.variable[ord]
        if (combine.data) data <- data[ord,,drop=FALSE]
        na.idx <- which(is.na(id.variable))
    } else {
        input.variable <- as.matrix(input.variable)
    }
    if (fast) {
        resl <- lapply(funs, function(fun_) {
            env <- environment()
            env$fun_ <- fun_
            if (length(args)>0) {                
                ff <- function(...) do.call(fun_,c(list(...),args))
                expr <- quote(ff(x_))
            } else {
                expr <- quote(fun_(x_))
            }
            .Call("mets_ApplyBy2",
                  idata=input.variable,
                  icluster=id.variable,
                  F=expr,
                  Env=env)
        })
    } else {
        resl <- lapply(funs, function(f) {
            f2 <- function(...) do.call(f,c(list(...),args))
            ApplyBy(input.variable,id.variable,f2)
        })
    }
    res <- Reduce(cbind,resl)
    if (missing(names)) names <- base::names(resl)
    try(colnames(res) <- names,silent=TRUE)
    if (!no.check && length(na.idx)>0) {
        res[na.idx,] <- NA
    }
    if (combine.data) {
        return(cbind(data,res))
    }
    return(res)
}

