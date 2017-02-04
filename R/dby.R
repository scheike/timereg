##' Calculate summary statistics grouped by variable
##'
##' Calculate summary statistics grouped by 
##' @title Calculate summary statistics grouped by 
##' @param data Data.frame
##' @param INPUT Input variables (character or formula)
##' @param ... functions
##' @param ID id variable
##' @param ORDER (optional) order variable
##' @param SORT sort order (id+order variable)
##' @param COMBINE If TRUE result is appended to data
##' @param NOCHECK No sorting or check for missing data
##' @param ARGS Optional list of arguments to functions (...)
##' @param NAMES Optional vector of column names
##' @param FAST if FALSE fallback to slower (safer) function evaluation
##' @param COLUMN If TRUE do the calculations for each column
##' @param REDUCE Reduce number of redundant rows
##' @param REGEX Allow regular expressions
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
##' dby(d,y~id|num, dlag, ARGS=list(k=1:2))
##'
##' dby(d, y~id|num, dlag)
##' dby(d, y~id|num, y1=dlag, ARGS=list(k=1:2), NAMES=c("y1","y2"))
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
dby <- function(data,INPUT,...,ID=NULL,ORDER=NULL,SORT=0,COMBINE=!REDUCE,NOCHECK=FALSE,ARGS=NULL,NAMES,FAST=TRUE,COLUMN=FALSE,REDUCE=FALSE,REGEX=FALSE) {
    if (inherits(INPUT,"formula")) {
        INPUT <- procformdata(INPUT,sep="\\|",data=data,na.action=na.pass,do.filter=FALSE,regex=REGEX)
        if (is.null(ID)) ID <- INPUT$predictor
        if (is.null(ORDER)) {
            if (length(INPUT$group)==0) ORDER <- NULL
            else ORDER <- INPUT$group[[1]]
        }
        INPUT <- INPUT$response
    }
    funs <- list(...)
    if (length(funs)==0) stop("Please specify function")
    if (inherits(ID,"formula")) ID <- model.frame(ID,data=data,na.action=na.pass)
    if (inherits(ORDER,"formula")) ORDER <- model.frame(ORDER,data=data,na.action=na.pass)
    if (length(ID)==0) stop("id variable needed")
    id.orig <- NULL
    if (!COMBINE) id.orig <- ID
    if (NCOL(ID)>1) ID <- interaction(ID)
    if (is.data.frame(ID)) {
        ID <- as.integer(ID[,1,drop=TRUE])
    } else {
        ID <- as.integer(ID)
    }
    if (is.data.frame(ORDER)) {
        ORDER <- ORDER[,1,drop=TRUE]
    }
    if (!NOCHECK) {
        if (length(ORDER)>0) {
            ord <- order(ID,ORDER,decreasing=SORT,method="radix")
        } else {
            ord <- order(ID,decreasing=SORT,method="radix")
        }
        INPUT <- as.matrix(INPUT[ord,,drop=FALSE])
        ID <- ID[ord]
        if (COMBINE) data <- data[ord,,drop=FALSE]
        na.idx <- which(is.na(ID))
    } else {
        INPUT <- as.matrix(INPUT)
    }
    if (FAST) {
        resl <- lapply(funs, function(fun_) {
            env <- new.env()
            env$fun_ <- fun_
            if (length(ARGS)>0) {                
                ff <- function(...) do.call(fun_,c(list(...),ARGS))
                expr <- quote(ff(x_))
            } else {
                expr <- quote(fun_(x_))
            }
            .ApplyBy2(INPUT,ID,F=expr,Env=env,Argument="x_",Columnwise=COLUMN)
        })
    } else {
        resl <- lapply(funs, function(f) {
            f2 <- function(...) do.call(f,c(list(...),ARGS))
            .ApplyBy(INPUT,ID,f2)
        })
    }
    res <- Reduce(cbind,resl)
    if (missing(NAMES)) NAMES <- base::names(resl)
    try(colnames(res) <- NAMES,silent=TRUE)
    if (!NOCHECK && length(na.idx)>0) {
        res[na.idx,] <- NA
    }
    if (COMBINE) {
        res <- cbind(data,res)
    } else {
        res <- cbind(id.orig,res)
    }
    if (REDUCE!=0) {
        cl <- attr(resl[[1]],"clustersize")
        if (REDUCE<0) { # Grab the last observation
            idx <- cumsum(cl)
        } else { # Grab the first
            idx <- cumsum(c(1,cl[-length(cl)]))
        }       
        res <- res[unique(idx),,drop=FALSE]
        rownames(res) <- NULL
    }
    return(res)
}

