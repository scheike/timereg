##' Calculate summary statistics grouped by variable
##'
##' Calculate summary statistics grouped by
##' @title Calculate summary statistics grouped by
##' @param data Data.frame
##' @param INPUT Input variables (character or formula)
##' @param ... functions
##' @param ID id variable
##' @param ORDER (optional) order variable
##' @param SUBSET (optional) subset expression
##' @param SORT sort order (id+order variable)
##' @param COMBINE If TRUE result is appended to data
##' @param NOCHECK No sorting or check for missing data
##' @param ARGS Optional list of arguments to functions (...)
##' @param NAMES Optional vector of column names
##' @param COLUMN If TRUE do the calculations for each column
##' @param REDUCE Reduce number of redundant rows
##' @param REGEX Allow regular expressions
##' @param ALL if FALSE only the subset will be returned
##' @export
##' @author Klaus K. Holst and Thomas Scheike
##' @details
##' dby2 for column-wise calculations
##' @aliases dby dby<- dby2 dby2<- dbyr
##' @examples
##' n <- 4
##' k <- c(3,rbinom(n-1,3,0.5)+1)
##' N <- sum(k)
##' d <- data.frame(y=rnorm(N),x=rnorm(N),id=rep(seq(n),k),num=unlist(sapply(k,seq)))
##' d2 <- d[sample(nrow(d)),]
##'
##' dby(d2, y~id, mean)
##' dby(d2, y~id + order(num), cumsum)
##'
##' dby(d,y ~ id + order(num), dlag)
##' dby(d,y ~ id + order(num), dlag, ARGS=list(k=1:2))
##' dby(d,y ~ id + order(num), dlag, ARGS=list(k=1:2), NAMES=c("l1","l2"))
##'
##' dby(d, y~id + order(num), mean=mean, csum=cumsum, n=length)
##' dby(d2, y~id + order(num), a=cumsum, b=mean, N=length, l1=function(x) c(NA,x)[-length(x)])
##'
##' dby(d, y~id + order(num), nn=seq_along, n=length)
##' dby(d, y~id + order(num), nn=seq_along, n=length)
##'
##' d <- d[,1:4]
##' dby(d, x<0) <- list(z=mean)
##' d <- dby(d, is.na(z), z=1)
##'
##' f <- function(x) apply(x,1,min)
##' dby(d, y+x~id, min=f)
##'
##' dby(d,y+x~id+order(num), function(x) x)
##'
##' f <- function(x) { cbind(cumsum(x[,1]),cumsum(x[,2]))/sum(x)}
##' dby(d, y+x~id, f)
##'
##' ## column-wise
##' a <- d
##' dby2(a, mean, median, REGEX=TRUE) <- '^[y|x]'~id
##' a
##' ## wildcards 
##' dby2(a,'y*'+'x*'~id,mean) 
##'
##'
##' ## subset
##' dby(d, x<0) <- list(z=NA)
##' d
##' dby(d, y~id|x>-1, v=mean,z=1)
##' dby(d, y+x~id|x>-1, mean, median, COLUMN=TRUE)
##'
##' dby2(d, y+x~id|x>0, mean, REDUCE=TRUE)
##'
##' dby(d,y~id|x<0,mean,ALL=FALSE)
##'
##' a <- iris
##' a <- dby(a,y=1)
##' dby(a,Species=="versicolor") <- list(y=2)
dby <- function(data,INPUT,...,ID=NULL,ORDER=NULL,SUBSET=NULL,SORT=0,COMBINE=!REDUCE,NOCHECK=FALSE,ARGS=NULL,NAMES,COLUMN=FALSE,REDUCE=FALSE,REGEX=mets.options()$regex,ALL=TRUE) {
    if (missing(INPUT)) INPUT <- .~1
    val <- substitute(INPUT)
    INPUT <- try(eval(val),silent=TRUE)
    funs <- substitute(list(...))[-1]
    if (inherits(INPUT,"try-error")) {
        INPUT <- as.formula(paste0(".~1|",deparse(val)))        
    }
    if (inherits(INPUT,c("formula","character"))) {
        INPUT <- procformdata(INPUT,sep="\\|",data=data,na.action=na.pass,do.filter=FALSE,regex=REGEX,specials="order")
        if (is.null(ID)) ID <- INPUT$predictor
        if (is.null(SUBSET)) {
            if (length(INPUT$group)!=0)
                SUBSET <- INPUT$group[[1]]
        }
        if (is.null(ORDER)) {
            if (length(INPUT$specials$order)!=0)
                ORDER <- INPUT$specials$order[[1]]
        }
        INPUT <- INPUT$response
    }

    if (inherits(ID,"formula")) ID <- model.frame(ID,data=data,na.action=na.pass)
    if (inherits(ORDER,"formula")) ORDER <- model.frame(ORDER,data=data,na.action=na.pass)
    if (inherits(SUBSET,"formula")) SUBSET <- model.frame(SUBSET,data=data,na.action=na.pass)
    if (is.null(INPUT)) {
        INPUT <- ID
        ID <- NULL
    }
    noID <- FALSE
    if (length(ID)==0) {
        noID <- TRUE
        ID = rep(1,NROW(INPUT))
    }
    group <- NULL
    if (!COMBINE || length(funs)==0) group <- ID
    if (NCOL(ID)>1) ID <- interaction(ID)
    if (is.data.frame(ID)) {
        ID <- as.numeric(ID[,1,drop=TRUE])
    } else {
        ID <- as.numeric(ID)
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
        numerics <- unlist(lapply(INPUT[1,],is.numeric))
        INPUT <- as.matrix(INPUT[ord,which(numerics),drop=FALSE])
        ID <- ID[ord]
        if (is.null(group)) data <- data[ord,,drop=FALSE]
        else {
            if (is.vector(group)) group <- group[ord]
            else group <- group[ord,,drop=FALSE]
        }
        na.idx <- which(is.na(ID))
    } else {
        INPUT <- as.matrix(INPUT)
    }

    if (length(SUBSET)>0) {
        SUBSET <- which(as.matrix(SUBSET))
        INPUT <- INPUT[SUBSET,,drop=FALSE]
        ID <- ID[SUBSET]
    }
    if (length(funs)==0) {
        if (noID) return(INPUT)
        return(cbind(INPUT,group))
    }

    if (NROW(INPUT)>0) {
    resl <- lapply(funs, function(fun_) {
        env <- new.env()
        isfun <- tryCatch(is.function(eval(fun_)),error=function(...) FALSE)
        if (isfun) {
            env$fun_ <- eval(fun_)
            if (length(ARGS)>0) {
                fun_ <- eval(fun_)
                ff <- function(...) do.call(fun_,c(list(...),ARGS))
                expr <- quote(ff(x))
            } else {
                expr <- quote(fun_(x))
            }
        } else {
            expr <- fun_
        }
        .ApplyBy2(INPUT,ID,F=expr,Env=env,Argument="x",Columnwise=COLUMN)
    })
    res <- Reduce(cbind,resl)
    } else {
        resl <- NULL
        res <- NULL
    }

    dn <- NULL
    ## Setting column names
    if (missing(NAMES)) {
        a <- match.call(expand.dots=FALSE)
        ff <- lapply(a[["..."]],deparse)
        fn <- lapply(ff,deparse)
        dn <- names(ff)
        if (is.null(dn)) dn <- rep("",length(ff))
        idx <- which(dn=="")
        if (length(idx)>0) {
            newn <- unlist(ff[idx])
            dn[idx] <- newn[seq_along(idx)]
            fcall <- grepl("^function",dn)
            if (any(fcall))
                dn[which(fcall)] <- which(fcall)
        }
        numcolf <- unlist(lapply(resl, NCOL))
        nn <- c()
        if (COLUMN) {
            colnames(INPUT)
            dn <- unlist(lapply(dn, function(x) paste(x,colnames(INPUT),sep=".")))
        }
        for (i in seq_along(resl)) {
            nc <- numcolf[i]
            curnam <- dn[i]
            if (COLUMN) {
                nc <- nc/NCOL(INPUT)
                pos <- NCOL(INPUT)*(i-1)+1
                pos <- seq(pos,pos+NCOL(INPUT)-1)
                curnam <- dn[pos]
            }
            if (nc>1) {
                nn <- c(nn, unlist(lapply(curnam, function(x) paste0(x,seq(nc)))))
            }
            else nn <- c(nn,curnam)
        }
        NAMES <- nn
    }

    try(colnames(res) <- NAMES,silent=TRUE)
    numidx <- grep("^[0-9]",colnames(res)) ## column names starting with digit
    if (length(numidx)>0)
       colnames(res)[numidx] <- paste0("_",colnames(res)[numidx])

    if (!NOCHECK && length(na.idx)>0) {
        res[na.idx,] <- NA
    }
    cl <- attr(resl[[1L]],"clustersize")
    if (COMBINE) {
        if (is.null(NAMES) & is.null(res) & length(dn)>0) {
            res <- matrix(NA,nrow(data),length(dn))
            colnames(res) <- dn
        }        
        nn <- colnames(res)
        nc <- ifelse(is.null(res), 0, ncol(res))
        res0 <- matrix(NA,nrow(data),nc)
        colnames(res0) <- nn
        didx <- which(colnames(data)%in%nn)
        if (length(didx)>0) {
            ridx <- match(colnames(data)[didx],nn)
            res0[,ridx] <- as.matrix(data[,didx])
            data <- data[,-didx,drop=FALSE]
        }
        if (length(SUBSET)>0) {
            res0[SUBSET,] <- res
            res <- cbind(data,res0)
        } else {
            if (!is.null(res))
                res <- cbind(data,res)
            else
                res <- data
        }

    } else {
        if (!noID) {
            if (length(SUBSET)>0) group <- group[SUBSET,,drop=FALSE]
            res <- cbind(group,res)
        }
    }
    if (REDUCE!=0) {
        if (REDUCE<0) { # Grab the last observation
            idx <- cumsum(cl)
        } else { # Grab the first
            idx <- cumsum(c(1,cl[-length(cl)]))
        }
        res <- res[unique(idx),,drop=FALSE]
        if (NROW(res)==1) {
            rownames(res) <- ""
        } else {
            rownames(res) <- NULL
        }
    }
    if (!ALL) {
        return(res[SUBSET,,drop=FALSE])
    }
    return(res)
}


##' @export
"dby<-" <- function(data,INPUT,...,value) {
    cl <- match.call()
    cl[[1L]] <- substitute(dby)
    a <- substitute(value)
    if (inherits(value,"formula")) {
        cl["value"] <- NULL
        names(cl)[names(cl)=="INPUT"] <- ""
        cl[["INPUT"]] <- value
    } else {
        if (is.list(value)) {
            cl[which(names(cl)=="value")] <- NULL
            start <- length(cl)
            for (i in seq_along(value)) {
                cl[start+i] <- value[i]
            }
            if (length(names(value))>0)
                names(cl)[start+seq_along(value)] <- names(value)
        } else {
            names(cl)[which(names(cl)=="value")] <- ""
        }
    }
    eval.parent(cl)
}

##' @export
dby2 <- function(data,INPUT,...) {
    cl <- match.call()
    cl[[1L]] <- substitute(dby)
    cl[["COLUMN"]] <- TRUE
    eval.parent(cl)
}

##' @export
"dby2<-" <- function(data,INPUT,...,value) {
    cl <- match.call()
    cl[[1L]] <- substitute(`dby<-`)
    cl[["COLUMN"]] <- TRUE
    eval.parent(cl)
}

##' @export
dbyr <- function(data,INPUT,...,COLUMN=FALSE) {
    cl <- match.call()
    cl[[1L]] <- substitute(`dby`)
    cl[["REDUCE"]] <- TRUE
    cl[["COLUMN"]] <- COLUMN
    eval.parent(cl)
}

