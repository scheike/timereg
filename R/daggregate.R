by2mat <- function(x,nam,...) {
    nulls <- which(unlist(lapply(x,is.null)))
    nonnulls <- setdiff(seq_along(x),nulls)
    nn <- do.call("expand.grid",attributes(x)$dimnames)
    if (length(nulls)>0) nn <- nn[-nulls,,drop=FALSE]
    res <- Reduce("rbind",x)
    if (is.null(colnames(res)) && !missing(nam)) {
        colnames(res) <- nam[seq(length(ncol(res)))]
    }
    suppressWarnings(res <- cbind(nn,res)) ## no warnings on row-names
    for (i in seq(ncol(res)-1)+1) {
        if (is.list(res[,i])) {
            if (!is.null(nn <- names(res[,i][[1]])))
                colnames(res)[i] <- paste0(colnames(res)[i],"(",paste0(nn,collapse=","),")")
        }
    }
    a <- rownames(x[[1]])
    res$"_var" <- a
    rownames(res) <- seq(nrow(res))
    return(res)
}


##' aggregating for for data frames
##'
##' aggregating for for data frames
##' @examples
##' data("sTRACE",package="timereg")
##' daggregate(iris, "^.e.al", x="Species", fun=cor, regex=TRUE)
##' daggregate(iris, Sepal.Length+Petal.Length ~Species, fun=summary)
##' daggregate(iris, log(Sepal.Length)+I(Petal.Length>1.5) ~ Species,
##'                  fun=summary)
##' daggregate(iris, "*Length*", x="Species", fun=head)
##' daggregate(iris, "^.e.al", x="Species", fun=tail, regex=TRUE)
##' daggregate(sTRACE, status~ diabetes, fun=table)
##' daggregate(sTRACE, status~ diabetes+sex, fun=table)
##' daggregate(sTRACE, status + diabetes+sex ~ vf+I(wmi>1.4), fun=table)
##' daggregate(iris, "^.e.al", x="Species",regex=TRUE)
##' dlist(iris,Petal.Length+Sepal.Length ~ Species |Petal.Length>1.3 & Sepal.Length>5,
##'             n=list(1:3,-(3:1)))
##' daggregate(iris, I(Sepal.Length>7)~Species | I(Petal.Length>1.5))
##' daggregate(iris, I(Sepal.Length>7)~Species | I(Petal.Length>1.5),
##'                  fun=table)
##'
##' dsum(iris, .~Species, matrix=TRUE, missing=TRUE)
##'
##' par(mfrow=c(1,2))
##' data(iris)
##' drename(iris) <- ~.
##' daggregate(iris,'sepal*'~species|species!="virginica",fun=plot)
##' daggregate(iris,'sepal*'~I(as.numeric(species))|I(as.numeric(species))!=1,fun=summary)
##'
##' dnumeric(iris) <- ~species
##' daggregate(iris,'sepal*'~species.n|species.n!=1,fun=summary)
##'
##' @export
##' @param data data.frame
##' @param y name of variable, or formula, or names of variables on data frame.
##' @param x name of variable, or formula, or names of variables on data frame.
##' @param subset subset expression
##' @param ... additional arguments to lower level functions
##' @param fun function defining aggregation
##' @param regex interpret x,y as regular expressions
##' @param missing Missing used in groups (x)
##' @param remove.empty remove empty groups from output
##' @param matrix if TRUE a matrix is returned instead of an array
##' @param silent suppress messages
##' @param na.action How model.frame deals with 'NA's
##' @aliases daggr
daggregate <- function(data,y=NULL,x=NULL,subset,...,fun="summary",regex=mets.options()$regex, missing=FALSE, remove.empty=FALSE, matrix=FALSE, silent=FALSE, na.action=na.pass, convert=NULL)
{# {{{
    if (is.vector(data)) data <- data.frame(data)
    subs <- substitute(subset)
    if (!base::missing(subs)) {
        expr <- suppressWarnings(inherits(try(subset,silent=TRUE),"try-error"))
        if (expr) data <- data[which(eval(subs,envir=data)),,drop=FALSE]
        else data[subset,,drop=FALSE]
    }
    if (is.null(y)) y <- colnames(data)
    if (inherits(y,"formula")) {
        yx <- procformdata(y,sep="\\|",data=data,na.action=na.action,regex=regex,...)
        y <- yx$response
        x0 <- yx$predictor
        if (is.null(x) && length(y)>0) x <- x0
        if (NCOL(x)==0) x <- NULL
        if (length(y)==0) {
            y <- x0
        }
    } else {
        yy <- c()
        for (y0 in y) {
            if (!regex) y0 <- glob2rx(y0)
            n <- grep(y0,names(data),perl=mets.options()$regex.perl)
            yy <- union(yy,names(data)[n])
        }
        y <- data[,yy,drop=FALSE]
    }
    if (is.character(x) && length(x)<NROW(data)) {
        xx <- c()
        for (x0 in x) {
            if (!regex) x0 <- glob2rx(x0)
            n <- grep(x0,names(data),perl=mets.options()$regex.perl)
            xx <- union(xx,names(data)[n])
        }
        x <- data[,xx,drop=FALSE]
    }
    if (inherits(x,"formula")) {
        x <- model.frame(x,data=data,na.action=na.action)
    }
    if (!is.null(x)) {
        xidx <- na.omit(match(colnames(x),colnames(y)))
        if (length(xidx)>0) y <- y[,-xidx,drop=FALSE]
    }
    if (is.character(fun)) fun <- get(fun)
    if (!is.null(convert) && is.logical(convert)) {
        if (convert) convert <- as.matrix
        else convert <- NULL
    }
    if (!is.null(convert)) {
        fun_ <- fun
        fun <- function(x,...) fun_(convert(x,...))
    }
    if (!is.null(x)) {
        if (missing) {
            x[is.na(x)] <- 'NA'
        }
        if (silent) {
                capture.output(res <- by(y,x,fun,...))
        } else {
            res <- by(y,x,fun,...)
        }
        if (remove.empty) {
            # ... need to abandon 'by'?
        }
        if (matrix) {
            res <- by2mat(res,colnames(y))
        }
        return(structure(res,ngroupvar=NCOL(x),class=c("daggregate",class(res))))
    }
    if (silent)
        capture.output(res <- do.call(fun, c(list(y),list(...))))
    else
        res <- do.call(fun, c(list(y),list(...)))
    res
    structure(res, ngroupvar=0, class=c("daggregate",class(res)))
}# }}}

##' @export
daggr <- function(data,...,convert=as.matrix) daggregate(data,...,convert=convert)

##' @export
print.daggregate <- function(x,quote=FALSE,...) {
    attr(x,c("ngroupvar")) <- NULL
    class(x) <- class(x)[-1]
    print(x,quote=quote,...)
}


##' @export
dhead <- function(data,y=NULL,x=NULL,...) daggregate(data,y,x,fun=function(z) utils::head(z,...))

##' @export
dtail <- function(data,y=NULL,x=NULL,...) daggregate(data,y,x,fun=function(z) utils::tail(z,...))

##' @export
dsummary <- function(data=NULL,y=NULL,x=NULL,...) daggregate(data,y,x,fun=function(z) base::summary(z,...))

##' @export
dstr <- function(data,y=NULL,x=NULL,...) invisible(daggregate(data,y,x,fun=function(z) utils::str(z,...)))

##' @export
dunique <- function(data,y=NULL,x=NULL,...) invisible(daggregate(data,y,x,fun=function(z) base::unique(z,...)))


##' summary, tables, and correlations for data frames
##'
##' summary, tables, and correlations for data frames
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param y name of variable, or fomula, or names of variables on data frame.
##' @param x possible group variable
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike
##' @examples
##' data("sTRACE",package="timereg")
##' dt<- sTRACE
##' dt$time2 <- dt$time^2
##' dt$wmi2 <- dt$wmi^2
##' head(dt)
##'
##' dcor(dt)
##'
##' dcor(dt,~time+wmi)
##' dcor(dt,~time+wmi,~vf+chf)
##' dcor(dt,time+wmi~vf+chf)
##'
##' dcor(dt,c("time*","wmi*"),~vf+chf)
##' @aliases dsummary dstr dcor dsubset dquantile dcount dmean dscalar deval deval2 dsum dsd
##' @export
dcor <- function(data,y=NULL,x=NULL,...) daggregate(data,y,x,...,fun=function(z,...) stats::cor(z,...))

##' @export
dscalar <- function(data,y=NULL,x=NULL,...,na.rm=TRUE,matrix=TRUE,fun=base::mean) {
    daggregate(data,y,x,matrix=matrix,...,
               fun=function(z,...) {
                   if (is.matrix(z)) {
                       apply(z,2,function(x) 
                           suppressWarnings(tryCatch(fun(x,na.rm=na.rm,...),error=function(e) return(NA))))
                   } else {                           
                       unlist(lapply(z,function(x) {
                           suppressWarnings(tryCatch(fun(x,na.rm=na.rm,...),error=function(e) return(NA)))
                       }))
                   }
               })
}


Summary <- function(object,na.rm=TRUE,...) {
    if (is.numeric(object)) {
        x <- c(summary(object,...),sd=sd(object,na.rm=TRUE))
    } else {
        x <- summary(object,...)
    }
    ## Formatting
    xx <- x
    if (is.numeric(x) || is.complex(x)) {
        finite <- is.finite(x)
        xx[finite] <- zapsmall(x[finite])
      }
    m <- match("NA's", names(xx), 0)
    if (inherits(x, "Date") || inherits(x, "POSIXct")) {
        xx <- if (length(a <- attr(x, "NAs"))) 
                 c(format(xx), `NA's` = as.character(a))
             else format(xx)
    }
    else if (m && !is.character(x)) 
        xx <- c(format(xx[-m]), `NA's` = as.character(xx[m]))
    xx
}

##' @export
deval2 <- function(data,...,matrix=simplify,simplify=TRUE)  deval(data,matrix=TRUE,simplify=TRUE,...)

##' @export
deval <- function(data,y=NULL,x=NULL,...,matrix=FALSE,fun=Summary,simplify=FALSE) {
    if (is.list(fun)) {
        newf <- function(x,...) {
            unlist(lapply(fun,function(f) f(x,...), ...))
        }
    }
    else newf <- fun
    res <- daggregate(data,y,x,matrix=matrix,...,
                     fun=function(z) lapply(z,function(x) {
                         suppressWarnings(tryCatch(newf(x,...),error=function(e) return(NA)))
                     }))
    if (simplify) {
        for (i in seq_len(ncol(res))) {
            if (is.list(res[,i])) res[,i] <- unlist(res[,i])
        }
    ##     Dim <- function(x) {
    ##         val <- dim(x)
    ##         if (is.null(val)) val <- c(1,length(x))
    ##         val
    ##     }
    ##     dm <- Dim(res[[1]])
    ##     dims <- unlist(lapply(res,function(x) identical(Dim(x),dm)))
    ##     if (all(dims)) {
    ##         Res <- res
    ##         n <- length(res)
    ##         res <- array(NA,dim=c(n,dm))
    ##         for (i in seq(n)) {
    ##             browser()
    ##         }
                
    ##     }
        ## }
    }
    res
}


##' @export
dmean <- function(data,...) dscalar(data,fun=base::mean,...)

##' @export
dsum <- function(data,...) dscalar(data,fun=base::sum,...)

##' @export
dsd <- function(data,...) dscalar(data,fun=stats::sd,...)


##' @export
dcount <- function(data,x=NULL,...,na.rm=TRUE) {
    res <- rbind(daggregate(data,x=x,matrix=TRUE,...,fun=function(z,...) NROW(z)))
    res[is.na(res)] <- 0
    rownames(res) <- seq(nrow(res))
    colnames(res)[ncol(res)] <- "count"
    res
}


##' @export
dsubset <- function(data,...) {
    daggregate(data,...,fun=function(z) z)
}


##' @export
dquantile <- function(data,y=NULL,x=NULL,probs=seq(0,1,by=1/breaks),breaks=4,matrix=TRUE,reshape=FALSE,...,na.rm=TRUE) {
    a <- daggregate(data,y,x,matrix=FALSE,...,fun=function(z,...) apply(z,2,function(x,...) quantile(x,probs=probs,na.rm=na.rm,...)))
    if (matrix) {
        res <- by2mat(a)
        xidx <- seq_len(attr(a, "ngroupvar"))
        if (!reshape || is.null(res[,"_var"]) || length(xidx)==0) return(res)
        res <- dreshape(res, id=colnames(res)[xidx], num="_var",sep="_")
        return(res)
    }
    return(a)
}

