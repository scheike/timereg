gsub2 <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
      pattern <- rep(pattern, length.out=length(replacement))
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- sub(pattern[i], replacement[i], result, ...)
  }
  result
}

procform <- function(formula, sep="\\|", nsep=1, return.formula=FALSE, data=NULL,
             no.match=TRUE, regex=FALSE, return.list=TRUE, ...) {
    res <- NULL
    if (is.null(formula)) {
        res <- colnames(data)
    } else if (is.character(formula)) {
        if (is.null(data)) {
            res <- unique(formula)
        } else {
            yy <-c()
            for (y0 in formula) {
                y0orig <- y0
                if (!regex) y0 <- glob2rx(y0)
                npos <- grep(y0,names(data),perl=mets.options()$regex.perl)
                if (no.match && length(npos)==0) {
                    yy <- union(yy, y0orig)
                } else {
                    yy <- union(yy,names(data)[npos])
                }
            }
            res <- unique(yy)
        }
    }
    if (is.numeric(formula)) res <- colnames(data)[formula]
    if (is.character(res)) {
        if (!return.list) return(res)
        if (return.formula) return(as.formula(paste("~",paste(res,collapse="+"))))
        return(list(response=res,predictor=NULL,filter=NULL))
    }

    
    ## Add parantheses around quotes if it is not a function call
    if (inherits(formula,"formula")) {

        st <- Reduce(paste,deparse(formula))
        strsplit(st,"\"")
        quotepos <- gregexpr("[\"']",st)[[1]]
        if (quotepos[1]>0) {
            sts <- strsplit(st,"[\"']")[[1]]
            foundsep <- any(grepl("|", sts, fixed=TRUE))
            p <- length(quotepos)
            ##repl <- rep(c("(\"","\")"),p)
            for (i in seq(p/2)*2-1) {
                sts[i] <- paste0(sts[i],"(\"")
                sts[i+1] <- paste0(sts[i+1],"\")")
            }
            ## To handle regular expression entered as strings in the formula, we add a 'filter' expression at the end of the formula
            if (!foundsep) sts <- c(sts,"|1")
            formula <- as.formula(paste(sts,collapse=""))
        }
    }

    aa <- attributes(terms(formula,data=data,specials="regex"))
    if (aa$response == 0) {
        res <- NULL
    } else {
        res <- paste(deparse(formula[[2]]), collapse = "")
    }
    filter.expression <- NULL
    foundsep <- FALSE
    pred <- filter <- c()
    if (!is.null(sep) && length(aa$term.labels) > 0) {
        foundsep <- any(grepl(sep,aa$term.labels))
        if (foundsep) {
            if (nsep>1) {
                xc <- gsub(" ","",unlist(lapply(aa$term.labels, function(z) strsplit(z,sep)[[1]])))
                pred <- xc[1]
                filter <- xc[-1]
            } else {
                xc <- gsub(" ","",unlist(lapply(aa$term.labels, function(z) {
                    spl <- regexpr(sep,z) ## first appearance
                    pred <- substr(z,1,spl-1)
                    filter <- substr(z,spl+1,nchar(z))
                    return(c(pred,filter))
                })))
                pred <- xc[1]
                filter <- xc[2]
            }
            if (any(pred==".")) {
                f <- as.formula(paste0(paste0(c(res,filter),collapse="+"),"~."))
                x <- attributes(terms(f,data=data))$term.labels
                pred <- x
            }
            if (filter%in%c("1","0","-1")) { 
                filter <- list()
                filter.expression <- NULL
            } else {
                filter.expression <- parse(text=filter)
                filter <- as.list(filter)
            }
        }
    }
    if (!foundsep) pred <- aa$term.labels

    expandst <- function(st) {
        st <- res <- unlist(strsplit(gsub(" ","",st),"\\+"))
        if (any(unlist(lapply(st, function(x) grepl("^\\(",x))))) {
            res <- c()
            for (x in st) {
                if (grepl("^\\(",x)) {
                    x <- gsub('\\"',"",x)
                    x <- gsub('^\\(',"",x)
                    x <- gsub('\\)$',"",x)
                    res <- c(res,unlist(procform(x,data=data,regex=regex, no.match=FALSE)$response))
                } else {
                    res <- c(res,x)
                }
                res <- unique(res)
            }
        }        
        return(res)
    }
    res <- expandst(res)
    pred <- expandst(pred)
    filter <- lapply(filter, expandst)

    if (return.formula) {
        if (foundsep && !is.null(filter)) {
            filter <- lapply(filter, function(z) as.formula(paste0(c("~", paste0(z,collapse="+")))))
        }
        if (length(pred)>0)
            pred <- as.formula(paste0(c("~", paste0(pred,collapse="+"))))
        if (length(res)>0)
            res <- as.formula(paste0(c("~", paste0(res,collapse="+"))))
    }
    res <- list(response=res, predictor=pred, filter=filter, filter.expression=filter.expression)
    if (!return.list) return(unlist(unique(res)))
    return(res)
}

procformdata <- function(formula,data,sep="\\|", na.action=na.pass, do.filter=TRUE, ...) {
    res <- procform(formula,sep=sep,data=data,return.formula=TRUE,...)
    if (inherits(res,"formula")) {
        res <- list(response=res)
    }
    y <- x <- NULL
    filter <- res$filter.expression
    if (!do.filter) {
        filter <- NULL
    }
    if (length(res$response)>0) {
        if (is.null(filter)) y <- model.frame(res$response,data=data,na.action=na.action)
        else y <- model.frame(res$response,data=subset(data,eval(filter)),na.action=na.action)
    }
    if (length(res$predictor)>0) {
        if (is.null(filter)) x <- model.frame(res$predictor,data=data,na.action=na.action)
        else x <- model.frame(res$predictor,data=subset(data,eval(filter)),na.action=na.action)

    }
    if (!do.filter) {
        group <- lapply(res$filter, function(x) model.frame(x,data=data,na.action=na.action))
        return(list(response=y,predictor=x,group=group))
    }
    return(list(response=y,predictor=x))
}



## Specials <- function(f,spec,split2="+",...) {
##   tt <- terms(f,spec)
##   pos <- attributes(tt)$specials[[spec]]
##   if (is.null(pos)) return(NULL)
##   x <- rownames(attributes(tt)$factors)[pos]
##   st <- gsub(" ","",x)
##   res <- unlist(strsplit(st,"[()]"))[2]
##   if (is.null(split2)) return(res)
##   unlist(strsplit(res,"+",fixed=TRUE))
## }

## f <- Surv(lefttime,time,status)~x1+id(~1+z,cluster)
## spec <- "id"
## split1=","
## split2="+"


## myspecials <- c("id","strata","f")
## f <- Event(leftime,time,cause) ~ id(~1+z+z2,cluster) + strata(~s1+s2) + f(a) + z*x
## ff <- Specials(f,"id",split2=",")

Specials <- function(f,spec,split1=",",split2=NULL,...) {
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
  if (is.null(pos)) return(NULL)
  x <- rownames(attributes(tt)$factors)[pos]
  st <- gsub(" ","",x) ## trim
##  res <- unlist(strsplit(st,"[()]"))
  spec <- unlist(strsplit(st,"[()]"))[[1]]
  res <- substr(st,nchar(spec)+2,nchar(st)-1) 
  if (!is.null(split1))    
    res <- unlist(strsplit(res,split1))
  res <- as.list(res)
  for (i in seq(length(res))) {
    if (length(grep("~",res[[i]]))>0)
      res[[i]] <- as.formula(res[[i]])
  }
  return(res)
##  if (is.null(split2)) return(res)
##  strsplit(res,"+",fixed=TRUE)
}

decomp.specials <- function (x, pattern = "[()]", sep = ",", ...) 
  {
    st <- gsub(" ", "", x)
    if (!is.null(pattern)) 
      st <- rev(unlist(strsplit(st, pattern, ...)))[1]
    unlist(strsplit(st, sep, ...))
  }



