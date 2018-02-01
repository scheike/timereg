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

##' @export
procform <- function(formula=NULL, sep="\\|", nsep=1, return.formula=FALSE, data=NULL,
             no.match=TRUE, regex=FALSE, return.list=TRUE, specials=NULL,...) {# {{{
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
            if (filter%in%c("0","-1")) {
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
    if (any(res==".")) {
        diffset <- c(".",setdiff(pred,res))
        res <- setdiff(union(res,colnames(data)),diffset)
    }
    filter <- lapply(filter, expandst)
    if (!is.null(specials)) {
        foundspec <- replicate(length(specials),c())
        names(foundspec) <- specials
        rmidx <- c()
        spec <- paste0("^",specials,"\\(")
        val <- lapply(spec, function(x) which(grepl(x,pred)))
        for (i in seq_along(val)) {
            if (length(val[[i]])>0) { # special function found
                rmidx <- c(rmidx,val[[i]])
                cleanpred <- gsub("\\)$","",gsub(spec[i],"",pred[val[[i]]]))
                foundspec[[i]] <- c(foundspec[[i]],cleanpred)
            }
        }
        if (length(rmidx)>0)
            pred <- pred[-rmidx]
        if (length(pred)==0) pred <- NULL
        specials <- foundspec
        for (i in seq_along(specials)) if (is.null(specials[[i]])) specials[i] <- NULL
        if (length(specials)==0) specials <- NULL
    }

    if (return.formula) {
        if (foundsep && !is.null(filter)) {
            filter <- lapply(filter, function(z) as.formula(paste0(c("~", paste0(z,collapse="+")))))
        }
        if (length(pred)>0)
            pred <- as.formula(paste0(c("~", paste0(pred,collapse="+"))))
        if (length(res)>0)
            res <- as.formula(paste0(c("~", paste0(res,collapse="+"))))
        if (!is.null(specials)) {
            specials <- lapply(specials,function(x)
                              as.formula(paste0(c("~", paste0(x,collapse="+")))))
        }
    }
    res <- list(response=res, predictor=pred, filter=filter, filter.expression=filter.expression, specials=specials)
    if (!return.list) return(unlist(unique(res)))
    return(res)
}

##' @export
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
    ### print(filter); print(missing(filter)); print(is.null(filter)); print(filter[[1]])
    ### when filter.expression is expression(1) then also no filter, ts 
    if ((!missing(filter))) if (!is.null(filter)) if (as.character(filter)=="1") filter <- NULL

    if (length(res$response)>0) {
        if (is.null(filter)) y <- model.frame(res$response,data=data,na.action=na.action)
        else y <- model.frame(res$response,data=subset(data,eval(filter)),na.action=na.action)
    }
    if (length(res$predictor)>0) {
        if (is.null(filter)) x <- model.frame(res$predictor,data=data,na.action=na.action)
        else x <- model.frame(res$predictor,data=subset(data,eval(filter)),na.action=na.action)

    }

    specials <- NULL
    if (!is.null(res$specials)) {
        specials <- lapply(res$specials,
                      function(x) {
                          if (is.null(filter)) model.frame(x,data=data,na.action=na.action)
                          else model.frame(x,data=subset(data,eval(filter)),na.action=na.action)
                      })
    }

    if (!do.filter) {
        group <- lapply(res$filter, function(x) model.frame(x,data=data,na.action=na.action))
        return(list(response=y,predictor=x,group=group,specials=specials))
    }
    return(list(response=y,predictor=x,specials=specials))
}# }}}

procform2 <- function(y,x=NULL,z=NULL,...) {# {{{
    yx <- procform(y,return.formula=FALSE,...)
    y <- yx$response
    x0 <- yx$predictor
    z0 <- NULL
    if (length(yx$filter)>0) z0 <- yx$filter[[1]]
    if (is.null(x) && length(y)>0) x <- x0
    if (NCOL(x)==0) x <- NULL
    if (length(y)==0) y <- x0
    if (!is.null(x)) {
        x <- unlist(procform(x,return.formula=FALSE,...))
    }
    if (is.null(z)) z <- z0
    if (!is.null(z)) {
        zz <- unlist(procform(z,return.formula=FALSE,...))
    }
    if (is.null(z)) z <- z0
    return(list(y=y,x=x,z=z))
}# }}}

##' @export
procform3 <- function(y,x=NULL,z=NULL,...) {# {{{
    yx <- procform(y,return.formula=FALSE,...)
    x0 <- yx$predictor
    y  <- yx$response
    if (is.null(yx$predictor)) { x0 <- yx$response ; y <- NULL} 

    if (is.null(y)) 
    if (!is.null(x)) {
        x <- procform(x,return.formula=FALSE,...)
        y <- c(x$predictor,x$response)
    }
    return(list(y=y,x=x0,z=z))
}# }}}


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

Specials <- function(f,spec,split1=",",split2=NULL,...) {# {{{
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
}# }}}

decomp.specials <- function (x, pattern = "[()]", sep = ",", ...)
  {
    st <- gsub(" ", "", x)
    if (!is.null(pattern))
      st <- rev(unlist(strsplit(st, pattern, ...)))[1]
    unlist(strsplit(st, sep, ...))
  }

