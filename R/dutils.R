##' Cutting, sorting, rm (removing), rename for data frames
##'
##' Cut variables, if breaks are given these are used, otherwise cuts into 
##' using group size given by probs, or equispace groups on range. Default 
##' is equally sized groups if possible
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param cuts vector of number of groups, 4 is default and gives quartiles.
##' @param probs groups defined from quantiles
##' @param breaks  possible breaks for cutting.
##' @param equi.space for equi-spaced breaks  
##' @param regex for regular expressions.
##' @param sep seperator for naming of cut names.
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike
##' @examples
##' data("sTRACE",package="timereg")
##' sTRACE$age2 <- sTRACE$age^2
##' sTRACE$age3 <- sTRACE$age^3
##'
##' mm <- dcut(sTRACE,~age+wmi)
##' head(mm)
##'
##' mm <- dcut(sTRACE,~age+wmi,cuts=c(2,4))
##' head(mm)
##'
##' mm <- dcut(sTRACE,c("age","wmi"))
##' head(mm)
##'
##' mm <- dcut(sTRACE,~.)
##' head(mm)
#
##' mm <- dcut(sTRACE,c("age","wmi"),cuts=c(2,4))
##' head(mm)
##'
##' gx <- dcut(sTRACE$age)
##' head(gx)
##'
##'
##' ## Removes all cuts variables with these names wildcards
##' mm1 <- drm(mm,c("*.2","*.4"))
##' head(mm1)
##'
##' ## wildcards, for age, age2, age4 and wmi
##' head(dcut(mm,c("a*","?m*")))
##'
##' ## with direct asignment
##' drm(mm) <- c("*.2","*.4")
##' head(mm)
##'
##' dcut(mm) <- c("age","*m*")
##' head(mm)
##'
##' ## renaming
##' ############################
##'
##' head(mm)
##' drename(mm, ~age+wmi) <- c("Wmi","Age")
##' head(mm)
##' mm1 <- mm
##'
##' ## all names to lower
##' drename(mm1) <- ~.
##' head(mm1)
##'
##' ## A* to lower
##' mm2 <-  drename(mm,c("A*","W*"))
##' head(mm2)
##' drename(mm,"A*") <- ~.
##' head(mm)
##'
##' dd <- data.frame(A_1=1:2,B_1=1:2)
##' funn <- function(x) gsub("_",".",x)
##' drename(dd) <- ~.
##' drename(dd,fun=funn) <- ~.
##' names(dd)
##' @aliases dcut dcut<- dunique drm drm<- dnames dnames<- drename drename<- dkeep dkeep<- ddrop ddrop<- dreshape
##' @export
dcut <- function(data,x,cuts=4,probs=NULL,equi=FALSE,breaks=NULL,regex=FALSE,sep=NULL,...)
{# {{{
    if (is.vector(data)) {
        if (is.null(breaks)) { 
             if (!is.null(probs))
	     {
                breaks <- quantile(data, probs, ...)
	     } else {
	   	if (!equi) { 
			probs <- seq(0, 1, length.out = cuts + 1)
			breaks <- quantile(data, probs, ...)
		} 
		if (equi) { 
			rr <- range(data)
			breaks <-  seq(rr[1],rr[2],length.out=cuts+1)
		}
	     }
	}

        if (sum(duplicated(breaks))==0)
             gx <- cut(data, breaks = breaks, include.lowest = TRUE,...)
        else {
	      wd <- which(duplicated(breaks))
              mb <- min(diff(breaks[-wd]))
	      breaks[wd] <- breaks[wd] +  (mb/2)*seq(length(wd))/length(wd)
              gx  <- cut(data,breaks=breaks,include.lowest=TRUE,...)
              warning(paste("breaks duplicated"))
        }
        return(gx)
    }

if (is.data.frame(data)) {

 if (is.null(sep) & is.null(breaks))    sep <- "."
 if (is.null(sep) & (!is.null(breaks))) sep <- "b"

if (missing(x)) x<- ~.

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data)
     xnames <- x
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }


  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) {
			           if (length(x)==1 && is.character(x)) x <- data[,x]
			           x
	       })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)
  if (ll>1 & length(cuts)==1) cuts <- rep(cuts,ll)
  if (length(x)!=length(cuts)) stop("length of variables not consistent with cuts")


for (k in 1:ll)
{
  xx <- x[[k]]
  print(xnames[k])
  if (is.numeric(xx)) {

      if (is.null(breaks)) { 
             if (!is.null(probs))
	     {
                bb <- quantile(xx, probs, ...)
	     } else {
	   	if (!equi) { 
			probs <- seq(0, 1, length.out = cuts[k] + 1)
			bb <- quantile(xx, probs, ...)
		} 
		if (equi) { 
			rr <- range(xx,na.rm=TRUE)
			bb <-  seq(rr[1],rr[2],length.out=cuts[k]+1)
		}
	     }
          name<-paste(xnames[k],cuts[k],sep=sep)
      } else { bb <- breaks; name<-paste(xnames[k],breaks[2],sep=sep) }


      if (sum(duplicated(bb))==0)
	     data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE,...)
      else {
	      wd <- which(duplicated(bb))
              mb <- min(diff(bb[-wd]))
	      bb[wd] <- bb[wd] +  (mb/2)*seq(length(wd))/length(wd)
          data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE,...)
          warning(paste("breaks duplicated for=",xnames[k]))
       }
   }
}

return(data)
}

}# }}}

##' @export
"dcut<-" <- function(data,...,value) dcut(data,value,...)

##' relev levels for data frames
##'
##' levels shows levels for variables in data frame, relevel relevels a factor in data.frame 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param ref new reference variable 
##' @param regex for regular expressions.
##' @param sep seperator for naming of cut names.
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike
##' @examples
##'
##' data(mena)
##' dstr(mena)
##' dfactor(mena) <- ~twinnum
##' 
##' dstr(mena)
##' 
##' mena2 <- drelevel(mena,"cohort","(1980,1982]")
##' mena2 <- drelevel(mena,~cohort,"(1980,1982]")
##' dlevels(mena)
##' dlevels(mena2)
##' drelevel(mena,ref="(1975,1977]")  <-  ~cohort
##' drelevel(mena,ref="(1980,1982]")  <-  ~cohort
##' dlevels(mena,"coh*")
##' dtable(mena,"coh*",level=1)
##' 
##' @aliases dlevels drelevel drelevel<- dfactor
##' @export
drelevel <- function(data,x,ref=NULL,regex=FALSE,sep=NULL,...)
{# {{{

 if (missing(x) & is.data.frame(data))  stop("specify factor to relevel for data frame\n")
 if (is.null(ref)) stop("specify baseline-reference level \n")

 if (is.null(sep))  sep <- "."

 if (is.vector(data) | inherits(data,"factor")) {
      if (is.vector(data)) data <- factor(data)
      gx <- relevel(data,ref=ref)
      return(gx)
 } else {

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data)
     xnames <- x
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }


  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) {
			           if (length(x)==1 && is.character(x)) x <- data[,x]
			           x
	       })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)
  if (ll>1 & length(ref)==1) ref <- rep(ref,ll)
  if (length(x)!=length(ref)) stop("length of baseline reference 'ref' not consistent with variables")

for (k in 1:ll)
{
  xx <- x[[k]]
  if (!is.factor(xx)) xx <- factor(xx)
  name<- paste(xnames[k],ref[k],sep=sep)

  data[,name] <- relevel(xx,ref=ref[k])
}

return(data)
}

}# }}}

##' @export
"drelevel<-" <- function(data,...,value) drelevel(data,value,...)

##' @export
dlevels <- function(data,x,regex=FALSE,max.levels=20,...)
{# {{{

 if (is.factor(data)) {
	 print(base::levels(data))
 } else {

 if (missing(x)) x <-  ~.

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data)
     xnames <- x
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }


  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) {
			           if (length(x)==1 && is.character(x)) x <- data[,x]
			           x
	       })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)

for (k in 1:ll)
{
  xx <- x[[k]]
  if (is.factor(xx))  {
	  cat(paste(xnames[k],":",sep=" #levels="));
	  nxx <- nlevels(xx) 
	  cat(paste(nxx,"\n")); 
 if (is.null(max.levels) || ((!is.null(max.levels)) & (nxx<max.levels)))
	  print(base::levels(xx)) 
   }
}

}

}# }}}


##' @export
drm <- function(data,x,regex=FALSE)
{# {{{
 if (inherits(x,"formula")) {
     xnames <- all.vars(x)
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }

  data[,xnames] <- NULL

return(data)
}# }}}

##' @export
"drm<-" <- function(data,...,value) drm(data,value,...)

##' @export
drename <- function(data,var=NULL,value=NULL,fun=base::tolower,...)
{  # {{{
    if (is.null(var)) {
        var <- colnames(data)
        varpos <- seq(ncol(data))
    } else {
        var <- procform(var,data=data,return.list=FALSE,...)
        varpos <- match(var,colnames(data))
    }

    if (is.null(value)) value <- do.call(fun,list(var))

    if (is.function(value)) {
        value <- do.call(value,list(var))
    } else  if (inherits(value,"formula")) {
        value <- all.vars(value)
        if (value[1]==".") value <- do.call(fun,list(var))
    } else {
        if (length(value)==1 && (value=="." || value=="*")) value <- do.call(fun,list(var))
    }

    if (length(varpos)!= length(value)) stop("length of old and new variables must match")
    colnames(data)[varpos] <- value
    return(data)
} # }}}


##' @export
"drename<-" <- function(data,...,value) drename(data,value=value,...)

##' @export
dnames <- function(data,...) drename(data,...)

##' @export
"dnames<-" <- function(data,...,value) drename(data,value=value,...)


##' @export
dkeep <- function(data,var,keep=TRUE,regex=FALSE)
{  # {{{

 if (inherits(var,"formula")) { var <- all.vars(var) }

if (is.character(var)) {
        varxnames <- var
        varnnames <- c()
        xxx<-c()
     for (xx in varxnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
        varnnames <- c(varnnames,n)
     }
     varxnames <- xxx[!duplicated(xxx)]
     varnnames <- varnnames[!duplicated(varnnames)]
  }


if (keep) data <- data[,varnnames] else data <- data[,-1*varnnames]
  return(data)
} # }}}


##' @export
"dkeep<-" <- function(x,...,value) dkeep(x,value,...)

##' @export
ddrop <- function(data,var,keep=FALSE) dkeep(data,var,keep=FALSE)

##' @export
"ddrop<-" <- function(x,...,value) dkeep(x,value,keep=FALSE,...)

##' Sort data according to columns in data frame
##'
##' @title Sort data frame
##' @param data Data frame
##' @param x variable to order by
##' @param ... additional variables to order by
##' @param decreasing sort order (vector of length x)
##' @return data.frame
##' @export
##' @examples
##' data(data="hubble",package="lava")
##' dsort(hubble, "sigma")
##' dsort(hubble, hubble$sigma,"v")
##' dsort(hubble,~sigma+v)
##' dsort(hubble,~sigma-v)
##'
##' ## with direct asignment
##' dsort(hubble) <- ~sigma-v
##' @aliases dsort dsort<-
##' @export
dsort <- function(data,x,...,decreasing=FALSE)
{# {{{
    if (missing(x)) return(data)
    if (inherits(x,"formula")) {
        xx <- lava::procformula(value=x)$res
        decreasing <- unlist(lapply(xx,function(x) substr(trim(x),1,1)=="-"))
        x <- all.vars(x)
    }
    if (is.character(x) && length(x)<nrow(data)) x <- lapply(x,function(z) data[,z])
    dots <- list(...)
    args <- lapply(dots, function(x) {
        if (length(x)==1 && is.character(x)) x <- data[,x]
        x
    })
    if (!is.list(x)) x <- list(x)
    data[do.call("order",c(c(x,args),list(decreasing=decreasing,method="radix"))),]
}# }}}

##' @export
"dsort<-" <- function(data,...,value)  dsort(data,value,...)


##' tables for data frames
##'
##' tables for data frames
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param y name of variable, or fomula, or names of variables on data frame.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param ... Optional additional arguments
##' @param level  1 for all marginal tables, 2 for all 2 by 2 tables, and null for the full table, possible versus group variable
##' @param response For level=2, only produce tables with columns given by 'response' (index)
##' @param flat  produce flat tables
##' @param total add total counts/proportions
##' @param prop Proportions instead of counts (vector of margins)
##' @param summary summary function
##' @author Klaus K. Holst and Thomas Scheike
##' @examples
##' data("sTRACE",package="timereg")
##'
##' dtable(sTRACE,~status)
##' dtable(sTRACE,~status+vf)
##' dtable(sTRACE,~status+vf,level=1)
##' dtable(sTRACE,~status+vf,~chf+diabetes)
#'
##' dtable(sTRACE,c("*f*","status"),~diabetes)
##' dtable(sTRACE,c("*f*","status"),~diabetes, level=2)
##' dtable(sTRACE,c("*f*","status"),level=1)
##'
##' dtable(sTRACE,~status+vf+sex|age>60)
##' dtable(sTRACE,status+vf+sex~+1|age>60, level=2)
##' dtable(sTRACE,.~status+vf+sex|age>60,level=1)
##' dtable(sTRACE,status+vf+sex~diabetes|age>60)
##' dtable(sTRACE,status+vf+sex~diabetes|age>60, flat=FALSE)
##'
##' dtable(sTRACE,status+vf+sex~diabetes|age>60, level=1)
##' dtable(sTRACE,status+vf+sex~diabetes|age>60, level=2)
##'
##' dtable(sTRACE,status+vf+sex~diabetes|age>60, level=2, prop=1, total=TRUE)
##' dtable(sTRACE,status+vf+sex~diabetes|age>60, level=2, prop=2, total=TRUE)
##' dtable(sTRACE,status+vf+sex~diabetes|age>60, level=2, prop=1:2, summary=summary)
##'
##' @aliases dtables
##' @export
dtable <- function(data,y=NULL,x=NULL,...,level=-1,response=NULL,flat=TRUE,total=FALSE,prop=FALSE,summary=NULL) {
    daggregate(data,y,x,...,
               fun=function(z) {
                   res <- sum <- c()
                   if (level==1 || ncol(z)==1) {
                       for (i in seq_len(ncol(z))) {
                           nn <- colnames(z)[i]
                           val <- table(z[,i],...)
                           names(attr(val,"dimnames")) <- nn
                           val <- list(val)
                           names(val) <- nn
                           res <- c(res, val)
                           if (!is.null(summary)) {
                               sval <- list(do.call(summary,list(val[[1]])))
                               names(sval) <- nn
                               c(sum, sval)
                           }
                       }
                       res <- list(table=res,summary=sum,...)
                       class(res) <- "dtable"
                       return(res)
                   }
                   if (level>1) {
                       if (level>2) response <- ncol(z)
                       idx1 <- seq(1,ncol(z)-1)
                       if (level>2 || !is.null(response)) {
                           idx1 <- response
                           idx2 <- setdiff(seq(ncol(z)),idx1)
                       }
                       for (i in idx1)  {
                           if (!(level>2 || !is.null(response))) {
                               idx2 <- seq(i+1,ncol(z))
                           }
                           for (j in idx2) {
                               n1 <- colnames(z)[j]
                               n2 <- colnames(z)[i]
                               val <- table(z[,c(j,i)],...)
                               if (prop[1]>0) {
                                   if (all(1:2 %in% prop)) {
                                       val <- prop.table(val)
                                   } else {
                                       val <- prop.table(val,prop)
                                   }
                               }
                               if (total) {
                                   tot <- prop
                                   if (length(prop)==1) tot <- setdiff(1:2,prop)
                                       val <- addmargins(val,tot)
                               }
                               val <- list(val)
                               names(val) <- paste0(n1,", ",n2)
                               res <- c(res, val)
                               if (!is.null(summary)) {
                                   sval <- list(do.call(summary,list(val[[1]])))
                                   #names(sval) <- names(val)
                                   sum <- c(sum, sval)
                               }
                           }
                       }
                       res <- list(table=res,summary=sum)
                       class(res) <- "dtable"
                       return(res)
                   }

                   res <- table(z,...)
                   if (!is.null(summary)) {
                       sum <- do.call(summary,c(list(res),list(...)))
                   }
                   if (prop[1]>0) res <- prop.table(res,prop)
                   if (total>0) res <- addmargins(res,prop)
                   if (flat) res <- ftable(res,...)
                   res <- list(table=res,summary=sum)
                   class(res) <- "dtable"
                   return(res)
               })
}

##' @export
print.dtable <- function(x,sep="\n",...) {
    cat(sep)
    if (inherits(x$table, c("table","ftable"))) {
        print(x$table)
        if (!is.null(x$summary)) print(x$summary)
        return(invisible(x))
    }

    for (i in seq_along(x$table)) {
        print(x$table[[i]],...)
        if (!is.null(x$summary))
            print(x$summary[[i]],...)
        cat(sep)
    }

}


procform <- function(formula, sep, nsep=1, return.formula=FALSE, data=NULL, regex=FALSE, return.list=TRUE, ...) {
    res <- NULL
    if (is.null(formula)) {
        res <- colnames(data)
    } else if (is.character(formula)) {
        if (is.null(data)) {
            res <- unique(formula)
        } else {
            yy <- c()
            for (y0 in formula) {
                if (!regex) y0 <- glob2rx(y0)
                n <- grep(y0,names(data))
                yy <- union(yy,names(data)[n])
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
    aa <- attributes(terms(formula,data=data))
    if (aa$response == 0) {
        res <- NULL
    } else {
        res <- paste(deparse(formula[[2]]), collapse = "")
    }
    filter.expression <- NULL
    foundsep <- FALSE
    pred <- filter <- c()
    if (!missing(sep) && length(aa$term.labels) > 0) {
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
            filter.expression <- parse(text=filter)
            filter <- as.list(filter)
        }
    }
    if (!foundsep) pred <- aa$term.labels
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

procformdata <- function(formula,data,sep="\\|", na.action=na.pass, ...) {
    res <- procform(formula,sep=sep,data=data,return.formula=TRUE,...)
    y <- x <- NULL
    filter <- res$filter.expression

    if (length(res$response)>0) {
        if (is.null(filter)) y <- model.frame(res$response,data=data,na.action=na.action)
        else y <- model.frame(res$response,data=subset(data,eval(filter)),na.action=na.action)
    }
    if (length(res$predictor)>0) {
        if (is.null(filter)) x <- model.frame(res$predictor,data=data,na.action=na.action)
        else x <- model.frame(res$predictor,data=subset(data,eval(filter)),na.action=na.action)

    }
    ## if (!is.null(res$group)) group <- lapply(res$,function(x) model.frame(x,data=data,...))
    return(list(response=y,predictor=x))
}


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
##' daggregate(iris, log(Sepal.Length)+I(Petal.Length>1.5) ~ Species, fun=summary)
##' daggregate(iris, "*Length*", x="Species", fun=head)
##' daggregate(iris, "^.e.al", x="Species", fun=tail, regex=TRUE)
##' daggregate(sTRACE, status~ diabetes, fun=table)
##' daggregate(sTRACE, status~ diabetes+sex, fun=table)
##' daggregate(sTRACE, status + diabetes+sex ~ vf+I(wmi>1.4), fun=table)
##' daggregate(iris, "^.e.al", x="Species",regex=TRUE)
##' dprint(iris,Petal.Length+Sepal.Length ~ Species |Petal.Length>1.3 & Sepal.Length>5, n=list(1:3,-(3:1)))
##' daggregate(iris, I(Sepal.Length>7)~Species | I(Petal.Length>1.5))
##' daggregate(iris, I(Sepal.Length>7)~Species | I(Petal.Length>1.5), fun=table)
##'
##' dsum(iris, .~Species, matrix=TRUE, missing=TRUE)
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
daggregate <- function(data,y=NULL,x=NULL,subset,...,fun="summary",regex=FALSE, missing=FALSE, remove.empty=FALSE, matrix=FALSE, silent=FALSE, na.action=na.pass)
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
        yx <- procformdata(y,sep="\\|",data=data,na.action=na.action,...)
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
            n <- grep(y0,names(data))
            yy <- union(yy,names(data)[n])
        }
        y <- data[,yy,drop=FALSE]
    }
    if (is.character(x) && length(x)<NROW(data)) {
        xx <- c()
        for (x0 in x) {
            if (!regex) x0 <- glob2rx(x0)
            n <- grep(x0,names(data))
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

    if (!is.null(x)) {
        if (missing) {
            x[is.na(x)] <- 'NA'
        }
        if (silent)
            capture.output(res <- by(y,x,fun,...))
        else
            res <- by(y,x,fun,...)
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
print.daggregate <- function(x,quote=FALSE,...) {
    attr(x,c("ngroupvar")) <- NULL
    class(x) <- class(x)[-1]
    print(x,quote=quote,...)
}


##' @export
dfactor <- function(data,x,regex=FALSE,sep=NULL,all=FALSE,...)
{# {{{

 if (is.null(sep))  sep <- ".f"

 if (is.vector(data)) {
	 if (!is.factor(data)) gx <- as.factor(data) else gx <- data
      return(gx)
 } else {

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data)
     xnames <- x
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }
 print(xnames)


  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) {
			           if (length(x)==1 && is.character(x)) x <- data[,x]
			           x
	       })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)
###  if (ll>1 & length(ref)==1) ref <- rep(ref,ll)
###  if (length(x)!=length(ref)) stop("length of baseline reference 'ref' not consistent with variables")

for (k in 1:ll)
{
  xx <- x[[k]]
  name<- paste(xnames[k],sep,sep="")
  if (!is.factor(xx) || all==TRUE) { 
	  gx <- as.factor(xx); 
          data[,name] <- gx } else data[,name] <- as.factor(xx) 
}

return(data)
}

}# }}}

##' @export
"dfactor<-" <- function(data,k=1,combine=TRUE,...,value) {
    dfactor(data,value,...)
}

##' @export
dnumeric <- function(data,x,regex=FALSE,sep=NULL,all=FALSE,...)
{# {{{

 if (is.null(sep))  sep <- ".n"

 if (is.factor(data)) {
      gx <- as.numeric(data) 
      return(gx)
 } else {

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data)
     xnames <- x
  } else if  (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        if (!regex) xx <- glob2rx(xx)
        n <- grep(xx,names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }
 print(xnames)


  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) {
			           if (length(x)==1 && is.character(x)) x <- data[,x]
			           x
	       })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)
###  if (ll>1 & length(ref)==1) ref <- rep(ref,ll)
###  if (length(x)!=length(ref)) stop("length of baseline reference 'ref' not consistent with variables")

for (k in 1:ll)
{
  xx <- x[[k]]
  name<- paste(xnames[k],sep,sep="")
  if (!is.numeric(xx) || all==TRUE) { 
	  gx <- as.numeric(xx); 
          data[,name] <- gx } else data[,name] <- as.numeric(xx) 
}

return(data)
}

}# }}}

##' @export
"dnumeric<-" <- function(data,...,value) {
    dnumeric(data,value,...)
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

###d <- data.frame(date=as.Date(1:30,origin="1960-01-01"),
###               g=factor(rep(letters[1:3],each=10)),
###               y=rbinom(30,1,0.5),z=rnorm(30),x=rep(c(0,1),15))


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
    ##x
    xx
}

##' @export
deval2 <- function(data,...,matrix=simplify,simplify=TRUE)  deval(data,matrix=TRUE,simplify=TRUE,...)

##' @export
deval <- function(data,y=NULL,x=NULL,...,matrix=FALSE,fun=Summary,simplify=FALSE) {
    res <- daggregate(data,y,x,matrix=matrix,...,
                     fun=function(z) lapply(z,function(x) {
                         suppressWarnings(tryCatch(fun(x,...),error=function(e) return(NA)))
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


##' list, head, print, tail 
##'
##' listing for data frames
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param y name of variable, or fomula, or names of variables on data frame.
##' @param n Index of observations to print (default c(1:nfirst, n-nlast:nlast)
##' @param ... Optional additional arguments (nfirst,nlast, and print options)
##' @param x possible group variable
##' @author Klaus K. Holst and Thomas Scheike
##' @examples
##' n <- 20
##' m <- lava::lvm(letters)
##' d <- lava::sim(m,n)
##'
##' dprint(d,~a+b+c)
##' dprint(d,~a+b+c|a<0 & b>0)
##' ## listing all : 
##' dprint(d,~a+b+c|a<0 & b>0,n=0)
##' dprint(d,a+b+c~I(d>0)|a<0 & b>0)
##' dprint(d,.~I(d>0)|a<0 & b>0)
##' dprint(d,~a+b+c|a<0 & b>0, nlast=0)
##' dprint(d,~a+b+c|a<0 & b>0, nfirst=3, nlast=3)
##' dprint(d,~a+b+c|a<0 & b>0, 1:5)
##' dprint(d,~a+b+c|a<0 & b>0, -(5:1))
##' dprint(d,~a+b+c|a<0 & b>0, list(1:5,50:55,-(5:1)))
##' dlist(d,a+b+c ~ I(d>0) |a<0 & b>0, list(1:5,50:55,-(5:1)))
##' @aliases dprint dlist dhead dtail 
##' @export
dprint <- function(data,y=NULL,n=NULL,...,x=NULL) daggregate(data,y,x,...,fun=function(z,...) Print(z,n=n,...),silent=FALSE)

##' @export
dlist <- function(data,...) dprint(data,...)

##' @export
dreshape <- function(data,...) fast.reshape(data,...)


##' Lag operator
##'
##' Lag operator
##' @examples
##' d <- data.frame(y=1:10,x=c(10:1))
##' dlag(d,k=1:2)
##' dlag(d,~x,k=0:1)
##' dlag(d$x,k=1)
##' dlag(d$x,k=-1:2, names=letters[1:4])
##' @export
##' @param data data.frame or vector
##' @param x optional column names or formula
##' @param k lag (vector of integers)
##' @param combine combine results with original data.frame
##' @param simplify Return vector if possible
##' @param names optional new column names
##' @param ... additional arguments to lower level functions
##' @aliases dlag dlag<-
dlag <- function(data,x,k=1,combine=TRUE,simplify=TRUE,names,...) {
    isvec <- FALSE
    if (!is.data.frame(data)) {
        isvec <- is.vector(data)
        data <- as.data.frame(data)
        combine <- FALSE
    }
    if (missing(x)) x <- base::names(data)
    if (is.character(x)) x <- data[,x,drop=FALSE]
    if (inherits(x,"formula")) {
        ##x <- as.data.frame(model.matrix(update(x,~.-1), model.frame(~.,data=data, na.action=na.pass)))
        x <- model.frame(x,data=data, na.action=na.pass)
    }
    kmin0 <- abs(min(k,0))
    kmax0 <- max(k,0)
    kmax <- kmin0+kmax0
    pad <- function(x) c(rep(NA,kmax0),x,rep(NA,kmin0))
    kidx <- match(k,seq(min(k,0),max(k,0)))
    val <- lapply(x,function(y) embed(pad(y),dimension=kmax+1)[,kidx,drop=FALSE])
    dval <- as.data.frame(val)    
    if (!missing(names)) {
        base::names(dval) <- names
    } else {
        nn <- as.vector(sapply(colnames(x),function(x) paste0(x,paste0(".",k))))        
        if (length(nn)==ncol(dval)) {
            nn <- gsub("-","_",nn)
            base::names(dval) <- nn
        }
    }
    if (combine) {        
        res <- cbind(data,dval)
        names(res) <- make.unique(base::names(res))
        return(res)
    }
    if (length(k)==1 && simplify && isvec) return(as.vector(val[[1]]))
    names(dval) <- base::make.unique(base::names(dval))
    return(dval)
}

##' @export
"dlag<-" <- function(data,k=1,combine=TRUE,...,value) {
    dlag(data,value,k=k,combine=combine,...)
}



Print <- function(x,n=NULL,nfirst=5,nlast=nfirst,digits=max(3,getOption("digits")-3),...) {
    mat <- !is.null(dim(x))
    if (!mat) {
        x <- cbind(x)
        colnames(x) <- ""
    }
    if (is.null(n)) {
        if (NROW(x)<=(nfirst+nlast)) n <- list(seq(NROW(x)))
        else {
            n <- c()
            if (nfirst>0)
                n <- c(n,list(seq(nfirst)))
            if (nlast>0)
                n <- c(n,list(-rev(seq(nlast))))
        }
    }
    if (is.null(n) || !is.list(n) && length(n)==1 && n==0) return(x)
    if (!is.list(n)) n <- list(n)
    d <- lapply(n,function(idx) {
        N <- NROW(x)
        idx <- idx[idx!=0 & abs(idx)<=N]
        idx[idx<0] <- N+idx[idx<0]+1
        base::format(x[idx,,drop=FALSE],digits=digits,...)
    })
    val <- c()
    sep <- rbind("---"=rep('',ncol(x)))
    for (i in seq_along(d)) {
        if (i>1) val <- rbind(val,sep)
        val <- rbind(val,base::as.matrix(d[[i]]))

    }
    return(structure(val,class=c("Print",class(val))))
}

##' @export
print.Print <- function(x,quote=FALSE,...) {
    class(x) <- class(x)[-1]
    print(x,quote=quote,...)
}



