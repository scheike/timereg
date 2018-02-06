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
##'
##' dtable(sTRACE,c("*f*","status"),~diabetes)
##' dtable(sTRACE,c("*f*","status"),~diabetes, level=2)
##' dtable(sTRACE,c("*f*","status"),level=1)
##'
##' dtable(sTRACE,~"*f*"+status,level=1)
##' dtable(sTRACE,~"*f*"+status+I(wmi>1.4)|age>60,level=2)
##' dtable(sTRACE,"*f*"+status~I(wmi>0.5)|age>60,level=1)
##' dtable(sTRACE,status~dcut(age))
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
##' @aliases dtable dtab
##' @export
dtable <- function(data,y=NULL,x=NULL,...,level=-1,response=NULL,flat=TRUE,total=FALSE,prop=FALSE,summary=NULL) {
       	daggregate(data,y,x,...,
               fun=function(z) {
                   res <- sum <- c()
                   if (level==1 || ncol(z)==1) {
                       for (i in seq_len(ncol(z))) {
                           nn <- colnames(z)[i]
                           val <- table(z[,i],...)
			   if (prop[1]>0) val <- prop.table(val)
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
dtab <- function(data,y=NULL,x=NULL,...) dtable(data,y=NULL,x=NULL,...)

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

