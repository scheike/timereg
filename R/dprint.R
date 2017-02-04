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
##' dlist(d,~a+b+c)
##' dlist(d,~a+b+c|a<0 & b>0)
##' ## listing all : 
##' dlist(d,~a+b+c|a<0 & b>0,n=0)
##' dlist(d,a+b+c~I(d>0)|a<0 & b>0)
##' dlist(d,.~I(d>0)|a<0 & b>0)
##' dlist(d,~a+b+c|a<0 & b>0, nlast=0)
##' dlist(d,~a+b+c|a<0 & b>0, nfirst=3, nlast=3)
##' dlist(d,~a+b+c|a<0 & b>0, 1:5)
##' dlist(d,~a+b+c|a<0 & b>0, -(5:1))
##' dlist(d,~a+b+c|a<0 & b>0, list(1:5,50:55,-(5:1)))
##' dprint(d,a+b+c ~ I(d>0) |a<0 & b>0, list(1:5,50:55,-(5:1)))
##' @aliases dprint dlist dhead dtail 
##' @export
dprint <- function(data,y=NULL,n=0,...,x=NULL) daggregate(data,y,x,...,fun=function(z,...) Print(z,n=n,...),silent=FALSE)

##' @export
dlist <- function(data,y=NULL,n=NULL,...) dprint(data,y=y,n=n,...)
