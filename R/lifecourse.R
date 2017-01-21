##' Life-course plot for event life data  with recurrent events 
##'
##' @title Life course plot
##' @param formula Formula (Event(start,slut,status) ~ ) 
##' @param data data.frame
##' @param id Id variable 
##' @param group group variable
##' @param type Type (line 'l', stair 's', ...)
##' @param lty Line type
##' @param col Colour
##' @param alpha transparency (0-1)
##' @param lwd Line width
##' @param recurrent.col col of recurrence type 
##' @param by make separate plot for each level in 'by' (formula, name of column, or vector)
##' @param xlab Label of X-axis
##' @param ylab Label of Y-axis
##' @param add Add to existing device
##' @param ... Additional arguments to lower level arguments
##' @author Thomas Scheike Klaus K. Holst
##' @examples
##' data = data.frame(id=c(10,10,10,2,2),start=c(0,1,2,3,4),slut=c(1,2,4,4,7),
##'                   type=c(1,2,3,2,3),status=c(0,1,2,1,2),group=c(1,1,1,2,2))
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id")
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id",recurrent.col="type")
##'
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id",group=~group,col=1:2)
##' par(mfrow=c(1,2))
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id",by=~group)
##' 
##' data$gid = data$g*10+ data$id/5
##' ll = lifecourse(Event(start,slut,status)~+gid,data,id="id")
##' ll = lifecourse(Event(start,slut,status)~+1,data,id="id")
##' ll = lifecourse(Event(start,slut,status)~+1,data,id="id",group=~group,col=1:2)
##' 
##' @export
lifecourse <- function(formula,data,id="id",group=NULL,
                      type="l",lty=1,col=1:10,alpha=0.3,lwd=1,
		      recurrent.col=NULL,
                      legend=NULL, by=NULL,
                      xlab="Time",ylab="Id",add=FALSE,...) 
{# {{{
    if (!lava.options()$cluster.index) stop("mets not available? Check 'lava.options()cluster.index'.")
    if (!is.null(by)) {
        if (is.character(by) && length(by==1)) {
            by <- data[,by]
        } else if (inherits(by,"formula")) {
            by <- model.frame(by,data,na.action=na.pass)
        }
        cl <- match.call(expand.dots=TRUE)
        cl$by <- NULL
        datasets <- split(data,by)
        res <- c()
        for (d in datasets) {
            cl$data <- d
            res <- c(res, eval(cl,parent.frame()))
        }
        return(invisible(res))
    }
    if (!is.null(group)) {
        if (is.character(group) && length(group==1)) {
            M <- data[,group]
        } else if (inherits(group,"formula")) {
            M <- model.frame(group,data,na.action=na.pass)
        } else {
            M <- group
        }
###        if (!add) plot(formula,data=data,xlab=xlab,ylab=ylab,...,type="n")        
        if (!add) {
	 if (inherits(id,"formula")) id <- all.vars(id)
	 if (inherits(group,"formula")) group <- all.vars(group)
	 if (is.character(id) && length(id)==1) Id <- 
	 y <- getoutcome(formula)
	 x <- attributes(y)$x

    if (length(x)==0) {# {{{

        y <- response <- all.vars( update(formula,.~+1))
        ###
	ccid <- cluster.index(data[,id])
	ccm <- ccid$idclustmat+1
        ccm <- ccm[!is.na(ccm)]
        ###
        x <- data[ccm,id]

	if (length(response)==3) { 
	t1 <-    data[ccm,response[1]] 
	t2 <-    data[ccm,response[2]] 
	tstat <- data[ccm,response[3]]
	} else {
	t1 <-    rep(0,length(ccm))
	t2 <-    data[ccm,response[1]] 
	tstat <- data[ccm,response[2]]
	}

	X <- c(c(t1),c(t2))
	Y <- rep(x,each=2)

# }}}
    } else {# {{{
        y <- response <- all.vars( update(formula,.~+1))
        ###
	ccid <- cluster.index(data[,id])
	ccm <- ccid$idclustmat+1
	ccm <- ccm[!is.na(ccm)]

	if (length(response)==3) { 
	t1 <-    data[ccm,response[1]] 
	t2 <-    data[ccm,response[2]] 
	tstat <- data[ccm,response[3]]
	} else {
	t1 <-    rep(0,length(ccm))
	t2 <-    data[ccm,response[1]] 
	tstat <- data[ccm,response[2]]
	}

        x <- data[ccm,x] 


	X <- c(c(t1),c(t2))
	Y <- rep(x,each=2)

    }# }}}
	 plot(X,Y,xlab=xlab,ylab=ylab,...,type="n")        
	}

        dd <- split(data,M)
        K <- length(dd)
        if (length(type)<K)        type <- rep(type,K)
        if (length(col)<K)         col <- rep(col,K)
        if (length(lty)<K)         lty <- rep(lty,K)
        if (length(lwd)<K)         lwd <- rep(lwd,K)
        if (length(alpha)<K)       alpha <- rep(alpha,K)
	res <- list()
        for (i in seq_len(K)) {            
           res[[i]] <- lifecourse(formula,data=dd[[i]],id=id,type=type[i],
                     lty=lty[i],col=col[i],lwd=lwd[i],
                     alpha=alpha[i],
		     recurrent.col=recurrent.col,
                     group=NULL,
                     add=TRUE,...)
        }
        if (!is.null(legend)) {
            graphics::legend(legend,names(dd),lwd=lwd,col=col,lty=lty)
        }
        return(invisible(res))
    }
    
    if (inherits(id,"formula")) id <- all.vars(id)
    if (inherits(group,"formula")) group <- all.vars(group)
    if (is.character(id) && length(id)==1) Id <- 
    y <- getoutcome(formula)
    x <- attributes(y)$x

    if (!requireNamespace("mets",quietly=TRUE)) stop("'mets' package required")

    if (length(x)==0) {# {{{

        y <- response <- all.vars( update(formula,.~+1))
        ###
	ccid <- cluster.index(data[,id])
	ccm <- ccid$idclustmat+1
        ccm <- ccm[!is.na(ccm)]
        ###
        x <- data[ccm,id]

	if (length(response)==3) { 
	t1 <-    data[ccm,response[1]] 
	t2 <-    data[ccm,response[2]] 
	tstat <- data[ccm,response[3]]
	} else {
	t1 <-    rep(0,length(ccm))
	t2 <-    data[ccm,response[1]] 
	tstat <- data[ccm,response[2]]
	}

	X <- cbind(c(t1),c(t2))
	Y <- matrix(rep(x,each=2),ncol=2,byrow=TRUE)
	status <- tstat

###	print(dim(X)); print(dim(Y)); print(dim(status)); print(summary(X)); print(summary(Y)); 

	if (is.null(recurrent.col))  {
            matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,
	          col=col[1],xlab=xlab,ylab=ylab,add=add,..)
	} else {
	   cn <- data[ccm,recurrent.col]
	   matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,
	       col=cn,xlab=xlab,ylab=ylab,add=add,...)
		}
         points(t2,x,pch=status)
# }}}
    } else {# {{{
        y <- response <- all.vars( update(formula,.~+1))
        ###
	ccid <- cluster.index(data[,id])
	ccm <- ccid$idclustmat+1
	ccm <- ccm[!is.na(ccm)]

	if (length(response)==3) { 
	t1 <-    data[ccm,response[1]] 
	t2 <-    data[ccm,response[2]] 
	tstat <- data[ccm,response[3]]
	} else {
	t1 <-    rep(0,length(ccm))
	t2 <-    data[ccm,response[1]] 
	tstat <- data[ccm,response[2]]
	}

        x <- data[ccm,x] 

	X <- cbind(c(t1),c(t2))
	Y <- matrix(rep(x,each=2),ncol=2,byrow=TRUE)
	status <- tstat

###	print(dim(X)); print(dim(Y)); print(dim(status)); print(summary(X)); print(summary(Y)); 

	if (is.null(recurrent.col))  {
         matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,
###	          col=Col(col[1],alpha[1]),xlab=xlab,ylab=ylab,...)
	          col=col[1],xlab=xlab,ylab=ylab,add=add,...)
	} else {
       cn <- data[ccm,recurrent.col]
       matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,
	       col=cn,xlab=xlab,ylab=ylab,add=add,...)
	}
         points(t2,x,pch=status)
    }# }}}

###    return(invisible(data.frame(Y=Y,X=X,status=status)))
    return(invisible(list(Y,X,status)))
}# }}}
