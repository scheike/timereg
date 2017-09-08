##' Life-course plot for event life data with recurrent events 
##'
##' @title Life-course plot
##' @param formula Formula (Event(start,slut,status) ~ ...) 
##' @param data data.frame
##' @param id Id variable 
##' @param group group variable
##' @param type Type (line 'l', stair 's', ...)
##' @param lty Line type
##' @param col Colour
##' @param alpha transparency (0-1)
##' @param lwd Line width
##' @param recurrent.col col of recurrence type
##' @param legend position of optional id legend
##' @param by make separate plot for each level in 'by' (formula, name of column, or vector)
##' @param status.legend Status legend
##' @param place.sl Placement of status legend
##' @param xlab Label of X-axis
##' @param ylab Label of Y-axis
##' @param add Add to existing device
##' @param ... Additional arguments to lower level arguments
##' @author Thomas Scheike Klaus K. Holst
##' @export
##' @examples
##' data = data.frame(id=c(1,1,1,2,2),start=c(0,1,2,3,4),slut=c(1,2,4,4,7),
##'                   type=c(1,2,3,2,3),status=c(0,1,2,1,2),group=c(1,1,1,2,2))
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id")
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id",recurrent.col="type")
##'
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id",group=~group,col=1:2)
##' op <- par(mfrow=c(1,2))
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id",by=~group)
##' par(op)
##' legends=c("censored","pregnant","married")
##' ll = lifecourse(Event(start,slut,status)~id,data,id="id",group=~group,col=1:2,status.legend=legends)
##'
lifecourse <- function(formula,data,id="id",group=NULL,
               type="l",lty=1,col=1:10,alpha=0.3,lwd=1,
               recurrent.col=NULL, recurrent.lty=NULL,
               legend=NULL,pchlegend=NULL,
               by=NULL,
               status.legend=NULL,place.sl="bottomright",
               xlab="Time",ylab="",add=FALSE,...) 
{# {{{
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
	status <- tstat

	X <- c(c(t1),c(t2))
	Y <- rep(x,each=2)

    }# }}}
	 plot(X,Y,xlab=xlab,ylab=ylab,...,type="n")        

       if (!is.null(status.legend)) {# {{{
	    if (is.null(pchlegend)) {
	       if (length(status.legend)!=length(unique(status)))
		   {
		           warning("Not all legends represented, legends could be wrong give pchlegend\n"); 
			   print(cbind(status.legend,sort(unique(status))))
		   }
                    points(t2,x,pch=status)
		    graphics::legend(place.sl,legend=status.legend,pch=sort(unique(status)))
	    } else { 
            points(t2,x,pch=pchlegend[status])
            graphics::legend(place.sl,legend=status.legend,pch=pchlegend)
       }
       }# }}}

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
		     recurrent.lty=recurrent.lty,
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
	    if (!is.null(recurrent.lty))   {
            cl <- data[ccm,recurrent.lty]
            matplot(t(X),t(Y),type=type,lty=cl,lwd=lwd,col=col[1],xlab=xlab,ylab=ylab,add=add,...)
	    } else matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,col=col[1],xlab=xlab,ylab=ylab,add=add,...)
	} else {
	   cn <- data[ccm,recurrent.col]
	   matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,col=cn,xlab=xlab,ylab=ylab,add=add,...)
	}

       if (!is.null(status.legend)) {# {{{
	    if (is.null(pchlegend)) {
	       if (length(status.legend)!=length(unique(status)))
		   {
		           warning("Not all legends represented, legends could be wrong give pchlegend\n"); 
			   print(cbind(status.legend,sort(unique(status))))
		   }
                    points(t2,x,pch=status)
		    graphics::legend(place.sl,legend=status.legend,pch=sort(unique(status)))
	    } else { 
            points(t2,x,pch=pchlegend[status])
            graphics::legend(place.sl,legend=status.legend,pch=pchlegend)
       }
       }# }}}


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



	if (is.null(recurrent.col))  cols <- col[1] else { 
            cols <- data[ccm,recurrent.lty]
	}
        if (is.null(recurrent.lty))  ltys <- lty else { 
            ltys <- data[ccm,recurrent.lty]
	}

###	print("hej")
###	print(t(Y))
###	print(t(X))
###	print(type)
###        matplot(t(X),t(Y),type="l",lwd=lwd,col=cols,xlab=xlab,ylab=ylab,add=add,...)

	if (is.null(recurrent.col))  {
###         matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,
###	          col=Col(col[1],alpha[1]),xlab=xlab,ylab=ylab,...)
###	          col=col[1],xlab=xlab,ylab=ylab,add=add,...)
	    if (!is.null(recurrent.lty))   {
            cl <- data[ccm,recurrent.lty]
            matplot(t(X),t(Y),type=type,lty=cl,lwd=lwd,col=col[1],xlab=xlab,ylab=ylab,add=add,...)
	    } else matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,col=col[1],xlab=xlab,ylab=ylab,add=add,...)
	} else {
            cn <- data[ccm,recurrent.col]
	    if (!is.null(recurrent.lty))   {
            cl <- data[ccm,recurrent.lty]
            matplot(t(X),t(Y),type=type,lty=cl,lwd=lwd,col=cn,xlab=xlab,ylab=ylab,add=add,...)
	    } else matplot(t(X),t(Y),type=type,lty=lty,lwd=lwd,col=cn,xlab=xlab,ylab=ylab,add=add,...)
       }

       if (!is.null(status.legend)) {# {{{
	    if (is.null(pchlegend)) {
	       if (length(status.legend)!=length(unique(status)))
		   {
		           warning("Not all legends represented, legends could be wrong give pchlegend\n"); 
			   print(cbind(status.legend,sort(unique(status))))
		   }
                    points(t2,x,pch=status)
		    graphics::legend(place.sl,legend=status.legend,pch=sort(unique(status)))
	    } else { 
            points(t2,x,pch=pchlegend[status])
            graphics::legend(place.sl,legend=status.legend,pch=pchlegend)
       }
       }# }}}

    }# }}}

###    return(invisible(data.frame(Y=Y,X=X,status=status)))
    return(invisible(list(Y,X,status)))
}# }}}
