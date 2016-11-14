
qcut <- function(x,cuts=4,breaks=NULL,...)
{# {{{
	if (is.null(breaks)) {
		   probs <- seq(0,1,length.out=cuts+1)
	   bb <- quantile(x,probs,...)
	} else bb <- breaks
	gx<- cut(x,breaks=bb,include.lowest=TRUE)
	return(gx)
}# }}}

dcut <- function(x,data=NULL,cuts=4,breaks=NULL,...)
{# {{{
  if (inherits(x,"formula")) {
     x <- all.vars(x)
     xnames <- x
  } 
  if (is.character(x) & is.null(data)) stop("when x character then data frame must be given\n"); 

  if (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        n <- grep(glob2rx(xx),names(data))
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

###  data[do.call("order",c(c(x,args),list(decreasing=decreasing,method="radix"))),]

for (k in 1:ll) 
{
  xx <- x[[k]]
  if (is.null(breaks)) {
     probs <- seq(0,1,length.out=cuts[k]+1)
     bb <- quantile(xx,probs,...)
     } else bb <- breaks


if (is.null(data)) gx <- cut(xx,breaks=bb,include.lowest=TRUE)
else { ### returns data frame with new argument 
   name<-paste(xnames[k],cuts[k],sep=".")
   data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE)
}
}

if (is.null(data)) return(gx) else return(data)
}# }}}

drm <- function(x,data)
{# {{{
  if (inherits(x,"formula")) {
     x <- all.vars(x)
     xnames <- x
  } 

  if (is.character(x)) {
     xnames <- x
     xxx<-c()
     for (xx in xnames)
     {
        n <- grep(glob2rx(xx),names(data))
        xxx <- c(xxx,names(data)[n])
     }
     xnames <- xxx[!duplicated(xxx)]
  }

  data[,xnames] <- NULL

return(data)
}# }}}



