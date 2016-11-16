##' Cutting, sorting, rm (removing), rename for data frames 
##' 
##' Cut variables, rm or rename 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param cuts vector of number of groups, 4 is default and gives quartiles.
##' @param breaks  possible breaks for cutting.
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
##' mm <- dcut(sTRACE,c("age","wmi"),cuts=c(2,4))
##' head(mm)
##'
##' gx <- qcut(sTRACE$age)
##' head(gx)
##'
##'
##' ## Removes all cuts variables with these names wildcards
##' mm <- drm(mm,c("*.2","*.4"))
##' head(mm)
##'
##' ## wildcards, for age, age2, age4 and wmi 
##' mm <- dcut(mm,c("a*","?m*"))
##' head(mm)
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
##' ## to lower 
##' head(drename(mm))
##' mm1 <-  drename(mm)
##' mm1 <-  drename(mm, ~age+wmi,~wmi+age)
##' head(mm1)
##'
##' drename(mm, ~age+wmi) <- c("wmi","age")
##' head(mm)
##' drename(mm, ~wmi+age) <- c("age","wmi")
##' head(mm)
##' #mm <-  drenname(mm,"age*","alder*")
##' #head(mm)
##' @export
dcut <- function(data,x,cuts=4,breaks=NULL,...)
{# {{{
 if (inherits(x,"formula")) {
     x <- all.vars(x)
     xnames <- x
     formular <- 1
  } else if  (is.character(x)) {
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
   name<-paste(xnames[k],cuts[k],sep=".")
   data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE)
}

###if (is.null(data)) gx <- cut(xx,breaks=bb,include.lowest=TRUE)
###else { ### returns data frame with new argument 
###   name<-paste(xnames[k],cuts[k],sep=".")
###   data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE)
###}


###if (is.null(data)) return(gx) else 
return(data)
}# }}}

##' @export
"dcut<-" <- function(data,...,value)
{# {{{
 dcut(data,value,...)
}# }}}

##' @export
drm <- function(data,x)
{# {{{
 if (inherits(x,"formula")) {
     xnames <- all.vars(x)
  } else if  (is.character(x)) {
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

##' @export
"drm<-" <- function(data,...,value)
{# {{{
 drm(data,value,...)
}# }}}

##' @export
drename <- function(data,var,value) 
{  # {{{

 if (missing(var)) var <- names(data)

 if (inherits(var,"formula")) {
      var <- all.vars(var)
   }

if (is.character(var)) {
        varxnames <- var 
        varnnames <- c()
        xxx<-c()
     for (xx in varxnames)
     {
        n <- grep(glob2rx(xx),names(data))
        xxx <- c(xxx,names(data)[n])
        varnnames <- c(varnnames,n) 
     }
     varxnames <- xxx[!duplicated(xxx)]
     varnnames <- varnnames[!duplicated(varnnames)]
  }

 if (missing(value)) value <- tolower(varxnames)

 if (inherits(value,"formula")) {
     value <- all.vars(value)
 }

   if (is.character(value)) {
        xnames <- value  
###        xxx<-c()
###     for (xx in xnames)
###     {
###        n <- grep(glob2rx(xx),names(data))
###        xxx <- c(xxx,names(data)[n])
###     }
###         xnames <- xxx[!duplicated(xxx)]
###         nnames <- nnames[!duplicated(nnames)]
  }

 print(xnames)
 print(varxnames)

  if (length(xnames)!= length(varxnames)) stop("length of old and new variables must mach")

  names(data)[varnnames] <- xnames; 
  return(data) 
} # }}}

##' @export
"drename<-" <- function(x, var, value) 
{ # {{{
  drename(x,var,value)
}# }}}




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
dsort <- function(data,x,...,decreasing=FALSE) 
{# {{{
    if (missing(x)) return(data)
    if (inherits(x,"formula")) {
        xx <- procformula(value=x)$res
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
"dsort<-" <- function(data,...,value) 
{ # {{{
  dsort(data,value,...)
}# }}}


