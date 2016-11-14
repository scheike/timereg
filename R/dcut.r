##' Cutting, rm (removing), renaming  for data frames 
##' 
##' Cut variables, rm or rename 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param cuts vector of number of groups, 4 is default and gives quartiles.
##' @param breaks  possible breaks for cutting.
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike 
##' @examples
##' data(sTRACE)
##' sTRACE$age2 <- sTRACE$age^2
##' sTRACE$age3 <- sTRACE$age^3
##'
##' mm <- dcut(sTRACE,~age+wmi,sTRACE)
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
##' dren(data, ~x+y) <- c("x","y")
##' dren(data, ~x+y) <- c("x","y")
##' drm(data) <- "*.2"
##'
##' "dren<-" <- function(x, var, value) { browser(); names(x)[var] <- value; return(x) }
##' dren <- function(x,var,value) {  } 
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
##'
##'
##' @export
dcut <- function(data,x,cuts=4,breaks=NULL,...)
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

##' @export
"dcut<-" <- function(data,x,cuts=4,breaks=NULL,...)
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

##' @export
drm <- function(data,x)
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

##' @export
"drm<-"  <- function(data,x)
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


 library(mets)

 data(sTRACE)
 sTRACE$age2 <- sTRACE$age^2
 sTRACE$age3 <- sTRACE$age^3
 names(sTRACE)

 mm <- dcut(sTRACE,~age+wmi)
 head(mm)

 mm <- dcut(sTRACE,~age+wmi,cuts=c(2,4))
 head(mm)

 mm <- dcut(sTRACE,c("age","wmi"))
 head(mm)

 mm <- dcut(sTRACE,c("age","wmi"),cuts=c(2,4))
 head(mm)

### gx <- dcut(sTRACE$age)
### head(gx)


##'## dren(data, ~x+y) <- c("x","y")
##'## dren(data, ~x+y) <- c("x","y")
##'## drm(data) <- "*.2"
##'##
##'## "dren<-" <- function(x, var, value) { browser(); names(x)[var] <- value; return(x) }
##'## dren <- function(x,var,value) {  } 

 ## Removes all cuts variables with these names wildcards
 mm <- drm(mm,c("*.2","*.4"))
 head(mm)

 ## wildcards, for age, age2, age4 and wmi 
 mm <- dcut(mm,c("a*","?m*"))
 head(mm)

 mm <- drm(mm,c("*.2","*.4"))
 head(mm)

 dcut(mm) <- c("a*","?m*")
 drm(mm) <- c("*.2","*.4")


