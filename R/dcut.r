##' Cutting, rm (removing), rename for data frames 
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
##' drename(mm, ~age+wmi) <- c("wmi","age")
##' head(mm)
##' drename(mm, ~wmi+age) <- c("age","wmi")
##' head(mm)
##' mm1 <-  drename(mm, ~age+wmi,~wmi+age)
##' head(mm1)
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


if (is.null(data)) return(gx) else return(data)
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
    if (inherits(var,"formula")) {
       var <- all.vars(var)
    }


    if (inherits(value,"formula")) {
       value <- all.vars(value)
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

 if (missing(value)) xnames <- tolower(varxnames) else if (is.character(value)) {
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

### library(mets)
###
### data(sTRACE)
### sTRACE$age2 <- sTRACE$age^2
### sTRACE$age3 <- sTRACE$age^3
### names(sTRACE)
###
### mm <- sTRACE
### mm <- dcut(sTRACE,~age+wmi)
### head(mm)
###
### ## Removes all cuts variables with these names wildcards
### mm <- drm(mm,c("*.2","*.4"))
### head(mm)
###
### ## wildcards, for age, age2, age4 and wmi 
### mm <- dcut(mm,c("a*","?m*"))
### head(mm)
###
### mm <- drm(mm,c("*.2","*.4"))
### head(mm)
###
### dcut(mm) <- c("a*","?m*")
### head(mm)
### drm(mm) <- c("*.2","*.4")
### head(mm)


###library(mets)
### data(sTRACE)
### sTRACE$age2 <- sTRACE$age^2
### sTRACE$age3 <- sTRACE$age^3
### names(sTRACE)
###
### mm <- sTRACE


##' dren(data, ~x+y) <- c("x","y")
###head(drename(mm, ~age+wmi,~ageny+wminy))
###mm$Age <- mm$age
###mm$Wmi <- mm$wmi
###head(drename(mm,~Age+Wmi))
###
###head(drename(mm, ~age+wmi,~ageny+wminy))
###head( drename(mm, ~age+wmi,c("wmi","age")))
###head( drename(mm, c("age","wmi"),c("wmi","age")))
###head( dren(mm, ~age+wmi,c("wminy","ageny")))
###
###drename(mm, ~age+wmi) <- c("wmi","age")
###drename(mm, ~age+wmi) <- c("wmi","age")

