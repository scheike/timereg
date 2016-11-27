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
##' ## all to lower 
##' head(drename(mm))
##'
##' ## A* to lower 
##' mm2 <-  drename(mm,c("A*","W*"))
##' head(mm2)
##' drename(mm,"A*") <- ~.
##' head(mm)
##' drename(mm) <- ~.
##' head(mm)
##' @export
dcut <- function(data,x,cuts=4,breaks=NULL,sep=NULL,...)
{# {{{
 if (is.null(sep) & is.null(breaks))    sep <- "."
 if (is.null(sep) & (!is.null(breaks))) sep <- "b"

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
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
     name<-paste(xnames[k],cuts[k],sep=sep)
     bb <- quantile(xx,probs)
     } else { bb <- breaks ; name<-paste(xnames[k],breaks[2],sep=sep) }

  if (sum(duplicated(bb))==0)
   data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE,...)
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
     if (value[1]==".") value <- tolower(varxnames)
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

### print(xnames)
### print(varxnames)

  if (length(xnames)!= length(varxnames)) stop("length of old and new variables must mach")

  names(data)[varnnames] <- xnames; 
  return(data) 
} # }}}

##' @export
"drename<-" <- function(x, var, value) 
{ # {{{
  drename(x,var,value)
}# }}}


##' @export
dkeep <- function(data,var,keep=TRUE) 
{  # {{{

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


if (keep) data <- data[,varnnames] else data <- data[,-1*varnnames]
###  names(data)[varnnames] <- xnames; 
  return(data) 
} # }}}

##' @export
"dkeep<-" <- function(x,...,value) 
{ # {{{
  dkeep(x,value,...)
}# }}}

##' @export
ddrop <- function(data,var,keep=FALSE) 
{  # {{{
	dkeep(data,var,keep=FALSE)
} # }}}

##' @export
"ddrop<-" <- function(x,...,value) 
{ # {{{
  dkeep(x,value,keep=FALSE,...)
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
"dsort<-" <- function(data,...,value) 
{ # {{{
  dsort(data,value,...)
}# }}}


### måske og x+y ~ group1+group2 to get correlations against groups defined by group1*group2


##' summary, tables, and correlations for data frames 
##' 
##' summary, tables, and correlations for data frames 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param g possible group variable
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike 
##' @examples
##' data("sTRACE",package="timereg")
##' dt<- sTRACE
##' dt$time2 <- dt$time^2
##' dt$wmi2 <- dt$wmi^2
##' head(dt)
##'
##' 
##' dcor(dt,~time+wmi)
##' dcor(dt,~time+wmi,~vf+chf)
##' dcor(dt,time+wmi~vf+chf)
##' 
##' dcor(dt,c("time*","wmi*"),~vf+chf)
##' 
##' @export
dcor <- function(data,x,g,...)
{# {{{

### x+y ~ group1+group2 to get correlations against groups defined by group1*group2
 rhs <- NULL
 if (inherits(x,"formula")) {
         rhs  <-  all.vars(update(x, 0~.))
	 lhs <-  all.vars(update(x, .~0))
	 if (lhs[1]!=".") { g <- rhs; x <- lhs;}
 } 

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
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

 if (!missing(g))
 if (inherits(g,"formula")) {
     g <- all.vars(g)
     if (g[1]==".") g <- names(data) 
     gnames <- g
     formular <- 1
  } else if  (is.character(g)) {
     gnames <- g
     xxx<-c()
     for (xx in gnames)
     {
        n <- grep(glob2rx(xx),names(data))
        xxx <- c(xxx,names(data)[n])
     }
     gnames <- xxx[!duplicated(xxx)]
  }


 if (!missing(g)) return(by(data[,xnames],data[,gnames],cor,...))
 if (missing(g)) return(cor(data[,xnames],...))

}# }}}

##' @export
dsummary <- function(data,x,g,...)
{# {{{

### x+y ~ group1+group2 to get correlations against groups defined by group1*group2
 rhs <- NULL
 if (inherits(x,"formula")) {
         rhs  <-  all.vars(update(x, 0~.))
	 lhs <-  all.vars(update(x, .~0))
	 if (lhs[1]!=".") { g <- rhs; x <- lhs;}
 } 


 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
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

 if (!missing(g))
 if (inherits(g,"formula")) {
     g <- all.vars(g)
     if (g[1]==".") g <- names(data) 
     gnames <- g
     formular <- 1
  } else if  (is.character(g)) {
     gnames <- g
     xxx<-c()
     for (xx in gnames)
     {
        n <- grep(glob2rx(xx),names(data))
        xxx <- c(xxx,names(data)[n])
     }
     gnames <- xxx[!duplicated(xxx)]
  }

 if (!missing(g)) return(by(data[,xnames],data[,gnames],summary,...))
 if (missing(g)) return(summary(data[,xnames],...))

}# }}}


##' tables for data frames 
##' 
##' tables for data frames 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param g possible group variable
##' @param ... Optional additional arguments
##' @author Klaus K. Holst and Thomas Scheike 
##' @examples
##' data("sTRACE",package="timereg")
##' dt<- sTRACE
##'
##' dtable(dt,~status+vf,~chf+diabetes)
##'
##' dtable(dt,status+vf~chf+diabetes)
##' 
##' dtable(dt,c("*f*","status"),~diabetes)
##' dtable(dt,c("*f*","status"),~diabetes,all2by2=FALSE)
##' 
##' @export
dtable<- function(data,x,g,all2by2=TRUE,...)
{# {{{

 rhs <- NULL
 if (inherits(x,"formula")) {
         rhs  <-  all.vars(update(x, 0~.))
	 lhs <-  all.vars(update(x, .~0))
	 if (lhs[1]!=".") { g <- rhs; x <- lhs;}
 } 

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
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

 if (!missing(g))
 if (inherits(g,"formula")) {
     g <- all.vars(g)
     if (g[1]==".") g <- names(data) 
     gnames <- g
     formular <- 1
  } else if  (is.character(g)) {
     gnames <- g
     xxx<-c()
     for (xx in gnames)
     {
        n <- grep(glob2rx(xx),names(data))
        xxx <- c(xxx,names(data)[n])
     }
     gnames <- xxx[!duplicated(xxx)]
  }

 if (all2by2==TRUE) {
 ### all 2 by 2 tables from xnames over g 
 ll<-list()
 k<-1
 nn <- length(xnames)
 for (i in seq(1,(nn-1)))
 for (j in seq((i+1),nn)) {
	 x1<-xnames[i]
	 x2<-xnames[j]
        llk<-table(data[,x1],data[,x2],...)
        ll[[k]]<- list(name=paste(x1,"x",x2,sep=""),table=llk)
	k<-k+1
 }
 } else { ## one big table 
     ll<-table(data[,xnames],...)
 }

 return(ll)

### if (!missing(g)) return(by(data[,xnames],data[,gnames],table,...)) else return(by(data[,xnames],data[,gnames],table,...))

}# }}}



##' @export
dhead <- function(data,x,...)
{# {{{

 if (inherits(x,"formula")) {
     x <- all.vars(x)
     if (x[1]==".") x <- names(data) 
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


return(head(data[,xnames],...))

}# }}}


##' @export
dreshape <- function(data,...)
{# {{{
    fast.reshape(data,...)
}# }}}



