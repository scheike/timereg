##' Cutting, sorting, rm (removing), rename for data frames
##'
##' Cut variables, if breaks are given these are used, otherwise cuts into 
##' using group size given by probs, or equispace groups on range. Default 
##' is equally sized groups if possible
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param probs groups defined from quantiles
##' @param breaks  number of breaks, for variables or vector of break points,
##' @param equi for equi-spaced breaks  
##' @param regex for regular expressions.
##' @param sep seperator for naming of cut names.
##' @param na.rm to remove NA for grouping variables.
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
##' mm <- dcut(sTRACE,~age+wmi,breaks=c(2,4))
##' head(mm)
##'
##' mm <- dcut(sTRACE,c("age","wmi"))
##' head(mm)
##'
##' mm <- dcut(sTRACE,~.)
##' head(mm)
##'
##' mm <- dcut(sTRACE,c("age","wmi"),breaks=c(2,4))
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
##' ############################
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
##' drename(mm) <- "A*"
##' head(mm)
##'
##' dd <- data.frame(A_1=1:2,B_1=1:2)
##' funn <- function(x) gsub("_",".",x)
##' drename(dd) <- ~.
##' drename(dd,fun=funn) <- ~.
##' names(dd)
##' @aliases dcut dcut<- dunique drm drm<- dnames dnames<- drename drename<- dkeep dkeep<- ddrop ddrop<- 
##' @export
dcut <- function(data,x,breaks=4,probs=NULL,equi=FALSE,regex=mets.options()$regex,sep=NULL,na.rm=TRUE,...)
{# {{{
    if (is.vector(data)) {# {{{
	if (is.list(breaks)) breaks <- unlist(breaks)

        if (length(breaks)==1) { 
             if (!is.null(probs))
	     {
                breaks <- quantile(data, probs, na.rm=na.rm, ...)
	     } else {
	   	if (!equi) { 
			probs <- seq(0, 1, length.out = breaks + 1)
			breaks <- quantile(data, probs,na.rm=na.rm, ...)
		} 
		if (equi) { 
			rr <- range(data,na.rm=na.rm)
			breaks <-  seq(rr[1],rr[2],length.out=breaks+1)
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
    }# }}}

if (is.data.frame(data)) {
  usernames <- FALSE

  if (inherits(x,"formula")) {
     vars <-mets:::procform(x,data=data,...)

      usernames <- FALSE
	 if (!is.null(vars$response)) usernames<-TRUE
	 if (usernames) {
            newnames <- vars$response
	    if (length(vars$response)!=length(vars$predictor)) { 
	    warning("length of new names not consistent with length of cut variables, uses default naming\n"); 
	   usernames <- FALSE
	   }
	 }
 }
 
 if (is.null(sep)) sep <- "."

if (missing(x)) x<- ~.

 if (inherits(x,"formula")) {
     x <- mets:::procform(x,data=data)$predictor
###     x <- all.vars(x)
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

  if (ll==1 & !is.list(breaks) & length(breaks)>1) breaks <- list(breaks)

  break.points <- FALSE
  if (is.list(breaks)) {
     break.points <- TRUE
     if (length(x)!=length(breaks) & length(breaks)!=1) 
	     warning("length of variables not consistent with list of breaks"); 
     if (length(breaks)!=ll) breaks <- rep(breaks[[1]],ll)
  }


  if (!break.points) {
     if (length(x)!=length(breaks) & length(breaks)!=1) 
	     warning("length of variables not consistent with breaks"); 
     if (length(breaks)!=ll) breaks<- rep(breaks[1],ll)
  }


for (k in 1:ll)
{
  xx <- x[[k]]
  if (is.numeric(xx)) {

      if (!is.list(breaks))
      {
          if (!is.null(probs))
	  {
                bb <- quantile(xx, probs,na.rm=na.rm, ...)
	  } else {
	   	if (!equi) { 
			probs <- seq(0, 1, length.out = breaks[k] + 1)
			bb <- quantile(xx, probs, na.rm=na.rm,...)
		} 
		if (equi) { 
			rr <- range(xx,na.rm=na.rm)
			bb <-  seq(rr[1],rr[2],length.out=breaks[k]+1)
		}
	     }
          name<-paste(xnames[k],breaks[k],sep=sep)
      } else { bb <- breaks[[1]]; name<-paste(xnames[k],breaks[[k]][1],sep=sep) }

      if (usernames) name <- newnames[k]

      print(breaks)


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
##' dfactor(mena)  <- ~twinnum
##' dnumeric(mena) <- ~twinnum.f
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
##' ### level 1 of zyg as baseline for new variable
##' drelevel(mena,ref=1) <- ~zyg
##' drelevel(mena,ref=c("DZ","[1973,1975]")) <- ~ zyg+cohort
##' ### level 2 of zyg and cohort as baseline for new variables
##' drelevel(mena,ref=2) <- ~ zyg+cohort
##' dlevels(mena)
##' 
##' @aliases dlevels drelevel drelevel<- dfactor dfactor<- dnumeric dnumeric<-
##' @export
drelevel <- function(data,x,ref=NULL,regex=mets.options()$regex,sep=NULL,...)
{# {{{

 if (missing(x) & is.data.frame(data))  stop("specify factor to relevel for data frame\n")
 if (is.null(ref)) stop("specify baseline-reference level \n")

 if (is.null(sep))  sep <- "."

 if (is.vector(data) | inherits(data,"factor")) {

      if (is.vector(data)) data <- factor(data)
      if (is.numeric(ref)) ref <-  levels(data)[ref]
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

  if (is.numeric(ref[k])) refk <- levels(xx)[ref[k]] else refk <- ref[k]
  data[,name] <- relevel(xx,ref=refk)
}

return(data)
}

}# }}}

##' @export
"drelevel<-" <- function(data,...,value) drelevel(data,value,...)

##' @export
dlevels <- function(data,x,regex=mets.options()$regex,max.levels=20,cols=FALSE,...)
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

  if (cols==TRUE) {
        antfactor <- 0
        namesfac <- c()
        nlev     <- c()
        lll <- list()
	maxl <- 0
	m <- 0
  }
                    
for (k in 1:ll)
{
  xx <- x[[k]]
  if (is.factor(xx))  {
	  cat(paste(xnames[k],":",sep=" #levels="));
	  nxx <- nlevels(xx) 
	  cat(paste(nxx,"\n")); 
     if (is.null(max.levels) || ((!is.null(max.levels)) & (nxx<max.levels))) {
        if (cols==FALSE)   print(base::levels(xx)) 
        if (cols==TRUE)  { maxl <- ifelse(base::nlevels(xx) > maxl, base::nlevels(xx), maxl); 
	                   antfactor <- antfactor+1; 
	                   namesfac <- c(namesfac,xnames[k])
	                   nlev <- c(nlev,base::nlevels(xx))
			   m <- m+1
                    	   lll[[m]] <- base::levels(xx)
                         }

     }
     if (cols==FALSE)   cat("-----------------------------------------\n")
   }
}

if (cols==TRUE) { 
	mout <- matrix("",maxl,antfactor)
	for (k in 1:antfactor) {
		mout[1:nlev[k],k] <- lll[[k]]
	}
	colnames(mout) <- namesfac
	rownames(mout) <- rep("  ",nrow(mout))
	prmatrix(mout,quote=FALSE)
}

}

}# }}}


##' @export
drm <- function(data,x,regex=mets.options()$regex)
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

	if (!is.null(var))    {
		var <- procform(var,data=data,return.list=TRUE,...)
	        varargs <-   c(!is.null(var$predictor),!is.null(var$response))
	} else varargs <- c(0,0)
	if (!is.null(value)) {
		value <- procform(value,data=data,return.list=TRUE,...)
	        valueargs <- c(!is.null(value$predictor),!is.null(value$response))
	} else valueargs <- c(0,0)

	vargs <- 1*varargs+1*valueargs

	if (sum(varargs)+sum(valueargs)>=3) 
		stop("arguments specified multiple times \n")

	if (sum(varargs)==2)   {
		value <- var$response;
		var <- var$predictor;   
	} else if (sum(valueargs)==2) {
		var <-   value$predictor; 
		value <- value$response;
	} else if (sum(vargs)==2) {
		value   <- c(value$predictor,value$response)
		var     <- c(var$predictor,var$response)
	} else if (sum(vargs)==1) {## only one set of variables , so  use fun on these vars 
		var <- c(value$predictor,value$response,var$response,var$predictor)
	        value <- do.call(fun,list(var))
        } else { ## nothing given so 
		var <- colnames(data); varpos <- seq(ncol(data))
	        value <- do.call(fun,list(var))
        }

        varpos <- match(var,colnames(data))

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
dkeep <- function(data,var,keep=TRUE,regex=mets.options()$regex)
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


##' @export
dfactor <- function(data,x,regex=mets.options()$regex,sep=NULL,levels,labels,...)
{# {{{

 if (is.null(sep))  sep <- ".f"

 if (is.vector(data)) {
	 if (!is.factor(data)) {
		 args <- list(data)
		 if (!missing(levels)) args <- c(args,list(levels=levels))
		 if (!missing(labels)) args <- c(args,list(labels=labels))
	         gx <- do.call(factor,args)
	 } else gx <- data
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

 if (!missing(levels) & !is.list(levels) )   levels <- list(levels)
 if (!missing(labels) & !is.list(labels) )   labels <- list(labels)

### if (!missing(levels) ) print(levels)
### if (!missing(labels) ) print(labels)

  misslabel <- TRUE
  if (!missing(labels))  {
	  misslabel <- FALSE
  if ((length(x)!=length(labels))) {
	  warning("length of label list not consistent with variables")
	  labels <- rep(list(labels[[1]]),ll)
###	  print(labels)
  }
  }


  misslevel <- TRUE
  if (!missing(levels))  { 
	  misslevel <- FALSE
  if ((length(x)!=length(levels))) {
	  warning("length of levels list not consistent with variables")
	  levels <- rep(list(levels[[1]]),ll)
  }
  }
 




for (k in 1:ll)
{
  xx <- x[[k]]
  name<- paste(xnames[k],sep,sep="")
  if (!is.factor(xx) || all==TRUE) { 
           args <- list(xx)
           if (!misslevel) args <- c(args,list(levels=levels[[k]]))
	   if (!misslabel) args <- c(args,list(labels=labels[[k]]))
	   gx <- do.call(factor,args)
           data[,name] <- gx  
   }
}

return(data)
}

}# }}}

##' @export
"dfactor<-" <- function(data,k=1,combine=TRUE,...,value) {
    dfactor(data,value,...)
}

##' @export
dnumeric <- function(data,x,regex=mets.options()$regex,sep=NULL,all=FALSE,...)
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
          data[,name] <- gx 
  } 
}

return(data)
}

}# }}}

##' @export
"dnumeric<-" <- function(data,...,value) {
    dnumeric(data,value,...)
}



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
    return(as.matrix(dval))
}

##' @export
"dlag<-" <- function(data,k=1,combine=TRUE,...,value) {
    dlag(data,value,k=k,combine=combine,...)
}

