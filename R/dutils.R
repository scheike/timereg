##' Cutting, sorting, rm (removing), rename for data frames
##'
##' Cut variables, if breaks are given these are used, otherwise cuts into 
##' using group size given by probs, or equispace groups on range. Default 
##' is equally sized groups if possible
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param y name of variable, or fomula, or names of variables on data frame.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param probs groups defined from quantiles
##' @param breaks  number of breaks, for variables or vector of break points,
##' @param equi for equi-spaced breaks  
##' @param regex for regular expressions.
##' @param sep seperator for naming of cut names.
##' @param na.rm to remove NA for grouping variables.
##' @param labels to use for cut groups 
##' @param all to do all variables, even when breaks are not unique 
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
##' mm <- dcut(sTRACE,catage4+wmi4~age+wmi)
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
##' dcut(mm) <- ageg1+wmig1~age+wmi
##' head(mm)
##'
##' ############################
##' ## renaming
##' ############################
##'
##' head(mm)
##' drename(mm, ~Age+Wmi) <- c("wmi","age")
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
dcut <- function(data,y=NULL,x=NULL,breaks=4,probs=NULL,equi=FALSE,regex=mets.options()$regex,sep=NULL,na.rm=TRUE,labels=NULL,all=FALSE,...)
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
             gx <- cut(data, breaks = breaks, include.lowest = TRUE,labels=labels,...)
        else {
	      wd <- which(duplicated(breaks))
              mb <- min(diff(breaks[-wd]))
	      breaks[wd] <- breaks[wd] +  (mb/2)*seq(length(wd))/length(wd)
              gx  <- cut(data,breaks=breaks,include.lowest=TRUE,labels=labels,...)
              warning(paste("breaks duplicated"))
        }
        return(gx)
    }# }}}

if (is.data.frame(data)) {# {{{

 if (is.null(sep)) sep <- "cat."

 usernames <- FALSE# {{{

     vars <-mets::procform3(y,x,data=data,regex=regex,...)
     x <-  xnames <- vars$x

     if (!is.null(vars$y)) {
         usernames<-TRUE
         newnames <- vars$y
	 if (length(vars$y)!=length(vars$x)) { 
	   warning("length of new names not consistent with length of cut variables, uses default naming\n"); 
	   usernames <- FALSE
      }
     }
# }}}

  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) { if (length(x)==1 && is.character(x)) x <- data[,x]; x })

  if (!is.list(x)) x <- list(x)
  ll <- length(x)
  if (ll==1 & !is.list(breaks) & length(breaks)>1) breaks <- list(breaks)

  break.points <- FALSE
  if (is.list(breaks)) {
     break.points <- TRUE
     if (length(x)!=length(breaks) & length(breaks)!=1) 
	     warning("length of variables not consistent with list of breaks"); 
     if (length(breaks)!=ll) breaks <- rep(list(breaks[[1]]),ll)
  }

  if (!break.points) {
     if (length(x)!=length(breaks) & length(breaks)!=1) 
	     warning("length of variables not consistent with breaks"); 
     if (length(breaks)!=ll) breaks<- rep(breaks[1],ll)
  }

  if (ll==1 & !is.list(labels)) labels <- list(labels)
  if (!is.list(labels)) labels <- list(labels); 
  if (length(labels)!=ll ) labels <- rep(list(labels[[1]]),ll)
  if (!is.list(labels)) stop("labels should be given as list"); 


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
      } else { bb <- breaks[[k]]; name<-paste(xnames[k],breaks[[k]][1],sep=sep) }

      if (usernames) name <- newnames[k]

      if (sum(duplicated(bb))==0)
	     data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE,labels=labels[[k]],...)
      else { 
	   if (all==TRUE) {
	      wd <- which(duplicated(bb))
              mb <- min(diff(bb[-wd]))
	      bb[wd] <- bb[wd] +  (mb/2)*seq(length(wd))/length(wd)
          data[,name] <- cut(xx,breaks=bb,include.lowest=TRUE,labels=labels[[k]],...)
          warning(paste("breaks duplicated for=",xnames[k]))
	   }
      }
   }
}

return(data)
}# }}}

}# }}}

##' @export
"dcut<-" <- function(data,...,value) dcut(data,y=value,...)

##' relev levels for data frames
##'
##' levels shows levels for variables in data frame, relevel relevels a factor in data.frame 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param y name of variable, or fomula, or names of variables on data frame.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param ref new reference variable 
##' @param newlevels to combine levels of factor in data frame
##' @param regex for regular expressions.
##' @param sep seperator for naming of cut names.
##' @param overwrite to overwrite variable 
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
##' mena2 <- drelevel(mena,"cohort",ref="(1980,1982]")
##' mena2 <- drelevel(mena,~cohort,ref="(1980,1982]")
##' mena2 <- drelevel(mena,cohortII~cohort,ref="(1980,1982]")
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
##' drelevel(mena,ref=c("DZ","[1973,1975]")) <- zygdz+cohort.early~ zyg+cohort
##' ### level 2 of zyg and cohort as baseline for new variables
##' drelevel(mena,ref=2) <- ~ zyg+cohort
##' dlevels(mena)
##' 
##' ##################### combining factor levels with newlevels argument
##' 
##' dcut(mena,labels=c("I","II","III","IV")) <- cat4~agemena
##' dlevels(drelevel(mena,~cat4,newlevels=1:3))
##' dlevels(drelevel(mena,ncat4~cat4,newlevels=3:2))
##' drelevel(mena,newlevels=3:2) <- ncat4~cat4
##' dlevels(mena)
##' 
##' dlevels(drelevel(mena,nca4~cat4,newlevels=list(c(1,4),2:3)))
##' 
##' drelevel(mena,newlevels=list(c(1,4),2:3)) <- nca4..2 ~ cat4
##' dlevels(mena)
##' 
##' drelevel(mena,newlevels=list("I-III"=c("I","II","III"),"IV"="IV")) <- nca4..3 ~ cat4
##' dlevels(mena)
##' 
##' drelevel(mena,newlevels=list("I-III"=c("I","II","III"))) <- nca4..4 ~ cat4
##' dlevels(mena)
##' 
##' drelevel(mena,newlevels=list(group1=c("I","II","III"))) <- nca4..5 ~ cat4
##' dlevels(mena)
##' 
##' drelevel(mena,newlevels=list(g1=c("I","II","III"),g2="IV")) <- nca4..6 ~ cat4
##' dlevels(mena)
##' 
##' @aliases dlevels dlevel dlev drelevel drelev dlev<- dlevel<- drelev<- drelevel<- dfactor dfactor<- dnumeric dnumeric<-
##' @export
drelevel <- function(data,y=NULL,x=NULL,ref=NULL,newlevels=NULL,regex=mets.options()$regex,sep=NULL,overwrite=FALSE,...)
{# {{{

 if (is.null(ref) && is.null(newlevels)) stop("specify baseline-reference level or new levels \n")

 if (!is.null(ref) & !is.null(newlevels)) { 
	 warning("can only either change ref or combine old levels, will change ref")
	 newlevels <- NULL
 }

 if (is.null(sep))  sep <- "."

 if (is.vector(data) | inherits(data,"factor")) {# {{{

      if (is.vector(data)) data <- factor(data)
      if (!is.null(ref)) {
	      if (is.numeric(ref)) ref <-  levels(data)[ref]
              gx <- relevel(data,ref=ref,...)
      return(gx)
      }
      if (!is.null(newlevels)) {
	      pnewlevels <- levlev(data,newlevels)
	      levels(data,...) <- pnewlevels
	      return(data)
      }
 } # }}}

if (is.data.frame(data)) {# {{{

 usernames <- FALSE# {{{

     vars <-mets::procform3(y,x,data=data,regex=regex,...)
     x <-  xnames <- vars$x

     if (!is.null(vars$y)) {
         usernames<-TRUE
         newnames <- vars$y
	 if (length(vars$y)!=length(vars$x)) { 
	   warning("length of new names not consistent with length of cut variables, uses default naming\n"); 
	   usernames <- FALSE
      }
     }
# }}}

  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) { if (length(x)==1 && is.character(x)) x <- data[,x]; x })

  if (!is.list(x)) x <- list(x)
  ll <- length(x)

  if (!is.null(ref)) 
  {  if (ll>1 & length(ref)==1) ref <- rep(ref,ll)
      if (length(x)!=length(ref)) stop("length of baseline reference 'ref' not consistent with variables")
  }

  if (!is.null(newlevels)) {
      if (ll==1 & !is.list(newlevels)) newlevels <- list(newlevels)
      if (is.list(newlevels) && !is.list(newlevels[[1]])) newlevels <- list(newlevels)
      if (length(x)!=length(newlevels)) 
          warning("length of variables not consistent with list of breaks"); 
      if (length(newlevels)!=ll) newlevels <- rep(list(newlevels[[1]]),ll)
  }


for (k in 1:ll)
{
  xx <- x[[k]]
  if (!is.factor(xx)) xx <- factor(xx)

  if (usernames) name <- newnames[k]
  if (!is.null(ref)) {
  name<- paste(xnames[k],ref[k],sep=sep)
  if (usernames) name <- newnames[k]
  if (overwrite) name<-xnames[k]
      if (is.numeric(ref[k])) refk <-  levels(xx)[ref[k]] else refk <- ref[k]
      gx <- relevel(xx,ref=refk,...)
      data[,name] <- gx
  }
  if (!is.null(newlevels)) {
  name<- paste(xnames[k],newlevels[[k]][1],sep=sep)
  if (usernames) name <- newnames[k]
  if (overwrite) name<-xnames[k]
      pnewlevels <- levlev(xx,newlevels[[k]])
      levels(xx,...) <- pnewlevels
      data[,name] <- xx
  }
}

return(data)
}# }}}

}# }}}

##' @export
"drelev<-" <- function(data,x=NULL,...,value) drelevel(data,y=value,x=x,...)

##' @export
drelev <- function(data,y=NULL,x=NULL,...) drelevel(data,y=y,x=x,...)

##' @export
"drelevel<-" <- function(data,x=NULL,...,value) drelevel(data,y=value,x=x,...)

tsglob2rx <- function(x) {
 glob2rx(gsub("\\+","\\\\+",x))
}

levlev <- function(fac,ref,regex=FALSE)
{# {{{
if (!is.list(ref)) ref <- list(ref)
lf <- levels(fac)
lfr <- lf
listnames <- names(ref)

newreflist <- list()
for (k in 1:length(ref)) {
	if (!is.null(listnames)) nn <- listnames[k] else nn <- NULL
	ln <- length(ref[[k]])
	if (is.numeric(ref[[k]])) refs <- lf[ref[[k]]] else refs <- ref[[k]]
	xxx <- c()
	for (xx in refs)
	{
	   if (!regex) xx <- tsglob2rx(xx)
           n <- grep(xx,lf)
	   xxx <- c(xxx,lf[n])
	}
	xxx<- xxx[!duplicated(xxx)]
	refs <- xxx

	ln <- length(refs)
	if (is.null(nn) || nn=="") 
	{
	if (length(refs)>1) nn <- paste(refs[1],refs[ln],sep="-") else nn <- refs[1]
	}
	newreflist <- c(newreflist,setNames(list(refs),nn))
	mm <- match(refs,lfr)

	lfr <- lfr[-mm]
}


if (length(lfr)>=1)
{
     for (k in 1:length(lfr)) 
     {
		nn <- paste(lfr[k])
###		newreflist <- c(newreflist,list(nn1=lfr[k]))
	newreflist <- c(newreflist,setNames(list(lfr[k]),nn))
     }
}

return(newreflist)
}# }}}

##' @export
dlevels <- function(data,y=NULL,x=NULL,regex=mets.options()$regex,max.levels=20,cols=FALSE,...)
{# {{{

 if (is.factor(data)) {
	 return(base::levels(data,...)) 
 } 

 if (is.data.frame(data)) {
   usernames <- FALSE# {{{

     vars <-mets::procform3(y,x,data=data,regex=regex,...)
     x <-  xnames <- vars$x

     if (!is.null(vars$y)) {
         usernames<-TRUE
         newnames <- vars$y
	 if (length(vars$y)!=length(vars$x)) { 
	   warning("length of new names not consistent with length of cut variables, uses default naming\n"); 
	   usernames <- FALSE
      }
     }
# }}}

  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) { if (length(x)==1 && is.character(x)) x <- data[,x]; x })
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
dlevel <- function(data,y=NULL,x=NULL,...) dlevels(data,y=y,x=x,...)

##' @export
dlev <- function(data,y=NULL,x=NULL,...) dlevels(data,y=y,x=x,...)

##' @export
drename <- function(data,y=NULL,x=NULL,fun=base::tolower,...)
{  # {{{

    vars <-mets::procform3(y,x,data=data,...)
    x <-  xnames <- vars$x

    if (!is.null(vars$y)) {
         newnames <- vars$y
	 if (length(vars$y)!=length(vars$x)) { 
	   stop("length of new names not consistent with length of cut variables, uses default naming\n"); 
      }
    } else { ## if newnames not given then use fun 
       newnames <- do.call(fun,list(x))
    }

    varpos <- match(x,colnames(data))

    if (length(varpos)!= length(newnames)) stop("length of old and new variables must match")
    colnames(data)[varpos] <- newnames
    return(data)
} # }}}

##' @export
"drename<-" <- function(data,x=NULL,...,value) drename(data,y=value,x=x,...)

##' @export
dfactor <- function(data,y=NULL,x=NULL,regex=mets.options()$regex,sep=NULL,usernames=NULL,levels,labels,...)
{# {{{

 if (is.null(sep))  sep <- ".f"

 if (!is.data.frame(data)) {
       if (!is.factor(data) || !missing(levels) || !missing(labels)) {
		 args <- list(data)
		 if (!missing(levels))  {
			 if (is.numeric(levels) & is.factor(data)) 
				 levels <-  levels(data)[levels]
			 args <- c(args,list(levels=levels,...))
		 }
		 if (!missing(labels)) {
			 args <- c(args,list(labels=labels,...))
		 }
	         gx <- do.call(factor,args)
      return(gx)
      } 
 }

if (is.data.frame(data)) {

   usernames <- FALSE# {{{

   vars <-mets::procform3(y,x,data=data,regex=regex,...)
   x <-  xnames <- vars$x

   if (!is.null(vars$y)) {
	usernames<-TRUE
	newnames <- vars$y
	if (length(vars$y)!=length(vars$x)) { 
		usernames <- FALSE
	}
   }
# }}}

  if (is.character(x) && length(x)<=ncol(data)) {
	  x <- lapply(xnames,function(z) data[,z])
  }
  dots <- list()
  args <- lapply(dots, function(x) { if (length(x)==1 && is.character(x)) x <- data[,x]; x })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)

 if (!missing(levels)) if (!is.list(levels))   levels <- list(levels)
 if (!missing(labels)) if (!is.list(labels) )   labels <- list(labels)

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

  if (usernames) name <- newnames[k]

  if (!is.factor(xx) || !missing(levels) || !missing(labels)) {
           args <- list(xx)
           if (!misslevel) {
		 if (is.numeric(levels[[k]]) & is.factor(xx)) 
				 llevels <-  levels(xx)[levels[[k]]]
			 else llevels <- levels[[k]]
		   args <- c(args,list(levels=llevels,...))
	   }
	   if (!misslabel) args <- c(args,list(labels=labels[[k]],...))
	   gx <- do.call(factor,args)
           data[,name] <- gx  
 }
}

return(data)
}

}# }}}

##' @export
"dfactor<-" <- function(data,x=NULL,...,value) dfactor(data,y=value,x=x,...)

#####' @export
###dfactor <- function(data,y=NULL,x=NULL,regex=mets.options()$regex,sep=NULL,usernames=NULL,levels,labels,...)
###{# {{{
###
### if (is.null(sep))  sep <- ".f"
###
### if (is.vector(data)) {
###	 if (!is.factor(data)) {
###		 args <- list(data)
###		 if (!missing(levels)) args <- c(args,list(levels=levels,...))
###		 if (!missing(labels)) args <- c(args,list(labels=labels,...))
###	         gx <- do.call(factor,args)
###	 } else gx <- data
###      return(gx)
### } 
###
###if (is.data.frame(data)) {
###
###   usernames <- FALSE# {{{
###
###vars <-mets::procform3(y,x,data=data,regex=regex,...)
###x <-  xnames <- vars$x
###
###if (!is.null(vars$y)) {
###	usernames<-TRUE
###	newnames <- vars$y
###	if (length(vars$y)!=length(vars$x)) { 
###		usernames <- FALSE
###	}
###}
#### }}}
###
###  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
###  dots <- list()
###  args <- lapply(dots, function(x) { if (length(x)==1 && is.character(x)) x <- data[,x]; x })
###  if (!is.list(x)) x <- list(x)
###  ll <- length(x)
###
###
### if (!missing(levels)) if (!is.list(levels))   levels <- list(levels)
### if (!missing(labels)) if (!is.list(labels) )   labels <- list(labels)
###
###### if (!missing(levels) ) print(levels)
###### if (!missing(labels) ) print(labels)
###
###  misslabel <- TRUE
###  if (!missing(labels))  {
###	  misslabel <- FALSE
###  if ((length(x)!=length(labels))) {
###	  warning("length of label list not consistent with variables")
###	  labels <- rep(list(labels[[1]]),ll)
######	  print(labels)
###  }
###  }
###
###
###  misslevel <- TRUE
###  if (!missing(levels))  { 
###	  misslevel <- FALSE
###  if ((length(x)!=length(levels))) {
###	  warning("length of levels list not consistent with variables")
###	  levels <- rep(list(levels[[1]]),ll)
###  }
###  }
### 
###
###for (k in 1:ll)
###{
###  xx <- x[[k]]
###  name<- paste(xnames[k],sep,sep="")
###
###  if (usernames) name <- newnames[k]
###
###  if (!is.factor(xx) ) { 
###           args <- list(xx)
###           if (!misslevel) args <- c(args,list(levels=levels[[k]],...))
###	   if (!misslabel) args <- c(args,list(labels=labels[[k]],...))
###	   gx <- do.call(factor,args)
###           data[,name] <- gx  
###   }
###}
###
###return(data)
###}
###
###}# }}}


##' @export
dnumeric <- function(data,y=NULL,x=NULL,regex=mets.options()$regex,sep=NULL,all=FALSE,...)
{# {{{

 if (is.null(sep))  sep <- ".n"

 if (is.factor(data)) {
      gx <- as.numeric(data) 
      return(gx)
 } 

 if (is.data.frame(data)) {

   usernames <- FALSE# {{{

     vars <-mets::procform3(y,x,data=data,regex=regex,...)
     x <-  xnames <- vars$x

     if (!is.null(vars$y)) {
         usernames<-TRUE
         newnames <- vars$y
	 if (length(vars$y)!=length(vars$x)) { 
	   warning("length of new names not consistent with length of cut variables, uses default naming\n"); 
	   usernames <- FALSE
      }
     }
# }}}

  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) { if (length(x)==1 && is.character(x)) x <- data[,x]; x })
  if (!is.list(x)) x <- list(x)
  ll <- length(x)


for (k in 1:ll)
{
  xx <- x[[k]]
  name<- paste(xnames[k],sep,sep="")
  if (usernames) name <- newnames[k]
  if (!is.numeric(xx) || all==TRUE) { 
	  gx <- as.numeric(xx,...); 
          data[,name] <- gx 
  } 
}

return(data)
}

}# }}}

##' @export
"dnumeric<-" <- function(data,x=NULL,...,value) dnumeric(data,y=value,x=x,...)

##' @export
drm <- function(data,x=NULL,regex=mets.options()$regex,...)
{# {{{
  vars <-mets::procform(x,data=data,no.match=FALSE,regex=regex,...)
  xnames <- c(vars$predictor,vars$response)

  data[,xnames] <- NULL
  return(data)
}# }}}

##' @export
"drm<-" <- function(data,...,value) drm(data,x=value,...)

##' @export
dkeep <- function(data,x=NULL,keep=TRUE,regex=mets.options()$regex,...)
{  # {{{
  vars <-mets::procform(x,data=data,no.match=FALSE,regex=regex,...)
  xnames <- c(vars$predictor,vars$response)
  nnames <- match(xnames,names(data))

  if (keep) data <- data[,nnames] else data <- data[,-1*nnames]
  return(data)
} # }}}

##' @export
"dkeep<-" <- function(data,...,value) dkeep(data,x=value,...)

##' @export
ddrop <- function(data,x=NULL,keep=FALSE,...) dkeep(data,x,keep=FALSE,...)

##' @export
"ddrop<-" <- function(data,...,value) dkeep(data,x=value,keep=FALSE,...)

##' @export
dnames <- function(data,...) drename(data,...)

##' @export
"dnames<-" <- function(data,...,value) drename(data,value=value,...)


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
dlag <- function(data,x,k=1,combine=TRUE,simplify=TRUE,names,...) {# {{{
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
        environment(x) <- environment()
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
}# }}}

##' @export
"dlag<-" <- function(data,k=1,combine=TRUE,...,value) {
    dlag(data,value,k=k,combine=combine,...)
}

##' Simple linear spline 
##' 
##' Constructs simple linear spline  on a data frame using the formula syntax of dutils
##' 
##' @param data if x is formula or names for data frame then data frame is needed.
##' @param y name of variable, or fomula, or names of variables on data frame.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param probs groups defined from quantiles
##' @param breaks  number of breaks, for variables or vector of break points,
##' @param equi for equi-spaced breaks  
##' @param regex for regular expressions.
##' @param sep seperator for naming of cut names.
##' @param na.rm to remove NA for grouping variables.
##' @param labels to use for cut groups 
##' @param all to do all variables, even when breaks are not unique 
##' @param ... Optional additional arguments
##' @author Thomas Scheike
##' @examples
##' data(TRACE)
##' TRACE <- dspline(TRACE,~wmi,breaks=c(1,1.3,1.7))
##' cca <- coxph(Surv(time,status==9)~age+vf+chf+wmi,data=TRACE)
##' cca2 <- coxph(Surv(time,status==9)~age+wmi+vf+chf+wmi.spline1+wmi.spline2+wmi.spline3,data=TRACE)
##' anova(cca,cca2)
##' 
##' nd=data.frame(age=50,vf=0,chf=0,wmi=seq(0.4,3,by=0.01))
##' nd <- dspline(nd,~wmi,breaks=c(1,1.3,1.7))
##' pl <- predict(cca2,newdata=nd)
##' plot(nd$wmi,pl,type="l")
##'
##' @export
##' @aliases dspline<-
dspline <- function(data,y=NULL,x=NULL,breaks=4,probs=NULL,equi=FALSE,regex=mets.options()$regex,sep=NULL,na.rm=TRUE,labels=NULL,all=FALSE,...)
{# {{{
    if (is.vector(data)) {# {{{
	if (is.list(breaks)) breaks <- unlist(breaks)

        if (length(breaks)==1) { 
             if (!is.null(probs))
	     {
                breaks <- quantile(data, probs, na.rm=na.rm, ...)
	        breaks <- breaks[-c(1,length(breaks))]
	     } else {
	   	if (!equi) { 
			probs <- seq(0, 1, length.out = breaks + 1)
			breaks <- quantile(data, probs,na.rm=na.rm, ...)
	                breaks <- breaks[-c(1,length(breaks))]
		} 
		if (equi) { 
			rr <- range(data,na.rm=na.rm)
			breaks <-  seq(rr[1],rr[2],length.out=breaks+1)
	                breaks <- breaks[-c(1,length(breaks))]
		}
	     }
	}

        if (sum(duplicated(breaks))==0) {
             gx <- LinSpline(data, breaks, ...)
	     attr(gx,"breaks") <- breaks
	} else {
	      wd <- which(duplicated(breaks))
              mb <- min(diff(breaks[-wd]))
	      breaks[wd] <- breaks[wd] +  (mb/2)*seq(length(wd))/length(wd)
              gx <- LinSpline(data, breaks,...)
	      attr(gx,"breaks") <- breaks
              warning(paste("breaks duplicated"))
        }
        return(gx)
    }# }}}

if (is.data.frame(data)) {# {{{

 if (is.null(sep)) sep <- "."

 usernames <- FALSE# {{{

     vars <-mets::procform3(y,x,data=data,regex=regex,...)
     x <-  xnames <- vars$x

     if (!is.null(vars$y)) {
         usernames<-TRUE
         newnames <- vars$y
	 if (length(vars$y)!=length(vars$x)) { 
	   warning("length of new names not consistent with length of cut variables, uses default naming\n"); 
	   usernames <- FALSE
      }
     }
# }}}

  if (is.character(x) && length(x)<nrow(data)) x <- lapply(xnames,function(z) data[,z])
  dots <- list()
  args <- lapply(dots, function(x) { if (length(x)==1 && is.character(x)) x <- data[,x]; x })

  if (!is.list(x)) x <- list(x)
  ll <- length(x)
  if (ll==1 & !is.list(breaks) & length(breaks)>1) breaks <- list(breaks)

  break.points <- FALSE
  if (is.list(breaks)) {
     break.points <- TRUE
     if (length(x)!=length(breaks) & length(breaks)!=1) 
	     warning("length of variables not consistent with list of breaks"); 
     if (length(breaks)!=ll) breaks <- rep(list(breaks[[1]]),ll)
  }

  if (!break.points) {
     if (length(x)!=length(breaks) & length(breaks)!=1) 
	     warning("length of variables not consistent with breaks"); 
     if (length(breaks)!=ll) breaks<- rep(breaks[1],ll)
  }

  if (ll==1 & !is.list(labels)) labels <- list(labels)
  if (!is.list(labels)) labels <- list(labels); 
  if (length(labels)!=ll ) labels <- rep(list(labels[[1]]),ll)
  if (!is.list(labels)) stop("labels should be given as list"); 


for (k in 1:ll)
{
  xx <- x[[k]]
  if (is.numeric(xx)) {

      if (!is.list(breaks))
      {
          if (!is.null(probs))
	  {
                bb <- quantile(xx, probs,na.rm=na.rm, ...)
	        bb <- bb[-c(1,length(bb))]
	  } else {
	   	if (!equi) { 
			probs <- seq(0, 1, length.out = breaks[k] + 1)
			bb <- quantile(xx, probs, na.rm=na.rm,...)
	                bb <- bb[-c(1,length(bb))]
		} 
		if (equi) { 
			rr <- range(xx,na.rm=na.rm)
			bb <-  seq(rr[1],rr[2],length.out=breaks[k]+1)
	                bb <- bb[-c(1,length(bb))]
		}
	     }
###          name<-paste(xnames[k],breaks[k],sep=sep)
          name<-xnames[k]
      } else { bb <- breaks[[k]]; 
###               name<-paste(xnames[k],breaks[[k]][1],sep=sep) 
          name<-xnames[k]
      }

      if (usernames) name <- newnames[k]

      if (sum(duplicated(bb))==0) {
	attr(data,paste(name,"spline.breaks",sep="")) <- bb
	for (i in seq_along(c(bb)))
	{
           namei <- paste(name,".spline",i,sep="")
	   data[,namei] <-  (xx-bb[i])*(xx>bb[i])
	}
      }
      else { 
	   if (all==TRUE) {
	      wd <- which(duplicated(bb))
              mb <- min(diff(bb[-wd]))
	      bb[wd] <- bb[wd] +  (mb/2)*seq(length(wd))/length(wd)
	      attr(data,paste(name,"spline.breaks",sep="")) <- bb
              for (i in seq_along(c(breaks)))
              {
		     namei <- paste(name,".spline",i,sep="")
		     data[,namei] <-  (xx-bb[i])*(xx>bb[i])
	      }
             warning(paste("breaks duplicated for=",xnames[k]))
	   }
      }
   }
}

return(data)
}# }}}

}# }}}

##' @export
"dspline<-" <- function(data,...,value) dspline(data,y=value,...)

##' Simple linear spline 
##' 
##' Simple linear spline 
##' 
##' @param x variable to make into spline
##' @param knots cut points 
##' @param num to give names x1 x2 and so forth
##' @param name name of spline expansion name.1 name.2 and so forth
##' @author Thomas Scheike
##' @keywords survival
##' @export
LinSpline <- function(x,knots,num=TRUE,name="Spline")
{# {{{

lspline <- matrix(0,length(c(x)),length(c(knots)))
for (i in seq_along(c(knots)))
{
    lspline[,i] <- (x-knots[i])*(x>knots[i])
}

lspline <- as.data.frame(lspline)
if (num==TRUE) names(lspline) <- paste(name,seq_along(c(knots)),sep="")
else if (!is.nulll(signif)) names(lspline) <- paste(name,round(c(knots),signif),sep="")

return(lspline)
}# }}}

