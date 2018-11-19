##' Regression for data frames with dutility call 
##'
##' Regression for data frames with dutility call 
##' @param data data frame
##' @param y name of variable, or fomula, or names of variables on data frame.
##' @param x name of variable, or fomula, or names of variables on data frame.
##' @param z name of variable, or fomula, or names of variables on data frame.
##' @param x.oneatatime x's one at a time
##' @param x.base.names base covarirates 
##' @param z.arg what is Z, c("clever","base","group","condition"), clever decides based on type of Z, base means that Z is used as fixed baseline covaraites for all X, group means the analyses is done based on groups of Z, and condition means that Z specifies a condition on the data
##' @param fun. function  lm is default
##' @param summary. summary to use 
##' @param regex regex
##' @param convert convert
##' @param special special's 
##' @param equal to do pairwise stuff
##' @param test development argument
##' @param ... Additional arguments for fun
##' @author Klaus K. Holst, Thomas Scheike
##' @examples##' 
##' data(iris)
##' data <- iris
##' drename(iris) <- ~.
##' names(iris)
##' set.seed(1)
##' iris$time <- runif(nrow(iris))
##' iris$time1 <- runif(nrow(iris))
##' iris$status <- rbinom(nrow(iris),1,0.5)
##' iris$S1 <- with(iris,Surv(time,status))
##' iris$S2 <- with(iris,Surv(time1,status))
##' iris$id <- 1:nrow(iris)
##' 
##' mm <- dreg(iris,"*.length"~"*.width"|I(species=="setosa" & status==1))
##' mm <- dreg(iris,"*.length"~"*.width"|species+status)
##' mm <- dreg(iris,"*.length"~"*.width"|species)
##' mm <- dreg(iris,"*.length"~"*.width"|species+status,z.arg="group")
##' 
##' \donttest{ ## Reduce Ex.Timings
##' y <- "S*"~"*.width"
##' xs <- dreg(iris,y,fun.=phreg)
##' xs <- dreg(iris,y,fun.=survdiff)
##' 
##' y <- "S*"~"*.width"
##' xs <- dreg(iris,y,x.oneatatime=FALSE,fun.=phreg)
##' 
##' ## under condition
##' y <- S1~"*.width"|I(species=="setosa" & sepal.width>3)
##' xs <- dreg(iris,y,z.arg="condition",fun.=phreg)
##' xs <- dreg(iris,y,fun.=phreg)
##' 
##' ## under condition
##' y <- S1~"*.width"|species=="setosa"
##' xs <- dreg(iris,y,z.arg="condition",fun.=phreg)
##' xs <- dreg(iris,y,fun.=phreg)
##' 
##' ## with baseline  after |
##' y <- S1~"*.width"|sepal.length
##' xs <- dreg(iris,y,fun.=phreg)
##' 
##' ## by group by species, not working
##' y <- S1~"*.width"|species
##' ss <- split(iris,paste(iris$species,iris$status))
##' 
##' xs <- dreg(iris,y,fun.=phreg)
##' 
##' ## species as base, species is factor so assumes that this is grouping
##' y <- S1~"*.width"|species
##' xs <- dreg(iris,y,z.arg="base",fun.=phreg)
##' 
##' ##  background var after | and then one of x's at at time
##' y <- S1~"*.width"|status+"sepal*"
##' xs <- dreg(iris,y,fun.=phreg)
##' 
##' ##  background var after | and then one of x's at at time
##' ##y <- S1~"*.width"|status+"sepal*"
##' ##xs <- dreg(iris,y,x.oneatatime=FALSE,fun.=phreg)
##' ##xs <- dreg(iris,y,fun.=phreg)
##' 
##' ##  background var after | and then one of x's at at time
##' ##y <- S1~"*.width"+factor(species)
##' ##xs <- dreg(iris,y,fun.=phreg)
##' ##xs <- dreg(iris,y,fun.=phreg,x.oneatatime=FALSE)
##' 
##' y <- S1~"*.width"|factor(species)
##' xs <- dreg(iris,y,z.arg="base",fun.=phreg)
##' 
##' y <- S1~"*.width"|cluster(id)+factor(species)
##' xs <- dreg(iris,y,z.arg="base",fun.=phreg)
##' xs <- dreg(iris,y,z.arg="base",fun.=coxph)
##' 
##' ## under condition with groups
##' y <- S1~"*.width"|I(sepal.length>4)
##' xs <- dreg(subset(iris,species=="setosa"),y,z.arg="group",fun.=phreg)
##' 
##' ## under condition with groups
##' y <- S1~"*.width"+I(log(sepal.length))|I(sepal.length>4)
##' xs <- dreg(subset(iris,species=="setosa"),y,z.arg="group",fun.=phreg)
##' 
##' y <- S1~"*.width"+I(dcut(sepal.length))|I(sepal.length>4)
##' xs <- dreg(subset(iris,species=="setosa"),y,z.arg="group",fun.=phreg)
##' 
##' ff <- function(formula,data,...) {
##'  ss <- survfit(formula,data,...)
##'  kmplot(ss,...)
##'  return(ss)
##' }
##'
##' if (interactive()) {
##' dcut(iris) <- ~"*.width"
##' y <- S1~"*.4"|I(sepal.length>4)
##' par(mfrow=c(1,2))
##' xs <- dreg(iris,y,fun.=ff)
##' }
##' }
##' 
##' @export
dreg <- function(data,y,x=NULL,z=NULL,x.oneatatime=TRUE,
	 x.base.names=NULL,z.arg=c("clever","base","group","condition"),
         fun.=lm,summary.=summary,regex=FALSE,convert=NULL,
	 special=NULL,equal=TRUE,test=1,...) {# {{{
### z.arg=clever,  if z is logical then condition
###                if z is factor  then group variable 
###                if z is numeric then baseline covariate 
### ... further arguments to fun 

###  fun <- as.character(substitute(fun))
###    if (is.character(fun))
###        fun <- get(fun)
###    if (!is.null(convert) && is.logical(convert)) {
###        if (convert)
###            convert <- as.matrix
###        else convert <- NULL
###    }
###    if (!is.null(convert)) {
###        fun_ <- fun
###        fun <- function(x, ...) fun_(convert(x, ...))
###    }

###    print(fun)
###    print(str(fun))


 yxzf <- procform(y,x=x,z=z,data=data,do.filter=FALSE,regex=regex)
 yxz <- procformdata(y,x=x,z=z,data=data,do.filter=FALSE,regex=regex)

 ## remove blank, to able to use also 	+1 on right hand side
 if (any(yxzf$predictor=="")) 
    yxzf$predictor <- yxzf$predictor[-which(yxzf$predictor=="")]


 yy <- yxz$response
 xx <- yxz$predictor
 ### group is list, so zz is data.frame 
 if ((length(yxzf$filter))==0) zz <- NULL else if ((length(yxzf$filter[[1]])==1 & yxzf$filter[[1]][1]=="1")) 
      zz <- NULL else   zz <- yxz$group[[1]]

 if (!is.null(zz)) {# {{{
 if (z.arg[1]=="clever")
 {
   if ((ncol(zz)==1) & is.logical(zz[1,1])) z.arg[1] <- "condition" 
       else if ((ncol(zz)==1) & is.factor(zz[,1])) z.arg[1] <- "group" 
       else  z.arg[1] <- "base"
  }
  }# }}}
### print(z.arg)


 basen <- NULL
 if (z.arg[1]=="base")
     basen <- yxzf$filter[[1]]
 if (z.arg[1]=="condition")
     data <- subset(data,eval(yxzf$filter.expression))
 if (z.arg[1]=="group")
     group <- interaction(zz) else group <- rep(1,nrow(data))
 if (z.arg[1]=="group")  levell <- levels(group) else levell <-1 


 res <- sum <- list()
 if (test==1) {
 if (is.null(summary)) sum <- NULL
 for (g in levell) {# {{{
 if (equal==TRUE) datal <- subset(data,group==g)
 else datal <- subset(data,group!=g)
 for (y in yxzf$response) {# {{{
	 if (x.oneatatime)  {
		 for (x in yxzf$predictor) {
	             if (length(c(x,basen))>1) 
                     basel <- paste(c(x,basen),collapse="+")
		     else basel <- c(x,basen)
	             form <- as.formula(paste(y,"~",basel))
		     if (!is.null(special)) form <- timereg::timereg.formula(form,special=special)
###		     val <- with(data,do.call(fun,c(list(formula=form),list(...))))
	     capture.output(
             val <- do.call(fun.,c(list(formula=form),list(data=datal),list(...))))
	     val$call <- paste(y,"~",basel)
             val <- list(val)
	     nn <- paste(y,"~",basel)
	     if (z.arg[1]=="group") {
		     if (equal==TRUE) nn <- paste(nn,"|",g)  else nn <- paste(nn,"| not",g);
	     }
	     names(val) <- nn
             res <- c(res, val)
	     if (!is.null(summary)) {
	        sval <- list(do.call(summary.,list(val[[1]])))
                names(sval) <- nn
	        sum <- c(sum, sval)
	    }
	 }
	 } else {
             basel <- paste(c(yxzf$predictor,basen),collapse="+")
             form <- as.formula(paste(y,"~",basel))
	     if (!is.null(special)) form <- timereg::timereg.formula(form,special=special)
	     capture.output(
             val <- do.call(fun.,c(list(formula=form),list(data=datal),list(...))))
	     nn <- paste(y,"~",basel)
	     if (z.arg[1]=="group") {
		     if (equal==TRUE) nn <- paste(nn,"|",g)  else nn <- paste(nn,"| not",g);
	     }
	     val$call <- nn
             val <- list(val)
	     names(val) <- paste(y,"~",basel)
             res <- c(res, val)
	     if (!is.null(summary)) {
	       sval <- list(do.call(summary.,list(val[[1]])))
	       names(sval) <- nn
	       sum <- c(sum, sval)
	     }

	 }
 }# }}}
 }# }}}
 }

   res <- list(reg=res,summary=sum)
###       res <- list(setNames(res,funn),summary=sum,...)
   class(res) <- "dreg"
###   structure(res,ngrouvar=0,class="dreg")
   return(res)

}# }}}

##' @export
print.dreg <- function(x,sep="-",...) {# {{{
    sep <- paste(rep(sep,50,sep=""),collapse="")
    sep <-  paste(sep,"\n")
    nn <-  names(x$reg)
    for (i in seq_along(x$reg)) {
        cat(paste("Model=",nn[i],"\n"))
            print(x$reg[[i]],...)
        cat(sep)
    }

}# }}}

##' @export
summary.dreg <- function(object,sep="-",...) {# {{{
    x <- object
    sep <- paste(rep(sep,50,sep=""),collapse="")
    sep <-  paste(sep,"\n")
###    cat(sep)
###    if (inherits(x$lm, c("lm"))) {
###        print(x$lm)
###        if (!is.null(x$summary)) print(x$summary)
###        return(invisible(x))
###    }
if (!is.null(x$summary)) {
    nn <-  names(x$summary)
    for (i in seq_along(x$summary)) {
        cat(paste("Model=",nn[i],"\n"))
        if (!is.null(x$summary))
            print(x$summary[[i]],...)
        else print(x$reg[[i]],...)
        cat(sep)
    }
}

}# }}}

