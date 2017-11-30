
dreg <- function(data,y,x=NULL,z=NULL,...,x.oneatatime=TRUE,
	 x.base.names=NULL,z.arg=c("clever","base","group","condition"),
         fun=lm,summ=TRUE,regex=FALSE,convert=NULL) {# {{{
### z.arg=clever,  if z is logical then condition
###                if z is factor  then group variable 
###                if z is numeric then baseline covariate 
### ... further arguments to fun 

### funn <- as.character(substitute(fun))
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

 yxzf <- mets:::procform(y,x=x,z=z,data=data,do.filter=FALSE,regex=regex)
 yxz <- mets:::procformdata(y,x=x,z=z,data=data,do.filter=FALSE,regex=regex)
### print(yxzf)

 yy <- yxz$response
 xx <- yxz$predictor
 ### group is list, so zz is data.frame 
 if ((length(yxzf$filter[[1]])==1 & yxzf$filter[[1]][1]=="1")) zz <- NULL else   zz <- yxz$group[[1]]

 if (!is.null(zz)) {# {{{
 if (z.arg[1]=="clever")
 {
   if ((ncol(zz)==1) & is.logical(zz[1,1])) z.arg[1] <- "condition" 
       else if ((ncol(zz)==1) & is.factor(zz[,1])) z.arg[1] <- "group" 
       else  z.arg[1] <- "base"
  }
  }# }}}
### print(z.arg)

###  ## some of the x covariates also baseline variables 
###  if (!is.null(x.base.names)) {
###          base.x <- xx[,x.base.names]
###	  xx <- xx[,-base.names]
###	  ### combined baseline, possibly also some in zz 
###	  base <- cbind(base,base.x)
###  } 
###  if (!is.null(condition)) 
###  {
###      yy <- subset(yy,condition==TRUE)
###      xx <- subset(xx,condition==TRUE)
###      if (!is.null(base)) base <- subset(base,condition==TRUE)
###  }

 basen <- NULL
 if (z.arg[1]=="base")
     basen <- yxzf$filter[[1]]
 if (z.arg[1]=="condition")
     data <- subset(data,eval(yxzf$filter.expression))

 res <- sum <- list()
 for (y in yxzf$response) {# {{{
	 if (x.oneatatime)  {
		 for (x in yxzf$predictor) {
	             if (length(c(x,basen))>1) 
                     basel <- paste(c(x,basen),collapse="+")
		     else basel <- c(x,basen)
	             form <- as.formula(paste(y,"~",basel))
###		     val <- with(data,do.call(fun,c(list(formula=form),list(...))))
	     capture.output(
             val <- do.call(fun,c(list(formula=form),list(data=data),list(...))))
	     val$call <- paste(y,"~",basel)
                     val <- list(val)
	             names(val) <- paste(y,"~",basel)
                    res <- c(res, val)
	   if (!is.null(summary)) {
	       sval <- list(do.call(summary,list(val[[1]])))
               names(sval) <- paste(y,"~",basel)
	       sum <- c(sum, sval)
	   }
	}
	 } else {
             basel <- paste(c(yxzf$predictor,basen),collapse="+")
             form <- as.formula(paste(y,"~",basel))
	     capture.output(
             val <- do.call(fun,c(list(formula=form),list(data=data),list(...))))
	     val$call <- paste(y,"~",basel)
             val <- list(val)
	     names(val) <- paste(y,"~",basel)
             res <- c(res, val)
	   if (!is.null(summary)) {
	       sval <- list(do.call(summary,list(val[[1]])))
	     names(sval) <- paste(y,"~",basel)
	       sum <- c(sum, sval)
	   }

	 }
 }# }}}

###  res <- sum <- list()
###  for (resp in names(yy)) { ## {{{ 
###     resp <- "S1"
###     yr <- yy[,resp]
###
###      if (x.oneatatime & (length(x.names)>1)) { ## over arguments in X# {{{
###        for (i in seq_len(ncol(xx))) {# {{{
######	   xn <- x.names(x)
######		   nn <- paste(ny,"~",nn,sep="")
######		   y <- z[,1]
###	   x <- xx[,i]
###	   if (!is.null(base)) x <- cbind(x,base)
######		   val <- do.call(funn,list(y~x,...))
###	###	   names(attr(val,"dimnames")) <- nn
######		   val <- lm(y~x,...)
###           xs <- names(x)
######           form <- as.formula(paste(resp,"~",paste(xs,collapse="+")))
###	   nn <- names(x)
######	   nn <- c(paste(resp,"~"),names(x))
######	   val <- coxph(form,data=data)
######	   x <- as.matrix(x)
######	   if (!is.null(group)) val <- by(group,x,fun)
###	   if (!is.null(group)) val <- by(data,group,function(x) coxph(yr~x))
###	   else val <- coxph(yr~x)
###	   val <- list(val)
######	   names(val) <- nn
###	   res <- c(res, val)
###	   if (!is.null(summary)) {
###	       sval <- list(do.call(summary,list(val[[1]])))
######	       names(sval) <- nn
###	       sum <- c(sum, sval)
###	   }
###       }# }}}
###   }# }}}
###   else { ## all x's at once 
###        x <- xx
###        if (!is.null(base)) x <- cbind(x,base)
###	   nn <- paste(resp,"~",names(x))
###           x <- as.matrix(x)
###	   val <- coxph(yr~x)
###	   val <- list(val)
######	   names(val) <- nn
###	   res <- c(res, val)
###	   if (!is.null(summary)) {
###	       sval <- list(do.call(summary,list(val[[1]])))
######	       names(sval) <- nn
###	       sum <- c(sum, sval)
###	   }
###   }
###
###  } ## }}}
###
   res <- list(reg=res,summary=sum)
###       res <- list(setNames(res,funn),summary=sum,...)
   class(res) <- "dreg"
   res
###   structure(res,ngrouvar=0,class="dreg")
###   return(res)

}# }}}

print.dreg <- function(x,sep="----------------------------\n",...) {# {{{
    nn <-  names(x$reg)
    for (i in seq_along(x$reg)) {
        cat(paste("Model=",nn[i],"\n"))
            print(x$reg[[i]],...)
        cat(sep)
    }

}# }}}

summary.dreg <- function(x,sep="----------------------------\n",...) {# {{{
###    cat(sep)
###    if (inherits(x$lm, c("lm"))) {
###        print(x$lm)
###        if (!is.null(x$summary)) print(x$summary)
###        return(invisible(x))
###    }
    nn <-  names(x$summary)
    for (i in seq_along(x$summary)) {
        cat(paste("Model=",nn[i],"\n"))
        if (!is.null(x$summary))
            print(x$summary[[i]],...)
        else print(x$reg[[i]],...)
        cat(sep)
    }

}# }}}


stop()

library(mets)
data(iris)
data=iris
drename(iris) <- ~.
names(iris)
iris$time <- runif(nrow(iris))
iris$time1 <- runif(nrow(iris))
iris$status <- rbinom(nrow(iris),1,0.5)
iris$S1 <- with(iris,Surv(time,status))
iris$S2 <- with(iris,Surv(time1,status))
iris$id <- 1:nrow(iris)

mm <- dreg(iris,"*.length"~"*.width"|I(species=="setosa" & status==1))
mm <- dreg(iris,"*.length"~"*.width"|species+status)
mm <- dreg(iris,"*.length"~"*.width")

x=NULL;z=NULL;level=1;
base=NULL;
z.arg=c("clever","base","group","condition");
fun=lm;summary=summary;regex=FALSE
###
data=iris
x=NULL;z=NULL;
x.oneatatime=TRUE
x.base.names=NULL;z.arg=c("clever","base","group","condition");
fun=lm;summ=TRUE;regex=FALSE
regex=FALSE
z.arg="base"
z.arg="clever"
convert = NULL


### testing forskellige calls
y <- "S*"~"*.width"
xs <- dreg(iris,y,fun=phreg)
xs <- dreg(iris,y,fun=survdiff)

### testing forskellige calls
y <- "S*"~"*.width"
xs <- dreg(iris,y,x.oneatatime=FALSE,fun=phreg)

## under condition 
y <- S1~"*.width"|I(species=="setosa" & sepal.width>3)
xs <- dreg(iris,y,z.arg="condition",fun=phreg)
xs <- dreg(iris,y,fun=phreg)

## under condition 
y <- S1~"*.width"|species=="setosa"
xs <- dreg(iris,y,z.arg="condition",fun=phreg)
xs <- dreg(iris,y,fun=phreg)

## with baseline  after | 
y <- S1~"*.width"|sepal.length
xs <- dreg(iris,y,fun=phreg)

## by group by species, not working 
y <- S1~"*.width"|species
xs <- dreg(iris,y,fun=phreg)

## species as base, species is factor so assumes that this is grouping 
y <- S1~"*.width"|species
xs <- dreg(iris,y,z.arg="base",fun=phreg)


##  background var after | and then one of x's at at time 
y <- S1~"*.width"|status+"sepal*"
xs <- dreg(iris,y,fun=phreg)

##  background var after | and then one of x's at at time 
y <- S1~"*.width"|status+"sepal*"
xs <- dreg(iris,y,x.oneatatime=FALSE,fun=phreg)
xs <- dreg(iris,y,fun=phreg)

##  background var after | and then one of x's at at time 
y <- S1~"*.width"+factor(species)
xs <- dreg(iris,y,fun=phreg)
xs <- dreg(iris,y,fun=phreg,x.oneatatime=FALSE)

y <- S1~"*.width"|factor(species)
xs <- dreg(iris,y,z.arg="base",fun=phreg)


y <- S1~"*.width"|cluster(id)+factor(species)
xs <- dreg(iris,y,z.arg="base",fun=phreg)
xs <- dreg(iris,y,z.arg="base",fun=coxph)


## under condition with groups  
y <- S1~"*.width"|I(sepal.length>4)
xs <- dreg(subset(iris,species=="setosa"),y,z.arg="group",fun=phreg)

## under condition with groups  
y <- S1~"*.width"+I(log(sepal.length))|I(sepal.length>4)
xs <- dreg(subset(iris,species=="setosa"),y,z.arg="group",fun=phreg)

y <- S1~"*.width"+I(dcut(sepal.length))|I(sepal.length>4)
xs <- dreg(subset(iris,species=="setosa"),y,z.arg="group",fun=phreg)

ff <- function(formula,data,...) {
   ss <- survfit(formula,data,...)
   kmplot(ss,...)
   return(ss)
}

dcut(iris) <- ~"*.width"
y <- S1~"*.4"|I(sepal.length>4)
par(mfrow=c(1,2))
xs <- dreg(iris,y,fun=ff)


#########################################

