
library(mets)

procformdata

dreg <- function(data,y,x=NULL,z=NULL,...,level=1,
		 fun=lm,summary=summary,regex=FALSE) {# {{{

  funn <- as.character(substitute(fun))
  pp2 <- mets:::procform2(y,x=x,z=z)

  z <- mets:::procformdata(pp2$y,data)$response
  xx <- mets:::procformdata(pp2$x,data)$response
  xx <- as.matrix(xx)
  group <- NULL
  if (!is.null(pp2$z)) group <- mets:::procformdata(pp2$z,data)$response

###  print("hej")
###  if (level==0) {
###     ll <- by(data, group,do.call(fun,list(formula=y~xx)))
###  }
###  print("hh hej")

###  ll <- by(data,group, function(x) lm(z[,1] ~ xx[,1] ))
  
###   pp <- mets:::procformdata(formula,data,regex=regex)
###   z <- pp$response

   res <- sum <- c()
   if (level==1) {
       for (i in seq_len(ncol(z))) {
	   nn <- colnames(z)[i]
           y <- z[,i]
	   x <- xx
           val <- do.call(funn,c(list(y~x),list(...)))
	   val <- list(val)
###	   names(attr(val,"dimnames")) <- nn
	   names(val) <- nn
###	   val <- setNames(val,nn)
	   res <- c(res, val)
	   if (!is.null(summary)) {
	       sval <- list(do.call(summary,list(val[[1]])))
	       names(sval) <- nn
	       sum <- c(sum, sval)
	   }
       }
       res <- list(reg=res,summary=sum,...)
###       res <- list(setNames(res,funn),summary=sum,...)
       class(res) <- "dreg"
       return(res)
   }
   if (level==-1)  {
	   ny <- colnames(z)[1]
           for (i in seq_len(ncol(xx))) {
	   nn <- colnames(xx)[i]
	   nn <- paste(ny,"~",nn,sep="")
           y <- z[,1]
	   x <- xx[,i]
           val <- do.call(funn,list(y~x,...))
###	   names(attr(val,"dimnames")) <- nn
	   val <- list(val)
	   names(val) <- nn
	   res <- c(res, val)
	   if (!is.null(summary)) {
	       sval <- list(do.call(summary,list(val[[1]])))
	       names(sval) <- nn
	       sum <- c(sum, sval)
	   }
       }
       res <- list(reg=res,summary=sum)
###       res <- list(setNames(res,funn),summary=sum,...)
       class(res) <- "dreg"
       return(res)
   }

}# }}}


print.dreg <- function(x,sep="----------------------------\n",...) {# {{{
###    cat(sep)
###    if (inherits(x$lm, c("lm"))) {
###        print(x$lm)
###        if (!is.null(x$summary)) print(x$summary)
###        return(invisible(x))
###    }
    nn <-  names(x$reg)
    for (i in seq_along(x$reg)) {
        cat(paste("Response=",nn[i],"\n"))
        if (!is.null(x$summary))
            print(x$summary[[i]],...)
        else print(x$reg[[i]],...)
        cat(sep)
    }

}# }}}


procform2 <- function(y,x=NULL,z=NULL,...) {# {{{
    yx <- procform(y,return.formula=FALSE,...)
    y <- yx$response
    x0 <- yx$predictor
    z0 <- NULL
    if (length(yx$filter)>0) z0 <- yx$filter[[1]]
    if (is.null(x) && length(y)>0) x <- x0
    if (NCOL(x)==0) x <- NULL
    if (length(y)==0) y <- x0
    if (!is.null(x)) {
        x <- unlist(procform(x,return.formula=FALSE,...))
    }
    if (is.null(z)) z <- z0
    if (!is.null(z)) {
        zz <- unlist(procform(z,return.formula=FALSE,...))
    }
    if (is.null(z)) z <- z0
    return(list(y=y,x=x,z=z))
}# }}}
ll <- procform2data(iris,"*.length"~"*.width"|species)
ll <- procform2data(iris,"*.length"~"*.width")

  z <- mets:::procformdata("*.length"~"*.width"|species,data=iris,do.filter=FALSE)

data(iris)
drename(iris) <- ~.
names(iris)
iris$time <- runif(nrow(iris))
iris$time1 <- runif(nrow(iris))
iris$status <- rbinom(nrow(iris),1,0.5)
iris$S1 <- with(iris,Surv(time,status))
iris$S2 <- with(iris,Surv(time1,status))

cc=dreg(iris,"S*"~"*.width"|species,fun=coxph,level=1)
cc
cc=dreg(iris,S1~"*.width"|species,fun=coxph,level=-1)
cc

cc=dreg(iris,Surv(time,status)+Surv(time1,status)~"*.width"|species,
	fun=coxph,level=-1)
cc$reg
cc$summary[[1]]


cc=dreg(iris,Surv(time,status)+Surv(time1,status)~"*.width",fun=coxph,level=1,
	summary=NULL)
cc
###
cc=dreg(iris,Surv(time,status)+Surv(time1,status)~"*.width",fun="coxph",level=-1)
cc


dreg(iris,"*.length"~"*.width"|species,level=1)
dreg(iris,"*.length"~"*.width"|species,level=-1)

dreg(iris,"*.length"~"*.width",level=1,summary=NULL)
dreg(iris,"*.length"~"*.width",level=-1)



cc=dreg(iris,I(sepal.length>6)~"*.width"|species,level=1,fun=glm)
cc=dreg(iris,I(sepal.length>6)~"*.width"|species,level=1,fun=lm)
cc
cc$reg

cc=dreg(iris,I(sepal.length>6)~"*.width"|species,level=1,fun="glm",family=binomial())
cc$reg


dsummary(iris,~"*.length"+"*.width")

dreg(iris,"*.length"~"*.width",level=1)
cc <- dreg(iris,"*.length"~"*.width",level=1,fun=lm)
names(cc)
cc
cc <- dreg(iris,"*.length"~"*.width",level=1,fun="lm")
names(cc)
cc[[1]]


  data=iris
  y <- "*.length"~"*.width"|species
  pp2 <- mets:::procform2(y)

  z <- mets:::procformdata(pp2$y,data)$response
  xx <- mets:::procformdata(pp2$x,data)$response

  xx <- as.matrix(xx)
  group <- NULL
  if (!is.null(pp2$z)) group <- mets:::procformdata(pp2$z,data)$response


  by(z,group,summary)

  by(z,group, do.call(summary, c(list(z)))

  do.call("lm",list(list(z[,1]~xx),list(z[,2]~xx)))


  by(data,group, function(x) lm(z[,1] ~ xx ))

  lapply(1:2, function(x) lm(z[,x] ~ xx ))
  lapply(1:2, function(x) summary(lm(z[,x] ~ xx )))
  by(data,group, function(x) summary(lm(z[,1] ~ xx )))


  z <- mets:::procformdata(pp2$y,data)$response
  xx <- mets:::procformdata(pp2$x,data)$response
  xx <- as.matrix(xx)
  group <- NULL
  if (!is.null(pp2$z)) group <- mets:::procformdata(pp2$z,data)$response

   ll <- by(data,group, function(x) lm(z[,1] ~ xx ))
   ll <- by(data,group, function(x) lm(z[,1] ~ xx ))
  print(ll)

  
#
