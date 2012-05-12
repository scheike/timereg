##' .. content for description (no empty lines) ..
##'
##' .. content for details ..
##' @title Calculate 
##' @param formula 
##' @param data 
##' @param ... 
##' @param time 
##' @param timestrata 
##' @param cumulative 
##' @param silent 
##' @author Klaus K. Holst and Thomas Scheike
##' @export
cumh <- function(formula,data,...,time,
                 timestrata=quantile(data[,time],c(0.25,0.5,0.75,1)),
                 cumulative=TRUE,
                 silent=FALSE) {
  
  time.  <- substitute(time)
  if (!is.character(time.)) time. <- deparse(time.)
  time <- time.
  
  res <- list(); i <- 0
  ht <- c()
  outcome <- as.character(terms(formula)[[2]])
  y0 <- data[,outcome]
  for (i in seq(length(timestrata))) {
    t <- timestrata[i]
    data[,outcome] <- y0
    newdata <- data
    if (!cumulative) {
      if (i==1) {
        idx <- data[,time]<t
      } else {
        idx <- (timestrata[i-1]<=data[,time] & data[,time]<t)
      }
    } else {
      data[,outcome] <- data[,outcome]*(data[,time]<t)
    }
    if (!silent) {
      message(t)##," ",sum(data[,outcome]))
    }
    if (!cumulative)
      res[[i]] <-bptwin(formula,data=data[idx,],...)
    else
      res[[i]] <-bptwin(formula,data=data,...)
    ht <- rbind(ht,c(t,summary(res[[i]])$h[1,]))
  }
  rownames(ht) <- timestrata
  colnames(ht) <- c("time","Heritability","Std.Err","2.5%","97.5%")
  res <- (list(ht=ht,models=res))
  class(res) <- "cumh"
  res
}

##' @S3method summary cumh
summary.cumh <- function(object,...) object 

##' @S3method print cumh
print.cumh <- function(x,...) {
  print(x$ht)
  invisible(x)
}

Col <- function (col, alpha = 0.2) {
    sapply(col, function(x) do.call(rgb, as.list(c(col2rgb(x)/255, 
        alpha))))
}

##' @S3method plot cumh
plot.cumh <- function(x,...,idx=seq(nrow(x$ht)),lwd=2,col,fillcol,alpha=0.2,ylim=c(0,1),xlab="Time",ylab="Heritability",add=FALSE) {

  if (missing(col)) col <- "darkblue"
  if (alpha>0 & missing(fillcol)) fillcol <- Col(col,alpha)
  if (!add) {
    plot(x$ht[idx,1:2,drop=FALSE],type="l",ylim=ylim,lwd=lwd,
         ylab=ylab,xlab=xlab,col=col,...)
  }
  xx <- with(x, c(ht[idx,1],rev(ht[idx,1])))
  yy <- with(x, c(ht[idx,4],rev(ht[idx,5])))           
  polygon(xx,yy,col=fillcol)
  lines(x$ht[idx,1:2,drop=FALSE],lwd=lwd,col=col,...)
  invisible(x)
}
