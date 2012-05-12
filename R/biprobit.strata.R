do.biprobit.strata <- function(x,fun,...) {
  res <- lapply(x$model,function(m) do.call(fun,c(list(m),list(...))))
  names(res) <- names(x$model)
  class(res) <- "do.biprobit.strata"
  return(res)
}

##' @S3method print do.biprobit.strata
print.do.biprobit.strata <- function(x,...) {
  for (i in seq_len(length(x))) {    
    message(rep("-",60),sep="")
    message("Strata '",names(x)[i],"'",sep="")
    print(x[[i]])
  }
  return(invisible(x))
}


##' @S3method plot biprobit.strata
plot.biprobit.strata <- function(x,...)
  suppressMessages(do.biprobit.strata(x,"plot",...))

##' @S3method print biprobit.strata
print.biprobit.strata <- function(x,...)
  print.do.biprobit.strata(x$model,...)

##' @S3method summary biprobit.strata
summary.biprobit.strata <- function(object,...)
  do.biprobit.strata(object,"summary",...)

##' @S3method coef biprobit.strata
coef.biprobit.strata <- function(object,...) object$coef

##' @S3method logLik biprobit.strata
logLik.biprobit.strata <- function(object,indiv=FALSE,list=FALSE,...) {
  ll <- lapply(object$model,function(x) logLik(x,indiv=indiv,...))
  if (list) return(ll)
  if (!indiv) {
    res <- structure(sum(unlist(ll)),df=0,nall=0)
    for (i in seq(length(ll))) {
      attributes(res)$nall <- attributes(res)$nall+attributes(ll[[i]])$nall
      attributes(res)$df <- attributes(res)$df+attributes(ll[[i]])$df
    }
    attributes(res)$nobs <- attributes(res)$nall-attributes(res)$df
    class(res) <- "logLik"
    return(res)
  }
  return(unlist(ll))
}

##' @S3method score biprobit.strata
score.biprobit.strata <- function(x,...) {
  ss <- lapply(x$model,function(m) score(m,indiv=FALSE,...))
  return(unlist(ss))
}
