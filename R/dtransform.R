##' Transform that allows condition 
##'
##' Defines new variables under condition for data frame 
##' @param data is data frame 
##' @param ... new variable definitions including possible if condition
##' @examples
##' data(mena)
##'
##' xx <- dtransform(mena,ll=log(agemena)+twinnum)
##'
##' xx <- dtransform(mena,ll=log(agemena)+twinnum,agemena<15)
##' xx <- dtransform(xx  ,ll=100+agemena,ll2=1000,agemena>15)
##' dsummary(xx,ll+ll2~I(agemena>15))
##' @aliases dtransform dtransform<- dtrans dtrans<- 
##' @export
dtransform <- function(data,...)
{# {{{
    if (is.vector(data)) data <- data.frame(data)

###    if (is.list(...)) e <- eval(substitute(...), data, parent.frame()) else 
###    if (!missing(EXPRLIST)) {
###	   e <- eval(substitute(c(list(...),EXPRLIST)), data, parent.frame())
###    } else
    e <- eval(substitute(list(...)), data, parent.frame())

    tags <- names(e)
    condn  <- match("",tags) 

    if (!is.na(condn)) {
	    condition <- TRUE
	    cond <- e[[condn[1]]]; 
            whereT <- which(cond) 
	    e[[condn]] <- NULL
	    tags <- tags[-condn]
    } else condition <- FALSE
    inx <- match(tags, names(data))
    matched <- !is.na(inx)
    matchedtags <- seq(length(e))[matched]

    if (any(matched)) {
      ### new values replaces old values 
      k <- 1
      if (condition==TRUE) {
	   for (i in inx[matched]) 
	   {
            mk <- matchedtags[k]
	    if (length(e[[mk]])==1) data[whereT,i] <- e[[mk]] else 
		                    data[whereT,i] <- e[[mk]][whereT]
	   k <- k+1
	   }
      } else data[inx[matched]] <- e[matched]
      data <- data.frame(data)
    }

    ### no matched, all new variables 
    if (!all(matched)) {
	    if (condition) 
	    for (i in seq(length(e))[!matched])  {
              if (length(e[[i]])==1) e[[i]] <- rep(e[[i]],nrow(data)) 
	      e[[i]][!cond] <- NA
	    } 
	    data <- cbind(data,data.frame(e[!matched]))
    } 

return(data)
}# }}}


##' @export
dtrans <- function(data,...) dtransform(data,...)

##' @export
`dtrans<-` <- function(data,...,value) {
    dtransform(data,...) <- value
    return(data)
}


##' @export
`dtransform<-` <- function(data,...,value) {
    cl <- match.call()
    cl[[1L]] <- substitute(dtransform)
    a <- substitute(value)
    if (inherits(value,"function")) {
        cl["value"] <- NULL
        names(cl)[names(cl)=="INPUT"] <- ""
        cl[["INPUT"]] <- value
    } else {
        if (is.list(value)) {
            cl[which(names(cl)=="value")] <- NULL
            start <- length(cl)
            for (i in seq_along(value)) {
                cl[start+i] <- value[i]
            }
            if (length(names(value))>0)
                names(cl)[start+seq_along(value)] <- names(value)
        } else {
            names(cl)[which(names(cl)=="value")] <- ""
        }
    }
    eval.parent(cl)

}
