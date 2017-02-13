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
"dtransform<-" <- function(data,value) dtransform(data,value)

##' @export
"dtrans<-" <- function(data,value) dtransform(data,value)



#####' @export
###dtransform <- function(data,...,subset)
###{# {{{
###    if (is.vector(data)) data <- data.frame(data)
###    subs <- substitute(subset)
###    if (!base::missing(subs)) {
###        expr <- suppressWarnings(inherits(try(subset,silent=TRUE),"try-error"))
###        if (expr) data <- data[which(eval(subs,envir=data)),,drop=FALSE]
###        else data[subset,,drop=FALSE]
###    }
###
###    datanew <- transform(data,...)
###return(datanew)
###}# }}}



