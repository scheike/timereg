
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


##' @export
dtransform <- function(data,...)
{# {{{
    if (is.vector(data)) data <- data.frame(data)

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





###xx <- dtransform(mena,ll=log(agemena)+twinnum,zyg=="DZ")
###xx <- dtransform(mena,ll=log(agemena)+twinnum,I(agemena)<15)

##' @export
"dtransform<-" function(data,subset,value) dtransform(data,value,subet)


