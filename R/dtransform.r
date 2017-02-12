
##' @export
dtransform <- function(data,...,subset)
{# {{{
    if (is.vector(data)) data <- data.frame(data)
    subs <- substitute(subset)
    if (!base::missing(subs)) {
        expr <- suppressWarnings(inherits(try(subset,silent=TRUE),"try-error"))
        if (expr) data <- data[which(eval(subs,envir=data)),,drop=FALSE]
        else data[subset,,drop=FALSE]
    }

    datanew <- transform(data,...)
return(datanew)
}# }}}

###xx <- dtransform(mena,ll=log(agemena)+twinnum,zyg=="DZ")
###xx <- dtransform(mena,ll=log(agemena)+twinnum,I(agemena)<15)

##' @export
"dtransform<-" function(data,subset,value) dtransform(data,value,subet)


