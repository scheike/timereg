Event <- function(time2,cause,cens.code=0,...) {
out <- cbind(time2,cause)
class(out) <- "Event"
attr(out,"cens.code") <- cens.code
return(out)
}
