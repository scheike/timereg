
#' @export
plotConfregion <- function(sfit,add=TRUE,polygon=TRUE,cols=1,ltys=1,...)
{# {{{
  nn <- names(sfit$strata)
  if (is.null(nn)) ll <- 1 else ll <- length(nn)

  if (length(cols)!=ll) cols <- seq(cols[1],ll) else cols <- cols
  if (length(cols)!=ll) ltys <- seq(ltys[1],ll)  else ltys <- ltys

ss <- 1
for (i in seq(ll))
{
   if (ll>1 & i>1) ss <- ss+sfit$strata[i-1]
   if (ll==1) index <- 1:length(sfit$time)
   else index <- ss:(ss+sfit$strata[i]-1)

   nl <- cbind(sfit$time[index],sfit$lower[index])
   ul <- cbind(sfit$time[index],sfit$upper[index])

   ## if call is from survfit type with hazard plot "fun="cumhaz"))
   ## checks if ... contains "fun="cumhaz"
   if (hasArg("fun")) {
      nl <- cbind(sfit$time[index],-log(sfit$lower[index]))
      ul <- cbind(sfit$time[index],-log(sfit$upper[index]))
   }

  if (!polygon) {
      lines(nl,type="s",col=cols[i],lty=ltys[i]+1,lwd=3,...)
      lines(ul,type="s",col=cols[i],lty=ltys[i]+1,lwd=3,...)
  } else {
	 ll <- length(nl[,1])
         timess <- nl[,1]
         ttp <- c(timess[1],rep(timess[-c(1,ll)],each=2),timess[ll])
         tt <- c(ttp,rev(ttp))
         yy <- c(rep(nl[-ll,2],each=rep(2)),rep(rev(ul[-ll,2]),each=2))
         col.alpha<-0.1
         col.ci<-cols[i]
         col.trans <- sapply(col.ci, FUN=function(x) 
           do.call(grDevices::rgb,as.list(c(grDevices::col2rgb(x)/255,col.alpha))))
	 polygon(tt,yy,lty=0,col=col.trans)
      }
}

}# }}}

#' @export
kmplot<- function(x,add=FALSE,loc=NULL,col=NULL,lty=NULL,conf.int=TRUE,polygon=TRUE,add.legend=TRUE,...)
{ ## {{{ 
	### default location if loc not given 
	if (is.null(loc)) { 
	    if (min(x$surv)>0.7) loc <- "bl" else loc <- "bl"
            if (hasArg("fun")) loc <- "tl"
	}
	if (loc=="bl") loc <- "bottomleft"
	else if (loc=="br") loc <- "bottomright"
	else if (loc=="tr") loc <- "topright"
	else if (loc=="tl") loc <- "topleft"
	else loc <- "bottomleft"
	nn <- names(x$strata)
	if (is.null(nn)) ll <- 1 else ll <- length(nn)
        if (is.null(col)) cols <- seq(ll)  else cols <- col
        if (is.null(lty)) ltys <- seq(ll)  else ltys <- lty
	plot(x,col=cols,lty=ltys,conf.int=FALSE,...)
	if (!is.null(nn) & add.legend) legend(loc,legend=names(x$strata),col=cols,lty=ltys)
	if (conf.int) {
	    plotConfregion(x,add=TRUE,ltys=ltys,cols=cols,polygon=polygon,...)
	}

} ## }}}


###  library(mets)
###  library(lava)  
###
###  data(melanoma)
###  dhead(melanoma)
###  dtable(melanoma,~status)
###  melanoma <- dtransform(melanoma,status=0,status==2)
###  melanoma <- dtransform(melanoma,stat1=(status==1)*1)
###  melanoma <- dtransform(melanoma,time=days/365.25)
###
###
### fit=survfit(Surv(time,status!=0)~sex,data=melanoma)
### par(mfrow=c(1,2))
### kmplot(fit)
### kmplot(fit,fun="cumhaz")
