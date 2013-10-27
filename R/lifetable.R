##' Creat simple life table 
##'
##' 
##' @title Life table
##' @param time Time (or matrix with columns time,status or entry,time,status)
##' @param status status (event=TRUE,censoring=FALSE)
##' @param entry Entry time
##' @param strata Strata
##' @param breaks Time intervals
##' @param confint If TRUE 95% confidence limits are calculated
##' @author Klaus K. Holst
##' @examples
##' library(timereg)
##' data(TRACE)
##' \donttest{
##' with(TRACE, lifetable(cbind(time,status==9),breaks=c(0.2,0.5),confint=TRUE))
##' }
##' d <- with(TRACE, lifetable(time,status==9,strata=list(sex=sex,vf=vf),breaks=c(0.2,0.5)))
##' summary(glm(events ~ offset(log(atrisk))+interval+sex+vf,data=d,poisson))
##' 
##' @export
lifetable <- function(time,status,entry,strata=list(),breaks=c(),confint=FALSE) {
    if (missing(entry)) entry <- rep(0,NROW(time))
    if ((is.matrix(time) || is.data.frame(time)) && ncol(time)>1) {
        if (ncol(time)==3) {
            status <- time[,3]
            entry <- time[,1]
            time <- time[,2]
        } else {
            status <- time[,2]
            time <- time[,1]
            entry <- rep(0,length(time))
        }
    }
    if (length(strata)>0) {
        a <- by(cbind(entry,time,status), strata,
                FUN=lifetable, breaks=breaks, confint=confint)
        nn <- do.call("expand.grid",attributes(a)$dimnames)
        nam <- nn[rep(seq(NROW(nn)),each=NROW(a[[1]])),]
        ##colnames(nam) <- colnames(nn)
        res <- cbind(Reduce("rbind",a),nam)
        return(res)
    }    
    if (length(breaks)>0) breaks <- sort(unique(breaks))
    en <- matrix(unlist(lapply(c(-Inf, breaks),function(x) pmax(x,entry))),
                  ncol=length(breaks)+1)
    ex <- matrix(unlist(lapply(c(breaks, Inf),function(x) pmin(x,time))),
                  ncol=length(breaks)+1)    
    dur <- ex-en
    endur <- rbind(c(breaks,Inf))%x%cbind(rep(1,nrow(en)))-en
    dur[dur<0] <- NA
    enter <- colSums(!is.na(dur))
    atrisk <- colSums(dur,na.rm=TRUE)
    eventcens <- dur<endur
    eventcens <- rbind(apply(dur<endur,2,function(x) x*(status+1)))
    lost <- colSums(eventcens==1,na.rm=TRUE)
    events <- colSums(eventcens==2,na.rm=TRUE)
    res <- data.frame(enter=enter,
                      atrisk=atrisk,
                      lost=lost,
                      events=events,
                      interval=as.factor(c(breaks,Inf)),
                      rate=events/atrisk)
    if (confint) {
        f <- events ~ offset(log(atrisk))
        if (length(breaks)>0) f <- update(f,.~.+factor(interval)-1)
        g <- glm(f,data=res,poisson)
        ci <- rbind(exp(stats::confint(g)))
        res[,"2.5%"] <- ci[,1]
        res[,"97.5%"] <- ci[,2]
    }
    res
}
