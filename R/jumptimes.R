##' @export
jumptimes <- function(time, status=TRUE, 
                      id,cause,
                      sample,
                      strata=NULL,num=NULL, ...) {
    if (missing(id)) {
        time <- if (missing(cause)) time[status>0]
                else time[status==cause]
    } else {
        ww <- na.omit(fast.reshape(cbind(time=time,status=status),id=id,num=num))
        statusvar <- grep("status",colnames(ww))
        timevar <- grep("time",colnames(ww))
        if (missing(cause)) {
            idx <- which(rowSums(ww[,statusvar]>0)==length(statusvar))
        } else {
            idx <- which(apply(as.matrix(ww[,statusvar]),1,function(x) all(x==cause)))
        }
            time <- na.omit(do.call(pmax,as.list(ww[idx,timevar])))
    }

    if (!missing(sample)) {
        time <- sort(time)
        t0 <- seq(min(time),max(time),length.out=min(sample,length(time)))
        ii <- fast.approx(time,t0)
        dup <- duplicated(ii)
        ii[dup] <- ii[dup]-1
        time <- unique(time[ii])
    }
    return(time)
}

## with(prt, jumptimes(time,status==2,id=id,sample=10))
