
event.split <- function(data,
			time="time",status="status",cuts="cuts",
			name.start="start",name.id="id",
			entry=NULL, cens.code=0,num="num")
{
## {{{ 
    n <- nrow(data)
    new.time <- data[,time]
    new.status <- data[,status]
    new.cuts <- data[,cuts]
    if (!is.null(entry)) {
      new.start <- data[,entry]
    } else new.start <- rep(0,n)
    if (num %in% names(data)) new.num <- data[,num] else 
	    new.num <- rep(0,n)
    id <- 1:n

    splits <- (new.cuts<new.time & new.start<new.cuts)
    id <- c(id,id[splits])
    new.time <- c(new.time,new.time[splits])
    new.start <- c(new.start,new.cuts[splits])
    new.status <- c(new.status,new.status[splits])
    new.num <- c(new.num,new.num[splits]+1)
    splits <- c(splits,rep(FALSE,sum(splits)))
    new.time[splits] <- new.cuts[splits]
    new.status[splits] <- cens.code

    data <- data[id,]
    data[,time] <- new.time
    data[,status] <- new.status
    data[,name.start] <- new.start
    data[,name.id] <- id
    data[,num] <- new.num
    data <- dsort(data,id)

    return(data)
    ## }}} 
} 


