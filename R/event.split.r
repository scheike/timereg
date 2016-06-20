event.split <- function(data,
			time="time",status="status",cuts="cuts",
			name.id="id",
			name.start="start", 
			cens.code=0, 
			order.id=TRUE,time.group=TRUE
			)
{
## {{{ 
    n <- nrow(data)
    new.time <- data[,time]
    new.status <- data[,status]

    if (is.numeric(cuts)) {
	    cutname <- paste("cut",cuts,sep=".")
            data[,cutname] <- cuts
    } else cutname <- cuts
    new.cuts <- data[,cutname]

    if (is.numeric(name.start)) {
	    start0 <- name.start
	    name.start <- paste("start",name.start,sep=".")
            data[,name.start] <- start0
    }  

  

    if ((name.start %in% names(data))) {
      new.start <- data[,name.start]
    } else new.start <- rep(0,n)

    if (any(new.start>= new.time)) cat("any(new.start>= new.time) is TRUE\n"); 

    if ((name.id %in% names(data))) idl <- data[,name.id] else {
	    idl <- 1:n
	    data[,name.id] <- idl 
    }

    splits <- which(new.cuts<new.time & new.start<new.cuts)

    if (length(splits)) {
	    rows  <- c(1:n,splits)
	    new.time <-   c(new.time,new.time[splits])
	    new.start <-  c(new.start,new.cuts[splits])
	    new.status <- c(new.status,new.status[splits])
	    new.ccc <-    c(new.cuts,new.cuts[splits])
	    idl <- c(idl,idl[splits])
	    new.time[splits] <- new.cuts[splits]
	    new.status[splits] <- cens.code
	    data <- data[rows,]
	    data[,time] <- new.time
	    data[,status] <- new.status
	    data[,name.start] <- new.start
	    data[,name.id] <- idl
###    if (num %in% names(data))
###        data[,num] <- data[,num] + new.num else data[,num] <- new.num

    }

    if (time.group) {
      group.time <- paste("before",cutname,sep=".")
      data[,group.time] <- c(rep(1,n),rep(0,length(splits)))
    } 

    if (order.id) data <- data[order(idl,new.start),] 

    return(data)
    ## }}} 
} 

