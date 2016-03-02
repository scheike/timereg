
event.split <- function(data,
			time="time",status="status",cuts="cuts",
			name.start="start",id=NULL,
			name.id="id",
			entry=NULL,cens.code=0,num="num")
{
## {{{ 
    n <- nrow(data)
    new.time <- data[,time]
    new.status <- data[,status]
    new.cuts <- data[,cuts]
    if (!is.null(entry)) {
      new.start <- data[,entry]
    } else new.start <- rep(0,n)
    idl <- 1:n
    if (any(new.start>= new.time)) cat("any(new.start>= new.time) is TRUE\n"); 

   if (name.id %in% names(data)) 
	   name.id <-paste(id,".split",sep="")

###    splits <- (new.cuts<new.time & new.start<new.cuts)
###    print(sum(splits))
    splits <- which(new.cuts<new.time & new.start<new.cuts)
    idl <- c(idl,idl[splits])
    new.time <- c(new.time,new.time[splits])
    new.start <- c(new.start,new.cuts[splits])
    new.status <- c(new.status,new.status[splits])
    new.ccc <- c(new.cuts,new.cuts[splits])
    ###
    new.time[splits] <- new.cuts[splits]
    new.status[splits] <- cens.code

    new.num <- c(rep(0,n),rep(1,length(splits)))
    new.num <- (new.start>=new.ccc)*1 
    data <- data[idl,]
    data[,time] <- new.time
    data[,status] <- new.status
    data[,name.start] <- new.start
    data[,name.id] <- idl
    data[,num] <- new.num
### data[,"ccc"] <- new.ccc

    if (is.null(id)) data <- data[order(idl),] else data <- data[order(data[,id],new.start),]

    return(data)
    ## }}} 
} 

###d <- data.frame(event=round(5*runif(5),2),
###		start=round(5*runif(5),2),
###		time=1:5,
###		status=rbinom(5,1,0.5),
###		x=1:5)
###d$fcut <- 2
###d
###event.split(d,cuts="event",num="nume")
###event.split(d,cuts="fcut",num="numf")
###
###### successive splitting 
###de <- event.split(d,cuts="event",num="nume")
###de
###def <- event.split(de,cuts="fcut",entry="start",id="id",num="numf")
######
###### successive splitting 
###df <- event.split(d,cuts="fcut",num="numf")
###df
###dfe <- event.split(df,cuts="event",entry="start",id="id",num="nume")
###dfe <- dfe[,names(def)]
###dfe
###dfe==def


