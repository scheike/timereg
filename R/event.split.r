#' EventSplit (SurvSplit).
#' 
#' contstructs start stop formulation of event time data after a variable in
#' the data.set. Similar to SurvSplit of the survival package but can also
#' split after random time given in data frame.
#' 
#' 
#' @param data data to be split
#' @param time time variable.
#' @param status status variable.
#' @param cuts cuts variable or numeric cut (only one value)
#' @param name.id name of id variable.
#' @param name.start name of start variable in data, start can also be numeric "0"
#' @param cens.code code for the censoring.
#' @param order.id order data after id and start.
#' @param time.group make variable "before"."cut" that keeps track of wether start,stop is before (1) or after cut (0).
#' @author Thomas Scheike
#' @keywords survival
#' @examples
#' 
#' set.seed(1)
#' d <- data.frame(event=round(5*runif(5),2),start=1:5,time=2*1:5,
#' 		status=rbinom(5,1,0.5),x=1:5)
#' d
#' 
#' d0 <- event.split(d,cuts="event",name.start=0)
#' d0
#' 
#' dd <- event.split(d,cuts="event")
#' dd
#' ddd <- event.split(dd,cuts=3.5)
#' ddd
#' event.split(ddd,cuts=5.5)
#' 
#' ### successive cutting for many values 
#' dd <- d
#' for  (cuts in seq(2,3,by=0.3)) dd <- event.split(dd,cuts=cuts)
#' dd
#' 
#' ###########################################################################
#' ### same but for situation with multiple events along the time-axis
#' ###########################################################################
#' d <- data.frame(event1=1:5+runif(5)*0.5,start=1:5,time=2*1:5,
#' 		status=rbinom(5,1,0.5),x=1:5,start0=0)
#' d$event2 <- d$event1+0.2
#' d$event2[4:5] <- NA 
#' d
#' 
#' d0 <- event.split(d,cuts="event1",name.start="start",time="time",status="status")
#' d0
#' ###
#' d00 <- event.split(d0,cuts="event2",name.start="start",time="time",status="status")
#' d00
#' 
#' @export
event.split <- function(data,
		time="time",status="status",cuts="cuts",name.id="id",
		name.start="start", cens.code=0,order.id=TRUE, time.group=FALSE)
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

###    if (newrow)  new.row <- rep(0,nrow(data))

    splits <- which(new.cuts<new.time & new.start<new.cuts)

    if (length(splits)) {
	    rows  <- c(1:n,splits)
	    new.time <-   c(new.time,new.time[splits])
	    new.start <-  c(new.start,new.cuts[splits])
	    new.status <- c(new.status,new.status[splits])
	    new.ccc <-    c(new.cuts,new.cuts[splits])
###	    new.row <- c(new.row,rep(1,length(splits)))
	    idl <- c(idl,idl[splits])
	    new.time[splits] <- new.cuts[splits]
	    new.status[splits] <- cens.code
	    data <- data[rows,]
	    data[,time] <- new.time
	    data[,status] <- new.status
	    data[,name.start] <- new.start
	    data[,name.id] <- idl
###	    if (newrow) data[,newrow.name]  <-  new.row 
###    if (num %in% names(data))
###        data[,num] <- data[,num] + new.num else data[,num] <- new.num

    }

    if (time.group) {
      group.time <- paste("before",cutname,sep=".")
      data[,group.time] <- 1*(data[,name.start]<data[,cutname]) ## sc(rep(1,n),rep(0,length(splits)))
    } 

    if (order.id) data <- data[order(idl,new.start),] 

    return(data)
    ## }}} 
} 


