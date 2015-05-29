
##' @export
force.same.cens <- function(data,id="id", time="time",cause="cause",entrytime=NULL,cens.code=0)
{ ## {{{ 
  ### no missing values

   ### handle time-variables separately
   data <- dsort(data,id)
   w <- which(names(data) %in% c(time,cause,entrytime,id) )
   datao <- data[,-w]

   data <- data[,c(time,cause,entrytime,id)]

  if (is.null(entrytime)) entry <- rep(0,nrow(data)) else entry <- data[,entrytime]

   censo <- (data[,cause]==cens.code)
   Wide <- fast.reshape(data,id=id)
   time1 <- paste(time,1,sep="")
   time2 <- paste(time,2,sep="")
   stat1 <- paste(cause,1,sep="")
   stat2 <- paste(cause,2,sep="")

   ### enforce same censoring ## {{{ 

   mintime <- pmin(Wide[,time1],Wide[,time2])
   statmin <- ifelse(Wide[,time1]<Wide[,time2],Wide[,stat1],Wide[,stat2])
###
   cens.first <- which(statmin==cens.code)
   Wide[cens.first,time1] <- mintime[cens.first]
   Wide[cens.first,time2] <- mintime[cens.first]
   Wide[cens.first,stat1] <- cens.code
   Wide[cens.first,stat2] <- cens.code

###   idn <- table(data[,id]); rep(mintime,idn)
   
   long <-   fast.reshape(Wide)
   long <- long[!is.na(long[,time]) & !is.na(long[,cause]),]

   ## }}} 

   if (!is.null(entrytime)) { ## {{{ enforce same truncation
      stop("not implemented yet\n")
   entry1 <- paste(entrytime,1,sep="")
   entry2 <- paste(entrytime,2,sep="")
###   trunc.first<- which(statmin==cens.code)
   trunc.max <- pmax(Wide[,entry1],Wide[,entry2])
   Wide[,entry1] <- Wide[,entry2] <- trunc.max
### drop those that enter later
  enter.after <- ( Wide[,entry1] < Wide[,time1]) & (Wide[,entry2] < Wide[,time2])
  Wide <- Wide[enter.after,] 
  } ## }}} 

 long <- cbind(long,datao)

 return(long)
} ## }}} 


