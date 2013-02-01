##' @export
surv.boxarea <- function(left.trunc,right.cens,data,timevar="time",status="status",id="id",covars=NULL,num=NULL,silent=1,boxtimevar="boxtime")
{ ## {{{

  if (is.null(data[,id])) stop("Wrong cluster variable")
  if (is.null(data[,timevar])) stop("Wrong time variable")
  if (is.null(data[,status])) stop("Wrong status variable")
  data <- data[order(data[,id]),]
###  num <- NULL
  if (is.null(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,c(num)] <- unlist(lapply(idtab,seq_len))
  }

  timevar2 <- paste(timevar,1:2,sep="")
  status2 <- paste(status,1:2,sep="")
  num2 <- paste(num,1:2,sep="")
  covars2 <- NULL; 
  if (length(covars)>0) covars2 <- paste(covars,1:2,sep="")

  ww0 <- fast.reshape(data[,c(timevar,status,covars,id,num)],id=id,num=num)[,c(timevar2,status2,covars2,id)] 
  mleft <-  (ww0[,timevar2[1]]>left.trunc[1]) & (ww0[,timevar2[2]]>left.trunc[2])  ## Both not-truncated

  if (length(na.idx <- which(is.na(mleft)))>0) {
    ##    warning("Removing incomplete cases", na.idx)
    mleft <- mleft[-na.idx]
    ww0 <- ww0[-na.idx,,drop=FALSE]
  }
  if (sum(mleft)==0) stop("No data selected\n"); 
  ww0 <- ww0[which(mleft),,drop=FALSE]

  right1 <- which(ww0[,timevar2[1]] > right.cens[1])
  right2 <- which(ww0[,timevar2[2]] > right.cens[2])
  ww0[,timevar2[1]][right1] <- right.cens[1]
  ww0[,timevar2[2]][right2] <- right.cens[2]
  ww0[,status2[1]][right1] <- 0
  ww0[,status2[2]][right2] <- 0
  truncvar2 <- c("left1","left2")
  ww0 <- cbind(ww0,left.trunc[1])
  ww0 <- cbind(ww0,left.trunc[2])
  colnames(ww0)[c(-1,0) + ncol(ww0)] <- truncvar2

  browser()

  if (silent<=0) message(paste("  Number of joint events:",sum(apply(ww0[,status2],1,sum)==2),"of ",nrow(ww0)),"\n");
  varying <- c(list(timevar2),list(status2),list(truncvar2),lapply(covars,function(x) paste(x,1:2,sep="")))
  lr.data <- data.frame(fast.reshape(ww0,var=c(timevar,status,"left")))
  ### ,v.names=c(timevar,status,"left",covars))
  lr.data[,boxtimevar] <- lr.data[,timevar]-lr.data[,"left"]
  return(structure(lr.data,num=num,time=boxtimevar,status=status,covars=covars,id=id))
}
