###{{{ bicomprisk0
bicomprisk0 <- function(formula, data, cause=c(1,1), cens=0, causes, indiv, strata=NULL, id,num,prodlim=FALSE,messages=FALSE,model,...) {
  mycall <- match.call()
  formulaId <- Specials(formula,"id")
  formulaIndiv <- Specials(formula,"indiv")
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- paste("~.-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")",
                     "-indiv(",paste(formulaIndiv,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  indiv <- formulaIndiv
  if (!is.null(formulaIndiv)) {
    mycall$indiv <- indiv
  } 
  if (missing(id)) stop("Missing 'id' variable")
  
  timevar <- terms(formula)[[2]]
  ##  hh <- with(data,eval(timevar))
  if (is.call(timevar)) {
    causes <- timevar[[3]]
    timevar <- timevar[[2]]
  }  
  timevar <- as.character(timevar)
  causes <- as.character(causes)
  
  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq_len(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        mycall$data <- dd[[i]]
        eval(mycall)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("bicomprisk.strata","biprobit.strata")
      res$N <- length(dd)
      return(res)
    }
  }

  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  indiv2 <- covars2 <- NULL 
  
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  which(data[,id]%in%names(idtab==2))
  data <- data[which(data[,id]%in%names(idtab==2)),]
  if (missing(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,num] <- unlist(lapply(idtab,seq_len))
  }

  timevar2 <- paste(timevar,1:2,sep=".")
  causes2 <- paste(causes,1:2,sep=".")
  if (length(covars)>0)
    covars2 <- paste(covars,1,sep=".")
  for (i in seq_len(length(indiv)))
    indiv2 <- c(indiv2, paste(indiv[i],1:2,sep="."))
  
  ww0 <- reshape(data[,c(timevar,causes,covars,indiv,id,num)],
                 direction="wide",idvar=id,timevar=num)[,c(timevar2,causes2,indiv2,covars2)]
 
  switchers <- which(ww0[,timevar2[1]]>ww0[,timevar2[2]])
  switchpos <- 1:4
  if (length(indiv)>0)
    switchpos <- c(switchpos, seq_len(2*length(indiv))+4)
  newpos <- switchpos
  for (i in seq_len(length(newpos)))   
    newpos[i] <- newpos[i] + ifelse(i%%2==1,1,-1)

  ww0[switchers,switchpos] <- ww0[switchers,newpos]
  ww0 <- na.omit(ww0)
 
  status <- rep(0,nrow(ww0))
  time <- ww0[,timevar2[2]]
  mycauses <- setdiff(unique(data[,causes]),0)

  time <- status <- rep(0,nrow(ww0))
  time <- ww0[,"time.1"]
  
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cause[2])
  status[idx2] <- 1
  time[idx2] <- ww0[idx2,timevar2[2]]
  
  idx2 <- which(ww0[,causes2[1]]!=cause[1] & ww0[,causes2[2]]!=cens)
  status[idx2] <- 2
  ##  time[idx2] <- ww0[idx2,"time.1"]
  
  idx1 <- which(ww0[,causes2[1]]==cens)
  ##status[idx1] <- 0
  ##  time[idx1] <- ww0[idx1,"time.1"]
  
  idx1 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cens)
  ##  status[idx1] <- 0
  time[idx1] <- ww0[idx1,timevar2[2]]
  
  idx1 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]!=cause[2] & ww0[,causes2[2]]!=cens)
  status[idx1] <- 2
  time[idx1] <- ww0[idx1,timevar2[2]]

  
  mydata0 <- mydata <- data.frame(time,status,ww0[,covars2],ww0[,indiv2])
  names(mydata) <- c(timevar,causes,covars,indiv2)
  
  if (!prodlim) {
    ff <- paste("Surv(",timevar,",",causes,"!=",cens,") ~ 1",sep="")
    if (length(c(covars,indiv))>0) {
      xx <- c(covars,indiv2)
      for (i in seq_len(length(xx)))
        xx[i] <- paste("const(",xx[i],")",sep="")
      ff <- paste(c(ff,xx),collapse="+")
      if (missing(model)) model <- "fg"      
    }
    if (missing(model)) model <- "additive"
    add<-comp.risk(as.formula(ff),data=mydata,
                   status,causeS=1,n.sim=0,resample.iid=0,model=model)
    padd <- predict(add,X=1,se=0,uniform=0)
  } else {
    ff <- as.formula(paste("Hist(",timevar,",",causes,")~",paste(c("1",covars,indiv2),collapse="+")))
    padd <- prodlim(ff,data=mydata)
  }
  class(padd) <- c("bicomprisk",class(padd))
  return(padd)  
}
###}}} bicomprisk0

bicomprisk <- function(formula, data, cause=c(1,1), cens=0, causes, indiv, strata=NULL, id,num,prodlim=FALSE,messages=TRUE,model,return.data=0,uniform=0,robust=0,...) {
  mycall <- match.call()
  formulaId <- Specials(formula,"id")
  formulaIndiv <- Specials(formula,"indiv")
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- paste("~.-id(",formulaId,")",
                     "-strata(",paste(formulaStrata,collapse="+"),")",
                     "-indiv(",paste(formulaIndiv,collapse="+"),")")
  formula <- update(formula,formulaSt)
  if (!is.null(formulaId)) {
    id <- formulaId
    mycall$id <- id
  }
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  indiv <- formulaIndiv
  if (!is.null(formulaIndiv)) {
    mycall$indiv <- indiv
  } 
  if (missing(id)) stop("Missing 'id' variable")
  
  timevar <- terms(formula)[[2]]
  ##  hh <- with(data,eval(timevar))
  if (is.call(timevar)) {
    causes <- timevar[[3]]
    timevar <- timevar[[2]]
  }  
  timevar <- as.character(timevar)
  causes <- as.character(causes)
  
  if (!is.null(strata)) {
    dd <- split(data,interaction(data[,strata]))
    nn <- unlist(lapply(dd,nrow))
    dd[which(nn==0)] <- NULL
    if (length(dd)>1) {
      fit <- lapply(seq_len(length(dd)),function(i) {
        if (messages>0) message("Strata '",names(dd)[i],"'")
        mycall$data <- dd[[i]]
        eval(mycall)
      })
      res <- list(model=fit)
      res$strata <- names(res$model) <- names(dd)
      class(res) <- c("bicomprisk.strata","biprobit.strata")
      res$N <- length(dd)
      return(res)
    }
  }

  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  indiv2 <- covars2 <- NULL 
  
  data <- data[order(data[,id]),]
  idtab <- table(data[,id])
  which(data[,id]%in%names(idtab==2))
  data <- data[which(data[,id]%in%names(idtab==2)),]
  if (missing(num)) {
    idtab <- table(data[,id])
    num <- "num"
    while (num%in%names(data)) num <- paste(num,"_",sep="")
    data[,num] <- unlist(lapply(idtab,seq_len))
  }

  timevar2 <- paste(timevar,1:2,sep=".")
  causes2 <- paste(causes,1:2,sep=".")
  if (length(covars)>0)
    covars2 <- paste(covars,1,sep=".")
  for (i in seq_len(length(indiv)))
    indiv2 <- c(indiv2, paste(indiv[i],1:2,sep="."))
  
  ww0 <- reshape(data[,c(timevar,causes,covars,indiv,id,num)],
                 direction="wide",idvar=id,timevar=num)[,c(timevar2,causes2,indiv2,covars2)]
 
  ## switchers <- which(ww0[,timevar2[1]]>ww0[,timevar2[2]])
  ## switchpos <- 1:4
  ## if (length(indiv)>0)
  ##   switchpos <- c(switchpos, seq_len(2*length(indiv))+4)
  ## newpos <- switchpos
  ## for (i in seq_len(length(newpos)))   
  ##   newpos[i] <- newpos[i] + ifelse(i%%2==1,1,-1)
  ##  ww0[switchers,switchpos] <- ww0[switchers,newpos]
  ww0 <- na.omit(ww0)
 
  status <- rep(0,nrow(ww0))
  time <- ww0[,timevar2[2]]
  mycauses <- setdiff(unique(data[,causes]),0)

  time <- status <- rep(0,nrow(ww0))
  time <- ww0[,"time.1"]

  ##  suppressMessages(browser())  

  ##(i,j)
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cause[2])
  status[idx2] <- 1
  time[idx2] <- apply(ww0[idx2,timevar2[1:2]],1,max)

  ##(0,0), (0,j)
  idx2 <- which(ww0[,causes2[1]]==cens & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  status[idx2] <- 0
  time[idx2] <- ww0[idx2,timevar2[1]]

  ##(ic,0), (ic,j)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]==cens | ww0[,causes2[2]]==cause[2]))
  status[idx2] <- 2
  time[idx2] <- ww0[idx2,timevar2[1]]

  ##(i,0)
  idx2 <- which(ww0[,causes2[1]]==cause[1] & ww0[,causes2[2]]==cens)
  status[idx2] <- 0
  time[idx2] <- ww0[idx2,timevar2[2]]

  ##(ic,jc)
  idx2 <- which(ww0[,causes2[1]]!=cens & ww0[,causes2[1]]!=cause[1] & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  status[idx2] <- 2
  time[idx2] <- apply(ww0[idx2,timevar2[1:2]],1,min)

  ##(0,jc),(i,jc)
  idx2 <- which((ww0[,causes2[1]]==cens | ww0[,causes2[1]]==cause[1]) & (ww0[,causes2[2]]!=cens & ww0[,causes2[2]]!=cause[2]))
  status[idx2] <- 2
  time[idx2] <- ww0[idx2,timevar2[2]]
  
  mydata0 <- mydata <- data.frame(time,status,ww0[,covars2],ww0[,indiv2])
  names(mydata) <- c(timevar,causes,covars,indiv2)
  
  if (return.data==2) return(list(data=mydata)) else {
  if (!prodlim) {
    ff <- paste("Surv(",timevar,",",causes,"!=",cens,") ~ 1",sep="")
    if (length(c(covars,indiv))>0) {
      xx <- c(covars,indiv2)
      for (i in seq_len(length(xx)))
        xx[i] <- paste("const(",xx[i],")",sep="")
      ff <- paste(c(ff,xx),collapse="+")
      if (missing(model)) model <- "fg"      
    }
    if (missing(model)) model <- "additive"
    add<-comp.risk(as.formula(ff),data=mydata,
                   status,causeS=1,n.sim=0,resample.iid=1,model=model,conservative=robust)
    padd <- predict(add,X=1,se=1,uniform=uniform,resample.iid=1)
  } else {
    ff <- as.formula(paste("Hist(",timevar,",",causes,")~",paste(c("1",covars,indiv2),collapse="+")))
    padd <- prodlim(ff,data=mydata)
  }
###  class(padd) <- c("bicomprisk",class(padd))
 if (return.data==1) return(list(comp.risk=padd,data=mydata)) else return(padd)  
  }
}

plot.bicomprisk <- function(x,add=FALSE,...) {
  if ("predict.timereg"%in%class(x)) {    
    if (!add) { plot.predict.timereg(x,...) }
    else {
      with(x,lines(time,P1,...))
    }
    
  } else {
    plot(x,...)
  }
  return(invisible(x))
}

