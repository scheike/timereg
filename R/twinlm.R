##' Fits a classical twin model for quantitative traits.
##'
##' @title Classic twin model for quantitative traits
##' @return   Returns an object of class \code{twinlm}.
##' @author Klaus K. Holst
##' @seealso \code{\link{bptwin}}, \code{\link{twinsim}}
##' @export
##' @examples
##' ## Simulate data
##' set.seed(1)
##' d <- twinsim(1000,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
##' ## E(y|z1,z2) = z1 - z2. var(A) = var(C) = var(E) = 1
##' 
##' ## E.g to fit the data to an ACE-model without any confounders we simply write
##' ace <- twinlm(y ~ 1, data=d, DZ="DZ", zyg="zyg", id="id")
##' ace
##' ## An AE-model could be fitted as
##' ae <- twinlm(y ~ 1, data=d, DZ="DZ", zyg="zyg", id="id", type="ae")
##' ## LRT:
##' compare(ae,ace)
##' ## AIC
##' AIC(ae)-AIC(ace)
##' ## To adjust for the covariates we simply alter the formula statement
##' ace2 <- twinlm(y ~ x11+x12, data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ## Summary/GOF
##' summary(ace2)
##' ## An interaction could be analyzed as:
##' ace3 <- twinlm(y ~ x11+x12 + x11:I(x12<0), data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ## Categorical variables are also supported
##' d2 <- transform(d,x12cat=cut(x12,3,labels=c("Low","Med","High")))
##' ace4 <- twinlm(y ~ x11+x12cat, data=d2, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ## plot the model structure
##' \dontrun{
##' plot(ace4)
##' }
##' @keywords models
##' @keywords regression
##' @param formula Formula specifying effects of covariates on the response.
##' @param data \code{data.frame} with one observation pr row. In
##'     addition a column with the zygosity (DZ or MZ given as a factor) of
##'     each individual much be
##'     specified as well as a twin id variable giving a unique pair of
##'     numbers/factors to each twin pair.
##' @param id The name of the column in the dataset containing the twin-id variable.
##' @param zyg The name of the column in the dataset containing the
##'     zygosity variable.
##' @param DZ Character defining the level in the zyg variable
##'     corresponding to the dyzogitic twins. If this argument is missing,
##'     the reference level (i.e. the first level) will be interpreted as
##'     the dyzogitic twins.
##' @param OS Optional. Character defining the level in the zyg variable
##'     corresponding to the oppposite sex dyzogitic twins.
##' @param weight Weight matrix if needed by the chosen estimator. For use
##'     with Inverse Probability Weights
##' @param type Character defining the type of analysis to be
##'     performed. Should be a subset of "aced" (additive genetic factors, common
##'     environmental factors, unique environmental factors, dominant
##'     genetic factors).
##' @param twinnum The name of the column in the dataset numbering the
##'     twins (1,2). If it does not exist in \code{data} it will
##'     automatically be created.
##' @param binary If \code{TRUE} a liability model is fitted. Note that if the right-hand-side of the formula is a factor, character vector, og logical variable, then the liability model is automatically chosen (wrapper of the \code{bptwin} function).
##' @param keep Vector of variables from \code{data} that are not
##'     specified in \code{formula}, to be added to data.frame of the SEM
##' @param estimator Choice of estimator/model
##' @param control Control argument parsed on to the optimization routine
##' @param constrain Development argument
##' @param ... Additional arguments parsed on to lower-level functions
twinlm <- function(formula, data, id, zyg, DZ, OS, weight=NULL, type=c("ace"), twinnum="twinnum", binary=FALSE,keep=weight,estimator="gaussian",constrain=TRUE,control=list(),...) {

  cl <- match.call(expand.dots=TRUE)
  opt <- options(na.action="na.pass")
  mf <- model.frame(formula,data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  formula <- update(formula, ~ . + 1)
  yvar <- getoutcome(formula)
  if (missing(zyg)) stop("Zygosity variable not specified")
  if (missing(id)) stop("Twin-pair variable not specified")

  if (binary | is.factor(data[,yvar]) | is.character(data[,yvar]) | is.logical(data[,yvar])) {
    args <- as.list(cl)
    args[[1]] <- NULL
    return(do.call("bptwin",args))
  }
  
  type <- tolower(type)
  ## if ("u" %in% type) type <- c("ue")
  
  varnames <- all.vars(formula)
  latentnames <- c("a1","a2","c1","c2","d1","d2","e1","e2")
  if (any(latentnames%in%varnames))
    stop(paste(paste(latentnames,collapse=",")," reserved for names of latent variables.",sep=""))
  
  mm <- model.matrix(formula,mf)
  options(opt)
  
  covars <- colnames(mm)
  if (attr(terms(formula),"intercept")==1)
    covars <- covars[-1]
  if(length(covars)<1) covars <- NULL

  zygstat <- data[,zyg]
  if(!is.factor(zygstat)) {
    zygstat <- as.factor(zygstat)
  }
  zyglev <- levels(zygstat)
  if (length(zyglev)>2 & missing(OS)) stop("More than two zygosity levels and no opposite sex (OS) group specified")
  
  ## Get data on wide format and divide into two groups by zygosity
  if (!twinnum%in%names(data)) {
    mynum <- rep(1,nrow(data))
    mynum[zygstat==zyglev[1]][duplicated(data[zygstat==zyglev[1],id])] <- 2
    mynum[zygstat==zyglev[2]][duplicated(data[zygstat==zyglev[2],id])] <- 2
    if (!missing(OS)) {
      mynum[zygstat==zyglev[3]][duplicated(data[zygstat==zyglev[3],id])] <- 2
    }
    data[,twinnum] <- mynum
  }
  
  cur <- cbind(data[,c(yvar,keep)],as.data.frame(cbind(as.numeric(data[,twinnum]),as.numeric(data[,id]),as.numeric(zygstat))));
  colnames(cur) <- c(yvar,keep,twinnum,id,zyg)
  mydata <- cbind(cur,mm)
  
  if (missing(DZ)) {
    warning("Using first level, `",zyglev[1],"', in status variable as indicator for 'dizygotic'", sep="")
    DZ <- zyglev[1]
  }
  myDZ <- which(levels(zygstat)==DZ)
  myOS <- NULL
  if (!missing(OS)) myOS <- which(levels(zygstat)==OS)
  numlev <- seq(2+!missing(OS))
  myMZ <- setdiff(numlev,c(myDZ,myOS))

 
  data1 <- mydata[mydata[,zyg]==myMZ,,drop=FALSE]
  if (nrow(data1)==0) stop("No MZ twins found")
  data2 <- mydata[mydata[,zyg]==myDZ,,drop=FALSE]
  if (nrow(data2)==0) stop("No DZ twins found")
  data3 <- NULL
  if (!missing(OS)) {
    data3 <- mydata[mydata[,zyg]==myOS,,drop=FALSE]
    if (nrow(data3)==0) stop("No opposite-sex twins found")
  }
  
  ## Monozygotic
  data1.1 <- data1[which(data1[,twinnum]==1),c(id,zyg,yvar,keep,covars)]; colnames(data1.1) <- c(id,zyg,paste(colnames(data1.1)[-c(1,2)],".1",sep=""))  
  data1.2 <- data1[which(data1[,twinnum]==2),c(id,yvar,keep,covars),drop=FALSE]; colnames(data1.2) <- c(id, paste(colnames(data1.2)[-1],".2",sep=""))

  ##Missing data?
  id1 <- data1.1[,id]
  id2 <- data1.2[,id]
  d1.mis <- setdiff(id2,id1) # Id's not in data1.1
  d2.mis <- setdiff(id1,id2) # Id's not in data1.2
  if (length(d1.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data1.1), nrow=length(d1.mis))
    d.temp[,1] <- d1.mis; colnames(d.temp) <- colnames(data1.1)
    data1.1 <- rbind(data1.1,d.temp)
  }
  if (length(d2.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data1.2), nrow=length(d2.mis))
    d.temp[,1] <- d2.mis; colnames(d.temp) <- colnames(data1.2)
    data1.2 <- rbind(data1.2,d.temp)
  } 


  ## Dizygotic 
  data2.1 <- data2[data2[,twinnum]==1,c(id,zyg,yvar,keep,covars)]; colnames(data2.1) <- c(id,zyg,paste(colnames(data2.1)[-c(1,2)],".1",sep=""))
  data2.2 <- data2[data2[,twinnum]==2,c(id,yvar,keep,covars),drop=FALSE]; colnames(data2.2) <- c(id, paste(colnames(data2.2)[-1],".2",sep=""))

  id1 <- data2.1[,id]
  id2 <- data2.2[,id]
  d1.mis <- setdiff(id2,id1) # Id's not in data1.1
  d2.mis <- setdiff(id1,id2) # Id's not in data1.2
  if (length(d1.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data2.1), nrow=length(d1.mis))
    d.temp[,1] <- d1.mis; colnames(d.temp) <- colnames(data2.1)
    data2.1 <- rbind(data2.1,d.temp)
  }
  if (length(d2.mis)>0) {
    d.temp <- matrix(NA,ncol= ncol(data2.2), nrow=length(d2.mis))
    d.temp[,1] <- d2.mis; colnames(d.temp) <- colnames(data2.2)
    data2.2 <- rbind(data2.2,d.temp)
  }

  wide1 <- merge(x=data1.1,y=data1.2, by=id); wide1[,zyg] <- myMZ
  wide2 <- merge(x=data2.1,y=data2.2, by=id); wide2[,zyg] <- myDZ
  wide3 <- NULL

  if (!missing(OS)) {
    ## Opposite sex
    data3.1 <- data3[which(data3[,twinnum]==1),c(id,zyg,yvar,keep,covars)]; colnames(data3.1) <- c(id,zyg,paste(colnames(data3.1)[-c(1,2)],".1",sep=""))
    data3.2 <- data3[which(data3[,twinnum]==2),c(id,yvar,keep,covars),drop=FALSE]; colnames(data3.2) <- c(id, paste(colnames(data3.2)[-1],".2",sep=""))

    id1 <- data3.1[,id]
    id2 <- data3.2[,id]
    d1.mis <- setdiff(id2,id1) # Id's not in data3.1
    d2.mis <- setdiff(id1,id2) # Id's not in data3.2
    if (length(d1.mis)>0) {
      d.temp <- matrix(NA,ncol= ncol(data3.1), nrow=length(d1.mis))
      d.temp[,1] <- d1.mis; colnames(d.temp) <- colnames(data3.1)
      data3.1 <- rbind(data3.1,d.temp)
    }
    if (length(d2.mis)>0) {
      d.temp <- matrix(NA,ncol= ncol(data3.2), nrow=length(d2.mis))
      d.temp[,1] <- d2.mis; colnames(d.temp) <- colnames(data3.2)
      datta3.2 <- rbind(data3.2,d.temp)
    }
    wide3 <- merge(x=data3.1,y=data3.2, by=id); wide3[,zyg] <- myOS

  }

  ## ###### The SEM
  outcomes <- paste(yvar,".",1:2,sep="")
  model1<-lvm(outcomes,silent=TRUE)
  f1 <- as.formula(paste(outcomes[1]," ~ ."))
  f2 <- as.formula(paste(outcomes[2]," ~ ."))
  regression(model1,silent=TRUE) <- update(f1, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e1,lambda[e])) ## + f(u,sdu1))
  regression(model1,silent=TRUE) <- update(f2, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e2,lambda[e]))## + f(u,sdu1))
  latent(model1) <- ~ a1+c1+d1+e1+e2##+u
  intercept(model1,latent(model1)) <- 0
  if (!is.null(covars))
    for (i in 1:length(covars)) {
      lava:::regfix(model1, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta[",i,"]",sep="")
      lava:::regfix(model1, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta[",i,"]",sep="")
    }
  covariance(model1) <- update(f1, . ~  v(0))
  covariance(model1) <- update(f2, . ~  v(0))
  lava:::covfix(model1, latent(model1), var2=NULL) <- 1
  if (!type%in%c("sat","flex")) {    
    intercept(model1,outcomes) <- "mu"
  }
  if (type%in%c("u","flex","sat")) {
    kill(model1) <- ~e1+e2
    covariance(model1,outcomes) <- "v1"
  }
  
  model2 <- model1
  cancel(model2) <- update(f2, . ~ a1)
  cancel(model2) <- update(f2, . ~ d1)
  regression(model2,silent=TRUE) <- update(f2, . ~ f(a2,lambda[a]))
  regression(model2,silent=TRUE) <- update(f2, . ~ f(d2,lambda[d]))
  
  covariance(model2) <- a1 ~ f(a2,0.5)
  covariance(model2) <- d1 ~ f(d2,0.25)
  latent(model2) <- ~ a2+d2
  intercept(model2, ~ a2+d2) <- 0
  covariance(model2) <- c(a2,d2) ~ v(1)
  model3 <- model2
  covariance(model3) <- a1 ~ f(a2,r1)
  suppressMessages(constrain(model3, r1 ~ ra) <- lava:::Range.lvm())
  covariance(model3) <- d1 ~ f(d2,r2)
  suppressMessages(constrain(model3, r2 ~ rd) <- lava:::Range.lvm())

  if (type=="flex") {
     intercept(model1,outcomes) <- "mu1"
     intercept(model2,outcomes) <- "mu2"
     intercept(model3,outcomes) <- "mu3"
     covariance(model1,outcomes) <- "var(MZ)"
     covariance(model2,outcomes) <- "var(DZ)"
     covariance(model3,outcomes) <- "var(OS)"
   }
  if (type=="sat") {
     covariance(model1,outcomes) <- c("var(MZ)1","var(MZ)2")
     covariance(model2,outcomes) <- c("var(DZ)1","var(DZ)2")
     covariance(model3,outcomes) <- c("var(OS)1","var(OS)2")
  }
  if (type%in%c("u","flex","sat")) {
    if (constrain) {
      if (type=="sat") {
        model1 <- covariance(model1,outcomes,constrain=TRUE,rname="atanh(rhoMZ)",cname="covMZ",lname="log(var(MZ)).1",l2name="log(var(MZ)).2")
        model2 <- covariance(model2,outcomes,constrain=TRUE,rname="atanh(rhoDZ)",cname="covDZ",lname="log(var(DZ)).1",l2name="log(var(DZ)).2")
        model3 <- covariance(model3,outcomes,constrain=TRUE,rname="atanh(rhoOS)",cname="covOS",lname="log(var(OS)).1",l2name="log(var(OS)).2")
      } else {
        if (type=="flex") {
          model1 <- covariance(model1,outcomes,constrain=TRUE,rname="atanh(rhoMZ)",cname="covMZ",lname="log(var(MZ))")
          model2 <- covariance(model2,outcomes,constrain=TRUE,rname="atanh(rhoDZ)",cname="covDZ",lname="log(var(DZ))")
          model3 <- covariance(model3,outcomes,constrain=TRUE,rname="atanh(rhoOS)",cname="covOS",lname="log(var(OS))")
        }  else {
          model1 <- covariance(model1,outcomes,constrain=TRUE,rname="atanh(rhoMZ)",cname="covMZ",lname="log(var)")
          model2 <- covariance(model2,outcomes,constrain=TRUE,rname="atanh(rhoDZ)",cname="covDZ",lname="log(var)")
          model3 <- covariance(model3,outcomes,constrain=TRUE,rname="atanh(rhoOS)",cname="covOS",lname="log(var)")          
        }        
      }     
    } else {
      covariance(model1,outcomes[1],outcomes[2]) <- "covMZ"
      covariance(model2,outcomes[1],outcomes[2]) <- "covDZ"
      covariance(model3,outcomes[1],outcomes[2]) <- "covOS"
    }
  }
  if (!is.null(covars) & type%in%c("flex","sat")) {
    sta <- ""
    if (type=="sat") sta <- "b"
       for (i in 1:length(covars)) {
         lava:::regfix(model1, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta1[",i,"]",sep="")         
         lava:::regfix(model1, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta1",sta,"[",i,"]",sep="")
         lava:::regfix(model2, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta2[",i,"]",sep="")
         lava:::regfix(model2, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta2",sta,"[",i,"]",sep="")
         lava:::regfix(model3, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta3[",i,"]",sep="")
         lava:::regfix(model3, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta3",sta,"[",i,"]",sep="")
       }
  }

  
  full <- list(MZ=model1,DZ=model2,OS=model3)
  isA <- length(grep("a",type))>0 & type!="sat"
  isC <- length(grep("c",type))>0
  isD <- length(grep("d",type))>0
  isE <- length(grep("e",type))>0 | type=="sat" | type=="u"
  if (!isA) {
    kill(model1) <- ~ a1 + a2
    kill(model2) <- ~ a1 + a2
    kill(model3) <- ~ a1 + a2 + ra
    constrain(model3,r1~1) <- NULL
  }
  if (!isD) {
    kill(model1) <- ~ d1 + d2
    kill(model2) <- ~ d1 + d2
    kill(model3) <- ~ d1 + d2 + rd
    constrain(model3,r2~1) <- NULL
  }
  if (!isC) {
    kill(model1) <- ~ c1 + c2
    kill(model2) <- ~ c1 + c2
    kill(model3) <- ~ c1 + c2
  }
  if (!isE) {
    kill(model1) <- ~ e1 + e2
    kill(model2) <- ~ e1 + e2
    kill(model3) <- ~ e1 + e2
  }

  ## Full rank covariate/design matrix?
  for (i in covars) {
    myvars <- paste(i,c(1,2),sep=".")
    dif <- wide1[,myvars[1]]-wide1[,myvars[2]]   
    mykeep <- myvars
    if (all(na.omit(dif)==00)) {
      mykeep <- mykeep[-2]
    }   
    trash <- setdiff(myvars,mykeep)
    if (length(mykeep)==1) {
      regression(model1, to=outcomes[2], from=mykeep) <- lava:::regfix(model1)$label[trash,outcomes[2]]
      kill(model1) <- trash
    }

    dif <- wide2[,myvars[1]]-wide2[,myvars[2]]   
    mykeep <- myvars
    if (all(na.omit(dif)==00)) {
      mykeep <- mykeep[-2]
    }  
    trash <- setdiff(myvars,mykeep)
    if (length(mykeep)==1) {
      regression(model2, to=outcomes[2], from=mykeep) <- lava:::regfix(model2)$label[trash,outcomes[2]]
      kill(model2) <- trash
    }

    if (!missing(OS)) {
      dif <- wide3[,myvars[1]]-wide2[,myvars[2]]   
      mykeep <- myvars
      if (all(na.omit(dif)==00)) {
        mykeep <- mykeep[-2]
      }  
      trash <- setdiff(myvars,mykeep)
      if (length(mykeep)==1) {
        regression(model3, to=outcomes[2], from=mykeep) <- lava:::regfix(model3)$label[trash,outcomes[2]]
        kill(model3) <- trash
      }      
    }
  }
  if (!is.null(weight)) {
    weight <- paste(weight,1:2,sep=".")
    estimator <- "weighted"
  }

  newkeep <- unlist(sapply(keep, function(x) paste(x, 1:2, 
                                                   sep = ".")))
  mm <- list(MZ = model1, DZ = model2)
  dd <- list(wide1, wide2)
  if (!missing(OS)) {
    mm <- c(mm, OS = list(model3))
    dd <- c(dd, list(wide3))
  }

  names(dd) <- names(mm)

  ## suppressWarnings(mg <- multigroup(mm, dd, missing=TRUE,fix=FALSE,keep=newkeep,type=2))

  if (is.null(estimator)) return(multigroup(mm, dd, missing=TRUE,fix=FALSE,keep=newkeep,type=2))

  optim <- list(method="nlminb2",refit=FALSE,gamma=1,start=rep(0.1,length(coef(mm[[1]]))*length(mm)))
  ##                                                       with(mg,npar+npar.mean)))
  if (length(control)>0) {
    optim[names(control)] <- control
  }

  if (is.Surv(data[,yvar])) {
    require("lava.tobit")
    if (is.null(optim$method))
       optim$method <- "nlminb1"
    e <- estimate(mm,dd,control=optim,...)
  } else {
    e <- estimate(mm,dd,weight=weight,estimator=estimator,fix=FALSE,control=optim,...)
  }
  if (!is.null(optim$refit) && optim$refit) {
    optim$method <- "NR"
    optim$start <- pars(e)
    if (is.Surv(data[,yvar])) {
      e <- estimate(mm,dd,estimator=estimator,fix=FALSE,control=optim,...)      
    } else {
      e <- estimate(mm,dd,weight=weight,estimator=estimator,fix=FALSE,control=optim,...)
    }
  }
  res <- list(coefficients=e$opt$estimate, vcov=Inverse(information(e)), estimate=e, model=mm, full=full, call=cl, data=data, zyg=zyg, id=id, twinnum=twinnum, type=type, model.mz=model1, model.dz=model2, model.dzos=model3, data.mz=wide1, data.dz=wide2, data.os=wide3, OS=!missing(OS), constrain=constrain)
  class(res) <- "twinlm"
  return(res)
}
