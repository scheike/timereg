##' Fits a classical twin model for quantitative traits.
##'
##' @title Classic twin model for quantitative traits
##' @return   Returns an object of class \code{twinlm}.
##' @author Klaus K. Holst
##' @seealso \code{\link{twinsim}}
##' @export
##' @examples
##' ## Simulate data
##' set.seed(1)
##' d <- twinsim(2000,b1=c(1,-1),b2=c(),acde=c(1,1,0,1))
##' ## E(y|z1,z2) = z1 - z2. var(A) = var(C) = var(E) = 1
##' 
##' ## E.g to fit the data to an ACE-model without any confounders we simply write
##' ace <- twinlm(y1 ~ 1, data=d, DZ="DZ", zyg="zyg", id="id")
##' ace
##' ## An AE-model could be fitted as
##' ae <- twinlm(y1 ~ 1, data=d, DZ="DZ", zyg="zyg", id="id", type="ae")
##' ## LRT:
##' compare(ae,ace)
##' ## AIC
##' AIC(ae)-AIC(ace)
##' ## To adjust for the covariates we simply alter the formula statement
##' ace2 <- twinlm(y1 ~ x11+x12, data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ## Summary/GOF
##' summary(ace2)
##' ## An interaction could be analyzed as:
##' ace3 <- twinlm(y1 ~ x11+x12 + x11:I(x12<0), data=d, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ## Categorical variables are also supported
##' d2 <- transform(d,x12cat=cut(x12,3,labels=c("Low","Med","High")))
##' ace4 <- twinlm(y1 ~ x11+x12cat, data=d2, DZ="DZ", zyg="zyg", id="id", type="ace")
##' ## plot the model structure
##' \donttest{
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
##' @param DZos Optional. Character defining the level in the zyg variable
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
##' @param binary If \code{TRUE} a liability model is fitted
##' @param keep Vector of variables from \code{data} that are not
##'     specified in \code{formula}, to be added to data.frame of the SEM
##' @param estimator Choice of estimator/model.
##' @param ... Additional arguments parsed on to lower-level functions
twinlm <- function(formula, data, id, zyg, DZ, DZos, weight=NULL, type=c("ace"), twinnum="twinnum", binary=FALSE,keep=weight,estimator="gaussian",...) {

  cl <- match.call(expand.dots=TRUE)
  mf <- model.frame(formula,data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  formula <- update(formula, ~ . + 1)
  yvar <- getoutcome(formula)

  if (binary | is.factor(data[,yvar]) | is.character(data[,yvar]) | is.logical(data[,yvar])) {
    args <- as.list(cl)
    args[[1]] <- NULL
    return(do.call("bptwin",args))
  }
  
  type <- tolower(type)
  if ("u" %in% type) type <- c("ue")
  
  varnames <- all.vars(formula)
  latentnames <- c("a1","a2","c1","c2","d1","d2","e1","e2")
  if (any(latentnames%in%varnames))
    stop(paste(paste(latentnames,collapse=",")," reserved for names of latent variables.",sep=""))
  
  
  opt <- options(na.action="na.pass")
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
  if (length(zyglev)>2) stop("Only support for two zygosity levels")
  ## Get data on wide format and divide into two groups by zygosity
  if (!twinnum%in%names(data)) {
    mynum <- rep(1,nrow(data))
    mynum[zygstat==zyglev[1]][duplicated(data[zygstat==zyglev[1],id])] <- 2
    mynum[zygstat==zyglev[2]][duplicated(data[zygstat==zyglev[2],id])] <- 2
    data[,twinnum] <- mynum
  }
  
  cur <- cbind(data[,c(yvar,keep)],as.numeric(data[,twinnum]),as.numeric(data[,id]),as.numeric(zygstat));
  colnames(cur) <- c(yvar,keep,twinnum,id,zyg)
  mydata <- cbind(cur,mm)
  if (missing(DZ)) {
    warning("Using first level, `",zyglev[1],"', in status variable as indicator for 'dizygotic'", sep="")
    DZ <- zyglev[1]
  }
  myDZ <- which(levels(zygstat)==DZ)
  myMZ <- setdiff(1:2,myDZ)

  data1 <- mydata[mydata[,zyg]==myMZ,,drop=FALSE]
  data2 <- mydata[mydata[,zyg]==myDZ,,drop=FALSE]
    
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
    
  ## ###### The SEM
  outcomes <- paste(yvar,".",1:2,sep="")
  model1<-lvm(outcomes,silent=TRUE)
  f1 <- as.formula(paste(outcomes[1]," ~ ."))
  f2 <- as.formula(paste(outcomes[2]," ~ ."))
##  parameter(model1) <- ~sdu1+sdu2
  regression(model1,silent=TRUE) <- update(f1, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e1,lambda[e]) + f(u,sdu1))
  regression(model1,silent=TRUE) <- update(f2, . ~ f(a1,lambda[a])+f(c1,lambda[c])+f(d1,lambda[d]) + f(e2,lambda[e]) + f(u,sdu1))
  latent(model1) <- ~ a1+c1+d1+e1+e2+u
  intercept(model1,latent(model1)) <- 0
  if (!is.null(covars))
    for (i in 1:length(covars)) {
      regfix(model1, from=paste(covars[i],".1",sep=""), to=outcomes[1],silent=TRUE) <- paste("beta[",i,"]",sep="")
      regfix(model1, from=paste(covars[i],".2",sep=""), to=outcomes[2],silent=TRUE) <- paste("beta[",i,"]",sep="")
    }
  covariance(model1) <- update(f1, . ~  v(0))
  covariance(model1) <- update(f2, . ~  v(0))
  covfix(model1, latent(model1), var2=NULL) <- 1
  intfix(model1,outcomes) <- "mu1"
  model2 <- model1
##  parameter(model2) <- ~sdu1+sdu2
  regression(model2) <- update(f1,.~f(u,sdu2))
  regression(model2) <- update(f2,.~f(u,sdu2))
  cancel(model2) <- update(f2, . ~ a1)
  cancel(model2) <- update(f2, . ~ d1)
  regression(model2,silent=TRUE) <- update(f2, . ~ f(a2,lambda[a]))
  regression(model2,silent=TRUE) <- update(f2, . ~ f(d2,lambda[d]))
  
  covariance(model2) <- a1 ~ f(a2,0.5)
  covariance(model2) <- d1 ~ f(d2,0.25)
  latent(model2) <- ~ a2+d2
  intercept(model2, ~ a2+d2) <- 0
  covariance(model2) <- c(a2,d2) ~ v(1)
  full <- list(model1,model2)
  ## #######
  isA <- length(grep("a",type))>0
  isC <- length(grep("c",type))>0
  isD <- length(grep("d",type))>0
  isE <- length(grep("e",type))>0
  isU <- length(grep("u",type))>0
  if (!isA) {
    kill(model1) <- ~ a1 + a2
    kill(model2) <- ~ a1 + a2
  }
  if (!isD) {
    kill(model1) <- ~ d1 + d2
    kill(model2) <- ~ d1 + d2
  }
  if (!isC) {
    kill(model1) <- ~ c1 + c2
    kill(model2) <- ~ c1 + c2
  }
  if (!isE) {
    kill(model1) <- ~ e1 + e2
    kill(model2) <- ~ e1 + e2
  }
  if (!isU) {
    kill(model1) <- ~ u##+sdu2
    kill(model2) <- ~ u##+sdu1
  }
  if (isU & isE) {
    regression(model1,outcomes[1],"e1") <- "lambda[e2]"
    regression(model1,outcomes[2],"e2") <- "lambda[e2]"
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
      regression(model1, to=outcomes[2], from=mykeep) <- regfix(model1)$label[trash,outcomes[2]]
      kill(model1) <- trash
    }

    dif <- wide2[,myvars[1]]-wide2[,myvars[2]]   
    mykeep <- myvars
    if (all(na.omit(dif)==00)) {
      mykeep <- mykeep[-2]
    }  
    trash <- setdiff(myvars,mykeep)
    if (length(mykeep)==1) {
      regression(model2, to=outcomes[2], from=mykeep) <- regfix(model2)$label[trash,outcomes[2]]
      kill(model2) <- trash
    }
  }

  if (!is.null(weight)) {
    weight <- paste(weight,1:2,sep=".")
    ## wide1[which(is.na(wide1[,weight[1]])),weight[1]] <- 0
    ## wide1[which(is.na(wide1[,weight[2]])),weight[2]] <- 0
    ## wide2[which(is.na(wide2[,weight[1]])),weight[1]] <- 0
    ## wide2[which(is.na(wide2[,weight[2]])),weight[2]] <- 0    
    ## wide1[which(is.na(wide1[,outcomes[1]])),c(outcomes[1],weight[1])] <- 0
    ## wide1[which(is.na(wide1[,outcomes[2]])),c(outcomes[2],weight[2])] <- 0
    ## wide2[which(is.na(wide2[,outcomes[1]])),c(outcomes[1],weight[1])] <- 0
    ## wide2[which(is.na(wide2[,outcomes[2]])),c(outcomes[2],weight[2])] <- 0
    estimator <- "weighted"
  }


  ## Estimate
  newkeep <- unlist(sapply(keep, function(x) paste(x,1:2,sep=".")))
  suppressWarnings(mg <- multigroup(list(model1,model2), list(wide1,wide2), missing=TRUE,fix=FALSE,keep=newkeep))
  if (is.null(estimator)) return(mg)

  e <- estimate(mg,weight=weight,estimator=estimator,fix=FALSE,...)
  res <- list(coefficients=e$opt$estimate, vcov=e$vcov, estimate=e, model=mg, full=full, call=cl, data=data, zyg=zyg, id=id, twinnum=twinnum, type=type, model.mz=model1, model.dz=model2, data.mz=wide1, data.dz=wide2)
  class(res) <- "twinlm"
  return(res)
}
