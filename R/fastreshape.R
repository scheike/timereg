##' Simple reshape/tranpose of data
##'
##' @title Fast reshape
##' @param data data.frame or matrix
##' @param id id-variable. If omitted then reshape from wide to long. 
##' @param varying Vector of prefix-names of the time varying
##' variables. Optional for long->wide reshaping.
##' @param num Optional number/time variable
##' @param sep String seperating prefix-name with number/time
##' @param ... Optional additional arguments to the \code{reshape} function used in the wide->long reshape.
##' @author Thomas Scheike, Klaus K. Holst
##' @export
##' @examples
##' m <- lvm(c(y1,y2,y3,y4)~x)
##' d <- sim(m,1e1)
##' 
##' dd <- fast.reshape(d,var="y")
##' d1 <- fast.reshape(dd,"id")
##'
##' ## From wide-format
##' d1 <- fast.reshape(dd,"id")
##' d2 <- fast.reshape(dd,"id",var="y")
##' d3 <- fast.reshape(dd,"id",var="y",num="time")
##' 
##' d4 <- fast.reshape(data.matrix(dd),"id",var="y")
##' 
##' ## From long-format
##' fast.reshape(d,var="y",idvar="a",timevar="b")
##' fast.reshape(d,var=list(c("y1","y2","y3","y4")),idvar="a",timevar="b")
##'
##' data(prt)
##' head(fast.reshape(prt,"id",var="cancer"))
fast.reshape <- function(data,id,varying,num,sep="",drop=NULL,
                         idname="id",timename="num",...) {
  if (NCOL(data)==1) data <- cbind(data)
'.onAttach' <- function(lib, pkg="mets")
  {    
    desc <- packageDescription(pkg)
    packageStartupMessage("\nLoading '", desc$Package, "' package...\n",
                          "Version    : ", desc$Version, "\n",
                          "Overview: help(package=", desc$Package, ")\n");
  }
  
  if (missing(id)) {
    ## reshape from wide to long format. Fall-back to stats::reshape
    nn <- colnames(data)
    nsep <- nchar(sep)
    vnames <- NULL
    if (missing(varying)) stop("Prefix of time-varying variables needed")    
    ncvar <- sapply(varying,nchar)
    newlist <- c()
    if (!is.list(varying)) {
      for (i in seq_len(length(varying))) {
        ii <- which(varying[i]==substr(nn,1,ncvar[i]))
        tt <- as.numeric(substring(nn[ii],ncvar[i]+1+nsep))      
        newlist <- c(newlist,list(nn[ii[order(tt)]]))
      }
      vnames <- varying      
      varying <- newlist
    }
    is_df <- is.data.frame(data)
    if (is_df) data <- data.matrix(data)
    
    fixed <- setdiff(nn,unlist(c(varying,drop,idname,timename)))
    nfixed <- length(fixed)
    nvarying <- length(varying)
    nclusts <- unlist(lapply(varying,length))
    nclust <- length(varying[[1]])
    if (any(nclusts!=nclust)) stop("Different length of varying vectors!")
    data <- data[,c(fixed,unlist(varying))]
    
    long <- as.data.frame(.Call("FastLong",idata=data,inclust=as.integer(nclust),
                                as.integer(nfixed),as.integer(nvarying)));
    colnames(long) <- c(fixed,vnames,idname,timename)

    if (is_df) { ## Recreate classes 
      vars.orig <- c(fixed,unlist(lapply(varying,function(x) x[1])))
    }   
    
    ##    if (test) return(reshape(as.data.frame(data),varying=varying,direction="long",v.names=vnames,...))
    return(long)
  }

  numvar <- idvar <- NULL 
  if (is.character(id) || is.factor(id)) {
    if (length(id)>1) stop("Expecting column name or vector of id's")
    idvar <- id
    id <- as.integer(data[,id,drop=TRUE])
  } else {
    if (length(id)!=nrow(data)) stop("Length of ids and data-set does not agree")
  }    
  if (!missing(num)) {
    if (is.character(num) || is.factor(num)) {
      numvar <- num
      num <- as.integer(data[,num,drop=TRUE])
    } else {
      if (length(num)!=nrow(data)) stop("Length of time and data-set does not agree")
    }
  } else {
    num <- NULL
  }  

  nn <- colnames(data)
  if (any(nn=="")) data <- data.frame(data)
 
  clustud <- cluster.index(id)
  maxclust <- clustud$maxclust
  idclust <- clustud$idclust  
  
  if (!is.null(numvar)) {
    ii <- which(colnames(data)==numvar)
    data <- data[,-ii,drop=FALSE]
  }
  if (missing(varying)) varying <- setdiff(colnames(data),c(idvar))
  vidx <- match(varying,colnames(data))
  N <- nrow(idclust)
  p <- length(varying)

  if (is.matrix(data) || all(apply(data[1,],2,is.numeric))) {
    ## Everything numeric - we can work with matrices
    dataw <- matrix(NA, nrow = N, ncol = p * (maxclust-1) + ncol(data))
    for (i in seq_len(maxclust)) {
      if (i==1) {
        dataw[, seq(ncol(data))] <- as.matrix(data[idclust[, i] + 1,])
        mnames <- colnames(data);
        mnames[vidx] <- paste(mnames[vidx],i,sep=sep)
      } else {
        dataw[, seq(p) + (ncol(data)-p) + (i - 1) * p] <- as.matrix(data[idclust[, i] + 1,varying])
        ##        mnames <- c(mnames,paste(varying,i,sep=sep))
      }
    }
    mnames <- c(mnames,as.vector(t(outer(varying,seq_len(maxclust-1)+1,function(...) paste(...,sep=sep)))))
    colnames(dataw) <- mnames
    return(dataw)
  } ## Potentially slower with data.frame where we use cbind
  
  for (i in seq_len(maxclust)) {
    if (i==1) {
      dataw <- data[idclust[,i]+1,,drop=FALSE]
      mnames <- names(data);
      mnames[vidx] <- paste(mnames[vidx],sep,i,sep="")
    } else {
      dataw <- cbind(dataw,data[idclust[,i]+1,varying,drop=FALSE])
      mnames <- c(mnames,paste(varying,sep,i,sep=""))
    }
  }
  names(dataw) <- mnames
  
  return(dataw)
 } 


simple.reshape <- function (data, id = "id", num = NULL) {
    cud <- cluster.index(data[, c(id)], num = num, Rindex = 1)
    N <- nrow(cud$idclust)
    p <- ncol(data)
    dataw <- matrix(NA, nrow = N, ncol = p * cud$maxclust)
    for (i in seq_len(cud$maxclust)) {
           dataw[, seq(p) + (i - 1) * p] <- as.matrix(data[cud$idclust[, i] + 1, ])
    }
   colnames(dataw) <- paste(names(data), rep(seq_len(cud$maxclust), each = p), sep = ".")
   return(dataw)
}



