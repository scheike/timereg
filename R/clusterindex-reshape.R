##' Finds subjects related to same cluster
##' 
##' @aliases cluster.index 
##' @references
##' Cluster indeces 
##' @examples
##' i<-c(1,1,2,2,1,3)
##' d<- cluster.index(i)
##' print(d)
##' @keywords cluster indeces 
##' @author Klaus Holst, Thomas Scheike
##' @export
cluster.index <- function(clusters,index.type=FALSE,num=NULL,Rindex=0,mat=NULL,return.all=FALSE,code.na=NA)
{ ## {{{
  n <- length(clusters)

  if (index.type==FALSE)  {
    if (is.numeric(clusters)) clusters <-  fast.approx(unique(clusters),clusters)-1 else  {
      max.clust <- length(unique(clusters))
      clusters <- as.integer(factor(clusters, labels = seq(max.clust)))-1
    }
  }
  
  if ((!is.null(num))) { ### different types in different columns
    mednum <- 1
    if (is.numeric(num)) numnum <-  fast.approx(unique(num),num)-1
    else {
      numnum <- as.integer(factor(num, labels = seq(length(unique(clusters))))) -1
    }
  } else { numnum <- 0; mednum <- 0; }

  clustud <- .Call("clusterindexM",as.integer(clusters),as.integer(mednum), as.integer(numnum),mat,return.all)
  if (!is.null(mat) && !return.all) return(clustud)
  
  if (Rindex==1) clustud$idclust <- clustud$idclustmat+1
  if (Rindex==1) clustud$firstclustid <- clustud$firstclustid +1 
  ### avoid NA's for C call
  if (Rindex==0 & !is.na(code.na)) clustud$idclust[is.na(clustud$idclust)] <- code.na
  
  clustud
} ## }}}

##' Finds all pairs within a cluster (family)
##' 
##' @aliases familycluster.index 
##' @references
##' Cluster indeces 
##' @examples
##' i<-c(1,1,2,2,1,3)
##' d<- familycluster.index(i)
##' print(d)
##' @keywords cluster indeces 
##' @author Klaus Holst, Thomas Scheike
##' @export
familycluster.index <- function(clusters,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
  clusters <- cluster.index(clusters,Rindex=Rindex)
  totpairs <- sum(clusters$cluster.size*(clusters$cluster.size-1)/2)
  clustud <- .Call("familypairindex",clusters$idclust,clusters$cluster.size,as.integer(2*totpairs))
  clustud$pairs <- matrix(clustud$familypairindex,ncol=2,byrow=TRUE)

  invisible(clustud)
} ## }}}

##' Finds all pairs within a cluster (famly)  with the proband (case/control) 
##' 
##' second column of pairs are the probands and the first column the related subjects
##' 
##' 
##' @aliases familycluster.index 
##' @references
##' Cluster indeces 
##' @examples
##' i<-c(1,1,2,2,1,3)
##' p<-c(1,0,0,1,0,1)
##' d<- familyclusterWithProbands.index(i,p)
##' print(d)
##' @keywords cluster indeces 
##' @author Klaus Holst, Thomas Scheike
##' @export
familyclusterWithProbands.index <- function(clusters,probands,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
    famc <-familycluster.index(clusters,index.type=index.type,num=num,Rindex=Rindex)
    if (length(probands)!=length(clusters)) stop("clusters and probands not same length\n"); 
    index.probs <- (1:length(clusters))[probands==1]
    subfamsWprobands <-famc$subfamilyindex[ famc$familypairindex %in% index.probs ]
    indexWproband <- famc$subfamilyindex %in% subfamsWprobands 
    famc$subfamilyindex <- famc$subfamilyindex[indexWproband]
    famc$familypairindex <- famc$familypairindex[indexWproband]
    pairs <- matrix(famc$familypairindex,ncol=2,byrow=TRUE)
    ipi1 <- pairs[,1] %in% index.probs
    gem2 <- pairs[,2]
    pairs[ipi1,2] <- pairs[ipi1,1]
    pairs[ipi1,1] <- gem2[ipi1]
    famc$pairs <- pairs

    famc$familypairindex <- c(t(pairs))
    invisible(famc)
} ## }}}


###library(mets)
###clusters <-   c(1,1,2,2,1,3,3,3,4,4)
###probands <-   c(0,1,0,1,0,1,0,0,0,0)
###index.type=FALSE;num=NULL;Rindex=1
###ilusters <- cluster.index(clusters,Rindex=1)
###ud <- familycluster.index(clusters)
###ud1 <- familyclusterWithProbands.index(clusters,probands)

##' @export
coarse.clust <- function(clusters,max.clust=100)
{ ## {{{ 

if (is.numeric(clusters)) 
   clusters <-  sindex.prodlim(unique(clusters),clusters)
cluster.size <- length(unique(clusters))

qq <- unique(quantile(clusters, probs = seq(0, 1, by = 1/max.clust)))
qqc <- cut(clusters, breaks = qq, include.lowest = TRUE)    
cclusters <-  as.integer(qqc)-1

return(cclusters)
} ## }}} 

##' @export
faster.reshape <- function(data,clusters,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
  if (NCOL(data)==1) data <- cbind(data)
 ### uses data.matrix 
  if (!is.matrix(data)) data <- data.matrix(data)
  if (is.character(clusters)) clusters <- data[,clusters]
  n <- length(clusters)

  if (nrow(data)!=n)  stop("nrow(data) and clusters of different lengths\n"); 

  if (index.type==FALSE)  {
    max.clust <- length(unique(clusters))
    if (is.numeric(clusters)) clusters <-  fast.approx(unique(clusters),clusters)-1 else 
    {
      max.clust <- length(unique(clusters))
      clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1
    }
  }

  if ((!is.null(num))) { ### different types in different columns
    if (length(num)!=n)  stop("clusters and num of different lengths\n"); 
    mednum <- 1
    if (is.numeric(num)) num <-  fast.approx(unique(num),num)-1
    else num <- as.integer(factor(num, labels = seq(length(unique(clusters))))) -1
  } else { num <- 0; mednum <- 0; }

  clustud <- .Call("clusterindexdata",as.integer(clusters),as.integer(mednum), as.integer(num),iddata=data)

  if (Rindex==1) clustud$idclust  <- clustud$idclust+1
###  if(Rindex==1) idclust[idclust==0] <- NA 
  maxclust <- clustud$maxclust

  xny <- clustud$iddata
  xnames <- colnames(data); 
  missingname <- (colnames(data)=="")
  xnames[missingname] <- paste(seq_len(maxclust))[missingname]
  xny <- data.frame(xny)
  mm <- as.vector(outer(xnames,seq_len(maxclust),function(...) paste(...,sep=".")))
  names(xny) <- mm

  return(xny); 
} ## }}}
