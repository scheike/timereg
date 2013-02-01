##' @export
cluster.index <- function(clusters,index.type=FALSE,num=NULL,Rindex=0)
{ ## {{{
  n <- length(clusters)

  if (index.type==FALSE)  {
    if (is.numeric(clusters)) clusters <-  timereg:::sindex.prodlim(unique(clusters),clusters)-1 else  {
      max.clust <- length(unique(clusters))
      clusters <- as.integer(factor(clusters, labels = seq(max.clust)))-1
    }
  }
  
  if ((!is.null(num))) { ### different types in different columns
    mednum <- 1
    if (is.numeric(num)) numnum <-  timereg:::sindex.prodlim(unique(num),num)-1
    else {
      numnum <- as.integer(factor(num, labels = seq(length(unique(clusters))))) -1
    }
  } else { numnum <- 0; mednum <- 0; }

  clustud <- .Call("clusterindexM",as.integer(clusters),as.integer(mednum), as.integer(numnum))
  
  if (Rindex==1) clustud$idclust <- clustud$idclustmat+1
  
  invisible(clustud)
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
    if (is.numeric(clusters)) clusters <-  timereg:::sindex.prodlim(unique(clusters),clusters)-1 else 
    {
      max.clust <- length(unique(clusters))
      clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1
    }
  }

  if ((!is.null(num))) { ### different types in different columns
    if (length(num)!=n)  stop("clusters and num of different lengths\n"); 
    mednum <- 1
    if (is.numeric(num)) num <-  timereg:::sindex.prodlim(unique(num),num)-1
    else num <- as.integer(factor(num, labels = seq(length(unique(clusters))))) -1
  } else { num <- 0; mednum <- 0; }

  clustud <- .Call("clusterindexdata",as.integer(clusters),as.integer(mednum), as.integer(num),iddata=data,DUP=FALSE)

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

