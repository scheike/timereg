
##' @export
cluster.index <- function(clusters,index.type=FALSE,num=NULL,Rindex=0)
{ ## {{{
n <- length(clusters)

if (index.type==FALSE)  {
	if (is.numeric(clusters)) clusters <-  sindex.prodlim(unique(clusters),clusters)-1 else  {
	   max.clust <- length(unique(clusters))
	   clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1
	}
}

 nclust <- .Call("nclust", as.integer(n), as.integer(clusters))
 maxclust <- nclust$maxclust
 antclust <- nclust$uniqueclust
 nclust <-   nclust$nclust[1:nclust$uniqueclust]

if ((!is.null(num))) { ### different types in different columns
   mednum <- 1
if (is.numeric(num)) numnum <-  sindex.prodlim(unique(num),num)-1
else
   numnum <- as.integer(factor(num, labels = 1:maxclust)) -1
} else { numnum <- 0; mednum <- 0; }

clustud <- .Call("clusterindexM",as.integer(n),as.integer(clusters), 
                as.integer(maxclust), as.integer(antclust),
		as.integer(mednum), as.integer(numnum))

if (Rindex==1) idclust  <- clustud$idclustmat+1 else idclust <- clustud$idclustmat+1
if(Rindex==1) idclust[idclust==0] <- NA 
out <- list(clusters=clusters,maxclust=maxclust,antclust=antclust,idclust=idclust,cluster.size=clustud$clustsize)
} ## }}}

komud<-function(){  ## {{{ 
library(mets)
clusters <- c(1,1,2,2,1,3)
if (is.numeric(clusters)) clusters <-  sindex.prodlim(unique(clusters),clusters)-1 
clusters
n <- length(clusters)
nclust <- .Call("nclust", as.integer(n), as.integer(clusters))
num <- numnum <- 0; mednum <- 0;
  maxclust <- nclust$maxclust
  antclust <- nclust$uniqueclust
  nclust <-   nclust$nclust[1:nclust$uniqueclust]

###
clustud <- .Call("clusterindexM",as.integer(n),as.integer(clusters), 
                as.integer(maxclust), as.integer(antclust),
		as.integer(mednum), as.integer(numnum))

out=cluster.index(clusters,Rindex=1)

data <- x <- matrix(1:12,6,2)
clusters <- c(1,1,2,2,1,3)

out=faster.reshape(x,clusters)

out=faster.reshapeM(x,clusters)

clustud <- .Call("clusterindexdata",as.integer(n),as.integer(clusters), 
                as.integer(maxclust), as.integer(antclust),
		as.integer(mednum), as.integer(num),idata=data,DUP=FALSE)
} ## }}}

##' @export
faster.reshape <- function(data,clusters,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
if (NCOL(data)==1) data <- cbind(data)
if (!is.matrix(data)) data <- as.matrix(data)

n <- length(clusters)
if (index.type==FALSE)  {
	max.clust <- length(unique(clusters))
	if (is.numeric(clusters)) clusters <-  sindex.prodlim(unique(clusters),clusters)-1 else 
	{
	max.clust <- length(unique(clusters))
	clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1
	}
}

 nclust <- .Call("nclust", as.integer(n), as.integer(clusters))
 maxclust <- nclust$maxclust
 antclust <- nclust$uniqueclust
 nclust <-   nclust$nclust[1:nclust$uniqueclust]

if ((!is.null(num))) { ### different types in different columns
   mednum <- 1
   if (is.numeric(num)) num <-  sindex.prodlim(unique(num),num)-1
   else num <- as.integer(factor(num, labels = 1:maxclust)) -1
} else { num <- 0; mednum <- 0; }

p <- ncol(data); 
init <- -1*Rindex;


###RcppExport SEXP clusterindexdata(SEXP in,SEXP iclusters,SEXP imaxclust,SEXP inclust,SEXP imednum,SEXP inum,
###		SEXP idata) 

clustud <- .Call("clusterindexdata",as.integer(n),as.integer(clusters), 
                as.integer(maxclust), as.integer(antclust),
		as.integer(mednum), as.integer(num),iddata=data,DUP=FALSE)

###print(clustud)
###clustud <- .C("clusterindexdata",
###	        as.integer(clusters), as.integer(antclust),as.integer(antpers),
###                as.integer(rep(init,antclust*maxclust)),as.integer(rep(0,antclust)), as.integer(mednum), 
###		as.integer(num), as.double(data) )

if (Rindex==1) idclust  <- clustud$idclustmat+1 else idclust <- clustud$idclustmat+1
if(Rindex==1) idclust[idclust==0] <- NA 
xny <- clustud$iddata

###mnames <- c()
###names(xny) <- mnames
return(xny); 
} ## }}}
