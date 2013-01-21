faster.reshape <- function(data,clusters,index.type=FALSE,num=NULL,Rindex=1)
{ ## {{{
if (NCOL(data)==1) data <- cbind(data)
if (!is.matrix(data)) data <- as.matrix(data)

antpers <- length(clusters)
if (index.type==FALSE)  {
	max.clust <- length(unique(clusters))
	if (is.numeric(clusters)) clusters <-  sindex.prodlim(unique(clusters),clusters)-1 else 
	{
	max.clust <- length(unique(clusters))
	clusters <- as.integer(factor(clusters, labels = 1:max.clust))-1
	}
}

 nclust <- .C("nclusters",
	as.integer(antpers), as.integer(clusters), as.integer(rep(0,antpers)), 
	as.integer(0), as.integer(0), package="timereg")
  maxclust <- nclust[[5]]
  antclust <- nclust[[4]]
  cluster.size <- nclust[[3]][1:antclust]

if ((!is.null(num))) { ### different types in different columns
   mednum <- 1
   if (is.numeric(num)) num <-  sindex.prodlim(unique(num),num)-1
   else num <- as.integer(factor(num, labels = 1:maxclust)) -1
} else {num<-0; mednum<-0;}
p <- ncol(data); 
init <- -1*Rindex;

clustud <- .C("clusterindexdata",
	        as.integer(clusters), as.integer(antclust),as.integer(antpers),
                as.integer(rep(init,antclust*maxclust)),as.integer(rep(0,antclust)), as.integer(mednum), 
		as.integer(num), as.double(c(data)), 
		as.integer(p),   as.double(rep(init,antclust*maxclust*p)), package="timereg")
xny <- matrix(clustud[[10]],antclust,maxclust*p)
###if(Rindex==1) xny[idclust==-1] <- NA 
###if(Rindex==1) xny[idclust==-1] <- NA 
if (Rindex==1) idclust  <- matrix(clustud[[4]],antclust,maxclust)+1
else idclust <- matrix(clustud[[4]],antclust,maxclust)
if(Rindex==1) idclust[idclust==0] <- NA 

xnames <- colnames(data); 
missingname <- (colnames(data)=="")
xnames[missingname] <- paste(seq_len(maxclust))[missingname]
xny <- data.frame(xny)
mm <- as.vector(outer(xnames,seq_len(maxclust),function(...) paste(...,sep=".")))
names(xny) <- mm
out <- xny; 
} ## }}}

