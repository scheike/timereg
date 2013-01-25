#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

/* how many are there of the different clusters, similar to table(clusters) */ 
RcppExport SEXP nclust(SEXP npers,SEXP clusters) {
  try {

  int i,maxclust=0;
  for (i=0;i<*npers;i++){
      if (nclust[clusters[i]]==0) uniqueclust[0]+=1; 
      nclust[clusters[i]]+=1; 
      if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 

  mclust[0]=maxclust; 

  }
return(Rcpp::List::create( mclust, nclust)); 
}

/* organize indeces to different clusters in matrix  of size nclust x maxclust */ 
RcppExport SEXP clusterindex(SEXP npers,SEXP clusters,SEXP nclust,SEXP mednum,SEXP num) 
{
  int i;

  if (*mednum==0) {
     for (i=0;i<*npers;i++){
         idclust[(clustsize[clusters[i]])*(*nclust)+clusters[i]]=i; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (i=0;i<*npers;i++){
        idclust[num[i]*(*nclust)+clusters[i]]=i; 
        clustsize[clusters[i]]+=1; 
     } 
  }
return(Rcpp::List::create( idclust, clustsize)); 
}

RcppExport SEXP clusterindexdata(SEXP npers,SEXP clusters,SEXP nclust,SEXP mednum,SEXP num,SEXP data) 
{

  int i,j;
  if (*mednum==0) {
     for (i=0;i<*npers;i++){
         idclust[(clustsize[clusters[i]])*(*nclust)+clusters[i]]=i; 
     for (j=0;j<*p;j++) nydata[(clustsize[clusters[i]]*(*p)+j)*(*nclust)+clusters[i]]=data[(*npers)*j+i]; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (i=0;i<*npers;i++){
        idclust[num[i]*(*nclust)+clusters[i]]=i; 
        for (j=0;j<*p;j++) nydata[(num[i]*(*p)+j)*(*nclust)+clusters[i]]=data[(*npers)*j+i]; 
        clustsize[clusters[i]]+=1; 
     } 
  }

return(Rcpp::List::create( idclust, clustsize,nydata)); 
}


