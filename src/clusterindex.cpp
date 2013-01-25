#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

/* how many are there of the different clusters, similar to table(clusters) */ 
RcppExport SEXP nclust(SEXP in,SEXP iclusters) {

  uvec clusters = Rcpp::as<uvec>(iclusters); 
  int  n = Rcpp::as<int>(in);
  int uniqueclust=0; 

  uvec nclust = zeros(n); 

  int i,maxclust=0;
  for (i=0;i<n;i++){
      if (nclust[clusters[i]]==0) uniqueclust+=1; 
      nclust[clusters[i]]+=1; 
      if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 

return(Rcpp::List::create(Rccp::Named("maxclust")=mclust,Rccp::Named("nclust")=nclust,
			  Rccp::Named("uniqueclust")=uniqueclust)); 
}

/* organize indeces to different clusters in matrix  of size nclust x maxclust */ 
RcppExport SEXP clusterindex(SEXP in,SEXP iclusters,SEXP imacxlust,SEXP inclust,SEXP imednum,SEXP inum) 
{

  int i;
  int  n = Rcpp::as<int>(in);
  uvec clusters = Rcpp::as<uvec>(iclusters); 
  int maxclust = Rcpp::as<int>(imaxclust); 
  int nclust = Rcpp::as<int>(inclust); 
  uvec num = Rcpp::as<uvec>(inum); 
  int  mednum = Rcpp::as<int>(imednum);

  umat idclust = zeros(nclust,maxclust); 
  uvec clustsize = zeros(nclust); 

  if (mednum==0) {
     for (i=0;i<n;i++){
         idclust[(clustsize[clusters[i]])*(nclust)+clusters[i]]=i; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (i=0;i<n;i++){
        idclust[num[i]*(nclust)+clusters[i]]=i; 
        clustsize[clusters[i]]+=1; 
     } 
  }
return(Rcpp::List::create(Rccp::Named("idclustmat")=idclust,Rccp::Named("clustsize")=clustsize)); 
}

//RcppExport SEXP clusterindexdata(SEXP npers,SEXP clusters,SEXP nclust,SEXP mednum,SEXP num,SEXP data) 
//{ // {{{
//
//  int i,j;
//  if (*mednum==0) {
//     for (i=0;i<*npers;i++){
//         idclust[(clustsize[clusters[i]])*(*nclust)+clusters[i]]=i; 
//     for (j=0;j<*p;j++) nydata[(clustsize[clusters[i]]*(*p)+j)*(*nclust)+clusters[i]]=data[(*npers)*j+i]; 
//         clustsize[clusters[i]]+=1; 
//      } 
//  } else {
//    for (i=0;i<*npers;i++){
//        idclust[num[i]*(*nclust)+clusters[i]]=i; 
//        for (j=0;j<*p;j++) nydata[(num[i]*(*p)+j)*(*nclust)+clusters[i]]=data[(*npers)*j+i]; 
//        clustsize[clusters[i]]+=1; 
//     } 
//  }
//
//return(Rcpp::List::create( idclust, clustsize,nydata)); 
//} // }}}
//


