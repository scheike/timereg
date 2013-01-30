#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

/* how many are there of the different clusters, similar to table(clusters) */ 
RcppExport SEXP nclust(SEXP iclusters) {
// {{{

  uvec clusters = Rcpp::as<uvec>(iclusters);
  unsigned n = clusters.n_elem;
 
  unsigned uniqueclust=0; 
  unsigned maxclust=0;
  uvec nclust(n); nclust.fill(0);
 
  for (unsigned i=0;i<n;i++){
    if (nclust[clusters[i]]==0) uniqueclust+=1; 
    nclust[clusters[i]]+=1; 
    if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 
  
  return(List::create(Named("maxclust")=maxclust,
		       Named("nclust")=nclust,
		       Named("uniqueclust")=uniqueclust)); 
 } // }}}



/* organize indeces to different clusters in matrix  of size nclust x maxclust */ 
RcppExport SEXP clusterindexM(SEXP iclusters, SEXP imednum, SEXP inum) {
// {{{

  uvec clusters = Rcpp::as<uvec>(iclusters);
  unsigned n = clusters.n_elem;
 
  unsigned uniqueclust=0; 
  unsigned maxclust=0;
  uvec nclust(n); nclust.fill(0);
 
  for (unsigned i=0;i<n;i++){
    if (nclust[clusters[i]]==0) uniqueclust+=1; 
    nclust[clusters[i]]+=1; 
    if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 

  uvec num = Rcpp::as<uvec>(inum); 
  int  mednum = Rcpp::as<int>(imednum);
  
  imat idclust = imat(uniqueclust,maxclust); idclust.fill(NA_INTEGER);
  uvec clustsize(uniqueclust); clustsize.fill(0);

  if (mednum==0) {
    for (unsigned i=0;i<n;i++){
      idclust[(clustsize[clusters[i]])*(uniqueclust)+clusters[i]]=i; 
      clustsize[clusters[i]]+=1; 
    } 
  } else {
    for (unsigned i=0;i<n;i++){
      idclust[num[i]*(uniqueclust)+clusters[i]]=i; 
      clustsize[clusters[i]]+=1; 
    } 
  }
  return(List::create(Named("clusters")=clusters,
		      Named("maxclust")=maxclust,
		      Named("idclustmat")=idclust,
		      Named("cluster.size")=clustsize,
		      Named("antclust")=clustsize
		      )); 
} // }}}




RcppExport SEXP clusterindexdata(SEXP iclusters,SEXP imaxclust,SEXP inclust,SEXP imednum,SEXP inum,
		SEXP idata) 
{ // {{{
  int i,j;
  uvec clusters = Rcpp::as<uvec>(iclusters); 
  int n = clusters.n_elem; 
  int maxclust = Rcpp::as<int>(imaxclust); 
  int nclust = Rcpp::as<int>(inclust); 
  uvec num = Rcpp::as<uvec>(inum); 
  int  mednum = Rcpp::as<int>(imednum);

  imat idclust = imat(nclust,maxclust); idclust.fill(NA_INTEGER);
  uvec clustsize(nclust); clustsize.fill(0);

  mat data = Rcpp::as<mat>(idata);
  int p= data.n_cols; 
  mat nydata(nclust,maxclust*p); nydata.fill(0);

  if (mednum==0) {
     for (i=0;i<n;i++){
         idclust[(clustsize[clusters[i]])*(nclust)+clusters[i]]=i; 
     for (j=0;j<p;j++) nydata[(clustsize[clusters[i]]*(p)+j)*(nclust)+clusters[i]]=data[(n)*j+i]; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (i=0;i<n;i++){
        idclust[num[i]*(nclust)+clusters[i]]=i; 
        for (j=0;j<p;j++) nydata[(num[i]*(p)+j)*(nclust)+clusters[i]]=data[(n)*j+i]; 
        clustsize[clusters[i]]+=1; 
     } 
  }

  return(List::create(Named("idclustmat")=idclust,Named("clustsize")=clustsize,Named("iddata")=nydata)); 
} // }}}

