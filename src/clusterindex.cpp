#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

/* how many are there of the different clusters, similar to table(clusters) */ 
RcppExport SEXP nclust(SEXP in, SEXP iclusters) {
// {{{
  uvec clusters = Rcpp::as<uvec>(iclusters); 
  int  n = Rcpp::as<int>(in);

  int uniqueclust=0; 
  uvec nclust(n); nclust.fill(0);
  int i,maxclust=0;

  for (i=0;i<n;i++){
    if (nclust[clusters[i]]==0) uniqueclust+=1; 
    nclust[clusters[i]]+=1; 
    if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 

  return(List::create(Named("maxclust")=maxclust,
		      Named("nclust")=nclust,
		      Named("uniqueclust")=uniqueclust)); 
} // }}}

/* organize indeces to different clusters in matrix  of size nclust x maxclust */ 
RcppExport SEXP clusterindexM(SEXP in, SEXP iclusters, SEXP imaxclust,
			     SEXP inclust, SEXP imednum, SEXP inum) {
// {{{
  int i;
  int  n = Rcpp::as<int>(in);
  uvec clusters = Rcpp::as<uvec>(iclusters); 
  int maxclust = Rcpp::as<int>(imaxclust); 
  int nclust = Rcpp::as<int>(inclust); 
  uvec num = Rcpp::as<uvec>(inum); 
  int  mednum = Rcpp::as<int>(imednum);

  umat idclust = umat(nclust,maxclust); idclust.fill(0);
  uvec clustsize(nclust); clustsize.fill(0);

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
  return(List::create(Named("idclustmat")=idclust,
		      Named("clustsize")=clustsize)); 
} // }}}


RcppExport SEXP clusterindexdata(SEXP in,SEXP iclusters,SEXP imaxclust,SEXP inclust,SEXP imednum,SEXP inum,
		SEXP idata) 
{ // {{{
  int i,j;
  int  n = Rcpp::as<int>(in);
  uvec clusters = Rcpp::as<uvec>(iclusters); 
  int maxclust = Rcpp::as<int>(imaxclust); 
  int nclust = Rcpp::as<int>(inclust); 
  uvec num = Rcpp::as<uvec>(inum); 
  int  mednum = Rcpp::as<int>(imednum);

  umat idclust = umat(nclust,maxclust); idclust.fill(0);
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

