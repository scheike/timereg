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
		      Named("antclust")=clustsize,
		      Named("uniqueclust")=uniqueclust
		      )); 
} // }}}

RcppExport SEXP familypairindex(SEXP iclustmat,SEXP iclustsize,SEXP inumfamindex) {
// {{{
  uvec clustsize = Rcpp::as<uvec>(iclustsize);
  umat clustmat;
  clustmat  = Rcpp::as<umat>(iclustmat);
  int uniqueclust=clustmat.n_rows;
  int  numfamindex= Rcpp::as<int>(inumfamindex);
 
  uvec famclustindex(numfamindex); famclustindex.fill(0);
  uvec subfamilyindex(numfamindex); subfamilyindex.fill(0);
  
  int i,j,k,v=0,h=0;

   for (i=0;i<uniqueclust;i++)
   {
	 if (clustsize(i)>=2)
         for (j=0;j<(clustsize(i)-1);j++)
         for (k=j+1;k<clustsize(i);k++)
         {
//	  Rprintf(" %d %d %d %d %d n=%d c=%d \n",i,j,k,h,v,numfamindex,clustsize(i)); 
	    famclustindex(v)=clustmat(i,j);
	    subfamilyindex(v)=h;
	    v+=1;
	    famclustindex(v)=clustmat(i,k);
	    subfamilyindex(v)=h;
	    v+=1;
	    h+=1;
	 }
  }

  return(List::create(Named("familypairindex")=famclustindex,
		       Named("subfamilyindex")=subfamilyindex
		       )); 
} // }}}

RcppExport SEXP clusterindexdata(SEXP iclusters, SEXP imednum,SEXP inum, SEXP idata) 
{ // {{{
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

  mat data = Rcpp::as<mat>(idata);
  unsigned p= data.n_cols; 
  mat nydata(uniqueclust,maxclust*p); nydata.fill(NA_REAL);

  if (mednum==0) {
     for (unsigned i=0;i<n;i++){
//         idclust[(clustsize[clusters[i]])*(uniqueclust)+clusters[i]]=i; 
     for (unsigned j=0;j<p;j++) nydata[(clustsize[clusters[i]]*(p)+j)*(uniqueclust)+clusters[i]]=data[(n)*j+i]; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (unsigned i=0;i<n;i++){
//        idclust[num[i]*(uniqueclust)+clusters[i]]=i; 
        for (unsigned j=0;j<p;j++) nydata[(num[i]*(p)+j)*(uniqueclust)+clusters[i]]=data[(n)*j+i]; 
        clustsize[clusters[i]]+=1; 
     } 
  }

  return(List::create(
//	              Named("idclust")=idclust,
		      Named("maxclust")=maxclust,
   		      Named("clustsize")=clustsize,Named("iddata")=nydata)); 
} // }}}

