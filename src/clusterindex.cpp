#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

/* how many are there of the different clusters, similar to table(clusters) */ 
RcppExport SEXP nclust(SEXP iclusters) {
// {{{
BEGIN_RCPP
  ivec clusters = Rcpp::as<ivec>(iclusters);
  int  n = clusters.n_elem;
 
  int  uniqueclust=0; 
  int  maxclust=0;
  ivec nclust(n); nclust.fill(0);
 
  for (int  i=0;i<n;i++){
    if (nclust[clusters[i]]==0) uniqueclust+=1; 
    nclust[clusters[i]]+=1; 
    if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 
  
  return(List::create(Named("maxclust")=maxclust,
		       Named("nclust")=nclust,
		       Named("uniqueclust")=uniqueclust)); 
END_RCPP
 } // }}}

/* organize indeces to different clusters in matrix  of size nclust x maxclust */ 
/* If optionally matrix 'mat' is supplied, the rows of mat corresponding to clusters is returned */
RcppExport SEXP clusterindexM(SEXP iclusters, SEXP imednum, SEXP inum, SEXP x, SEXP all) {
// {{{
BEGIN_RCPP
  ivec clusters = Rcpp::as<ivec>(iclusters);
  int n = clusters.n_elem;
  bool All = Rcpp::as<bool>(all);
  int  uniqueclust=0; 
  int  maxclust=0;
  ivec nclust(n); nclust.fill(0);
  bool hasX = !((Rf_isNull)(x));

    
  for (int  i=0;i<n;i++) {
    if (nclust[clusters[i]]==0) uniqueclust+=1; 
    nclust[clusters[i]]+=1; 
    if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 
 
  ivec num = Rcpp::as<ivec>(inum); 
  int  mednum = Rcpp::as<int>(imednum);
  ivec unum;
  if (mednum==1) {
    // unum = unique(num);
    // maxclust = unum.n_elem;
    maxclust = max(num)+1;
  }     
  
  imat idclust = imat(uniqueclust,maxclust); idclust.fill(NA_INTEGER);
  ivec clustsize(uniqueclust); clustsize.fill(0);
  ivec firstclustid(uniqueclust); firstclustid.fill(0);

  if (mednum==0) {
    for (int  i=0;i<n;i++){
      idclust[(clustsize[clusters[i]])*(uniqueclust)+clusters[i]]=i; 
      clustsize[clusters[i]]+=1; 
      if (clustsize[clusters[i]]==1) firstclustid[clusters[i]]=i; 
    } 
  } else {
    for (int  i=0;i<n;i++){
      idclust[num[i]*(uniqueclust)+clusters[i]]=i; 
      clustsize[clusters[i]]+=1; 
      if (clustsize[clusters[i]]==1) firstclustid[clusters[i]]=i; 
    } 
  }
 
 if (hasX) {   
   mat X = Rcpp::as<mat>(x); 
   //if (X.n_rows!=n) 
   mat res(uniqueclust,X.n_cols); res.fill(0);
   for (unsigned i=0; i<idclust.n_rows; i++) {
     for (int k=0; k<clustsize[i]; k++) {
       res.row(i) += X.row(idclust(i,k));
     }
   }
   if (!All) { return(Rcpp::wrap(res)); }
   return(List::create(Named("clusters")=clusters,
		       Named("maxclust")=maxclust,
		       Named("idclustmat")=idclust,
		       Named("cluster.size")=clustsize,
		       Named("antclust")=clustsize,
		       Named("uniqueclust")=uniqueclust,
		       Named("firstclustid")=firstclustid,
		       Named("X")=res
		       ));    
 }


  return(List::create(Named("clusters")=clusters,
		      Named("maxclust")=maxclust,
		      Named("idclustmat")=idclust,
		      Named("cluster.size")=clustsize,
		      Named("antclust")=clustsize,
		      Named("uniqueclust")=uniqueclust,
		      Named("firstclustid")=firstclustid
		      )); 
END_RCPP
} // }}}

RcppExport SEXP familypairindex(SEXP iclustmat,SEXP iclustsize,SEXP inumfamindex) {
// {{{
  ivec clustsize = Rcpp::as<ivec>(iclustsize);
  imat clustmat;
  clustmat  = Rcpp::as<imat>(iclustmat);
  int uniqueclust=clustmat.n_rows;
  int  numfamindex= Rcpp::as<int>(inumfamindex);
 
  ivec famclustindex(numfamindex); famclustindex.fill(0);
  ivec subfamilyindex(numfamindex); subfamilyindex.fill(0);
  
  int i,j,v=0,h=0;

   for (i=0;i<uniqueclust;i++)
   {
	 if (clustsize(i)>=2)
         for (j=0;j<(clustsize(i)-1);j++)
         for (int k=j+1;k<clustsize(i);k++)
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
BEGIN_RCPP
  ivec clusters = Rcpp::as<ivec>(iclusters);
  int  n = clusters.n_elem;
 
  int  uniqueclust=0; 
  int  maxclust=0;
  ivec nclust(n); nclust.fill(0);
 
  for (int  i=0;i<n;i++){
    if (nclust[clusters[i]]==0) uniqueclust+=1; 
    nclust[clusters[i]]+=1; 
    if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 

  ivec num = Rcpp::as<ivec>(inum); 
  int  mednum = Rcpp::as<int>(imednum);
  
  imat idclust = imat(uniqueclust,maxclust); idclust.fill(NA_INTEGER);
  ivec clustsize(uniqueclust); clustsize.fill(0);

  mat data = Rcpp::as<mat>(idata);
  int  p= data.n_cols; 
  mat nydata(uniqueclust,maxclust*p); nydata.fill(NA_REAL);

  if (mednum==0) {
     for (int  i=0;i<n;i++){
//         idclust[(clustsize[clusters[i]])*(uniqueclust)+clusters[i]]=i; 
     for (int  j=0;j<p;j++) nydata[(clustsize[clusters[i]]*(p)+j)*(uniqueclust)+clusters[i]]=data[(n)*j+i]; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (int  i=0;i<n;i++){
//        idclust[num[i]*(uniqueclust)+clusters[i]]=i; 
        for (int  j=0;j<p;j++) nydata[(num[i]*(p)+j)*(uniqueclust)+clusters[i]]=data[(n)*j+i]; 
        clustsize[clusters[i]]+=1; 
     } 
  }

  return(List::create(
//	              Named("idclust")=idclust,
		      Named("maxclust")=maxclust,
   		      Named("clustsize")=clustsize,Named("iddata")=nydata)); 
END_RCPP
} // }}}

