// *[[Rcpp::depends(RcppArmadillo)]]
// *[[Rcpp::interfaces(r, cpp)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
#include <vector>

// const double epsilon=1.0e-16;

using namespace std;
using namespace Rcpp;

// [[Rcpp::export(name = ".ApplyBy2")]]
NumericMatrix ApplyBy2(NumericMatrix idata,
		       NumericVector icluster,
		       SEXP F,
		       Environment Env,
		       std::string Argument="x",
		       int Columnwise=0,
		       int Reduce=0,
		       double epsilon=1.0e-16
		       ) {
BEGIN_RCPP
   unsigned n = idata.nrow();
   unsigned posstart=0;
   unsigned p = idata.ncol();
   if (Columnwise==0) {
     Env[Argument] = idata( Range(0,0), Range(0, idata.ncol()-1) );
   } else {
     Env[Argument] = 1;
   }
   NumericVector res1(Rf_eval(F,Env));
   unsigned nf = res1.size();
   unsigned P=nf;
   if (Columnwise!=0) {
     P *= idata.ncol();
   }   

   // First we get the number of clusters
   double curcluster=icluster[0];
   double prevcluster=icluster[0];
   unsigned nclusters=1;
   for (unsigned i=0; i<n; i++) {
     curcluster=icluster[i];
     if (fabs(curcluster-prevcluster)>epsilon) {
       nclusters++;
     }
     prevcluster=curcluster;
   }

   unsigned clpos=0;
   NumericVector clustersize(nclusters);
   curcluster=icluster[0];
   prevcluster=icluster[0];
   NumericMatrix res(n,P);
   for (unsigned i=0; i<=n; i++) {
     if (i<n) curcluster=icluster[i];    
     if (fabs(curcluster-prevcluster)>epsilon || i==n) {
       if (i<n) clpos++;
       unsigned pos=i-1;
       unsigned nr=pos-posstart+1;
       NumericMatrix tmp = idata(Range(posstart,pos),_);
       NumericVector val;
       unsigned valsize=1;
       bool valmany=1;
       if (Columnwise==0 || p==1) {
	 // Apply function on matrix
	 Env[Argument] = tmp;
	 SEXP fres = Rf_eval(F,Env);
	 val = NumericVector(fres);
	 valsize = val.size();
	 valmany = (valsize>=(nr*nf));
       } else {
	 // Apply function on each column 
	 val = NumericVector(nr*P);
	 valsize = 1;
	 for (unsigned k=0; k<p; k++) { // column-wise
	   Env[Argument] = tmp(_,k);
	   SEXP fres = Rf_eval(F,Env);
	   NumericVector val0(fres);
	   valsize = val0.size();
	   if (valsize<nr) { // Results is just one row
	     valmany = 0;
	     for (unsigned j=0; j<valsize; j++) {
	       val[(j+k*nf)*nr] = val0[j];
	     }  
	   } else for (unsigned j=0; j<valsize; j++) {
	     val[j + k*(nf*nr)] = val0[j];
	   }
	 }
	 valsize = valsize*nf;
       }
       double myval=0;
       for (unsigned j=0; j<nr; j++) {
	 for (unsigned c=0; c<P; c++) {	   
	   if (valmany) {
	     myval = val[j + nr*c];
	   } else {
	     unsigned pos=c;
	     if (Columnwise!=0) pos *=nr;
	     myval = val[pos];
	   }
	   res(j+posstart,c) = myval;	   
	 }
       }       
       posstart=i;
     }
     if (i<n) clustersize[clpos]++;
     prevcluster=curcluster;
   }

   res.attr("clustersize") = clustersize;
   return(res);
END_RCPP
}  



// [[Rcpp::export(name = ".ApplyBy")]]
NumericMatrix ApplyBy(NumericMatrix idata,
		      IntegerVector icluster,
		      Function f) {
  /*
    This version uses Rcpp::Function but is much slow due to repeatedly calling tryCatch
  */
BEGIN_RCPP
   unsigned n = idata.nrow();
   NumericMatrix row1 = idata( Range(0,0), Range(0, idata.ncol()-1) );
   Function asMatrix("as.matrix");
   NumericMatrix res1 = asMatrix(f(row1));
   unsigned p = res1.ncol()*res1.nrow();
   NumericMatrix res(n,p);
   unsigned posstart=0;
   unsigned curcluster=icluster[0];
   unsigned prevcluster=icluster[0];
   for (unsigned i=0; i<=n; i++) {
     if (i<n) curcluster=icluster[i];
     if (curcluster!=prevcluster || i==n) {
       unsigned pos=i-1;
       unsigned nr=pos-posstart+1;
       NumericMatrix tmp = idata(Range(posstart,pos),_);
       // Rcout << "tmp=" << tmp << std::endl << std::endl;       
       // NumericMatrix val = asmatrix(f(tmp));
       // unsigned valsize = val.nrow()*val.ncol();       
       NumericVector val = Rcpp::as<NumericVector>(f(tmp));
       unsigned valsize = val.size();
       // Rcout << "val=" << val << std::endl << std::endl;	      
       bool valmany = (valsize==(nr*p));
       double myval=0;
       for (unsigned j=0; j<nr; j++) {
	 for (unsigned c=0; c<p; c++) {	   
	   if (valmany) {
	     myval = val[j + nr*c];
	   } else {
	     myval = val[c];
	   }
	   res(j+posstart,c) = myval;	   
	 }
       }       
       posstart=i;
     }
     prevcluster=curcluster;
   }
      
   return(res);
END_RCPP
}  
