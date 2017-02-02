// *[[Rcpp::depends(RcppArmadillo)]]
// *[[Rcpp::interfaces(r, cpp)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export(name = ".ApplyBy2")]]
NumericMatrix ApplyBy2(NumericMatrix idata,
		       IntegerVector icluster,
		       SEXP F,
		       Environment Env
		       ) {
BEGIN_RCPP
   unsigned n = idata.nrow();
   NumericMatrix row1 = idata( Range(0,0), Range(0, idata.ncol()-1) );
   unsigned posstart=0;
   unsigned curcluster=icluster[0];
   unsigned prevcluster=icluster[0];

   Env["x_"] = row1;
   NumericVector res1(Rf_eval(F,Env));
   unsigned p = res1.size();;
   NumericMatrix res(n,p);
   for (unsigned i=0; i<=n; i++) {
     if (i<n) curcluster=icluster[i];
     if (curcluster!=prevcluster || i==n) {       
       unsigned pos=i-1;
       unsigned nr=pos-posstart+1;
       NumericMatrix tmp = idata(Range(posstart,pos),_);
       Env["x_"] = tmp;
       SEXP fres = Rf_eval(F,Env);
       NumericVector val(fres);
       unsigned valsize = val.size();
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

//* [[Rcpp::export(name = ".ApplyBy")]]
// [[Rcpp::export]]
NumericMatrix ApplyBy(NumericMatrix idata,
		      IntegerVector icluster,
		      Function f) {
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
