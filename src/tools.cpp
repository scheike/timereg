#include "tools.h"
#include <vector>

void fastpattern(const umat &y, umat &pattern, uvec &group)  {
  unsigned n = y.n_rows;
  unsigned k = y.n_cols;  
  uvec mygroup(n);
  unsigned npattern = (unsigned) pow((double) 2,(double) k);
  umat mypattern(npattern,k);
  mypattern.fill(1);
  unsigned K=0;

  for (unsigned i=0; i<n; i++) {
    urowvec ri = y.row(i);
    bool found = FALSE;
    for (unsigned j=0; j<K; j++) {      
      if (sum(ri!=mypattern.row(j))==0) {
	found = true;
	mygroup(i) = j;	
	break;
      }
    }
    if (!found) {
      K++;
      mypattern.row(K-1) = ri;
      mygroup(i) = K-1;
    }
  } 
  pattern = mypattern.rows(0,K-1);
  group = mygroup;
}

SEXP FastPattern(SEXP y1,SEXP y2)  {
BEGIN_RCPP
  mat Y1 = Rcpp::as<mat>(y1);
  unsigned n = Y1.n_rows;
  unsigned k = Y1.n_cols;
  umat stat;
  if (!Rf_isNull(y2)) {
    mat Y2 = Rcpp::as<mat>(y2);   
    if ((Y2.n_rows!=n) || (Y2.n_cols!=k)) 
       throw(Rcpp::exception("Dimension did not agree","tools.cpp",1));
    stat = (Y1==Y2);
  } else {
    stat = conv_to<umat>::from(Y1);
  }
  uvec group(n);
  //  unsigned npattern = (unsigned) pow(2,(double) k);
  umat pattern; //(npattern,k);  
  fastpattern(stat,pattern,group);
  
  return(Rcpp::List::create(
			    Rcpp::Named("pattern")=pattern,
			    Rcpp::Named("group")=wrap(group)
			    ));
END_RCPP
} 



RcppExport SEXP FastApprox(const SEXP a,
			   const SEXP t,
			   const SEXP z) {
  NumericVector A(a);
  NumericVector T(t);
  NumericVector Z(z);
  vector<unsigned> idx(Z.size());
  vector<double> newT(Z.size());

  NumericVector::iterator it;  
  double lower,upper; int pos=200;
  for (int i=0; i<Z.size(); i++) {
    it = lower_bound(A.begin(), A.end(), Z[i]);
    if (it == A.begin()) { 
      pos = 0; 
      // upper = *it; // no smaller value  than val in vector
    } 
    else if (int(it-A.end())==0) {
      pos = A.size()-1;
      //lower = *(it-1); // no bigger value than val in vector
    } else {
      lower = *(it-1);
      upper = *it;
      pos = int(it- A.begin());
      if (abs(Z[i]-lower)<abs(Z[i]-upper)) {
	pos = int(it- A.begin());
      }
    }
    idx[i] = pos;
    newT[i] = T[pos];
  }
  List ans;
  ans["t"] = newT;
  ans["pos"] = idx;
  return(ans);
}
