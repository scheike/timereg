#include "mvn.h"

RcppExport SEXP cdf(SEXP a, SEXP b, SEXP r) { 
  double u1=-Rcpp::as<double>(a);
  double u2=-Rcpp::as<double>(b);
  double cr=Rcpp::as<double>(r);  
  double val = bvnd_(&u1, &u2, &cr);
  NumericVector res(1); res[0] = val;
  return(res);
}

double dbvnorm(double y1, double y2, double R) {
  double detR = 1-R*R;
  // inv(R) = [1 -r; -r 1]/detR (prove by gauss elim.)
  double res = 1/(2*M_PI*sqrt(detR))*exp(-0.5/detR*(y1*y1+y2*y2-2*R*y1*y2));
  return(res);
}

vecmat Dbvn(double y1, double y2, double R) {
  vec DP(2);
  double R2 = R*R;
  DP(0) = Rf_dnorm4(y1,0.0,1.0,0)*Rf_pnorm5(y2-R*y1,0.0,sqrt(1-R2),1,0);
  DP(1) = Rf_dnorm4(y2,0.0,1.0,0)*Rf_pnorm5(y1-R*y2,0.0,sqrt(1-R2),1,0);
  mat HP(2,2);
  HP(1,0) = HP(0,1) = dbvnorm(y1,y2,R);
  HP(0,0) = -y1*DP(0) - R*HP(1,0);
  HP(1,1) = -y2*DP(1) - R*HP(1,0);  
  vecmat res;
  res.V = DP;
  res.M= HP;
  return(res);
}

double Sbvn(double &l1, double &l2,double &r) {     
  double val = bvnd_(&l1, &l2, &r);
  return(val);
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
