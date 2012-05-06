#ifndef MVN_H
#define MVN_H

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
/* #include <iostream> */
/* #include <sstream> */
/* #include <cstring> */


using namespace std;
using namespace Rcpp;
using namespace arma;

const double log2pi = log(2*math::pi());

struct vecmat
{
  vec V;
  mat M;
};

extern "C" double bvnd_(const double *dh, const double *dk, const double *r);
RcppExport SEXP cdf(SEXP a, SEXP b, SEXP r);
double Sbvn(double &l1, double &l2,double &r);
inline double Fbvn(double u1, double u2, double r) { 
  u1 *= -1; u2 *= -1; 
  return(Sbvn(u1,u2,r)); 
};
vecmat Dbvn(double y1, double y2, double R);
double dbvnorm(double y1, double y2, double R);
RcppExport SEXP FastApprox(const SEXP a, const SEXP t, const SEXP z);

/* template <class T> */
/* string numStr(T x) { */
/*   ostringstream nmbstr; nmbstr << x; */
/*   string ss = nmbstr.str();     */
/*   return ss; */
/* } */





#endif /* MVN_H */
