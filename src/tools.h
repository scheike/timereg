#ifndef TOOLS_H
#define TOOLS_H

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstring>

using namespace std;
using namespace Rcpp;
using namespace arma;

RcppExport SEXP FastApprox(const SEXP a, const SEXP t, const SEXP z);
RcppExport SEXP FastPattern(SEXP y1,SEXP y2);

void fastpattern(const umat &y, umat &pattern, uvec &group);

template <class T>
string numStr(T x) {
  ostringstream nmbstr; nmbstr << x;
  string ss = nmbstr.str();    
  return ss;
}


#endif /* TOOLS_H */
