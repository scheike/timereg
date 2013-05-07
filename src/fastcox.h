#ifndef FASTCOX_H
#define FASTCOX_H

#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;


RcppExport SEXP FastCoxPrep( SEXP entry, SEXP exit, SEXP status, SEXP X, SEXP id);
RcppExport SEXP FastCoxPL( SEXP b, SEXP x, SEXP xx, SEXP sgn, SEXP jumps);

#endif /* FASTCOX_H */
