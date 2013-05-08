#include "fastcox.h"

using namespace arma;

RcppExport SEXP FastCoxPrep( SEXP entry, SEXP exit, SEXP status, SEXP x, SEXP id) {
BEGIN_RCPP

  vec Exit = Rcpp::as<vec>(exit);  
  ivec Status = Rcpp::as<ivec>(status);
  mat X = Rcpp::as<mat>(x);
  bool haveId = (Rf_isNull)(id);
  bool Truncation = !((Rf_isNull)(entry));
  unsigned n = X.n_rows;
  uvec idx;

  mat XX(X.n_rows, X.n_cols*X.n_cols); // Calculate XX' at each time-point
  for (unsigned i=0; i<X.n_rows; i++) {
    rowvec Xi = X.row(i);
    XX.row(i) = reshape(Xi.t()*Xi,1,XX.n_cols);
  }

  if (Truncation) {
    vec Entry = Rcpp::as<vec>(entry);  
    Exit.insert_rows(0,Entry);
    XX.insert_rows(0,XX);
    X.insert_rows(0,X);
    Status.insert_rows(0,Status);
  } 
  idx = sort_index(Exit);
  ivec Sign;
  if (Truncation) {
    Sign.reshape(2*n,1); Sign.fill(1);
    for (unsigned i=0; i<n; i++) Sign(i) = -1;
    Status = Status%Sign;
    Sign = Sign.elem(idx);  
  }
  XX = XX.rows(idx);
  X = X.rows(idx);  
  Status = Status.elem(idx);
  uvec jumps = find(Status==1);
  
  uvec newId;
  if (!haveId) {
    uvec Id = Rcpp::as<uvec>(id);
    if (Truncation) {
      Id.insert_rows(0,Id);
    }
    newId = Id.elem(idx);
  }

  return(Rcpp::List::create(Rcpp::Named("XX")=XX,
			    Rcpp::Named("X")=X,
			    Rcpp::Named("jumps")=jumps,
			    Rcpp::Named("sign")=Sign,
			    Rcpp::Named("ord")=idx,
			    Rcpp::Named("time")=Exit,
			    Rcpp::Named("id")=newId
			    ));
END_RCPP
  }



RcppExport SEXP FastCoxPL( SEXP b, SEXP x, SEXP xx, SEXP sgn, SEXP jumps ) {
BEGIN_RCPP
  colvec beta = Rcpp::as<colvec>(b);
  mat X = Rcpp::as<mat>(x);
  mat XX = Rcpp::as<mat>(xx);
  uvec Jumps = Rcpp::as<uvec>(jumps);
  ivec Sign = Rcpp::as<ivec>(sgn);
  unsigned n = X.n_rows;

  colvec a = X*beta;
  colvec ea = exp(a);  
  if (Sign.n_cols==n) { // Truncation
    ea = ea%Sign;
  }

  colvec b1 = flipud(cumsum(flipud(ea)));
  // colvec b1(n); b1(0) = ea(n-1);
  // for (unsigned i=1; i<n; i++) {
  //   b1(i) = b1(i-1)+ea(n-i);
  // }
  vec val = a-log(b1); // Partial log-likelihood
  mat b = X;
  for (unsigned j=0; j<b.n_cols; j++) {
    b.col(j) = flipud(cumsum(flipud(b.col(j)%ea))) / b1;
  }
  mat grad = (X-b); // Score
  for (unsigned j=0; j<XX.n_cols; j++) {
    XX.col(j) = flipud(cumsum(flipud(XX.col(j)%ea))) / b1;
  }
  XX = XX.rows(Jumps);
  grad = grad.rows(Jumps);
  b = b.rows(Jumps);
  val = val.elem(Jumps);
  mat hess = -(reshape(sum(XX),X.n_cols,X.n_cols)-b.t()*b);

  return(Rcpp::List::create(Rcpp::Named("jumps")=jumps,
			    Rcpp::Named("ploglik")=sum(val),
			    Rcpp::Named("U")=grad,
			    Rcpp::Named("gradient")=sum(grad),
			    Rcpp::Named("hessian")=hess
			    ));
END_RCPP
  }
