#include "fastcox.h"

using namespace arma;

RcppExport SEXP FastCoxPrep( SEXP entry, SEXP exit, SEXP status, SEXP x, SEXP id) {
BEGIN_RCPP
  vec Entry = Rcpp::as<vec>(entry);
  vec Exit = Rcpp::as<vec>(exit);  
  uvec Status = Rcpp::as<uvec>(status);
  mat X = Rcpp::as<mat>(x);
  bool haveId = (Rf_isNull)(id);

  unsigned n = X.n_rows;
  Entry.insert_rows(0,Exit);
  uvec idx = sort_index(Entry);

  mat XX(X.n_rows, X.n_cols*X.n_cols); // Calculate XX' at each time-point
  for (unsigned i=0; i<X.n_rows; i++) {
    rowvec Xi = X.row(i);
    XX.row(i) = reshape(Xi.t()*Xi,1,XX.n_cols);
  }
  XX.insert_rows(0,XX);
  XX = XX.rows(idx);
  X.insert_rows(0,X);
  X = X.rows(idx);
  
  uvec newId;
  if (!haveId) {
    uvec Id = Rcpp::as<uvec>(id);
    Id.insert_rows(0,Id);
    newId = Id.elem(idx);
  }

  Status.insert_rows(0,Status);
  ivec Sign(2*n); Sign.fill(-1);
  for (unsigned i=0; i<n; i++) Sign(i) = 1; 
  Sign = Sign.elem(idx);
  Status = Status.elem(idx);
  uvec jumps = find(Status%Sign==1);

  return(Rcpp::List::create(Rcpp::Named("XX")=XX,
			    Rcpp::Named("X")=X,
			    Rcpp::Named("jumps")=jumps,
			    Rcpp::Named("sign")=Sign,
			    Rcpp::Named("ord")=idx,
			    Rcpp::Named("time")=Entry,
			    Rcpp::Named("id")=newId
			    ));
END_RCPP
  }


RcppExport SEXP FastCoxPL( SEXP b, SEXP x, SEXP xx, SEXP sgn, SEXP jumps ) {
BEGIN_RCPP
  colvec beta = Rcpp::as<colvec>(b);
  mat X = Rcpp::as<mat>(x);
  mat XX = Rcpp::as<mat>(xx);
  ivec Sign = Rcpp::as<ivec>(sgn);
  uvec Jumps = Rcpp::as<uvec>(jumps);

  colvec a = X*beta;
  colvec ea = Sign%exp(a);
  colvec b1 = flipud(cumsum(flipud(ea)));
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
