#include "fastcox.h"

using namespace arma;

RcppExport SEXP FastCoxPrep( SEXP entry, SEXP exit, SEXP status, SEXP x, SEXP id) {
BEGIN_RCPP

  vec Exit = Rcpp::as<vec>(exit);  
  ivec Status = Rcpp::as<ivec>(status);
  mat X = Rcpp::as<mat>(x);
  bool haveId = (Rf_isNull)(id);
  bool Truncation = !((Rf_isNull)(entry));
  unsigned n = Exit.n_elem;

  mat XX(X.n_rows, X.n_cols*X.n_cols); // Calculate XX' at each time-point
  for (unsigned i=0; i<X.n_rows; i++) {
    rowvec Xi = X.row(i);
    XX.row(i) = reshape(Xi.t()*Xi,1,XX.n_cols);
  }

  ivec Sign;
  if (Truncation) {
    vec Entry = Rcpp::as<vec>(entry);  
    Exit.insert_rows(0,Entry);
    XX.insert_rows(0,XX);
    X.insert_rows(0,X);
    Status.insert_rows(0,Status);
    Sign.reshape(2*n,1); Sign.fill(1);
    for (unsigned i=0; i<n; i++) Sign(i) = -1;
    Status = Status%(1+Sign);
  }
  uvec idx0 = sort_index(Status,1); 
  uvec idx = stable_sort_index(Exit.elem(idx0),0);
  idx = idx0.elem(idx);
  if (Truncation) {
    Sign = Sign.elem(idx);  
  }

  if (X.n_rows>0) {
    XX = XX.rows(idx);
    X = X.rows(idx);  
  }
  Status = Status.elem(idx);
  uvec jumps = find(Status>0);
  
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


// colvec revcumsum(const colvec &a) {
//   return(flipud(cumsum(flipud(a))));
// }
colvec revcumsum(const colvec &a) {
  unsigned n = a.n_rows;
  colvec res = a; double prev=0;  
  for (unsigned i=0; i<n; i++) {
    prev += a(n-i-1);
    res(n-i-1) = prev;
  }  
  return(res);
}
colvec revcumsum(const colvec &a, const colvec &v1, const colvec &v2) {
  return(revcumsum(a%v1)/v2);
}


RcppExport SEXP FastCoxPL( SEXP b, SEXP x, SEXP xx, SEXP sgn, SEXP jumps ) {
BEGIN_RCPP
  colvec beta = Rcpp::as<colvec>(b);
  mat X = Rcpp::as<mat>(x);
  mat XX = Rcpp::as<mat>(xx);
  uvec Jumps = Rcpp::as<uvec>(jumps);
  ivec Sign = Rcpp::as<ivec>(sgn);
  unsigned n = X.n_rows;

  colvec Xb = X*beta;
  colvec eXb = exp(Xb);  
  if (Sign.n_rows==eXb.n_rows) { // Truncation
    eXb = Sign%eXb;
  }

  colvec S0 = revcumsum(eXb);
  // mat S1(X.n_rows,X.n_cols);
  // for (unsigned j=0; j<X.n_cols; j++) {
  //   S1.col(j) = revcumsum(X.col(j),eXb);
  // }
  // mat S2(X.n_rows,XX.n_cols);
  // for (unsigned j=0; j<X.n_cols; j++) {
  //   S2.col(j) = revcumsum(XX.col(j),eXb);
  // }
  mat E(X.n_rows,X.n_cols); // S1/S0(s)
  for (unsigned j=0; j<X.n_cols; j++) {
    E.col(j) = revcumsum(X.col(j),eXb,S0);
  }
  for (unsigned j=0; j<XX.n_cols; j++) { // int S2/S0(s)
    XX.col(j) = revcumsum(XX.col(j),eXb,S0);
  }

  XX = XX.rows(Jumps);
  X = X.rows(Jumps);
  E = E.rows(Jumps);
  S0 = S0.elem(Jumps);
  mat grad = (X-E); // Score
  vec val = Xb.elem(Jumps)-log(S0); // Partial log-likelihood
  mat hess = -(reshape(sum(XX),X.n_cols,X.n_cols)-E.t()*E);

  return(Rcpp::List::create(Rcpp::Named("jumps")=jumps,
			    Rcpp::Named("ploglik")=sum(val),
			    Rcpp::Named("U")=grad,
			    Rcpp::Named("gradient")=sum(grad),
			    Rcpp::Named("hessian")=hess,
			    Rcpp::Named("S2S0")=XX,
			    Rcpp::Named("E")=E,
			    Rcpp::Named("S0")=S0
			    ));
END_RCPP
  }
