// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <Rmath.h>

//#include "fastcox.h"
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List FastCoxPrep( arma::vec Entry,
		  arma::vec Exit,
		  arma::ivec Status,
		  arma::mat X,
		  arma::uvec Id,
		  bool haveId,
		  bool Truncation) {
  // vec Exit = Rcpp::as<vec>(exit);  
  // ivec Status = Rcpp::as<ivec>(status);
  // mat X = Rcpp::as<mat>(x);
  // bool haveId = (Rf_isNull)(id);
  // bool Truncation = !((Rf_isNull)(entry));
  // bool Truncation = Entry.n_elem>0;
  // bool haveId = Id.n_elem>0;

  unsigned p = X.n_cols;
  unsigned n = Exit.n_elem;
  if (Truncation) n *= 2;

 Rcout << "Hallo" << std::endl;
 Rcout << "n=" << X.n_rows << ", p=" << X.n_cols << std::endl;

 mat XX(n, X.n_cols*X.n_cols); // Calculate XX' at each time-point
  for (unsigned i=0; i<X.n_rows; i++) {
    rowvec Xi = X.row(i);
    //    XX.row(i) = reshape(Xi.t()*Xi,1,XX.n_cols);
    XX.row(i) = vectorise(Xi.t()*Xi,1);
    if (Truncation) XX.row(i+n/2) = XX.row(i);
  }
  

 Rprintf("Hallo");
  ivec Sign;
  if (Truncation) {
    // vec Entry = Rcpp::as<vec>(entry);  
    Exit.insert_rows(0,Entry);
    Rcout << "X\n";
    X.insert_rows(0,X);
    Rcout << "Status\n";
    Status.insert_rows(0,Status);
    Sign.reshape(n,1); Sign.fill(1);
    for (unsigned i=0; i<n; i++) Sign(i) = -1;
    Status = Status%(1+Sign);
  }
 Rprintf("Hallo");
  uvec idx0 = sort_index(Status,1); 
  uvec idx = stable_sort_index(Exit.elem(idx0),0);
  idx = idx0.elem(idx);
  if (Truncation) {
    Sign = Sign.elem(idx);  
  }
 Rprintf("Hallo");
  if (X.n_rows>0) {
    XX = XX.rows(idx);
    X = X.rows(idx);  
  }
  Status = Status.elem(idx);
  uvec jumps = find(Status>0);
   Rprintf("Hallo");
  uvec newId;
  // if (haveId) {
  //   // uvec Id = Rcpp::as<uvec>(id);
  //   if (Truncation) {
  //     Id.insert_rows(0,Id);
  //   }
  //   newId = Id.elem(idx);
  // }

  return(Rcpp::List::create(Rcpp::Named("XX")=XX,
			    Rcpp::Named("X")=X,
			    Rcpp::Named("jumps")=jumps,
			    Rcpp::Named("sign")=Sign,
			    Rcpp::Named("ord")=idx,
			    Rcpp::Named("time")=Exit,
			    Rcpp::Named("id")=newId
			    ));
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


// [[Rcpp::export]]
List FastCoxPL(const arma::colvec& beta, 
	       const arma::colvec& X, 
	       const arma::colvec& XX, 
	       const arma::ivec& Sign, 
	       const arma::uvec& Jumps) {
BEGIN_RCPP
  // colvec beta = Rcpp::as<colvec>(b);
  // mat X = Rcpp::as<mat>(x);
  // mat XX = Rcpp::as<mat>(xx);
  // uvec Jumps = Rcpp::as<uvec>(jumps);
  // ivec Sign = Rcpp::as<ivec>(sgn);
  // unsigned n = X.n_rows;
  unsigned p = X.n_cols;

  colvec Xb = X*beta;
  colvec eXb = exp(Xb);  
  if (Sign.n_rows==eXb.n_rows) { // Truncation
    eXb = Sign%eXb;
  }

  colvec S0 = revcumsum(eXb);
  // mat S1(X.n_rows,p);
  // for (unsigned j=0; j<X.n_cols; j++) {
  //   S1.col(j) = revcumsum(X.col(j),eXb);
  // }
  // mat S2(X.n_rows,XX.n_cols);
  // for (unsigned j=0; j<p; j++) {
  //   S2.col(j) = revcumsum(XX.col(j),eXb);
  // }
  mat E(X.n_rows,p); // S1/S0(s)
  for (unsigned j=0; j<p; j++) {
    E.col(j) = revcumsum(X.col(j),eXb,S0);
  }
  mat XX2 = XX;
  for (unsigned j=0; j<XX2.n_cols; j++) { // int S2/S0(s)
    XX2.col(j) = revcumsum(XX2.col(j),eXb,S0);
  }

  XX2 = XX2.rows(Jumps);
  //  X = X.rows(Jumps);
  E = E.rows(Jumps);
  S0 = S0.elem(Jumps);
  mat grad = (X.rows(Jumps)-E); // Score
  vec val = Xb.elem(Jumps)-log(S0); // Partial log-likelihood
  mat hess = -(reshape(sum(XX2),p,p)-E.t()*E);

  return(Rcpp::List::create(Rcpp::Named("jumps")=Jumps,
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
