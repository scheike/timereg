// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <Rmath.h>

//#include "fastcox.h"
using namespace Rcpp;
using namespace arma;


RcppExport SEXP FastCoxPrep(SEXP EntrySEXP,
			    SEXP ExitSEXP,
			    SEXP StatusSEXP,
			    SEXP XSEXP,
			    SEXP IdSEXP,
			    SEXP TruncationSEXP) {
BEGIN_RCPP/*{{{*/
  arma::vec Entry = Rcpp::as<arma::vec>(EntrySEXP);
  arma::vec  Exit  = Rcpp::as<arma::vec>(ExitSEXP);
  arma::Col<int> Status= Rcpp::as<arma::Col<int> >(StatusSEXP);
  arma::mat  X     = Rcpp::as<arma::mat>(XSEXP);
  try {
    arma::Col<unsigned> Id    = Rcpp::as<arma::Col<unsigned> >(IdSEXP);
  }
  catch(...) {}

  //bool haveId = Rcpp::as<bool>(haveIdSEXP);
  bool Truncation = Rcpp::as<bool>(TruncationSEXP);
  // vec Exit = Rcpp::as<vec>(exit);  
  // ivec Status = Rcpp::as<ivec>(status);
  // mat X = Rcpp::as<mat>(x);
  // bool haveId = (Rf_isNull)(id);
  // bool Truncation = !((Rf_isNull)(entry));
  // bool Truncation = Entry.n_elem>0;
  // bool haveId = Id.n_elem>0;

  //  unsigned p = X.n_cols;
  unsigned n = Exit.n_elem;
  if (Truncation) n *= 2;

  //Rcout << "n=" << X.n_rows << ", p=" << X.n_cols << std::endl;

 mat XX(n, X.n_cols*X.n_cols); // Calculate XX' at each time-point
  for (unsigned i=0; i<X.n_rows; i++) {
    rowvec Xi = X.row(i);
    //    XX.row(i) = reshape(Xi.t()*Xi,1,XX.n_cols);
    XX.row(i) = vectorise(Xi.t()*Xi,1);
    if (Truncation) XX.row(i+n/2) = XX.row(i);
  }
  

  arma::Col<int> Sign;
  if (Truncation) {
    // vec Entry = Rcpp::as<vec>(entry);  
    Exit.insert_rows(0,Entry);
    X.insert_rows(0,X);
    Status.insert_rows(0,Status);
    Sign.reshape(n,1); Sign.fill(1);
    for (unsigned i=0; i<(n/2); i++) Sign(i) = -1;
    Status = Status%(1+Sign);
  }
  //Rcout << "Status=" << Status << std::endl;
  arma::uvec idx0 = sort_index(Status,"descend"); 
  arma::uvec idx = stable_sort_index(Exit.elem(idx0),"ascend");
  idx = idx0.elem(idx);
  //Rcout << "idx=" << idx << std::endl;
  if (Truncation) {
    Sign = Sign.elem(idx);  
  }
  if (X.n_rows>0) {
    XX = XX.rows(idx);
    X = X.rows(idx);  
  }
  Status = Status.elem(idx);
  arma::uvec jumps = find(Status>0);
  //Rprintf("jumps");
  arma::Col<unsigned> newId;
  // if (haveId) {
  //   // uvec Id = Rcpp::as<uvec>(id);
  //   if (Truncation) {
  //     Id.insert_rows(0,Id);
  //   }
  //   newId = Id.elem(idx);
  // }

  return(Rcpp::wrap(Rcpp::List::create(Rcpp::Named("XX")=XX,
				       Rcpp::Named("X")=X,
				       Rcpp::Named("jumps")=jumps,
				       Rcpp::Named("sign")=Sign,
				       Rcpp::Named("ord")=idx,
				       Rcpp::Named("time")=Exit,
				       Rcpp::Named("id")=newId				       
				       )));
END_RCPP
}/*}}}*/

RcppExport SEXP FastCoxPrepStrata(SEXP EntrySEXP,
			    SEXP ExitSEXP,
			    SEXP StatusSEXP,
			    SEXP XSEXP,
			    SEXP IdSEXP,
			    SEXP TruncationSEXP,
			    SEXP strataSEXP,
               		    SEXP weightsSEXP,
			    SEXP offsetsSEXP
			    ) {
BEGIN_RCPP/*{{{*/
  arma::vec Entry = Rcpp::as<arma::vec>(EntrySEXP);
  arma::vec  Exit  = Rcpp::as<arma::vec>(ExitSEXP);
  arma::Col<int> Status= Rcpp::as<arma::Col<int> >(StatusSEXP);
  arma::mat  X     = Rcpp::as<arma::mat>(XSEXP);
  arma::Col<int> strata= Rcpp::as<arma::Col<int> >(strataSEXP);
  try {
    arma::Col<unsigned> Id    = Rcpp::as<arma::Col<unsigned> >(IdSEXP);
  }
  catch(...) {}

  colvec weights = Rcpp::as<colvec>(weightsSEXP);
  colvec offsets = Rcpp::as<colvec>(offsetsSEXP);
  //bool haveId = Rcpp::as<bool>(haveIdSEXP);
  bool Truncation = Rcpp::as<bool>(TruncationSEXP);
  // vec Exit = Rcpp::as<vec>(exit);  
  // ivec Status = Rcpp::as<ivec>(status);
  // mat X = Rcpp::as<mat>(x);
  // bool haveId = (Rf_isNull)(id);
  // bool Truncation = !((Rf_isNull)(entry));
  // bool Truncation = Entry.n_elem>0;
  // bool haveId = Id.n_elem>0;

  //  unsigned p = X.n_cols;
  unsigned n = Exit.n_elem;
  if (Truncation) n *= 2;

  //Rcout << "n=" << X.n_rows << ", p=" << X.n_cols << std::endl;

 mat XX(n, X.n_cols*X.n_cols); // Calculate XX' at each time-point
  for (unsigned i=0; i<X.n_rows; i++) {
    rowvec Xi = X.row(i);
    //    XX.row(i) = reshape(Xi.t()*Xi,1,XX.n_cols);
    XX.row(i) = vectorise(Xi.t()*Xi,1);
    if (Truncation) XX.row(i+n/2) = XX.row(i);
  }
  

  arma::Col<int> Sign;
  if (Truncation) {
    // vec Entry = Rcpp::as<vec>(entry);  
    Exit.insert_rows(0,Entry);
    X.insert_rows(0,X);
    Status.insert_rows(0,Status);
    strata.insert_rows(0,strata);
//    weights.insert_rows(0,weights);
//    offsets.insert_rows(0,offsets);
    Sign.reshape(n,1); Sign.fill(1);
    for (unsigned i=0; i<(n/2); i++) Sign(i) = -1;
    Status = Status%(1+Sign);
  }
  //Rcout << "Status=" << Status << std::endl;
  arma::uvec idx0 = sort_index(Status,"descend"); 
  arma::uvec idx = stable_sort_index(Exit.elem(idx0),"ascend");
  idx = idx0.elem(idx);
  //Rcout << "idx=" << idx << std::endl;
  if (Truncation) {
    Sign = Sign.elem(idx);  
  }
  if (X.n_rows>0) {
    XX = XX.rows(idx);
    X = X.rows(idx);  
  }
  strata = strata.elem(idx); 
  weights = weights.elem(idx); 
  offsets = offsets.elem(idx); 
  Status = Status.elem(idx);
  arma::uvec jumps = find(Status>0);
  //Rprintf("jumps");
  arma::Col<unsigned> newId;
  // if (haveId) {
  //   // uvec Id = Rcpp::as<uvec>(id);
  //   if (Truncation) {
  //     Id.insert_rows(0,Id);
  //   }
  //   newId = Id.elem(idx);
  // }

  return(Rcpp::wrap(Rcpp::List::create(Rcpp::Named("XX")=XX,
				       Rcpp::Named("X")=X,
				       Rcpp::Named("jumps")=jumps,
				       Rcpp::Named("sign")=Sign,
				       Rcpp::Named("ord")=idx,
				       Rcpp::Named("time")=Exit,
				       Rcpp::Named("id")=newId,				       
				       Rcpp::Named("weights")=weights,
				       Rcpp::Named("offset")=offsets,				       
				       Rcpp::Named("strata")=strata				       
				       )));
END_RCPP
}/*}}}*/


mat  vecmatrow(const colvec &a, const mat b) {
  unsigned n = b.n_cols;
  mat res=b; 
  for (unsigned i=0; i<n; i++) {
    res.col(i)=a%b.col(i); 
  }  
  return(res);
} 


// colvec revcumsum(const colvec &a) {
//   return(flipud(cumsum(flipud(a))));
// }
colvec revcumsum(const colvec &a) {/*{{{*/
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
}/*}}}*/

RcppExport SEXP revcumsumR(SEXP ia) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  unsigned n = a.n_rows;
  colvec res = a; double prev=0;  
  for (unsigned i=0; i<n; i++) {
    prev += a(n-i-1);
    res(n-i-1) = prev;
  }  
List rres; 
rres["res"]=res; 

return(rres);
}/*}}}*/


RcppExport SEXP cumsumstrataR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); 
  tmpsum=tmpsum*0; 
  colvec res = a; 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(i); 
    tmpsum(ss) += a(i); 
    res(i) = tmpsum(ss);
  }  

  List rres; 
  rres["res"]=res; 
  return(rres);
} /*}}}*/

RcppExport SEXP revcumsumstrataR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); 
  tmpsum=tmpsum*0; 
  colvec res = a; 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(n-i-1); 
    tmpsum(ss) += a(n-i-1); 
    res(n-i-1) = tmpsum(ss);
  }  

  List rres; 
  rres["res"]=res; 
  return(rres);
} /*}}}*/

colvec revcumsumstrata(const colvec &a,IntegerVector strata,int nstrata) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); 
  tmpsum=tmpsum*0; 
  colvec res = a; 

  for (unsigned i=0; i<n; i++) {
    int ss=strata(n-i-1); 
    tmpsum(ss) += a(n-i-1); 
    res(n-i-1) = tmpsum(ss);
//    printf("%d %d %d %lf %lf \n",i,ss,strata(n-i-1),tmpsum(ss),a(n-i-1)); 
  }  
//  printf("===========================\n"); 
  return(res);
}
colvec revcumsumstrata1(const colvec &a, const colvec &v1, const colvec &v2,IntegerVector strata,int nstrata) {
//  unsigned n = a.n_rows;
//  colvec dd= revcumsumstrata(a%v1,strata,nstrata); 
//  for (unsigned i=0; i<n; i++) {
//	  if (v2(i)>0) dd(i)=dd(i)/v2(i); else dd(i)=0; 
//  }  
//  return(dd); 
  return(revcumsumstrata(a%v1,strata,nstrata)/v2);
}/*}}}*/


RcppExport SEXP FastCoxPL(SEXP betaSEXP,
			  SEXP XSEXP,
			  SEXP XXSEXP,
			  SEXP SignSEXP,
			  SEXP JumpsSEXP) {
BEGIN_RCPP/*{{{*/
  colvec beta = Rcpp::as<colvec>(betaSEXP);
  mat X = Rcpp::as<mat>(XSEXP);
  mat XX = Rcpp::as<mat>(XXSEXP);
  arma::uvec Jumps = Rcpp::as<uvec >(JumpsSEXP);
  arma::Col<int> Sign = Rcpp::as<arma::Col<int> >(SignSEXP);
  // unsigned n = X.n_rows;
  unsigned p = X.n_cols;

  colvec Xb = X*beta;
  colvec eXb = exp(Xb);
  if (Sign.n_rows==eXb.n_rows) { // Truncation
    eXb = Sign%eXb;
  }

  colvec S0 = revcumsum(eXb);
  //Rcout << "S0\n" << S0;

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

  E = E.rows(Jumps);
  mat E2(E.n_rows, E.n_cols*E.n_cols); // Calculate E' E at each time-point
  for (unsigned i=0; i<E.n_rows; i++) {
    rowvec Xi = E.row(i);
    E2.row(i) = vectorise(Xi.t()*Xi,1);
  }

  mat XX2 = XX;
  for (unsigned j=0; j<XX2.n_cols; j++) { // int S2/S0(s)
    XX2.col(j) = revcumsum(XX2.col(j),eXb,S0);
  }

  XX2 = XX2.rows(Jumps);
  //  X = X.rows(Jumps);
  S0 = S0.elem(Jumps);
  mat grad = (X.rows(Jumps)-E); // Score
  vec val = Xb.elem(Jumps)-log(S0); // Partial log-likelihood
  mat hesst = -(XX2-E2);
  mat hess = reshape(sum(hesst),p,p);

//  hess.print("hess"); 
//  S0.print("S0"); 
//  grad.print("grad"); 
//  E.print("E"); 
//  XX2.print("XX"); 
//  printf("============================ \n"); 


  return(Rcpp::List::create(Rcpp::Named("jumps")=Jumps,
			    Rcpp::Named("ploglik")=sum(val),
			    Rcpp::Named("U")=grad,
			    Rcpp::Named("gradient")=sum(grad),
			    Rcpp::Named("hessian")=hess,
			    Rcpp::Named("hessianttime")=hesst,
			    Rcpp::Named("S2S0")=XX,
			    Rcpp::Named("E")=E,
			    Rcpp::Named("S0")=S0
			    ));
END_RCPP
  }/*}}}*/

RcppExport SEXP FastCoxPLstrata(SEXP betaSEXP,
			  SEXP XSEXP,
			  SEXP XXSEXP,
			  SEXP SignSEXP,
			  SEXP JumpsSEXP, 
			  SEXP strataSEXP, 
			  SEXP nstrataSEXP,
			  SEXP weightsSEXP,
			  SEXP offsetsSEXP
			  ) {
BEGIN_RCPP/*{{{*/
  colvec beta = Rcpp::as<colvec>(betaSEXP);
  mat X = Rcpp::as<mat>(XSEXP);
  mat XX = Rcpp::as<mat>(XXSEXP);
  arma::uvec Jumps = Rcpp::as<uvec >(JumpsSEXP);
  arma::Col<int> Sign = Rcpp::as<arma::Col<int> >(SignSEXP);
  IntegerVector strata(strataSEXP);
  int nstrata = Rcpp::as<int>(nstrataSEXP);
  // unsigned n = X.n_rows;
  unsigned p = X.n_cols;
  colvec weights = Rcpp::as<colvec>(weightsSEXP);
  colvec offsets = Rcpp::as<colvec>(offsetsSEXP);


  colvec Xb = X*beta;
  colvec eXb = exp(Xb)%weights%offsets;
//  colvec eXb = exp(Xb);
  if (Sign.n_rows==eXb.n_rows) { // Truncation
    eXb = Sign%eXb;
  }

  colvec S0 = revcumsumstrata(eXb,strata,nstrata);

  //Rcout << "S0\n" << S0;

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
    E.col(j) = revcumsumstrata1(X.col(j),eXb,S0,strata,nstrata);
  }

  E = E.rows(Jumps);
  mat E2(E.n_rows, E.n_cols*E.n_cols); // Calculate E' E at each time-point
  for (unsigned i=0; i<E.n_rows; i++) {
    rowvec Xi = E.row(i);
    E2.row(i) = vectorise(Xi.t()*Xi,1);
  }

  mat XX2 = XX;
  for (unsigned j=0; j<XX2.n_cols; j++) { // int S2/S0(s)
    XX2.col(j) = revcumsumstrata1(XX2.col(j),eXb,S0,strata,nstrata);
  }

  XX2 = XX2.rows(Jumps);
  //  X = X.rows(Jumps);
  colvec weightsJ=weights.elem(Jumps);  
  S0 = S0.elem(Jumps);
  mat grad = (X.rows(Jumps)-E); // Score
  vec val =  (Xb.elem(Jumps)-log(S0)); // Partial log-likelihood
  colvec S02 = weightsJ%S0;
  mat grad2=vecmatrow(weightsJ,grad);   // score 
  vec val2 = weightsJ%(Xb.elem(Jumps)-log(S0)); // Partial log-likelihood
  mat hesst = -(XX2-E2);

  mat hesst2 = -vecmatrow(weightsJ,XX2-E2);
  mat hess  = reshape(sum(hesst),p,p);
  mat hess2 = reshape(sum(hesst2),p,p);

  if (hess.has_nan()) {
	printf("============================ \n"); 
	S0.print("S0"); 
	grad.print("grad"); 
	E.print("E"); 
	XX2.print("XX"); 
	printf("============================ \n"); 
	}

  return(Rcpp::List::create(Rcpp::Named("jumps")=Jumps,
			    Rcpp::Named("ploglik")=sum(val2),
			    Rcpp::Named("U")=grad2,
			    Rcpp::Named("gradient")=sum(grad2),
			    Rcpp::Named("hessian")=hess2,
			    Rcpp::Named("hessianttime")=hesst,
			    Rcpp::Named("S2S0")=XX,
			    Rcpp::Named("E")=E,
			    Rcpp::Named("S0")=S0
			    ));
END_RCPP
  }/*}}}*/


RcppExport SEXP CubeVec(SEXP XXSEXP, SEXP betaSEXP)
		  {
BEGIN_RCPP/*{{{*/
  colvec beta = Rcpp::as<colvec>(betaSEXP);
  mat XX = Rcpp::as<mat>(XXSEXP);
  unsigned p = beta.n_rows;
  unsigned n = XX.n_rows;

//  XX.print("XX"); 

  mat XXbeta(n,p);
  for (unsigned j=0; j<n; j++)  {
	  XXbeta.row(j)=(reshape(XX.row(j),p,p)*beta).t();
  }

  return(Rcpp::List::create(Rcpp::Named("XXbeta")=XXbeta));
END_RCPP
}/*}}}*/

RcppExport SEXP CubeMat( SEXP XXSEXP, SEXP XSEXP)
		  {
BEGIN_RCPP/*{{{*/
  mat XX = Rcpp::as<mat>(XXSEXP);
  mat X  = Rcpp::as<mat>(XSEXP);
  unsigned p = X.n_cols;
  unsigned n = XX.n_rows;

  mat XXX(n,p*p);
  for (unsigned j=0; j<n; j++)  {
	  XXX.row(j)=vectorise(reshape(XX.row(j),p,p)*X,1);
  }

  return(Rcpp::List::create(Rcpp::Named("XXX")=XXX));
END_RCPP
}/*}}}*/


RcppExport SEXP PropTestCox(SEXP iU, SEXP idUt, SEXP insim, SEXP iobssup
			  ) {
BEGIN_RCPP/*{{{*/
  mat U      = Rcpp::as<mat>(iU);
  mat dUt = Rcpp::as<mat>(idUt);
  arma::vec osup = Rcpp::as<arma::vec>(iobssup);
  int nsim = Rcpp::as<int>(insim);
  unsigned p = U.n_cols;
  unsigned n = U.n_rows;

  vec pval(p); 
  mat Uthati(n,p); 
  mat Uti(n,p); 
  mat sup(nsim,p); 
  mat simUti(n,50*p); 


  GetRNGstate();  /* to use R random normals */

  for (unsigned j=0; j<nsim; j++) {
     int thissiml=0; 
     for (unsigned k=0; k<n; k++)  Uti.row(k)=norm_rand()*U.row(k); 
     for (unsigned k=0; k<p; k++)  Uti.col(k) = cumsum(Uti.col(k));

     for (unsigned k=0; k<n; k++)  {
	  Uthati.row(k)=(reshape(dUt.row(k),p,p)*(Uti.row(n-1)).t()).t();
     }
     Uthati=Uti-Uthati; 

//     if(j==0) Uthati.print("one sim"); 

     for (unsigned k=0;k<p;k++)  {
//	     printf(" %lf \n",sup(j,k)); 
        sup(j,k)=max(abs(Uthati.col(k))); 
//      count if sup for this realization is larger than supObs
        if (sup(j,k)>=osup(k) & thissiml==0) { pval(k)++; thissiml=1;}
        if (j<50) { simUti.col(j*p+k)=Uthati.col(k); }
     }
  }
  pval=pval/nsim; 

  PutRNGstate();  /* to use R random normals */

  return(Rcpp::List::create(Rcpp::Named("supUsim")=sup,
			    Rcpp::Named("simUt")=simUti,
			    Rcpp::Named("pval")=pval)); 
END_RCPP
  }/*}}}*/


