// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include "twostage.h"

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
			    SEXP offsetsSEXP, 
			    SEXP ZSEXP
			    ) {
BEGIN_RCPP/*{{{*/
  arma::vec Entry = Rcpp::as<arma::vec>(EntrySEXP);
  arma::vec  Exit  = Rcpp::as<arma::vec>(ExitSEXP);
  arma::Col<int> Status= Rcpp::as<arma::Col<int> >(StatusSEXP);
  arma::mat  X     = Rcpp::as<arma::mat>(XSEXP);
  arma::mat  Z     = Rcpp::as<arma::mat>(ZSEXP);
  arma::Col<int> strata= Rcpp::as<arma::Col<int> >(strataSEXP);
  arma::Col<int> Id= Rcpp::as<arma::Col<int> >(IdSEXP);
//  arma::uvec idx = Rcpp::as<arma::Col<int> >(sortid);
//  try {
//    arma::Col<unsigned> Id    = Rcpp::as<arma::Col<unsigned> >(IdSEXP);
//  }
//  catch(...) {}


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

  unsigned nZ = Z.n_rows;
  if (Truncation) nZ = 2*Z.n_rows;
  mat ZX(nZ , Z.n_cols * X.n_cols);
  if (Z.n_rows==X.n_rows) 
  for (unsigned i=0; i<X.n_rows; i++) {
    rowvec Xi = X.row(i);
    rowvec Zi = Z.row(i);
    ZX.row(i) = vectorise((Xi.t()*Zi),1); // to get back to right form with reshape
//    if (i==-1) {
//       rowvec zx=ZX.row(i) ; 
//       Xi.print("xi"); 
//       Zi.print("zi"); 
//       zx.print("zx"); 
//       mat mm=reshape(zx,Z.n_cols,X.n_cols);
//       mm.print("mm"); 
//    }
    if (Truncation) ZX.row(i+n/2) = ZX.row(i);
  }


  arma::Col<int> Sign;
  Sign.reshape(n,1); Sign.fill(1);
  if (Truncation) {
    // vec Entry = Rcpp::as<vec>(entry);  
    Exit.insert_rows(0,Entry);
    X.insert_rows(0,X);
    Status.insert_rows(0,Status);
    Id.insert_rows(0,Id);
    strata.insert_rows(0,strata);
    weights.insert_rows(0,weights);
    offsets.insert_rows(0,offsets);
    for (unsigned i=0; i<(n/2); i++) Sign(i) = -1;
    Status = Status%(1+Sign);
  }
  //Rcout << "Status=" << Status << std::endl;
  
  // also sorting after id to use multiple phregs together
  // ts 20/3-2018
  arma::uvec idx00 = sort_index(Id,"ascend"); 
  arma::uvec idx0 = stable_sort_index(Status.elem(idx00),"descend"); 
  idx0 = idx00.elem(idx0);
  arma::uvec idx = stable_sort_index(Exit.elem(idx0),"ascend");
  idx = idx0.elem(idx);

//  arma::uvec idx0 = stable_sort_index(Status.elem(idx00),"descend"); 

//  arma::uvec idx0 = sort_index(Status,"descend"); 
//  arma::uvec idx = stable_sort_index(Exit.elem(idx0),"ascend");
//  idx = idx0.elem(idx);

//  arma::uvec idx00 = stable_sort_index(Id.elem(idx),"ascend"); 
//  idx = idx00.elem(idx);


  //Rcout << "idx=" << idx << std::endl;
  if (Truncation) {
    Sign = Sign.elem(idx);  
  }
  if (X.n_rows>0) {
    XX = XX.rows(idx);
    X = X.rows(idx);  
  }
  if ((ZX.n_rows==XX.n_rows) & (XX.n_rows>0)) {
    ZX = ZX.rows(idx);  
  }
  Exit = Exit.elem(idx); 
  weights = weights.elem(idx); 
  offsets = offsets.elem(idx); 
  Status = Status.elem(idx);
  Id = Id.elem(idx); 
  strata = strata.elem(idx); 
  arma::uvec jumps = find(Status>0);
  //Rprintf("jumps");
//  arma::Col<unsigned> newId;
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
				       Rcpp::Named("id")=Id,				       
				       Rcpp::Named("weights")=weights,
				       Rcpp::Named("offset")=offsets,				       
				       Rcpp::Named("strata")=strata,				       
				       Rcpp::Named("ZX")=ZX				       
				       )));
END_RCPP
}/*}}}*/

colvec  whichi(IntegerVector a,int n, int j) {/*{{{*/
  colvec res(n); 
  for (int i=0; i<n; i++) {
	  res(i)=(a(i)==j); 
  }  
  return(res);
} /*}}}*/

mat  vecmatrow(const colvec &a, const mat &b) {/*{{{*/
  unsigned n = b.n_cols;
  mat res=b; 
  for (unsigned i=0; i<n; i++) {
    res.col(i)=a%b.col(i); 
  }  
  return(res);
} /*}}}*/

// colvec revcumsum(const colvec &a) {
//   return(flipud(cumsum(flipud(a))));
// }
RcppExport SEXP revcumsumR(SEXP ia) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  unsigned n = a.n_rows;
  colvec res = a; 
  double prev=0;  
  for (unsigned i=0; i<n; i++) {
    prev += a(n-i-1);
    res(n-i-1) = prev;
  }  
List rres; 
rres["res"]=res; 

return(rres);
}/*}}}*/

colvec revcumsum(const colvec &a) {/*{{{*/
  unsigned n = a.n_rows;
  colvec res = a; 
  double prev=0;  
  for (unsigned i=0; i<n; i++) {
    prev += a(n-i-1);
    res(n-i-1) = prev;
  }  
  return(res);
}
colvec revcumsum(const colvec &a, const colvec &v1, const colvec &v2) {
  return(revcumsum(a%v1)/v2);
}/*}}}*/

RcppExport SEXP sumstrataR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); tmpsum.zeros(); 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(i); 
    if ((ss< nstrata) & (ss>=0))
    tmpsum(ss) += a(i); 
  }  
  List rres; 
  rres["res"]=tmpsum; 
  return(rres);
} /*}}}*/

colvec  sumstrata(colvec a,IntegerVector strata,int nstrata) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); 
  tmpsum.zeros(); tmpsum.zeros(); 

  for (unsigned i=0; i<n; i++) {
    int ss=strata(i); 
    if ((ss< nstrata) & (ss>=0))
    tmpsum(ss) += a(i); 
  }  
  return(tmpsum);
} /*}}}*/

RcppExport SEXP cumsumstrataR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec res = a; 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(i); 
    if ((ss<nstrata) & (ss>=0))  {
       tmpsum(ss) += a(i); 
       res(i) = tmpsum(ss);
    }
  }  
  List rres; 
  rres["res"]=res; 
  return(rres);
} /*}}}*/

colvec  cumsumstrata(colvec a,IntegerVector strata,int nstrata) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); 
  tmpsum.zeros(); tmpsum.zeros(); 
  colvec res = a; 

  for (unsigned i=0; i<n; i++) {
    int ss=strata(i); 
    if ((ss<nstrata) & (ss>=0))  {
    tmpsum(ss) += a(i); 
    res(i) = tmpsum(ss);
    }
  }  

  return(res);
} /*}}}*/

colvec  cumsumstrataPO(colvec a,IntegerVector strata,int nstrata,double propodds,colvec exb) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); 
  tmpsum.zeros(); 
  colvec res = a; 
  colvec pow = a; 

  for (unsigned i=0; i<n; i++) {
    int ss=strata(i); 
    if ((ss<nstrata) & (ss>=0))  {
    if (propodds>0)  pow(i)=(1+propodds*exb(i)*tmpsum(ss)); 
    tmpsum(ss) += pow(i)/a(i); 
    res(i) = tmpsum(ss);
    }
  }  

  return(pow);
} /*}}}*/


mat DLambeta(colvec weights,colvec S0,mat E,mat Xi,IntegerVector strata,int nstrata,double propodds,colvec exb) {/*{{{*/
  unsigned n = S0.n_rows;
  unsigned p = E.n_cols;
  colvec tmpsum(nstrata); tmpsum.zeros(); 
  mat dLbetatminus(nstrata,p); dLbetatminus.zeros(); 
  colvec res = S0; 
  colvec pow = S0; 
  mat dLbeta(n,p); dLbeta.zeros(); 

  for (unsigned i=0; i<n; i++) {
    int ss=strata(i); 
//    if (propodds>0)  
    pow(i)=(1+propodds*exb(i)*tmpsum(ss)); 
    dLbeta.row(i) = dLbetatminus.row(ss)+weights(i)*
   ( (dLbetatminus.row(ss)*exb(i)+Xi.row(i)*(pow(i)-1))/S0(i)-E.row(i)*pow(i)/S0(i));
    tmpsum(ss) += weights(i)*pow(i)/S0(i); 
    res(i) = tmpsum(ss);
    dLbetatminus.row(ss) = dLbeta.row(i); 
  }  

  return(dLbeta);
} /*}}}*/

colvec  cumsumstrataAddGam(colvec a,IntegerVector strata,int nstrata,
		colvec exb,colvec etheta,cube thetades,cube rv,mat ags,
		uvec Jumps) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); 
  tmpsum.zeros(); tmpsum.zeros(); 
  colvec res = a; 
  colvec pow = a; 
  vec allvec(6); 
  vec DthetaS(etheta.n_elem),DthetaDtS(etheta.n_elem); 
  colvec exbs(2); 

  for (unsigned i=0; i<n; i++) {
    int ss=strata(i); 
    mat thetadesv=thetades.slice(i); 
    mat rv1=rv.slice(i); 
    // ordered after time and comes two and two 
    if (strata(i)==0) { exbs(0)=exb(Jumps(i)); exbs(1)=exb(Jumps(i)+1); }
    if (strata(i)==1) { exbs(1)=exb(Jumps(i)); exbs(0)=exb(Jumps(i)+1); }
//    exbs(0)=exb(Jumps(i)); exbs(1)=exb(0); 
//    exbs(1)=exb(i); 
//    printf(" %d %d \n",(int) JumpsCauses(i,0),JumpsCauses(i,1)); 
//    exbs.print("exbs"); 
    double ll=survivalRVCmarg(etheta,thetadesv,ags,strata(i)+1,exbs%tmpsum,rv1,DthetaS,DthetaDtS,allvec);
//    Rprintf(" %d %d %lf %lf %lf \n",i,strata(i),ll,allvec(0),1/a(i)); 
//    etheta.print("etheta"); 
//    thetadesv.print("thetades"); 
//    tmpsum.print("tmpsum"); 
//    rv1.print("rv1"); 
    pow(i)=allvec(0)/ll; //   S / D_1 S
    if ((ss<nstrata) & (ss>=0))  {
       tmpsum(ss) += pow(i)/a(i); 
       res(i) = tmpsum(ss);
    }
  }  

  return(pow);
} /*}}}*/

RcppExport SEXP revcumsumstrataR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec res = a; 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(n-i-1); 
    if ((ss<nstrata) & (ss>=0))  {
       tmpsum(ss) += a(n-i-1); 
       res(n-i-1) = tmpsum(ss);
    }
  }  

  List rres; 
  rres["res"]=res; 
  return(rres);
} /*}}}*/

RcppExport SEXP riskstrataR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); tmpsum.zeros(); 
//  colvec res = a; 
  mat res(n,nstrata); res.zeros(); 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(n-i-1); 
       tmpsum(ss) += a(n-i-1); 
       res(n-i-1,ss) = a(n-i-1);
  }  

  List rres; 
  rres["risk"]=res; 
  return(rres);
} /*}}}*/

//
//RcppExport SEXP riskidstrataR(SEXP ia,SEXP iid, SEXP inid, SEXP istrata, SEXP instrata) {/*{{{*/
//  colvec a = Rcpp::as<colvec>(ia);
//  IntegerVector intstrata(istrata); 
//  int nstrata = Rcpp::as<int>(instrata);
//  unsigned n = a.n_rows;
//  int nid = Rcpp::as<int>(inid);
//  IntegerVector id(iid); 
//  int ss,lid; 
//
//  mat tmpsuma(nstrata,nid); tmpsuma.zeros(); 
////  colvec res = a; 
//  mat res(n,nid); res.zeros(); 
//  for (unsigned i=0; i<n; i++) {
//       ss=intstrata(n-i-1); lid=id(n-i-1); 
//       tmpsuma(ss,lid) += a(n-i-1); 
//       res(n-i-1,lid) = a(n-i-1); 
//  }
//
//  List rres; 
//  rres["risk"]=res; 
//  return(rres);
//} /*}}}*/
//


colvec revcumsumstrata(const colvec &a,IntegerVector strata,int nstrata) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec res = a; 

  for (unsigned i=0; i<n; i++) {
    int ss=strata(n-i-1); 
    if ((ss<nstrata) & (ss>=0))  {
       tmpsum(ss) += a(n-i-1); 
       res(n-i-1) = tmpsum(ss);
    }
  }  
  return(res);
}// }}}

colvec revcumsumstrata1(const colvec &a,const  colvec &v1,const  colvec &v2, IntegerVector strata,int nstrata) {
  return(revcumsumstrata(a%v1,strata,nstrata)/v2);
}

colvec revcumsumstratalag(const colvec &a,IntegerVector strata,int nstrata) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec res = a; 

  for (unsigned i=0; i<n; i++) {
       int ss=strata(n-i-1); 
       res(n-i-1) = tmpsum(ss);
       tmpsum(ss) += a(n-i-1); 
  }  
  return(res);
}/*}}}*/


mat  revcumsumstrataMatCols(const mat  &a,const  colvec &v1,const  colvec &v2, IntegerVector strata,int nstrata) { // {{{
  mat res =a; 
  unsigned p=a.n_cols; 
  for (unsigned j=0; j<p; j++) {
    res.col(j) = revcumsumstrata1(a.col(j),v1,v2,strata,nstrata);
  }
  return(res); 
}/*}}}*/


RcppExport SEXP cumsumstratasumR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec ressum = a; 
  colvec lagressum = a; 
  colvec ressqu = a; 
  double cumsum=0; 
  int first=0,ss; 
  for (unsigned i=0; i<n; i++) {
    ss=intstrata(i); 
    // valid strata update 
    if ((first>0.1) & (i>=1)& (ss<nstrata) & (ss>=0))
       ressqu(i)=ressqu(i-1)+pow(a(i),2)+2*a(i)*tmpsum(ss);
    if ((ss<nstrata) & (ss>=0))  {
       lagressum(i)=tmpsum(ss); 
       tmpsum(ss) += a(i); 
       cumsum+=a(i); 
       if (first<0.1) ressqu(i) = pow(a(i),2); 
       first=1; 
    }
    ressum(i) = cumsum;
  }  

  List rres; 
  rres["sumsquare"]=ressqu; 
  rres["sum"]=ressum; 
  rres["lagsum"]=lagressum; 
  return(rres);
} /*}}}*/

RcppExport SEXP revcumsumstratasumR(SEXP ia,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec tmpsqr(nstrata);  tmpsqr.zeros(); 
  colvec cumsum(nstrata);  cumsum.zeros(); 
  colvec ressum = a; 
  colvec lagressum = a; 
  colvec ressqu = a; 
  colvec lagressqu(n); 
  int ss; 
  for (unsigned i=0; i<n; i++) {
       ss=intstrata(n-i-1); 
       lagressqu(n-i-1)=tmpsqr(ss); 
       lagressum(n-i-1)=cumsum(ss); 
       ressqu(n-i-1)=tmpsqr(ss)+pow(a(n-i-1),2)+2*a(n-i-1)*tmpsum(ss);
       tmpsum(ss) += a(n-i-1); 
       cumsum(ss)+=a(n-i-1); 
       ressum(n-i-1) = cumsum(ss);
       tmpsqr(ss)=ressqu(n-i-1); 
  }  

  List rres; 
  rres["sumsquare"]=ressqu; 
  rres["lagsumsquare"]=lagressqu; 
  rres["sum"]=ressum; 
  rres["lagsum"]=lagressum; 
  return(rres);
} /*}}}*/

RcppExport SEXP covrfR(SEXP ia,SEXP ib, SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  colvec b = Rcpp::as<colvec>(ib);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  colvec tmpsumrev(nstrata); tmpsumrev.zeros(); 
  colvec ressqu = a; 
  int ss; 
  for (unsigned i=0; i<n; i++) {
    ss=intstrata(n-i-1); 
    if ((ss<nstrata) & (ss>=0)) tmpsumrev(ss) += b(n-i-1); 
  }  

  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec tmpsqr(nstrata);  tmpsqr.zeros(); 
//  colvec first(nstrata);  first.zeros(); 
  for (unsigned i=0; i<n; i++) {
    ss=intstrata(i); 
    // valid strata update 
    if (((ss<nstrata) & (ss>=0))) {
       ressqu(i)=tmpsqr(ss)-a(i)*tmpsumrev(ss)+b(i)*tmpsum(ss)+a(i)*b(i);
       tmpsumrev(ss) -= b(i); 
       tmpsum(ss)    += a(i); 
       tmpsqr(ss)=ressqu(i); 
    }
    }

  List rres; rres["covs"]=ressqu; 
  return(rres);
} /*}}}*/

RcppExport SEXP cumsumidstratasumCovR(SEXP ia,SEXP ib,SEXP iid,SEXP inid,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  colvec b = Rcpp::as<colvec>(ib);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector id(iid); 
  int nid = Rcpp::as<int>(inid);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  int lid,ss; 

  mat tmpsuma(nstrata,nid); tmpsuma.zeros(); 
  mat tmpsumb(nstrata,nid); tmpsumb.zeros(); 
  colvec tmpsqr(nstrata);  tmpsqr.zeros(); 
//  colvec ressqu = a; colvec ressumu = a; 
  colvec ressuma = a; 
  colvec ressumb = b; 
  colvec ressqu = a; 
  colvec cumsuma(nstrata); cumsuma.zeros(); 
  colvec cumsumb(nstrata); cumsumb.zeros(); 
  colvec first(nstrata);  first.zeros(); 

  for (unsigned i=0; i<n; i++) {
    ss=intstrata(i); lid=id(i); 
    // valid strata update 
    if ((ss<nstrata) & (ss>=0)) {
       ressqu(i)=tmpsqr(ss)+a(i)*b(i)+a(i)*tmpsumb(ss,lid)+b(i)*tmpsuma(ss,lid);
       tmpsuma(ss,lid) += a(i); 
       tmpsumb(ss,lid) += b(i); 
       cumsuma(ss) += a(i); 
       cumsumb(ss) += b(i); 
       ressuma(i) = cumsum(ss);
       ressumb(i) = cumsum(ss);
       tmpsqr(ss)=ressqu(i); 
    }
  }  

  List rres; 
  rres["sumsquare"]=ressqu; 
  rres["suma"]=ressuma; 
  rres["sumb"]=ressumb; 
  return(rres);
} /*}}}*/

RcppExport SEXP cumsumidstratasumR(SEXP ia,SEXP iid,SEXP inid,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector id(iid); 
  int nid = Rcpp::as<int>(inid);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;
  int lid,ss; 

  mat tmpsum(nstrata,nid); tmpsum.zeros(); 
  colvec tmpsqr(nstrata);  tmpsqr.zeros(); 
//  colvec ressqu = a; colvec ressumu = a; 
  colvec ressum = a; 
  colvec ressumid = a; 
  colvec lagressum = a; 
  colvec ressqu = a; 
  colvec cumsum(nstrata); cumsum.zeros(); 
  colvec first(nstrata);  first.zeros(); 

  for (unsigned i=0; i<n; i++) {
    ss=intstrata(i); lid=id(i); 
    // valid strata update 
	ressqu(i)=tmpsqr(ss)+pow(a(i),2)+2*a(i)*tmpsum(ss,lid);
	tmpsum(ss,lid) += a(i); 
	lagressum(i)=cumsum(ss); 
	cumsum(ss) += a(i); 
        ressum(i) = cumsum(ss);
        ressumid(i) = tmpsum(ss,lid);
        tmpsqr(ss)=ressqu(i); 
  }  

  List rres; 
  rres["sumsquare"]=ressqu; 
  rres["sum"]=ressum; 
  rres["sumidstrata"]=ressumid; 
  rres["lagsum"]=lagressum; 
  return(rres);
} /*}}}*/


colvec cumsumAS(const colvec &a,IntegerVector strata,int nstrata) {/*{{{*/
  unsigned n = a.n_rows;
  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec ressum = a; 
  double sums=0; 
  for (unsigned i=0; i<n; i++) {
    int ss=strata(i); 
    ressum(i) = sums+a(i)-tmpsum(ss); 
    tmpsum(ss)=a(i); 
    sums=ressum(i); 
  }  
  return(ressum);
} /*}}}*/

RcppExport SEXP cumsumASR(SEXP ia, SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;

  a.print("a"); 

  colvec tmpsum(nstrata); tmpsum.zeros(); 
  colvec ressum = a; 
  double sums=0; 
  for (unsigned i=0; i<n; i++) {
    int ss=intstrata(i); 
	  Rprintf(" %ld %ld %ld \n",n,ss,i); 
    ressum(i) = sums+a(i)-tmpsum(ss); 
    tmpsum(ss)=a(i); 
    sums=ressum(i); 
  }  

  List rres; 
  rres["sum"]=ressum; 
  return(rres);
} /*}}}*/


RcppExport SEXP meanriskR(SEXP ia, SEXP iid,SEXP inid,SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
//  IntegerVector a = Rcpp::as<colvec>(ia);
//  IntegerVector a(ia); 
//  int n = a1.n_rows;
  int n = a.n_rows;
  IntegerVector id(iid); 
  int nid = Rcpp::as<int>(inid);
  IntegerVector strata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);

  colvec res=a;  
  colvec mean(n); mean.zeros(); 
  colvec tot(n);  tot.zeros();  

  for (int i=0; i<nid; i++) {
   res=revcumsumstrata(a%whichi(id,n,i),strata,nstrata);
   if (i>0) mean=mean+i*res; 
   tot=tot+res; 
  }  
  mean=mean/tot; 

  List rres; 
  rres["meanrisk"]=mean; 
  rres["risk"]=tot; 
  return(rres);
} /*}}}*/

RcppExport SEXP revcumsumidstratasumR(SEXP ia,SEXP iid, SEXP inid, SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;
  IntegerVector id(iid); 
  int nid = Rcpp::as<int>(inid);
  int lid,ss; 

  mat tmpsum(nstrata,nid); tmpsum.zeros(); 
  colvec tmpsqr(nstrata);  tmpsqr.zeros(); 
//  colvec ressqu = a; colvec ressumu = a; 
  colvec ressum = a; 
  colvec ressumid = a; 
  colvec lagressum(n); 
  colvec ressqu = a; 
  colvec lagressqu(n); // lagressqu.zeros
  colvec cumsum(nstrata); cumsum.zeros(); 
  colvec first(nstrata);  first.zeros(); 

  for (unsigned i=0; i<n; i++) {
       ss=intstrata(n-i-1); lid=id(n-i-1); 
       lagressqu(n-i-1)=tmpsqr(ss);  // previous from  
       ressqu(n-i-1)=tmpsqr(ss)+pow(a(n-i-1),2)+2*a(n-i-1)*tmpsum(ss,lid);
       tmpsum(ss,lid) += a(n-i-1); 
       lagressum(n-i-1)=cumsum(ss);  // previous sum from  
       cumsum(ss)+=a(n-i-1); 
       ressum(n-i-1)=cumsum(ss);
       ressumid(n-i-1)=tmpsum(ss,lid); 
       tmpsqr(ss)=ressqu(n-i-1); 
  }  

  List rres; 
  rres["sumsquare"]=ressqu; 
  rres["lagsumsquare"]=lagressqu; 
  rres["lagsum"]=lagressum; 
  rres["sum"]=ressum; 
  rres["sumidstrata"]=ressumid; 
  return(rres);
} /*}}}*/

RcppExport SEXP revcumsumidstratasumCovR(SEXP ia,SEXP ib,SEXP iid, SEXP inid, SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  colvec b = Rcpp::as<colvec>(ib);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;
  IntegerVector id(iid); 
  int nid = Rcpp::as<int>(inid);
  int lid,ss; 

  mat tmpsuma(nstrata,nid); tmpsuma.zeros(); 
  mat tmpsumb(nstrata,nid); tmpsumb.zeros(); 
  colvec cumsuma(nstrata); cumsuma.zeros(); 
  colvec cumsumb(nstrata); cumsumb.zeros(); 
  colvec tmpsqr(nstrata);  tmpsqr.zeros(); 
//  colvec ressqu = a; colvec ressumu = a; 
  colvec ressum = a; 
  colvec lagressum(n); 
  colvec ressqu = a; 
  colvec lagressqu(n); // lagressqu.zeros
  colvec first(nstrata);  first.zeros(); 

  for (unsigned i=0; i<n; i++) {
    ss=intstrata(n-i-1); lid=id(n-i-1); 
    // valid strata update 
     if ((ss<nstrata) & (ss>=0)) {
//    if ((first(ss)>0.1)) {
       lagressqu(n-i-1)=tmpsqr(ss);  // previous from  
       lagressum(n-i-1)=cumsum(ss);  // previous sum from  
       ressqu(n-i-1)=tmpsqr(ss)+a(n-i-1)*b(n-i-1)+a(n-i-1)*tmpsumb(ss,lid)+b(n-i-1)*tmpsuma(ss,lid);
//    }
    tmpsuma(ss,lid) += a(n-i-1); 
    tmpsumb(ss,lid) += b(n-i-1); 
//    cumsuma(ss)+=a(n-i-1); 
//    cumsumb(ss)+=b(n-i-1); 
    // first 
//    if (first(ss)<0.1) {ressqu(n-i-1)=a(n-i-1)*b(n-i-1); first(ss)=1; lagressqu(n-i-1)=0; lagressum(n-i-1)=0; }
//    ressqu(n-i-1) = sum(tmpsum%tmpsum);
//    tmpsum.print("pp"); 
//    ressuma(n-i-1)=cumsum(ss);
//    ressumb(n-i-1)=cumsum(ss);
    tmpsqr(ss)=ressqu(n-i-1); 
    }
  }  

  List rres; 
  rres["sumsquare"]=ressqu; 
  rres["lagsumsquare"]=lagressqu; 
  return(rres);
} /*}}}*/

RcppExport SEXP covrfstrataR(SEXP ia,SEXP ib, SEXP iid,SEXP inid, SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  colvec b = Rcpp::as<colvec>(ib);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;
  IntegerVector id(iid); 
  int nid = Rcpp::as<int>(inid);
  int lid,ss; 

  mat tmpsumrev(nstrata,nid); tmpsumrev.zeros(); 
  mat tmpsum(nstrata,nid); tmpsum.zeros(); 
  colvec tmpcov(nstrata);  tmpcov.zeros(); 
  colvec ressum = a; colvec rescov = a; 
  colvec cumsum(nstrata); cumsum.zeros(); 
  colvec first(nstrata);  first.zeros(); 

  for (unsigned i=0; i<n; i++) {
    ss=intstrata(n-i-1); lid=id(n-i-1); 
    tmpsumrev(ss,lid) += b(n-i-1); 
  }

  for (unsigned i=0; i<n; i++) {
    ss=intstrata(i); lid=id(i); 
    rescov(i)= tmpcov(ss)+a(i)*tmpsumrev(ss,lid)-b(i)*tmpsum(ss,lid)-a(i)*b(i);
    tmpsum(ss,lid)    += a(i); 
    tmpsumrev(ss,lid) -= b(i); 
    tmpcov(ss)=rescov(i); 
  }  

  List rres; rres["covs"]=rescov; 
  return(rres);
} /*}}}*/

RcppExport SEXP covrfstrataCovR(SEXP ia,SEXP ib,SEXP ia2,SEXP ib2,SEXP iid,SEXP inid, SEXP istrata, SEXP instrata) {/*{{{*/
  colvec a = Rcpp::as<colvec>(ia);
  colvec b = Rcpp::as<colvec>(ib);
  colvec a2 = Rcpp::as<colvec>(ia2);
  colvec b2 = Rcpp::as<colvec>(ib2);
//  mat b = Rcpp::as<mat>(ib);
  IntegerVector intstrata(istrata); 
  int nstrata = Rcpp::as<int>(instrata);
  unsigned n = a.n_rows;
  IntegerVector id(iid); 
  int nid = Rcpp::as<int>(inid);
  int lid,ss; 

  mat tmpsumrev(nstrata,nid); tmpsumrev.zeros(); 
  mat tmpsumrev2(nstrata,nid); tmpsumrev2.zeros(); 
  mat tmpsum(nstrata,nid); tmpsum.zeros(); 
  mat tmpsum2(nstrata,nid); tmpsum2.zeros(); 
  colvec tmpsqr(nstrata);  tmpsqr.zeros(); 
//  colvec ressqu = a; colvec ressumu = a; 
  colvec ressum = a; 
  colvec ressqu = a; 
  colvec cumsum(nstrata); cumsum.zeros(); 
  colvec first(nstrata);  first.zeros(); 

  for (unsigned i=0; i<n; i++) {
    ss=intstrata(n-i-1); lid=id(n-i-1); 
    tmpsumrev(ss,lid) += b(n-i-1); 
    tmpsumrev2(ss,lid) += b2(n-i-1); 
  }

  for (unsigned i=0; i<n; i++) {
    ss=intstrata(i); lid=id(i); 
       ressqu(i)=tmpsqr(ss)-a(i)*tmpsumrev2(ss,lid)+b2(i)*tmpsum(ss,lid)+a(i)*b2(i);
       tmpsumrev2(ss,lid) -= b2(i); 
       tmpsum(ss,lid)    += a(i); 
       tmpsqr(ss)=ressqu(i); 
  }  

  List rres; rres["covs"]=ressqu; 
  return(rres);
} /*}}}*/


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
				SEXP offsetsSEXP,
				SEXP ZXSEXP) {
BEGIN_RCPP/*{{{*/
  colvec beta = Rcpp::as<colvec>(betaSEXP);
  mat X = Rcpp::as<mat>(XSEXP);
  mat XX = Rcpp::as<mat>(XXSEXP);
  mat ZX = Rcpp::as<mat>(ZXSEXP);
  arma::uvec Jumps = Rcpp::as<uvec >(JumpsSEXP);
  arma::Col<int> Sign = Rcpp::as<arma::Col<int> >(SignSEXP);
  IntegerVector strata(strataSEXP);
  int nstrata = Rcpp::as<int>(nstrataSEXP);
  // unsigned n = X.n_rows;
  unsigned p = X.n_cols;
  colvec weights = Rcpp::as<colvec>(weightsSEXP);
  colvec offsets = Rcpp::as<colvec>(offsetsSEXP);

  colvec Xb = X*beta+offsets;
  colvec eXb = exp(Xb)%weights;
  if (Sign.n_rows==eXb.n_rows) { // Truncation
    eXb = Sign%eXb;
  }

  colvec S0 = revcumsumstrata(eXb,strata,nstrata);
  mat E=revcumsumstrataMatCols(X,eXb,S0,strata,nstrata); 

//  for (unsigned j=0; j<p; j++) {
//    E.col(j) = revcumsumstrata1(X.col(j),eXb,S0,strata,nstrata);
//  }

  E = E.rows(Jumps);
  mat E2(E.n_rows, E.n_cols*E.n_cols); // Calculate E' E at each time-point
  for (unsigned i=0; i<E.n_rows; i++) {
    rowvec Xi = E.row(i);
    E2.row(i) = vectorise(Xi.t()*Xi,1);
  }

  mat XX2=revcumsumstrataMatCols(XX,eXb,S0,strata,nstrata); 
//  mat XX2 = XX;
//  for (unsigned j=0; j<XX2.n_cols; j++) { // int S2/S0(s)
//    XX2.col(j) = revcumsumstrata1(XX2.col(j),eXb,S0,strata,nstrata);
//  }

  mat ZX2 = ZX;
  if (ZX.n_rows==X.n_rows) {
     ZX2=revcumsumstrataMatCols(ZX,eXb,S0,strata,nstrata); 
  } 

  XX2 = XX2.rows(Jumps);
  colvec weightsJ=weights.elem(Jumps);  
  S0 = S0.elem(Jumps);
  mat grad = (X.rows(Jumps)-E);        // Score
  vec val =  (Xb.elem(Jumps)-log(S0)); // Partial log-likelihood

//  colvec iweightsJ=1/weightsJ; 
  colvec S02 = S0/weightsJ;            // S0 with weights to estimate baseline 
  mat grad2= vecmatrow(weightsJ,grad); // score  with weights
  vec val2 = weightsJ%val;             // Partial log-likelihood with weights

  mat hesst = -(XX2-E2);               // hessian contributions in jump times 
  mat hess  = reshape(sum(hesst),p,p);
  if (ZX.n_rows==X.n_rows) {
     ZX2 = ZX2.rows(Jumps);
  }

  mat hesst2 = vecmatrow(weightsJ,hesst); // hessian over time with weights 
  mat hess2 = reshape(sum(hesst2),p,p);  // hessian with weights 

//  hesst2.print("hessiantime"); hess2.print("hessian");
//  if (hess.has_nan()) {
//	printf("============================ \n"); 
//	S0.print("S0"); exb.print("exb"); grad.print("grad"); e.print("e"); xx2.print("xx"); X.print("X"); 
//	printf("============================ \n"); 
//	}

  return(Rcpp::List::create(Rcpp::Named("jumps")=Jumps,
			    Rcpp::Named("ploglik")=sum(val2),
			    Rcpp::Named("U")=grad2,
			    Rcpp::Named("gradient")=sum(grad2),
			    Rcpp::Named("hessian")=hess2,
			    Rcpp::Named("hessianttime")=hesst2,
			    Rcpp::Named("S2S0")=XX2,
			    Rcpp::Named("E")=E,
			    Rcpp::Named("S0")=S02,
			    Rcpp::Named("ZXeXb")=ZX2,
			    Rcpp::Named("weights")=weightsJ
			    ));
END_RCPP
  }/*}}}*/

RcppExport SEXP FastCoxPLstrataPO(SEXP betaSEXP,
				SEXP XSEXP,
				SEXP XXSEXP,
				SEXP SignSEXP,
				SEXP JumpsSEXP, 
				SEXP strataSEXP, 
				SEXP nstrataSEXP,
				SEXP weightsSEXP,
				SEXP offsetsSEXP,
				SEXP ZXSEXP,
				SEXP propoddsSEXP) {
BEGIN_RCPP/*{{{*/
  colvec beta = Rcpp::as<colvec>(betaSEXP);
  mat X = Rcpp::as<mat>(XSEXP);
  mat XX = Rcpp::as<mat>(XXSEXP);
  mat ZX = Rcpp::as<mat>(ZXSEXP);
  arma::uvec Jumps = Rcpp::as<uvec >(JumpsSEXP);
  arma::Col<int> Sign = Rcpp::as<arma::Col<int> >(SignSEXP);
  IntegerVector strata(strataSEXP);
  int nstrata = Rcpp::as<int>(nstrataSEXP);
  double propodds = Rcpp::as<double>(propoddsSEXP);
  // unsigned n = X.n_rows;
  unsigned p = X.n_cols;
  colvec weights = Rcpp::as<colvec>(weightsSEXP);
  colvec offsets = Rcpp::as<colvec>(offsetsSEXP);

  colvec Xb = X*beta+offsets;
  colvec eXb = exp(Xb)%weights;
  if (Sign.n_rows==eXb.n_rows) { // Truncation
    eXb = Sign%eXb;
  }

  colvec S0 = revcumsumstrata(eXb,strata,nstrata);
  mat E=revcumsumstrataMatCols(X,eXb,S0,strata,nstrata); 

  E = E.rows(Jumps);
  mat E2(E.n_rows, E.n_cols*E.n_cols); // Calculate E' E at each time-point
  for (unsigned i=0; i<E.n_rows; i++) {
    rowvec Xi = E.row(i);
    E2.row(i) = vectorise(Xi.t()*Xi,1);
  }

  mat XX2=revcumsumstrataMatCols(XX,eXb,S0,strata,nstrata); 

  mat ZX2 = ZX;
  if (ZX.n_rows==X.n_rows) {
     ZX2=revcumsumstrataMatCols(ZX,eXb,S0,strata,nstrata); 
  } 

  XX2 = XX2.rows(Jumps);
  colvec weightsJ=weights.elem(Jumps);  
  S0 = S0.elem(Jumps);
  IntegerVector strataJ = seq_len(Jumps.n_rows);  
  for (unsigned i=0; i<Jumps.n_rows; i++) {
	  strataJ(i)=strata(Jumps(i)); 
//  printf("%d %d %d \n",Jumps.n_rows,strataJ(i),Jumps(i)); 
  }

  colvec pow=cumsumstrataPO(S0/weightsJ,strataJ,nstrata,propodds,eXb.elem(Jumps)); 
  mat DLam =DLambeta(weightsJ,S0,E,X.rows(Jumps),strataJ,nstrata,propodds,eXb.elem(Jumps)); 
//  DLam.print(" er den der"); 

  mat grad = (X.rows(Jumps)-E);        // Score
  vec val =  (Xb.elem(Jumps)-log(S0)); // Partial log-likelihood

  colvec S02 = S0/(pow%weightsJ);            // S0 with weights to estimate baseline 
  mat grad2  = vecmatrow(pow%weightsJ,grad); // score  with weights
  vec val2   = pow%weightsJ%val;             // Partial log-likelihood with weights

  // derivative adjustment for PO model 
  mat gradAdjt(Jumps.n_rows,p*p); 
  vec eXBj= eXb.elem(Jumps); 
  mat XI= X.rows(Jumps); 
  for (unsigned i=0; i<Jumps.n_rows; i++) {
    rowvec Xi = grad.row(i);
    gradAdjt.row(i)= vectorise((DLam.row(i)*eXBj(i)+(pow(i)-1)*XI.row(i)).t()*Xi,1); 
  }

  mat hesst = -(XX2-E2);               // hessian contributions in jump times 
  hesst = vecmatrow(pow,hesst); 
  hesst = hesst+ gradAdjt; 
//  mat hess  = reshape(sum(hesst),p,p);
  if (ZX.n_rows==X.n_rows) {
     ZX2 = ZX2.rows(Jumps);
  }

//  mat hesst2 = vecmatrow(pow%weightsJ,hesst); // hessian over time with weights 
  mat hesst2 = vecmatrow(weightsJ,hesst); // hessian over time with weights 
  mat hess2 = reshape(sum(hesst2),p,p);  // hessian with weights 

  return(Rcpp::List::create(Rcpp::Named("jumps")=Jumps,
			    Rcpp::Named("ploglik")=sum(val2),
			    Rcpp::Named("U")=grad2,
			    Rcpp::Named("gradient")=sum(grad2),
			    Rcpp::Named("hessian")=hess2,
			    Rcpp::Named("hessianttime")=hesst2,
			    Rcpp::Named("S2S0")=XX2,
			    Rcpp::Named("E")=E,
			    Rcpp::Named("S0")=S02,
			    Rcpp::Named("ZXeXb")=ZX2,
			    Rcpp::Named("weights")=weightsJ,
			    Rcpp::Named("propoddsW")=pow,
			    Rcpp::Named("DLam")=DLam 
			    ));
END_RCPP
  }/*}}}*/


RcppExport SEXP FastCoxPLstrataAddGam(SEXP betaSEXP,
				SEXP XSEXP,
				SEXP XXSEXP,
				SEXP SignSEXP,
				SEXP JumpsSEXP, 
				SEXP strataSEXP, 
				SEXP nstrataSEXP,
				SEXP weightsSEXP,
				SEXP offsetsSEXP,
				SEXP ZXSEXP,
		SEXP itheta, SEXP idimthetades,SEXP ithetades, 
		SEXP iags, SEXP ivarlink, SEXP idimjumprv,SEXP ijumprv,
		SEXP iJumpsCauses
		 ) {
BEGIN_RCPP/*{{{*/
  colvec beta = Rcpp::as<colvec>(betaSEXP);
  mat X = Rcpp::as<mat>(XSEXP);
  mat XX = Rcpp::as<mat>(XXSEXP);
  mat ZX = Rcpp::as<mat>(ZXSEXP);
  arma::uvec Jumps = Rcpp::as<uvec >(JumpsSEXP);
  arma::Col<int> Sign = Rcpp::as<arma::Col<int> >(SignSEXP);
  IntegerVector strata(strataSEXP);
  int nstrata = Rcpp::as<int>(nstrataSEXP);
//  double propodds = Rcpp::as<int>(propoddsSEXP);
  // unsigned n = X.n_rows;
  unsigned p = X.n_cols;
  colvec weights = Rcpp::as<colvec>(weightsSEXP);
  colvec offsets = Rcpp::as<colvec>(offsetsSEXP);


// {{{ reading in matrices and cubes for AddGam  cause is in strata 
    vec                  theta = Rcpp::as<vec>(itheta); 
    mat                    ags = Rcpp::as<mat>(iags);
    int                varlink = Rcpp::as<int>(ivarlink);

// array for xjump covariates of jump subject, for all causes 
// NumericVector vxjump(ixjump);
// IntegerVector arrayDims(idimxjump);
// arma::cube xjump(vxjump.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vecthetades(ithetades);
 IntegerVector arrayDims1(idimthetades);
 arma::cube thetades(vecthetades.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vrv(ijumprv);
 IntegerVector arrayDims2(idimjumprv);
 arma::cube rv(vrv.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);

 // indeces of the causes relatd to the two jumps
 arma::umat JumpsCauses = Rcpp::as<umat >(iJumpsCauses);
// arma::uvec ij1 = JumpsCauses.col(0); 
// arma::uvec ij2 = JumpsCauses.col(1); 
 // }}}
 
  vec etheta=theta; 
  if (varlink==1) etheta=exp(theta); 

  colvec Xb = X*beta+offsets;
  colvec eXb = exp(Xb)%weights;
  if (Sign.n_rows==eXb.n_rows) { // Truncation
    eXb = Sign%eXb;
  }

  colvec S0 = revcumsumstrata(eXb,strata,nstrata);
  mat E=revcumsumstrataMatCols(X,eXb,S0,strata,nstrata); 

//  for (unsigned j=0; j<p; j++) {
//    E.col(j) = revcumsumstrata1(X.col(j),eXb,S0,strata,nstrata);
//  }

  E = E.rows(Jumps);
  mat E2(E.n_rows, E.n_cols*E.n_cols); // Calculate E' E at each time-point
  for (unsigned i=0; i<E.n_rows; i++) {
    rowvec Xi = E.row(i);
    E2.row(i) = vectorise(Xi.t()*Xi,1);
  }

  mat XX2=revcumsumstrataMatCols(XX,eXb,S0,strata,nstrata); 
//  mat XX2 = XX;
//  for (unsigned j=0; j<XX2.n_cols; j++) { // int S2/S0(s)
//    XX2.col(j) = revcumsumstrata1(XX2.col(j),eXb,S0,strata,nstrata);
//  }

  mat ZX2 = ZX;
  if (ZX.n_rows==X.n_rows) {
     ZX2=revcumsumstrataMatCols(ZX,eXb,S0,strata,nstrata); 
  } 

  XX2 = XX2.rows(Jumps);
  colvec weightsJ=weights.elem(Jumps);  
  S0 = S0.elem(Jumps);

  IntegerVector strataJ = seq_len(Jumps.n_rows);  
  for (unsigned i=0; i<Jumps.n_rows; i++) {
	  strataJ(i)=strata(Jumps(i)); 
  }
//  colvec pow=cumsumstrataPO(S0,strataJ,nstrata,propodds,eXb.elem(Jumps)); 

  // for now use that covariates are the same for the two causes 
colvec  pow=cumsumstrataAddGam(S0,strataJ,nstrata,eXb,etheta,thetades,rv,ags,Jumps); 


  mat grad = (X.rows(Jumps)-E);        // Score
  vec val =  (Xb.elem(Jumps)-log(S0)); // Partial log-likelihood

  colvec S02 = S0/(pow%weightsJ);            // S0 with weights to estimate baseline 
  mat grad2  = vecmatrow(pow%weightsJ,grad); // score  with weights
  vec val2   = pow%weightsJ%val;             // Partial log-likelihood with weights

  mat hesst = -(XX2-E2);               // hessian contributions in jump times 
  mat hess  = reshape(sum(hesst),p,p);
  if (ZX.n_rows==X.n_rows) {
     ZX2 = ZX2.rows(Jumps);
  }

//  mat hesst2 = vecmatrow(pow%weightsJ,hesst); // hessian over time with weights 
  mat hesst2 = vecmatrow(weightsJ,hesst); // hessian over time with weights 
  mat hess2 = reshape(sum(hesst2),p,p);  // hessian with weights 

//  if (hess.has_nan()) {
//	printf("============================ \n"); 
//	S0.print("S0"); exb.print("exb"); grad.print("grad"); e.print("e"); xx2.print("xx"); X.print("X"); 
//	printf("============================ \n"); 
//	}

  return(Rcpp::List::create(Rcpp::Named("jumps")=Jumps,
			    Rcpp::Named("ploglik")=sum(val2),
			    Rcpp::Named("U")=grad2,
			    Rcpp::Named("gradient")=sum(grad2),
			    Rcpp::Named("hessian")=hess2,
			    Rcpp::Named("hessianttime")=hesst2,
			    Rcpp::Named("S2S0")=XX2,
			    Rcpp::Named("E")=E,
			    Rcpp::Named("S0")=S02,
			    Rcpp::Named("ZXeXb")=ZX2,
			    Rcpp::Named("weights")=weightsJ
			    ));
END_RCPP
  }/*}}}*/


RcppExport SEXP Matdoubleindex(SEXP im,SEXP irow,SEXP icols,SEXP ilength) {/*{{{*/
  mat m = Rcpp::as<mat>(im);
  IntegerVector rows(irow); 
  IntegerVector cols(icols); 
  unsigned l = Rcpp::as<int>(ilength);
  vec res(l); 
  for (unsigned i=0; i<l; i++) res(i)=m(rows(i),cols(i)); 
  List rres; 
  rres["mat"]=res; 
  return(rres);
} /*}}}*/

mat CubeVecC(mat XX, vec beta,int dim1) {/*{{{*/
  unsigned p = beta.n_rows;
  unsigned n = XX.n_rows;

  mat XXbeta(n,dim1);
  for (unsigned j=0; j<n; j++)  {
	  XXbeta.row(j)=(reshape(XX.row(j),dim1,p)*beta).t();
  }
  return(XXbeta); 
}/*}}}*/

RcppExport SEXP CubeVec(SEXP XXSEXP, SEXP betaSEXP)
		  {
BEGIN_RCPP/*{{{*/
  colvec beta = Rcpp::as<colvec>(betaSEXP);
  mat XX = Rcpp::as<mat>(XXSEXP);
  unsigned p = beta.n_rows;
  unsigned n = XX.n_rows;

  mat XXbeta(n,p);
  for (unsigned j=0; j<n; j++)  {
	  XXbeta.row(j)=(reshape(XX.row(j),p,p)*beta).t();
  }

  return(Rcpp::List::create(Rcpp::Named("XXbeta")=XXbeta));
END_RCPP
}/*}}}*/

RcppExport SEXP CubeMat(SEXP XXSEXP,SEXP XSEXP)
		  {
BEGIN_RCPP/*{{{*/
  mat XX = Rcpp::as<mat>(XXSEXP);
  mat X  = Rcpp::as<mat>(XSEXP);
  unsigned p = X.n_cols;
  unsigned n = XX.n_rows;

  mat XXX(n,p*p);
  for (unsigned j=0; j<n; j++)  {
	  XXX.row(j)=(vectorise(reshape(XX.row(j),p,p)*X)).t();
  }

  return(Rcpp::List::create(Rcpp::Named("XXX")=XXX));
END_RCPP
}/*}}}*/

mat  vecmatmat(mat a,mat b) 
{/*{{{*/
  unsigned n = b.n_rows;
  unsigned p1 = a.n_cols;
  unsigned p2 = b.n_cols;

  mat res(n,p1*p2); 
  for (unsigned i=0; i<n; i++) {
     rowvec ai = a.row(i);
     rowvec bi = b.row(i);
     res.row(i)=vectorise(bi.t()*ai,1);
//     mat tt=reshape(res.row(i),p1,p2);  // now tt is ai.t() * bi
//     tt.print("tt"); 
  }  
  return(res);
} /*}}}*/

RcppExport SEXP  vecMatMat(SEXP iX,SEXP iZ) {
BEGIN_RCPP/*{{{*/
  arma::mat X = Rcpp::as<arma::mat>(iX);
  arma::mat Z = Rcpp::as<arma::mat>(iZ);

//  unsigned n =Z.n_rows; 
//  unsigned p1=X.n_cols; 
//  unsigned p2=Z.n_cols; 

  mat res=vecmatmat(X,Z); 
 return(Rcpp::List::create(Rcpp::Named("vXZ")=res)); 
END_RCPP
} /*}}}*/

RcppExport SEXP PropTestCox(SEXP iU, SEXP idUt, SEXP insim, SEXP iobssup) {
BEGIN_RCPP/*{{{*/
  mat U      = Rcpp::as<mat>(iU);
  mat dUt = Rcpp::as<mat>(idUt);
  arma::vec osup = Rcpp::as<arma::vec>(iobssup);
  unsigned nsim = Rcpp::as<int>(insim);
  unsigned p = U.n_cols;
  unsigned n = U.n_rows;

  vec pval(p); pval.zeros(); 
  mat Uti(n,p); 
  mat sup(nsim,p); 
  mat simUti(n,50*p); 

  GetRNGstate();  /* to use R random normals */

  for (unsigned j=0; j<nsim; j++) {
     vec nr=rnorm(n); 
     Uti=vecmatrow(nr,U); 
     for (unsigned k=0; k<p; k++)  Uti.col(k) = cumsum(Uti.col(k));
     mat Uthati= CubeVecC(dUt,(Uti.row(n-1)).t(),p); 
     Uthati=Uti-Uthati; 

     for (unsigned k=0;k<p;k++)  {
        sup(j,k)=max(abs(Uthati.col(k))); 
        if ((sup(j,k)>=osup(k))) { pval(k)++; }
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

RcppExport SEXP PropTestCoxClust(SEXP iU, SEXP idUt, SEXP iwrr, SEXP iZ,
		SEXP iLam, SEXP iEdLam, SEXP insim, SEXP iobssup,
		SEXP inclust,SEXP iid, SEXP istrata,SEXP instrata,
		SEXP iJumps) {
BEGIN_RCPP/*{{{*/
  mat U      = Rcpp::as<mat>(iU);
  mat dUt = Rcpp::as<mat>(idUt);
  arma::vec wrr = Rcpp::as<arma::vec>(iwrr);
  arma::vec Lam = Rcpp::as<arma::vec>(iLam);
  mat EdLam = Rcpp::as<mat>(iEdLam);
  mat Z  = Rcpp::as<mat>(iZ);
  IntegerVector id(iid); 
  arma::vec osup = Rcpp::as<arma::vec>(iobssup);
  unsigned nsim = Rcpp::as<int>(insim);
  unsigned nclust = Rcpp::as<int>(inclust);
  unsigned p = U.n_cols;
  unsigned n = U.n_rows;
  IntegerVector strata(istrata); 
  arma::uvec Jumps = Rcpp::as<uvec >(iJumps);
  int nstrata = Rcpp::as<int>(instrata);

  vec pval(p); pval.zeros(); 
  mat Uti(n,p); 
  mat E1(n,p); 
  mat E2(n,p); 
  mat sup(nsim,p); 
  int nj=dUt.n_rows; 
  mat simUti(nj,50*p); 
  vec nr(n); 

  GetRNGstate();  /* to use R random normals */

//     vec vvv=revcumsumstratalag(wrr,strata,nstrata); 
//     vvv.print("vvv"); 

  for (unsigned j=0; j<nsim; j++) {
     vec nnr=rnorm(nclust); 
     for (unsigned k=0; k<n; k++) nr(k) = nnr(id(k)); 
     Uti=vecmatrow(nr,U); 
     for (unsigned k=0; k<p; k++)  Uti.col(k) = cumsum(Uti.col(k));
//     Uti.print("Uti"); 

     vec nrwrr=wrr%nr;
     vec vvv=revcumsumstratalag(nrwrr,strata,nstrata); 
//     vvv.print("vvv"); 

     for (unsigned k=0; k<p; k++) {
         E1.col(k)=revcumsumstratalag(Z.col(k)%nrwrr,strata,nstrata)%Lam;
	 E1.col(k)=cumsumAS(E1.col(k),strata,nstrata); 
     }
     E2=vecmatrow(vvv,EdLam); 
     for (unsigned k=0; k<p; k++) {
	 E2.col(k)=cumsumAS(E2.col(k),strata,nstrata); 
     }

     mat Uthati= CubeVecC(dUt,(Uti.row(n-1)).t(),p); 
//     E1.print("E1"); E2.print("E2"); Uthati.print("Pt"); 
     Uti=Uti-E1+E2; 
     mat UU =Uti.rows(Jumps); 
//     UU.print("uu"); 
     Uthati=UU-Uthati; 

     for (unsigned k=0;k<p;k++)  {
        sup(j,k)=max(abs(Uthati.col(k))); 
        if ((sup(j,k)>=osup(k))) { pval(k)++; }
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

RcppExport SEXP ModelMatrixTestCox(SEXP iU, SEXP idUt,SEXP ibetaiid, SEXP insim, SEXP iobssup) {
BEGIN_RCPP/*{{{*/
  mat U   = Rcpp::as<mat>(iU);
  mat dUt = Rcpp::as<mat>(idUt);
  mat betaiid = Rcpp::as<mat>(ibetaiid);
  arma::vec osup = Rcpp::as<arma::vec>(iobssup);
  unsigned nsim = Rcpp::as<int>(insim);
  unsigned mp = U.n_cols;
  unsigned p = betaiid.n_cols;
  unsigned n = U.n_rows;

  vec pval(mp); pval.zeros(); 
//  mat Uthati(n,mp); 
  mat Uti(n,mp); 
  mat betati(n,p); 
  mat sup(nsim,mp); 
  mat last(nsim,mp); 
  mat simUti(n,50*mp); 

  GetRNGstate();  /* to use R random normals */

  colvec nr(Uti.n_rows);

  for (unsigned j=0;j<nsim; j++) {
     nr=rnorm(n); 
     Uti=vecmatrow(nr,U); 
     betati=vecmatrow(nr,betaiid); 
     for (unsigned k=0;k<mp;k++)  Uti.col(k)=cumsum(Uti.col(k));
     colvec betatti=(sum(betati)).t();
     mat Uthati= CubeVecC(dUt,betatti,mp); 
     Uthati=Uti-Uthati; //     if(j==0) Uthati.print("one sim"); 

     for (unsigned k=0;k<mp;k++)  {
        last(j,k)=Uthati(n-1,k); 
        sup(j,k)=max(abs(Uthati.col(k))); 
        if ((sup(j,k)>=osup(k))) {pval(k)++;}
        if (j<50) { simUti.col(j*mp+k)=Uthati.col(k); }
     }
  }
  pval=pval/nsim; 

  PutRNGstate();  /* to use R random normals */

  return(Rcpp::List::create(Rcpp::Named("supUsim")=sup,
			    Rcpp::Named("last")=last,
			    Rcpp::Named("simUt")=simUti,
			    Rcpp::Named("pval")=pval)); 
END_RCPP
  }/*}}}*/

//RcppExport SEXP simBandCumHazCox(SEXP iU, SEXP idUt,SEXP ibetaiid, SEXP insim,  SEXP isecum) {
//BEGIN_RCPP/*{{{*/
//  mat U   = Rcpp::as<mat>(iU);
//  mat dUt = Rcpp::as<mat>(idUt);
//  mat betaiid = Rcpp::as<mat>(ibetaiid);
//  arma::vec secum = Rcpp::as<arma::vec>(isecum);
//  unsigned nsim = Rcpp::as<int>(insim);
//  unsigned mp = U.n_cols;
//  unsigned p = betaiid.n_cols;
//  unsigned n = U.n_rows;
//
////  vec pval(mp); pval.zeros(); mat Uthati(n,mp); mat simUti(n,50*mp); 
//  mat Uti(n,mp); 
//  mat betati(n,p); 
//  mat sup(nsim,mp); 
//  mat simUti(n,50*mp); 
//
//  GetRNGstate();  /* to use R random normals */
//  colvec nr(Uti.n_rows);
//
//  for (unsigned j=0;j<nsim; j++) {
//     nr=rnorm(n); 
//     Uti=vecmatrow(nr,U); 
//     betati=vecmatrow(nr,betaiid); 
//     colvec betatti=(sum(betati)).t();
//     mat Uthati= dUt * betatti; 
//     for (unsigned k=0;k<mp;k++)  Uti.col(k)=cumsum(Uti.col(k));
//     Uthati=Uti-Uthati; //     if(j==0) Uthati.print("one sim"); 
//     for (unsigned k=0;k<mp;k++)  {
//        sup(j,k)=max(abs(Uthati.col(k)/secum)); 
//        if (j<50) { simUti.col(j*mp+k)=Uthati.col(k); }
//     }
//  }
//
//  PutRNGstate();  /* to use R random normals */
//
//  return(Rcpp::List::create(Rcpp::Named("supUsim")=sup,
//			    Rcpp::Named("simUt")=simUti)); 
//END_RCPP
//  }/*}}}*/
//
