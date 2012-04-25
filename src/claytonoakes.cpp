#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace arma;
using namespace Rcpp;

RcppExport SEXP claytonoakes( SEXP ds, SEXP ts, SEXP es, SEXP allcs, 
			 SEXP cs, SEXP cuts, SEXP hs, SEXP mulths, SEXP var ) {
  // try {
  colvec event = Rcpp::as<colvec>(ds);
  colvec time = Rcpp::as<colvec>(ts);
  colvec entry = Rcpp::as<colvec>(es);
  colvec uniqueclusters = Rcpp::as<colvec>(cs);
  colvec clusters = Rcpp::as<colvec>(allcs);
  NumericVector cut = Rcpp::as<NumericVector>(cuts);
  colvec basehaz = Rcpp::as<colvec>(hs);
  colvec mhaz = Rcpp::as<colvec>(mulths);
  colvec theta0 = Rcpp::as<colvec>(var);
  
  unsigned ncluster = uniqueclusters.n_rows;
  unsigned n = time.n_rows;
  unsigned ncuts = cut.size();    
  colvec dt(ncuts-1);
  for (unsigned i=1; i<ncuts; i++) dt(i-1) = cut(i)-cut(i-1);
  colvec L0 = cumsum(basehaz%dt);
  colvec L00(L0.n_rows+1); 
  for (unsigned i=0; i<L0.n_rows; i++) L00(i+1)=L0(i);
  colvec L1(L0.n_rows); L1(0)=1;
  
  double curtheta0;
  int pos,epos;
  NumericVector::iterator it;   
  double logLik=0;
  unsigned j=0;
  for (unsigned i=0; i<ncluster; i++) {
    curtheta0 = theta0(j);
    int curcluster = uniqueclusters(i);
    unsigned nevent = 0;    
    double Haz1=0;
    double Haz0=0;
    double lgammaratio=0;
    do {
      double multhaz = mhaz(j);      
      {
	L1(0)=1;
	for (unsigned k=1; k<L0.n_rows; k++)
	  L1(k)=exp(curtheta0*multhaz*L0(k-1));
	for (unsigned k=0; k<L0.n_rows; k++)
	  L1(k)*=(exp(curtheta0*multhaz*basehaz(k)*dt(k))-1)/(curtheta0);
	L1=cumsum(L1);
      }

      double t = time(j);
      double e = entry(j);
      it = std::lower_bound(cut.begin(), cut.end(), t);
      pos = int(it-cut.begin())-1;
      if (event(j)==1) {
	nevent += event(j);
	// hazard contribution
	logLik += log(multhaz*basehaz(pos)*exp(curtheta0*multhaz*(L00(pos)+basehaz(pos)*(t-cut(pos)))));
	// Ratio of gamma functions
	lgammaratio += log(1/curtheta0+nevent-1);
      }

      Haz1 += 1/(curtheta0)*exp(curtheta0*multhaz*L00(pos))*
	(exp(curtheta0*multhaz*(basehaz(pos)*(t-cut(pos))))-1);
      if (pos>0) Haz1 += L1(pos-1);     

      if (e>0) { // Truncation
	it = std::lower_bound(cut.begin(), cut.end(), e);
	epos= int(it-cut.begin())-1;	
	Haz0 += 1/(curtheta0)*exp(curtheta0*multhaz*L00(epos))*
	  (exp(curtheta0*multhaz*(basehaz(epos)*(e-cut(epos))))-1);      
	if (epos>0) Haz0 += L1(epos-1);     
      }

      j++;
      if (j==n) break;
    } while (clusters(j)==curcluster);
    logLik += -(1/curtheta0+nevent)*log(1+curtheta0*Haz1); // Survival
    logLik += (1/curtheta0)*log(1+curtheta0*Haz0); // Truncation
    logLik += nevent*log(curtheta0);
    logLik += lgammaratio;
  }
  
  return(Rcpp::List::create(Rcpp::Named("logLik")=logLik));  
//  } catch( std::exception &ex ) {
//   forward_exception_to_r( ex );
//  } catch(...) { 
//   ::Rf_error( "C++ exception (unknown reason)" ); 
//  }
// return R_NilValue; // -Wall

}

