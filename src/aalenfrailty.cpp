#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
#include "twostage.h"

using namespace arma;

RcppExport SEXP Bhat(SEXP ds,
		     SEXP Xs, 
		     SEXP theta, 		     
		     SEXP id, 
		     SEXP ididx,
		     SEXP idsize) { // {{{ 
  try {

    uvec           event = Rcpp::as<uvec>(ds); 
    mat                X = Rcpp::as<mat>(Xs);
    mat               Xc = zeros<mat>(X.n_cols,X.n_cols);
    double      thetahat = Rcpp::as<double>(theta);
    unsigned  stop,start = X.n_rows;
    uvec        eventpos = find(event==1);
    mat               dB = zeros(eventpos.n_elem,X.n_cols);
    uvec         cluster = Rcpp::as<uvec>(id);

    uvec clustersize, clustpos;
    umat clusterindex;
    bool doclust = (Rf_isNull)(idsize);
    if (!doclust)  {
      clustersize        = Rcpp::as<uvec>(idsize);
      clusterindex       = Rcpp::as<umat>(ididx);
    }

    // Obtain usual estimates of increments, dB, of the
    // cumulative time-varying effects in Aalens Additive Model
    for (unsigned ij=0; ij<eventpos.n_elem; ij++) {
      unsigned rij = eventpos.n_rows-ij-1;
      stop  = start-1;
      start = eventpos(rij);
      mat Xij = X(span(start,stop), span::all);
      Xc = Xc + Xij.st()*Xij;
      mat U, V; vec s; 
      svd(U,s,V,Xc);
      mat Xci = U*diagmat(1/s)*V.st();
      dB.row(rij) = trans(Xci*trans(X.row(start)));
    }
    if (thetahat==0) {
      return(Rcpp::List::create(Rcpp::Named("dB")=dB)); //Increments of marg. aalenn
    }


    vec       Hij(X.n_rows); // Vector to hold cumulative hazard; updated as t increases
    Hij.fill(0);
    mat        B2 = zeros(eventpos.n_elem,X.n_cols); // Cond. cumulative coef 
    for (unsigned k=0; k<eventpos.n_elem; k++) { // Iterate over events
      unsigned ij = eventpos(k);
      unsigned i = cluster(ij);  // cluster
      if (doclust) {
         clustpos = find(cluster==i); // position of within cluster observations
      } else {
        unsigned csize = clustersize(i);
        clustpos  = conv_to<uvec>::from(clusterindex.submat(i,0,i,csize-1));
      }
      uvec posL = find(clustpos<ij); // earlier events/censoring within cluster
      uvec posR = find(clustpos>=ij); // later/current events within cluster
      unsigned Ni = 0; // Number of events in cluster before current event,time t-
      double Hi = 0 ; // Sum of cum.haz. within cluster up to time t-
      if (posL.n_elem>0) {
	Ni = sum(event.elem(clustpos.elem(posL)));
	Hi = sum(Hij.elem(clustpos.elem(posL)));
      }
      uvec pos;
      if (posR.n_elem>0 && k>0) {
	pos = clustpos.elem(posR);
	mat Xi = X.rows(pos);
	Hij.elem(pos) = Xi*trans(B2.row(k-1));
	Hi += sum(Hij.elem(pos));
      }
      double psi = (1/thetahat+Ni)/(1/thetahat+fmax(0,Hi));
      B2.row(k) = dB.row(k)/psi;
      if (k>0) { B2.row(k) += B2.row(k-1); }
    }
    return(Rcpp::List::create(Rcpp::Named("dB")=dB, //Increments of marg. aalen
			      Rcpp::Named("B2")=B2  // Cum.coef of frailty aalen
			      ));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}

RcppExport SEXP Uhat(SEXP ds, SEXP H, SEXP theta, SEXP id, SEXP idsize) {
  try { // {{{ 
    uvec                event = Rcpp::as<uvec>(ds);
    vec                   Hij = Rcpp::as<vec>(H);
    double           thetahat = Rcpp::as<double>(theta);
    umat              cluster = Rcpp::as<umat>(id);
    uvec clustersize, ucluster, clustpos;
    unsigned nclust;
    bool doclust = (Rf_isNull)(idsize);
    //bool doclust = (cluster.n_cols==1);
    if (doclust)  {
      ucluster              = unique(cluster);
      nclust                = ucluster.n_elem;
    } else {
      clustersize           = Rcpp::as<uvec>(idsize);
      nclust = cluster.n_rows;
    }
    vec res(nclust);
    for (unsigned i=0; i<nclust; i++) {
      if (doclust) {
        unsigned ic = ucluster(i);  // cluster 
        clustpos = find(cluster==ic); // position of within cluster observations
      } else {
        unsigned csize = clustersize(i);
        clustpos  = conv_to<uvec>::from(cluster.submat(i,0,i,csize-1));
	//cluster(span(i,i),span(0,csize-1));
      }
      double Ni = sum(event.elem(clustpos));
      double Hi = sum(Hij.elem(clustpos));
      double thetaH = thetahat*Hi+1;
      double R = (log(thetaH)/thetahat + (Ni-Hi)/(thetaH));
      for (unsigned h=0; h<Ni; h++) R -= 1/(1+thetahat*h);
      res(i) = R/thetahat;
    }
    return(Rcpp::wrap(res));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}} 


//RcppExport SEXP Bhatprobandpairs( SEXP ds, SEXP Xs, SEXP theta, 		     SEXP id, 
//		     SEXP ididx, SEXP idsize, 
//		     SEXP iweights
////		     SEXP istatusproband, SEXP itimeproband, SEXP rvdes, SEXP thetades, SEXP antrvs
//		     ) 
//{ // {{{ 
//  try {
//
//    uvec           event = Rcpp::as<uvec>(ds); 
//    mat                X = Rcpp::as<mat>(Xs);
//    mat               Xc = zeros<mat>(X.n_cols,X.n_cols);
//    double      thetahat = Rcpp::as<double>(theta);
//    unsigned  stop,start = X.n_rows;
//    uvec        eventpos = find(event==1);
//    mat               dB = zeros(eventpos.n_elem,X.n_cols);
//    uvec         cluster = Rcpp::as<uvec>(id);
//    vec          weights = Rcpp::as<vec>(iweights);
////    vec          timeproband   = Rcpp::as<vec>(itimeproband);
////    uvec         statusproband = Rcpp::as<uvec>(istatusproband);
//
//    uvec clustersize, clustpos;
//    umat clusterindex;
//    bool doclust = (Rf_isNull)(idsize);
//    if (!doclust)  {
//      clustersize        = Rcpp::as<uvec>(idsize);
//      clusterindex       = Rcpp::as<umat>(ididx);
//    }
//
//    // Obtain usual estimates of increments, dB, of the
//    // cumulative time-varying effects in Aalens Additive Model
//    for (unsigned ij=0; ij<eventpos.n_elem; ij++) {
//      unsigned rij = eventpos.n_rows-ij-1;
//      stop  = start-1;
//      start = eventpos(rij);
//      mat Xij = X(span(start,stop), span::all);
//      Xc = Xc + Xij.st()*Xij;
//      mat U, V; vec s; 
//      svd(U,s,V,Xc);
//      mat Xci = U*diagmat(1/s)*V.st();
//      dB.row(rij) = trans(Xci*trans(X.row(start)));
//    }
//    if (thetahat==0) {
//      return(Rcpp::List::create(Rcpp::Named("dB")=dB)); //Increments of marg. aalenn
//    }
//
//
//    vec       Hij(X.n_rows); // Vector to hold cumulative hazard; updated as t increases
//    Hij.fill(0);
//    mat        B2 = zeros(eventpos.n_elem,X.n_cols); // Cond. cumulative coef 
//    for (unsigned k=0; k<eventpos.n_elem; k++) { // Iterate over events
//      unsigned ij = eventpos(k);
//      unsigned i = cluster(ij);  // cluster
//      if (doclust) {
//         clustpos = find(cluster==i); // position of within cluster observations
//      } else {
//        unsigned csize = clustersize(i);
//        clustpos  = conv_to<uvec>::from(clusterindex.submat(i,0,i,csize-1));
//      }
//      uvec posL = find(clustpos<ij); // earlier events/censoring within cluster
//      uvec posR = find(clustpos>=ij); // later/current events within cluster
//      unsigned Ni = 0; // Number of events in cluster before current event,time t-
//      double Hi = 0 ; // Sum of cum.haz. within cluster up to time t-
//
////      if (posL.n_elem>0) {
////	Ni = sum(event.elem(clustpos.elem(posL)));
////	Hi = sum(Hij.elem(clustpos.elem(posL)));
////      }
////      uvec pos;
////      if (posR.n_elem>0 && k>0) {
////	pos = clustpos.elem(posR);
////	mat Xi = X.rows(pos);
////	Hij.elem(pos) = Xi*trans(B2.row(k-1));
////	Hi += sum(Hij.elem(pos));
////      }
//
//      B2.row(k) = dB.row(k)*weights(k);
//      if (k>0) { B2.row(k) += B2.row(k-1); }
//    }
//    return(Rcpp::List::create(Rcpp::Named("dB")=dB, //Increments of marg. aalen
//			      Rcpp::Named("B2")=B2  // Cum.coef of frailty aalen
//			      ));
//  } catch( std::exception &ex ) {
//    forward_exception_to_r( ex );
//  } catch(...) {  
//    ::Rf_error( "c++ exception (unknown reason)" ); 
//  }
//  return R_NilValue; // -Wall
//} // }}}

RcppExport SEXP BhatAddGam(SEXP idBaalen,SEXP icause,
		SEXP idimxjump,SEXP ixjump, // cube 
		SEXP itheta,
		SEXP idimthetades,SEXP ithetades, // cube 
		SEXP iags, SEXP ivarlink, 
		SEXP idimjumprv,SEXP ijumprv)  // cube 
{ // {{{ 
  try {

//   wall_clock timer; 
//   timer.tic(); 

// {{{ reading in matrices and cubes 
    mat                dBaalen = Rcpp::as<mat>(idBaalen);
    uvec                 cause = Rcpp::as<uvec>(icause); 
    vec                  theta = Rcpp::as<vec>(itheta); 
    mat                    ags = Rcpp::as<mat>(iags);
    int                varlink = Rcpp::as<int>(ivarlink);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vxjump(ixjump);
 IntegerVector arrayDims(idimxjump);
 arma::cube xjump(vxjump.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vecthetades(ithetades);
 IntegerVector arrayDims1(idimthetades);
 arma::cube thetades(vecthetades.begin(), arrayDims1[0], arrayDims1[1], arrayDims1[2], false);

// array for xjump covariates of jump subject, for all causes 
 NumericVector vrv(ijumprv);
 IntegerVector arrayDims2(idimjumprv);
 arma::cube rv(vrv.begin(), arrayDims2[0], arrayDims2[1], arrayDims2[2], false);

 // }}}
 
// double nt = timer.toc();
// printf("timer-ind %lf \n",nt); 

  vec casev(cause.n_elem); 
  vec etheta=theta; 
  if (varlink==1) etheta=exp(theta); 

//  xjump for each jump contains matrix of covariates such that vec cumhaz1= xjump.slice(s) * Bhat 

    mat  Bhat(dBaalen.n_rows, dBaalen.n_cols); 
//    cube  DthetaBhat(theta.n_elem, dBaalen.n_cols,dBaalen.n_rows); 
//    vec dBB(theta.n_elem); 

//    Bhat.fill(0); // initialize 
//
    vec DthetaS(theta.n_elem),DthetaDtS(theta.n_elem),DthetaW(theta.n_elem); 
    vec allvec(6); 
    int ncr=rv.n_rows; 
    vec cumhaz1(ncr); cumhaz1.fill(0); 
    vec Dcumhaz1(ncr); 
//    vec cumhaz2(ncr); cumhaz2.fill(0); 

    double  caseweight=1,ll; 
//    mat rv2=0*rv.slice(0); 
    mat rv1=rv.slice(0); 

//	rv1.print("rv1"); 
//	cumhaz1.print("ch1"); 
//	ags.print("ags"); 

 wall_clock timer; 
 timer.tic(); 

    for (unsigned k=0; k<cause.n_elem; k++) { // Iterate over events
	    // {{{ 

        // computes weights based on additive gamma model 
        mat thetadesv=thetades.slice(k); 
	rv1=rv.slice(k); 
//	thetadesv.print("thetades"); 
//	rv1.print("rv1"); 
        ll=survivalRVCmarg(etheta,thetadesv,ags,(int) cause(k),cumhaz1,rv1,DthetaS,DthetaDtS,allvec,Dcumhaz1);
        caseweight=allvec(0)/ll; //   S / D_1 S
	casev(k)=caseweight; 
	DthetaW=(ll*DthetaDtS-allvec(0)*DthetaS)/(ll*ll);

        //  increments 
        Bhat.row(k)=dBaalen.row(k)*caseweight;
//        DthetaBhat.slice(k)= DthetaW * dBaalen.row(k);

	// derivative of baseline wrt theta
//	mat dBthetamat=xjump.slice(k) * trans(DthetaBhat.slice(k)); 
//	dBB = trans(dBthetamat) *Dcumhaz1; 

//  cumulative  for all causes
        if (k>0) { Bhat.row(k) += Bhat.row(k-1); }
//        if (k>0) { DthetaBhat.slice(k)+= DthetaBhat.slice(k-1)+ 
//                      (dBB * dBaalen.row(k)); 
//	}
	mat xj=xjump.slice(k); 
//	vec bb=trans(Bhat.row(k)); 
//	bb.print("bb"); 
//	cumulative hazard at time t- for all causes 
	cumhaz1=xjump.slice(k) * trans(Bhat.row(k)); 
//	mat pp=DthetaBhat.slice(k); 
//	pp.print("pp"); Dcumhaz1.print("Dcumhaz1");  

//		trans(sum(dBthetamat,0)); 
    } // }}}

 double nt2 = timer.toc();
 printf("Bhat-profile timer-loop %lf \n",nt2); 

    return(Rcpp::List::create(Rcpp::Named("B")=Bhat, 
			      Rcpp::Named("caseweights")=casev)
//			      Rcpp::Named("DthetaBhat")=-1*DthetaBhat)
		    );
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}


RcppExport SEXP MatxCube(
		SEXP imat,
		SEXP idim,SEXP iDBhat 
		)  
{ // {{{ 
  try {

 mat  xmat = Rcpp::as<mat>(imat);

 NumericVector vDBhat(iDBhat);
 IntegerVector arrayDims(idim);
 arma::cube DBhat(vDBhat.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

 mat X(arrayDims[2],arrayDims[0]); 

 for (int k=0; k<arrayDims[2]; k++) { // Iterate over events
	    X.row(k)= xmat.row(k) * trans(DBhat.slice(k)); 
 } 

    return(Rcpp::List::create(Rcpp::Named("X")=X));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
} // }}}

