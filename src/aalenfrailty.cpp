#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

using namespace arma;

RcppExport SEXP Bhat(SEXP ds,
		     SEXP Xs, 
		     SEXP theta, 		     
		     SEXP id, 
		     SEXP ididx,
		     SEXP idsize) {
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


    vec       Hij(X.n_rows); // Vector to hold cumulative hazard; updated as t increases
    Hij.fill(0);
    mat        B2 = zeros(eventpos.n_elem,X.n_cols); // Conditional cumulative hazard 
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
    return(Rcpp::List::create(Rcpp::Named("dB")=dB,
			      Rcpp::Named("B2")=B2
			      ));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
}


RcppExport SEXP Uhat(SEXP ds, SEXP H, SEXP theta, SEXP id, SEXP idsize) {
  try {
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
      for (unsigned h=0; h<=fmax(0,Ni-1); h++) R -= 1/(1+thetahat*h);
      res(i) = R/thetahat;
    }
    return(Rcpp::wrap(res));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
}

