#include "mvn.h"

using namespace std;
using namespace Rcpp;
using namespace arma;


// {{{ uniprobit

RcppExport SEXP uniprobit(SEXP m, SEXP dm,
			  SEXP s, SEXP ds, SEXP y, SEXP w, SEXP std,
			  SEXP eqmarg) {
BEGIN_RCPP
  NumericMatrix mm(m); 
  NumericMatrix dss(ds); 
  NumericMatrix dmm(dm);
  NumericMatrix yy(y);
  double S = Rcpp::as<double>(s);
  mat mu(mm.begin(), mm.nrow(), mm.ncol(), false);
  mat dS(dss.begin(), dss.nrow(), dss.ncol(), false);
  mat tdmu(dmm.begin(), dmm.nrow(), dmm.ncol(), false);
  mat Y(yy.begin(), yy.nrow(), yy.ncol(), false);
  bool weights = Rcpp::as<bool>(std);

  unsigned n = mm.nrow();
  NumericVector W;
  if (weights) {
    NumericVector W0(w);
    W=W0;
  }

  double sigma = sqrt(S);  
  colvec alpha(n), alpha0(n), M(n);
  for (unsigned i=0; i<n; i++) {
    alpha[i] = alpha0[i] = Rf_pnorm5(mu[i],0,sigma,1,0);
    if (Y[i]==0) alpha[i] = 1-alpha[i];  
    M[i] = -sigma*Rf_dnorm4(mu[i]/sigma,0.0,1.0,0);
  }
  mat U1(dmm.ncol(),n);
  for (unsigned i=0; i<n; i++) {
    U1.col(i) = -trans(tdmu.row(i)*1/S*M[i]);
  }
  if (dS.n_cols>0) {
    colvec V = S*alpha0+mu%M;  
    mat U2 = 0.5*trans(dS)*trans(-alpha0+V/S)/S;
    U1 = join_cols(U1,U2);    
  }

  mat logLik = log(alpha);
  for (unsigned i=0; i<n; i++) {
    double wi = 1; 
    if (weights) wi = W[i];
    logLik[i] *= wi;
    if (Y[i]==0) wi *= -1;;
    U1.col(i) *= wi/alpha[i];
  }
  
  List res;
  res["score"] = trans(U1);
  res["loglik"] = logLik;  
  return(res); 
END_RCPP
}

// }}} uniprobit


// {{{ logLik

mat logLik(mat &y, mat &mu, mat &Sigma) {
  unsigned n = y.n_rows;
  colvec res(n);
  double l = Sigma(0,0);
  double R = Sigma(0,1)/l;
  for (unsigned i=0;i<n;i++) {
    double mysign = 1;
    rowvec lo = mu.row(i)/sqrt(l);
    if (y(i,0)==1)  lo(0) *= -1;
    if (y(i,1)==1) lo(1) *= -1;
    if (y(i,0)!=y(i,1)) mysign = -1;
    double R0 = R*mysign;
    res(i) = log(Sbvn(lo(0),lo(1),R0)); // calculates lower tail (mult.by -1).
  }
  return(res);
}

// }}} logLik

// {{{ score
vecmat score(mat &y, mat &mu, mat &x, mat &Sigma, mat &dS0, mat &w) {
  unsigned n = y.n_rows;
  mat iSigma = inv(Sigma);  
  double l = Sigma(0,0);
  vec ll = sqrt(diagvec(Sigma)); 
  mat Lambda = diagmat(ll);
  mat iLambda = diagmat(1/ll);
  mat R = iLambda*Sigma*iLambda;
  mat R0 = R;
  double r = Sigma(0,1)/l;
  
  mat LR = Lambda*R;
  mat LR0 = LR;
  mat S= Sigma;
  mat iS = iSigma;
  mat iSxiS1 = kron(iS,iS);
  iS(1,0) = iS(0,1) = -iSigma(1,0);
  mat iSxiS2 = kron(iS,iS);

  vec Lik(n);
  mat Score(n,x.n_cols/2+dS0.n_rows);
  mat U;
  for (unsigned i=0;i<n;i++) {
    int mysign = 1;
    double r0 = r;
    rowvec lo = mu.row(i)%(1/trans(ll));
    mat L = eye(2,2);
    if (y(i,0)==1) { lo(0) *= -1; L(0,0)= -1; } 
    if (y(i,1)==1) { lo(1) *= -1; L(1,1)= -1; }
    mat dS = dS0;    
    mat iSxiS = iSxiS1;    
    if (y(i,0)!=y(i,1)) { 
      dS.col(1) *= -1; dS.col(2) *= -1;
      iSxiS = iSxiS2;
      mysign = -1;
    } 
    S = Sigma;
    S(1,0) = mysign*Sigma(1,0); S(0,1) = S(1,0);
    iS = iSigma;
    iS(1,0) = mysign*iSigma(1,0); iS(0,1) = iS(1,0);
    LR = LR0;
    LR(0,1) = mysign*LR0(1,0); LR(1,0) = LR(0,1);
    r0 *= mysign;
    rowvec up = -lo;
    double alpha = Sbvn(lo(0),lo(1),r0);    
    vecmat D = Dbvn(up(0),up(1),r0);
    mat M = -LR*D.V;    
    mat V = alpha*S + LR*D.M*trans(LR);
    if (!(w.is_empty())) {  
      mat W = reshape(w.row(i),2,2); //diagmat(w.row(i));
      V = W%V;     
      iS = W%iS;
    }
    mat vecV = reshape(V,4,1);
    mat veciS = reshape(iS,4,1);
    mat dmu = L*trans(reshape(x.row(i),x.n_cols/2,2));     
    mat d1 = trans(dmu)*iS*M/alpha;
    mat d2 = -0.5*(dS)*(veciS - iSxiS*vecV/alpha);
    rowvec val = trans(join_cols(d1,d2));
    Score.row(i) = val;
    Lik(i) = alpha;
  }
  vecmat res; 
  res.V = Lik; res.M = Score;
  return(res);
}
// }}} score 

// {{{ biprobit 0

RcppExport SEXP biprobit0(SEXP m, 
			  SEXP s, SEXP ds, SEXP y, SEXP x, SEXP w, SEXP std,
			  SEXP sameweight) {

  NumericMatrix mm(m); 
  NumericMatrix ss(s); 
  NumericMatrix dss(ds); 
  NumericMatrix xx(x);
  NumericMatrix yy(y);

  mat mu(mm.begin(), mm.nrow(), mm.ncol(), false);
  mat S(ss.begin(), ss.nrow(), ss.ncol(), false);
  mat dS(dss.begin(), dss.nrow(), dss.ncol(), false);
  mat X(xx.begin(), xx.nrow(), xx.ncol(), false);
  mat Y(yy.begin(), yy.nrow(), yy.ncol(), false);
  bool weights = Rcpp::as<bool>(std);
  bool OneWeight = Rcpp::as<bool>(sameweight); 

  mat U;
  mat W;
  List res;
  if (weights) {
    NumericMatrix ww(w);
    mat W0(ww.begin(), ww.nrow(), ww.ncol(), false);
    W = W0;
    //    WW = trans(reshape(W1,2,N/2));    
  }
  if (OneWeight & weights) {
    vecmat U0 = score(Y,mu,X,S,dS,U);
    res["loglik"] = log(U0.V)%W.col(0);
    for (unsigned i=0; i<U0.M.n_rows; i++) U0.M.row(i) *= W(i,0);
    res["score"] = U0.M;
  } else {
    vecmat U0 = score(Y,mu,X,S,dS,W);
    res["loglik"] = log(U0.V);
    res["score"] = U0.M;
  }

  return(res);
}

// }}} biprobit0


vecmat scorecor(mat &y, mat &mu, mat &x, mat &r, mat &dr, mat &w, bool eqmarg) {

  unsigned n = y.n_rows;  
  mat I = eye(2,2); 
  mat dS(dr.n_cols,4);
  dS.fill(0);
  
  vec Lik(n);
  unsigned nparmean = x.n_cols;
  if (eqmarg) nparmean *= 0.5;
  mat Score(n,nparmean+dr.n_cols);
  mat U;
  for (unsigned i=0;i<n;i++) {
    int mysign = 1;
    double r0 = r(i);
    mat S = I; 
    rowvec lo = mu.row(i);
    mat L = I;
    if (y(i,0)==1) { lo(0) *= -1; L(0,0)= -1; } 
    if (y(i,1)==1) { lo(1) *= -1; L(1,1)= -1; }
    if (y(i,0)!=y(i,1)) { 
      mysign = -1;
    } 
    dS.col(1) = mysign*trans(dr.row(i)); dS.col(2) = dS.col(1);
    r0 *= mysign;
    S(1,0) = S(0,1) = r0;
    double detS = 1-r0*r0;
    mat iS = S/detS; iS(1,0) = iS(0,1) *= -1;
    mat iSxiS = kron(iS,iS);
    rowvec up = -lo;
    double alpha = Sbvn(lo(0),lo(1),r0);    
    vecmat D = Dbvn(up(0),up(1),r0);
    mat M = -S*D.V;    
    mat V = alpha*S + S*D.M*trans(S);

    if (!(w.is_empty())) {  
      mat W = diagmat(w.row(i));
      V = W*V;     
      iS = W*iS;
    }
    mat vecV = reshape(V,4,1);
    mat veciS = reshape(iS,4,1);
    mat dmu;
    if (eqmarg) {
      dmu = L*trans(reshape(x.row(i),x.n_cols/2,2));
    }
    else {
      dmu = join_cols(x.row(i),x.row(i));
      for (unsigned j=0; j<x.n_cols/2; j++) {
	dmu(1,j) = 0;
	dmu(0,j+x.n_cols/2) = 0;
      }
      dmu = L*dmu;
    }
    // Rcpp::Rcout << "dmu" << dmu << endl;  
    // Rcpp::Rcout << "iS" << iS << endl;  
    //    Rcpp::Rcout << "dS" << dS << endl;  
    // Rcpp::Rcout << "M" << M << endl;  
    // Rcpp::Rcout << "iSxiS" << iSxiS << endl;  
    // Rcpp::Rcout << "vecV" << vecV << endl;      
    mat d1 = trans(dmu)*iS*M/alpha;
    mat d2 = -0.5*(dS)*(veciS - iSxiS*vecV/alpha);

    rowvec val = trans(join_cols(d1,d2));
    Score.row(i) = val;
    Lik(i) = alpha;
  }
  vecmat res; 
  res.V = Lik; res.M = Score;
  return(res);
}












// {{{ score 2
/// allows different marginals for individual 1 and 2
vecmat score2(mat &y, mat &mu, mat &x, mat &Sigma, mat &dS0, mat &w, bool eqmarg) {
  unsigned n = y.n_rows;
  mat iSigma = inv(Sigma);
  double l = Sigma(0,0);
  vec ll = sqrt(diagvec(Sigma)); 
  mat Lambda = diagmat(ll);
  mat iLambda = diagmat(1/ll);
  mat R = iLambda*Sigma*iLambda;
  mat R0 = R;
  double r = Sigma(0,1)/l;
  
  mat LR = Lambda*R;
  mat LR0 = LR;
  mat S= Sigma;
  mat iS = iSigma;
  mat iSxiS1 = kron(iS,iS);
  iS(1,0) = iS(0,1) = -iSigma(1,0);
  mat iSxiS2 = kron(iS,iS);
  
  vec Lik(n);
  unsigned nparmean = x.n_cols;
  if (eqmarg) nparmean *= 0.5;
  mat Score(n,nparmean+dS0.n_rows);
  mat U;
  for (unsigned i=0;i<n;i++) {
    int mysign = 1;
    double r0 = r;
    rowvec lo = mu.row(i)%(1/trans(ll));
     mat L = eye(2,2);
    if (y(i,0)==1) { lo(0) *= -1; L(0,0)= -1; } 
    if (y(i,1)==1) { lo(1) *= -1; L(1,1)= -1; }
    mat dS = dS0;    
    mat iSxiS = iSxiS1;    
    if (y(i,0)!=y(i,1)) { 
      dS.col(1) *= -1; dS.col(2) *= -1;
      iSxiS = iSxiS2;
      mysign = -1;
    } 
    S = Sigma;
    S(1,0) = mysign*Sigma(1,0); S(0,1) = S(1,0);
    iS = iSigma;
    iS(1,0) = mysign*iSigma(1,0); iS(0,1) = iS(1,0);
    LR = LR0;
    LR(0,1) = mysign*LR0(1,0); LR(1,0) = LR(0,1);
    r0 *= mysign;
    rowvec up = -lo;
    double alpha = Sbvn(lo(0),lo(1),r0);    
    vecmat D = Dbvn(up(0),up(1),r0);
    mat M = -LR*D.V;    
    mat V = alpha*S + LR*D.M*trans(LR);
    if (!(w.is_empty())) {  
      mat W = diagmat(w.row(i));
      V = W*V;     
      iS = W*iS;
    }
    mat vecV = reshape(V,4,1);
    mat veciS = reshape(iS,4,1);
    mat dmu;
    if (eqmarg) {
      dmu = L*trans(reshape(x.row(i),x.n_cols/2,2));
    }
    else {
      dmu = join_cols(x.row(i),x.row(i));
      for (unsigned j=0; j<x.n_cols/2; j++) {
	dmu(1,j) = 0;
	dmu(0,j+x.n_cols/2) = 0;
      }
      dmu = L*dmu;
    }
    mat d1 = trans(dmu)*iS*M/alpha;
    mat d2 = -0.5*(dS)*(veciS - iSxiS*vecV/alpha);
    rowvec val = trans(join_cols(d1,d2));
    Score.row(i) = val;
    Lik(i) = alpha;
  }
  vecmat res; 
  res.V = Lik; res.M = Score;
  return(res);
}
// }}} score 2


// {{{ biprobit2

RcppExport SEXP biprobit2(SEXP m, SEXP dm,
			  SEXP s, SEXP ds, SEXP y, SEXP w, SEXP std,
			  SEXP eqweight, SEXP eqmarg, SEXP correlation) {
BEGIN_RCPP
  NumericMatrix mm(m); 
  NumericMatrix ss(s); 
  NumericMatrix dss(ds); 
  NumericMatrix dmm(dm);
  NumericMatrix yy(y);

  mat mu(mm.begin(), mm.nrow(), mm.ncol(), false);
  mat S(ss.begin(), ss.nrow(), ss.ncol(), false);
  mat dS(dss.begin(), dss.nrow(), dss.ncol(), false);
  mat dmu(dmm.begin(), dmm.nrow(), dmm.ncol(), false);
  mat Y(yy.begin(), yy.nrow(), yy.ncol(), false);
  bool weights = Rcpp::as<bool>(std);
  bool OneWeight = Rcpp::as<bool>(eqweight); 
  bool EqMarginal = Rcpp::as<bool>(eqmarg);
  bool Cor = Rcpp::as<bool>(correlation);

  vecmat U0;
  mat W,U;
  List res;
  if (weights) {
    NumericMatrix ww(w);
    mat W0(ww.begin(), ww.nrow(), ww.ncol(), false);
    W = W0;
    //    WW = trans(reshape(W1,2,N/2));    
  }
  if (OneWeight & weights) {
    if (Cor) { 
      U0 = scorecor(Y,mu,dmu,S,dS,U,EqMarginal);
    } else {
      U0 = score2(Y,mu,dmu,S,dS,U,EqMarginal);
    }
    res["loglik"] = log(U0.V)%W.col(0);
    for (unsigned i=0; i<U0.M.n_rows; i++) U0.M.row(i) *= W(i,0);
    res["score"] = U0.M;
  } else {
    if (Cor) {
      U0 = scorecor(Y,mu,dmu,S,dS,W,EqMarginal);
    } else {
      U0 = score2(Y,mu,dmu,S,dS,W,EqMarginal);
    }
    res["loglik"] = log(U0.V);
    res["score"] = U0.M;
  }

  return(res);
END_RCPP
}

// }}} biprobit 2


