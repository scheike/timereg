#include "mvn.h"
#include "tools.h"
#include <math.h>

int _mvt_maxpts=25000; 
double _mvt_abseps=0.001; 
double _mvt_releps=0;  
int _mvt_df = 0;
int _mvt_inform;
double _mvt_error[3]; 


double dmvn(const vec &y, const mat &W, 
	    bool log=true, double logdet=datum::inf) {
  int n = W.n_rows;  
  double res = -0.5*n*log2pi;
  if (logdet!=datum::inf) {
    res += -0.5*(logdet + as_scalar(trans(y)*W*y));    
  } else {
    double sign=0;
    mat iW = inv(W);
    log_det(logdet,sign,W);
    res += -0.5*(logdet + as_scalar(trans(y)*iW*y));  
  }
  if (!log) res = exp(res);
  return(res);  
}
	    

double cdfmvn(mat &upper, mat &cor) {  
  double val=0;
  int n = cor.n_cols;
  NumericVector _mvt_delta(n);   
  unsigned ncor = n*(n-1)/2;
  rowvec Cor(ncor);
  int j = 0;
  for (int r=0; r<n; r++) {
    for (int c=r+1; c<n; c++) {
      Cor(j) = cor(r,c); 
      j++;
    }
  }
  irowvec infin(n); infin.fill(0); //
  mvtdst_(&n, &_mvt_df,
  	  &upper[0], // Lower, ignored
  	  &upper[0], 
  	  &infin[0], // Infinity argument (all 0 since CDF)
  	  &Cor[0],
  	  &_mvt_delta[0], &_mvt_maxpts,
  	  &_mvt_abseps, &_mvt_releps,
  	  &_mvt_error[0], &val, &_mvt_inform);  
  return(val);
}



RcppExport SEXP Dpmvn (SEXP lower, SEXP upper, SEXP mu, SEXP sigma, SEXP std) { 
BEGIN_RCPP
  colvec y = Rcpp::as<colvec>(upper);
  NumericVector Lower(lower);
  bool Std = Rcpp::as<bool>(std);
  mat S = Rcpp::as<mat>(sigma);  
  colvec Mu = Rcpp::as<colvec>(mu);  
  int n = S.n_cols;
  vec L(n);
  for (int j=0; j<n; j++) {    
    if (y(j)==datum::inf) L(j) = 1; // (lower,+Inf)
    if (Lower(j)==-datum::inf) L(j) = 0; // (-Inf,upper)
    // =2 (low,up)   For now we actually don't consider the interval-censored case
    // <0 (-Inf,Inf)
  }  
  mat LL = diagmat(L);
  mat iLambda;
  if (!Std) {
    iLambda = diagmat(1/sqrt(diagvec(S)));
    S = iLambda*S*iLambda;    
    y = iLambda*(y-Mu);
  }
  
  List res;

  vec D(n);  
  mat H(n,n);
  if (n==1) {
    List res;    
    D[0] = Rf_dnorm4(y[0],0.0,1.0,0);
    H[0] = -y[0]*D[0];
    if (!Std) {
      D[0] *= iLambda[0];
      H[0] *= iLambda[0]*iLambda[0];
    }
    res["gradient"] = D;
    res["hessian"] = H;
    return(res);    
  }

  for (unsigned j=0; j<n; j++) {
    mat Sj = S; Sj.shed_row(j); 
    mat Sj0 = Sj.col(j);
    Sj.shed_col(j);
    Sj -= Sj0*Sj0.t();
    mat iL = diagmat(1/sqrt(diagvec(Sj)));    
    Sj = iL*Sj*iL;
    vec muj = y;  
    muj.shed_row(j);
    muj = iL*(muj - Sj0*y[j]);
    D[j] = Rf_dnorm4(y[j],0.0,1.0,0)*
      cdfmvn(muj,Sj);
  }
  
  if (n==2) {
    H(0,1) = H(1,0) = dmvn(y, S, false);
    H.diag() = -y%D - H(0,1)*S(0,1);
  } else {
    mat phis(n,n); phis.fill(0);
    mat Phis = phis;
    uvec idx1(n-2);
    uvec idx2(2);
    for (unsigned i=0; i<(n-1); i++) {
      for (unsigned j=(i+1); j<n; j++) {
	idx2(0) = i; idx2(1) = j;
	unsigned pos = 0;
	for (unsigned k=0; k<n; k++) {
	  if (k!=i && k!=j) {
	    idx1(pos) = k; 
	    pos++;
	  }
	}
	mat Snij = S.submat(idx1,idx2);
	mat S2 = S.submat(idx2,idx2);
	mat B = Snij*inv(S2);
	mat Sij = S.submat(idx1,idx1) - B*trans(Snij);
	vec muij = y.elem(idx1) - B*y.elem(idx2);

	mat iL = diagmat(1/sqrt(diagvec(Sij)));    
	Sij = iL*Sij*iL;
	muij = iL*muij;
	Phis(i,j) = Phis(j,i) = cdfmvn(muij,Sij);
	phis(i,j) = phis(j,i) = dmvn(y.elem(idx2),S2,false);
      }
    }
    H = Phis%phis;    
    colvec ones(n); ones.fill(1);
    H.diag() = -y%D - (S%H)*ones;
  }

  //double val = cdfmvn(y,S);  
  //return(Rcpp::wrap(val));

  if (!Std) {
    H = iLambda*H*iLambda;    
    D = iLambda*D;
  }
  res["gradient"] = D;
  res["hessian"] = H;

  return(res);
END_RCPP
}





vec loglikmvn(mat &Yl, mat &Yu, uvec &Status,
	      mat &Mu, mat &S, 
	      mat &Z,  mat &Su,
	      mat &Threshold) {
  
  int k = Yl.n_cols;
  int n = Yl.n_rows;
  uvec Obs = find(Status==0);
  uvec Cens = find(Status==1);
  uvec Ord = find(Status>1);
  uvec NonObs = find(Status>0);
  int nObs = Obs.size();
  int nNonObs = NonObs.size();
  int nOrd = Ord.size();
  int nu = Su.n_cols;
  bool nonconstvar = (nu>0); 

  double sign, logdetS0;
  mat Se,S0,iS0;
    
  vec loglik(n); loglik.fill(0);
  if (nObs>0) {
    mat Y0 = Yl.cols(Obs)-Mu.cols(Obs);
    Se = S0 = S.submat(Obs,Obs); 
    iS0 = inv(S0);
    log_det(logdetS0,sign,S0);
    double normconst = -0.5*nObs*log2pi;  

    for (int i=0; i<n; i++) { // Iterate over subjects

      if (nonconstvar) {
	mat Z0 = reshape(Z.row(i),k,nu);	
	S0 = Se+Z0.rows(Obs)*Su*trans(Z0.rows(Obs));
	iS0 = inv(S0);
	log_det(logdetS0,sign,S0);
      }
      
      loglik(i) = -0.5*(logdetS0 + as_scalar(Y0.row(i)*iS0*trans(Y0.row(i))));
      loglik(i) += normconst;
    }
  }


  if (nNonObs>0) {
    NumericVector _mvt_delta(nNonObs); 

    mat MuNonObs = Mu.cols(NonObs);
    mat SNonObs = S.submat(NonObs,NonObs);

    if ((nObs>0) & (!nonconstvar)) { // Calculate conditional on observed
      mat S01 = S.submat(NonObs,Obs);
      mat iS1 = inv(S.submat(Obs,Obs));
      MuNonObs = MuNonObs + 
	trans(S01*iS1*trans(Yl.cols(Obs)-Mu.cols(Obs)));   
      //      MuNonObs.each_row() += Mu.(NonObs);
      SNonObs = SNonObs - S01*iS1*trans(S01);

    }
    vec il = 1/sqrt(diagvec(SNonObs));
    mat iL = diagmat(il);
    
    Se = S0 = iL*SNonObs*iL; // Correlation matrix
    int ncor = nNonObs*(nNonObs-1)/2;
    rowvec Cor(ncor);
    if (ncor>0) {
      int j = 0;
      for (int r=0; r<nNonObs; r++) {
	for (int c=r+1; c<nNonObs; c++) {
	  Cor(j) = S0(r,c); 
	  j++;
	}
      }
    }
    int nthresmax = Threshold.n_cols;
    uvec OrdNonObs;

    mat Thres;
    if (nOrd>0) { 
      OrdNonObs = find(Status.elem(NonObs)>1);
      Thres = Threshold.rows(Ord);
    }
   
    rowvec lower(nNonObs); 
    rowvec upper(nNonObs);
    irowvec infin(nNonObs); // 0: right, 1: left, 2: interval

    uvec currow(1); 
    for (int i=0; i<n; i++) { // Iterate over subjects

      currow(0) = i;
      rowvec Mi = MuNonObs.row(i);
      lower = Yl.submat(currow,NonObs);
      upper = Yu.submat(currow,NonObs);      
	
      if (nonconstvar) {
	mat Z0 = reshape(Z.row(i),k,nu);	
	mat SS = S+Z0*Su*trans(Z0);

	if (nObs>0) {
	  S0 = SS.submat(NonObs,NonObs);
	  mat S01 = SS.submat(NonObs,Obs);
	  mat iS1 = inv(SS.submat(Obs,Obs));	  
	  Mi = Mi + 
	    trans(S01*iS1*trans(Yl.submat(currow,Obs)-Mu.submat(currow,Obs)));   

	  SNonObs = S0 - S01*iS1*trans(S01);
	} else {
	  SNonObs = SS;
	}

	il = 1/sqrt(diagvec(SNonObs));
	iL = diagmat(il);    
	Se = S0 = iL*SNonObs*iL; // Correlation matrix
	if (ncor>0) {
	  int j = 0;
	  for (int r=0; r<nNonObs; r++) {
	    for (int c=r+1; c<nNonObs; c++) {
	      Cor(j) = S0(r,c); 
	      j++;
	    }
	  }
	}
      }

      umat StatusNonObs = Status.elem(NonObs);
      infin.fill(2);
      
      for (int j=0; j<nNonObs; j++) {
	if (upper(j)==datum::inf && StatusNonObs(j)==1) infin(j) = 1;
	if (lower(j)==-datum::inf && StatusNonObs(j)==1) infin(j) = 0;
      }
      // uvec infplus = find(upper==datum::inf);
      // uvec infminus = find(lower==-datum::inf);
      // if (infplus.size()>0) { infin.elem(infplus) -= 1; }
      // if (infminus.size()>0) { infin.elem(infminus) -= 2; }
          
      if (nOrd>0) {	
	for (int j=0; j<nOrd; j++) {
	  int jj = OrdNonObs(j);
	  int yval = (int) lower(jj);
	  if (yval<1) { // if Y=0
	    infin(jj) = 0; // Integrate over left tail
	    upper(jj) = Thres(j,0);
	  } else { // Y>1
	    if (yval>=nthresmax) { // Y=k (last)
	      double val = Thres(j,yval-1);
	      infin(jj) = 1; // Integrate over right tail
	      lower(jj) = val;
	    } else {
	      double val = Thres(j,yval-1);
	      double val2 = Thres(j,yval);
	      if (val>=val2) { // Also Y=k (last)
		infin(jj) = 1; 
		lower(jj) = val;
	      } else { //Y=k-i (in between)	      
		lower(jj) = val;
		upper(jj) = val2;
	      } 
	    }
	  }
	}
      }

      lower = (lower-Mi)%trans(il);
      upper = (upper-Mi)%trans(il);

      double val = 0;
      mvtdst_(&nNonObs, &_mvt_df,
	      &lower[0], &upper[0], 
	      &infin[0], &Cor[0],
	      &_mvt_delta[0], &_mvt_maxpts,
	      &_mvt_abseps, &_mvt_releps,
	      &_mvt_error[0], &val, &_mvt_inform);        

      // if (isnan(val)) {
      // 	cerr << "***i=" << i << endl;
      // 	cerr << "Threshold=" << Threshold << endl;
      // 	cerr << "Thres=" << Thres << endl;
      // 	cerr << "status=" << Status;
      // 	cerr << "yl=" << Yl.row(i);
      // 	cerr << "yu=" << Yu.row(i);
      // 	cerr << "lower=" << lower;
      // 	cerr << "upper=" << upper;
      // 	cerr << "infin=" << infin;
      // 	cerr << "Cor=\n" << Cor;
      // 	cerr << "   val=" << val << endl;
      // }
      
      loglik(i) += log(val);
    }
  }
  
  return(loglik);
}

RcppExport SEXP loglikMVN(SEXP yl, SEXP yu, 
			  SEXP status,
			  SEXP mu, SEXP dmu,
			  SEXP s, SEXP ds,
			  SEXP z, SEXP su, SEXP dsu,
			  SEXP threshold, SEXP dthreshold) {
BEGIN_RCPP
  mat Yl = Rcpp::as<mat>(yl);
  uvec Status = Rcpp::as<uvec>(status);
  mat Yu = Rcpp::as<mat>(yu);
  mat Mu = Rcpp::as<mat>(mu);  
  mat S = Rcpp::as<mat>(s);  

  if ((Mu.n_cols!=Yl.n_cols) || (Mu.n_rows!=Yl.n_rows))
      throw(Rcpp::exception("Dimension of 'mu' and 'yl' did not agree","mvn.cpp",1));
  if (Status.size()!=Yl.n_cols)
      throw(Rcpp::exception("Dimension of 'status' and 'yl' did not agree","mvn.cpp",1));

  uvec Cens = find(Status==1);
  uvec Obs = find(Status==0);
  uvec NonObs = find(Status>0);
  uvec Ord = find(Status>1);
  // int nObs = Obs.size();
  // int nNonObs = NonObs.size();
  int nOrd = Ord.size();
  int nCens = Cens.size();
  unsigned n = Yl.n_rows;
  mat Z,Zsub;
  mat Su;
  if (!Rf_isNull(z)) {
    Z = Rcpp::as<mat>(z);
    Su = Rcpp::as<mat>(su);
    if (Z.n_cols!=(Yl.n_cols/Su.n_cols))
      throw(Rcpp::exception("Dimension of 'z' and 'su' did not agree","mvn.cpp",1));
    if (Z.n_rows!=n) 
      throw(Rcpp::exception("Dimension of 'z' and 'yl' did not agree","mvn.cpp",1));	   
  }
  mat Threshold;
  if (nOrd>0) {
    Threshold = Rcpp::as<mat>(threshold);  
    if (Threshold.n_rows!=Yl.n_cols) 
      throw(Rcpp::exception("Dimension of 'threshold' and 'yl' did not agree","mvn.cpp",1));
  }

  vec loglik(n); loglik.fill(0);
  if (nCens>0) {
    if ((Yl.n_cols!=Yu.n_cols) || (Yl.n_rows!=Yu.n_rows)) 
      throw(Rcpp::exception("Dimension of 'yl' and 'yu' did not agree","mvn.cpp",1));

    umat stat = (Yl.cols(Cens)==Yu.cols(Cens));
    uvec group(n);
    umat pattern; 
    fastpattern(stat,pattern,group);    
    uvec NewStatus = Status;        
    unsigned K = pattern.n_rows;
    for (unsigned i=0; i<K; i++) { // Iterate over patterns
      uvec idx = find(group==i);
      NewStatus.elem(Cens) = 1-pattern.row(i);      
      mat Ylsub = Yl.rows(idx);
      mat Yusub = Yu.rows(idx);
      mat Musub = Mu.rows(idx);
      if (!Rf_isNull(z)) Zsub = Z.rows(idx);
      // cerr << "i=" << i << endl;
      // cerr << "idx =" << idx << endl; 
      // cerr << "Yl =" << Ylsub << endl;
      // cerr << "Yu =" << Yusub << endl;
      // cerr << "Musub =" << Musub << endl;
      // cerr << "Zsub =" << Zzub << endl;
      // cerr << "Status =" << NewStatus << endl;
      
      vec ll = loglikmvn(Ylsub, Yusub, NewStatus,
				   Musub, S, 
				   Zsub,  Su,
				   Threshold);      
      loglik.elem(idx) = ll;
    }   
    

  } else {
    loglik = loglikmvn(Yl,Yu,Status,
		       Mu, S, 
		       Z,  Su,
		       Threshold);
  }

  return (wrap(loglik));
END_RCPP
}
//   return(Rcpp::List::create(
// 			    Rcpp::Named("loglik")=loglik,
// 			    Rcpp::Named("k")=k,
// 			    Rcpp::Named("n")=n
// 			    ));
// }





RcppExport SEXP pmvn(SEXP upper, SEXP cor) { 
  NumericVector Upper(upper); int n = Upper.size();
  NumericVector Lower(n);
  IntegerVector infin(n); 
            // INFIN(I) < 0, Ith limits are (-infinity, infinity);
            // INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
            // INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
            // INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
  NumericVector Cor(cor);
  NumericVector _mvt_delta(n); 

  double val;
  mvtdst_(&n, &_mvt_df,
	  &Lower[0], &Upper[0], 
	  &infin[0], &Cor[0],
	  &_mvt_delta[0], &_mvt_maxpts,
	  &_mvt_abseps, &_mvt_releps,
	  &_mvt_error[0], &val, &_mvt_inform);  
  List res;
  res["val"] = val;
  return(res);
}


RcppExport SEXP bvncdf(SEXP a, SEXP b, SEXP r) { 
  double u1=-Rcpp::as<double>(a);
  double u2=-Rcpp::as<double>(b);
  double cr=Rcpp::as<double>(r);  
  double val = bvnd_(&u1, &u2, &cr);
  NumericVector res(1); res[0] = val;
  return(res);
}

double dbvnorm(double y1, double y2, double R) {
  double detR = 1-R*R;
  // inv(R) = [1 -r; -r 1]/detR (prove by gauss elim.)
  double res = 1/(2*M_PI*sqrt(detR))*exp(-0.5/detR*(y1*y1+y2*y2-2*R*y1*y2));
  return(res);
}

vecmat Dbvn(double y1, double y2, double R) {
  vec DP(2);
  double R2 = R*R;
  DP(0) = Rf_dnorm4(y1,0.0,1.0,0)*Rf_pnorm5(y2-R*y1,0.0,sqrt(1-R2),1,0);
  DP(1) = Rf_dnorm4(y2,0.0,1.0,0)*Rf_pnorm5(y1-R*y2,0.0,sqrt(1-R2),1,0);
  mat HP(2,2);
  HP(1,0) = HP(0,1) = dbvnorm(y1,y2,R);
  HP(0,0) = -y1*DP(0) - R*HP(1,0);
  HP(1,1) = -y2*DP(1) - R*HP(1,0);  
  vecmat res;
  res.V = DP;
  res.M= HP;
  return(res);
}

double Sbvn(double &l1, double &l2,double &r) {     
  double val = bvnd_(&l1, &l2, &r);
  return(val);
}

