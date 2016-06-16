#include "tools.h"
#include <vector>

SEXP FastLong2(SEXP idata, SEXP inclust,
	       SEXP infixed, SEXP invarying) {
BEGIN_RCPP
  unsigned nvarying = Rcpp::as<unsigned>(invarying); // Number of varying var
  unsigned nfixed = Rcpp::as<unsigned>(infixed); // Number of non-varying
  unsigned nclust = Rcpp::as<unsigned>(inclust); // Number within cluster      
  Rcpp::DataFrame d = Rcpp::as<Rcpp::DataFrame>(idata); // Wide data
  //bool Missing = Rcpp::as<bool>(missing);
  unsigned n= d.nrows();
  unsigned M = nclust*n;
  unsigned K = nvarying+nfixed+2;

  Row<unsigned> idx(nvarying);
  uvec mis(M); mis.fill(0);
  for (unsigned k=0; k<nvarying; k++)
    idx[k] = nfixed+k*nclust;
  Rcpp::List myList(K);
     
  for (unsigned k=0; k<nfixed; k++) {
    if (Rf_isLogical(d[k])) {
      LogicalVector mm = Rcpp::as<LogicalVector>(d[k]);
      myList[k] = Rcpp::rep_each(mm,nclust);
    } else if (Rf_isInteger(d[k])) {
      IntegerVector mm = Rcpp::as<IntegerVector>(d[k]);	  
      myList[k] = Rcpp::rep_each(mm,nclust);
    } else if (Rf_isNumeric(d[k])) {
      NumericVector mm = Rcpp::as<NumericVector>(d[k]);	  
      myList[k] = Rcpp::rep_each(mm,nclust);
    } else if (Rf_isComplex(d[k])) {
      ComplexVector mm = Rcpp::as<ComplexVector>(d[k]);	  
      myList[k] = Rcpp::rep_each(mm,nclust);
    } else if (Rf_isString(d[k]) || Rf_isFactor(d[k])) {
      CharacterVector mm = Rcpp::as<CharacterVector>(d[k]);	  
      myList[k] = Rcpp::rep_each(mm,nclust);
    }
  }

  for (unsigned k=0; k<nvarying; k++) {
    unsigned type=0; // 1:logical,2:integer,3:numeric,4:complex,5:character
    for (unsigned i=0; i<nclust; i++) {
      bool assigned = false;
      if (Rf_isLogical(d[idx[k]+i]) & (type<2)) {
	type=1;
	assigned=true;
      }
      if (Rf_isInteger(d[idx[k]+i]) & (type<3)) {
	type=2;
	assigned=true;
      }
      if (Rf_isNumeric(d[idx[k]+i]) & (type<4)) {
	type=3;
	assigned=true;
      }
      if (Rf_isComplex(d[idx[k]+i]) & (type<5)) {
	type=4;
	assigned=true;
      }
      if (!assigned) {
	//Rf_isString(d[idx[k]]) || Rf_isFactor(d[idx[k]])
	type=5;
      }
    }
    if (type==1) {
      LogicalMatrix mm(nclust,n);
      for (unsigned i=0; i<nclust; i++)
	mm(i,_) = Rcpp::as<LogicalVector>(d[idx[k]+i]);	    
      myList[nfixed+k] = Rcpp::LogicalVector(mm.begin(),mm.end());
    } else if (type==2) {
      IntegerMatrix mm(nclust,n);
      for (unsigned i=0; i<nclust; i++)
	mm(i,_) = Rcpp::as<IntegerVector>(d[idx[k]+i]);	    
      myList[nfixed+k] = Rcpp::IntegerVector(mm.begin(),mm.end());
    } else if (type==3) {
      NumericMatrix mm(nclust,n);
      for (unsigned i=0; i<nclust; i++)
	mm(i,_) = Rcpp::as<NumericVector>(d[idx[k]+i]);	    
      //myList[nfixed+k] = Rcpp::NumericMatrix(M,1,mm.begin());      
      myList[nfixed+k] = Rcpp::NumericVector(mm.begin(),mm.end());
    } else if (type==4) {
      ComplexMatrix mm(nclust,n);
      for (unsigned i=0; i<nclust; i++)
	mm(i,_) = Rcpp::as<ComplexVector>(d[idx[k]+i]);	    
      myList[nfixed+k] = Rcpp::ComplexVector(mm.begin(),mm.end());
    } else if (type==5) {
      CharacterMatrix mm(nclust,n);
      for (unsigned i=0; i<nclust; i++)
	mm(i,_) = Rcpp::as<CharacterVector>(d[idx[k]+i]);	    
      myList[nfixed+k] = Rcpp::CharacterVector(mm.begin(),mm.end());
    } 
  }

  //  IntegerVector Id = Rcpp::seq_len(n);
  IntegerVector Id = Rcpp::rep_each(Rcpp::seq_len(n),nclust);
  IntegerVector Num = Rcpp::rep(Rcpp::seq_len(nclust),n);  
  myList[K-2] = Id;
  myList[K-1] = Num;
  myList.attr("names") = Rcpp::seq_len(K);
  myList.attr("row.names") = Rcpp::seq_len(M);
  myList.attr("class") = "data.frame";

  return(Rcpp::wrap(myList));
END_RCPP
}


SEXP FastLong(SEXP idata, SEXP inclust,
	      SEXP infixed, SEXP invarying, SEXP missing) {
BEGIN_RCPP
      unsigned nvarying = Rcpp::as<unsigned>(invarying); // Number of varying var
      unsigned nfixed = Rcpp::as<unsigned>(infixed); // Number of non-varying
      unsigned nclust = Rcpp::as<unsigned>(inclust); // Number within cluster
      mat d = Rcpp::as<mat>(idata); // Wide data
      bool Missing = Rcpp::as<bool>(missing);
      unsigned M = nclust*d.n_rows;
      unsigned K = nvarying+nfixed+2;
      mat dd(M,K); // Long data      
      uvec idx(nvarying);
      uvec mis(M); mis.fill(0); // NA_INTEGER
      for (unsigned k=0; k<nvarying; k++)
         idx[k] = nfixed+k*nclust;
      rowvec xx(K); xx.fill(NA_REAL);      
      unsigned count = 0;
      for (unsigned i=0; i<d.n_rows; i++) {
         rowvec xx0 = xx;
         rowvec d0 = d.row(i);
         for (unsigned k=0; k<nfixed; k++) {
	   xx0[k] = d(i,k);
	 }
         for (unsigned j=0; j<nclust; j++) {
	   uvec idx0 = idx+j;
            xx0.subvec(nfixed,K-3) = trans(d0.elem(idx0));
            xx0[K-2] = i+1; xx0[K-1] = j+1;
	    bool allmiss = Missing;
	    if (Missing) 
	      for (unsigned k=nfixed*(i>0); k<nvarying+nfixed; k++) {
		if (!R_IsNA(xx0[k])) {
		  allmiss = false;		 
		}
		if (!allmiss) break;
	      }
	    if (allmiss) {
	      mis[count] = 1;
	    } else {
	      dd.row(count) = xx0;
	    }
            count++;
         }
      }
      if (Missing) dd = dd.rows(find(mis==0));
      return(Rcpp::wrap(dd));
END_RCPP
}

void fastpattern(const umat &y, umat &pattern, uvec &group, unsigned categories /*=2*/) {
  unsigned n = y.n_rows;
  unsigned k = y.n_cols;

  uvec mygroup(n);
  unsigned npattern = (unsigned) pow((double) categories,(double) k);
  umat mypattern(npattern,k);
  mypattern.fill(1);
  unsigned K=0;

  for (unsigned i=0; i<n; i++) {
    urowvec ri = y.row(i);
    bool found = FALSE;
    for (unsigned j=0; j<K; j++) {      
      if (sum(ri!=mypattern.row(j))==0) {
	found = true;
	mygroup(i) = j;	
	break;
      }

    }
    if (!found) {
      K++;
      mypattern.row(K-1) = ri;
      mygroup(i) = K-1;
    }
  } 
  pattern = mypattern.rows(0,K-1);
  group = mygroup;
}

SEXP FastPattern(SEXP y1,SEXP y2,SEXP cat)  {
BEGIN_RCPP
  mat Y1 = Rcpp::as<mat>(y1);
  unsigned n = Y1.n_rows;
  unsigned k = Y1.n_cols;
  unsigned Cat = Rcpp::as<unsigned>(cat);
  umat stat;
  if (!Rf_isNull(y2)) {
    mat Y2 = Rcpp::as<mat>(y2);   
    if ((Y2.n_rows!=n) || (Y2.n_cols!=k)) 
       throw(Rcpp::exception("Dimension did not agree","tools.cpp",1));
    stat = (Y1==Y2);
  } else {
    stat = conv_to<umat>::from(Y1);
  }
  uvec group(n);
  //  unsigned npattern = (unsigned) pow(2,(double) k);
  umat pattern; //(npattern,k);  
  fastpattern(stat,pattern,group,Cat);
  
  return(Rcpp::List::create(
			    Rcpp::Named("pattern")=pattern,
			    Rcpp::Named("group")=wrap(group)
			    ));
END_RCPP
} 


RcppExport SEXP FastCluster(const SEXP x) {
BEGIN_RCPP
  arma::Row<unsigned> xx = as<arma::Row<unsigned> >(x);
  vector<unsigned> cpos, csize; 
  unsigned val=0, prev=0, cursize=0;  
  for (unsigned i=0; i<xx.n_elem; i++) {
    if (xx[i]!=0) {
      if (prev==0) {
	cpos.push_back(i+1);
	val++;
      }
      xx[i] = val;
      cursize++;
    } else if (prev!=0) {
      csize.push_back(cursize);
      cursize = 0;
    }
    prev = xx[i];
  }
  if (prev!=0) csize.push_back(cursize);
  
  return(Rcpp::List::create(
			    Rcpp::Named("cluster")=conv_to< vector<unsigned> >::from(xx),
			    Rcpp::Named("cluster.first")=cpos,
			    Rcpp::Named("cluster.size")=csize
			    ));  
  // return( wrap(list(xx,) );
END_RCPP
}
		     
RcppExport SEXP FastApprox(const SEXP time,
			   const SEXP newtime,
			   const SEXP equal,
			   const SEXP type // (0: nearest, 1: right, 2: left)
			   ) {
BEGIN_RCPP
  unsigned Type = Rcpp::as<unsigned>(type);
  NumericVector NewTime(newtime);
  NumericVector Sorted(time);
  // IntegerVector Order;
  // NumericVector Sorted = Time;
  // std::sort(Sorted.begin(), Sorted.end());
  // //    .sort();
  // IntegerVector Order = match(Sorted, Time);
  // return(Rcpp::wrap(Time));  
  bool Equal = Rcpp::as<bool>(equal);
  vector<int> idx(NewTime.size());
  vector<int> eq(NewTime.size());

  double vmax = Sorted[Sorted.size()-1];
  NumericVector::iterator it;  
  double upper=0.0; int pos=0;
  for (int i=0; i<NewTime.size(); i++) {    
    eq[i] = 0;
    if (NewTime[i]>vmax) {
      pos = Sorted.size()-1;
    } else {
      it = lower_bound(Sorted.begin(), Sorted.end(), NewTime[i]);
      upper = *it;
      if (it == Sorted.begin()) { 
	pos = 0; 
	if (Equal && (NewTime[i]==upper)) { eq[i] = 1; }
      }
      // else if (int(it-Sorted.end())==0) {
      // 	pos = Sorted.size()-1;
      // }
      else {
	pos = int(it-Sorted.begin());
	if (Type==0 && fabs(NewTime[i]-Sorted[pos-1])<fabs(NewTime[i]-Sorted[pos])) pos -= 1;
	if (Equal && (NewTime[i]==upper)) { eq[i] = pos+1; }
      }
    }
    if (Type==2 && NewTime[i]<upper) pos--;
    idx[i] = pos+1;
  }
  if (Equal) {
    List ans;
    ans["idx"] = idx;
    ans["eq"] = eq;    
    return(ans);   
  } 
  return(Rcpp::wrap(idx));
END_RCPP
}
