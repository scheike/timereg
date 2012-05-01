#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace arma;
using namespace Rcpp;

double plack(double theta,double cif1,double cif2,vec &dp) 
{ // {{{
double valr,valn,val1,cifs,
       thetad,val1d,valnd,valrd,d,
       cif1d, cif2d, cifsd;

cifs=cif1+cif2; // {{{
if (theta!=1) {
valn=2*(theta-1); 
val1=(1+(theta-1)*(cifs))-pow( pow((1+(theta-1)*cifs),2)-4*cif1*cif2*theta*(theta-1),0.5); 
valr=val1/valn; 
} else {
valr=cif1*cif2;
} // }}}

d=0.000001; thetad=theta+d; // {{{
if (thetad!=1) {
valnd=2*(thetad-1); 
val1d=(1+(thetad-1)*(cifs))-pow( pow((1+(thetad-1)*cifs),2)-4*cif1*cif2*thetad*(thetad-1),0.5); 
valrd=val1d/valnd; 
} else {
valrd=cif1*cif2;
} // }}}
dp[0]=(valrd-valr)/d;  

cif1d=cif1+d; cifsd=cif1d+cif2; // {{{
if (theta!=1) {
valnd=2*(theta-1); 
val1d=(1+(theta-1)*(cifsd))-pow( pow((1+(theta-1)*cifsd),2)-4*cif1d*cif2*theta*(theta-1),0.5); 
valrd=val1d/valnd; 
} else {
valrd=cif1d*cif2;
} // }}}
dp[1]=(valrd-valr)/d;  

cif2d=cif2+d; cifsd=cif1+cif2d; // {{{
if (theta!=1) {
valnd=2*(theta-1); 
val1d=(1+(theta-1)*(cifsd))-pow( pow((1+(theta-1)*cifsd),2)-4*cif1d*cif2*theta*(theta-1),0.5); 
valrd=val1d/valnd; 
} else {
valrd=cif1d*cif2;
} // }}}
dp[2]=(valrd-valr)/d;  

//if (theta!=1) {
//dval1= cifs-(2*(1+(theta-1)*cifs)*cifs-4*2*cif1*cif2*theta+4*cif1*cif2)/
//	(2*pow( pow((1+(theta-1)*cifs),2)-4*cif1*cif2*theta*(theta-1),0.5)); 
//val=valn*dval1-val1*2; 
//dp[0]= val/pow(valn,2); 
//dp[0]=(valrd-valr)/0.000001; 
//} else {
//dp[0]=1; 
//}

return(valr); 
} // }}}
               
double min(double a, double b) { if (a<b) return(a); else return(b); }
double max(double a, double b) { if (a>b) return(a); else return(b); }

RcppExport SEXP cor(SEXP itimes,SEXP iy,SEXP icause, SEXP iCA1, SEXP iKMc,
		SEXP iz, SEXP iest,SEXP iZgamma, SEXP isemi,SEXP izsem,
//		detail,biid,gamiid,timepow,theta,vartheta,
		SEXP itheta, SEXP iXtheta, 
		SEXP ithetades,
		SEXP icluster,SEXP iclustsize,SEXP iclusterindex,
		SEXP iinverse,SEXP iCA2,
		SEXP ix2, // SEXP iz2,
		SEXP isemi2, SEXP iest2,SEXP iZ2gamma2,
//		b2iid, gam2iid, SEXP htheta,SEXP dhtheta,SEXP rhoR,
	        SEXP iflexfunc,
		SEXP iiid, 
		SEXP isym,SEXP  iweights, 
                SEXP isamecens, SEXP istabcens,SEXP iKMtimes,SEXP isilent,SEXP icifmodel,
		SEXP idepmodel,
		SEXP iestimator, SEXP ientryage,SEXP icif1entry,SEXP icif2entry,SEXP itrunkp
) // {{{
{
// {{{ setting matrices and vectors, and exporting to armadillo matrices
//      NumericMatrix zz(z); 
//	mat mz(zz.begin(),zz.nrow(),zz.ncol(),false); 
//mat z2 = Rcpp::as<mat>(iz2);
mat est = Rcpp::as<mat>(iest);
mat est2 = Rcpp::as<mat>(iest2);
mat z = Rcpp::as<mat>(iz);
mat zsem = Rcpp::as<mat>(izsem);
mat z2 = Rcpp::as<mat>(ix2);
mat thetades = Rcpp::as<mat>(ithetades);
mat clusterindex = Rcpp::as<mat>(iclusterindex);
colvec theta = Rcpp::as<colvec>(itheta);
colvec Xtheta = Rcpp::as<colvec>(iXtheta);
colvec y = Rcpp::as<colvec>(iy);
colvec clustsize = Rcpp::as<colvec>(iclustsize);
int antclust = clusterindex.n_rows; 
 colvec times = Rcpp::as<colvec>(itimes);
 int Ntimes=times.n_rows; 
 colvec cause = Rcpp::as<colvec>(icause);
 colvec cluster = Rcpp::as<colvec>(icluster);
 colvec Zgamma = Rcpp::as<colvec>(iZgamma);
 colvec KMtimes  = Rcpp::as<colvec>(iKMtimes );
 colvec Z2gamma2 = Rcpp::as<colvec>(iZ2gamma2);
 colvec KMc= Rcpp::as<colvec>(iKMc);
 colvec weights = Rcpp::as<colvec>(iweights);
 colvec entryage = Rcpp::as<colvec>(ientryage);
 colvec cif1entry = Rcpp::as<colvec>(icif1entry);
 colvec cif2entry = Rcpp::as<colvec>(icif2entry);
 colvec trunkp = Rcpp::as<colvec>(itrunkp);

  int samecens = Rcpp::as<int>(isamecens);
  int inverse= Rcpp::as<int>(iinverse);
  int semi = Rcpp::as<int>(isemi);
  int semi2 = Rcpp::as<int>(isemi2);
  int flexfunc = Rcpp::as<int>(iflexfunc);
  int stabcens = Rcpp::as<int>(istabcens);
  int silent = Rcpp::as<int>(isilent);
  int cifmodel = Rcpp::as<int>(icifmodel);
  int CA1 = Rcpp::as<int>(iCA1); 
  int CA2 = Rcpp::as<int>(iCA2);
  int sym = Rcpp::as<int>(isym); 
  int depmodel= Rcpp::as<int>(idepmodel); 
  int estimator= Rcpp::as<int>(iestimator); 
  int iid= Rcpp::as<int>(iiid); 

  int ci,ck,i,j,c,s,k,v,naprint; 
  double Li,Lk,weight,p11t,ormarg,sdj,diff,cweight1,cweight2,resp3,time,resp1,resp2,response,thetak,respst; 
//  double plack(); 
  vec dplack(4); 
  dplack=0*dplack; 
  int pt=theta.n_rows; 

  mat thetiid(antclust,pt); 
  if (iid==1) thetiid=0*thetiid; 

  vec Utheta(pt); 
  vec vthetascore(pt); 
  rowvec pthetavec(pt); 
  vec vtheta2(pt); 
  mat DUtheta(pt,pt); 
  DUtheta=0*DUtheta; 

  rowvec bhatt2 = est.row(est2.n_cols); 
  colvec pbhat2(z.n_rows); 

//  Rprintf("semi2 pt depm %d %d %d %d %d \n",semi2,pt,depmodel,CA1,CA2); 
//  est2.print("est2"); 
//  est.print("est"); 
//  z2.print("z2"); 
//  thetades.print("pp"); 
//  Rprintf(" pt depm %d %d \n",pt,depmodel); 

  // }}}

  for (s=0;s<Ntimes;s++) //   if (KMtimes[s]>0) 
  {
	  time=times(s); 
          rowvec bhatt = est.row(s); 
          vec pbhat = z * bhatt; 
	  if ((semi==1) & (cifmodel==1)) pbhat = pbhat + Zgamma*time;
	  if ((semi==1) & (cifmodel==2))  pbhat=pbhat%exp(Zgamma); 
//	  bhatt.print("bhatt");  pbhat.print("pbhatt"); 

	  if ((depmodel!=1) & (CA1!=CA2)) {
		  bhatt2 = est2.row(s); 
//		  bhatt2.print("est2"); 
		  pbhat2 = z2 * bhatt2; 
		  if ((semi2==1) & (cifmodel==1)) pbhat2 = pbhat2 + Z2gamma2*time;
	     if ((semi2==1) & (cifmodel==2)) pbhat2=pbhat2%exp(Z2gamma2); 
	  }

    for (j=0;j<antclust;j++) if (clustsize(j)>=2) { 
	  diff=0; sdj=0; 

          for (c=0;c<clustsize[j];c++) for (v=0;v<clustsize[j];v++) // {{{
	  if ((sym==1 && c!=v) || (sym==0 && c<v)) { 
	    i=clusterindex(j,c); k=clusterindex(j,v); 

     if ((entryage[i] < time) && (entryage[k]< time)) 
     if ((KMc[i] > 0) && (KMc[k] > 0)) {
	 ci=cause(i); ck=cause(k); 
	 resp1= ((y[k]<=time) && (ck==CA2));
	 resp2= ((y[i]<=time) && (ci==CA1))* ((y[k]<=time) && (ck==CA2));
         
	respst=((y[i]<=entryage[i]) && (ci==CA1))* ((y[k]<=time) && (ck==CA2)) + 
	       ((y[i]<=time) && (ci==CA1))* ((y[k]<=entryage[k]) && (ck==CA2)) ;
//	  Rprintf(" %d %d %d %d %d \n",j,i,k,ci,ck); 

         if (flexfunc==0) thetak= Xtheta(i); else thetak=Xtheta(s,i); 
	 pthetavec= thetades.row(i); 
//	 else { thetak=evalh(vtheta1,vtime,pthetavec,htheta,rhoR); }
         Li=pbhat(i); 
         Lk=pbhat(k); 
         if (CA1!=CA2) Lk=pbhat2(k); 

	 if (depmodel==1) ormarg=(1-exp(-Li))/exp(-Li); 
	 else if (depmodel==2) ormarg=(1-exp(-Li))*(1-exp(-Lk)); 
	 else if (depmodel==3) ormarg=(1-exp(-Li))*(1-exp(-Lk)); 

if (j<0) printf("mmmm  %lf %lf %lf %lf %lf %lf \n",pbhat(i),pbhat2(k),pbhat(k),Li,Lk,ormarg); 

	 int nocens= (ci!=0)+(ck!=0); 
	 nocens=min(nocens,2); 

	if (depmodel==1) { // cor model  // {{{

	 if (stabcens==0) { // {{{ responses 
	    if (samecens==1) { resp2=resp2/min(KMc[i],KMc[k]); 
		                respst=respst/min(KMc[i],KMc[k]); 
	    } 
	    else { resp2=resp2/(KMc[i]*KMc[k]); respst=respst/(KMc[i]*KMc[k]); 
	    }
	    resp1=resp1/KMtimes[k];
	 } else {
//	    cweight1= max(KMtimes[s],KMc[k]); 
	    cweight2= max(KMtimes[s],KMc[i]);  
	    if (samecens==1) {
		    resp2=resp2/min(cweight2,cweight2); 
		    respst=respst/min(cweight2,cweight2); 
	    }
	    resp1=resp1/KMtimes[k];
	 } // }}}

           if (trunkp[i]<1) {	
		   // DENNE skal tilpasse COR modellen
	      response=weight*(resp2- exp(thetak)*(ormarg+cif1entry[i]*cif2entry[k]
	    	   -(1-exp(-Li))*cif2entry[k]-(1-exp(-Lk))*cif1entry[i])/trunkp[i]);
	       diff=diff+response; 
	       sdj=sdj- weight*exp(thetak)*(ormarg+cif1entry[i]*cif2entry[k]- 
	 (1-exp(-Li))*cif2entry[k]- (1-exp(-Lk))*cif1entry[i])/trunkp[i];
	       resp3=-exp(thetak);
	   } else {
	     response=exp(thetak)*(exp(thetak)*ormarg*(resp1-resp2)-resp2); 
	     sdj=sdj+2*exp(2*thetak)*ormarg*(resp1-resp2)-exp(thetak)*resp2;
	     diff=diff+response; 
	     resp3=exp(2*thetak)*(resp1-resp2)*exp(Li);
	    } // }}}
	} else if (depmodel==2) { // RR model  // {{{
	 if (estimator==1) {
	    if (samecens==1) { resp2=resp2/min(KMc[i],KMc[k]); respst=respst/min(KMc[i],KMc[k]); } 
	    else { resp2=resp2/(KMc[i]*KMc[k]); respst=respst/(KMc[i]*KMc[k]); }
	    weight=1; 
	 } else if (estimator==0){ 
	    if (samecens==1) weight=1/min(KMc[i],KMc[k]); else weight=1/(KMc[i]*KMc[k]);
	 } else { weight=(time<KMc[i])*(time<KMc[k])*1; }

        if (trunkp[i]<1) {	
//           stpart=respst/trunkp[i]; 
	   response=weight*(resp2- exp(thetak)*(ormarg+cif1entry[i]*cif2entry[k]
		   -(1-exp(-Li))*cif2entry[k]-(1-exp(-Lk))*cif1entry[i])/trunkp[i]);
	   diff=diff+response; 
	   sdj=sdj- weight*exp(thetak)*(ormarg+cif1entry[i]*cif2entry[k]- 
	 (1-exp(-Li))*cif2entry[k]- (1-exp(-Lk))*cif1entry[i])/trunkp[i];
	   resp3=-exp(thetak);
	} else {
	   response=weight*(resp2-exp(thetak)*ormarg); 
	   diff=diff+response; 
	   sdj=sdj-exp(thetak)*ormarg*weight;
	   resp3=-exp(thetak);
	}

	} // }}}
	else if (depmodel==3) { // OR model  // {{{
	 if (estimator==1) {
	    if (samecens==1) { resp2=resp2/min(KMc[i],KMc[k]); 
		                 respst=respst/min(KMc[i],KMc[k]); 
	    } 
	    else { resp2=resp2/(KMc[i]*KMc[k]); respst=respst/(KMc[i]*KMc[k]); 
	    }
	    weight=1; 
	 } else if (estimator==0){ 
	    if (samecens==1) weight=1/min(KMc[i],KMc[k]); else weight=1/(KMc[i]*KMc[k]);
	 } else { weight=(time<KMc[i])*(time<KMc[k])*1; }

        if (trunkp[i]<1) {	
//           stpart=respst/trunkp[i]; 
           p11t=plack(exp(thetak),(1-exp(-Li)),(1-exp(-Lk)),dplack);
	   response=weight*dplack[0]*exp(thetak)*(resp2-p11t);
	   diff=diff+response; 
	   sdj=sdj-weight*dplack[0]*dplack[0]*exp(2*thetak); 
	   resp3=0;
	} else {
           p11t=plack(exp(thetak),(1-exp(-Li)),(1-exp(-Lk)),dplack);
	   response=weight*dplack[0]*exp(thetak)*(resp2-p11t);
	   diff=diff+response; 
	   sdj=sdj-2*weight*exp(2*thetak)*pow(dplack[0],2);
	   resp3=-dplack[0]*exp(thetak);
	}
	} // }}}

if (j<0) printf("uu %d %d %d %lf %lf %lf %lf %lf %lf \n",j,i,k,time,resp1,resp2,y[i],y[k],response); 
if (j<0) printf("uu2 %lf %lf %lf %lf %lf %lf %d %d \n",pbhat(i),pbhat(k),pbhat2(k),response,thetak,ormarg,ci,ck); 
}

        } /* for (c=0....... */   // }}}
	
//        if (flexfunc==0) vtheta2= pthetavec;  
//        else  evaldh(vtheta1,vtime,pthetavec,vtheta2,dhtheta,rhoR);

//       printf("%ld %ld %lf %lf\n",s,j,diff,weights[j]); 
       DUtheta =DUtheta+ sdj*weights[j]* (pthetavec * trans(pthetavec));
       vthetascore = (weights[j]*diff)*pthetavec; 
       Utheta=Utheta+vthetascore; 
//       Utheta.print("Ut"); 
       

       if (iid==1) for (c=0;c<pt;c++) thetiid(j,c)+=vthetascore(c); 
   } /* j in antclust */ 
 } // s < Ntimes

List res; 

res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) res["theta.iid"]=thetiid; 

return(res); 
} // }}}


