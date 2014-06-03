#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <complex.h>

using namespace arma;
using namespace Rcpp;

double claytonoakesP(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
{ // {{{
double valr=1,x,y,z;
double p00,p10,p01,p11; 

//double cifs=cif1+cif2; //double S=1+(cifs*(theta-1)); 
x=theta; y=cif1; z=cif2; 

valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);

p11=valr; 
p10=x-p11; 
p01=y-p11; 
p00=1-x-y+p11; 

if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 

if (status1==0 && status2==0) { // {{{
valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
} // }}}

if (status1==1 && status2==0) { // {{{
valr=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
dp(0)=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
} // }}}

if (status1==0 && status2==1) { // {{{
valr=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
dp(0)=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
} // }}}

if (status1==1 && status2==1) { // {{{
valr= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
dp(0)=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
} // }}}

return(valr); 
} // }}}
 
double placklikeP(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
{ // {{{
//double S,S2,a;
//S=1+cifs*(theta-1); S2=4*cif1*cif2*theta*(theta-1);
double x,y,z,valr=1,p11,p10,p01,p00; 
//double cifs=cif1+cif2; 
//a=(1+(theta-1)*(cifs)); 
x=theta; y=cif1; z=cif2; 

dp(0)=0; 

if (theta!=1) {
p11=(1+(y+z)*(x-1)-sqrt(pow(1+(y+z)*(x-1),2)-4*x*(x-1)*y*z))/(2*(x-1));
dp(0)= (y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (  -1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z  ),2)))/(2.*pow(-1 + x,2));
} else p11=cif1*cif2;

p11=p11;
p10=y-p11; 
p01=z-p11; 
p00=1-y-z+p11; 

if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 

return(valr); 
} // }}}

// double CclaytonoakesP(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
//{ // {{{
//double valr=1,x,y,z;
//double p00,p10,p01,p11; 
//
////double cifs=cif1+cif2; //double S=1+(cifs*(theta-1)); 
//x=theta; y=cif1; z=cif2; 
//
//valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
//dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
//
//p11=valr; 
//p10=x-p11; 
//p01=y-p11; 
//p00=1-x-y+p11; 
//
//double epsilon=1E-20; 
//cx_double Ctheta,Cvalr,Cy,Cz; 
//Ctheta=cx_double(theta,epsilon); 
//Cy=cx_double(y,0); 
//Cz=cx_double(z,0); 
//
//printf(" mig \n"); 
////Cvalr=  pow((1/pow(Cy,1/Ctheta) + 1/pow(Cz,1/Ctheta)) - 1,-Ctheta);
//double dd=imag(Cvalr)/epsilon; 
//
//printf("complex  %lf ",dd); 
//
//if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
//if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
//if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
//if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 
//
//if (status1==0 && status2==0) { // {{{
//valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
//dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
//} // }}}
//
//if (status1==1 && status2==0) { // {{{
//valr=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==0 && status2==1) { // {{{
//valr=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==1 && status2==1) { // {{{
//valr= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
//dp(0)=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
//} // }}}
//
//printf(" %lf \n",dp(0)); 
//return(valr); 
//} // }}}

cx_double Cpij(cx_double x, cx_double y, cx_double z,int status1, int status2)
{ // {{{ 
cx_double p11,one,two,four; 
one=cx_double(1,0); two=cx_double(2,0); four=cx_double(4,0); 
//p11=(1+(y+z)*(x-1)-sqrt(exp(ln(1+(y+z)*(x-1))*2)-4*x*(x-1)*y*z))/(2*(x-1));
p11=(one+(y+z)*(x-one)-sqrt(pow((one+(y+z)*(x-one)),two)-four*x*(x-one)*y*z));
//p11=one+(y+z)*(x-one); 
p11=p11/(two*(x-one));

// calculates probs depending on status 
if (status1==1 && status2==0) p11=y-p11;
if (status1==0 && status2==1) p11=z-p11; 
if (status1==0 && status2==0) p11=one-y-z+p11; 

return(p11); 
} // }}}

double CplacklikeP(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
{ // {{{
//double S,S2,a;
//S=1+cifs*(theta-1); S2=4*cif1*cif2*theta*(theta-1);
double x,y,z,valr=1,p11,p10,p01,p00; 
//double cifs=cif1+cif2; 
//a=(1+(theta-1)*(cifs)); 
x=theta; y=cif1; z=cif2; 

dp(0)=0; 

if (theta!=1) {
p11=(1+(y+z)*(x-1)-sqrt(pow(1+(y+z)*(x-1),2)-4*x*(x-1)*y*z))/(2*(x-1));
dp(0)= (y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (  -1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z  ),2)))/(2.*pow(-1 + x,2));
} else p11=cif1*cif2;

//det komplekse trick, derivative wrt y og zi, dvs D_1 P(y,z,theta) og D_2 P
cx_double CCp11,Ctheta,Cy,Cz; 
Ctheta=cx_double(theta,0); 
Cy=cx_double(y,1E-20); 
Cz=(cx_double) z; 
CCp11=Cpij(Ctheta,Cy,Cz,status1,status2); 
dp(1)=imag(CCp11)/1E-20; 
Cz=cx_double(z,1E-20); 
Cy=(cx_double) y; 
CCp11=Cpij(Ctheta,Cy,Cz,status1,status2); 
dp(2)=imag(CCp11)/1E-20; 
//printf(" %lf  ",imag(CCp11)/1E-20); 

p11=p11;
p10=y-p11; 
p01=z-p11; 
p00=1-y-z+p11; 

if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 
//printf(" %lf \n",dp(0)); 

return(valr); 
} // }}}
               
//double min(double a, double b) { if (a<b) return(a); else return(b); }
//double max(double a, double b) { if (a>b) return(a); else return(b); }

RcppExport SEXP twostageloglikebin( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, SEXP ithetades,
		SEXP icluster,SEXP iclustsize,SEXP iclusterindex, SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid 
) // {{{
{
// {{{ setting matrices and vectors, and exporting to armadillo matrices
 mat thetades = Rcpp::as<mat>(ithetades); 
 mat clusterindex = Rcpp::as<mat>(iclusterindex);
 colvec theta = Rcpp::as<colvec>(itheta);
 colvec clustsize = Rcpp::as<colvec>(iclustsize);
 int antclust = clusterindex.n_rows; 
 colvec cause = Rcpp::as<colvec>(icause);
 colvec pmargsurv = Rcpp::as<colvec>(ipmargsurv);
 colvec cluster = Rcpp::as<colvec>(icluster);
 colvec weights = Rcpp::as<colvec>(iweights);
// colvec entryage = Rcpp::as<colvec>(ientryage);
 colvec trunkp = Rcpp::as<colvec>(itrunkp);
 colvec secluster = Rcpp::as<colvec>(isecluster);

// array for derivative of flexible design
 NumericVector DXthetavec(iDXtheta);
 IntegerVector arrayDims(idimDX);
 arma::cube DXtheta(DXthetavec.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
 IntegerVector strata(istrata);

 int varlink= Rcpp::as<int>(ivarlink);
 int silent = Rcpp::as<int>(isilent);
 int depmodel= Rcpp::as<int>(idepmodel); 
 int iid= Rcpp::as<int>(iiid); 
 int antiid = Rcpp::as<int>(iantiid);

 mat Xtheta = Rcpp::as<mat>(iXtheta);

  int udtest=0; 
  if (udtest==1) { // {{{
//  Rprintf(" %d %d %d %d %d %d %d \n",samecens,inverse,semi,semi2,flexfunc,stabcens,silent); 
//  Rprintf(" %d %d %d %d %d %d %d \n",cifmodel,CA1,CA2,sym,depmodel,estimator,iid); 
//        est.print("est"); 
//	est2.print("est2"); 
//        z.print("z"); 
//	zsem.print("zsemi"); 
//	z2.print("z2"); 
        thetades.print("theta.des"); 
        clusterindex.print("clusterindex"); 
//        rvdes.print("rvdes"); 
	theta.print("theta"); 
	Xtheta.print("Xtheta"); 
//	  y.print("y-times"); 
	  clustsize.print("clustsize"); 
	  pmargsurv.print("margsurv"); 
	  cause.print("cause"); 
	  cluster.print("cluster"); 
//	  Zgamma.print("zgam"); 
//	  Z2gamma2.print("zgam2"); 
//	  KMtimes.print("KMtimes"); 
//	  KMc.print("KMc"); 
	  weights.print("weights"); 
//	  entryage.print("entryage"); 
//	  cif1entry.print("cif1entry"); 
//	  cif2entry.print("cif2entry"); 
	  trunkp.print("trunkp"); 
  } else if (udtest==2) 
  {
//  Rprintf(" %d %d %d %d %d %d %d \n",samecens,inverse,semi,semi2,flexfunc,stabcens,silent); 
//     Rprintf(" %d %d %d %d %d %d %d \n",cifmodel,CA1,CA2,sym,depmodel,estimator,iid); 
//      Rprintf("est %lf \n",mean(mean(est))); 
//      Rprintf("est2 %lf \n",mean(mean(est2))); 
//      Rprintf("z %lf \n",mean(mean(z))); 
//      Rprintf("zsem %lf \n",mean(mean(zsem))); 
//      Rprintf("z2 %lf \n",mean(mean(z2))); 
      mat mt=mean(thetades); 
      mt.print("meancol thetades"); 
//      Rprintf("theatdes %lf \n",mean(mean(thetades))); 
      Rprintf("ci %lf \n",mean(mean(clusterindex))); 
//      Rprintf("rvdes %lf \n",mean(mean(rvdes))); 
      Rprintf("theta %lf \n",mean(theta)); 
      Rprintf("Xtheta %lf \n",mean(mean(Xtheta))); 
//      Rprintf("y %lf \n",mean(y)); 
      Rprintf("ci %lf \n",mean(clustsize)); 
//      Rprintf("times %lf \n",mean(times)); 
      Rprintf("cause %lf \n",mean(cause)); 
      Rprintf("cluster %lf \n",mean(cluster)); 
//      Rprintf("Zgamma %lf \n",mean(Zgamma)); 
//      Rprintf("Z2gamma2 %lf \n",mean(Z2gamma2)); 
//      Rprintf("KMtimes %lf \n",mean(KMtimes)); 
//      Rprintf("KMc %lf \n",mean(KMc)); 
      Rprintf("weights %lf \n",mean(weights)); 
//      Rprintf("entry %lf \n",mean(entryage)); 
//      Rprintf("cif1entry %lf \n",mean(cif1entry)); 
//      Rprintf("cif2entry %lf \n",mean(cif2entry)); 
      Rprintf("trunkp %lf \n",mean(trunkp)); 
  } // }}}

  int ci,ck,i,j,c,s=0,k,v,c1,v1; 
  double ll=1,Li,Lk,diff=0;
  //double sdj=0;
  //  double Lit=1,Lkt=1,llt=1;
  double deppar=1,ssf=0,thetak=0; 
//  double plack(); 
  vec dplack(4); dplack.fill(0);
  vec dplackt(4); dplackt.fill(0);
  int pt=theta.n_rows; 
  vec ckij(4),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  if (iid==1) thetiid.fill(0); 

  colvec p11tvec(antclust); 
//  p11tvec=0; 
//  Rprintf(" %d \n",pt); 
  colvec Utheta(pt); 
  colvec vthetascore(pt); 
  colvec pthetavec(pt); 
  vec vtheta2(pt); 
  mat DUtheta(pt,pt); 
  DUtheta.fill(0); 
  Utheta.fill(0); 
//  if (!Utheta.is_finite()) {  Rprintf(" NA's i def U\n"); Utheta.print("U"); }
//  if (!DUtheta.is_finite()) { Rprintf(" NA's i def DU\n"); DUtheta.print("DU"); }
  // }}}

for (j=0;j<antclust;j++) if (clustsize(j)>=2) { 

    R_CheckUserInterrupt(); diff=0; //sdj=0; 

  for (c=0;c<clustsize(j)-1;c++) for (v=c+1;v<clustsize(j);v++) // {{{ //  if ((c<v)) 
  { 
     i=clusterindex(j,c); k=clusterindex(j,v); 
     if (strata(i)==strata(k)) { // {{{

     ci=cause(i); ck=cause(k); Li=pmargsurv(i); Lk=pmargsurv(k); 
         
     int flexfunc=0; 
      if (flexfunc==0) {
          thetak=Xtheta(i,0);  
	  pthetavec= trans(thetades.row(i)); 
	  vthetascore=1*pthetavec; 
      } else { // {{{ 
	  thetak=Xtheta(i,s); 
	  pthetavec = DXtheta(span(s),span(i),span::all); 
      } // }}} 

  if (depmodel==1){ if (varlink==1) deppar=1/exp(thetak); else deppar=1/thetak;}
  if (depmodel==2){ if (varlink==1) deppar=exp(thetak); else deppar=thetak; }

	if (depmodel==1) { // clayton-oakes  // {{{
	   ll=claytonoakesP(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*log(ll); 
	   diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	   if (varlink==1) diff=pow(deppar,1)*diff;  
	   if (varlink==0) diff=1*pow(deppar,2)*diff; 
	   //sdj=-pow(diff,2); 
	   // }}}
	} else if (depmodel==2) { // plackett model  // {{{
           ll=placklikeP(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*log(ll); 
	   // sdj=pow(dplack(0)/ll,2); 
	   if (varlink==1) diff=-deppar*dplack(0)/ll; 
	   if (varlink==0) diff=-dplack(0)/ll;
	} // }}}

     for (c1=0;c1<pt;c1++) 
     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)-=weights(i)*pow(diff,2)*vthetascore(c1)*vthetascore(v1);
     vthetascore=weights(i)*diff*vthetascore; 
     Utheta=Utheta+vthetascore; 

     if (iid==1) for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
     } // }}} strata(i)==strata(k) indenfor strata

  } /* for (c=0....... */   // }}}

} /* j in antclust */ 

//printf("Sum of squares %lf \n",ssf);theta.print("theta");Utheta.print("Utheta");DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) res["theta.iid"]=thetiid; 

return(res); 
} // }}}

