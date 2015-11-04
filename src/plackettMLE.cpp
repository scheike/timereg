#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace arma;
using namespace Rcpp;

//extern "C" double 
RcppExport SEXP claytonoakesR(SEXP itheta,SEXP  iistatus1,SEXP  iistatus2,SEXP  icif1,SEXP  icif2,SEXP iags) 
{ // {{{
 colvec theta = Rcpp::as<colvec>(itheta);
 colvec cif1 = Rcpp::as<colvec>(icif1);
 colvec cif2 = Rcpp::as<colvec>(icif2);
 colvec istatus1 = Rcpp::as<colvec>(iistatus1);
 colvec istatus2 = Rcpp::as<colvec>(iistatus2);

 colvec L=theta; 
 colvec dL=theta; 
 int n=cif1.size(); 
 double valr=1,dp=0;
 double x,y,z;
 int status1,status2; 

//  theta.print("theta"); 
//  istatus1.print("theta"); istatus2.print("theta"); cif1.print("cif1 "); cif2.print("cif2 "); 

  for (int i=0;i<n;i++)
  { // {{{
  x=theta(i); y=cif1(i); z=cif2(i); 
  status1=istatus1(i); status2=istatus2(i); 

// {{{ 

if (status1==0 && status2==0) { // {{{
valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
dp= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
} // }}}

if (status1==1 && status2==0) { // {{{
valr=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
dp=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
} // }}}

if (status1==0 && status2==1) { // {{{
valr=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
dp=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
} // }}}

if (status1==1 && status2==1) { // {{{
valr= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
dp=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
} // }}}

// }}} 

  L(i)=valr;
  dL(i)=dp; 
  } // }}} 

List res; 
res["like"]=L; 
res["dlike"]=dL; 

return(res);  
} // }}}
 
double claytonoakes(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
{ // {{{
  double valr=1,x,y,z;
  //double cifs=cif1+cif2; 
  //double S=1+(cifs*(theta-1)); 
  x=theta; y=cif1; z=cif2; 

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
 
double placklike(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
{ // {{{
  double valr=1,x,y,z; 
  // double cifs=cif1+cif2; 
  // double S=1+cifs*(theta-1); 
  // double S2=4*cif1*cif2*theta*(theta-1);
  // double a=(1+(theta-1)*(cifs)); 
  x=theta; y=cif1; z=cif2; 

dp(0)=0; 

if (status1==0 && status2==0) { // {{{
if (theta!=1) {
valr=(1+(y+z)*(x-1)-sqrt(pow(1+(y+z)*(x-1),2)-4*x*(x-1)*y*z))/(2*(x-1));
dp(0)= (y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (  -1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z  ),2)))/(2.*pow(-1 + x,2));
} else valr=cif1*cif2;
} // }}}

if (status1==1 && status2==0) { // {{{
if (theta!=1) {
valr= (-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x));
dp(0)= (1 + ((-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 +   x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*z - 4*x*z + 2*(-1   + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1   + x)) - (-1 + x - (-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (  -1 + x)*(y + z),2))))/(2.*pow(-1 + x,2));
} else valr=cif2;
} // }}}

if (status1==0 && status2==1) { // {{{
if (theta!=1) {
valr= (-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x)) ;
dp(0)= (1 + ((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z  ),2))))/(2.*(-1 + x)) - (-1 + x - (-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(  -1 + x,2));
} else valr=cif2;
} // }}}

if (status1==1 && status2==1) { // {{{
if (theta!=1) {
valr=(((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*(-1 + x));
dp(0)= ((-3*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(8.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),2.5)) + ((-4*(-1 + x)*z - 4*x*z + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((-4*(-1 + x)*y - 4*x*y + 2*(-1 + x)*(y + z) + 2*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + ((2*pow(-1 + x,2) - 4*(-1 + x)*x)*(-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) + (2*x)/sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2)))/(2.*(-1 + x)) - (((-4*(-1 + x)*x*y + 2*(-1 + x)*(1 + (-1 + x)*(y + z)))*(-4*(-1 + x)*x*z + 2*(-1 + x)*(1 + (-1 + x)*(y + z))))/(4.*pow(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2),1.5)) - (2*pow(-1 + x,2) - 4*(-1 + x)*x)/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z),2))))/(2.*pow(-1 + x,2));
} else  valr=1; 
} // }}}

return(valr); 
} // }}}
               
//double min(double a, double b) { if (a<b) return(a); else return(b); }
//double max(double a, double b) { if (a>b) return(a); else return(b); }

// {{{ laplace and derivatives
double ilapsf(double y, double x, double z)
{
return( exp((y*log(x)-log(z))/y)-x );
}
double lapsf(double y,double x, double z) 
{
return( pow(x,y)/pow((z + x),y));
}
vec Dlapsf(double y, double x, double z)
{ 
vec dL(3); 
dL(0) =(pow(z+x,y)*log(x)*pow(x,y) - pow(x,y)*log(x+z)*pow(z+x,y))/pow((z + x),(2*y)); 
dL(1) =(pow(z+x,y)*(y/x)*pow(x,y) - pow(x,y)*(y/(x+z))*pow(z+x,y))/pow((z+x),(2*y)); 
dL(2) =(- pow(x,y)*(y/(x+z))*pow(z+x,y))/pow(z+x,(2*y)); 
return(dL);
} 

vec D2lapsf(double y, double x, double z) 
{ 
vec dL(6); 
dL(0)= pow(x,y)* pow(x+z,(-y-1))* (y* log(x+z)-y* log(x)-1) ;
dL(1)= y* pow(x,(y-1))* pow(x+z,(-y-2))*(x-y* z) ;
dL(2)= y* (y+1)* pow(x,y)*pow((x+z),(-y-2));
dL(4)= pow(y,2)* (y+1)* pow(x,(y-1))* pow(x+z,(-y-2))+(-y-2)* y* (y+1)* pow(x,y)* pow(x+z,(-y-3));
dL(3)= y* pow(x,y)* pow(x+z,(-y-2))+(y+1)*pow(x,y)* 
	pow(x+z,(-y-2))+y* (y+1)* pow(x,y) *log(x)* 
	pow(x+z,(-y-2))-y *(y+1)* pow(x,y)* pow(x+z,(-y-2))* log(x+z);
dL(5)= y* (y+1)* (y+2)* (-pow(x,y))* pow(x+z,(-y-3));
return(dL); 
} 
vec Dilapsf(double y, double x, double z) 
{ 
vec dL(3); 
dL(0) = (x* pow(z,(-1/y))* log(z))/pow(y,2);
dL(1) = pow(z,(-1/y))-1;
dL(2) = -(x* pow(z,(-(y+1)/y)))/y;
return(dL);
} 
vec D2ilapsf(double y, double x, double z) 
{ 
vec dL(6); 
dL(0)=(x* pow(z,(-(y+1)/y))* (y-log(z)))/pow(y,3);
dL(1)= - pow(z,(-(y+1)/y))/y;
dL(2)= (x* (y+1)* pow(z,(-1/y-2)))/pow(y,2);
dL(4)= ((y+1)* pow(z,(-1/y-2)))/pow(y,2);
dL(3)= -(x* pow(z,(-1/y-2))* (y* (y+2)-(y+1)* log(z)))/pow(y,4);
dL(5)= -(x*(y+1)* (2* y+1)* pow(z,(-1/y-3)))/pow(y,3);
return(dL);
} 

// }}} 

double claytonoakesRVC(vec theta,mat thetades,mat ags, int status1,int status2,
		double cif1,double cif2,vec x1, vec x2, vec &dp) 
{ // {{{
  double valr=1;
  //double cifs=cif1+cif2; 
  //double S=1+(cifs*(theta-1)); 

// colvec theta = Rcpp::as<colvec>(itheta);
// mat thetades = Rcpp::as<mat>(ithetades);
// colvec x1= Rcpp::as<colvec>(irv1);
// colvec x2= Rcpp::as<colvec>(irv2);
// vec x1= irv1; vec x2= irv2;
// double cif1 = Rcpp::as<double>(icif1);
// double cif2 = Rcpp::as<double>(icif2);
// int status1 = Rcpp::as<int>(istatus1);
// int status2 = Rcpp::as<int>(istatus2);
// double cif1 = *icif1; double cif2 = *icif2;
// int status1 = *istatus1; int status2 = *istatus2;


 colvec dL=theta; dL.fill(0); 
 double f1,f2;
 //double cifs=cif1+cif2; 
 //double S=1+(cifs*(theta-1)); 
 f1=cif1; f2=cif2; 

 colvec par = thetades * theta; 

int nn=thetades.n_rows; 
int lpar=thetades.n_cols; 

 // {{{ first basic laplace derivatives
double lamtot1=  sum(x1 % par);
double ii1 = ilapsf(lamtot1,lamtot1,f1);
double ii2 = ilapsf(lamtot1,lamtot1,f2);

vec resv(nn); resv.fill(0); 
//resv <- rep(0,length(par))
//iresv <- rep(0,length(par))
vec iresv(nn); iresv.fill(0); 

double like=1,iisum; 
int i; 
for (i=0;i<nn;i++) 
{
iisum = x1(i)*ii1+x2(i)*ii2;
resv(i) = lapsf(par(i),lamtot1,iisum);
iresv(i) = iisum;
like=like*resv(i); 
}


vec D1(nn),D2(nn),D3(nn); 
vec D13(nn),D23(nn),D33(nn),D133(nn),D233(nn),D333(nn);


vec res(6),res0(6); 
for (i=0;i<nn;i++) 
{ // {{{ 
iisum   = x1(i)*ii1+x2(i)*ii2;
lamtot1=sum(trans(ags.row(i)) % theta); 
res     = D2lapsf(par(i),lamtot1,iisum);
res0    = Dlapsf( par(i),lamtot1,iisum);
  D1(i)   = res0(0);
  D2(i)   = res0(1);
  D3(i)   = res0(2);
 D13(i)  =  res(0);
 D23(i)  =  res(1);
 D33(i)  =  res(2);
D133(i) =  res(3);
D233(i) =  res(4);
D333(i) =  res(5);
} // }}} 

// {{{  derivatives for ilap 
vec restd = D2ilapsf(lamtot1,lamtot1,f1);
vec rest= Dilapsf(lamtot1,lamtot1,f1);
double Di1t  = rest(0);
double Di2t  = rest(1);
double Di3t  = rest(2);
double Di13t =  restd(0);
double Di23t =  restd(1);
//double Di33t =  restd(2);
//double Di133t = restd(3);
//double Di233t = restd(4);
//
vec ressd   = D2ilapsf(lamtot1,lamtot1,f2);
vec ress  = Dilapsf(lamtot1,lamtot1,f2);
double Di1s   = ress(0);
double Di2s   = ress(1);
double Di3s   = ress(2);
double Di13s  =  ressd(0);
double Di23s  =  ressd(1);
//double Di33s  =  ressd(2);
//double Di133s =  ressd(3);
//double Di233s =  ressd(4);
// }}}

// }}} 

//x1.print("x1"); thetades.print("td");

//rowvec tx1=x1.t(); tx1.print("x1");

vec msum(lpar); //msum.fill(1);
//vec msum= trans(trans(x1) * thetades); 
vec mdesi(lpar);
//msum.print("mdesi");
//mdesi.print("mdesi");

// {{{ derivatives of like ,ds dt, dtheta

//ddd1 <- ddd2 <- ddd3 <- matrix(0,npar,length(msum))
//dtheta <- dtt <- dts <- dt <- ds <- dsdt <- dtdtds <- 0
vec indtheta0(lpar),indtheta0t(lpar),indtheta0s(lpar),indtheta0ts(lpar);
vec dthetaj(lpar);
//,dtj(lpar),dsj(lpar);
double dtj,dsj;

vec dtheta(lpar),dtt(lpar),
    dts(lpar),
//    dt(lpar), ds(lpar),dsdt(lpar),
    dtdtds(lpar); 
double dt=0,ds=0,dsdt=0;


dtheta.fill(0); dtt.fill(0); dts.fill(0); 
//dt.fill(0); ds.fill(0); dsdt.fill(0); 
dtdtds.fill(0); 
vec dthetad3(lpar), 
    dilapthetad3t(lpar), 
    dilapthetad3s(lpar), 
dthetadtj(lpar), dthetadsj(lpar), dthetad33(lpar), 
dthetadtdsj(lpar), led3(lpar), led2(lpar), dilapthetad3(lpar); 


double dsdtj, numdsdtj ;

for (i=0;i<nn;i++) 
{
msum=trans(ags.row(i)); 
indtheta0  =msum*x1(i)*(Di1t+Di2t)+ msum*x2(i)*(Di1s+Di2s);
indtheta0t =msum*x1(i)*(Di13t+Di23t);
indtheta0s =msum*x2(i)*(Di13s+Di23s);
indtheta0ts= indtheta0t+indtheta0s;
mdesi=trans(thetades.row(i)); //mdesi.print("mdesi");
dthetaj = (mdesi*D1(i)+msum*D2(i)+ D3(i)*indtheta0);
dtheta =dtheta+dthetaj/resv(i);
//
dtj = D3(i)*x1(i)*Di3t;
dt  = dt+dtj/resv(i);
dsj = D3(i)*x2(i)*Di3s;
ds  = ds+dsj/resv(i);

// {{{  2nd deriv 
dsdtj    = D33(i)*x2(i)*x1(i)*Di3t*Di3s;
numdsdtj = (dsdtj*resv(i)-dsj*dtj);
dsdt     = dsdt+(dsdtj*resv(i)-dsj*dtj)/pow(resv(i),2);

dthetad3 =  x1(i)*(mdesi*D13(i)+
		   msum*D23(i)+ D33(i)*indtheta0); 
dilapthetad3t =  x1(i)*msum*(Di13t+Di23t);
dthetadtj = Di3t*dthetad3+D3(i)*dilapthetad3t;
dtt = dtt+(dthetadtj*resv(i)-dthetaj*dtj)/pow(resv(i),2);

dthetad3=x2(i)*(mdesi*D13(i)+ msum*D23(i)+ D33(i)*indtheta0);
dilapthetad3s = x2(i)*msum*(Di13s+Di23s);
dthetadsj = Di3s*dthetad3+D3(i)*dilapthetad3s;
dts = dts+(dthetadsj*resv(i)-dthetaj*dsj)/pow(resv(i),2);
// }}} 

// {{{ 3rd deriv 
dthetad33 = (mdesi*D133(i)+
              msum*D233(i)+ D333(i)*indtheta0);
dilapthetad3s = x2(i)*msum*(Di13s+Di23s);
dilapthetad3t = x1(i)*msum*(Di13t+Di23t);
dthetadtdsj = x1(i)*x2(i)*( D33(i)*Di3t*dilapthetad3s+ D33(i)*Di3s*dilapthetad3t+ Di3s*Di3t*dthetad33);
led3 = pow(resv(i),2)*( 
     dthetaj*dsdtj+resv(i)*dthetadtdsj-
     dthetadtj*dsj-dtj*dthetadsj);
led2 = numdsdtj*2*dthetaj*resv(i);

dtdtds = dtdtds+(led3-led2)/pow(resv[i],4); 

// }}} 
}

vec dttheta(lpar),dstheta(lpar),d3(lpar); 
dttheta = dtheta*dt*like + like*dtt;
dstheta = dtheta*ds*like + like*dts;
//
d3 = dts*dt*like+ds*dtt*like+ds*dt*dtheta*like+
      like*dtdtds+dtheta*dsdt*like;
//
//dsdt1 = ds*dt*like;
//dsdt2 =  like*dsdt;
dsdt =  ds*dt*like+like*dsdt;
dtheta =  dtheta *like;
dt = like*dt;
ds = like*ds;

// }}} 


if (status1==0 && status2==0) { // {{{
	valr=like; dp=dtheta;
} // }}}
if (status1==0 && status2==1) { // {{{
	valr=ds; 
	dp=dstheta ;
} // }}}
if (status1==1 && status2==0) { // {{{
	valr=dt; 
	dp=dttheta ;
} // }}}
if (status1==1 && status2==1) { // {{{
	valr=dsdt; 
	dp=d3; 
} // }}}

return(valr); 
} // }}}

double survivalRVC(vec theta,mat thetades,mat ags,int cause1,int cause2,vec cif1,vec cif2,mat x1, mat x2, vec &dp, vec &alllike) 
{ // {{{
  double lamtot1=1,valr=1;
  // index variable som angiver cause, men hvis cause==0 er index -1 
  int icause1,icause2;
  icause1=cause1; icause2=cause2; 
  if (cause1==0) icause1=1; 
  if (cause2==0) icause2=1; 

  int test=0,itest=0;
  if (itest==1) {
	  theta.print("theta"); 
	  thetades.print("theta-des"); 
	  printf(" %d %d \n",cause1,cause2); 
	  cif1.print("ci1"); 
	  cif2.print("ci1"); 
	  x1.print("x1"); 
	  x2.print("x2"); 
	  ags.print("ags"); 
  }


 colvec dL=theta; dL.fill(0); 
 vec f1=cif1; 
 vec f2=cif2; 

 if (test==1) { cif1.print("c1"); x1.print("x1"); }

 colvec par = thetades * theta; 

 if (test==1) { theta.print("theta"); thetades.print("t-des "); par.print("pp"); }

int nn=thetades.n_rows; 
int lpar=thetades.n_cols; 

 // {{{ first basic laplace derivatives
vec x11=trans(x1.row(0)); 

if (test==1) { x11.print("x11"); par.print("par "); }

vec resv(nn); resv.fill(0); 

if (test==1) { x1.print("x1"); cif1.print("cif1"); }

vec x1f1, x2f2; 
x1f1= trans(x1) * cif1; 
x2f2= trans(x2) * cif2; 

if (test==1) { x1f1.print("x1f1"); } 
//ags.print("ags"); theta.print("par"); 

double like=1,iisum; 
int i; 
for (i=0;i<nn;i++) 
{
lamtot1=sum(trans(ags.row(i)) % theta); 
iisum = x1f1(i)+x2f2(i);
resv(i) = lapsf(par(i),lamtot1,iisum);
like=like*resv(i); 
}

if (test==1) { resv.print("resv"); printf(" like %lf \n",like); }

vec D1(nn),D2(nn),D3(nn); 
vec D13(nn),D23(nn),D33(nn),D133(nn),D233(nn),D333(nn);

vec res(6),res0(6); 
for (i=0;i<nn;i++) 
{ // {{{ 
iisum = x1f1(i)+x2f2(i);
lamtot1=sum(trans(ags.row(i)) % theta); 
res0    = Dlapsf( par(i),lamtot1,iisum);
  D1(i)   = res0(0);
  D2(i)   = res0(1);
  D3(i)   = res0(2);
res     = D2lapsf(par(i),lamtot1,iisum);
 D13(i)  =  res(0);
 D23(i)  =  res(1);
 D33(i)  =  res(2);
D133(i) =  res(3);
D233(i) =  res(4);
D333(i) =  res(5);
} // }}} 

// }}} 

if (test==1) { x11.print("x11"); thetades.print("td"); }

vec msum(lpar); //msum.fill(1);
vec mdesi(lpar);

if (test==1) { msum.print("msum"); mdesi.print("mdesi"); }

// {{{ derivatives of like ,ds dt, dtheta

double dtj,dsj;
double dt=0,ds=0,dsdt=0;

vec dthetaj(lpar);
vec dtheta(lpar),dtt(lpar), dts(lpar),dtdtds(lpar); 
dtheta.fill(0); dtt.fill(0); dts.fill(0); dtdtds.fill(0); 
vec dthetad3(lpar), dthetadtj(lpar), dthetadsj(lpar), dthetad33(lpar), 
    dthetadtdsj(lpar), led3(lpar), led2(lpar), numdj(lpar); 

double dsdtj, numdsdtj ;

for (i=0;i<nn;i++) 
{ // {{{ 
mdesi=trans(thetades.row(i)); //mdesi.print("mdesi");
msum=trans(ags.row(i)); 
dthetaj = (mdesi*D1(i)+msum*D2(i));
dtheta =dtheta+dthetaj/resv(i);
//
dtj = D3(i)*x1(icause1-1,i);
dt  = dt+dtj/resv(i);
dsj = D3(i)*x2(icause2-1,i);
ds  = ds+dsj/resv(i);

// {{{  2nd deriv 
dsdtj    = D33(i)*x2(icause2-1,i)* x1(icause1-1,i);
numdsdtj = (dsdtj*resv(i)-dsj*dtj);
dsdt     = dsdt+numdsdtj/pow(resv(i),2);

dthetadtj =  x1(icause1-1,i)*(mdesi*D13(i)+ msum*D23(i)); 
dtt = dtt+(dthetadtj*resv(i)-dthetaj*dtj)/pow(resv(i),2);

dthetadsj=x2(icause2-1,i)*(mdesi*D13(i)+ msum*D23(i));
numdj=(dthetadsj*resv(i)-dthetaj*dsj); 
dts = dts+numdj/pow(resv(i),2);
// }}} 

// {{{ 3rd deriv 
dthetad33 = (mdesi*D133(i)+msum*D233(i))*x2(icause2-1,i)*x1(icause1-1,i);
//dthetad33.print("d33"); numdj.print("numdj"); 
led2 = numdj*2*dtj*resv(i);
//led2.print("numdj"); dthetadsj.print("dthetadsj"); dthetadtj.print("dthetadsj"); 
led3 = pow(resv(i),2)*(dthetad33*resv(i)+dthetadsj*dtj-dthetadtj*dsj-dthetaj*dsdtj); 
dtdtds = dtdtds+(led3-led2)/pow(resv(i),4); 
// }}} 
} // }}} 


vec dttheta(lpar),dstheta(lpar),d3(lpar); 
dttheta = dtheta*dt*like + like*dtt;
dstheta = dtheta*ds*like + like*dts;
//
d3 = dts*dt*like+ds*dtt*like+ds*dt*dtheta*like+like*dtdtds+dtheta*dsdt*like;
//
//dsdt1 = ds*dt*like;
//dsdt2 =  like*dsdt;
dsdt =  ds*dt*like+like*dsdt;
dtheta =  dtheta *like;
dt = like*dt;
ds = like*ds;


//printf(" %lf %lf \n",dt,ds); 

// }}} 


if (test==1) {
printf("32 her \n"); printf(" like %lf  \n",like); 
dtheta.print("dtheta"); dstheta.print("dstheta"); 
dttheta.print("dttheta"); d3.print("dsttheta"); 
}


alllike(0)=like; alllike(1)=-1*dt; alllike(2)=-1*ds; alllike(3)=dsdt; 
alllike(4)=cause1; alllike(5)=cause2; 
//printf("%d %d %lf %lf %lf %lf \n",cause1,cause2,like,dt,ds,dsdt); 
//alllike.print("all-2"); 

if (cause1==0 && cause2==0) { // {{{
	valr=like; dp=-1*dtheta; } // }}}
if (cause1==0 && cause2!=0) { // {{{ 
	valr=-1*ds; dp=dstheta ; } // }}}
if (cause1!=0 && cause2==0) { // {{{
	valr=-1*dt; dp=dttheta ; } // }}}
if (cause1!=0 && cause2!=0) { // {{{
	valr=dsdt; dp=-1*d3; } // }}}

return(valr); 
} // }}}

RcppExport SEXP survivalRV(SEXP itheta,SEXP istatus1,SEXP istatus2,
	   	     SEXP icif1,SEXP icif2,
                     SEXP irv1, SEXP irv2,SEXP ithetades,SEXP iags, SEXP ivarlink)
{ // {{{

 colvec theta = Rcpp::as<colvec>(itheta);
 mat thetades = Rcpp::as<mat>(ithetades);
 mat x1= Rcpp::as<mat>(irv1);
 mat x2= Rcpp::as<mat>(irv2);
 mat ags= Rcpp::as<mat>(iags);
// colvec x2= Rcpp::as<colvec>(irv2);
 vec cif1 = Rcpp::as<vec>(icif1);
 vec cif2 = Rcpp::as<vec>(icif2);
 int varlink = Rcpp::as<int>(ivarlink);
 int status1 = Rcpp::as<int>(istatus1);
 int status2 = Rcpp::as<int>(istatus2);

 int test=0; 
 if (test==1) {
	 theta.print("the"); 
	 thetades.print("the"); 
	 ags.print("ags"); 
	 x1.print("x1 "); 
	 x2.print("x2 "); 
	 cif1.print("cif"); 
	 cif2.print("cif"); 
 }

List ressl; 
ressl["par"]=theta; 

if (varlink==1) theta=exp(theta); 

 colvec dL=theta; dL.fill(0); 
// double f1,f2;
 //double cifs=cif1+cif2; 
 //double S=1+(cifs*(theta-1)); 
 vec f1=cif1; vec f2=cif2; 

 colvec par = thetades * theta; 

//int nn=thetades.n_rows; 
int lpar=thetades.n_cols; 

vec dp(lpar); dp.fill(0); 
vec all(4); 

double like=0; 
if (status1==0 && status2==0) { // {{{
like=survivalRVC(theta,thetades,ags,0,0,f1,f2,x1,x2,dp,all) ;
} // }}}
if (status1==0 && status2!=0) { // {{{
like=survivalRVC(theta,thetades,ags,0,status2,f1,f2,x1,x2,dp,all) ;
} // }}}
if (status1!=0 && status2==0) { // {{{
like=survivalRVC(theta,thetades,ags,status1,0,f1,f2,x1,x2,dp,all) ;
} // }}}
if (status1!=0 && status2!=0) { // {{{
like=survivalRVC(theta,thetades,ags,status1,status2,f1,f2,x1,x2,dp,all) ;
} // }}}

ressl["like"]=like; 
if (varlink==1) dp=dp % theta;  
ressl["dlike"]=dp;

ressl["theta"]=theta; 
ressl["par.des"]=thetades; 
ressl["varlink"]=varlink; 
ressl["alllike"]=all; 

return(ressl);  
} // }}}

RcppExport SEXP claytonoakesRV(SEXP itheta,SEXP istatus1,SEXP istatus2,SEXP icif1,SEXP icif2,
                     SEXP irv1, SEXP irv2,SEXP ithetades,SEXP iags, SEXP ivarlink)
{ // {{{
 colvec theta = Rcpp::as<colvec>(itheta);
 mat thetades = Rcpp::as<mat>(ithetades);
 mat ags= Rcpp::as<mat>(iags);
 colvec x1= Rcpp::as<colvec>(irv1);
 colvec x2= Rcpp::as<colvec>(irv2);
 double cif1 = Rcpp::as<double>(icif1);
 double cif2 = Rcpp::as<double>(icif2);
 int varlink = Rcpp::as<int>(ivarlink);
 int status1 = Rcpp::as<int>(istatus1);
 int status2 = Rcpp::as<int>(istatus2);

List ressl; 
ressl["par"]=theta; 

if (varlink==1) theta=exp(theta); 

 colvec dL=theta; dL.fill(0); 
 double f1,f2;
 //double cifs=cif1+cif2; 
 //double S=1+(cifs*(theta-1)); 
 f1=cif1; f2=cif2; 

 colvec par = thetades * theta; 

//int nn=thetades.n_rows; 
int lpar=thetades.n_cols; 

vec dp(lpar); dp.fill(0); 

if (status1==0 && status2==0) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,0,0,f1,f2,x1,x2,dp) ;
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
} // }}}
if (status1==0 && status2==1) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,0,1,f1,f2,x1,x2,dp) ;
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
} // }}}
if (status1==1 && status2==0) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,1,0,f1,f2,x1,x2,dp) ;
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
} // }}}
if (status1==1 && status2==1) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,1,1,f1,f2,x1,x2,dp) ;
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
} // }}}
ressl["theta"]=theta; 
ressl["par.des"]=thetades; 

vec obs(4); 
obs(0)=status1; obs(1)=f1; obs(2)=status2; obs(3)=f2; 
ressl["obs"]=obs; 
ressl["varlink"]=varlink; 

return(ressl);  
} // }}}

// for binary case, additive gamma from twostage survival
double claytonoakesbinRVC(vec theta,mat thetades,mat ags,int status1,int status2,double cif1,double cif2,vec x1, vec x2, vec &dp) 
{ // {{{
	double f1,f2,valr=1;
	colvec dL=theta; dL.fill(0); 
	f1=cif1; f2=cif2; 

	colvec par = thetades * theta; 
	int nn=thetades.n_rows; 
	int lpar=thetades.n_cols; 

	// {{{ first basic laplace derivatives
	double lamtot1= sum(x1 % par);
	double ii1 = ilapsf(lamtot1,lamtot1,f1);
	double ii2 = ilapsf(lamtot1,lamtot1,f2);

	vec resv(nn); resv.fill(0); 
	//resv <- rep(0,length(par))
	//iresv <- rep(0,length(par))
	vec iresv(nn); iresv.fill(0); 

	double like=1,iisum; 
	int i; 
	for (i=0;i<nn;i++) 
	{
//              lamtot1=sum(trans(ags.row(i)) % theta); 
//	        ii1 = ilapsf(lamtot1,lamtot1,f1);
//	        ii2 = ilapsf(lamtot1,lamtot1,f2);
		iisum = x1(i)*ii1+x2(i)*ii2;
		resv(i) = lapsf(par(i),lamtot1,iisum);
		iresv(i) = iisum;
		like=like*resv(i); 
	}

	vec D1(nn),D2(nn),D3(nn); 
	vec D13(nn),D23(nn),D33(nn),D133(nn),D233(nn),D333(nn);

	vec res(6),res0(6); 
	for (i=0;i<nn;i++) 
	{ // {{{ 

//                lamtot1=sum(trans(ags.row(i)) % theta); 
//	        ii1 = ilapsf(lamtot1,lamtot1,f1);
//	        ii2 = ilapsf(lamtot1,lamtot1,f2);
		iisum   = x1(i)*ii1+x2(i)*ii2;
		res     = D2lapsf(par(i),lamtot1,iisum);
		res0    = Dlapsf( par(i),lamtot1,iisum);
		D1(i)   = res0(0);
		D2(i)   = res0(1);
		D3(i)   = res0(2);
		D13(i)  =  res(0);
		D23(i)  =  res(1);
		D33(i)  =  res(2);
		D133(i) =  res(3);
		D233(i) =  res(4);
		D333(i) =  res(5);
	} // }}} 


	// {{{  derivatives for ilap 
	vec restd = D2ilapsf(lamtot1,lamtot1,f1);
	vec rest= Dilapsf(lamtot1,lamtot1,f1);
	double Di1t  = rest(0);
	double Di2t  = rest(1);
	double Di3t  = rest(2);
	double Di13t =  restd(0);
	double Di23t =  restd(1);
	//double Di33t =  restd(2);
	//double Di133t = restd(3);
	//double Di233t = restd(4);
	//
	vec ressd   = D2ilapsf(lamtot1,lamtot1,f2);
	vec ress  = Dilapsf(lamtot1,lamtot1,f2);
	double Di1s   = ress(0);
	double Di2s   = ress(1);
	double Di3s   = ress(2);
	double Di13s  =  ressd(0);
	double Di23s  =  ressd(1);
	//double Di33s  =  ressd(2);
	//double Di133s =  ressd(3);
	//double Di233s =  ressd(4);
	// }}}
	// }}} 

	//x1.print("x1"); thetades.print("td");
	//rowvec tx1=x1.t(); tx1.print("x1");

	vec msum= trans(trans(x1) * thetades); 
	vec mdesi(lpar);
//	vec msum(lpar);
	//msum.print("mdesi");
	//mdesi.print("mdesi");

	// {{{ derivatives of like ,ds dt, dtheta

	//ddd1 <- ddd2 <- ddd3 <- matrix(0,npar,length(msum))
	//dtheta <- dtt <- dts <- dt <- ds <- dsdt <- dtdtds <- 0
	vec indtheta0(lpar),indtheta0t(lpar),indtheta0s(lpar),indtheta0ts(lpar);
	vec dthetaj(lpar);
	//,dtj(lpar),dsj(lpar);
	double dtj,dsj;

	vec dtheta(lpar),dtt(lpar),
	    dts(lpar),
	    //    dt(lpar), ds(lpar),dsdt(lpar),
	    dtdtds(lpar); 
	double dt=0,ds=0,dsdt=0;


	dtheta.fill(0); dtt.fill(0); dts.fill(0); 
	//dt.fill(0); ds.fill(0); dsdt.fill(0); 
	dtdtds.fill(0); 
	vec dthetad3(lpar), 
	    dilapthetad3t(lpar), 
	    dilapthetad3s(lpar), 
	    dthetadtj(lpar), dthetadsj(lpar), dthetad33(lpar), 
	    dthetadtdsj(lpar), led3(lpar), led2(lpar), dilapthetad3(lpar); 


	//double dsdtj, numdsdtj ;

	for (i=0;i<nn;i++) 
	{
		mdesi=trans(thetades.row(i)); //mdesi.print("mdesi");
//                msum=trans(ags.row(i)); 
		indtheta0  =msum*x1(i)*(Di1t+Di2t)+ msum*x2(i)*(Di1s+Di2s);
		indtheta0t =msum*x1(i)*(Di13t+Di23t);
		indtheta0s =msum*x2(i)*(Di13s+Di23s);
		indtheta0ts= indtheta0t+indtheta0s;
		dthetaj = (mdesi*D1(i)+msum*D2(i)+ D3(i)*indtheta0);
		dtheta =dtheta+dthetaj/resv(i);
		//
		dtj = D3(i)*x1(i)*Di3t;
		dt  = dt+dtj/resv(i);
		dsj = D3(i)*x2(i)*Di3s;
		ds  = ds+dsj/resv(i);
	}

	vec dttheta(lpar),dstheta(lpar),d3(lpar); 
	dttheta = dtheta*dt*like + like*dtt;
	dstheta = dtheta*ds*like + like*dts;
	//
	d3 = dts*dt*like+ds*dtt*like+ds*dt*dtheta*like+
		like*dtdtds+dtheta*dsdt*like;
	//
	//dsdt1 = ds*dt*like;
	//dsdt2 =  like*dsdt;
	dsdt =  ds*dt*like+like*dsdt;
	dtheta =  dtheta *like;
	dt = like*dt;
	ds = like*ds;

	// }}} 

	double p11=like; 
	double p10=f1-p11; 
	double p01=f2-p11; 
	double p00=1-f1-f2+p11; 
	//printf("%lf %lf  %lf %lf %lf %lf \n",f1,f2,p11,p10,p10,p00); 
	//d3.print("d3"); 

	if (status1==0 && status2==0) { // {{{
		valr=p00; dp=dtheta;
	} // }}}
	if (status1==0 && status2==1) { // {{{
		valr=p01; 
		dp=-1*dtheta; 
	} // }}}
	if (status1==1 && status2==0) { // {{{
		valr=p10; 
		dp=-1*dtheta; 
	} // }}}
	if (status1==1 && status2==1) { // {{{
		valr=p11; 
		dp=dtheta; 
	} // }}}


	return(valr); 
} // }}}

// for binary case, R version 
RcppExport SEXP claytonoakesbinRV(SEXP itheta,SEXP istatus1,SEXP istatus2,SEXP icif1,SEXP icif2,
		SEXP irv1, SEXP irv2,SEXP ithetades,SEXP iags, SEXP ivarlink)
{ // {{{
	colvec theta = Rcpp::as<colvec>(itheta);
	mat thetades = Rcpp::as<mat>(ithetades);
	mat ags = Rcpp::as<mat>(iags);
	colvec x1= Rcpp::as<colvec>(irv1);
	colvec x2= Rcpp::as<colvec>(irv2);
	double cif1 = Rcpp::as<double>(icif1);
	double cif2 = Rcpp::as<double>(icif2);
	int varlink = Rcpp::as<int>(ivarlink);
	int status1 = Rcpp::as<int>(istatus1);
	int status2 = Rcpp::as<int>(istatus2);

	List ressl; 
	ressl["par"]=theta; 

	if (varlink==1) theta=exp(theta); 

	colvec dL=theta; dL.fill(0); 
	double f1,f2;
	//double cifs=cif1+cif2; 
	//double S=1+(cifs*(theta-1)); 
	f1=cif1; f2=cif2; 

	colvec par = thetades * theta; 

	//int nn=thetades.n_rows; 
	int lpar=thetades.n_cols; 

	vec dp(lpar); dp.fill(0); 

	if (status1==0 && status2==0) { // {{{
		double like=claytonoakesbinRVC(theta,thetades,ags,0,0,f1,f2,x1,x2,dp) ;
		ressl["like"]=like; 
		if (varlink==1) dp=dp % theta;  
		ressl["dlike"]=dp;
	} // }}}
	if (status1==0 && status2==1) { // {{{
		double like=claytonoakesbinRVC(theta,thetades,ags,0,1,f1,f2,x1,x2,dp) ;
		ressl["like"]=like; 
		if (varlink==1) dp=dp % theta;  
		ressl["dlike"]=dp;
	} // }}}
	if (status1==1 && status2==0) { // {{{
		double like=claytonoakesbinRVC(theta,thetades,ags,1,0,f1,f2,x1,x2,dp) ;
		ressl["like"]=like; 
		if (varlink==1) dp=dp % theta;  
		ressl["dlike"]=dp;
	} // }}}
	if (status1==1 && status2==1) { // {{{
		double like=claytonoakesbinRVC(theta,thetades,ags,1,1,f1,f2,x1,x2,dp) ;
		ressl["like"]=like; 
		if (varlink==1) dp=dp % theta;  
		ressl["dlike"]=dp;
	} // }}}
	ressl["theta"]=theta; 
	ressl["par.des"]=thetades; 

	vec obs(4); 
	obs(0)=status1; obs(1)=f1; obs(2)=status2; obs(3)=f2; 
	ressl["obs"]=obs; 
	ressl["varlink"]=varlink; 

	return(ressl);  
} // }}}

RcppExport SEXP twostageloglike( 
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
  double ll=1,Li,Lk,sdj=0,diff=0,loglikecont=0;
  double Lit=1,Lkt=1,llt=1,deppar=1,ssf=0,thetak=0; 
//  double plack(); 
  int pt=theta.n_rows; 
  vec dplack(pt); dplack.fill(pt);
  vec dplackt(pt); dplackt.fill(pt);
//  vec ckij(pt),dckij(pt),ckijvv(pt),dckijvv(pt),ckijtv(pt),dckijtv(pt),ckijvt(pt),dckijvt(pt);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  colvec loglikeiid(antiid); 
  if (iid==1) { thetiid.fill(0); 
	        loglikeiid.fill(0); 
  }

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

//  rowvec bhatt2 = est.row(est2.n_cols); 
//  colvec pbhat2(z.n_rows); 
// depmodel=5 
//  rvdes.print("rvdes"); 
//  thetades.print("ttt"); 
//  nr=rvdes.n_cols; 
//  vec alphaj(nr),alphai(nr),alpha(nr),
//      rvvec(nr),rvvec1(nr),rvvec2vv(nr),rvvec2vt(nr),rvvec2tv(nr);
//  vec  rvvec2(nr); 

  // }}}

for (j=0;j<antclust;j++) if (clustsize(j)>=2) { 

   R_CheckUserInterrupt(); diff=0; sdj=0; 

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
      } else { 
	  thetak=Xtheta(i,s); 
	  pthetavec = DXtheta(span(s),span(i),span::all); 
      }

  if (depmodel==1){ if (varlink==1) deppar=1/exp(thetak); else deppar=1/thetak;}
  if (depmodel==2){ if (varlink==1) deppar=exp(thetak); else deppar=thetak; }

	if (depmodel==1) { // clayton-oakes  // {{{

           if (trunkp(i)<1 || trunkp(k)<1) {	
		   Lit=trunkp(i); Lkt=trunkp(k); 
		   llt=claytonoakes(deppar,0,0,Lit,Lkt,dplackt);
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));
	           diff=dplack(0)/ll-dplackt(0)/llt; 
	   } else {
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
	           diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	   } 
	   if (varlink==1) diff=-pow(deppar,1)*diff;  
	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
	   sdj=pow(diff,2); 
	   // }}}
	} else if (depmodel==2) { // plackett model  // {{{
        if (trunkp(i)<1 || trunkp(k)<1) {	
           Lit=trunkp(i); Lkt=trunkp(k); 
           llt=placklike(deppar,0,0,Lit,Lkt,dplackt);
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*(log(ll)-log(llt));
	   loglikecont=(log(ll)-log(llt));
	   diff=dplack(0)/ll-dplackt(0)/llt; 
	   sdj=pow(diff,2); 
	} else {
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*log(ll); 
	   loglikecont=log(ll);
//	printf(" %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf \n",j,ci,ck,thetak,deppar,Li,Lk,weights(i),ll,log(ll),diff); 
//	printf(" %d %lf \n",j,ll); 
	}
	   if (varlink==1) diff=deppar*dplack(0)/ll; 
	   if (varlink==0) diff=dplack(0)/ll; 
	   sdj=pow(diff,2); 
	} // }}}

     for (c1=0;c1<pt;c1++) 
     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)+=weights(i)*sdj*vthetascore(c1)*vthetascore(v1);
     vthetascore=weights(i)*diff*vthetascore; 
     Utheta=Utheta+vthetascore; 

     if (iid==1) { for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
	           loglikeiid(secluster(i))+=loglikecont; 
     }
     } // }}} strata(i)==strata(k) indenfor strata

  } /* for (c=0....... */   // }}}

} /* j in antclust */ 

//printf("Sum of squares %lf \n",ssf); theta.print("theta"); Utheta.print("Utheta"); DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) { res["theta.iid"]=thetiid; 
	      res["loglikeiid"]=loglikeiid; 
            }

return(res); 
} // }}}

RcppExport SEXP twostageloglikeRV( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, SEXP ithetades,
		SEXP icluster,SEXP iclustsize,SEXP iclusterindex, SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid, SEXP irvdes, SEXP iags
) 
{ // {{{
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
 mat rvdes= Rcpp::as<mat>(irvdes); 
 mat ags = Rcpp::as<mat>(iags);

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
  double ll=1,Li,Lk,sdj=0,diff=0,loglikecont=0;
  double Lit=1,Lkt=1,llt=1,deppar=1,ssf=0,thetak=0; 
//  double plack(); 
 
  int pt=theta.n_rows; 
  vec dplack(pt); dplack.fill(0);
  vec dplackt(pt); dplackt.fill(0);
//  vec ckij(pt),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  colvec loglikeiid(antiid); 
  if (iid==1) { thetiid.fill(0); 
	        loglikeiid.fill(0); 
  }

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

//  rowvec bhatt2 = est.row(est2.n_cols); 
//  colvec pbhat2(z.n_rows); 
// depmodel=5 
//  rvdes.print("rvdes"); 
//  thetades.print("ttt"); 
    int nr=rvdes.n_cols; 
    vec rv2(nr),rv1(nr);
//  vec  rvvec2(nr); 
//  rv1.print("rv1"); 
//  rv2.print("rv2"); 

  vec etheta=theta; 
  // }}}

for (j=0;j<antclust;j++) if (clustsize(j)>=2) { 

   R_CheckUserInterrupt(); diff=0; sdj=0; 

  for (c=0;c<clustsize(j)-1;c++) for (v=c+1;v<clustsize(j);v++) // {{{ //  if ((c<v)) 
  { 
     i=clusterindex(j,c); k=clusterindex(j,v); 
//	  printf("cci 2 \n"); 
     if (strata(i)==strata(k)) { // {{{

     ci=cause(i); ck=cause(k); Li=pmargsurv(i); Lk=pmargsurv(k); 
         
     int flexfunc=0; 
      if (flexfunc==0) {
	  if (depmodel!=3) {
             thetak=Xtheta(i,0);  
	     pthetavec= trans(thetades.row(i)); 
	     vthetascore=1*pthetavec; 
	  }
      } else { 
	  thetak=Xtheta(i,s); 
	  pthetavec = DXtheta(span(s),span(i),span::all); 
      }

  if (depmodel==1){ if (varlink==1) deppar=1/exp(thetak); else deppar=1/thetak;}
  if (depmodel==2){ if (varlink==1) deppar=exp(thetak); else deppar=thetak; }
  if (depmodel==3){ if (varlink==1) etheta=exp(theta); else etheta=theta; }

	if (depmodel==1) { // clayton-oakes  // {{{

           if (trunkp(i)<1 || trunkp(k)<1) {	
		   Lit=trunkp(i); Lkt=trunkp(k); 
		   llt=claytonoakes(deppar,0,0,Lit,Lkt,dplackt);
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));
	           diff=dplack(0)/ll-dplackt(0)/llt; 
	   } else {
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
	           diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
//	printf("%lf %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",deppar,j,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	   } 
	   if (varlink==1) diff=-pow(deppar,1)*diff;  
	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
	   sdj=pow(diff,2); 
	   // }}}
	} else if (depmodel==3) { //  additive random gamma clayton-oakes  // {{{

	   rv1=trans(rvdes.row(i)); rv2=trans(rvdes.row(k)); 
//	   printf(" %d %d %d \n",j,i,k);
//         rv1.print("rv1");    rv2.print("rv2"); 
//	   thetades.print("thet"); 

           if (trunkp(i)<1 || trunkp(k)<1) { // {{{ 
		   Lit=trunkp(i); Lkt=trunkp(k); 
//		   llt=claytonoakesRV(theta,thetades,0,0,Lit,Lkt,rv1,rv2,dplackt);
//		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   llt=claytonoakesRVC(etheta,thetades,ags,0,0,Lit,Lkt,rv1,rv2,dplackt);
		   ll=claytonoakesRVC(etheta,thetades,ags,ci,ck,Li,Lk,rv1,rv2,dplack);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));

	           if (varlink==1) { 
			   dplackt=dplackt % etheta;  
			   dplack=dplack % etheta;  
		   }
	           vthetascore=dplack/ll-dplackt/llt; 
		   // }}}
	   } else {
		   ll=claytonoakesRVC(etheta,thetades,ags,ci,ck,Li,Lk,rv1,rv2,dplackt);
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);

	           if (varlink==1) dplackt=dplackt % etheta;  
	           vthetascore=dplackt/ll; 
	   } 
//	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
//	   sdj=pow(diff,2); 
	   // }}}
	} else if (depmodel==2) { // plackett model  // {{{
        if (trunkp(i)<1 || trunkp(k)<1) {	
           Lit=trunkp(i); Lkt=trunkp(k); 
           llt=placklike(deppar,0,0,Lit,Lkt,dplackt);
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*(log(ll)-log(llt));
	   loglikecont=(log(ll)-log(llt));
	   diff=dplack(0)/ll-dplackt(0)/llt; 
	   sdj=pow(diff,2); 
	} else {
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*log(ll); 
	   loglikecont=log(ll);
//	printf(" %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf \n",j,ci,ck,thetak,deppar,Li,Lk,weights(i),ll,log(ll),diff); 
//	printf(" %d %lf \n",j,ll); 
	}
	   if (varlink==1) diff=deppar*dplack(0)/ll; 
	   if (varlink==0) diff=dplack(0)/ll; 
	   sdj=pow(diff,2); 
	} // }}}

     if (depmodel!=3) {
	     for (c1=0;c1<pt;c1++) 
	     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)+=weights(i)*sdj*vthetascore(c1)*vthetascore(v1);
	     vthetascore=weights(i)*diff*vthetascore; 
	     Utheta=Utheta+vthetascore; 
	} else  { // additive gamma structure 
//		printf(" mig 1\n"); 
//		vthetascore.print("vvv"); 
	     for (c1=0;c1<pt;c1++) 
	     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)+=weights(i)*vthetascore(c1)*vthetascore(v1);
//		printf(" mig 1\n"); 
	     vthetascore=weights(i)*vthetascore; 
	     Utheta=Utheta+vthetascore; 
//		vthetascore.print("vvv 2"); 
	}

     if (iid==1) { for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
	           loglikeiid(secluster(i))+=loglikecont; 
     }


     } // }}} strata(i)==strata(k) indenfor strata

  } /* for (c=0....... */   // }}}

} /* j in antclust */ 

//printf("Sum of squares %lf \n",ssf); theta.print("theta"); Utheta.print("Utheta"); DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) { res["theta.iid"]=thetiid; 
	      res["loglikeiid"]=loglikeiid; 
            }

return(res); 
} // }}}


//This version of the program specifies which pairs that should be considered 
//in contrast to the general programs that considers all pairs within a cluster
//
//in addition and more importantly the random effects specifiation is within each pair 
//therefore all needed quantities comes for each pair : 
//theta.des and the random effects vectors, cluster and pair id's
//cluster og secluster now also follows the pairs 
RcppExport SEXP twostageloglikeRVpairs( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, 
		SEXP ithetades,SEXP icluster,SEXP iclustsize,SEXP iclusterindex, 
		SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid, SEXP irvdes,
		SEXP idimthetades, SEXP idimrvdes, SEXP inrvs, SEXP iags 
)  
{ // {{{ 
// {{{ 
//  setting matrices and vectors, and exporting to armadillo matrices
//  // {{{
// colvec nrvs = Rcpp::as<colvec>(inrvs);
 IntegerVector nrvs(inrvs);
 mat clusterindex = Rcpp::as<mat>(iclusterindex);
 colvec theta = Rcpp::as<colvec>(itheta);
int pt=theta.n_rows; 

 mat ags = Rcpp::as<mat>(iags);
 colvec clustsize = Rcpp::as<colvec>(iclustsize);
 // this is number of pairs (rather than clusters)
 int antclust= clusterindex.n_rows; 
 colvec cause = Rcpp::as<colvec>(icause);
 colvec pmargsurv = Rcpp::as<colvec>(ipmargsurv);
 colvec cluster = Rcpp::as<colvec>(icluster);
 colvec weights = Rcpp::as<colvec>(iweights);
// colvec entryage = Rcpp::as<colvec>(ientryage);
 colvec trunkp = Rcpp::as<colvec>(itrunkp);
 colvec secluster = Rcpp::as<colvec>(isecluster);
// mat rvdes= Rcpp::as<mat>(irvdes); 
 int depmodel= Rcpp::as<int>(idepmodel); 
 IntegerVector strata(istrata);

// array for derivative of flexible design
 NumericVector DXthetavec(iDXtheta);
 IntegerVector arrayDims(idimDX);
 arma::cube DXtheta(DXthetavec.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

//printf(" mig thetades 222222\n"); 
// array for parameter restrictions (one for each pair) pairs * (ant random effects)* (ant par)
 NumericVector thetadesvec(ithetades);
 IntegerVector arrayDims1(idimthetades); 
 IntegerVector arrayDD(3); 
//printf(" mig thetades 222222\n"); 

//printf(" %d %d %d \n",arrayDims1[0],arrayDims1[1],arrayDims1[2]); 

 arrayDD[2]=1; 
 if (depmodel==3)  {
	 arrayDD[0]=arrayDims1[0]; arrayDD[1]=arrayDims1[1]; arrayDD[2]=arrayDims1[2];  }
 else { arrayDD[0]=1; arrayDD[1]=1; arrayDD[2]=1;  }
arma::cube thetadesi(thetadesvec.begin(), arrayDD[0], arrayDD[1], arrayDD[2], false);
//printf(" mig cube thetades 222222\n"); 

// mat thetades(arrayDD[0],arrayDD[1]); 
// if  (depmodel!=3) 
// mat thetades=mat(arrayDims1[0],arrayDims1[1]*arrayDD[2],thetadesvec.begin()); 
mat thetades=mat(thetadesvec.begin(),arrayDims1[0],arrayDims1[1]*arrayDD[2],false); 

//printf(" mig rvdes \n"); 
// array for parameter restrictions (one for each pair) pairs * (ant random effects)* (ant par)
// array for pairwise random effects (two vectors for each pair)  pairs * 2* (ant random effects)
// mat rvdes= Rcpp::as<mat>(irvdes); 
 NumericVector rvdesvec(irvdes);
 IntegerVector arrayDims2(idimrvdes); 
 if (depmodel==3)  {
	 arrayDD[0]=arrayDims2[0]; arrayDD[1]=arrayDims2[1]; arrayDD[2]=arrayDims2[2];  }
 else { arrayDD[0]=1; arrayDD[1]=1; arrayDD[2]=1;  }
// printf(" %d %d %d \n",arrayDD[0], arrayDD[1], arrayDD[2]); 
arma::cube rvdesC(rvdesvec.begin(), arrayDD[0], arrayDD[1], arrayDD[2], false);

// printf(" her er lidt knas\n"); 
mat rvdes=mat(rvdesvec.begin(),arrayDims2[0],arrayDims2[1]*arrayDD[2],false); 
// if  (depmodel!=3) {
//	 mat rvdes = Rcpp::as<mat>(irvdes); 
// } else  mat rvdes(arrayDims2[1],arrayDims2[2]); 
// printf(" not !\n"); 


//  thetades.fill(0); 


 int varlink= Rcpp::as<int>(ivarlink);
 int silent = Rcpp::as<int>(isilent);
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
//        thetades.print("theta.des"); 
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
//      mat mt=mean(thetades); 
//      mt.print("meancol thetades"); 
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
  }   // }}}


  int ci,ck,i,j,s=0,k,c1,v1; 
  double ll=1,Li,Lk,sdj=0,diff=0,loglikecont=0;
  double Lit=1,Lkt=1,llt=1,deppar=1,ssf=0,thetak=0; 
//  double plack(); 
 
  vec dplack(pt); dplack.fill(0);
  vec dplackt(pt); dplackt.fill(0);
//  vec ckij(pt),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  colvec loglikeiid(antiid); 
  if (iid==1) { thetiid.fill(0); 
	        loglikeiid.fill(0); 
  }
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

//  rowvec bhatt2 = est.row(est2.n_cols); 
//  colvec pbhat2(z.n_rows); 
// depmodel=5 
//  rvdes.print("rvdes"); 
//  thetades.print("ttt"); 
//    int nr=rvdes.n_cols; 
//    printf("--------------------- %d \n",nr); 
    int nr=1; 
    if  (depmodel==3) nr=arrayDD[2]; 
    vec rv2(nr),rv1(nr);
//  vec  rvvec2(nr); 
//  rv1.print("rv1"); 
//  rv2.print("rv2"); 

  vec etheta=theta; 
  // 
  // // }}}
  
colvec likepairs(antclust); 

for (j=0;j<antclust;j++) { 

   R_CheckUserInterrupt(); diff=0; sdj=0; 

// index of subject's in pair "j"
   i=clusterindex(j,0); k=clusterindex(j,1); 
//	  printf("cci 2 \n"); 
     if (strata(i)==strata(k)) { // 

     // basic survival status 
     ci=cause(i); ck=cause(k); Li=pmargsurv(i); Lk=pmargsurv(k); 
//     printf(" %d %d %lf %lf \n",ci,ck,Li,Lk); 
         
     int flexfunc=0; 
      if (flexfunc==0) {
	  if (depmodel!=3) {
//  printf("pthetavec"); 
             thetak=Xtheta(i,0);  
	     pthetavec= trans(thetades.row(i)); 
	     vthetascore=1*pthetavec; 
//  printf("1 pthetavec \n"); 
//	     pthetavec= thetadesi.subcube(span(j),span(i),span::all); 
	  }
      } else { 
	  thetak=Xtheta(i,s); 
	  pthetavec = DXtheta(span(s),span(i),span::all); 
      }

  if (depmodel==1){ if (varlink==1) deppar=1/exp(thetak); else deppar=1/thetak;}
  if (depmodel==2){ if (varlink==1) deppar=exp(thetak); else deppar=thetak; }
  if (depmodel==3){ if (varlink==1) etheta=exp(theta); else etheta=theta; }

	if (depmodel==1) { // clayton-oakes  // {{{ 

           if (trunkp(i)<1 || trunkp(k)<1) {	
		   Lit=trunkp(i); Lkt=trunkp(k); 
		   llt=claytonoakes(deppar,0,0,Lit,Lkt,dplackt);
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));
	           diff=dplack(0)/ll-dplackt(0)/llt; 
	   } else {
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
	           diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf("%lf %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",deppar,j,(int) secluster(i),i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	  } 
	   if (varlink==1) diff=-pow(deppar,1)*diff;  
	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
	   sdj=pow(diff,2); 
	   // // }}} 
	} else if (depmodel==3) { //  additive random gamma clayton-oakes  // {{{ 

	// takes random effects specification for each pair
//	etheta.print("etheta"); 
	   
	// 3-dimensional array pairs*(2xrandom effects)
        int lnrv= nrvs(j)-1; // number of random effects for this cluster 	

	mat rv=rvdesC.slice(j); 
//	rv.print("rv"); 
        vec rv1= trans(rv.row(0)); 
        vec rv2= trans(rv.row(1)); 
//        mat rv2= rv.rows(nnn/2,nnn-1);
//        rv1= rvdesC.subcube( span(0),span(0,lnrv),span(j));
//	rv1.print("rv1"); 
//        rv2= rvdesC.subcube( span(1),span(0,lnrv),span(j));
//	rv2.print("rv2"); 

	// takes parameter relations for each pair
	// 3-dimensional array pairs*(random effects* pars )
//	mat thetades=thetadesi.subcube( span(j),span(0,lnrv),span::all);
//
//	mat thetadesv=thetadesi.subcube( span(j),span(0,lnrv),span(0,pt-1));
//	mat thetades=mat(thetadesv.begin(),nrvs(j),pt); 
	mat thetadesv=thetadesi.slice(j); 

	if (j< -10)  {
	   Rprintf(" %d %d \n",lnrv,pt); 
           rv1.print("rv1");    rv2.print("rv2"); 
	   thetadesv.print("thetades "); 
	   etheta.print("e-theta"); 
//	   mat test=mat(thetades.begin(),3,1); 
//	   test.print("test"); 
	}

           if (trunkp(i)<1 || trunkp(k)<1) { //  
		   Lit=trunkp(i); Lkt=trunkp(k); 
//		   llt=claytonoakesRV(theta,thetades,0,0,Lit,Lkt,rv1,rv2,dplackt);
//		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   llt=claytonoakesRVC(etheta,thetadesv,ags,0,0,Lit,Lkt,rv1,rv2,dplackt);
		   ll=claytonoakesRVC(etheta,thetadesv,ags,ci,ck,Li,Lk,rv1,rv2,dplack);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));

	           if (varlink==1) { 
			   dplackt=dplackt % etheta;  
			   dplack=dplack % etheta;  
		   }
	           vthetascore=dplack/ll-dplackt/llt; 
		   // 
	   } else {
		   ll=claytonoakesRVC(etheta,thetadesv,ags,ci,ck,Li,Lk,rv1,rv2,dplackt);
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
//		   printf("%lf %lf \n",weights(i),ll); 

	           if (varlink==1) dplackt=dplackt % etheta;  
	           vthetascore=dplackt/ll; 
	   } 
//	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
//	   sdj=pow(diff,2); 
	   // }}} 
	} else if (depmodel==2) { // plackett model  //  // {{{ 
        if (trunkp(i)<1 || trunkp(k)<1) {	
           Lit=trunkp(i); Lkt=trunkp(k); 
           llt=placklike(deppar,0,0,Lit,Lkt,dplackt);
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(j)*(log(ll)-log(llt));
	   loglikecont=(log(ll)-log(llt));
	   diff=dplack(0)/ll-dplackt(0)/llt; 
	   sdj=pow(diff,2); 
	} else {
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*log(ll); 
	   loglikecont=log(ll);
//	printf(" %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf \n",j,ci,ck,thetak,deppar,Li,Lk,weights(i),ll,log(ll),diff); 
	   if (varlink==1) diff=deppar*dplack(0)/ll; 
	   if (varlink==0) diff=dplack(0)/ll; 
	   sdj=pow(diff,2); 
	} // // }}} 
	} // }}}


        if (depmodel!=3) {
	     for (c1=0;c1<pt;c1++) 
	     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)+=weights(i)*sdj*vthetascore(c1)*vthetascore(v1);
	     vthetascore=weights(i)*diff*vthetascore; 
	     Utheta=Utheta+vthetascore; 
	} else  { // additive gamma structure 
//		vthetascore.print("vvv"); 
	     for (c1=0;c1<pt;c1++) 
	     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)+=weights(i)*vthetascore(c1)*vthetascore(v1);
	     vthetascore=weights(i)*vthetascore; 
	     Utheta=Utheta+vthetascore; 
//		vthetascore.print("vvv 2"); 
	}

     if (iid==1) { for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
	           loglikeiid((int) secluster(i))+=loglikecont; 
     }


     } //  strata(i)==strata(k) indenfor strata

if (iid==1)  likepairs(j)=ll; 

} /* j in antpairs */ 

//printf("Sum of squares %lf \n",ssf); theta.print("theta"); Utheta.print("Utheta"); DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) { res["theta.iid"]=thetiid; 
	      res["loglikeiid"]=loglikeiid; 
              res["likepairs"]=likepairs; 
            }

return(res); 
} // }}}  


// here consider (V^T Z) hazard(t,X)
// where hazard(t,X) is given as pmargsurv 
//
//This version of the program specifies which pairs that should be considered 
//
//in addition and more importantly the random effects specifiation is within each pair 
//therefore all needed quantities comes for each pair : 
//theta.des and the random effects vectors, cluster and pair id's
//cluster og secluster now also follows the pairs 
//
// computes also the needed weights for fitting the hazard models 
// given theta 
// That is the pairwise intensities are
// weight(Hist) * lambda(t,X) 
// where weight(Hist) is computed via additive hazards model
// for all jumps (conditional on what is known), these can
// be expressed via joint model 
// Lam(T_1) are the cumulative marginal hazards in T_1 and T_2
// D_1 S_(Lam(T_1),Lam(T_2))/S_0, D_2 S_(T_1,T_2)/S_0, 
// D_1 D_t S_(T_1,T_2), D_s D_2 S_(T_1,T_2), 
// D_1 D_s D_t S_(T_1,T_2), D_t D_s D_2 S_(T_1,T_2), 
RcppExport SEXP survivalloglikeRVpairs( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, 
		SEXP ithetades,SEXP icluster,SEXP iclustsize,SEXP iclusterindex, 
		SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid, SEXP irvdes,
		SEXP idimthetades, SEXP idimrvdes, SEXP inrvs ,
		SEXP iags
)  
{ // {{{ 
// {{{ 
//  setting matrices and vectors, and exporting to armadillo matrices
//  // {{{
// colvec nrvs = Rcpp::as<colvec>(inrvs);
 IntegerVector nrvs(inrvs);
 mat clusterindex = Rcpp::as<mat>(iclusterindex);
 colvec theta = Rcpp::as<colvec>(itheta);
 int pt=theta.n_rows; 
 colvec clustsize = Rcpp::as<colvec>(iclustsize);
 // this is number of pairs (rather than clusters)
 int antclust= clusterindex.n_rows; 
 colvec cause = Rcpp::as<colvec>(icause);
 mat pmargsurv = Rcpp::as<mat>(ipmargsurv);
 mat ags = Rcpp::as<mat>(iags);
 colvec cluster = Rcpp::as<colvec>(icluster);
 colvec weights = Rcpp::as<colvec>(iweights);
// colvec entryage = Rcpp::as<colvec>(ientryage);
 colvec trunkp = Rcpp::as<colvec>(itrunkp);
 colvec secluster = Rcpp::as<colvec>(isecluster);
// mat rvdes= Rcpp::as<mat>(irvdes); 
 int depmodel= Rcpp::as<int>(idepmodel); 
 IntegerVector strata(istrata);
 vec all(4); 

// array for derivative of flexible design
 NumericVector DXthetavec(iDXtheta);
 IntegerVector arrayDims(idimDX);
 arma::cube DXtheta(DXthetavec.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

//printf(" mig thetades 222222\n"); 
// array for parameter restrictions (one for each pair) pairs * (ant random effects)* (ant par)
 NumericVector thetadesvec(ithetades);
 IntegerVector arrayDims1(idimthetades); 
 IntegerVector arrayDD(3); 
//printf(" mig thetades 222222\n"); 

//printf(" %d %d %d \n",arrayDims1[0],arrayDims1[1],arrayDims1[2]); 

 if (depmodel==3)  {
	 arrayDD[0]=arrayDims1[0]; arrayDD[1]=arrayDims1[1]; arrayDD[2]=arrayDims1[2];  }
 else { arrayDD[0]=1; arrayDD[1]=1; arrayDD[2]=1;  }
arma::cube thetadesi(thetadesvec.begin(), arrayDD[0], arrayDD[1], arrayDD[2], false);
//printf(" mig cube thetades 222222\n"); 

// mat thetades(arrayDD[0],arrayDD[1]); 
// if  (depmodel!=3) 
// mat thetades=mat(arrayDims1[0],arrayDims1[1]*arrayDD[2],thetadesvec.begin()); 
mat thetades=mat(thetadesvec.begin(),arrayDims1[0],arrayDims1[1]*arrayDD[2],false); 

//printf(" %d %d %d \n",arrayDD[0], arrayDD[1], arrayDD[2]); 

//printf(" mig rvdes \n"); 
// array for parameter restrictions (one for each pair) pairs * (ant random effects)* (ant par)
// array for pairwise random effects (two vectors for each pair)  pairs * 2* (ant random effects)
// mat rvdes= Rcpp::as<mat>(irvdes); 
 NumericVector rvdesvec(irvdes);
 IntegerVector arrayDims2(idimrvdes); 
 if (depmodel==3)  {
	 arrayDD[0]=arrayDims2[0]; arrayDD[1]=arrayDims2[1]; arrayDD[2]=arrayDims2[2];  }
 else { arrayDD[0]=1; arrayDD[1]=1; arrayDD[2]=1;  }
// printf(" %d %d %d \n",arrayDD[0], arrayDD[1], arrayDD[2]); 
arma::cube rvdesC(rvdesvec.begin(), arrayDD[0], arrayDD[1], arrayDD[2], false);

mat rvdes=mat(rvdesvec.begin(),arrayDims2[0],arrayDims2[1]*arrayDD[2],false); 

 int varlink= Rcpp::as<int>(ivarlink);
 int silent = Rcpp::as<int>(isilent);
 int iid= Rcpp::as<int>(iiid); 
 int antiid = Rcpp::as<int>(iantiid);
 mat Xtheta = Rcpp::as<mat>(iXtheta);
 int ci,ck,i,j,s=0,k,c1,v1; 
 double ll=1,sdj=0,diff=0,loglikecont=0;
 double llt=1,ssf=0,thetak=0; 
 
 vec dplack(pt); dplack.fill(0);
 vec dplackt(pt); dplackt.fill(0);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  colvec loglikeiid(antiid); 
  if (iid==1) { thetiid.fill(0); 
	        loglikeiid.fill(0); 
  }
  colvec p11tvec(antclust); 
//
  colvec Utheta(pt); 
  colvec vthetascore(pt); 
  colvec pthetavec(pt); 
  vec vtheta2(pt); 
  mat DUtheta(pt,pt); 
  DUtheta.fill(0); 
  Utheta.fill(0); 
// printf("2 her er lidt knas\n"); 

  vec etheta=theta; 
  double deppar=1; 
  // 
  // // }}}
  
colvec likepairs(antclust); 
mat matlikepairs(antclust,6); 
vec allvec(6); 

//printf(" hej 1 \n"); 
//pmargsurv.print("pmarg"); 

for (j=0;j<antclust;j++) { // {{{  

   R_CheckUserInterrupt(); diff=0; sdj=0; 

// index of subject's in pair "j"
   i=clusterindex(j,0); k=clusterindex(j,1); 
//	  printf("cci 2 \n"); 
     if (strata(i)==strata(k)) { //  // {{{ 

     // basic survival status 
     ci=cause(i); ck=cause(k); 
     vec Li=trans(pmargsurv.row(i)); 
     vec Lk=trans(pmargsurv.row(k)); 

     if (j<-10) {
     printf("%d %d %d %d  \n",i,k,ci,ck); 
	     Li.print("Lk"); 
	     Lk.print("Lk"); 
     }
         
     int flexfunc=0; 
      if (flexfunc==0) {
	  if (depmodel!=3) {
             thetak=Xtheta(i,0);  
	     pthetavec= trans(thetades.row(i)); 
	     vthetascore=1*pthetavec; 
//	     pthetavec= thetadesi.subcube(span(j),span(i),span::all); 
	  }
      } else { 
	  thetak=Xtheta(i,s); 
	  pthetavec = DXtheta(span(s),span(i),span::all); 
      }

  if (depmodel==1){ if (varlink==1) deppar=1/exp(thetak); else deppar=1/thetak;}
  if (depmodel==2){ if (varlink==1) deppar=exp(thetak); else deppar=thetak; }
  if (depmodel==3){ if (varlink==1) etheta=exp(theta); else etheta=theta; }

	if (depmodel==1) { 
	} else if (depmodel==3) { //  additive random gamma clayton-oakes  // {{{ 

	// takes random effects specification for each pair
	// 3-dimensional array pairs*(2xrandom effects)
        int lnrv= nrvs(j)-1; // number of random effects for this cluster 	

	mat rv=rvdesC.slice(j);
	int nnn=rv.n_rows; 
	// first half of rows for person1 and second half for subject 2
	// nn number of competing risks
        mat rv1= rv.rows(0,nnn/2-1);
        mat rv2= rv.rows(nnn/2,nnn-1);

	if (j<-1 ) {
//	rv.print("rv"); 
	rv1.print("rv1"); 
	rv2.print("rv2"); 
	}

	// takes parameter relations for each pair
	// 3-dimensional array pairs*(random effects* pars )
	mat thetadesv=thetadesi.slice(j);

	if (j<-1)  {
	   Rprintf(" %d %d \n",lnrv,pt); 
           rv1.print("rv1");    
	   rv2.print("rv2"); 
	   thetadesv.print("thetades "); 
	   etheta.print("e-theta"); 
//	   mat test=mat(thetades.begin(),3,1); 
//	   test.print("test"); 
	}

           if (trunkp(i)<1 || trunkp(k)<1) { //  
                 vec Lit=trunkp.row(i); vec Lkt=trunkp.row(k); 
	   llt=survivalRVC(etheta,thetadesv,ags,0,0,Lit,Lkt,rv1,rv2,dplackt,allvec);
	   ll=survivalRVC(etheta,thetadesv,ags,(int) ci,(int) ck,Li,Lk,rv1,rv2,dplack,allvec);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));

	           if (varlink==1) { 
			   dplackt=dplackt % etheta;  
			   dplack=dplack % etheta;  
		   }
	           vthetascore=dplack/ll-dplackt/llt; 
		   // 
	   } else {
//		   printf(" %d %d %d %d \n",i,k,(int) ci,(int) ck); 
//printf(" hej 1 %d \n",nnn); ags.print("ags"); ags.print("thetadesv"); 
		   ll=survivalRVC(etheta,thetadesv,ags,(int) ci,(int) ck,Li,Lk,rv1,rv2,dplackt,allvec);
//printf(" hej 1 %d \n",nnn); 
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
		   if (j<-10)    {
			   printf("%lf %lf \n",weights(i),ll); 
			   allvec.print("allvec"); 
			   dplackt.print("dll"); 
		   }
	           if (varlink==1) dplackt=dplackt % etheta;  
	           vthetascore=dplackt/ll; 
	   } 
//	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
//	   sdj=pow(diff,2); 
	   // }}} 
	} else if (depmodel==2) { } 


        if (depmodel!=3) {
	     for (c1=0;c1<pt;c1++) 
	     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)+=weights(i)*sdj*vthetascore(c1)*vthetascore(v1);
	     vthetascore=weights(i)*diff*vthetascore; 
	     Utheta=Utheta+vthetascore; 
	} else  { // additive gamma structure 
//		vthetascore.print("vvv"); 
	     for (c1=0;c1<pt;c1++) 
	     for (v1=0;v1<pt;v1++) DUtheta(c1,v1)-=weights(i)*vthetascore(c1)*vthetascore(v1);
	     vthetascore=weights(i)*vthetascore; 
	     Utheta=Utheta+vthetascore; 
//		vthetascore.print("vvv 2"); 
	}

     if (iid==1) { for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
	           loglikeiid((int) secluster(i))+=loglikecont; 
     }


     } // // }}}  strata(i)==strata(k) indenfor strata

if (iid==1)  { likepairs(j)=ll; matlikepairs.row(j)=trans(allvec); }

} // }}} /* j in antpairs */ 

//printf("Sum of squares %lf \n",ssf); theta.print("theta"); Utheta.print("Utheta"); DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) { res["theta.iid"]=thetiid; 
	      res["loglikeiid"]=loglikeiid; 
              res["likepairs"]=likepairs; 
              res["all.likepairs"]=matlikepairs; 
            }


return(res); 
} // }}}  


