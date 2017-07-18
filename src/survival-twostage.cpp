#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace arma;
using namespace Rcpp;


double claytonoakes(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
{ // {{{
  double valr=1,x,y,z;
  //double cifs=cif1+cif2; 
  //double S=1+(cifs*(theta-1)); 
  // theta is 1/variance, which is how this function is called  
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

//extern "C" double 
RcppExport SEXP claytonoakesR(SEXP itheta,SEXP  iistatus1,SEXP  iistatus2,SEXP  icif1,SEXP  icif2,SEXP ivarlink) 
{ // {{{
 colvec theta = Rcpp::as<colvec>(itheta);
 colvec cif1 = Rcpp::as<colvec>(icif1);
 colvec cif2 = Rcpp::as<colvec>(icif2);
 colvec istatus1 = Rcpp::as<colvec>(iistatus1);
 colvec istatus2 = Rcpp::as<colvec>(iistatus2);
 int varlink = Rcpp::as<int>(ivarlink);

 // parametrization of Clayton-Oakes model 
 if (varlink==1) theta=1/exp(theta); else theta=1/theta;

 colvec L=theta; 
 colvec logL=theta; 
 colvec dlogL=theta; 
 int n=cif1.size(); 
 double valr=1;
 double x,y,z;
 int status1,status2; 
 colvec vdp(1); 

//  theta.print("theta"); 
//  istatus1.print("theta"); istatus2.print("theta"); cif1.print("cif1 "); cif2.print("cif2 "); 

  for (int i=0;i<n;i++)
  { // {{{
  x=theta(i); y=cif1(i); z=cif2(i); 
  status1=istatus1(i); status2=istatus2(i); 

  valr=claytonoakes(x,status1,status2,y,z,vdp); 
  L(i)=valr;
  logL(i)=log(valr);
  if (varlink==1) dlogL(i)=-pow(x,1)*vdp(0)/logL(i);  
  if (varlink==0) dlogL(i)=-1*pow(x,2)*vdp(0)/logL(i); 

  } // }}} 

List res; 
res["like"]=L; 
res["loglike"]=logL; 
res["dloglike"]=dlogL; 

return(res);  
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

//double lapsf(double y,double x, double z) 
//{
//return( exp(log(x)*y)/exp(log(z + x)*y));
//}
vec DlapsfOrig(double y, double x, double z)
{ 
vec dL(3); 
dL(0) =(pow(z+x,y)*log(x)*pow(x,y) - pow(x,y)*log(x+z)*pow(z+x,y))/pow((z + x),(2*y)); 
dL(1) =(pow(z+x,y)*(y/x)*pow(x,y) - pow(x,y)*(y/(x+z))*pow(z+x,y))/pow((z+x),(2*y)); 
dL(2) =(- pow(x,y)*(y/(x+z))*pow(z+x,y))/pow(z+x,(2*y)); 
return(dL);
} 

vec D2lapsfOrig(double y, double x, double z) 
{ 
vec dL(6); 
// D13
dL(0)= pow(x,y)* pow(x+z,(-y-1))* (y* log(x+z)-y* log(x)-1) ;
// D23
dL(1)= y* pow(x,(y-1))* pow(x+z,(-y-2))*(x-y* z) ;
// D33
dL(2)= y* (y+1)* pow(x,y)*pow((x+z),(-y-2));
// D133
dL(4)= 
pow(y,2)* (y+1)* pow(x,(y-1))* pow(x+z,(-y-2))+(-y-2)* y* (y+1)* pow(x,y)* pow(x+z,(-y-3));
// D233
dL(3)= y* pow(x,y)* pow(x+z,(-y-2))+(y+1)*pow(x,y)* 
	pow(x+z,(-y-2))+y* (y+1)* pow(x,y) *log(x)* 
	pow(x+z,(-y-2))-y *(y+1)* pow(x,y)* pow(x+z,(-y-2))* log(x+z);
// D333
dL(5)= y* (y+1)* (y+2)* (-pow(x,y))* pow(x+z,(-y-3));
return(dL); 
} 


vec Dlapsf(double y, double x, double z)
{ 
vec dL(3); 
double zxiy=pow(z+x,y),xiy=pow(x,y),zxi2y=pow((z + x),(2*y));
double ydxz=y/(x+z);
//printf("%lf %lf %lf %lf %lf \n",y,x,z,ydxz,zxi2y); 
double p2=xiy*zxiy; 
//if (zxi2y>0.000000000001) {
dL(0) =p2*(log(x) - log(x+z))/zxi2y; 
dL(1) =p2*((y/x) - ydxz)/zxi2y; 
dL(2) =p2*(-ydxz)/zxi2y; 
//}
return(dL);
} 

vec D2lapsf(double y, double x, double z) 
{ 
vec dL(6); 
double zximym2=pow(z+x,-y-2), zximym1=pow(z+x,-y-1), zximym3=pow(z+x,-y-3),
       xiym1=pow(x,y-1), xiy=pow(x,y), lxz=log(x+z), lx=log(x),yp1=y+1,
       ygyp1=y*(y+1);
       double p3=ygyp1*xiy; 
// D13
dL(0)= xiy* zximym1* (y* lxz-y* lx-1) ;
// D23
dL(1)= y* xiym1* zximym2*(x-y* z) ;
// D33
dL(2)= y* yp1*xiy*zximym2;
// D133
dL(4)= pow(y,2)* yp1* xiym1* zximym2+(-y-2)* p3* zximym3;
// D233
dL(3)= zximym2*( y* xiy+(y+1)*xiy + p3 *lx - p3 * lxz);
// D333
dL(5)= ygyp1 * (y+2)* (-xiy* zximym3);
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


cube vcrossmat(vec d, mat x1x2)
{ // {{{ 
cube dd(d.n_elem,x1x2.n_rows,2); 
	dd.slice(0)=d * trans(x1x2.col(0)); 
	dd.slice(1)=d * trans(x1x2.col(1)); 
return(dd); 
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
  if (itest==1) { // {{{
	  theta.print("theta"); 
	  thetades.print("theta-des"); 
	  Rprintf(" %d %d \n",cause1,cause2); 
	  cif1.print("ci1"); 
	  cif2.print("ci1"); 
	  x1.print("x1"); 
	  x2.print("x2"); 
	  ags.print("ags"); 
  } // }}} 

 vec dL=theta; 
 dL.fill(0); 
 vec par = thetades * theta; 

 if (test==1) { cif1.print("c1"); x1.print("x1"); }
 if (test==1) { theta.print("theta"); thetades.print("t-des "); par.print("pp"); }

 // nn number of parameters, ncr number of competing risks, lpar number of parameters pars=thetades * theta
int nn=thetades.n_rows,lpar=thetades.n_cols; 

// x1 = ncr x nrvs , Dcif1=ntheta x ncr  

vec sumtheta=ags * theta; 

// test=3; 
// wall_clock timer; 
// timer.tic(); 

 // {{{ first basic laplace derivatives
vec resv(nn); // resv.fill(0); 

//if (test==1) { x1.print("x1"); cif1.print("cif1"); }

// x1 = ncr x nrvs , cif1=ncr  
vec x1f1, x2f2; 
x1f1= trans(x1) * cif1; 
x2f2= trans(x2) * cif2; 

//if (test==1) { x1f1.print("x1f1"); } 
//ags.print("ags"); theta.print("par"); 

vec D1(nn),D2(nn),D3(nn); 
vec D13(nn),D23(nn),D33(nn),D133(nn),D233(nn),D333(nn);

vec res(6),res0(6); res.fill(0); res0.fill(0); 
//double x,y,z; 

double like=1,iisum; 
int i; 
for (i=0;i<nn;i++) 
{
lamtot1= sumtheta(i); 
iisum = x1f1(i)+x2f2(i);
resv(i) = lapsf(par(i),lamtot1,iisum);
//resv(i) = pow(lamtot1,par(i))/pow((iisum + lamtot1),par(i));
like=like*resv(i); 
res0    = Dlapsf( par(i),lamtot1,iisum);
  D1(i)   = res0(0); D2(i)   = res0(1); D3(i)   = res0(2);
res     = D2lapsf(par(i),lamtot1,iisum);
 D13(i)  =  res(0); D23(i)  =  res(1); D33(i)  =  res(2);
D133(i)  =  res(3); D233(i) =  res(4); D333(i) =  res(5);
}

// }}} 

//if (test==3) {
// double nt2 = timer.toc();
// printf("timer-loop 2  %lf antal rv's %d \n",100000*nt2,nn); 
//}

// timer.tic(); 
//if (test==1) { x11.print("x11"); thetades.print("td"); }

vec msum(lpar),mdesi(lpar);

if (test==1) { msum.print("msum"); mdesi.print("mdesi"); }

// {{{ derivatives of like ,ds dt, dtheta

double dtj,dsj,dsdtj,numdsdtj ;
double dt=0,ds=0,dsdt=0;

vec dthetaj(lpar),dthetaj1(lpar);
vec dtheta(lpar),dtt(lpar), dts(lpar),dtdtds(lpar); 
dtheta.fill(0); dtt.fill(0); dts.fill(0); dtdtds.fill(0); 
vec dthetad3(lpar), dthetadtj(lpar), dthetadsj(lpar), dthetad33(lpar), 
    dthetadtdsj(lpar), led3(lpar), led2(lpar), numdj(lpar),num1(lpar); 


//int ncr=cif1.n_elem; 
//mat dbasel(ncr,2),dbase(ncr,2); 
//dbase.fill(0); 
//mat dbasedt=dbase, dbaseds=dbase, dbasedtds=dbase; 
//vec db2(ncr,2),db1(ncr,2), db3(ncr,2),db4(ncr,2), db5(ncr,2);
//vec dbds2(ncr,2),dbds1(ncr,2), dbdt1(ncr,2),dbdt2(ncr,2);
//vec dbdtds2(ncr,2),dbdtds1(ncr,2); 

//mat dbasedtheta1(lpar,ncr); 
//mat dbasedtheta2(lpar,ncr); 
//
//cube dbasedtheta(lpar,ncr,2); 

//dbasedtheta.fill(0); 
//cube dbasedthetadt  =dbasedtheta, 
//     dbasedthetads  =dbasedtheta, 
//     dbasedthetadtds=dbasedtheta,
//     ddd            =dbasedtheta, 
//     ddd1           =dbasedtheta, 
//     ddd2           =dbasedtheta; 


// setting up random effect covariates for derivative wrt baselines for the two subjects
//cube x1x2(ncr,2,nn); 
//for (i=0;i<nn;i++) {
//       dbasel.col(0)= x1.col(i); dbasel.col(1)= x2.col(i); x1x2.slice(i)=dbasel; 
//}

//double x,y,z; // ,d1333,d2333,d3333; 

for (i=0;i<nn;i++) 
{ // {{{ 
mdesi=trans(thetades.row(i)); //mdesi.print("mdesi");
msum=trans(ags.row(i)); 
lamtot1= sumtheta(i); 
iisum = x1f1(i)+x2f2(i);
//
// // {{{  derivs Dtheta, Dt, Ds
dthetaj = (mdesi*D1(i)+msum*D2(i));
dtheta =dtheta+dthetaj/resv(i);
//
dtj = D3(i)*x1(icause1-1,i);
dt  = dt+dtj/resv(i);
dsj = D3(i)*x2(icause2-1,i);
ds  = ds+dsj/resv(i);
// }}} 

//  derivs with 3rd argument for derivative with respect profiled baseline
// {{{ 

//ddd3=D3(i)*x1x2.slice(i); 
//dbase= dbase+ D3(i)*x1x2.slice(i)/resv(i); 
//
//dbdt1   = D33(i)*x1(icause1-1,i)*x1x2.slice(i);
//dbasedt = dbasedt+(dbdt1*resv(i)-ddd3*dtj)/pow(resv(i),2);
//
//dbds1    = D33(i)*x2(icause2-1,i)*x1x2.slice(i);
//dbaseds  = dbaseds+(dbds1*resv(i)-ddd3*dsj)/pow(resv(i),2);
//
//dthetaj1 = resv(i)*(mdesi*D13(i)+msum*D23(i))+D3(i)*dthetaj;
//dbasedtheta=dbasedtheta+vcrossmat(dthetaj1,x1x2.slice(i))/pow(resv(i),2); 

// }}} 

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


// {{{  2nd deriv + Dbase 
//dbdtds1    = D333(i)*x2(icause2-1,i)*x1(icause1-1,i)*x1x2.slice(i);
//db4        = (dbdtds1*resv(i)+dsdtj*D3(i)*x1x2.slice(i)-dbds1*dtj-dsj*dbdt1);
//dbasedtds  = dbasedtds+(pow(resv(i),2)*db4-2*resv(i)*numdsdtj*D3(i)*x1x2.slice(i))/pow(resv(i),4);

// num * x1x2.slice(i) giver  num*x1 and num*x2 dvs cube
//ddd1 =  x1(icause1-1,i)*(mdesi*D133(i)+ msum*D233(i))*x1x2.slice(i); 
//ddd2=   ddd1*resv(i) + dthetadj * ddd3 -  (mdesi*D13(i)+msum*D23(i))*x1x2.slice(i)*dtj - dthetaj *dbdt1; 
//dtt = dtt+(dthetadtj*resv(i)-dthetaj*dtj)/pow(resv(i),2);
//ddd1=crossvmat((dthetadtj*resv(i)-dthetaj*dtj),2*resv(i)*D3(i)*x1x2.slice(i))
//ddd2=vcrossmat(
//ddd1= vcrossmat(dthetaj1,x1x2.slice(i)); 
//dbasedthetadt =dbasdthetadt+(ddd2*pow(resv(i),2)-2*resv(i)*ddd3* )/pow(resv(i),4);

// }}} 


// {{{ 3rd deriv 
dthetad33 = (mdesi*D133(i)+msum*D233(i))*x2(icause2-1,i)*x1(icause1-1,i);
led2 = numdj*2*dtj*resv(i);
led3 = pow(resv(i),2)*(dthetad33*resv(i)+dthetadsj*dtj-dthetadtj*dsj-dthetaj*dsdtj); 
dtdtds = dtdtds+(led3-led2)/pow(resv(i),4); 
// }}} 
//
//// {{{ 3rd deriv + Dbase
//
//x=lamtot1 ; y=par(i) ; z=iisum; 
//
//// d/dz y* pow(x,y)* pow(x+z,(-y-2))+(y+1)*pow(x,y)* pow(x+z,(-y-2))+y* (y+1)* pow(x,y) *log(x)* pow(x+z,(-y-2))-y *(y+1)* pow(x,y)* pow(x+z,(-y-2))* log(x+z);
//d2333= 
//pow(x,y)* (-2 - y)* y* pow(x + z,-3 - y) + pow(x,y)* (-2 - y)* (1 + y)* pow(x + z,-3 - y) - pow(x,y)* y* (1 + y)* pow(x + z,-3 - y) + 
//pow(x,y)* (-2 - y)* y* (1 + y)* pow(x + z,-3 - y)* log(x) - pow(x,y)* (-2 - y)* y* (1 + y)* pow(x + z,-3 - y)* log(x + z);
//// d/dz pow(y,2)* (y+1)* pow(x,(y-1))* pow(x+z,(-y-2))+(-y-2)* y* (y+1)* pow(x,y)* pow(x+z,(-y-3));
//d1333= pow(x,(-1 + y))* y* (1 + y)* (2 + y)*pow(x + z,-4 - y)*(3* x - y* z); 
//// d/dz y* (y+1)* (y+2)* (-pow(x,y))* pow(x+z,(-y-3));
//d3333=  pow(x,y)* y* (1+y)* (2+y)* (3+y)* pow(x+z,(-4-y)); 
//
////dthetad33 = (mdesi*D133(i)+msum*D233(i))*x2(icause2-1,i)*x1(icause1-1,i);
//////dthetad33.print("d33"); numdj.print("numdj"); 
////led2 = numdj*2*dtj*resv(i);
//////led2.print("numdj"); dthetadsj.print("dthetadsj"); dthetadtj.print("dthetadsj"); 
////led3 = pow(resv(i),2)*(dthetad33*resv(i)+dthetadsj*dtj-dthetadtj*dsj-dthetaj*dsdtj); 
////dtdtds = dtdtds+(led3-led2)/pow(resv(i),4); 
//// }}} 
//
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


//if (test==3) {
// double nt3 = timer.toc();
// printf("timer-loop 3 %lf \n",100000*nt3); 
//}

if (test==1) {
Rprintf("32 her \n"); Rprintf(" like %lf  \n",like); 
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

double survivalRVC2(vec theta,mat thetades,mat ags,int cause1,int cause2,vec cif1,vec cif2,mat x1, mat x2, vec &dp, vec &alllike) 
{ // {{{
  double lamtot1=1,x,y,z,like=1,iisum; 
  // index variable som angiver cause, men hvis cause==0 er index -1 
  int icause1=cause1,icause2=cause2,nn=thetades.n_rows,lpar=thetades.n_cols,i; 
  if (cause1==0) icause1=1; 
  if (cause2==0) icause2=1; 

  int itest=0;
  if (itest==1) { // {{{
	  theta.print("theta"); 
	  thetades.print("theta-des"); 
	  printf(" %d %d \n",cause1,cause2); 
	  cif1.print("ci1"); 
	  cif2.print("ci1"); 
	  x1.print("x1"); 
	  x2.print("x2"); 
	  ags.print("ags"); 
  } // }}} 

int type=0; 
if (cause1==0 && cause2==0) type=0; 
if (cause1!=0 && cause2==0) type=1; 
if (cause1==0 && cause2!=0) type=2; 
if (cause1!=0 && cause2!=0) type=3; 

 colvec dL=theta; dL.fill(0); 
 colvec par = thetades * theta; 

// if (test==1) { cif1.print("c1"); x1.print("x1"); }
// if (test==1) { theta.print("theta"); thetades.print("t-des "); par.print("pp"); }

vec sumtheta=ags * theta; 

// test=3; 
// wall_clock timer; 
// timer.tic(); 

 // {{{ first basic laplace derivatives
vec resv(nn); // resv.fill(0); 

vec msum(lpar),mdesi(lpar);
double dtj=0,dsj=0,dsdtj=0,numdsdtj=0;
double dt=0,ds=0,dsdt=0;

vec dthetaj(lpar);
vec dtheta(lpar),dtt(lpar), dts(lpar),dtdtds(lpar); 
dtheta.fill(0); dtt.fill(0); dts.fill(0); dtdtds.fill(0); 
vec dthetad3(lpar), dthetadtj(lpar), dthetadsj(lpar), dthetad33(lpar), 
    dthetadtdsj(lpar), led3(lpar), led2(lpar), numdj(lpar); 


//if (test==1) { x1.print("x1"); cif1.print("cif1"); }

vec x1f1, x2f2; 
x1f1= trans(x1) * cif1; 
x2f2= trans(x2) * cif2; 

//printf(" hej 1  \n"); 
//if (test==1) printf(" hej\n"); 
//if (test==1) { x1f1.print("x1f1"); } 
//ags.print("ags"); theta.print("par"); 

vec D1(nn),D2(nn),D3(nn); 
vec D13(nn),D23(nn),D33(nn),D133(nn),D233(nn),D333(nn);
vec res(6),res0(6); res.fill(0); res0.fill(0); 

for (i=0;i<nn;i++) 
{ // {{{ 
lamtot1= sumtheta(i); 
iisum = x1f1(i)+x2f2(i);
resv(i) = lapsf(par(i),lamtot1,iisum);
//resv(i) = pow(lamtot1,par(i))/pow((iisum + lamtot1),par(i));
like=like*resv(i); 
res0    = Dlapsf( par(i),lamtot1,iisum);
  D1(i)   = res0(0); D2(i)   = res0(1); D3(i)   = res0(2);
  y=par(i); x=lamtot1; z=iisum; 
//res     = D2lapsf(par(i),lamtot1,iisum);
 if (type!=0) {
 D13(i)  =   pow(x,y)* pow(x+z,(-y-1))* (y* log(x+z)-y* log(x)-1) ;
 D23(i)  =   y* pow(x,(y-1))* pow(x+z,(-y-2))*(x-y* z) ;
 }

if (type==3)  {
D33(i)  =  y* (y+1)* pow(x,y)*pow((x+z),(-y-2));
D133(i)  =  y* pow(x,y)* pow(x+z,(-y-2))+(y+1)*pow(x,y)* 
            pow(x+z,(-y-2))+y* (y+1)* pow(x,y) *log(x)* 
            pow(x+z,(-y-2))-y *(y+1)* pow(x,y)* pow(x+z,(-y-2))* log(x+z);
D233(i) =  pow(y,2)* (y+1)* pow(x,(y-1))* pow(x+z,(-y-2))+(-y-2)* y* (y+1)* pow(x,y)* pow(x+z,(-y-3));
D333(i) = y* (y+1)* (y+2)* (-pow(x,y))* pow(x+z,(-y-3));
}

//printf(" hej  \n"); 


//if (test==3) {
// double nt2 = timer.toc();
// printf("timer-loop 2  %lf antal rv's %d \n",100000*nt2,nn); 
//}

// timer.tic(); 
//if (test==1) { x11.print("x11"); thetades.print("td"); }

//if (test==1) printf(" hej\n"); 


//if (test==1) { msum.print("msum"); mdesi.print("mdesi"); }

// {{{ derivatives of like ,ds dt, dtheta

//mat ddcif1=Dcif, ddcif2=Dcif, ddcif3=Dcif;
//Dcif.fill(0); 
//ddcif1.fill(0); ddcif2.fill(0); 
//
//vec dciflike(Dcif.n_rows), 
//    dcifdt(Dcif.n_rows), dcifds(Dcif.n_rows), 
//    d3like1(Dcif.n_rows), 
//    d3like2(Dcif.n_rows), 
//    d3like3(Dcif.n_rows), 
//    d3like(Dcif.n_rows), dcifdtds(Dcif.n_rows);
//dciflike.fill(0); dcifdt.fill(0); dcifds.fill(0); dcifdtds.fill(0); 

mdesi=trans(thetades.row(i)); //mdesi.print("mdesi");
msum=trans(ags.row(i)); 
dthetaj = (mdesi*D1(i)+msum*D2(i));
dtheta =  dtheta+dthetaj/resv(i);
//
//printf(" hej 2  \n"); 

if (type==1 || type==3) {
dtj = D3(i)*x1(icause1-1,i);
dt  = dt+dtj/resv(i);
}
if (type==2 || type==3) {
dsj = D3(i)*x2(icause2-1,i);
ds  = ds+dsj/resv(i);
}

//d3like1=D3(i)*x1.col(i); 
//dciflike=dciflike+d3like1/resv(i); 
//
//d3like2=D33(i)*x1(icause1-1,i)*x1.col(i); 
//dcifdt=dcifdt+(resv(i)*d3like2-dtj*d3like)/pow(resv(i),2);
//d3like3=D33(i)*x2(icause2-1,i)*x2.col(i); 

//printf(" hej 3  \n"); 

// {{{  2nd deriv 

if (type==1 || type==3) {
dthetadtj =  x1(icause1-1,i)*(mdesi*D13(i)+ msum*D23(i)); 
dtt = dtt+(dthetadtj*resv(i)-dthetaj*dtj)/pow(resv(i),2);
}

if (type==2 || type==3) {
dthetadsj=x2(icause2-1,i)*(mdesi*D13(i)+ msum*D23(i));
numdj=(dthetadsj*resv(i)-dthetaj*dsj); 
dts = dts+numdj/pow(resv(i),2);
}

if (type==3) {
dsdtj    = D33(i)*x2(icause2-1,i)* x1(icause1-1,i);
numdsdtj = (dsdtj*resv(i)-dsj*dtj);
dsdt     = dsdt+numdsdtj/pow(resv(i),2);
}
// }}} 

//printf(" hej 4  \n"); 

// {{{ 3rd deriv 
if (type==3) {
dthetad33 = (mdesi*D133(i)+msum*D233(i))*x2(icause2-1,i)*x1(icause1-1,i);
//dthetad33.print("d33"); numdj.print("numdj"); 
led2 = numdj*2*dtj*resv(i);
//led2.print("numdj"); dthetadsj.print("dthetadsj"); dthetadtj.print("dthetadsj"); 
led3 = pow(resv(i),2)*(dthetad33*resv(i)+dthetadsj*dtj-dthetadtj*dsj-dthetaj*dsdtj); 
dtdtds = dtdtds+(led3-led2)/pow(resv(i),4); 
}
// }}} 

} // }}} 

//printf(" hej 5  \n"); 

vec dttheta(lpar),dstheta(lpar),d3(lpar); 
if (type==1 ) dttheta = dtheta*dt*like + like*dtt;
if (type==2 ) dstheta = dtheta*ds*like + like*dts;
//
if (type==3)  {
d3 = dts*dt*like+ds*dtt*like+ds*dt*dtheta*like+like*dtdtds+dtheta*dsdt*like;
dsdt =  ds*dt*like+like*dsdt;
}
dtheta =  dtheta *like;
dt = like*dt;
ds = like*ds;

//printf(" hej 6  \n"); 
//printf(" %lf %lf \n",dt,ds); 
// }}} 

//if (test==3) {
// double nt3 = timer.toc();
// printf("timer-loop 3 %lf \n",100000*nt3); 
//}

//if (test==1) {
//printf("32 her \n"); printf(" like %lf  \n",like); 
//dtheta.print("dtheta"); dstheta.print("dstheta"); 
//dttheta.print("dttheta"); d3.print("dsttheta"); 
//}


alllike(0)=like; 
if (type==1) alllike(1)=-1*dt; 
if (type==2) alllike(2)=-1*ds; 
if (type==3) alllike(3)=dsdt; 
alllike(4)=cause1; 
alllike(5)=cause2; 
//printf("%d %d %lf %lf %lf %lf \n",cause1,cause2,like,dt,ds,dsdt); 
//alllike.print("all-2"); 


double valr=1; 
if (type==0) {  valr=like;  dp=-1*dtheta; } 
if (type==1) {  valr=-1*dt; dp=dttheta ; } 
if (type==2) {  valr=-1*ds; dp=dstheta ; } 
if (type==3) {  valr=dsdt;  dp=-1*d3; } 


return(valr); 
} // }}}
// }}}

double survivalRVC2all(vec theta,mat thetades,mat ags,int cause1,int cause2,vec cif1,vec cif2,mat x1, mat x2, vec &dp, vec &alllike) 
{ // {{{
  double lamtot1=1,x,y,z,like=1,iisum; 
  // index variable som angiver cause, men hvis cause==0 er index -1 
  int icause1=cause1,icause2=cause2,nn=thetades.n_rows,lpar=thetades.n_cols,i; 
  if (cause1==0) icause1=1; 
  if (cause2==0) icause2=1; 

//  int test=0,itest=0;
//  if (itest==1) { // {{{
//	  theta.print("theta"); 
//	  thetades.print("theta-des"); 
//	  printf(" %d %d \n",cause1,cause2); 
//	  cif1.print("ci1"); 
//	  cif2.print("ci1"); 
//	  x1.print("x1"); 
//	  x2.print("x2"); 
//	  ags.print("ags"); 
//  } // }}} 

int type=0; 
if (cause1==0 && cause2==0) type=0; 
if (cause1!=0 && cause2==0) type=1; 
if (cause1==0 && cause2!=0) type=2; 
if (cause1!=0 && cause2!=0) type=3; 

 colvec dL=theta; dL.fill(0); 
 colvec par = thetades * theta; 

// if (test==1) { cif1.print("c1"); x1.print("x1"); }
// if (test==1) { theta.print("theta"); thetades.print("t-des "); par.print("pp"); }

vec sumtheta=ags * theta; 

// test=3; 
// wall_clock timer; 
// timer.tic(); 

 // {{{ first basic laplace derivatives
vec resv(nn); // resv.fill(0); 

vec msum(lpar),mdesi(lpar);
double dtj=0,dsj=0,dsdtj=0,numdsdtj=0;
double dt=0,ds=0,dsdt=0;

vec dthetaj(lpar);
vec dtheta(lpar),dtt(lpar), dts(lpar),dtdtds(lpar); 
dtheta.fill(0); dtt.fill(0); dts.fill(0); dtdtds.fill(0); 
vec dthetad3(lpar), dthetadtj(lpar), dthetadsj(lpar), dthetad33(lpar), 
    dthetadtdsj(lpar), led3(lpar), led2(lpar), numdj(lpar); 


//if (test==1) { x1.print("x1"); cif1.print("cif1"); }

vec x1f1, x2f2; 
x1f1= trans(x1) * cif1; 
x2f2= trans(x2) * cif2; 

//if (test==1) printf(" hej\n"); 
//if (test==1) { x1f1.print("x1f1"); } 
//ags.print("ags"); theta.print("par"); 

vec D1(nn),D2(nn),D3(nn); 
vec D13(nn),D23(nn),D33(nn),D133(nn),D233(nn),D333(nn);
vec res(6),res0(6); res.fill(0); res0.fill(0); 

for (i=0;i<nn;i++) 
{ // {{{ 
lamtot1= sumtheta(i); 
iisum = x1f1(i)+x2f2(i);
resv(i) = lapsf(par(i),lamtot1,iisum);
//resv(i) = pow(lamtot1,par(i))/pow((iisum + lamtot1),par(i));
like=like*resv(i); 
res0    = Dlapsf( par(i),lamtot1,iisum);
  D1(i)   = res0(0); D2(i)   = res0(1); D3(i)   = res0(2);
  y=par(i); x=lamtot1; z=iisum; 
//res     = D2lapsf(par(i),lamtot1,iisum);
 D13(i)  =   pow(x,y)* pow(x+z,(-y-1))* (y* log(x+z)-y* log(x)-1) ;
 D23(i)  =   y* pow(x,(y-1))* pow(x+z,(-y-2))*(x-y* z) ;

D33(i)  =  y* (y+1)* pow(x,y)*pow((x+z),(-y-2));
D133(i)  =  y* pow(x,y)* pow(x+z,(-y-2))+(y+1)*pow(x,y)* 
            pow(x+z,(-y-2))+y* (y+1)* pow(x,y) *log(x)* 
            pow(x+z,(-y-2))-y *(y+1)* pow(x,y)* pow(x+z,(-y-2))* log(x+z);
D233(i) =  pow(y,2)* (y+1)* pow(x,(y-1))* pow(x+z,(-y-2))+(-y-2)* y* (y+1)* pow(x,y)* pow(x+z,(-y-3));
D333(i) = y* (y+1)* (y+2)* (-pow(x,y))* pow(x+z,(-y-3));


//if (test==3) {
// double nt2 = timer.toc();
// printf("timer-loop 2  %lf antal rv's %d \n",100000*nt2,nn); 
//}

// timer.tic(); 
//if (test==1) { x11.print("x11"); thetades.print("td"); }

//if (test==1) printf(" hej\n"); 


//if (test==1) { msum.print("msum"); mdesi.print("mdesi"); }

// {{{ derivatives of like ,ds dt, dtheta

//mat ddcif1=Dcif, ddcif2=Dcif, ddcif3=Dcif;
//Dcif.fill(0); 
//ddcif1.fill(0); ddcif2.fill(0); 
//
//vec dciflike(Dcif.n_rows), 
//    dcifdt(Dcif.n_rows), dcifds(Dcif.n_rows), 
//    d3like1(Dcif.n_rows), 
//    d3like2(Dcif.n_rows), 
//    d3like3(Dcif.n_rows), 
//    d3like(Dcif.n_rows), dcifdtds(Dcif.n_rows);
//dciflike.fill(0); dcifdt.fill(0); dcifds.fill(0); dcifdtds.fill(0); 

mdesi=trans(thetades.row(i)); //mdesi.print("mdesi");
msum=trans(ags.row(i)); 
dthetaj = (mdesi*D1(i)+msum*D2(i));
dtheta =  dtheta+dthetaj/resv(i);
//
dtj = D3(i)*x1(icause1-1,i);
dt  = dt+dtj/resv(i);
dsj = D3(i)*x2(icause2-1,i);
ds  = ds+dsj/resv(i);

//d3like1=D3(i)*x1.col(i); 
//dciflike=dciflike+d3like1/resv(i); 
//
//d3like2=D33(i)*x1(icause1-1,i)*x1.col(i); 
//dcifdt=dcifdt+(resv(i)*d3like2-dtj*d3like)/pow(resv(i),2);
//d3like3=D33(i)*x2(icause2-1,i)*x2.col(i); 


// {{{  2nd deriv 

dthetadtj =  x1(icause1-1,i)*(mdesi*D13(i)+ msum*D23(i)); 
dtt = dtt+(dthetadtj*resv(i)-dthetaj*dtj)/pow(resv(i),2);

dthetadsj=x2(icause2-1,i)*(mdesi*D13(i)+ msum*D23(i));
numdj=(dthetadsj*resv(i)-dthetaj*dsj); 
dts = dts+numdj/pow(resv(i),2);

dsdtj    = D33(i)*x2(icause2-1,i)* x1(icause1-1,i);
numdsdtj = (dsdtj*resv(i)-dsj*dtj);
dsdt     = dsdt+numdsdtj/pow(resv(i),2);
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
dsdt =  ds*dt*like+like*dsdt;
dtheta =  dtheta *like;
dt = like*dt;
ds = like*ds;
//printf(" %lf %lf \n",dt,ds); 
// }}} 

//if (test==3) {
// double nt3 = timer.toc();
// printf("timer-loop 3 %lf \n",100000*nt3); 
//}

//if (test==1) {
//printf("32 her \n"); printf(" like %lf  \n",like); 
//dtheta.print("dtheta"); dstheta.print("dstheta"); 
//dttheta.print("dttheta"); d3.print("dsttheta"); 
//}


 alllike(0)=like; 
 alllike(1)=-1*dt; 
 alllike(2)=-1*ds; 
 alllike(3)=dsdt; 
 alllike(4)=cause1; 
 alllike(5)=cause1; 
//printf("%d %d %lf %lf %lf %lf \n",cause1,cause2,like,dt,ds,dsdt); 
//alllike.print("all-2"); 

 // allike(4) 
double valr=1; 
if (type==0) {  valr=like;  dp=-1*dtheta; alllike(4)=1; } 
if (type==1) {  valr=-1*dt; dp=dttheta ;  alllike(4)=like;  } 
if (type==2) {  valr=-1*ds; dp=dstheta ;  alllike(4)=like;  } 
if (type==3) {  valr=dsdt;  dp=-1*d3;     alllike(4)=-1*ds; } 

return(valr); 
} // }}}
// }}}

double survivalRVCmarg(vec theta,mat thetades,mat ags,int cause1,vec cif1,mat x1, vec &dp, vec &ddp, vec &alllike) 
{ // {{{
  double lamtot1=1,valr=1,x,y,z;
  int icause1=cause1; // ,icause2;
  if (cause1==0) icause1=1; 

//  int test=0,itest=0;
//  if (itest==1) {
//	  theta.print("theta"); 
//	  thetades.print("theta-des"); 
//	  cif1.print("ci1"); 
//	  x2.print("x1"); 
//	  ags.print("ags"); 
//  }

 vec dL=theta; dL.fill(0); 
 vec par = thetades * theta; 

// if (test==1) { cif1.print("c1"); x1.print("x1"); }
//if (test==1) { theta.print("theta"); thetades.print("t-des "); par.print("pp"); }

int nn=thetades.n_rows; 
int lpar=thetades.n_cols; 
vec sumtheta=ags * theta; 

vec resv(nn); // resv.fill(0); 

//if (test==1) { x1.print("x1"); cif1.print("cif1"); }

vec x1f1; 
x1f1= trans(x1) * cif1; 

//if (test==1) { x1f1.print("x1f1"); } 
//ags.print("ags"); theta.print("par"); 

vec res0(6); 
vec D1(nn),D2(nn),D3(nn), D13(nn), D23(nn), D33(nn); 

vec msum(lpar); //msum.fill(1);
vec mdesi(lpar);

double dtj=0,dt=0;
vec dthetaj(lpar),dtheta(lpar),dtt(lpar),dthetadtj(lpar);
dtheta.fill(0); dtt.fill(0);  // dtdt.fill(0); 
double like=1,iisum; 
int i; 

for (i=0;i<nn;i++) 
{ // {{{ 
lamtot1= sumtheta(i); 
iisum = x1f1(i);
resv(i) = lapsf(par(i),lamtot1,iisum);
like=like*resv(i); 

//vec res(6),res0(6); 
//double y,x,z; 

res0    = Dlapsf( par(i),lamtot1,iisum);
 D1(i)   = res0(0); D2(i)   = res0(1); D3(i)   = res0(2);
//  D1(i)   = res0(0); D2(i)   = res0(1); D3(i)   = res0(2);
  y=par(i); x=lamtot1; z=iisum; 
// D3(i) =(- pow(x,y)*(y/(x+z))*pow(z+x,y))/pow(z+x,(2*y)); 
 if (cause1!=0) {
 D13(i) = pow(x,y)* pow(x+z,(-y-1))* (y* log(x+z)-y* log(x)-1) ;
 D23(i)  =   y* pow(x,(y-1))* pow(x+z,(-y-2))*(x-y* z) ;
// D33(i)  =  y* (y+1)* pow(x,y)*pow((x+z),(-y-2));
 }

//int ncr=Dcif.n_elem; 
//vec D2Dtcif1(ncr),ddcif1(ncr),ddcif2(ncr); 
//Dcif.fill(0); D2Dtcif1.fill(0); ddcif2.fill(0); 

mdesi=trans(thetades.row(i)); //mdesi.print("mdesi");
msum=trans(ags.row(i)); 
dtj = D3(i)*x1(icause1-1,i);
dt  = dt+dtj/resv(i);

dthetaj = (mdesi*D1(i)+msum*D2(i));
dtheta =dtheta+dthetaj/resv(i);

//ddcif1=D3(i)*x1.col(i);
//ddcif2  = ddcif2+ ddcif1/resv(i);
//D2Dtcif1= D2Dtcif1+
//	(resv(i)*D33(i)*x1(icause1-1,i)*x1.col(i)- ddcif1*dtj)/pow(resv(i),2); 

 if (cause1!=0) {
dthetadtj =  x1(icause1-1,i)*(mdesi*D13(i)+ msum*D23(i)); 
dtt = dtt+(dthetadtj*resv(i)-dthetaj*dtj)/pow(resv(i),2);
 }

//dtdtj    = D33(i)*x1(icause1-1,i)* x1(icause1-1,i);
//numdtdtj = (dtdtj*resv(i)-dtj*dtj);
//dtdt     = dtdt+numdtdtj/pow(resv(i),2);

} // }}} 

vec dttheta(lpar); 
if (cause1!=0) dttheta = dtheta*dt*like + like*dtt;
dtheta =  dtheta *like;

//D2Dtcif1= ddcif2*like*dt+ like*D2Dtcif1; 
//ddcif2=ddcif2*like; 
//dtdt =  dt*dt*like+like*dtdt;
dt = like*dt;

// Derivative of marginal caseweight in cif1 direction D_2 (like/ D_t like)
//Dcif= -1*(dt*ddcif2-like*D2Dtcif1)/pow(dt,2); 

alllike(0)=like; alllike(1)=-1*dt; 
//alllike(2)=dtdt;
//alllike(3)=(-dt*dt+like*dtdt)/(dt*dt);
alllike(4)=cause1; alllike(5)=0; 

if (cause1==0 ) { valr=like;  dp=-1*dtheta;   ddp=-1*dtheta; }
if (cause1!=0 ) { valr=-1*dt; dp=dttheta;   ddp=-1*dtheta; }

return(valr); 
} // }}}

RcppExport SEXP survivalRV(SEXP itheta,SEXP istatus1,SEXP istatus2,
	   	     SEXP icif1,SEXP icif2,
                     SEXP irv1, SEXP irv2,SEXP ithetades,
		     SEXP iags, SEXP ivarlink)
{ // {{{

	try {
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
 colvec par = thetades * theta; 

 int lpar=thetades.n_cols; 

vec dp(lpar); dp.fill(0); 
vec all(6); 

double like=0; 
if (status1==0 && status2==0) { // {{{
like=survivalRVC(theta,thetades,ags,0,0,cif1,cif2,x1,x2,dp,all) ;
} // }}}
if (status1==0 && status2!=0) { // {{{
like=survivalRVC(theta,thetades,ags,0,status2,cif1,cif2,x1,x2,dp,all) ;
} // }}}
if (status1!=0 && status2==0) { // {{{
like=survivalRVC(theta,thetades,ags,status1,0,cif1,cif2,x1,x2,dp,all) ;
} // }}}
if (status1!=0 && status2!=0) { // {{{
like=survivalRVC(theta,thetades,ags,status1,status2,cif1,cif2,x1,x2,dp,all);
} // }}}

ressl["like"]=like; 
if (varlink==1) dp=dp % theta;  
ressl["dlike"]=dp;

ressl["theta"]=theta; 
ressl["par.des"]=thetades; 
ressl["varlink"]=varlink; 
ressl["alllike"]=all; 

return(ressl);  
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall

} // }}}

RcppExport SEXP survivalRV2(SEXP itheta,SEXP istatus1,SEXP istatus2,
	   	     SEXP icif1,SEXP icif2,
                     SEXP irv1, SEXP irv2,SEXP ithetades,
		     SEXP iags, SEXP ivarlink)
{ // {{{

	try {
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
 colvec par = thetades * theta; 

 int lpar=thetades.n_cols; 

vec dp(lpar); dp.fill(0); 
vec all(6); 

double like=0; 
if (status1==0 && status2==0) { // {{{
like=survivalRVC2(theta,thetades,ags,0,0,cif1,cif2,x1,x2,dp,all) ;
} // }}}
if (status1==0 && status2!=0) { // {{{
like=survivalRVC2(theta,thetades,ags,0,status2,cif1,cif2,x1,x2,dp,all) ;
} // }}}
if (status1!=0 && status2==0) { // {{{
like=survivalRVC2(theta,thetades,ags,status1,0,cif1,cif2,x1,x2,dp,all) ;
} // }}}
if (status1!=0 && status2!=0) { // {{{
like=survivalRVC2(theta,thetades,ags,status1,status2,cif1,cif2,x1,x2,dp,all);
} // }}}

ressl["like"]=like; 
if (varlink==1) dp=dp % theta;  
ressl["dlike"]=dp;

ressl["theta"]=theta; 
ressl["par.des"]=thetades; 
ressl["varlink"]=varlink; 
ressl["alllike"]=all; 

return(ressl);  
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall

} // }}}


double claytonoakesRVC(vec theta,mat thetades,mat ags,
		       int status1,int status2,double cif1,double cif2,vec x1, vec x2, vec &dp,vec &ccw) 
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
vec iresv(nn); iresv.fill(0); 

double like=1,iisum; 
int i; 
for (i=0;i<nn;i++) 
{
iisum = x1(i)*ii1+x2(i)*ii2;
resv(i) = lapsf(par(i),lamtot1,iisum);
//printf("%lf %lf %lf %lf \n",par(i),lamtot1,iisum,resv(i)); 
iresv(i) = iisum;
like=like*resv(i); 
//printf("===%d  %lf %lf %lf %lf %lf \n",i,x1(i),ii1,x2(i),ii2,like); 
}

//resv.print("resv"); 
//iresv.print("iresv"); 

vec D1(nn),D2(nn),D3(nn); 
vec D13(nn),D23(nn),D33(nn),D133(nn),D233(nn),D333(nn);

vec res(6),res0(6); 
for (i=0;i<nn;i++) 
{ // {{{ 
iisum   = x1(i)*ii1+x2(i)*ii2;
lamtot1= sum( trans(ags.row(i)) % theta); 
res     = D2lapsf(par(i),lamtot1,iisum);
res0    = Dlapsf( par(i),lamtot1,iisum);
//printf(" %lf %lf %lf %lf \n",x1(i),ii1,x2(i),ii2); 
//printf(" %d %lf %lf %lf \n",i,par(i),lamtot1,iisum); 
//res.print("res"); 
//res0.print("res0"); 
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

//D3.print("D3"); 


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
//
//printf(" %lf %lf \n",lamtot1,f1,f2); 
//restd.print("td D2ilapsf"); 
//rest.print("rest D2ilapsf"); 
//
//ressd.print("sd D2ilapsf"); 
//ress.print("ress D2ilapsf"); 

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

vec dtheta(lpar),dtt(lpar), dts(lpar), dtdtds(lpar); 
//    dt(lpar), ds(lpar),dsdt(lpar),
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

dthetad3 =  x1(i)*(mdesi*D13(i)+ msum*D23(i)+ D33(i)*indtheta0); 
dilapthetad3t =  x1(i)*msum*(Di13t+Di23t);
dthetadtj = Di3t*dthetad3+D3(i)*dilapthetad3t;
dtt = dtt+(dthetadtj*resv(i)-dthetaj*dtj)/pow(resv(i),2);

dthetad3=x2(i)*(mdesi*D13(i)+ msum*D23(i)+ D33(i)*indtheta0);
dilapthetad3s = x2(i)*msum*(Di13s+Di23s);
dthetadsj = Di3s*dthetad3+D3(i)*dilapthetad3s;
dts = dts+(dthetadsj*resv(i)-dthetaj*dsj)/pow(resv(i),2);
// }}} 

if (status1==1 && status2==1) { // {{{ 3rd deriv 
dthetad33 = (mdesi*D133(i)+ msum*D233(i)+ D333(i)*indtheta0);
dilapthetad3s = x2(i)*msum*(Di13s+Di23s);
dilapthetad3t = x1(i)*msum*(Di13t+Di23t);
dthetadtdsj = x1(i)*x2(i)*( D33(i)*Di3t*dilapthetad3s+ D33(i)*Di3s*dilapthetad3t+ Di3s*Di3t*dthetad33);
led3 = pow(resv(i),2)*( dthetaj*dsdtj+resv(i)*dthetadtdsj- dthetadtj*dsj-dtj*dthetadsj);
led2 = numdsdtj*2*dthetaj*resv(i);

dtdtds = dtdtds+(led3-led2)/pow(resv[i],4); 
} // }}} 

}


vec dttheta(lpar),dstheta(lpar),d3(lpar); 
dttheta = dtheta*dt*like + like*dtt;
dstheta = dtheta*ds*like + like*dts;
//
d3 = dts*dt*like+ds*dtt*like+ds*dt*dtheta*like+ like*dtdtds+dtheta*dsdt*like;
//
dsdt =  ds*dt*like+like*dsdt;
dtheta =  dtheta *like;
dt = like*dt;
ds = like*ds;

// }}} 

if (status1==0 && status2==0) { // {{{
	valr=like; dp=dtheta;
	ccw(0)=1; 
} // }}}
if (status1==0 && status2==1) { // {{{
	valr=ds; 
	dp=dstheta ;
	ccw(0)=like; 
} // }}}
if (status1==1 && status2==0) { // {{{
	valr=dt; 
	dp=dttheta ;
	ccw(0)=like; 
} // }}}
if (status1==1 && status2==1) { // {{{
	valr=dsdt; 
	dp=d3; 
	ccw(0)=ds; 
} // }}}

return(valr); 
} // }}}


RcppExport SEXP claytonoakesRV(SEXP itheta,SEXP istatus1,SEXP istatus2,SEXP icif1,SEXP icif2,
        SEXP irv1, SEXP irv2,SEXP ithetades,SEXP iags, SEXP ivarlink,SEXP iccw)
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
 vec ccw= Rcpp::as<colvec>(iccw);

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

vec DbetaDtheta(2*lpar); 


if (status1==0 && status2==0) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,0,0,f1,f2,x1,x2,dp,ccw) ;
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
} // }}}
if (status1==0 && status2==1) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,0,1,f1,f2,x1,x2,dp,ccw) ;
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
} // }}}
if (status1==1 && status2==0) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,1,0,f1,f2,x1,x2,dp,ccw) ;
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
} // }}}
if (status1==1 && status2==1) { // {{{
double like=claytonoakesRVC(theta,thetades,ags,1,1,f1,f2,x1,x2,dp,ccw) ;
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
double claytonoakesbinRVC(vec theta,mat thetades,mat ags,int status1,int status2,double cif1,double cif2,vec x1, vec x2, vec &dp,vec &DbetaDtheta,vec &ps,vec &dp00) 
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
	d3 = dts*dt*like+ds*dtt*like+ds*dt*dtheta*like+like*dtdtds+dtheta*dsdt*like;
	//
	dsdt =  ds*dt*like+like*dsdt;
	dtheta =  dtheta *like;
	dt = like*dt;
	ds = like*ds;

	// }}} 

	//printf("%lf %lf  %lf %lf %lf %lf \n",f1,f2,p11,p10,p10,p00); 
	//d3.print("d3"); 
//	dttheta.print("dt"); 
//	dstheta.print("ds"); 
	double p11=like; double p10=f1-p11; double p01=f2-p11; double p00=1-f1-f2+p11; 
	
	DbetaDtheta.subvec(0,lpar-1)=dttheta; 
	DbetaDtheta.subvec(lpar,2*lpar-1)=dstheta; 

	ps(0)=p00; ps(1)=p10; ps(2)=p01; ps(3)=p11; 
	dp00=dtheta; 
	ps(6)=dt-1; ps(7)=ds-1; 

	if (status1==0 && status2==0) { // {{{
		valr=p00; dp=dtheta; ps(4)=dt-1; ps(5)=ds-1; 
	} // }}}
	if (status1==0 && status2==1) { // {{{
		valr=p01; dp=-1*dtheta; DbetaDtheta=-1*DbetaDtheta; 
	        ps(4)=-dt; ps(5)=1-ds; 
	} // }}}
	if (status1==1 && status2==0) { // {{{
		valr=p10; dp=-1*dtheta; DbetaDtheta=-1*DbetaDtheta; 
	        ps(4)=1-dt; ps(5)=-ds; 
	} // }}}
	if (status1==1 && status2==1) { // {{{
		valr=p11; dp=dtheta; ps(4)=dt; ps(5)=ds; 
	} // }}}


	return(valr); 
} // }}}

// for binary case, R version 
RcppExport SEXP claytonoakesbinRV(SEXP itheta,SEXP istatus1,SEXP istatus2,SEXP icif1,SEXP icif2,
	                  	  SEXP irv1, SEXP irv2,SEXP ithetades,SEXP iags, SEXP ivarlink)
{ // {{{
	try {
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
	vec dp00(lpar); dp00.fill(0); 
	vec ps(8); dp.fill(0); 
	vec DbetaDtheta(2*lpar); 
	double like; 

	if (status1==0 && status2==0) { // {{{
		like=claytonoakesbinRVC(theta,thetades,ags,0,0,f1,f2,x1,x2,dp,DbetaDtheta,ps,dp00) ;
	} // }}}
	if (status1==0 && status2==1) { // {{{
		like=claytonoakesbinRVC(theta,thetades,ags,0,1,f1,f2,x1,x2,dp,DbetaDtheta,ps,dp00) ;
	} // }}}
	if (status1==1 && status2==0) { // {{{
		like=claytonoakesbinRVC(theta,thetades,ags,1,0,f1,f2,x1,x2,dp,DbetaDtheta,ps,dp00) ;
	} // }}}
	if (status1==1 && status2==1) { // {{{
		like=claytonoakesbinRVC(theta,thetades,ags,1,1,f1,f2,x1,x2,dp,DbetaDtheta,ps,dp00) ;
	} // }}}
	ressl["like"]=like; 
	if (varlink==1) dp=dp % theta;  
	ressl["dlike"]=dp;
	ressl["ps"]=ps;
	ressl["dp00"]=dp00;
	ressl["theta"]=theta; 
	ressl["par.des"]=thetades; 
	vec obs(4); 
	obs(0)=status1; obs(1)=f1; obs(2)=status2; obs(3)=f2; 
	ressl["obs"]=obs; 
	ressl["varlink"]=varlink; 

	return(ressl);  
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall



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
  try {
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

  int ci,ck,i,j,c,s=0,k,v,c1; 
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

     DUtheta+=weights(i)*sdj*vthetascore*trans(vthetascore);
     vthetascore=weights(i)*diff*vthetascore; 
     Utheta+=vthetascore; 

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
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall


} // }}}

RcppExport SEXP twostageloglikeRV( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, SEXP ithetades,
		SEXP icluster,SEXP iclustsize,SEXP iclusterindex, SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid, SEXP irvdes, SEXP iags,
		SEXP iascertained
) 
{ // {{{
  try {
// {{{ setting matrices and vectors, and exporting to armadillo matrices
 mat thetades = Rcpp::as<mat>(ithetades); 
 mat clusterindex = Rcpp::as<mat>(iclusterindex);
 colvec theta = Rcpp::as<colvec>(itheta);
 colvec clustsize = Rcpp::as<colvec>(iclustsize);
 int antclust = clusterindex.n_rows; 
 IntegerVector cause(icause);
 colvec pmargsurv = Rcpp::as<colvec>(ipmargsurv);
 IntegerVector cluster(icluster);
 colvec weights = Rcpp::as<colvec>(iweights);
// colvec entryage = Rcpp::as<colvec>(ientryage);
 colvec trunkp = Rcpp::as<colvec>(itrunkp);
 colvec secluster = Rcpp::as<colvec>(isecluster);
 mat rvdes= Rcpp::as<mat>(irvdes); 
 mat ags = Rcpp::as<mat>(iags);
 int ascertained= Rcpp::as<int>(iascertained); 

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
//	  cause.print("cause"); 
//	  cluster.print("cluster"); 
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
//      Rprintf("cause %lf \n",mean(cause)); 
//      Rprintf("cluster %lf \n",mean(cluster)); 
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

  int ci,ck,i,j,c,s=0,k,v,c1; 
  double ll=1,Li,Lk,diff=0,loglikecont=0,sdj=0;
  double Lit=1,Lkt=1,llt=1,deppar=1,ssf=0,thetak=0; 
//  double plack(); 
 
  int pt=theta.n_rows; 
  vec dplack(pt); dplack.fill(0);
  vec dplackt(pt); dplackt.fill(0);
//  vec ckij(pt),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  colvec loglikeiid(antclust); 
  colvec trunclikeiid(antclust); 
  if (iid==1) { thetiid.fill(0); 
	  loglikeiid.fill(0); trunclikeiid.fill(0); 
  }

  vec p11tvec(antclust); 
  vec Utheta(pt); 
  vec vthetascore(pt); 
  vec pthetavec(pt); 
  vec vtheta2(pt); 
  mat DUtheta(pt,pt); 
  DUtheta.fill(0); Utheta.fill(0); 
//  if (!Utheta.is_finite()) {  Rprintf(" NA's i def U\n"); Utheta.print("U"); }
//  if (!DUtheta.is_finite()) { Rprintf(" NA's i def DU\n"); DUtheta.print("DU"); }

  int nr=rvdes.n_cols; 
  vec rv2(nr),rv1(nr);
  vec etheta=theta; 
  vec wwc(4); 

  // }}}

for (j=0;j<antclust;j++) if (clustsize(j)>=2) { 

   R_CheckUserInterrupt(); diff=0; 

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
//  if (depmodel==3){ if (varlink==1) etheta=exp(theta); else etheta=theta; }

	if (depmodel==1) { // clayton-oakes  // {{{

           if (trunkp(i)<1 || trunkp(k)<1) {	
		   Lit=trunkp(i); Lkt=trunkp(k); 
		   if ((ascertained==0) || (ascertained==2)) llt=claytonoakes(deppar,0,0,Lit,Lkt,dplackt);
                   if (ascertained==2) llt=1-llt;  // 1-p00, no censoring case
		   if (ascertained==1) llt=claytonoakes(deppar,0,1,Lit,Lkt,dplackt);
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));
	           if (ascertained==2) diff=dplack(0)/ll+dplackt(0)/llt; else diff=dplack(0)/ll-dplackt(0)/llt; 
	   } else {
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
	           diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
//	dplack.print("dp"); 
//	printf("%lf %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",deppar,j,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	   } 
	   if (varlink==1) diff=-pow(deppar,1)*diff;  
	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
	   sdj=pow(diff,2); 
	   // }}}
	} else if (depmodel==3) { //  additive random gamma clayton-oakes  // {{{

	   rv1=trans(rvdes.row(i)); rv2=trans(rvdes.row(k)); 
//	   printf(" %d %d %d %d %d \n",j,i,k,ci,ck);
//         rv1.print("rv1");    rv2.print("rv2"); 
//	   thetades.print("thet"); 

           if (trunkp(i)<1 || trunkp(k)<1) { // {{{ 
		   Lit=trunkp(i); Lkt=trunkp(k); 
//		   llt=claytonoakesRV(theta,thetades,0,0,Lit,Lkt,rv1,rv2,dplackt);
//		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
//
		   if ((ascertained==0) || (ascertained==2)) llt=claytonoakesRVC(etheta,thetades,ags,0,0,Lit,Lkt,rv1,rv2,dplackt,wwc);
                   if (ascertained==2) llt=1-llt;  // 1-p00, no censoring case
		   if (ascertained==1) llt=claytonoakesRVC(etheta,thetades,ags,0,1,Lit,Lkt,rv1,rv2,dplackt,wwc);
		   ll=claytonoakesRVC(etheta,thetades,ags,ci,ck,Li,Lk,rv1,rv2,dplack,wwc);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));

                   if (ascertained==2) vthetascore=dplack/ll+dplackt/llt; else  vthetascore=dplack/ll-dplackt/llt; 
		   // }}}
	   } else {
		   ll=claytonoakesRVC(etheta,thetades,ags,ci,ck,Li,Lk,rv1,rv2,dplackt,wwc);
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
//	dplackt.print("dp"); 
          	   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
//	           if (varlink==1) dplackt=dplackt % etheta;  
	           vthetascore=dplackt/ll; 
	   } 
//	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
//	   sdj=pow(diff,2); 
	   // }}}
	} else if (depmodel==2) { // plackett model  // {{{
        if (trunkp(i)<1 || trunkp(k)<1) {	
           Lit=trunkp(i); Lkt=trunkp(k); 
	   if ((ascertained==0) || (ascertained==2)) llt=placklike(deppar,0,0,Lit,Lkt,dplackt);
           if (ascertained==2) llt=1-llt;  // 1-p00, no censoring case
	   if (ascertained==1) llt=placklike(deppar,0,1,Lit,Lkt,dplackt);
           if (ascertained==1) llt=1-llt;
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(i)*(log(ll)-log(llt));
	   loglikecont=(log(ll)-log(llt));
	   if (ascertained==2) diff=dplack(0)/ll+dplackt(0)/llt; else diff=dplack(0)/ll-dplackt(0)/llt; 
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

//	if (j<10)  printf(" %d  %lf varlink \n",varlink,dplack(0)); 

     if (depmodel!=3) {
//	     printf(" %lf weights \n",weights(i)); 
//		vthetascore.print("vvv 2"); 
	    DUtheta+=weights(i)*diff*diff*vthetascore*trans(vthetascore);
	    vthetascore=weights(i)*diff*vthetascore; 
	    Utheta-=vthetascore; 
	} else  { // additive gamma structure 
	    DUtheta+=weights(i)*vthetascore*trans(vthetascore);
	    vthetascore=weights(i)*vthetascore; 
	    Utheta-=vthetascore; 
//		vthetascore.print("vvv 2"); 
	}

     if (iid==1) { for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)-=vthetascore(c1); 
	           loglikeiid(j)+=loglikecont; 
	           trunclikeiid(j)+=llt; 
     }


     } // }}} strata(i)==strata(k) indenfor strata

  } /* for (c=0....... */   // }}}

} /* j in antclust */ 

//printf("SSum of squares %lf \n",ssf); theta.print("theta"); Utheta.print("Utheta"); DUtheta.print("DUtheta"); 

List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) { res["theta.iid"]=thetiid; 
	      res["loglikeiid"]=loglikeiid; 
	      res["trunclikeiid"]=trunclikeiid; 
            }

return(res); 
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall


} // }}}


//This version of the program specifies which pairs that should be considered 
//in contrast to the general programs that considers all pairs within a cluster
//
//in addition and more importantly the random effects specifiation is within each pair 
//therefore all needed quantities comes for each pair : 
//theta.des and the random effects vectors, cluster and pair id's
//cluster og secluster now also follows the pairs 
//
// left-truncation by giving ptrunc such that truncation probability can be computed
// case control is also a matter of conditioning on marginal status of proband (second in pair)
// L(T1,T2,d1,d2)/L(0,T2,0,d2) propto L(T1,T2,d1,d2)
// ascertainment is equivalent except we here have delayed entry for first component so 
// T2 is subject that leads to ascertainment, d2=1 
// L(T1,T2,d1,d2)/L(T2,T2,0,1) 
// we handle this by giving the survival functions of ascertainment pairs/controls appropriately
RcppExport SEXP twostageloglikeRVpairs( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, 
		SEXP ithetades,SEXP icluster,SEXP iclustsize,SEXP iclusterindex, 
		SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid, SEXP irvdes,
		SEXP idimthetades, SEXP idimrvdes, SEXP inrvs, SEXP iags, 
	        SEXP iascertained	
)  
{ // {{{ 
  try {
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
 int ascertained= Rcpp::as<int>(iascertained); 
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

//  printf("rvdesC %d %d %d \n",arrayDD[0], arrayDD[1], arrayDD[2]); 
  arma::cube rvdesC(rvdesvec.begin(), arrayDD[0], arrayDD[1], arrayDD[2], false);

//  mat B=rvdesC.slice(1); B.print("rv.1"); 
//  mat A=rvdesC.slice(0); A.print("rv.1"); 

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


  int ci,ck,i,j,s=0,k,c1; 
  double ll=1,Li,Lk,sdj=0,diff=0,loglikecont=0;
  double Lit=1,Lkt=1,llt=1,deppar=1,ssf=0,thetak=0; 
//  double plack(); 
 
  vec dplack(pt); dplack.fill(0);
  vec dplackt(pt); dplackt.fill(0);
//  vec ckij(pt),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  colvec loglikeiid(antclust); 
  colvec trunclikeiid(antclust); 
  if (iid==1) { thetiid.fill(0); 
	        loglikeiid.fill(0); 
	        trunclikeiid.fill(0); 
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
  vec wwc(2); 
  // 
  // // }}}
  
colvec likepairs(antclust); 

for (j=0;j<antclust;j++) { 

   R_CheckUserInterrupt(); diff=0; sdj=0; 

// index of subject's in pair "j"
   i=clusterindex(j,0); k=clusterindex(j,1); 
//	  printf("cci 2 %d %d \n",i,k); 
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
//  if (depmodel==3){ if (varlink==1) etheta=exp(theta); else etheta=theta; }

	if (depmodel==1) { // clayton-oakes  // {{{ 

           if (trunkp(i)<1 || trunkp(k)<1) {	
		   Lit=trunkp(i); Lkt=trunkp(k); 
		   if ((ascertained==0) || (ascertained==2)) llt=claytonoakes(deppar,0,0,Lit,Lkt,dplackt);
                   if (ascertained==2) llt=1-llt;  // 1-p00, no censoring case
		   if (ascertained==1) llt=claytonoakes(deppar,0,1,Lit,Lkt,dplackt);
		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));
	           if (ascertained==2) diff=dplack(0)/ll+dplackt(0)/llt; else diff=dplack(0)/ll-dplackt(0)/llt; 
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
	// 3-dimensional array pairs*(2xrandom effects)
        int lnrv= nrvs(j)-1; // number of random effects for this cluster 	
	mat rv=rvdesC.slice(j); 

        vec rv1= trans(rv.submat(span(0),span(0,lnrv)));
        vec rv2= trans(rv.submat(span(1),span(0,lnrv)));


	if (j<-1 ) { rv.print("rv"); Rprintf(" %d \n",lnrv); rv1.print("rv1"); rv2.print("rv2"); }

	// takes parameter relations for each pair
	// 3-dimensional array pairs*(random effects* pars )
	mat thetadesvv=thetadesi.slice(j); 
	mat thetadesv=thetadesvv.rows(0,lnrv); 


	if (j< -10)  {
	   Rprintf("%d %d %d %d %d %lf %lf \n",j,i,k,ci,ck,Li,Lk); 
         rv1.print("rv1");    rv2.print("rv2");    thetadesv.print("thetades "); 
	   etheta.print("e-theta");    ags.print("ags"); 
	}

           if (trunkp(i)<1 || trunkp(k)<1) { 
		   
		   Lit=trunkp(i); Lkt=trunkp(k); 
//		   llt=claytonoakesRV(theta,thetades,0,0,Lit,Lkt,rv1,rv2,dplackt);
//		   ll=claytonoakes(deppar,ci,ck,Li,Lk,dplack);
		   if ((ascertained==0) || (ascertained==2)) llt=claytonoakesRVC(etheta,thetadesv,ags,0,0,Lit,Lkt,rv1,rv2,dplackt,wwc);
                   if (ascertained==2) llt=1-llt;  // 1-p00, no censoring case
		   if (ascertained==1) llt=claytonoakesRVC(etheta,thetadesv,ags,0,1,Lit,Lkt,rv1,rv2,dplackt,wwc);
		   ll=claytonoakesRVC(etheta,thetadesv,ags,ci,ck,Li,Lk,rv1,rv2,dplack,wwc);
//		   dplack.print("dp"); 
//		   dplackt.print("dp---------------t"); 
//	   Rprintf("%d %d %d %d %d %lf %lf %lf %lf %lf %lf  \n",j,i,k,ci,ck,Lit,Lkt,Li,Lk,log(ll),log(llt)); 
		   ssf+=weights(i)*(log(ll)-log(llt));
		   loglikecont=(log(ll)-log(llt));
                   if (ascertained==2) vthetascore=dplack/ll+dplackt/llt; else  vthetascore=dplack/ll-dplackt/llt; 
//		   vthetascore.print("vtheta-score"); 
	   } else {
		   ll=claytonoakesRVC(etheta,thetadesv,ags,ci,ck,Li,Lk,rv1,rv2,dplackt,wwc);
		   ssf+=weights(i)*log(ll); 
		   loglikecont=log(ll);
	           vthetascore=dplackt/ll; 
	   } 
	   // }}} 
	} else if (depmodel==2) { // plackett model  //  // {{{ 
        if (trunkp(i)<1 || trunkp(k)<1) {	
           Lit=trunkp(i); Lkt=trunkp(k); 
	   if ((ascertained==0) || (ascertained==2)) llt=placklike(deppar,0,0,Lit,Lkt,dplackt);
           if (ascertained==2) llt=1-llt;  // 1-p00, no censoring case
	   if (ascertained==1) llt=placklike(deppar,0,1,Lit,Lkt,dplackt);
           ll=placklike(deppar,ci,ck,Li,Lk,dplack);
	   ssf+=weights(j)*(log(ll)-log(llt));
	   loglikecont=(log(ll)-log(llt));
//	   diff=dplack(0)/ll-dplackt(0)/llt; 
	   if (ascertained==2) diff=dplack(0)/ll+dplackt(0)/llt; else diff=dplack(0)/ll-dplackt(0)/llt; 
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
	     DUtheta+=weights(i)*sdj*vthetascore*trans(vthetascore);
	     vthetascore=weights(i)*diff*vthetascore; 
	     Utheta-=vthetascore; 
	} else  { // additive gamma structure 
//		vthetascore.print("vvv"); 
	     DUtheta+=weights(i)*vthetascore*trans(vthetascore);
	     vthetascore=weights(i)*vthetascore; 
	     Utheta-=vthetascore; 
//		vthetascore.print("vvv 2"); 
	}

     if (iid==1) { for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
	           loglikeiid(j)+=loglikecont; 
		   trunclikeiid(j)+=llt; 
     }
     } //  strata(i)==strata(k) indenfor strata

if (iid==1)  likepairs(j)=ll; 

} /* j in antpairs */ 

//printf("Sum of squares %lf \n",ssf); theta.print("theta"); Utheta.print("Utheta"); DUtheta.print("DUtheta"); 

List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) { res["theta.iid"]   =thetiid; 
	      res["loglikeiid"]  =loglikeiid; 
              res["likepairs"]   =likepairs; 
              res["trunclikeiid"]=trunclikeiid; 
            }

return(res); 
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall

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
//
// left-truncation handled by giving ptrunc such that truncation probability can be computed
// case control is also a matter of conditioning on marginal status of proband (second in pair), so similar to normal left truncation
// L(T1,T2,d1,d2)/L(0,T2,0,d2)
// ascertainment is equivalent except we here have delayed entry for first component also a matter of conditioning on marginal status of proband (second in pair)
// L(T1,T2,d1,d2)/L(T2,T2,0,d2) (where by construction T1>T2, since T2 is the first jump)
// we handle this by giving the cumulatives of ascertainment pairs/controls appropriately
RcppExport SEXP survivalloglikeRVpairs( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, 
		SEXP ithetades,SEXP icluster,SEXP iclustsize,SEXP iclusterindex, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid, SEXP irvdes,
		SEXP idimthetades, SEXP idimrvdes, SEXP inrvs ,
		SEXP iags,SEXP ientrycause, SEXP iascertained	
)  
{ // {{{ 
  try {
// {{{ 
//  setting matrices and vectors, and exporting to armadillo matrices
//  // {{{
// colvec nrvs = Rcpp::as<colvec>(inrvs);
 int ascertained= Rcpp::as<int>(iascertained);  // 1= ascertained or casecontrol, controlled via trunkp 
 IntegerVector nrvs(inrvs);
 IntegerVector entrycause(ientrycause);
 IntegerVector cause(icause);
 mat clusterindex = Rcpp::as<mat>(iclusterindex);
 colvec theta = Rcpp::as<colvec>(itheta);
 int pt=theta.n_rows; 
 colvec clustsize = Rcpp::as<colvec>(iclustsize);
 // this is number of pairs (rather than clusters)
 int      antclust   = clusterindex.n_rows; 
 mat      pmargsurv  = Rcpp::as<mat>(ipmargsurv);
 mat      trunkp     = Rcpp::as<mat>(itrunkp);
 mat      ags        = Rcpp::as<mat>(iags);
 colvec   cluster    = Rcpp::as<colvec>(icluster);
 colvec   weights    = Rcpp::as<colvec>(iweights);
 colvec  secluster   = Rcpp::as<colvec>(isecluster);
// mat rvdes= Rcpp::as<mat>(irvdes); 
// colvec entryage = Rcpp::as<colvec>(ientryage);
 int     depmodel    = Rcpp::as<int>(idepmodel); 
 IntegerVector strata(istrata);
 vec all(4); 


// array for derivative of flexible design
 NumericVector DXthetavec(iDXtheta);
 IntegerVector arrayDims(idimDX);
 arma::cube DXtheta(DXthetavec.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

// array for parameter restrictions (one for each pair) pairs * (ant random effects)* (ant par)
 NumericVector thetadesvec(ithetades);
 IntegerVector arrayDims1(idimthetades); 
 IntegerVector arrayDD(3); 

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

 int silent = Rcpp::as<int>(isilent);
 int iid= Rcpp::as<int>(iiid); 
 int antiid = Rcpp::as<int>(iantiid);
 mat Xtheta = Rcpp::as<mat>(iXtheta);
 int ci,ck,i,j,s=0,k,c1; 
 double ll=1,loglikecont=0;
 double llt=1,ssf=0,thetak=0; 
 
 vec dplack(pt); dplack.fill(0);
 vec dplackt(pt); dplackt.fill(0);
  i=silent+1; 

  mat thetiid(antiid,pt); 
  colvec loglikeiid(antiid); 
  colvec trunclikeiid(antiid); 
  if (iid==1) { thetiid.fill(0); 
	        loglikeiid.fill(0); 
		trunclikeiid.fill(0); 
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
  mat Dcif(arrayDims2[0]/2,pt); 
//  Dcif.print("Dcif"); 
  
  // // }}}
  
colvec likepairs(antclust); 
mat matlikepairs(antclust,6); 
vec allvec(6); 

// wall_clock timer; timer.tic(); 

//printf(" thomas \n"); 

for (j=0;j<antclust;j++) { // {{{  

   R_CheckUserInterrupt(); 

// index of subject's in pair "j"
   i=clusterindex(j,0); k=clusterindex(j,1); 

   if (strata(i)==strata(k)) { //  // {{{ 

     // basic survival status 
     ci=cause(i); ck=cause(k); 
     vec Li=trans(pmargsurv.row(i));  // vector of cumulative hazards 
     vec Lk=trans(pmargsurv.row(k)); 

     if (j<-10) { Rprintf("%d %d %d %d  \n",i,k,ci,ck); Li.print("Lk"); Lk.print("Lk"); }
         
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

      etheta=theta; 

      if (depmodel==1) {  // only additive gamma model 
      } else if (depmodel==3) { //  additive random gamma clayton-oakes  // {{{ 

	// takes random effects specification for each pair
	// 3-dimensional array pairs*(2xrandom effects)
        int lnrv= nrvs(j)-1; // number of random effects for this cluster 	

	// first half of rows for person1 and second half for subject 2
	// nn number of competing risks
//        mat rv1= rv.rows(0,nnn/2-1); mat rv2= rv.rows(nnn/2,nnn-1);
      mat rv=rvdesC.slice(j);
      int nnn=rv.n_rows; 
      mat rv1= rv.submat(0,0,nnn/2-1,lnrv);
      mat rv2= rv.submat(nnn/2,0,nnn-1,lnrv);

	if (j<-1 ) { rv1.print("rv1"); rv2.print("rv2"); }

	// takes parameter relations for each pair
	// 3-dimensional array pairs*(random effects* pars )
	mat thetadesvv=thetadesi.slice(j); 
	mat thetadesv=thetadesvv.rows(0,lnrv); 

	if (j<-1)  {
	   Rprintf(" %d %d \n",lnrv,pt); 
           rv1.print("rv1");    rv2.print("rv2"); 
	   thetadesv.print("thetades "); etheta.print("e-theta"); 
//	   mat test=mat(thetades.begin(),3,1); 
//	   test.print("test"); 
	}

           if (any(trunkp.row(i)>0) || any(trunkp.row(k)>0)) { //  /*{{{*/
	      vec Lit=trans(trunkp.row(i)); 
	      vec Lkt=trans(trunkp.row(k)); 
              if (j<-10) {  
		   Rprintf(" så betinges der ! %d %d %d %d \n",entrycause(i),entrycause(k),cause(i),cause(k)); 
		   Lit.print("Lit"); Lkt.print("Lkt"); Lk.print("Lk"); 
              }

              if (ascertained==0) // standard left truncation, both surviving Lit, Lkt
	      llt=survivalRVC2(etheta,thetadesv,ags,0,0,Lit,Lkt,rv1,rv2,dplackt,allvec);
	      if (ascertained==1) // ascertainment correction depending on cumulative hazards of proband Lkt
		                  // case control adjustment  depending on cumulative hazards of proband Lkt
	      llt=survivalRVC2(etheta,thetadesv,ags,0,cause(k),Lit,Lkt,rv1,rv2,dplackt,allvec);
	      ll=survivalRVC2(etheta,thetadesv,ags,cause(i),cause(k),Li,Lk,rv1,rv2,dplack,allvec);

//	      printf("%d %lf %lf \n",j,ll,llt); 
//            printf(" så betinges der ! %d %d %d %d \n",entrycause(i),entrycause(k),cause(i),cause(k)); 
//	      Li.print("Li"); Lk.print("Lk");  Lit.print("Lit"); Lkt.print("Lkt"); 

	      ssf+=weights(i)*(log(ll)-log(llt));
	      loglikecont=(log(ll)-log(llt));
	      vthetascore=dplack/ll-dplackt/llt; 
	   /*}}}*/
	   } else {/*{{{*/
//	      printf("hhhhh %d  %d %d %d %d \n",ascertained,i,k,(int) ci,(int) ck); 
	      ll=survivalRVC2(etheta,thetadesv,ags,cause(i),cause(k),Li,Lk,rv1,rv2,dplack,allvec);
	      ssf+=weights(i)*log(ll); 
	      loglikecont=log(ll);
	      if (j<-10)    { 
		etheta.print("theta"); 
		Li.print("Li"); Lk.print("Lk"); 
		Rprintf("%d %d  %lf %lf \n",cause(i),cause(k),weights(i),ll); 
		allvec.print("allvec"); 
		dplack.print("dll"); 
	      }
	      vthetascore=dplack/ll; 
	   } /*}}}*/
	   // }}} 
	} else if (depmodel==2) { 
        // not possible only additive gamma model 
	} 


        if (depmodel!=3) {
	} else  { // additive gamma structure 
	     DUtheta+=weights(i)*vthetascore*trans(vthetascore);
	     vthetascore=weights(i)*vthetascore; 
	     Utheta=Utheta+vthetascore; 
	}

     if (iid==1) { for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=weights(i)*vthetascore(c1); 
	           loglikeiid(j)+=loglikecont; 
		   trunclikeiid(j)+=llt; 
     }

     } // // }}}  strata(i)==strata(k) indenfor strata

if (iid==1)  { likepairs(j)=ll; matlikepairs.row(j)=trans(allvec); }

} // }}} /* j in antpairs */ 

// double nt2 = timer.toc();
// printf("timer-twostage  %lf \n",nt2); 

//printf("Sum of squares %lf \n",ssf); theta.print("theta"); Utheta.print("Utheta"); DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) { res["theta.iid"]=thetiid; 
	      res["loglikeiid"]=loglikeiid; 
              res["likepairs"]=likepairs; 
	      res["trunclikeiid"]=trunclikeiid; 
              res["all.likepairs"]=matlikepairs; 
            }

return(res); 
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall

} // }}}  


