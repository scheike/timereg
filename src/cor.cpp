#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace arma;
using namespace Rcpp;

// {{{ laplace and derivatives for structured random cif
double lapgam(double alpha,double beta,double t)
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
return(val); 
}

double ilapgam(double alpha,double beta,double y)
{
double val; 
val=beta*(exp(-log(y)/alpha)-1); 
return(val); 
}

double Dilapgam(double alpha,double beta,double y)
{
double val; 
val=beta*exp(-log(y)/alpha)*(log(y)/(alpha*alpha)); 
return(val); 
}

double Dbetailapgam(double alpha,double beta,double y)
{
double val; 
val=(exp(-log(y)/alpha)-1);
return(val); 
}

double Dalphalapgam(double alpha,double beta,double t)
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
val=(log(beta)-log(beta+t))*val; 
return(val); 
}

double Dbetalapgam(double alpha,double beta,double t)
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
val=alpha*(1/beta-1/(beta+t))*val; 
return(val); 
}


double Dtlapgam(double alpha,double beta,double t)
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
val=-alpha*(1/(beta+t))*val; 
return(val); 
}
// }}}

void ckrvdes(vec &alphai,vec &alphak, // {{{
		double beta, double x,double y,
		vec &ckij, vec &dckij,vec &rvi,vec &rvk)
{ 
double val,val1,val2,val3,alphi=0,alphk=0,alph=0;
double test=1; // lapgam(),ilapgam(),Dtlapgam(),Dalphalapgam(),Dilapgam();
int prv,k; 

if (test<1) {
Rprintf("ckr \n"); 
//print_vec(dckij); print_vec(rvk); print_vec(rvi); print_vec(alphai); print_vec(alphak); 
} 

// fix denne i CPP version
//alphi=trans(rvi) * alphai; alphk= trans(rvk) * alphak; 

//if (test<1) Rprintf("=============================ckr %lf %lf \n",alphi,alphk); 

prv=rvi.n_rows; 
vec Dphi(prv),Dphk(prv);

val=1; 
for (k=0;k<prv;k++) if (rvi(k)+rvk(k)>0) 
{
val1=rvi(k)*ilapgam(alphi,beta,exp(-x))+
     rvk(k)*ilapgam(alphk,beta,exp(-y)); 
if (rvi(k)>0) alph=alphai(k); else alph=alphak(k); 
val1=lapgam(alph,beta,val1); 
val=val*val1; 
}
ckij[0]=1-exp(-x)-exp(-y)+val; 

if (test<1) Rprintf(" %lf ckij \n",ckij[0]); 

val1=0;
for (k=0;k<prv;k++) if (rvi(k)+rvk(k)>0) 
{
if (rvi(k)>0) alph=alphai(k); else alph=alphak(k); 

val2=rvi(k)*ilapgam(alphi,beta,exp(-x))+rvk(k)*ilapgam(alphk,beta,exp(-y)); 
val1= Dtlapgam(alph,beta,val2);
val3= lapgam(alph,beta,val2);

dckij(k)=dckij(k)+Dalphalapgam(alph,beta,val2)/val3;

Dphi=Dphi+(val1*rvi(k)*Dilapgam(alphi,beta,exp(-x))/val3)*rvi; 
//scl_vec_mult(val1*rvi(k)*Dilapgam(alphi,beta,exp(-x))/val3,rvi,Dphi); 
Dphk=Dphk+(val1*rvk(k)*Dilapgam(alphk,beta,exp(-y))/val3)*rvk; 

dckij=dckij+Dphi+Dphk; 
//vec_add(Dphi,dckij,dckij); vec_add(Dphk,dckij,dckij); 
};
dckij=val*dckij; 
//scl_vec_mult(val,dckij,dckij); 

val2=rvi[k]*ilapgam(alphi,beta,exp(-x))+rvk[k]*ilapgam(alphk,beta,exp(-y)); 
val1= Dtlapgam(alph,beta,val2);
val3= lapgam(alph,beta,val2);

dckij[k]=dckij[k]+Dalphalapgam(alph,beta,val2)/val3;

Dphi=val1*rvi*rvi[k]*Dilapgam(alphi,beta,exp(-x))/val3; 
Dphk=val1*rvk*rvk[k]*Dilapgam(alphk,beta,exp(-y))/val3; 

dckij=(dckij+Dphi+Dphk);
dckij=val*dckij; 

//if (test<1) print_vec(dckij); 
//free_vecs(&Dphi,&Dphk,NULL); 
if (test<1) Rprintf("=============================================== ude af cvrks \n");
} // }}}

void funkdes2(vec &alphai,vec &alphak, // {{{
		double beta, double x,double y,
		vec &ckij, vec &dckij,vec &rvi,vec &rvk)
{ 
double val,val1,alphi,alphk,alph,betai,betak;
double test=1;
//       lapgam(),ilapgam(),Dtlapgam(), Dalphalapgam(),Dilapgam(),Dbetalapgam(),Dbetailapgam();
int prv,k; 

if (test<1) {
Rprintf("ckr \n"); 
//print_vec(dckij); print_vec(rvk); print_vec(rvi); print_vec(alphai); print_vec(alphak); 
}

alphi=dot(rvi,alphai); alphk=dot(rvk,alphak); 
betai=alphi; betak=alphk;

if (test<1) Rprintf("=============================ckr %lf %lf \n",alphi,alphk); 

prv=rvk.n_rows; 

val=1; 
for (k=0;k<prv;k++) if (rvi(k)+rvk(k)>0) 
{
val1=rvi(k)*ilapgam(alphi,betai,exp(-x))+
     rvk(k)*ilapgam(alphk,betak,exp(-y)); 
if (rvi(k)>0) alph=alphai(k); else alph=alphak(k); 
val1=lapgam(alph,betai,val1); 
val=val*val1; 
}
ckij[0]=1-exp(-x)-exp(-y)+val; 
} // }}}

void ckrvdes2(vec &alphai,vec &alphak, // {{{
		double beta, double x,double y,
		vec &ckij, vec &dckij,vec &rvi,vec &rvk)
{
double val,val1,val2,val3,alphi,alphk,alph,betai,betak;
double test=1; //lapgam(),ilapgam(),Dtlapgam(), Dalphalapgam(),Dilapgam(),Dbetalapgam(),Dbetailapgam();
int prv,k; 
//void funkdes2(); 

if (test<1) {
Rprintf("ckr \n"); 
//print_vec(dckij); print_vec(rvk); print_vec(rvi); print_vec(alphai); print_vec(alphak); 
}

// fix denne i CPP version
//rvi.print("ckrv rvi"); 
//alphai.print("alph ckrv rvi"); 
alphi=sum(rvi % alphai); alphk=sum(rvk% alphak); 
betai=alphi; betak=alphk;

prv=rvi.n_rows;
vec Dphi(prv),Dphk(prv);
Dphi=0*Dphi; Dphk=0*Dphk; 

val=1; 
for (k=0;k<prv;k++) if (rvi(k)+rvk(k)>0) 
{
val1=rvi(k)*ilapgam(alphi,betai,exp(-x))+
     rvk(k)*ilapgam(alphk,betak,exp(-y)); 
if (rvi(k)>0) alph=alphai(k); else alph=alphak(k); 
val1=lapgam(alph,betai,val1); 
val=val*val1; 
}
ckij[0]=1-exp(-x)-exp(-y)+val; 

for (k=0;k<prv;k++) if (rvi(k)+rvk(k)>0) 
{
if (rvi(k)>0) alph=alphai(k); else alph=alphak(k); 
val2=rvi(k)*ilapgam(alphi,betai,exp(-x))+
     rvk(k)*ilapgam(alphk,betak,exp(-y)); 
val1= Dtlapgam(alph,betai,val2);
val3= lapgam(alph,betai,val2);

dckij(k)=dckij(k)+Dalphalapgam(alph,betai,val2)/val3;
Dphi=(val1*rvi(k)*(Dilapgam(alphi,betai,exp(-x))+Dbetailapgam(alphi,betai,exp(-x)))/val3)*rvi; 
Dphk=(val1*rvk(k)*(Dilapgam(alphk,betak,exp(-y))+Dbetailapgam(alphk,betai,exp(-y)))/val3)*rvk;
dckij=dckij+Dphi+Dphk;
//vec_add(Dphi,dckij,dckij); vec_add(Dphk,dckij,dckij); 

Dphi=(Dbetalapgam(alph,betai,val2)/val3)*rvi;
//vec_add(Dphi,dckij,dckij); 
dckij=dckij+Dphi; 
};
dckij=val*dckij;

} // }}}

// {{{ Laplace for random-cif model 
double laplace(double t,double x)
{ // {{{
double val,val1; 
val=(1+x*t); 
if (val<0) val=0; 
//if (fabs(t)< 0.000000000000001) val1=0; else val1=exp(-log(val)*(1/t)); 
val1=exp(-log(val)*(1/t)); 
// Rprintf("laplace %lf %lf  \n",val,val3); 
return(val1); 
} // }}}

double ilaplace(double t,double y)
{ // {{{
double val,laplace(); 
val=exp(-log(y)*t); val= (val-1)/t;  
// Rprintf("ilaplace y^(1/t)  %lf %lf  \n",exp(log(y)/t),pow(y,1.0/t)); 
// Rprintf("ilaplace  %lf %lf  \n",val,val1); 
return(val); 
} // }}}

double Dilaplace(double theta,double y)
{ // {{{
double val4,val2,val,val1;
val=exp(log(y)/theta); 
val2=-log(y)*val/(theta*theta); 
val4=(1-val)+theta*val2; 
val1=(val4+log(y)*(1-val)/theta)/val;  
return(val1); 
} // }}}

double Dlaplace(double theta,double t)
{ // {{{
double val,val1;

val=1+t/theta; val1=theta*val-log(val); 
val=val1*laplace(theta,t); 
return(val); 
} // }}}

double D2laplace(double theta,double t)
{ // {{{
double val,val1,val2,val3; 

val=1+t/theta; 
val1=theta*val-log(val); 
val3=(t/(theta*theta))/val+(val+t/theta)/(val*val); 
val2=Dlaplace(theta,t)*val1+laplace(theta,t)*val3; 
return(val2); 
} // }}}
// }}}

void ckf(double t,double x,double y,vec &ckij,vec &dckij)
{ // {{{
double val,val2,val3,val4;
//double laplace(),ilaplace(),Dilaplace(),Dlaplace(),D2laplace(); 
double t0;

if (x<0) x=0.0001; if (y<0) y=0.0001; 

val=ilaplace(t,exp(-x))+ilaplace(t,exp(-y)); val2=laplace(t,val); 
ckij[0]=1-exp(-x)-exp(-y)+val2; 

val3=exp(x*t)+exp(y*t)-1; 
val4=val3*log(val3)+exp(x*t)*(-x*t)+exp(y*t)*(-y*t); 

// t0 =exp(-log(t)*2)*exp(log(val3)*(-1/t-1))*val4; 
t0 =pow(1/t,2)*exp(log(val3)*(-1/t-1))*val4; 

dckij[0]=t0; 
} // }}}

void DUetagamma(double t, double x,double y,vec &xi,vec &xk) 
{ // {{{
double y1,y2,t1,val3,val4;

y1=exp(-x); y2=exp(-y); 
val3=exp(x*t)+exp(y*t)-1; 

val4 =exp(log(val3)*(-1/t)); 

t1=val4/(val3); 
//if (isnan(t1)) {
//Rprintf(" missing values in DUetagamma \n"); 
//Rprintf(" t x y val3=exp(x*t)+exp(y*t)-1 %lf %lf %lf %lf  \n",t,x,y,val3); 
////print_vec(xi); 
////print_vec(xk); 
//}; 

xi= (y1-t1*exp(t*x))*xi;
xk= (y2-t1*exp(t*y))*xk; 
xi=xi+xk; 
} // }}}
 
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
		SEXP itheta, SEXP iXtheta, SEXP ithetades,
		SEXP icluster,SEXP iclustsize,SEXP iclusterindex,
		SEXP iinverse,SEXP iCA2, SEXP ix2, // SEXP iz2,
		SEXP isemi2, SEXP iest2,SEXP iZ2gamma2,
//		b2iid, gam2iid, SEXP htheta,SEXP dhtheta,SEXP rhoR,
	        SEXP iflexfunc, SEXP iiid, SEXP isym,SEXP  iweights, 
                SEXP isamecens, SEXP istabcens,SEXP iKMtimes,SEXP isilent,SEXP icifmodel,
		SEXP idepmodel, SEXP iestimator, SEXP ientryage,SEXP icif1entry,SEXP icif2entry,SEXP itrunkp, SEXP irvdes
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
mat thetades = Rcpp::as<mat>(ithetades); // thetades muligvis anderledes for random-cif-des, check
mat clusterindex = Rcpp::as<mat>(iclusterindex);
mat rvdes= Rcpp::as<mat>(irvdes); 
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
 vec cif1lin=-log(1-cif1entry); 

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

  int nr,ci,ck,i,j,c,s,k,v; 
  double Li,Lk,weight=0,p11t,ormarg=0,sdj,diff,cweight2,resp3,time,resp1,resp2;
  double Dinverse=1,DDinverse=1,ddd,edd,ssf=0,response=0,thetak=0,respst=0; 
//  double plack(); 
  vec dplack(4); 
  dplack=0*dplack; 
  int pt=theta.n_rows; 
  vec ckij(4),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

  mat thetiid(antclust,pt); 
  if (iid==1) thetiid=0*thetiid; 

  colvec Utheta(pt); 
  colvec vthetascore(pt); 
  colvec pthetavec(pt); 
  vec vtheta2(pt); 
  mat DUtheta(pt,pt); 
  DUtheta=0*DUtheta; 

  rowvec bhatt2 = est.row(est2.n_cols); 
  colvec pbhat2(z.n_rows); 

  // depmodel=5 
  nr=rvdes.n_cols; 
  vec alphaj(nr),alphai(nr),alpha(nr),
      rvvec(nr),rvvec1(nr),rvvec2vv(nr),rvvec2vt(nr),rvvec2tv(nr);
  vec  rvvec2(nr); 

//  rvvec.print("rvv"); 
//  rvdes.print("rff"); 
//  Xtheta.print("XT"); 

//  Rprintf("semi2 pt depm %d %d %d %d %d \n",semi2,pt,depmodel,CA1,CA2); 
//  est2.print("est2"); 
//  est.print("est"); 
//  z2.print("z2"); 
//  thetades.print("pp"); 
//  Rprintf(" pt depm %d %d \n",pt,depmodel); 

  // }}}

  for (s=0;s<Ntimes;s++) //   if (KMtimes[s]>0) 
  {
          R_CheckUserInterrupt();
	  time=times(s); 
          rowvec bhatt = est.row(s); 
          vec pbhat = z * bhatt; 
	  if ((semi==1) & (cifmodel==1)) pbhat = pbhat + Zgamma*time;
	  if ((semi==1) & (cifmodel==2))  pbhat=pbhat%exp(Zgamma); 
//	  bhatt.print("bhatt");  pbhat.print("pbhatt"); 

	  if ((CA1!=CA2)) {
		  bhatt2 = est2.row(s); 
//		  bhatt2.print("est2"); 
		  pbhat2 = z2 * bhatt2; 
	     if ((semi2==1) & (cifmodel==1)) pbhat2 = pbhat2 + Z2gamma2*time;
	     if ((semi2==1) & (cifmodel==2)) pbhat2=pbhat2%exp(Z2gamma2); 
	  }

    for (j=0;j<antclust;j++) if (clustsize(j)>=2) { 
	  diff=0; sdj=0; 
//	 printf("0 mig \n"); 

	  if (depmodel==5) {
//	     if (j< 2) {  vtheta2.print("================================= theta"); }
             if (inverse==1)  vtheta2=exp(theta); else vtheta2=theta;
//	     if (j< 2) {  printf(" %d %d \n",depmodel,inverse); vtheta2.print(" theta"); }
             alphai= thetades * vtheta2;
	     alphaj= thetades * vtheta2;
	  }
//	 printf("00 mig \n"); 

          for (c=0;c<clustsize[j];c++) for (v=0;v<clustsize[j];v++) // {{{
	  if ((sym==1 && c!=v) || (sym==0 && c<v)) { 
	    i=clusterindex(j,c); k=clusterindex(j,v); 
	    response=0; 

     if ((entryage[i] < time) && (entryage[k]< time)) 
     if ((KMc[i] > 0) && (KMc[k] > 0)) {
	 ci=cause(i); ck=cause(k); 
	 resp1= ((y[k]<=time) && (ck==CA2));
	 resp2= ((y[i]<=time) && (ci==CA1))* ((y[k]<=time) && (ck==CA2));
         
	respst=((y[i]<=entryage[i]) && (ci==CA1))* ((y[k]<=time) && (ck==CA2)) + 
	       ((y[i]<=time) && (ci==CA1))* ((y[k]<=entryage[k]) && (ck==CA2)) ;
//	  Rprintf(" %d %d %d %d %d \n",j,i,k,ci,ck); 

	 if (depmodel!=5)  {
              if (flexfunc==0) thetak= Xtheta(i); else thetak=Xtheta(s,i); 
	      pthetavec= trans(thetades.row(i)); 
	 }
//	 printf("mig \n"); 
//	 else { thetak=evalh(vtheta1,vtime,pthetavec,htheta,rhoR); }
         Li=pbhat(i); Lk=pbhat(k); 
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
	      // DENNE skal tilpasse COR modellen FIXES
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
	     p11t=(1-exp(-Li))*(1-exp(-Lk))*exp(thetak)/(exp(-Li)+exp(thetak)*exp(-Li));
	     ssf+=pow(resp2-p11t,2); 
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
	   ssf+=weight*pow((resp2-exp(thetak)*ormarg),2); 
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
	   ssf+=weight*pow(resp2-p11t,2); 
	}
	} // }}}
	else if (depmodel==4) { // random effects model // {{{

	 if (samecens==1) response=resp2/min(KMc[i],KMc[k]); else
	 {response=resp2/(KMc[i]*KMc[k]);}

	 double ithetak=0; 
	 if (inverse==1) { ithetak=exp(thetak); Dinverse=ithetak; DDinverse=exp(2*thetak); }
       	 else ithetak=thetak; 
//	 if (j< 10) printf("%d  %lf %lf \n",inverse,ithetak,thetak); 

            ckf(ithetak,Li,Lk,ckij,dckij); 
//if (j<10) printf("aaa %d %d %d %d %lf %lf %lf %lf %lf %lf \n",s,j,i,k,thetak,ckij(0),dckij(0),Li,Lk,response); 
            if (trunkp[i]<1) { // {{{
               ckf(ithetak,cif1lin[i],cif1lin[k],ckijvv,dckijvv); 
               ckf(ithetak,Li,cif1lin[k],ckijtv,dckijtv); 
               ckf(ithetak,cif1lin[i],Lk,ckijvt,dckijvt); 
	       ddd=(dckij[0]+dckijvv[0]-dckijtv[0]-dckijvt[0])/trunkp[i]; 
	       edd=(ckij[0]+ckijvv[0]-ckijtv[0]-ckijvt[0])/trunkp[i]; 
               if (inverse==1) response=ddd*Dinverse*(response-edd);   
	       else  response=Dinverse*(response-edd); 
               diff=diff+response; 
               if (inverse==1) sdj=sdj+DDinverse*ddd*ddd; 
               else  sdj=sdj+DDinverse*ddd; 
	       ssf+=pow(response-ckij[0],2);  // }}}
            } else {
	     ssf=ssf+pow(response-ckij[0],2); 
            response=dckij[0]*Dinverse*(response-ckij[0]); 
//            else  response=Dinverse*(response-ckij[0]); 
             diff=diff+response; 
             sdj=sdj-DDinverse*dckij[0]*dckij[0]; 
//             else  sdj=sdj-DDinverse*dckij[0]; 
           }
       } // }}}
	else if (depmodel==5) { // structured random effects model // {{{

	  if (samecens==1) response=resp2/min(KMc[i],KMc[k]); else response=resp2/(KMc[i]*KMc[k]);
	   rvvec=trans(rvdes.row(i)); rvvec1=trans(rvdes.row(k)); 

	if (trunkp[i]<1) { // {{{
	       ckrvdes2(alphai,alphaj,1.0,Li,Lk,ckij,rvvec2,rvvec,rvvec1); 
	       ckrvdes2(alphai,alphaj,1.0,cif1lin[i],cif1lin[k],ckijvv,rvvec2vv,rvvec,rvvec1); 
	       ckrvdes2(alphai,alphaj,1.0,Li,cif1lin[k],ckijtv,rvvec2tv,rvvec,rvvec1); 
	       ckrvdes2(alphai,alphaj,1.0,cif1lin[i],Lk,ckijvt,rvvec2vt,rvvec,rvvec1); 
	       rvvec=0*rvvec; 
	//	  ddd=(dckij[0]+dckijvv[0]-dckijtv[0]-dckijvt[0])/trunkp[i]; 
		  edd=(ckij[0]+ckijvv[0]-ckijtv[0]-ckijvt[0])/trunkp[i]; 
		  rvvec2=rvvec2+rvvec2vv;
//		  vec_add(rvvec2vv,rvvec2,rvvec2); 
		  rvvec2=rvvec2+rvvec2tv;
//		  vec_subtr(rvvec2,rvvec2tv,rvvec2); 
//		  vec_subtr(rvvec2,rvvec2vt,rvvec2); 
		  rvvec2=rvvec2-rvvec2vt;
		  diff=(response-edd); 
		 vthetascore= thetades * rvvec2; 
//		 vM(pardes[i],rvvec2,vthetascore);  // }}}
		} else {
//		printf("2 %d \n",j); 
	       ckrvdes2(alphai,alphaj,1.0,Li,Lk,ckij,rvvec2,rvvec,rvvec1); 
//		printf("3 %d \n",j); 
	       diff=(response-ckij[0]); 
//	       rvvec2.print("rvv2"); 
//	       thetades.print("td"); 
	       vthetascore= trans(thetades) * rvvec2; 
//	       vthetascore.print("vthetascore"); 
//	       vM(pardes[i],rvvec2,vthetascore); 
	       ssf=ssf+pow(response-ckij[0],2); 
	       }
//	vthetascore.print("vt s"); 
//	vtheta2.print("vt 2"); 

              if (inverse==1)  vthetascore=vthetascore%vtheta2;

//              DUtheta=DUtheta+weights[j]*(vhetascore*trans(vthetascore));
	      vthetascore = (weights[j]*diff)*vthetascore; 
	      Utheta=Utheta+vthetascore; 

	      if (iid==1) for (c=0;c<pt;c++) thetiid(j,c)+=vthetascore(c); 
	} // }}}
	else if (depmodel==6) { // random effects model two causes // {{{
            ckf(thetak,Li,Lk,ckij,dckij); 
            if (trunkp[i]<1) {
               ckf(thetak,cif1lin[i],cif1lin[k],ckijvv,dckijvv); 
               ckf(thetak,Li,cif1lin[k],ckijtv,dckijtv); 
               ckf(thetak,cif1lin[i],Lk,ckijvt,dckijvt); 
	       ddd=(dckij[0]+dckijvv[0]-dckijtv[0]-dckijvt[0])/trunkp[i]; 
	       edd=(ckij[0]+ckijvv[0]-ckijtv[0]-ckijvt[0])/trunkp[i]; 
               if (inverse==1) response=ddd*Dinverse*(response-edd);   
	       else  response=Dinverse*(response-edd); 
               diff=diff+response; 
               if (inverse==1) sdj=sdj+DDinverse*ddd*ddd; 
               else  sdj=sdj+DDinverse*ddd; 
	       ssf+=pow(response-ckij[0],2); 
            } else {
            if (inverse==1) response=dckij[0]*Dinverse*(response-ckij[0]); 
            else  response=Dinverse*(response-ckij[0]); 
            diff=diff+response; 
             if (inverse==1) sdj=sdj+DDinverse*dckij[0]*dckij[0]; 
             else  sdj=sdj+DDinverse*dckij[0]; 
	     ssf=ssf+pow(response-ckij[0],2); 
           }
       } // }}}

//if (j<10) printf("uu %d %d %d %lf %lf %lf %lf %lf %lf \n",j,i,k,time,resp1,resp2,y[i],y[k],response); 
//if (j<10) printf("uu2 %lf %lf %lf %lf %lf %lf %d %d \n",pbhat(i),pbhat(k),0*pbhat2(k),response,thetak,ormarg,ci,ck);
}
        } /* for (c=0....... */   // }}}
	
//        if (flexfunc==0) vtheta2=pthetavec;  
//        else  evaldh(vtheta1,vtime,pthetavec,vtheta2,dhtheta,rhoR);

//if (j<10)  printf("%ld %ld %lf %lf\n",s,j,diff,sdj); 
//pthetavec.print("ptvec"); 
//Utheta.print("Utthetavec"); 
//DUtheta.print("DUtvec"); 
   if (depmodel!=5) {
       DUtheta=DUtheta+sdj*weights[j]*(pthetavec*trans(pthetavec));
       vthetascore = (weights[j]*diff)*pthetavec; 
       Utheta=Utheta+vthetascore; 

       if (iid==1) for (c=0;c<pt;c++) thetiid(j,c)+=vthetascore(c); 
   }

//DUtheta.print("DUtvec"); 
//vthetascore.print("vs"); 
//vthetascore.print("vs"); 
//Utheta.print("Utvec"); 
   } /* j in antclust */ 
 } // s < Ntimes

List res; 
//theta.print("theta"); 
//printf(" %lf \n",ssf); 

res["ssf"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
if (iid==1) res["theta.iid"]=thetiid; 

return(res); 
} // }}}


