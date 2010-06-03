#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <R.h>
#include <Rinternals.h>

void genoLogLike(int *oh, int *nphppp, int *nph, int *np, 
		 double *hf, double *rho, double *logLike, double *score, 
		 vector *personScoreLL[*np], double *d2l);

void genoLogLikeHp(int *oh, int *nphppp, int *nph, int *amountToReturn, 
		   int *np, double *hp, double *rho, double *logLike, 
		   double *score, vector *personScoreLL[*np], double *d2l);  

SEXP mkansv(double *x,int dim);

void designeval(double *x, double *h, SEXP f, double *res, SEXP rho, int dimx,int dimres);

void haplodesign(vector *xi,int *haplotype,vector *xih,SEXP f,SEXP rho); 

void fDUeta(int c2,int *oh,int *nph,int *hapdim,
	    double *haplofreq,matrix *hapDes,vector *etascore);

void Fdesigneval(double *x, double *z,double *h, SEXP f, double *res, SEXP rho, int dimx,int dimz, int dimres); 

void Fhaplodes(vector *xi,int *haplotype,vector *xih,SEXP f,SEXP rho) ;

void FhaplodesMM(vector *xi,vector *zi,int *haplotype,vector *xih,SEXP f,SEXP rho); 

void haplodesignMM(vector *xi,vector *zi,
		int *haplotypeD,int *haplotypeP, vector *xih, SEXP f,SEXP rho);

void designevalMM(double *x, double *z,double *hD, double *hP, 
	SEXP f, double *res, SEXP rho, int dimx,int dimz,int dimres);

SEXP mkansvM(double *x,int dim);

void designevalM(double *x, double *hD, double *hP, SEXP f, double *res, SEXP rho, int dimx,int dimres);

void haplodesignM(vector *xi,int *haplotypeD,int *haplotypeP, vector *xih, SEXP f,SEXP rho);

void semihaplo( double *times,int *Ntimes,double *x,int *delta,int *cause,double *KMc,double *z,int *antpers,int *px,int *Nit,
double *score,double *hess,double *est,double *var,int *sim,int *antsim,int *rani,double *test,double *testOBS,double *Ut,double *simUt,int *weighted,
double *gamma,double *vargamma,int *semi,double *zsem,int *pg,int *trans,double *gamma2,int *CA,int *line,int *detail,double *biid,
double *gamiid,int *resample,double *timepow,int *clusters,int *antclust,
double *haplofreq,double *alphaiid,int *hapdim,
int *nph,int *oh,int *nphpp,SEXP designfuncX,SEXP designfuncZ,SEXP rhoR,int *dimxih,int *dimzih,double *haplodes,int *fixhaplo,int *designtest);

//void semihaplo(
//times,Ntimes,x,delta,cause,KMc,z,antpers,px,Nit,
//score,hess,est,var,sim,antsim,rani,test,testOBS,Ut,simUt,weighted,
//gamma,vargamma,semi,zsem,pg,trans,gamma2,CA,line,detail,biid,
//gamiid,resample,timepow,clusters,antclust,
//haplofreq,alphaiid,hapdim,
//nph,oh,nphpp,designfuncX,designfuncZ,rhoR,dimxih,dimzih,haplodes,fixhaplo,
//designtest)
//double *times,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
//*Ut,*simUt,*gamma,*zsem,*vargamma,*gamma2,*biid,*gamiid,*timepow,
//*haplofreq,*alphaiid,*haplodes;
//int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*rani,*weighted,
//*semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,
//*nph,*oh,*nphpp,*hapdim,*dimxih,*dimzih,*fixhaplo,*designtest;
//SEXP designfuncX,designfuncZ,rhoR; 


void crfitsemimatch( double *times,int *Ntimes,double *x,int *delta,int *cause,
        double *KMc,double *z,int *antpers,int *px,int *Nit,
        double *score,double *hess,double *est,double *var,int *sim,
	int *antsim,int *rani,double *test,double *testOBS,double *Ut,
	double *simUt,int *weighted, double *gamma,double *vargamma,int *semi,
	double *zsem,int *pg,int *trans,double *gamma2,int *CA,
	int *line,int *detail,double *biid, double *gamiid,int *resample,
	double *timepow,int *clusters,int *antclust,double *haplofreq,double *alphaiid,
	int *hapdim, int *nph,int *oh,int *nphpp,SEXP designfuncX,
	SEXP designfuncZ,SEXP rhoR,int *dimxih,int *dimzih,double *haplodes,
	int *fixhaplo,int *designtest);
//double *times,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
//*Ut,*simUt,*gamma,*zsem,*vargamma,*gamma2,*biid,*gamiid,*timepow,
//*haplofreq,*alphaiid,*haplodes;
//int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*rani,*weighted,
//*semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,
//*nph,*oh,*nphpp,*hapdim,*dimxih,*dimzih,*fixhaplo,*designtest;


