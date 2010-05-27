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


SEXP mkansvL(double *x,int dim);

double evalh(vector *theta,double *t,vector *xih,SEXP f,SEXP rho); 

void funcevalh(double *theta,double *t,double *x, SEXP f, SEXP rho,double *res,int dimtheta,int dimx); 

void evaldh(vector *theta,double *t,vector *xih,vector *dh,SEXP f,SEXP rho); 

void funcevaldh(double *theta,double *t,double *x, double *res, SEXP f,SEXP rho,int dimtheta,int dimx,int dimres); 

double laplace(double t,double x); 

void ck(double t,double x,double y,double *ckij,double *dckij); 

double mypow(double x,double p); 

void DUetagamma(double t, double x,double y,vector *xi,vector *xk) ;

double laplace(double t,double x); 

double ilaplace(double t,double y); 

double Dilaplace(double theta,double y); 

double Dlaplace(double theta,double t);

double D2laplace(double theta,double t);

