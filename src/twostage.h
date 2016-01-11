#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

RcppExport SEXP claytonoakesR(SEXP itheta,SEXP  iistatus1,SEXP  iistatus2,SEXP  icif1,SEXP  icif2); 

RcppExport SEXP claytonoakesbinRV(SEXP itheta,SEXP istatus1,SEXP istatus2,SEXP icif1,SEXP icif2,
                     SEXP irv1, SEXP irv2,SEXP ithetades,SEXP iags, SEXP ivarlink); 

double claytonoakesbinRVC(vec theta,mat thetades,mat ags,int status1,int status2,double cif1,double cif2,vec x1, vec x2, vec &dp,vec &DbetaDtheta,vec &ps,vec &dp00) ; 

double survivalRVC(vec theta,mat thetades,mat ags,int cause1,int cause2,vec cif1,vec cif2,mat x1, mat x2, vec &dp, vec &alllike); 

double survivalRVCmarg(vec theta,mat thetades,mat ags,int cause1,vec cif1,mat x1, 
		vec &dp,vec &ddp, vec &alllike); 

double survivalRVC2all(vec theta,mat thetades,mat ags,int cause1,int cause2,vec cif1,vec cif2,mat x1, mat x2, vec &dp, vec &alllike); 

double claytonoakesRVC(vec theta,mat thetades,mat ags, int status1,int status2,
		double cif1,double cif2,vec x1, vec x2, vec &dp,vec  &ccw);

double ilapsf(double y, double x, double z);

double lapsf(double y,double x, double z) ;

vec Dlapsf(double y, double x, double z);

vec D2lapsf(double y, double x, double z) ;

vec Dilapsf(double y, double x, double z) ;

vec D2ilapsf(double y, double x, double z) ;
