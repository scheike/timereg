#include<R.h>
#include<Rmath.h>

double F77_SUB(phid)(double *x){
  return pnorm(*x, 0.0, 1.0, 1, 0);
  
}

double F77_SUB(studnt)(int *nu, double *x){
  return pt(*x, *nu, 1, 0);
}
