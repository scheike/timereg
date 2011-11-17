//#include <stdio.h>
#include <stdlib.h>     /* declares malloc() */
#include <math.h>
#include <R_ext/BLAS.h>
#include "matrix.h"

/* DGEMV - perform one of the matrix-vector operations */
/* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,  */

extern void F77_SUB(dgemv)(
                const char *trans, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx, const double *beta,
		double *y, const int *incy);

void confBandBasePredict (double *delta, int *nObs, int *nt, int *n,
			  double *se, double *mpt, int *nSims){    
  int nRowDelta = *nObs * *nt;
  // nColDelta = *n
  // se is a vector of length nRowDelta
  // mpt is a vector of length *n
  // pt is a vector of length nRowDelta
  int i,j,k; // dummy variables, counters
  // The next line does: double g[*n]; // vector of IID random normals
  double *g = (double *)malloc(*n * sizeof(double));
  // The next line does: double pt[nRowDelta];
  double *pt = (double *)malloc(nRowDelta * sizeof(double));
  double pt1, pt2; // temporary variables used while calculating maxima 
  // Some parameters to give to DGEMV in the BLAS library
  char trans = 'n';
  double alpha = 1.0;
  double beta = 0.0;
  int incx = 1;
  int incy = 1;
  double norm_rand(); 
  void GetRNGstate(),PutRNGstate();  

  GetRNGstate();

  for(i = 0; i < *nSims; i++){ // Number of draws
    // First generate IID random normal vector of length *n
    for(j = 0; j < *n; j++){
      g[j] = norm_rand();
    }
    // Matrix multiplication:
    // pt := delta %*% g
    
    F77_CALL(dgemv)(&trans, &nRowDelta, n, &alpha, 
                    delta, &nRowDelta, g, &incx, &beta, pt, &incy);
    
    for(k = 0; k < *nObs; k++){
      pt1 = -1.0e99; // initially set to -INF
      for(j = 0; j < *nt; j++){
	pt2 = fabs(pt[k * *nt + j])/se[k * *nt + j];
	if(pt1 < pt2){
	  pt1 = pt2;
	}
      }
      mpt[i * *nObs + k] = pt1;
    }
  }

  PutRNGstate();

  // prevent memory leaks by unallocating memory allocated by malloc
  free(g);
  free(pt);

}

