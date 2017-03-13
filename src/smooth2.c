//#include <stdio.h>
#include <math.h>
#include "matrix.h"

void smooth2B(designX,nx,p,bhat,nb,b,degree,coef)
double *designX,*bhat,*b;
int *coef,*nx,*p,*degree,*nb;
{
  matrix *mat1,*mat2,*I,*XWy,*Y,*sm1,*sm2,*sY,*RES;
  matrix *sm1sm2t; // not in original
  int med,j,k,s,count,starti=0,d;
  double x,w;

  malloc_mats(*nx,*degree+1,&mat1,&mat2,NULL);
  malloc_mats(*nx,*p-1,&Y,NULL);
  malloc_mats((*degree+1),*p-1,&XWy,&RES,NULL);
  malloc_mats((*degree+1),*degree+1,&I,NULL);

  for (s=0;s<*nb;s++){
    med=0; x=bhat[s]; count=0; 

    for (j=starti;((j<*nx) && (designX[j]<x+*b));j++) {
      if ((designX[j]>x-(*b)) && (med==0)) {med=1; starti=j;}

      if (fabs(designX[j]-x)<*b) {
	w=tukey(designX[j]-x,*b);/*Rprintf("%lf %lf \n",designX[j]-x,w);*/ 
	ME(mat1,count,0)=1.0; 
	ME(mat2,count,0)=w; 
	for (d=1;d<=*degree;d++) {
	  ME(mat1,count,d)=pow(designX[j]-x,d); 
	  ME(mat2,count,d)=w*ME(mat1,count,d);
	}

	for (k=1;k<*p;k++){
	  ME(Y,count,k-1)=w*designX[k*(*nx)+j];
	}
	count=count+1; 
      }
    }
    /* */
    malloc_mats(count,*degree+1,&sm1,&sm2,NULL);
    malloc_mats(count,*p-1,&sY,NULL);
    malloc_mat(count,count,sm1sm2t);
   
    mat_subsec(mat1,0,0,count-1,*degree,sm1);   
    mat_subsec(mat2,0,0,count-1,*degree,sm2);   
    mat_subsec(Y,0,0,count-1,*p-2,sY); 

    MtA(sm1,sm2,sm1sm2t);
    invert(sm1sm2t,I); 
    MtA(sm1,sY,XWy);
    MxA(I,XWy,RES); 

    for (k=1;k<*p;k++){
      bhat[k*(*nb)+s]=ME(RES,*coef,k-1); 
    }
    free_mats(&sm1,&sm2,&sY,sm1sm2t,NULL); 
  }
  free_mats(&mat1,&mat2,&Y,&XWy,&RES,&I,NULL); 
}
