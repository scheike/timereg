//#include <stdio.h>
#include <math.h>
#include "matrix.h"

double tukey(x,b) double x,b;
{
  return((1/b)*((cos(3.141592 *(x/b))+ 1)/2) * (fabs(x/b) < 1)); 
}

double dtukey(x,b) double x,b;
{
  return((-3.141592/b*b)*(sin(3.141592 *(x/b))/2)*(fabs(x/b) < 1));
}

void smoothB(designX,nx,p,bhat,nb,b,degree,coef)
double *designX,*bhat,*b;
int *coef,*nx,*p,*degree,*nb;
{ // {{{ 
  matrix *mat1,*mat2,*II,*I;
  vector *XWy,*Y,*RES,*sY;
  int count,j,k,s,d;
  int silent=1;
  double tukey(),x,w,band;
  matrix *sm1,*sm2;

  malloc_mat(*nx,(*degree)+1,mat1);
  malloc_mat(*nx,(*degree)+1,mat2);
  malloc_mat(*nx,(*degree)+1,sm1);
  malloc_mat(*nx,(*degree)+1,sm2);
  malloc_vec(*nx,Y);
  malloc_vec(*nx,sY);
  malloc_vec((*degree)+1,XWy);
  malloc_vec((*degree)+1,RES);
  malloc_mat((*degree)+1,(*degree)+1,II);
  malloc_mat((*degree)+1,(*degree)+1,I);

  for (s=0;s<*nb;s++){
    x=bhat[s]; 
    for (k=1;k<*p;k++)  { 
      vec_zeros(Y); 
      mat_zeros(mat1); 
      mat_zeros(mat2); 
      count=0; 
      vec_zeros(RES); 
      band=b[(k-1)*(*nb)+s]; 
      /* Rprintf("band %lf %ld \n",band,k);  */

      for (j=0;j<*nx;j++) {
	if (fabs(designX[j]-x)<band) {
	  /* Rprintf("smooth %lf %lf \n",band,designX[k*(*nx)+j]);  */
	  w=tukey(designX[j]-x,band); 
	  ME(mat1,count,0)=1.0; 
	  ME(mat2,count,0)=w; 
	  for (d=1;d<=*degree;d++) {
	    ME(mat1,count,d)=pow(designX[j]-x,d); 
	    ME(mat2,count,d)=w*ME(mat1,count,d);
	  }
   
	  VE(Y,count)=w*designX[k*(*nx)+j];
	  count=count+1; 
	}
      }
      /*
	mat1=m_resize(mat1,count,*degree+1);
	mat2=m_resize(mat2,count,*degree+1);
	Y=v_resize(Y,count); m_output(mat2); 
      */

      if (count>=4) {
	MtA(mat1,mat2,II); 
	invertS(II,I,silent); 
	vM(mat1,Y,XWy);
	vM(I,XWy,RES); 
      };
      bhat[k*(*nb)+s]=VE(RES,*coef);  
    } /* components */ 
  } /* times */ 
  free_mat(sm1); free_mat(sm2); free_mat(mat1); free_mat(mat2); 
  free_mat(I); free_mat(II); 
  free_vec(sY); free_vec(Y); free_vec(XWy); free_vec(RES); 
} // }}}


void localTimeReg(designX,nx,p,times,response,bhat,nb,b,lin,dens)
double *designX,*bhat,*b,*times,*response,*dens;
int *nx,*p,*nb,*lin;
{
  matrix *X,*AI,*A;
  vector *res,*Y,*XY;
  int c,j,k,s,silent=1;
  double band,tukey(),dtukey(),x,w,delta; 
  j=(*lin+1)*(*p);  
  malloc_mat(*nx,j,X);
  malloc_mat(j,j,A);
  malloc_mat(j,j,AI);
  malloc_vec(*nx,Y);
  malloc_vec(j,XY);
  malloc_vec(j,res);

  /* Rprintf("enters Local Time Regression \n");  */
  for (s=0;s<*nb;s++){
    x=bhat[s]; 
    for (c=0;c<*nx;c++){
      delta=times[c]-x;
      band=b[s]; 
      w=tukey(delta,band); 
      dens[s]=dens[s]+w;
      dens[(*nb)+s]=dens[(*nb)+s]+dtukey(delta,b[s]); 
      for(j=0;j<*p;j++) {
	ME(X,c,j)=designX[j*(*nx)+c]*sqrt(w);
	if (*lin>=1) ME(X,c,*p+j)=designX[j*(*nx)+c]*delta*sqrt(w);
	if (*lin>=2) ME(X,c,2*(*p)+j)=delta*ME(X,c,*p+j);
	if (*lin==3) ME(X,c,3*(*p)+j)=delta*ME(X,c,2*(*p)+j);
      }
      VE(Y,c)=response[c]*sqrt(w);  
    }
    dens[s]=dens[s]/(*nx); dens[(*nb)+s]=dens[(*nb)+s]/(*nx);

    MtA(X,X,A); 
    invertS(A,AI,silent); 
    if (ME(AI,0,0)==0.0){
      Rprintf("Non-invertible design in local smoothing at time %lf \n",x); 
    }
    vM(X,Y,XY); 
    Mv(AI,XY,res); 

    for (k=1;k<((*lin)+1)*(*p)+1;k++){ 
      bhat[k*(*nb)+s]=VE(res,k-1); 
    }
  }

  free_mat(A); free_mat(AI); free_mat(X); 
  free_vec(Y); free_vec(XY); free_vec(res); 
}
