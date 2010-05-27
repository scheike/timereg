#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "haplosurv.h" 

void genoMLEdes( // {{{
antpers,Nit,detail,nph,oh,
nphpp,haplopars,rho,score,d2lgeno, 
Iinv,varpar,step,lm,minlm,haplodes,alpha,dimhap,haplofreq,iid)
double *Iinv,*score,*haplopars,
       *rho,*d2lgeno,*varpar,*step,*lm,*minlm,*haplodes,*alpha,
       *haplofreq,*iid; 
int *antpers,*Nit,*detail,*nph,*oh,*nphpp,*dimhap;
// }}}
{
// {{{ variable defs
  matrix *AA,*S1f,*S1hap,*S2,*S1,*SI,*VU,*SIhapDes,*hapDes; 
  vector *Vhaplopars,*Uf,*U,*delta,*pars,*Uhap; 
  vector *scorei[*antpers]; 
  int nphm1,i,j,k,l,it,dimpar; 
  nphm1=nph[0]-1; dimpar=nph[0]-1;
  int amount[1]; 
  unsigned ii=nphm1*nphm1; 
  double maxdelt,LogLikegeno[1],sumscore=10.0,fabs(); 
  double *scoregenoC= calloc(nphm1, sizeof(double));
  double *d2lgenoC= calloc(ii, sizeof(double));
//  void genoLogLikeHp(); 
// }}}


// {{{ matrix allocation
  malloc_mats(dimpar,*dimhap,&SIhapDes,&hapDes,NULL); 
  malloc_mats(*dimhap,dimpar,NULL); 
  malloc_mat(dimpar,dimpar,S1f); 
  malloc_mat(*dimhap,*dimhap,S1hap); 
  malloc_vecs(dimpar,&Uf,&Vhaplopars,NULL); 
  malloc_vec(*dimhap,Uhap); 
  malloc_mat(*dimhap,nphm1,AA); 

  malloc_mats(*dimhap,*dimhap,&S2,&S1,&SI,&VU,NULL); 
  malloc_vecs(*dimhap,&U,&delta,&pars,NULL); 
  for(j=0;j<*antpers;j++)  {malloc_vec(dimpar,scorei[j]);
			    }
// }}}

  for (i=0;i<nphm1;i++)
  for(j=0;j<*dimhap;j++)  ME(hapDes,i,j)=haplodes[j*(nphm1)+i];
  for (j=0;j<*dimhap;j++)  VE(pars,j)=alpha[j];

  /* Main procedure ================================== */
  for (it=0;it<*Nit;it++) {
  maxdelt=0; 

  // initializing {{{ score and haplopars 
  vec_zeros(U); mat_zeros(S1);  

  for (j=0;j<*dimhap;j++) VE(Uhap,j)=VE(pars,j);
  Mv(hapDes,Uhap,Vhaplopars); 
  for (j=0;j<nphm1;j++) { 
      haplopars[j]=VE(Vhaplopars,j);
      haplofreq[j]=exp(haplopars[j]); 
      haplofreq[nphm1]=haplofreq[nphm1]+haplofreq[j]; }
      haplofreq[nphm1]=1+haplofreq[nphm1];
   for (j=0;j<nphm1;j++) { 
	 haplofreq[j]=haplofreq[j]/haplofreq[nphm1];} 
	 haplofreq[nphm1]=1/haplofreq[nphm1]; 
  
  // }}}
  
  // genologlike  {{{
  for (j=0;j<nphm1;j++) {scoregenoC[j]=0; 
    for (l=0;l<nphm1;l++) d2lgenoC[j*(nphm1)+l]=0; } 
    amount[0]=3; 
    genoLogLikeHp(oh,nphpp,nph,amount,antpers,haplopars,rho,LogLikegeno,scoregenoC,scorei,d2lgenoC); 

   for (j=0;j<nphm1;j++) VE(Uf,j)=scoregenoC[j];
   for (j=0;j<nphm1;j++) for (l=0;l<nphm1;l++) ME(S1f,j,l)=-d2lgenoC[j*(nphm1)+l]; 

    //  printf(" for ================== \n"); 
    //  print_vec(Uf); print_mat(hapDes); 

    vM(hapDes,Uf,U); MtA(hapDes,S1f,AA); MxA(AA,hapDes,S1);

    // print_vec(U); print_mat(S1); 
    //  printf("================== \n"); 
// }}}
    
      if (*lm>0) { /* Levenberg-Marquardt algorithm */ 
if (*detail>=1) printf("Levenberg-Marquardt steps, sumscore %d %lf \n",it,sumscore); 
         mat_transp(S1,S2); MxA(S1,S1,S2); 
         for (k=0;k<*dimhap;k++) ME(S2,k,k)=ME(S2,k,k)+lm[0]; // *ME(S2,k,k); 
         invert(S2,SI); MxA(SI,S1,S2); Mv(S2,U,delta); 
       }
      else {
      /* Newton-Raphson step */ 
     if (*detail>=1)  printf(" Newton-Raphson steps, sumscore %d %lf \n",it,sumscore); 
        invert(S1,SI); Mv(SI,U,delta); 
      }
       scl_vec_mult(step[0],delta,delta); 

      for (k=0;k<*dimhap;k++) {if (maxdelt<fabs(VE(delta,k))) maxdelt=fabs(VE(delta,k));}
     if (*detail>=1)  printf(" maximum delta value  %d %lf \n",it,maxdelt); 

       if (sumscore<2.05) lm[0]=minlm[0]; 

       if (*detail>=2) { 
         printf("====================Iteration %d ==================== \n",it);
         printf("Estimate beta \n"); print_vec(pars); 
         printf("Score D l\n"); print_vec(U); 
         printf("Information D^2 l\n"); print_mat(SI); 
         printf("simple D2 l\n");  print_mat(S1); 
         printf("delta \n"); print_vec(delta); 
       }
       vec_add(pars,delta,pars); 
       sumscore=0; 
       for (k=0;k<*dimhap;k++) {sumscore +=fabs(VE(U,k));}
       if ((fabs(sumscore)<0.000000001) & (it<*Nit-2)){ it=*Nit-2; }

  } /* it */

   for (j=0;j<nphm1;j++) {scoregenoC[j]=0; 
      for (l=0;l<nphm1;l++) d2lgenoC[j*(nphm1)+l]=0; }

   amount[0]=2;
   genoLogLikeHp(oh,nphpp,nph,amount,antpers,haplopars,rho,LogLikegeno,scoregenoC,scorei,d2lgenoC);
     //
    
   MxA(hapDes,SI,SIhapDes); 
   for (j=0;j<*antpers;j++) {
	   vM(SIhapDes,scorei[j],delta);  
           for(k=0;k<*dimhap;k++)  { iid[k*(*antpers)+j]=VE(delta,k); 
           for(i=0;i<*dimhap;i++) ME(VU,i,k)+= VE(delta,i)*VE(delta,k); }
   }

  // returning some arguments to R  {{{
    for(j=0;j<*dimhap;j++)  {
         score[j]=VE(U,j); alpha[j]=VE(pars,j);
      for (k=0;k<*dimhap;k++){Iinv[k*(*dimhap)+j]=ME(SI,j,k);
	      varpar[k*(*dimhap)+j]=ME(VU,j,k);
      }
   }
  // }}}

// {{{ freeing all memory 
  free_mats(&SIhapDes,&hapDes,&S1f, &S1hap, &AA, &S2,&S1,&SI,&VU,NULL); 
  free_vecs(&Uf, &Vhaplopars, &Uhap, &U,&delta,&pars,NULL); 

  for (j=0;j<*antpers;j++) { free_vec(scorei[j]); }
  free(scoregenoC); free(d2lgenoC); 
// }}}

}
