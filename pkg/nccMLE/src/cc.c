#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
void casecohort(designX,nx,px,T,status,cohort,betaS,Nit,w,loglike,Vbeta,RVbeta,detail,breslow,nbre)
double *designX,*betaS,*w,*loglike,*Vbeta,*RVbeta,*T,*breslow; 
int *nx,*px,*Nit,*detail,*status,*nbre,*cohort;
{
  matrix *dU,*UI,*UI2,*S2,*E2,*TMP,*dU2; 
  vector *xis,*U,*beta,*xi[*nx],*Ui,*Di,*S1,*Wpi,*pi,*Yi,*delta;
  int i,j,k,s,it,count=0;
  double xb,S0,ll,sumU,
         *xbhat=calloc(*nx,sizeof(double));

malloc_mats(*px,*px,&TMP,&S2,&dU,&UI,&UI2,&E2,&dU2,NULL);
malloc_vecs(*px,&xis,&S1,&beta,&U,&Ui,&Di,&delta,NULL); 
malloc_vecs(*nx,&Wpi,&pi,&Yi,NULL); 
for (j=0;j<*nx;j++) malloc_vec(*px,xi[j]);

for(j=0;j<*px;j++) VE(beta,j)=betaS[j]; 
for (j=0;j<*nx;j++) { for (i=0;i<*px;i++) VE(xi[j],i)=designX[i*(*nx)+j];}

for (it=0;it<*Nit;it++)
{ 
vec_zeros(U); mat_zeros(dU); ll=0; count=0; sumU=0; mat_zeros(dU2); 

for (j=0;j<*nx;j++) { 
vec_star(beta,xi[j],xis); xbhat[j]=vec_sum(xis);
}

for (s=0;s<*nx;s++)
{  
if (status[s]==1)
{
  S0=0; vec_zeros(S1); mat_zeros(S2); 

  for (j=s;j<*nx;j++)
  {
   if ((cohort[j]==1))  /* Prentice Self */
   {
   xb=xbhat[j]; S0=S0+exp(xb);
   scl_vec_mult(exp(xb),xi[j],xis); vec_add(S1,xis,S1);

   for(k=0;k<*px;k++) 
   for(i=0;i<*px;i++) ME(S2,k,i)=ME(S2,k,i)+
				  VE(xi[j],k)*VE(xi[j],i)*exp(xb);
   } /* cohort */
  }
   scl_vec_mult(1/S0,S1,S1); scl_mat_mult(1/S0,S2,S2); 

   for(k=0;k<*px;k++) for(i=0;i<*px;i++) ME(E2,k,i)=VE(S1,k)*VE(S1,i);

   vec_subtr(xi[s],S1,Ui); vec_add(U,Ui,U); 
   mat_subtr(S2,E2,TMP); mat_add(dU,TMP,dU); 

   for(k=0;k<*px;k++) 
   for(i=0;i<*px;i++) ME(dU2,k,i)=ME(dU2,k,i)-VE(Ui,k)*VE(Ui,i);
   breslow[0*(*nbre)+count]=T[s];
   breslow[1*(*nbre)+count]=breslow[1*(*nbre)+count-1]+(1/S0); 
/* */
   count=count+1; 
} /* if status[s]==1 */ 
} /* s nx*/ 

invert(dU,UI); Mv(UI,U,delta); /* print_mat(SI); print_vec(U); */ 
if (*Nit>1) vec_add(beta,delta,beta); 

for(j=0;j<*px;j++) sumU=sumU+fabs(VE(U,j)); 
/* if (sumU>0.00001) it=it-1; else it=*Nit-1; */

if (*detail==1) { 
printf("======================COX=========================\n"); 
print_mat(UI); print_mat(dU);  
printf(" score er "); print_vec(U); 
printf(" beta er "); print_vec(beta); };
} /* it */


/* Baseline beregnes i slut punkt */
count=0; 
for (j=0;j<*nx;j++) {vec_star(beta,xi[j],xis); xbhat[j]=vec_sum(xis); }

for (s=0;s<*nx;s++)
{  
if (status[s]==1)
{
  S0=0; vec_zeros(S1); 

  for (j=s;j<*nx;j++) {
   if ((cohort[j]==1)) {
   xb=xbhat[j]; S0=S0+exp(xb); } }

   breslow[0*(*nbre)+count]=T[s];
   breslow[1*(*nbre)+count]=breslow[1*(*nbre)+count-1]+(1/S0); 
   count=count+1; 
} /* if status[s]==1 */ 
} /* s nx*/ 


for(j=0;j<*px;j++) {betaS[j]= VE(beta,j); loglike[0]=ll;
for (k=0;k<*px;k++) {
     Vbeta[k*(*px)+j]=-ME(UI,j,k); RVbeta[k*(*px)+j]=ME(UI2,j,k); } 
}

free_mats(&E2,&S2,&UI,&UI2,&dU,&TMP,&dU2,NULL);
free_vecs(&xis,&U,&beta,&pi,&Yi,&Di,&Wpi,&S1,&Ui,&delta,NULL); 
for(j=0;j<*nx;j++) free_vec(xi[j]); 

free(xbhat); 
}

void chenlo(designX,nx,px,T,status,cohort,betaS,Nit,npop,ncohort,ncases,ncoca)
double *designX,*betaS,*T; 
int *nx,*px,*Nit,*status,*cohort,*npop,*ncases,*ncohort,*ncoca;
{
  matrix *dU,*UI,*UI2,*S2,*S2ca,*S2co,*E2,*TMP,*dU2; 
  vector *xis,*U,*beta,*xi[*nx],*Ui,*Di,*S1,*S1ca,*S1co,*Wpi,*pi,*Yi,*delta;
  int i,j,k,s,it,count=0,detail=0;
  double N,n,N1,N0,n1,n0, // xbhat[*nx]; 
  *xbhat=calloc(*nx,sizeof(double));
  double xb,S0,S0ca,S0co,ll,sumU;
  
N=*npop; n=*ncohort; N1=*ncases; n1=*ncoca; n0=n-n1; N0=N-N1; 

/* printf(" %lf %lf %lf %lf \n",N,n,n1,n0);  */

malloc_mats(*px,*px,&TMP,&S2,&S2co,&S2ca,&dU,&UI,&UI2,&E2,&dU2,NULL);
malloc_vecs(*px,&xis,&S1co,&S1,&S1ca,&beta,&U,&Ui,&Di,&delta,NULL); 
malloc_vecs(*nx,&Wpi,&pi,&Yi,NULL); 
for (j=0;j<*nx;j++) malloc_vec(*px,xi[j]);

for(j=0;j<*px;j++) VE(beta,j)=betaS[j]; 
for (j=0;j<*nx;j++) for (i=0;i<*px;i++) VE(xi[j],i)=designX[i*(*nx)+j];

for (it=0;it<*Nit;it++)
{ 
vec_zeros(U); mat_zeros(dU); ll=0; count=0; sumU=0; mat_zeros(dU2); 
for (j=0;j<*nx;j++) {vec_star(beta,xi[j],xis); xbhat[j]=vec_sum(xis);}

for (s=0;s<*nx;s++)
{  
if (status[s]==1)
{
  S0ca=0; vec_zeros(S1ca); mat_zeros(S2ca); 
  S0co=0; vec_zeros(S1co); mat_zeros(S2co); 

  for (j=s;j<*nx;j++)
  {
   if ( (cohort[j]==1) & (status[j]==0) ) {
   xb=xbhat[j]; S0co=S0co+exp(xb); 
   scl_vec_mult(exp(xb),xi[j],xis); vec_add(S1co,xis,S1co);
   for(k=0;k<*px;k++) for(i=0;i<*px;i++) ME(S2co,k,i)=ME(S2co,k,i)+
				  VE(xi[j],k)*VE(xi[j],i)*exp(xb);
   } /* non-cases cohort */
   if (status[j]==1) {
   xb=xbhat[j]; S0ca=S0ca+exp(xb); 
   scl_vec_mult(exp(xb),xi[j],xis); vec_add(S1ca,xis,S1ca);
   for(k=0;k<*px;k++) for(i=0;i<*px;i++) ME(S2ca,k,i)=ME(S2ca,k,i)+
				  VE(xi[j],k)*VE(xi[j],i)*exp(xb);
   } /* cases */
  }
   S0=(S0ca/N)+(N0/(N*n0))*S0co; 
   scl_vec_mult(1/N,S1ca,S1ca); 
   scl_vec_mult(N0/(N*n0),S1co,S1co); 
   vec_add(S1ca,S1co,S1); scl_vec_mult(1/S0,S1,S1); 

   scl_mat_mult(1/N,S2ca,S2ca); 
   scl_mat_mult(N0/(N*n0),S2co,S2co); 
   mat_add(S2ca,S2co,S2); 
   scl_mat_mult(1/S0,S2,S2); 

   for(k=0;k<*px;k++) for(i=0;i<*px;i++) ME(E2,k,i)=VE(S1,k)*VE(S1,i);

   vec_subtr(xi[s],S1,Ui); vec_add(U,Ui,U); 
   mat_subtr(S2,E2,TMP); mat_add(dU,TMP,dU); 

   for(k=0;k<*px;k++) 
   for(i=0;i<*px;i++) ME(dU2,k,i)=ME(dU2,k,i)-VE(Ui,k)*VE(Ui,i);
/* breslow[0*(*nbre)+count]=T[s];
   breslow[1*(*nbre)+count]=breslow[1*(*nbre)+count-1]+(1/S0); 
*/
   count=count+1; 
} /* if status[s]==1 */ 
} /* s nx*/ 
/*
//v_output(S1co); v_output(S1ca); 
print_vec(S1); print_vec(S1co); print_vec(S1ca); 
printf(" %lf %lf %lf \n",S0,S0ca,S0co); 
*/

invert(dU,UI); Mv(UI,U,delta); /* print_mat(SI); print_vec(U); */ 
vec_add(beta,delta,beta); 

for(j=0;j<*px;j++) sumU=sumU+fabs(VE(U,j)); 
if (sumU>0.000001) it=it-1; else it=*Nit-1; 

if (detail==1) { 
printf("======================COX=========================\n"); 
print_mat(UI); print_mat(dU);  
printf(" score er "); print_vec(U); 
printf(" beta er "); print_vec(beta); };
} /* it */


for(j=0;j<*px;j++) {betaS[j]= VE(beta,j); 
/*
loglike[0]=ll;
for (k=0;k<*px;k++) {
  //Vbeta[k*(*px)+j]=-UI->me[j][k]; RVbeta[k*(*px)+j]=UI2->me[j][k]; } 
  Vbeta[k*(*px)+j]=-ME(UI,j,k); RVbeta[k*(*px)+j]=ME(UI2,j,k); } 
*/
}

free_mats(&E2,&S2,&S2co,&S2ca,&UI,&UI2,&dU,&TMP,&dU2,NULL);
free_vecs(&xis,&U,&beta,&pi,&Yi,&Di,&Wpi,&S1,&S1ca,&S1co,&Ui,&delta,NULL); 
for(j=0;j<*nx;j++) free_vec(xi[j]); 
}


void emCC(designX,nx,px,T,status,
   	tau,nno,betaS,Nit,loglike,
        Uii,score,detail, breslow,nbre,
	p,W,antW,oWj,betait,
	delt,ntot)
double *designX,*betaS,*loglike,*Uii,*score,*T,*breslow,*p,*tau,*W,
*delt,*oWj; 
int *nx,*px,*Nit,*detail,*status,*nbre,*nno,*antW,*betait,*ntot;
{
  matrix *dU,*UI,*UI2,*S2,*E2,*TMP,*dU2; 
  vector *xis,*U,*beta,*xi[*nx],*Ui,*Di,*S1,*Wpi,*pi,*Yi,*delta;
  vector *Wj[*antW],*betagam;
  int oit,i,j,k,s,it,*imin=calloc(1,sizeof(int)),count=0;
  double dummy,xb,S0,ll,sumU=0; // ,Wbeta[*antW],xbhat[*nx];
  double rralpha,sumB=0,sumP=0,sumtot; 
  double  *Wbeta=calloc(antW[0],sizeof(double)),
          *xbhat=calloc(nx[0],sizeof(double)),
          *alphaij=calloc(antW[0],sizeof(double)),
          *gamp=calloc(antW[0],sizeof(double)), 
          *gambreslow=calloc(nbre[0],sizeof(double));

if (*detail==1) printf("i emcc programmet %d %d %d %d \n",*nbre,*nno,*ntot,*antW); 

for (j=0;j<*antW;j++) gamp[j]=0;  
for (j=0;j<*nbre;j++) gambreslow[j]=0;  

malloc_mats(*px,*px,&TMP,&S2,&dU,&UI,&UI2,&E2,&dU2,NULL);
malloc_vecs(*px,&xis,&S1,&betagam,&beta,&U,&Ui,&Di,&delta,NULL); 
malloc_vecs(*nx,&Wpi,&pi,&Yi,NULL); 
for (j=0;j<*nx;j++) malloc_vec(*px,xi[j]);
for (j=0;j<*antW;j++) malloc_vec(*px,Wj[j]);

for (j=0;j<*px;j++) {VE(beta,j)=betaS[j];VE(betagam,j)=betaS[j];}
for (j=0;j<*nx;j++) for (i=0;i<*px;i++) VE(xi[j],i)=designX[i*(*nx)+j];
for (j=0;j<*antW;j++) for (i=0;i<*px;i++) VE(Wj[j],i)=W[i*(*antW)+j];

dummy=0; 
for (j=0;j<*antW;j++) {
vec_star(beta,Wj[j],xis); Wbeta[j]=vec_sum(xis);
alphaij[j]=exp(-breslow[1*(*nbre)+*nbre-1]*exp(Wbeta[j]))*p[j]; 
dummy=dummy+alphaij[j]; }
for (j=0;j<*antW;j++) {alphaij[j]=alphaij[j]/dummy; } 
for (j=0;j<*antW;j++) {p[j]=((*nno)*alphaij[j]+oWj[j])/(*ntot);} 

for (oit=0;oit<*Nit;oit++) { 
for (it=0;it<*betait;it++) { 
vec_zeros(U); mat_zeros(dU); ll=0; count=0; sumU=0; mat_zeros(dU2); sumB=0; sumP=0; 
for (j=0;j<*nx;j++) {vec_star(beta,xi[j],xis); xbhat[j]=vec_sum(xis);}
for (j=0;j<*antW;j++) {vec_star(beta,Wj[j],xis); Wbeta[j]=vec_sum(xis);} 

for (s=0;s<*nx;s++) {  
  S0=0; vec_zeros(S1); mat_zeros(S2); 
if ((status[s]==1)) {

  for (j=s;j<*nx;j++)
  {
   xb=xbhat[j]; S0=S0+exp(xb);
   scl_vec_mult(exp(xb),xi[j],xis); vec_add(S1,xis,S1);

   for(k=0;k<*px;k++) 
   for(i=0;i<*px;i++) ME(S2,k,i)=ME(S2,k,i)+
				  VE(xi[j],k)*VE(xi[j],i)*exp(xb);
  } /* observerede :  cohort og case */

  for (j=0;j<*antW;j++) 
  { 
  rralpha=exp(Wbeta[j])*(*nno)*alphaij[j]; 
  S0=S0+rralpha; 
  scl_vec_mult(rralpha,Wj[j],xis); 
  vec_add(S1,xis,S1);

  for(k=0;k<*px;k++) 
  for(i=0;i<*px;i++) ME(S2,k,i)=ME(S2,k,i)+
		  VE(Wj[j],k)*VE(Wj[j],i)*rralpha; 
  } /* bidrag fra uobserverede */  

   scl_vec_mult(1/S0,S1,S1); scl_mat_mult(1/S0,S2,S2); 

   for(k=0;k<*px;k++) for(i=0;i<*px;i++) ME(E2,k,i)=VE(S1,k)*VE(S1,i);

   vec_subtr(xi[s],S1,Ui); vec_add(U,Ui,U); 
   mat_subtr(S2,E2,TMP); mat_add(dU,TMP,dU); 

   for(k=0;k<*px;k++) 
   for(i=0;i<*px;i++) ME(dU2,k,i)=ME(dU2,k,i)-VE(Ui,k)*VE(Ui,i);

   if (*betait>1) {
   breslow[0*(*nbre)+count]=T[s];
   if (count>0) breslow[1*(*nbre)+count]=breslow[1*(*nbre)+count-1]+(1/S0); 
   else breslow[1*(*nbre)]=(1/S0); 
   }; 

   count=count+1; 
} /* if status[s]==1 */ 
} /* s nx*/ 
invert(dU,UI); Mv(UI,U,delta); /* print_mat(SI); print_vec(U); */ 

if (*betait>1) vec_add(beta,delta,beta); 
for(j=0;j<*px;j++) sumU=sumU+fabs(VE(U,j)); 

if (*detail==1) { 
printf(" CC==== nested CC ====NPMLE ======================== \n"); 
printf(" Nit betait  %d %d \n",oit,it); 
print_mat(UI); print_mat(dU);  
printf(" score er "); print_vec(U); 
printf(" beta er "); print_vec(beta);  }

if ((sumU<0.001) & (*betait>1) ) break; 
} /* betait */

dummy=0; 
for (j=0;j<*antW;j++) {
vec_star(beta,Wj[j],xis); Wbeta[j]=vec_sum(xis);
alphaij[j]=exp(-breslow[1*(*nbre)+*nbre-1]*exp(Wbeta[j]))*p[j]; 
dummy=dummy+alphaij[j]; }
for (j=0;j<*antW;j++) {alphaij[j]=alphaij[j]/dummy; } 
for (j=0;j<*antW;j++) {p[j]=((*nno)*alphaij[j]+oWj[j])/(*ntot);} 

vec_subtr(betagam,beta,delta); 
for(j=0;j<*px;j++) sumU=sumU+fabs(VE(delta,j)); 
for(j=0;j<*antW;j++) sumP=sumP+fabs(gamp[j]-p[j]); 
for(j=0;j<*nbre;j++) sumB=sumB+fabs(breslow[1*(*nbre)+j]-gambreslow[j]); 
if (*betait>1) sumtot=sumU+sumP+sumB; else sumtot=sumP+sumB;

if (*detail==1) { 
printf(" sumU sumP sumB %lf %lf  %lf \n",sumU,sumP,sumB); 
};
if ((sumtot<0.000001) ) oit=*Nit-1; 

for(j=0;j<*px;j++) VE(betagam,j)=VE(beta,j);
for(j=0;j<*nbre;j++) gambreslow[j]=breslow[1*(*nbre)+j]; 
for(j=0;j<*antW;j++) gamp[j]=p[j]; 
/* printf(" %lf %lf %d \n",sumP,sumB,oit);  */
} /* Nit */

for(j=0;j<*px;j++) {
betaS[j]= VE(beta,j); loglike[0]=ll;
delt[j]= VE(delta,j); score[j]=VE(U,j); 
for (k=0;k<*px;k++) {Uii[k*(*px)+j]=-ME(UI,j,k);} 
}

free_mats(&E2,&S2,&UI,&UI2,&dU,&TMP,&dU2,NULL);
free_vecs(&xis,&U,&beta,&pi,&Yi,&Di,&Wpi,&S1,&Ui,&delta,NULL); 
for(j=0;j<*nx;j++) free_vec(xi[j]); for(j=0;j<*antW;j++) free_vec(Wj[j]); 

  free(Wbeta); free(xbhat); free(alphaij);free(gamp);free(gambreslow);free(imin); 

}

void varemCC(designX,nx,px,T,status,
	     tau,nno,betaS,Nit,loglike,
	     Vbeta,score,detail,breslow,nbre,
	     p,W,antW,oWj,betait,delt,
	     ntot,antstrat,stratum,nstrat,d,
	     indi,fullZ,pz,antbase,basestratum,
	     konv,juul)
double *designX,*betaS,*loglike,*Vbeta,*score,*T,*breslow,*p,*tau,*W,*d,
*delt,*oWj,*fullZ,*konv; 
int *nx,*px,*Nit,*detail,*status,*nbre,*nno,*antW,*betait,*ntot,*stratum,
*antstrat,*nstrat,*indi,*pz,*antbase,*basestratum,*juul;
{
double *betaj=calloc(*px,sizeof(double)),
*Uii=calloc((*px)*(*px),sizeof(double));
int j,k,*one=calloc(1,sizeof(int)),
    *nul=calloc(1,sizeof(int)); 
one[0]=1; nul[0]=0; 
void emCCindi(),emCC(); 

if (*detail>2) printf("i emcc programmet %d %d %d %d \n",*nbre,*nno,*ntot,*antW); 

for (j=0;j<*px;j++) { 
for (k=0;k<*px;k++) betaj[k]=betaS[k]+(k==j)*(*d); 

if (*detail>2) { 
   printf(" EM aided differentation \n"); 
   printf(" MLE of beta \n"); 
   printf(" direction %d out of %d \n",j,*px); 
   for (k=0;k<*px;k++) printf(" %lf \n",betaj[k]); 
   printf(" \n");  }

if (*indi==0)
emCC(designX,nx,px,T,status,tau,nno,betaj,Nit,loglike,Uii,
score,detail,breslow,nbre,p,W,antW,oWj,one,delt,ntot);
else if (*indi==1) 
emCCindi(designX,nx,px,T,status,tau,nno,betaj,Nit,loglike,Uii,score,detail,
breslow,nbre,p,W,antW,oWj,one,delt,ntot,antstrat,stratum,nstrat,fullZ,pz,
antbase,basestratum,konv,juul); 

for (k=0;k<*px;k++) { Vbeta[k*(*px)+j]=score[k]/(*d); }  
} 

if (*detail>2) {
for (j=0;j<*px;j++)  printf("score %lf ",score[j]); printf(" \n");   
for (j=0;j<*px;j++)  
{for (k=0;k<*px;k++) printf("score/d %lf ",Vbeta[k*(*px)+j]); printf(" \n"); }  
}

free(betaj); free(Uii); 
}  /*varem() ---------------------------------- */



void
emCCindi(designX,nx,px,T,status,tau,nno,betaS,Nit,loglike,
Uii,score,detail,breslow,nbre,p,W,antW,oWj,betait,delt,ntot,
antstrat,stratum,nstrat,
fullZ,pz,antbase,basestratum,konv,ju)
double *designX,*betaS,*loglike,*score,*T,*breslow,*p,*tau,*W,*delt,
*oWj,*fullZ,*konv,*Uii; 
int *nx,*px,*Nit,*detail,*status,*nbre,*nno,*antW,*betait,*ntot,
*antstrat,*stratum,*nstrat,*pz,*antbase,*basestratum,*ju;
{
if (*detail>2) printf(" i emcc programmet \n"); 
  matrix *dU,*UI,*UI2,*S2,*E2,*TMP,*dU2; 
  vector *xis,*U,*beta,*xi[*nx],*Ui,*Di,*S1,*Wpi,*pi,*Yi,*delta;
  vector *Wj[*antW],*betagam,*zi[*nno],*betaz,*zis,*betax;
  vector *alphaij[*nno],*ta,*xisuz; 
  int nobs,juul,bases=0,ppx,medx,ppz,medz,oit,i,j,k,s,m,it, difdim,
      *imin=calloc(2,sizeof(int));   
  double dummy,xb,S0,ll,sumU=0,rralpha,sumB=0,sumP=0,sumtot,Zuobeta=0; 
  int count,
      *timeind=calloc(*nno,sizeof(int)),
      *nnstrat=calloc(*antstrat,sizeof(int));
  double *Wbeta=calloc(antW[0],sizeof(double)),
         *Zbeta=calloc(*nno,sizeof(double)),
         *xbhat=calloc(nx[0],sizeof(double)),
         *zbhatz=calloc(*nno,sizeof(double)),
         *alphaj=calloc(antW[0]*antstrat[0],sizeof(double)),
         *sumalpha=calloc(antstrat[0],sizeof(double)),
         *gamp=calloc(antW[0]*antstrat[0],sizeof(double)), 
         *gambreslow=calloc(nbre[0]*antbase[0],sizeof(double));

if (*detail>2) printf(" %d %d %d %d  \n",*nbre,*antbase,*antW,*antstrat); 
if (*detail>2) printf(" %d %d %d %d %d \n",*px,*nx,*antW,*nno,*ntot); 

juul=*ju; nobs=ntot[0]-nno[0]; difdim=px[0]-pz[0]; 

/* if (juul==1) difdim=1;  */
if (*pz==0) {medz=0; ppz=1;} else {ppz=pz[0]; medz=1;} 
if (*px==0) {medx=0; px[0]=pz[0];difdim=1;} else {ppx=px[0]; medx=1;} 

malloc_vec(ppz,betaz); malloc_vec(ppz,zis); 
for (j=0;j<*antstrat;j++) nnstrat[j]=nstrat[j]; 

malloc_mats(*px,*px,&TMP,&S2,&dU,&UI,&UI2,&E2,&dU2,NULL);
malloc_vecs(*px,&xis,&S1,&betagam,&beta,&U,&Ui,&Di,&delta,NULL); 
malloc_vecs(*nx,&Wpi,&pi,&Yi,NULL); 
for (j=0;j<*nx;j++) malloc_vec(*px,xi[j]); 
for (j=0;j<*antW;j++) malloc_vec(difdim,Wj[j]);
for (j=0;j<*nno;j++) {malloc_vec(*antW,alphaij[j]); malloc_vec(ppz,zi[j]);} 

malloc_vec(*nbre,ta); 
malloc_vec(difdim,betax); malloc_vec(difdim,xisuz); 

juul=*ju; 

for (j=0;j<*px;j++) {VE(beta,j)=betaS[j];VE(betagam,j)=betaS[j];}
for (j=0;j<difdim;j++) VE(betax,j)=betaS[j];
if (*pz==0) VE(betaz,0)=0.0; else {
for (j=0;j<*pz;j++) VE(betaz,j)=betaS[difdim+j];}; 

for (j=0;j<*nx;j++) for (i=0;i<*px;i++) VE(xi[j],i)=designX[i*(*nx)+j];
for (j=0;j<*antW;j++) for (i=0;i<difdim;i++) VE(Wj[j],i)=W[i*(*antW)+j]; 

if (medz==1) {
for (j=0;j<*nno;j++) for (i=0;i<*pz;i++) VE(zi[j],i)=fullZ[i*(*nno)+j];}  

for (i=0;i<*nno;i++) {for(j=0;j<*nbre;j++) VE(ta,j)=fabs(breslow[j]-tau[i]); 
dummy=vec_min(ta,imin); timeind[i]=*imin; }

/* ========================================================== */
for (i=0;i<*antstrat;i++) for (j=0;j<*antW;j++) alphaj[i*(*antW)+j]=0; 
for (j=0;j<*nx;j++)  {vec_star(beta,xi[j],xis); xbhat[j]=vec_sum(xis); }
for (j=0;j<*antW;j++) {vec_star(betax,Wj[j],xisuz);Wbeta[j]=vec_sum(xisuz);}
for (j=0;j<*nno;j++)  {vec_star(betaz,zi[j],zis);zbhatz[j]=vec_sum(zis); }

for (i=0;i<*nno;i++) { dummy=0; 
if (juul==0) {bases=basestratum[nobs+i]; } 
for (j=0;j<*antW;j++) {
if (juul==0) {Zuobeta=Wbeta[j]; } 
if (juul==1) {bases=j+1; Zuobeta=xbhat[nobs+i];} 

if (i< -1987 ) 
//%lf %lf %d %d %d %lf %d %d \n",alphaij[i]->ve[j], 
printf(" %lf %lf %d %d %d %lf %d %d \n",VE(alphaij[i],j),
breslow[bases*(*nbre)+timeind[i]], bases,*nbre,timeind[i],
    p[(stratum[nobs+i]-1)*(*antW)+j],
    stratum[nobs+i],*antW); 

VE(alphaij[i],j)=exp(-breslow[bases*(*nbre)+timeind[i]]*
    exp(Zuobeta+zbhatz[i]))*p[(stratum[nobs+i]-1)*(*antW)+j]; 
dummy=dummy+VE(alphaij[i],j);}

for (j=0;j<*antW;j++) {VE(alphaij[i],j)= VE(alphaij[i],j)/dummy;
alphaj[(stratum[nobs+i]-1)*(*antW)+j] 
=alphaj[(stratum[nobs+i]-1)*(*antW)+j]+VE(alphaij[i],j);} }  

for (i=0;i<*antstrat;i++) for (j=0;j<*antW;j++) 
p[i*(*antW)+j]=(alphaj[i*(*antW)+j]+oWj[i*(*antW)+j])/nstrat[i]; 

/* ========================================================== */
for (oit=0;oit<*Nit;oit++)
{ 

for (it=0;it<*betait;it++)
{ 
vec_zeros(U); mat_zeros(dU); ll=0; count=0; sumU=0; mat_zeros(dU2); sumB=0; sumP=0; 
for (j=0;j<difdim;j++) {VE(betax,j)=VE(beta,j);}
if (medz==1) for (j=0;j<*pz;j++) VE(betaz,j)=VE(beta,difdim+j); 

for (j=0;j<*nx;j++) {vec_star(beta,xi[j],xis); xbhat[j]=vec_sum(xis); }
/* printf(" %lf ",xbhat[j]); printf(" \n"); */
for (j=0;j<*antW;j++) {vec_star(betax,Wj[j],xisuz); Wbeta[j]=vec_sum(xisuz);} 
for (j=0;j<*nno;j++)  {vec_star(betaz,zi[j],zis);   zbhatz[j]=vec_sum(zis); }

for (s=0;s<nobs;s++)
{
/* printf(" %d s \n",s);  */
if (status[s]==1) {
  S0=0; vec_zeros(S1); mat_zeros(S2); 

  for (j=s;j<nobs;j++)
  {if (basestratum[j]==basestratum[s]) {
   xb=xbhat[j]; S0=S0+exp(xb);
   scl_vec_mult(exp(xb),xi[j],xis); vec_add(S1,xis,S1);

   for(k=0;k<*px;k++) 
   for(i=0;i<*px;i++) ME(S2,k,i)=ME(S2,k,i)+
                                   //xi[j]->ve[k]*xi[j]->ve[i]*exp(xb);} 
                                   VE(xi[j],k)*VE(xi[j],i)*exp(xb);}
  } /* observerede :  cohort og case */

if (juul==0) { 
  for (i=0;i<*nno;i++) { 
  if ((basestratum[nobs+i]==basestratum[s]) & (tau[i]>T[s])) {
  for (j=0;j<*antW;j++) {
  rralpha=exp(Wbeta[j]+zbhatz[i])*VE(alphaij[i],j); 
  S0=S0+rralpha; 
  
  for (m=0;m<difdim;m++) VE(Di,m)=VE(Wj[j],m); 
  if (medz==1) for (m=difdim;m<*px;m++) VE(Di,m)=VE(zi[i],m-difdim); 
  scl_vec_mult(rralpha,Di,xis); vec_add(S1,xis,S1);

  for(k=0;k<*px;k++) 
  for(m=0;m<*px;m++) ME(S2,k,m)=ME(S2,k,m)+ VE(Di,k)*VE(Di,m)*rralpha; 
  } }
  } /* bidrag fra uobserverede */  
}  else { 
  for (i=0;i<*nno;i++) { 
  if (tau[i]>T[s]) {
  j=basestratum[s]-1; 
  rralpha=exp(xbhat[nobs+i]+zbhatz[i])*VE(alphaij[i],j); 
  S0=S0+rralpha; 
  
  for (m=0;m<difdim;m++) VE(Di,m)=VE(xi[nobs+i],m); 
  if (medz==1) for (m=difdim;m<*px;m++) VE(Di,m)=VE(zi[i],m-difdim); 
/*
printf(" %lf %lf %lf \n",xbhat[nobs+i],zbhatz[i],VE(alphaij[i],j)); 
print_vec(Di); 
*/
  scl_vec_mult(rralpha,Di,xis); vec_add(S1,xis,S1);

  for(k=0;k<*px;k++) 
  for(m=0;m<*px;m++) ME(S2,k,m)=ME(S2,k,m)+VE(Di,k)*VE(Di,m)*rralpha; 
  } }
} /* bidrag fra uobserverede */  

   scl_vec_mult(1/S0,S1,S1); scl_mat_mult(1/S0,S2,S2); 

   for(k=0;k<*px;k++) for(i=0;i<*px;i++) ME(E2,k,i)=VE(S1,k)*VE(S1,i);
   vec_subtr(xi[s],S1,Ui); vec_add(U,Ui,U); 
   mat_subtr(S2,E2,TMP); mat_add(dU,TMP,dU); 

   for(k=0;k<*px;k++) 
   for(i=0;i<*px;i++) ME(dU2,k,i)=ME(dU2,k,i)-VE(Ui,k)*VE(Ui,i);

   breslow[0*(*nbre)+count]=T[s];
   for (j=1;j<=(*antbase);j++) {
   if (count==0) breslow[j*(*nbre)]=(1/S0)*(j==basestratum[s]);
   else breslow[j*(*nbre)+count]=breslow[j*(*nbre)+count-1]+(1/S0)*(j==basestratum[s]);}
   count=count+1; 
} /* if status[s]==1 */ } /* s nx*/ 
invert(dU,UI); Mv(UI,U,delta); /* print_mat(SI); print_vec(U); */ 

if (*betait>1) vec_add(beta,delta,beta); 
for(j=0;j<*px;j++) sumU=sumU+fabs(VE(U,j)); 

if (*detail>0) { 
printf(" CC==== nested CC ====NPMLE ======================== \n"); 
printf(" Nit betait  %d %d \n",oit,it); 
//m_output(dU); 
print_mat(UI); print_mat(dU);  
printf(" score er "); print_vec(U); 
printf(" beta er "); print_vec(beta); };

if ((sumU<0.000000001) & (*betait>1) ) break;  
} /* betait */
 
for (i=0;i<*antstrat;i++) for (j=0;j<*antW;j++) alphaj[i*(*antW)+j]=0; 

for (i=0;i<*nno;i++) { dummy=0; 
for (j=0;j<*antW;j++) {
if (juul==0) { bases=basestratum[nobs+i]; Zuobeta=Wbeta[j]; } 
if (juul==1) { bases=j+1; Zuobeta=xbhat[nobs+i]; } 

VE(alphaij[i],j)=exp(-breslow[bases*(*nbre)+timeind[i]]*
    exp(Zuobeta+zbhatz[i]))*p[(stratum[nobs+i]-1)*(*antW)+j]; 
dummy=dummy+VE(alphaij[i],j);}
for (j=0;j<*antW;j++) {VE(alphaij[i],j)=VE(alphaij[i],j)/dummy;
/* alpharet[j*(*nno)+i]=VE(alphaij[i],j);  */
alphaj[(stratum[nobs+i]-1)*(*antW)+j] 
=alphaj[(stratum[nobs+i]-1)*(*antW)+j]+VE(alphaij[i],j);} }  

for (i=0;i<*antstrat;i++) for (j=0;j<*antW;j++) 
p[i*(*antW)+j]=(alphaj[i*(*antW)+j]+oWj[i*(*antW)+j])/nstrat[i]; 

vec_subtr(betagam,beta,delta); 
for(j=0;j<*px;j++) sumU=sumU+fabs(VE(delta,j)); 
for (i=0;i<*antstrat;i++) 
for(j=0;j<*antW;j++) sumP=sumP+fabs(gamp[i*(*antW)+j]-p[i*(*antW)+j]); 

for (i=0;i<*antbase;i++) for(j=0;j<*nbre;j++) 
sumB=sumB+fabs(breslow[(i+1)*(*nbre)+j]-gambreslow[i*(*nbre)+j]); 
 
if (*betait>1) sumtot=sumU+sumP+sumB; else sumtot=sumP+sumB;
konv[0]=sumU; konv[1]=sumP; konv[2]=sumB; 

if (*detail>2) printf("sumU sumP sumB %lf %lf %lf \n",sumU,sumP,sumB); 

for(j=0;j<*px;j++) VE(betagam,j)=VE(beta,j);
for (i=0;i<*antstrat;i++) for(j=0;j<*antW;j++) gamp[i*(*antW)+j]=p[i*(*antW)+j]; 
for (i=0;i<*antbase;i++) 
for(j=0;j<*nbre;j++) gambreslow[i*(*nbre)+j]=breslow[(i+1)*(*nbre)+j];

if (sumtot<0.00000001) break; 
} /* Nit */

for(j=0;j<*px;j++) {
betaS[j]= VE(beta,j); loglike[0]=ll; delt[j]= VE(delta,j); score[j]=VE(U,j); 
 for (k=0;k<*px;k++) { Uii[k*(*px)+j]=-ME(UI,j,k); }  
}

free_mats(&E2,&S2,&UI,&UI2,&dU,&TMP,&dU2,NULL);
free_vecs(&xisuz,&ta,&betagam,&zis,&betaz,&xis,&U,&beta,&pi,&Yi,&Di,&Wpi,&S1,&Ui,&delta,NULL); 
for(j=0;j<*nno;j++) {free_vec(zi[j]); free_vec(alphaij[j]); }
for(j=0;j<*nx;j++) free_vec(xi[j]); for(j=0;j<*antW;j++) free_vec(Wj[j]); 

 free(imin); free(timeind); free(nnstrat); free(Wbeta); free(Zbeta);
 free(xbhat); free(zbhatz); free(alphaj); free(sumalpha); free(gamp);
 free(gambreslow);
}
