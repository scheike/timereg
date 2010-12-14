#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "multicomp.h"
#include <Rdefines.h>
#include <R.h>

                 
/* ====================================================== */
void cor(times,Ntimes,x, delta,cause,CA1,
		KMc,z,antpers, px,Nit,score,
		hess,est,gamma, semi,zsem,pg,
		detail,biid,gamiid,timepow,theta,vartheta,
		thetades,ptheta,antclust, cluster,clustsize,clusterindex,
		maxclust,step,inverse,CA2,x2,px2,
		semi2,z2,pg2,est2,gamma2,b2iid,
		gam2iid, htheta,dhtheta,rhoR,dimpar,flexfunc,
		thetiid, sym, weights, notaylor
) // {{{
double *theta,*times,*x,*KMc,*z,*score,*hess,*est,*gamma,*zsem,*vartheta,*biid,*gamiid,*timepow,*thetades,*step,*x2,*z2,*est2,*gamma2,*b2iid,*gam2iid,*thetiid,*weights;
int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*semi,*pg,*CA1,*CA2,*detail,*ptheta,
*antclust,*cluster,*clustsize,*clusterindex,*maxclust,*inverse,
	*pg2,*px2,*semi2,*dimpar,*flexfunc,*sym,*notaylor;
SEXP htheta,dhtheta,rhoR; 
{
// {{{ allocation and def's
 matrix *ldesignX,*ldesignG,*X2,*Z2;
 matrix *DUeta[*Ntimes],*DUeta2[*Ntimes],*DUgamma2,*DUgamma;
 matrix *Biid[*antclust],*B2iid[*antclust]; 
 matrix *destheta,*d2UItheta,*d2Utheta,*varthetascore; //*Sthetaiid[*antclust]; 
 vector *gamma2iid[*antclust],
	*gammaiid[*antclust],*W2[*antclust],*W3[*antclust];
 vector *dB,*VdB,*bhatt,*pbhat,*plamt;
 vector *pghat0,*pghat,*gam,*pthetavec;
 vector *xk,*xi,*rowX,*rowZ,*difX,*zk,*zi;
 vector *Utheta,*vthetascore,*vtheta1,*vtheta2,*vtheta3,*dtheta;//*thetaiid[*antclust]; 
 vector *gam2,*bhatt2,*pbhat2,*pghat2,*pghat02,*rowX2,*xi2,*xk2,*zi2,*zk2,
	*rowZ2; 

 int pmax,v,itt,i,j,k,l,s,c;
 double Li,thetak,response=0,time,dtime,timem=0;
 double fabs(),diff;// Dinverse,DDinverse; 
 double sumscore,sdj,resp3=0,resp1,resp2,ormarg,
	*vtime=calloc(1,sizeof(double));
	
//if (*trans==1) for (j=0;j<*pg;j++) if (timepow[j]!= 1) {timem=1;break;}
//if (*trans==2) for (j=0;j<*pg;j++) if (timepow[j]!= 0) {timem=1;break;}

 if (*notaylor==0) 
  for (j=0;j<*Ntimes;j++) { 
    malloc_mat(*px,*dimpar,DUeta[j]); malloc_mat(*px2,*dimpar,DUeta2[j]); 
  }
  for (j=0;j<*antclust;j++) { 
  if (*notaylor==0) {
     malloc_vec(*pg,gammaiid[j]); malloc_mat(*Ntimes,*px,Biid[j]); 
     malloc_vec(*pg2,gamma2iid[j]); malloc_mat(*Ntimes,*px2,B2iid[j]); 
     malloc_vec(*dimpar,W3[j]); 
  }
     malloc_vec(*dimpar,W2[j]); 
     // malloc_vec(*dimpar,thetaiid[j]); malloc_mat(*dimpar,*dimpar,Sthetaiid[j]); 
  }

  malloc_mats(*antpers,*px,&ldesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,NULL); 
  malloc_mats(*antpers,*pg2,&Z2,NULL);
  malloc_mats(*antpers,*px2,&X2,NULL);
  malloc_mats(*pg,*dimpar,&DUgamma,NULL);
  malloc_mats(*pg2,*dimpar,&DUgamma2,NULL);
  malloc_mats(*dimpar,*dimpar,&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  malloc_mat(*antpers,*ptheta,destheta); 

  malloc_vecs(*dimpar,&Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,&vtheta3,NULL);
  malloc_vecs(*ptheta,&pthetavec,NULL);
  malloc_vecs(*antpers,&pbhat2,&pghat02,&pghat2,NULL);
  malloc_vecs(*px2,&bhatt2,&rowX2,&xi2,&xk2,NULL);
  malloc_vecs(*pg2,&gam2,&rowZ2,&zi2,&zk2,NULL);
  malloc_vecs(*px,&xk,&xi,&rowX,&difX,&dB,&VdB,&bhatt,NULL);
  malloc_vecs(*pg,&zk,&zi,&rowZ,&gam,NULL);
  malloc_vecs(*antpers,&pbhat,&pghat0,&pghat,&plamt,NULL);
  // }}}

// {{{ reading in variables 

  pmax=max(*px,*pg); pmax=max(pmax,*ptheta); 
  // int pmax2=max(*px2,*pg2); 

  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 
  for (j=0;j<*dimpar;j++) VE(vtheta1,j)=theta[j]; 

  if (*notaylor==0) {
  if (*semi==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++)
        VE(gammaiid[i],j)=gamiid[j*(*antclust)+i]; 
  /*
  if (*semi2==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg2;j++)
        VE(gamma2iid[i],j)=gam2iid[j*(*antclust)+i]; 
   */

   for (i=0;i<*antclust;i++) {
   for (s=0;s<*Ntimes;s++) 
   for (c=0;c<*px;c++) {l=i*(*px)+c; ME(Biid[i],s,c)=biid[l*(*Ntimes)+s]; }
// for (s=0;s<*Ntimes;s++) 
// for (c=0;c<*px2;c++) {l=i*(*px2)+c; ME(B2iid[i],s,c)=b2iid[l*(*Ntimes)+s];}
   }
  }


    for (c=0;c<*antpers;c++) {
      for(j=0;j<pmax;j++)  {
	if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c];
	if (j<*ptheta) ME(destheta,c,j)= thetades[j*(*antpers)+c];
        if (*semi==1) if (j<*pg) { ME(ldesignG,c,j)=zsem[j*(*antpers)+c];}} 
        //if (*CA1!=*CA2) {
        //  for(j=0;j<pmax2;j++)  {
        //    if (j<*px2) ME(X2,c,j)=x2[j*(*antpers)+c];
        //   if (*semi2==1) if (j<*pg2) { ME(Z2,c,j)=z2[j*(*antpers)+c];}}
        //  for (j=0;j<*pg2;j++) VE(gam2,j)=gamma2[j]; 
        //} 
    }

   if (*semi==1) Mv(ldesignG,gam,pghat0);
   if (*CA1!=*CA2 && *semi2==1 && 3==2) Mv(Z2,gam2,pghat02);

// }}}

  for (itt=0;itt<*Nit;itt++)
  {
    R_CheckUserInterrupt();
    sumscore=0; mat_zeros(d2Utheta); vec_zeros(Utheta); 

      for (s=0;s<*Ntimes;s++)
      {
	  time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 
	  vtime[0]=time; 
	  for(j=1;j<=*px;j++) VE(bhatt,j-1)=est[j*(*Ntimes)+s];
	  Mv(ldesignX,bhatt,pbhat); 
	  if (*semi==1) {scl_vec_mult(time,pghat0,pghat);vec_add(pbhat,pghat,pbhat);}

          if (*CA1!=*CA2 && 3==2) { // {{{
	     for(j=1;j<=*px2;j++) {VE(bhatt2,j-1)=est2[j*(*Ntimes)+s];}
	     Mv(X2,bhatt2,pbhat2); 
	     if (*semi2==1) {scl_vec_mult(time,pghat02,pghat2);vec_add(pbhat2,pghat2,pbhat2);}
	  } // }}}

    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) { 
         // if (itt==0) printf(" %d \n",j); 

          vec_zeros(vtheta2);diff=0;sdj=0; vec_zeros(rowX);vec_zeros(rowZ); 

          for (c=0;c<clustsize[j];c++) for (v=0;v<clustsize[j];v++) 
	  if ((*sym==1 && c!=v) || (*sym==0 && c<v)) { // {{{
	    i=clusterindex[c*(*antclust)+j]; k=clusterindex[v*(*antclust)+j];
	    resp1= ((x[k]<=time) && (cause[k]==*CA2));
	    resp2= ((x[i]<=time) && (cause[i]==*CA1))*
		   ((x[k]<=time) && (cause[k]==*CA2));
         
	   if (c==0 && v==1)  extract_row(destheta,clusterindex[j],pthetavec); 

	    if (*flexfunc==0) thetak=vec_prod(pthetavec,vtheta1); 
	    else { thetak=evalh(vtheta1,vtime,pthetavec,htheta,rhoR); }
         Li=VE(pbhat,i); ormarg=(1-exp(-Li))/(exp(-Li)); 

	 if (KMc[i]<0.001) { resp2=resp2/0.001; } else { resp2=resp2/KMc[i]; }
	 if (KMc[k]<0.001) { resp1=resp1/0.001; resp2=resp2/0.001; }
	 else { resp2=resp2/KMc[k]; resp1=resp1/KMc[k];}

	if (*inverse==0) {
	   response=exp(thetak)*ormarg*(resp1-resp2)-resp2; 
	   sdj=sdj+exp(thetak)*ormarg*(resp1-resp2);
	   diff=diff+response; 
	   resp3=exp(thetak)*(resp1-resp2)*exp(Li);
	}

	if (*inverse==1) {
	   response=exp(thetak)*(exp(thetak)*ormarg*(resp1-resp2)-resp2); 
	   sdj=sdj+2*exp(2*thetak)*ormarg*(resp1-resp2)-exp(thetak)*resp2;
	   diff=diff+response; 
	   resp3=exp(2*thetak)*(resp1-resp2)*exp(Li);
	}


if (isnan(response))   { // {{{ print diverse
printf(" %d %d %d %d %d \n",s,j,clustsize[j],i,k);  
printf(" %lf %lf %lf %lf %lf %lf \n",exp(thetak),resp1,resp2,resp3,ormarg,response); 
printf(" %d %d %d \n",clustsize[j],i,k); 
printf(" %lf %lf %lf \n",x[i],x[k],time); 
printf(" resp %lf  %lf %lf  \n",response,KMc[i],KMc[k]); 
printf(" %d %d  \n",i,k); 

printf("============================== \n"); 
} // }}}

        if (*notaylor==0) 
	if (itt==*Nit-1) { // {{{
	  extract_row(ldesignX,i,xi); 
	  scl_vec_mult(resp3,xi,xi); vec_add(xi,rowX,rowX); 

	  if (*semi==1)  { 
	     extract_row(ldesignG,i,zi); 
	     scl_vec_mult(resp3*time,zi,zi); vec_add(zi,rowZ,rowZ); 
	  }

        } // }}}

        } /* for (c=0....... */   // }}}

	
        if (*flexfunc==0) scl_vec_mult(1,pthetavec,vtheta2);  
        else { evaldh(vtheta1,vtime,pthetavec,vtheta2,dhtheta,rhoR);}

       for (k=0;k<*dimpar;k++) for (c=0;c<*dimpar;c++) 
	 ME(d2Utheta,k,c)= ME(d2Utheta,k,c)+sdj*weights[j]*VE(vtheta2,k)*VE(vtheta2,c);

       scl_vec_mult(diff,vtheta2,vthetascore);  
       scl_vec_mult(weights[j],vthetascore,vthetascore); 
       vec_add(vthetascore,Utheta,Utheta); 


      if (itt==*Nit-1) { // {{{
         if (*notaylor==0) {
         for (k=0;k<*px;k++) for (c=0;c<*dimpar;c++) 
         ME(DUeta[s],k,c)+= VE(rowX,k)*VE(vtheta2,c);  
         if (*semi==1) for (k=0;k<*pg;k++) for (c=0;c<*dimpar;c++) 
         ME(DUgamma,k,c)+= VE(rowZ,k)*VE(vtheta2,c);  
	 }
         vec_add(vthetascore,W2[j],W2[j]);
//          for (k=0;k<*dimpar;k++) for (c=0;c<*dimpar;c++) 
//          ME(Sthetaiid[j],k,c)+= sdj*VE(vtheta2,k)*VE(vtheta2,c); 
      } // }}}

   } /* j in antclust */ 

   if (*notaylor==0)
   if (itt==*Nit-1) for (j=0;j<*antclust;j++) {
      extract_row(Biid[j],s,rowX); 
      vM(DUeta[s],rowX,dtheta); vec_add(dtheta,W3[j],W3[j]); }

    } /* s=1,...Ntimes */


//LevenbergMarquardt(d2Utheta,d2UItheta,Utheta,dtheta,step,step);
invert(d2Utheta,d2UItheta); Mv(d2UItheta,Utheta,dtheta);
scl_vec_mult(step[0],dtheta,dtheta); 

if (*detail==1) { // {{{
printf("===============Iteration %d ==================== \n",itt);

printf("Theta \n"); print_vec(vtheta1);
printf("delta theta \n"); print_vec(dtheta);
printf("Score D l\n"); print_vec(Utheta);
printf("Information D^2 l\n"); print_mat(d2UItheta); 
} // }}}

for (k=0;k<*dimpar;k++) sumscore= sumscore+fabs(VE(Utheta,k));
if ((sumscore<0.00001) & (itt<*Nit-2)) itt=*Nit-2;
if (isnan(vec_sum(Utheta)))  itt=*Nit-1;

vec_subtr(vtheta1,dtheta,vtheta1);
} /*itt løkke */ 


//print_mat(DUgamma); 

vec_zeros(dtheta); 
for (j=0;j<*antclust;j++) 
{
//vec_add(dtheta,W2[j],dtheta); 
// printf("===============lomse lomse ==================== \n",j);
// print_vec(W2[j]); print_vec(W3[j]); 
 if (*notaylor==0) {
    vec_add(W2[j],W3[j],W2[j]); 
     if (*semi==1) { //printf(" =W4======== \n"); print_vec(W4[j]); 
         vM(DUgamma,gammaiid[j],vtheta3); vec_add(W2[j],vtheta3,W2[j]);
 // printf(" semi 1 %d \n",j); print_vec(dtheta); 
      }
  }

for (k=0;k<*dimpar;k++) for (c=0;c<*dimpar;c++) 
ME(varthetascore,k,c)=ME(varthetascore,k,c)+VE(W2[j],c)*VE(W2[j],k); 

Mv(d2UItheta,W2[j],vtheta2);
// print_mat(d2UItheta); print_vec(W2[j]); print_vec(vtheta2); 
for (k=0;k<*dimpar;k++) thetiid[k*(*antclust)+j]=VE(vtheta2,k);
}

//  printf("W2 sum ==== \n"); print_vec(dtheta); 
//  printf("second derivative \n"); print_mat(d2UItheta); 
//  printf("iid decomp variance \n"); print_mat(varthetascore); 

   MxA(varthetascore,d2UItheta,d2Utheta); 
   MxA(d2UItheta,d2Utheta,varthetascore);

  for (j=0;j<*dimpar;j++) {theta[j]=VE(vtheta1,j); score[j]=VE(Utheta,j);
    for (k=0;k<*dimpar;k++) {vartheta[k*(*dimpar)+j]=ME(varthetascore,j,k);
                             hess[k*(*dimpar)+j]=ME(d2UItheta,j,k);}}

   // {{{ freeing 
  free(vtime); 

  if (*notaylor==0) 
  for (j=0;j<*Ntimes;j++) { free_mat(DUeta[j]); free_mat(DUeta2[j]); }

  for (j=0;j<*antclust;j++) { 
	  if (*notaylor==0) {
               free_vec(gammaiid[j]); free_mat(Biid[j]); 
               free_vec(gamma2iid[j]); free_vec(W3[j]);free_mat(B2iid[j]); 
	  }
     free_vec(W2[j]);   // free_vec(W4[j]); 
     // free_vec(thetaiid[j]);  free_mat(Sthetaiid[j]); 
  }

  free_mats(&ldesignX,&ldesignG,&Z2,&X2,&DUgamma,&DUgamma2,NULL);
  free_mats(&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  free_mats(&destheta,NULL); 

  free_vecs(&pthetavec,&pbhat2,&pghat02,&pghat2,
  &Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,&vtheta3,
  &bhatt2,&rowX2,&xi2,&xk2,&gam2,&rowZ2,&zi2,&zk2,NULL); 


   free_vecs(&xk,&xi,&rowX,&difX,&dB,&VdB,&bhatt,
  	   &zk,&zi,&rowZ,&gam,&pbhat,&pghat0,&pghat,&plamt,NULL);

//
  // }}}
} // }}}


SEXP mkansvL(double *x,int dim) 
// {{{
{ 
   int i;
   SEXP ans;
   PROTECT(ans = allocVector(REALSXP, dim));
   for(i=0;i<dim;i++) {REAL(ans)[i] = x[i]; }
   UNPROTECT(1);
   return ans;
} 
// }}}

double evalh(vector *theta,double *t,vector *xih,SEXP f,SEXP rho) 
	// {{{
{
   int j,dimtheta,dimxih; 
   dimtheta=length_vector(theta); 
   dimxih=length_vector(xih); 
   double thet[dimtheta],xx[dimxih]; 
   for (j=0;j<dimxih;j++) xx[j]=VE(xih,j);
   for (j=0;j<dimtheta;j++) thet[j]=VE(theta,j); ; 
   double res[1]; 
   void funcevalh(); 
//   printf(" ============================ \n"); printf(" %lf \n",t[0]); 
//  print_vec(theta); print_vec(xih); 
   funcevalh(thet,t,xx, f ,rho,res, dimtheta, dimxih);
   return(res[0]); 
} 
// }}}

void funcevalh(double *theta,double *t,double *x, SEXP f, SEXP rho,double *res,int dimtheta,int dimx)
// {{{
{ 
   SEXP ans;
   defineVar(install("theta"), mkansvL(theta,dimtheta),rho);
   defineVar(install("t"), mkansvL(t,1),rho);
   defineVar(install("x"), mkansvL(x,dimx),rho);
   ans=eval(f,rho);
   PROTECT(ans);
     res[0]=REAL(ans)[0]; 
   UNPROTECT(1);
} // }}}

void evaldh(vector *theta,double *t,vector *xih,vector *dh,SEXP f,SEXP rho) 
 // {{{
{
   int i,j,dimx,dimxih; 
   dimx=length_vector(theta); 
   dimxih=length_vector(xih); 
   // double x[dimx]; 
   double res[dimx]; 
   double thet[dimx],xx[dimxih]; 
   void funcevaldh(); 
   for (j=0;j<dimxih;j++) xx[j]=VE(xih,j);
   for (j=0;j<dimx;j++) thet[j]=VE(theta,j); ; 
   for (j=0;j<dimx;j++) res[j]=0; 
 
   // for (j=0;j<dimx;j++) x[j]=VE(xi,j);
   // for (j=0;j<dimxih;j++) res[j]=0; 
   funcevaldh(thet,t,xx,res, f,rho, dimx, dimxih,dimx);
   for(i=0;i<dimx;i++) VE(dh,i)=res[i];
} 
// }}}

void funcevaldh(double *theta,double *t,double *x, double *res, SEXP f,SEXP rho,int dimtheta,int dimx,int dimres)
// {{{
{ 
   SEXP ans;
   int i;
   defineVar(install("theta"), mkansvL(theta,dimtheta),rho);
   defineVar(install("t"), mkansvL(t,1),rho);
   defineVar(install("x"), mkansvL(x,dimx),rho);
   ans=eval(f,rho);
   PROTECT(ans);
     for(i=0;i<dimres;i++) { res[i]=REAL(ans)[i]; }
   UNPROTECT(1);
} // }}}


void mcifrr(times,Ntimes,x, delta,cause,CA1,
		KMc,z,antpers, px,Nit,score,
		hess,est,gamma, semi,zsem,pg,
		detail,biid,gamiid,timepow,theta,vartheta,
		thetades,ptheta,antclust, cluster,clustsize,clusterindex,
		maxclust,step,inverse,CA2,x2,px2,
		semi2,z2,pg2,est2,gamma2,b2iid,
		gam2iid, htheta,dhtheta,rhoR,dimpar,flexfunc,
		thetiid, sym, weights, notaylor
) // {{{
double *theta,*times,*x,*KMc,*z,*score,*hess,*est,*gamma,*zsem,*vartheta,*biid,*gamiid,*timepow,*thetades,*step,*x2,*z2,*est2,*gamma2,*b2iid,*gam2iid,*thetiid,*weights;
int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*semi,*pg,*CA1,*CA2,*detail,*ptheta,
*antclust,*cluster,*clustsize,*clusterindex,*maxclust,*inverse,
	*pg2,*px2,*semi2,*dimpar,*flexfunc,*sym,*notaylor;
SEXP htheta,dhtheta,rhoR; 
{
// {{{ allocation and def's
 matrix *ldesignX,*ldesignG,*X2,*Z2;
 matrix *DUeta[*Ntimes],*DUeta2[*Ntimes],*DUgamma2,*DUgamma;
 matrix *Biid[*antclust],*B2iid[*antclust]; 
 matrix *destheta,*d2UItheta,*d2Utheta,*varthetascore; // *Sthetaiid[*antclust]; 
 vector *gamma2iid[*antclust],
	*gammaiid[*antclust],*W2[*antclust],*W3[*antclust];
 vector *dB,*VdB,*bhatt,*pbhat,*plamt;
 vector *pghat0,*pghat,*gam,*pthetavec;
 vector *xk,*xi,*rowX,*rowZ,*difX,*zk,*zi;
 vector *Utheta,*vthetascore,*vtheta1,*vtheta2,*vtheta3,*dtheta; // *thetaiid[*antclust]; 
 vector *gam2,*bhatt2,*pbhat2,*pghat2,*pghat02,*rowX2,*xi2,*xk2,*zi2,*zk2,
	*rowZ2; 

 int pmax,v,itt,i,j,k,l,s,c;
 double Li,Lk,thetak,response,time,dtime;
 double fabs(),diff;// Dinverse,DDinverse; 
 double sumscore,sdj,resp3,resp1=0,resp2,ormarg,
	*vtime=calloc(1,sizeof(double));
	
  if (*notaylor==0) 
  for (j=0;j<*Ntimes;j++) { 
	  malloc_mat(*px,*dimpar,DUeta[j]); malloc_mat(*px2,*dimpar,DUeta2[j]); 
  }
  for (j=0;j<*antclust;j++) { 
	  if (*notaylor==0) {
             malloc_vec(*pg,gammaiid[j]); malloc_mat(*Ntimes,*px,Biid[j]); 
             malloc_vec(*pg2,gamma2iid[j]); malloc_mat(*Ntimes,*px2,B2iid[j]); 
             malloc_vec(*dimpar,W3[j]); 
	  }
     malloc_vec(*dimpar,W2[j]); 
     // malloc_vec(*dimpar,thetaiid[j]); malloc_mat(*dimpar,*dimpar,Sthetaiid[j]); 
  }
  malloc_mat(*antpers,*ptheta,destheta); 
  malloc_mats(*dimpar,*dimpar,&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  malloc_mats(*antpers,*px,&ldesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,NULL); 
  malloc_mats(*antpers,*pg2,&Z2,NULL);
  malloc_mats(*antpers,*px2,&X2,NULL);
  malloc_mats(*pg,*dimpar,&DUgamma,NULL);
  malloc_mats(*pg2,*dimpar,&DUgamma2,NULL);

  malloc_vecs(*dimpar,&Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,&vtheta3,NULL);
  malloc_vecs(*ptheta,&pthetavec,NULL);
  malloc_vecs(*antpers,&pbhat2,&pghat02,&pghat2,NULL);
  malloc_vecs(*px2,&bhatt2,&rowX2,&xi2,&xk2,NULL);
  malloc_vecs(*pg2,&gam2,&rowZ2,&zi2,&zk2,NULL);
  malloc_vecs(*px,&xk,&xi,&rowX,&difX,&dB,&VdB,&bhatt,NULL);
  malloc_vecs(*pg,&zk,&zi,&rowZ,&gam,NULL);
  malloc_vecs(*antpers,&pbhat,&pghat0,&pghat,&plamt,NULL);
  // }}}

// {{{ reading in variables 

  pmax=max(*px,*pg); pmax=max(pmax,*ptheta); 
  int pmax2=max(*px2,*pg2); 

  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 
  for (j=0;j<*dimpar;j++) VE(vtheta1,j)=theta[j]; 

  if (*notaylor==0)  {
  if (*semi==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++)
        VE(gammaiid[i],j)=gamiid[j*(*antclust)+i]; 
  if (*semi2==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg2;j++)
        VE(gamma2iid[i],j)=gam2iid[j*(*antclust)+i]; 

   for (i=0;i<*antclust;i++) {
   for (s=0;s<*Ntimes;s++) 
   for (c=0;c<*px;c++) {l=i*(*px)+c; ME(Biid[i],s,c)=biid[l*(*Ntimes)+s]; }
   for (s=0;s<*Ntimes;s++) 
   for (c=0;c<*px2;c++) {l=i*(*px2)+c; ME(B2iid[i],s,c)=b2iid[l*(*Ntimes)+s]; }
   }
   }


    for (c=0;c<*antpers;c++) {
      for(j=0;j<pmax;j++)  {
	if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c];
	if (j<*ptheta) ME(destheta,c,j)= thetades[j*(*antpers)+c];
        if (*semi==1) if (j<*pg) { ME(ldesignG,c,j)=zsem[j*(*antpers)+c];}} 
        if (*CA1!=*CA2) {
          for(j=0;j<pmax2;j++)  {
            if (j<*px2) ME(X2,c,j)=x2[j*(*antpers)+c];
           if (*semi2==1) if (j<*pg2) { ME(Z2,c,j)=z2[j*(*antpers)+c];}}
          for (j=0;j<*pg2;j++) VE(gam2,j)=gamma2[j]; 
        } 
    }

   if (*semi==1) Mv(ldesignG,gam,pghat0);
   if (*CA1!=*CA2 && *semi2==1) Mv(Z2,gam2,pghat02);

// }}}

for (itt=0;itt<*Nit;itt++)
{
 R_CheckUserInterrupt();
 sumscore=0; mat_zeros(d2Utheta); vec_zeros(Utheta); 

for (s=0;s<*Ntimes;s++) // {{{
{
    time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 
    vtime[0]=time; 
    for(j=1;j<=*px;j++) VE(bhatt,j-1)=est[j*(*Ntimes)+s];
    Mv(ldesignX,bhatt,pbhat); 
    if (*semi==1) {scl_vec_mult(time,pghat0,pghat);vec_add(pbhat,pghat,pbhat);}
    if (*CA1!=*CA2) {
       for(j=1;j<=*px2;j++) {VE(bhatt2,j-1)=est2[j*(*Ntimes)+s];}
       Mv(X2,bhatt2,pbhat2); 
       if (*semi2==1) {scl_vec_mult(time,pghat02,pghat2);vec_add(pbhat2,pghat2,pbhat2);}
    }
	
    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) { // {{{
         // if (itt==0) printf(" %d \n",j); 

          vec_zeros(vtheta2);diff=0;sdj=0; vec_zeros(rowX);vec_zeros(rowZ); 
	  vec_zeros(rowX2);vec_zeros(rowZ2); 

          for (c=0;c<clustsize[j];c++) for (v=0;v<clustsize[j];v++) // {{{
	  if ((*sym==1 && c!=v) || (*sym==0 && c<v)) { 
	    i=clusterindex[c*(*antclust)+j]; k=clusterindex[v*(*antclust)+j];
	    resp2= ((x[i]<=time) && (cause[i]==*CA1))*
		   ((x[k]<=time) && (cause[k]==*CA2));
         
	   if (c==0 && v==1)  extract_row(destheta,clusterindex[j],pthetavec); 

	    if (*flexfunc==0) thetak=vec_prod(pthetavec,vtheta1); 
	    else { thetak=evalh(vtheta1,vtime,pthetavec,htheta,rhoR); }
         Li=VE(pbhat,i); 
         if (*CA1!=*CA2 ) Lk=VE(pbhat2,k); else Lk=VE(pbhat,k); 
	 ormarg=(1-exp(-Li))*(1-exp(-Lk)); 

	 if (KMc[i]<0.001) { resp2=resp2/0.001; } else { resp2=resp2/KMc[i]; }
	 if (KMc[k]<0.001) { resp1=resp1/0.001; resp2=resp2/0.001; }
	 else { resp2=resp2/KMc[k]; resp1=resp1/KMc[k];}

	response=resp2-exp(thetak)*ormarg; 
	diff=diff+response; 
	sdj=sdj-exp(thetak)*ormarg;
	resp3=-exp(thetak);

if (isnan(response))   { // {{{ print diverse
printf(" %d %d %d %d %d \n",s,j,clustsize[j],i,k);  
printf(" %lf %lf %lf %lf %lf %lf \n",exp(thetak),resp1,resp2,resp3,ormarg,response); 
printf(" %d %d %d \n",clustsize[j],i,k); 
printf(" %lf %lf %lf \n",x[i],x[k],time); 
printf(" resp %lf  %lf %lf  \n",response,KMc[i],KMc[k]); 
printf(" %d %d  \n",i,k); 

printf("============================== \n"); 
} // }}}

       if(*notaylor==0) { // {{{
	if (itt==*Nit-1) { 
	    extract_row(ldesignX,i,xi); 
	    scl_vec_mult(resp3*exp(-Li)*(1-exp(-Lk)),xi,xi); vec_add(xi,rowX,rowX); 

	if (*CA1!=*CA2 ) {
	    extract_row(X2,k,xk);
	    scl_vec_mult(resp3*exp(-Lk)*(1-exp(-Li)),xk,xk); vec_add(xk,rowX2,rowX2); 
	   } else {
	    extract_row(ldesignX,k,xk); 
	    scl_vec_mult(resp3*exp(-Lk)*(1-exp(-Li)),xk,xk); vec_add(xk,rowX,rowX); 
	   }

	if (*semi==1)  { // {{{
	   extract_row(ldesignG,i,zi); 
	   scl_vec_mult(resp3*exp(-Li)*(1-exp(-Lk))*time,zi,zi); vec_add(zi,rowZ,rowZ); 

	   if (*CA1!=*CA2 ) {
	      extract_row(Z2,k,zk);
	      scl_vec_mult(resp3*exp(-Lk)*(1-exp(-Li))*time,zk,zk); vec_add(zk,rowZ2,rowZ2); 
	   } else {
	      extract_row(ldesignG,k,zk);
	      scl_vec_mult(resp3*exp(-Lk)*(1-exp(-Li))*time,zk,zk); vec_add(zk,rowZ,rowZ); 
	   }
	    
        } 
	} // }}}


	} // }}}

        } /* for (c=0....... */   // }}}
	
        if (*flexfunc==0) scl_vec_mult(1,pthetavec,vtheta2);  
        else { evaldh(vtheta1,vtime,pthetavec,vtheta2,dhtheta,rhoR);}

       for (k=0;k<*dimpar;k++) for (c=0;c<*dimpar;c++) 
	 ME(d2Utheta,k,c)= ME(d2Utheta,k,c)+weights[j]*sdj*VE(vtheta2,k)*VE(vtheta2,c);

       scl_vec_mult(diff,vtheta2,vthetascore);  
       scl_vec_mult(weights[j],vthetascore,vthetascore);  
       vec_add(vthetascore,Utheta,Utheta); 

// scl_vec_mult(exp(thetak)*diff,vthetascore,vthetascore);  
// printf(" %lf %lf %lf \n",thetak,ormarg,diff); 
//  print_vec(Utheta); 

      if (itt==*Nit-1) { // {{{
	 if (*notaylor==0) {
         for (k=0;k<*px;k++) for (c=0;c<*dimpar;c++) 
            ME(DUeta[s],k,c)= ME(DUeta[s],k,c)+VE(rowX,k)*VE(vtheta2,c);  
         if (*semi==1) for (k=0;k<*pg;k++) for (c=0;c<*dimpar;c++) 
         ME(DUgamma,k,c)= ME(DUgamma,k,c)+VE(rowZ,k)*VE(vtheta2,c);  
      if (*CA1!=*CA2) {
         for (k=0;k<*px2;k++) for (c=0;c<*dimpar;c++) 
            ME(DUeta2[s],k,c)= ME(DUeta2[s],k,c)+VE(rowX2,k)*VE(vtheta2,c);  
         if (*semi==1) for (k=0;k<*pg2;k++) for (c=0;c<*dimpar;c++) 
         ME(DUgamma2,k,c)= ME(DUgamma2,k,c)+VE(rowZ2,k)*VE(vtheta2,c);  
      }
      }
      vec_add(vthetascore,W2[j],W2[j]);

//      for (k=0;k<*dimpar;k++) for (c=0;c<*dimpar;c++) 
//          ME(Sthetaiid[j],k,c)=ME(Sthetaiid[j],k,c)+sdj*VE(vtheta2,k)*VE(vtheta2,c); 
      } // }}}

   } /* j in antclust */  // }}}


   if (*notaylor==0) 
   if (itt==*Nit-1) for (j=0;j<*antclust;j++) {
      extract_row(Biid[j],s,rowX); 
      vM(DUeta[s],rowX,dtheta); vec_add(dtheta,W3[j],W3[j]); 

      if (*CA1!=*CA2) {
       extract_row(B2iid[j],s,rowX2); 
       vM(DUeta2[s],rowX2,dtheta); vec_add(dtheta,W3[j],W3[j]);
      }
   } 

} /* s=1,...Ntimes */ // }}}


//   printf(" %ld \n",*inverse); 
// print_mat(d2Utheta);  print_vec(Utheta); 

invert(d2Utheta,d2UItheta); Mv(d2UItheta,Utheta,dtheta);
scl_vec_mult(-1*step[0],dtheta,dtheta); 

if (*detail==1) { // {{{
printf("===============Iteration %d ==================== \n",itt);

printf("Theta \n"); print_vec(vtheta1);
printf("delta theta \n"); print_vec(dtheta);
printf("Score D l\n"); print_vec(Utheta);
printf("Information D^2 l\n"); print_mat(d2UItheta); 
} // }}}

for (k=0;k<*dimpar;k++) sumscore= sumscore+fabs(VE(Utheta,k));
if ((sumscore<0.000001) & (itt<*Nit-2)) itt=*Nit-2;
if (isnan(vec_sum(Utheta)))  itt=*Nit-1;

vec_add(vtheta1,dtheta,vtheta1);

} /*itt løkke */ 

//print_mat(DUgamma); 

for (j=0;j<*antclust;j++) 
{
//vec_add(dtheta,W2[j],dtheta); 
// printf("===============lomse lomse ==================== \n",j);
// print_vec(W2[j]); print_vec(W3[j]); 
if  (*notaylor==0) {
 vec_add(W2[j],W3[j],W2[j]); 
if (*semi==1) { //printf(" =W4======== \n"); print_vec(W4[j]); 
   vM(DUgamma,gammaiid[j],vtheta3); 
  vec_add(W2[j],vtheta3,W2[j]);
 // printf(" semi 1 %d \n",j); print_vec(dtheta); 
}

if (*CA1!=*CA2 && *semi2==1) { //printf(" =W4======== \n"); print_vec(W4[j]); 
  vM(DUgamma2,gamma2iid[j],vtheta3); 
  vec_add(W2[j],vtheta3,W2[j]);
 // printf(" semi 1 %d \n",j); print_vec(dtheta); 
}
}

for (k=0;k<*dimpar;k++) for (c=0;c<*dimpar;c++) 
ME(varthetascore,k,c)=ME(varthetascore,k,c)+VE(W2[j],c)*VE(W2[j],k); 

Mv(d2UItheta,W2[j],vtheta2);
// print_mat(d2UItheta); print_vec(W2[j]); print_vec(vtheta2); 
for (k=0;k<*dimpar;k++) thetiid[k*(*antclust)+j]=VE(vtheta2,k);

}

//  printf("W2 sum ==== \n"); print_vec(dtheta); 
//  printf("second derivative \n"); print_mat(d2UItheta); 
//  printf("iid decomp variance \n"); print_mat(varthetascore); 

   MxA(varthetascore,d2UItheta,d2Utheta); 
   MxA(d2UItheta,d2Utheta,varthetascore);

  for (j=0;j<*dimpar;j++) {theta[j]=VE(vtheta1,j); score[j]=VE(Utheta,j);
  for (k=0;k<*dimpar;k++) {vartheta[k*(*dimpar)+j]=ME(varthetascore,j,k);
                             hess[k*(*dimpar)+j]=ME(d2UItheta,j,k);}}

   // {{{ freeing 
  free(vtime); 
   
  if (*notaylor==0) 
  for (j=0;j<*Ntimes;j++) { free_mat(DUeta[j]); free_mat(DUeta2[j]); }

  for (j=0;j<*antclust;j++) { 
     if (*notaylor==0) {
        free_vec(gammaiid[j]); free_mat(Biid[j]); 
        free_vec(gamma2iid[j]); free_mat(B2iid[j]); 
        free_vec(W3[j]);  
     }
     free_vec(W2[j]); 
     // free_vec(thetaiid[j]);  free_mat(Sthetaiid[j]); 
  }

  free_mats(&ldesignX,&ldesignG,&Z2,&X2,&DUgamma,&DUgamma2,NULL);
  free_mats(&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  free_mats(&destheta,NULL); 

  free_vecs(&pthetavec,&pbhat2,&pghat02,&pghat2,
  &Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,&vtheta3,
  &bhatt2,&rowX2,&xi2,&xk2,&gam2,&rowZ2,&zi2,&zk2,NULL); 

 free_vecs(&xk,&xi,&rowX,&difX,&dB,&VdB,&bhatt,&zk,&zi,&rowZ,&gam,
           &pbhat,&pghat0,&pghat,&plamt,NULL);

  // }}}
}
