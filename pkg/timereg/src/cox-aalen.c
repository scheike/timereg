#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include <time.h>
#include <sys/types.h>
                 
void score(times,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,
betaS,Nit,cu,vcu,w,mw,loglike,Iinv,Vbeta,detail,offs,mof,sim,antsim,
rani,Rvcu,RVbeta,
test,testOBS,Ut,simUt,Uit,XligZ,aalen,nb,id,status,wscore,dNit,ratesim,score,dhatMit,gammaiid,dmgiid,
retur,robust,covariance,Vcovs,addresamp,addproc,
resample,gamiid,biid,clusters,antclust,vscore,betafixed,weights,entry,exactderiv,
timegroup,maxtimepoint,stratum,silent)
double
*designX,*designG,*times,*betaS,*start,*stop,*cu,*w,*loglike,*Vbeta,*RVbeta,*vcu,*offs,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*aalen,*score,*dhatMit,*gammaiid,*dmgiid,*Vcovs,*addproc,*gamiid,*biid,*vscore,*weights,*dNit;
int*covariance,*nx,*px,*ng,*pg,*antpers,*Ntimes,*mw,*Nit,*detail,*mof,*sim,*antsim,*rani,*XligZ,*nb,*id,*status,*wscore,*ratesim,*retur,*robust,*addresamp,*resample,*clusters,*antclust,*betafixed,*entry,*exactderiv,*timegroup,*maxtimepoint,*stratum,*silent;
{ 
  int basesim=0,mjump=0,timing=0; 
  clock_t c0,c1; 
  c0=clock();
  mjump=sim[2];   // multiple jumps in clusters, relevant for ratesim=0 simulering via cholesky simulering 
  basesim=sim[1]; // baseline is also simulated and variance estimated (can be omitted for some for models) 
  // basesim=0 no simulations but variance, basesim=1 (simul and variance), basesim=-1 (no simulations no variance) 
  
  if (*detail==2) Rprintf("Memory allocation starting %d %d %d  \n",*antpers,*antclust,*maxtimepoint); 
// {{{ setting up memory 
  matrix *X,*Z,*WX,*WZ,*cdesX,*cdesX2,*cdesX3,*CtVUCt,*A,*AI;
  matrix *Vcov,*dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *tmp2,*tmp3,*dS,*S1,*SI,*S2,*M1,*VU,*ZXAI,*VUI; 
  matrix *ZPZ,*RobVbeta,*Delta,*tmpM1,*Utt,*Delta2,*tmpM2;
//  matrix *St[*maxtimepoint],*M1M2[*Ntimes],*C[*maxtimepoint],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
//  matrix *St[*Ntimes],
//  matrix *M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*AIs[*Ntimes];
  matrix *Stg[*maxtimepoint],*Cg[*maxtimepoint]; 
  matrix *ZPX1,*ZPZ1,*ZPXo,*ZPZo; 
  vector *dA,*VdA,*MdA,*delta,*zav,*lamt,*lamtt;
  vector *xi,*zi,*U,*beta,*xtilde,*Gbeta,*zcol,*one,*difzzav;
  vector *offset,*weight,*varUthat[*maxtimepoint],*Uprofile;
//  vector *ZXdA[*Ntimes];
  vector *ta,*ahatt,*vrisk,*tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB; 
  vector *W2[*antclust],*W3[*antclust];
  matrix *W3t[*antclust],*W4t[*antclust],*W2t[*antclust],*Uti[*antclust]; 
  vector *Ui[*antclust]; 
  vector *reszpbeta,*res1dim;
  matrix *dAt; 
  int cin=0,ci=0,c,pers=0,i=0,j,k,l,s,s1,it,count,pmax,
      *imin=calloc(1,sizeof(int)),
      *cluster=calloc(*antpers,sizeof(int)),
      *strata=calloc(*antpers,sizeof(int)),
      *ipers=calloc(*Ntimes,sizeof(int)); 
  double S0,RR=1,time=0,ll,lle,llo;
  double tau,hati,random,scale,sumscore;
  double *cug=calloc((*maxtimepoint)*(*px+1),sizeof(double)),
         *timesg=calloc((*maxtimepoint),sizeof(double));
  double norm_rand();
  void GetRNGstate(),PutRNGstate();

  int stratpers=0,antstrat=stratum[1]; 
  double *S0strata=calloc(antstrat,sizeof(double)); 
  matrix *ZPZs[antstrat],*ZPXs[antstrat]; // ,*As[antstrat],*ZXs[antstrat]; 

//  for (j=0;j<*nx;j++) printf(" %d ",stratum[j+2]);
  if (*detail==1) Rprintf("antstrat %d \n",antstrat); 

  for (j=0;j<antstrat;j++) {
    S0strata[j]=0; 
    malloc_mat(*pg,*pg,ZPZs[j]); 
    malloc_mat(*px,*pg,ZPXs[j]); 
//    malloc_mat(*pg,*pg,As[j]); 
//    malloc_mat(*px,*pg,ZXs[j]); 
  }

  if (*detail==1) Rprintf("Memory allocation starting \n"); 

  /* float gasdev(),expdev(),ran1(); */ 

  GetRNGstate();  /* to use R random normals */

  if (*robust==1) { // {{{ 
    for (j=0;j<*antclust;j++) { 
     malloc_vec(*pg,W2[j]);
     if (basesim>=0) {
         malloc_mat(*maxtimepoint,*px,W3t[j]); 
         malloc_mat(*maxtimepoint,*px,W4t[j]); 
         malloc_vec(*px,W3[j]); 
     }
      malloc_mat(*maxtimepoint,*pg,W2t[j]); 
      malloc_mat(*maxtimepoint,*pg,Uti[j]);  
      malloc_vec(*pg,Ui[j]);  
    }
    for(j=0;j<*maxtimepoint;j++) malloc_vec(*pg,varUthat[j]);
  } // }}} 

  for (c=0;c<*nx;c++) cluster[id[c]]=clusters[c]; 

  if (sim[0]==1) {
    malloc_mat(*maxtimepoint,*pg,Delta2); malloc_mat(*maxtimepoint,*pg,tmpM2);  
  }
  if (basesim==1) {
    malloc_mat(*maxtimepoint,*px,Delta);  malloc_mat(*maxtimepoint,*px,tmpM1); 
  }

  malloc_mat(*maxtimepoint,*pg,Utt); 
  malloc_mats(*antpers,*px,&WX,&X,&cdesX,&cdesX2,&cdesX3,NULL); 
  malloc_mats(*antpers,*pg,&WZ,&ZP,&Z,NULL); 
  malloc_mats(*px,*px,&Vcov,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*pg,*pg,&RobVbeta,&ZPZ,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*px,*pg,&tmp3,&ZPX,&dYI,&Ct,NULL); 
  malloc_mats(*px,*pg,&ZPX1,NULL); malloc_mats(*pg,*pg,&ZPZ1,NULL); 
  malloc_mats(*px,*pg,&ZPXo,NULL); malloc_mats(*pg,*pg,&ZPZo,NULL); 
  malloc_mat(*Ntimes,*px,dAt); 

  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 
  malloc_vecs(*antpers,&weight,&lamtt,&lamt,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&xtilde,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,NULL); 
  malloc_vec(*nb,ta); 
  malloc_vec(*antpers,vrisk); 

  if (*detail==1) Rprintf("Memory allocation starting \n"); 

  for(j=0;j<*maxtimepoint;j++) { malloc_mat(*px,*pg,Cg[j]); malloc_mat(*pg,*pg,Stg[j]);}

  matrix *Cn,*M1M2n,*ZXAIn,*AIn; 
  if (basesim>=0) 
  {
     malloc_mat((*px)*(*Ntimes),*pg,Cn); 
     malloc_mat(*pg,(*px)*(*Ntimes),M1M2n); 
     malloc_mat((*px)*(*Ntimes),*px,AIn); 
  }
  malloc_mat(*pg,(*px)*(*Ntimes),ZXAIn); 

//  matrix *Uiclustert[*antclust]; 
//  matrix *Uicluster[*antclust]; 
//if (*ratesim==0 && mjump==1) {
//  for(j=0;j<*antclust;j++) { 
//	  malloc_mat((*pg)*(*maxtimepoint),*pg,Uiclustert[j]); 
//	  malloc_mat((*pg),(*pg),Uicluster[j]); 
//}
//}

  vector *ranvec,*vectmp;  
  malloc_vec(*pg,ranvec); 
  malloc_vec((*pg)*(*maxtimepoint),vectmp); 

//  for(j=0;j<*Ntimes;j++) { 
//    malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]); 
//    malloc_mat(*pg,*px,ZXAIs[j]);  malloc_vec(*px,dAt[j]); malloc_mat(*px,*pg,dYIt[j]); 
////    malloc_vec(*pg,ZXdA[j]); malloc_mat(*pg,*pg,St[j]); 
//  } 

  pmax=max(*px,*pg); ll=0; 
  for(j=0;j<*pg;j++) VE(beta,j)=betaS[j]; 
  for(j=0;j<*antpers;j++) {VE(weight,j)=1; VE(offset,j)=1;} 
  // }}}
  
  if (*detail==1) Rprintf("Memory allocation done \n"); 
  if (timing==2) { // {{{
  c1=clock();
  Rprintf ("\telapsed CPU time: setting up allocation  %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}

  R_CheckUserInterrupt();

  if (*detail==1) Rprintf("Iterations start \n"); 
  cu[0]=times[0]; 
  for (it=0;it<*Nit || (*Nit==0 && it==0);it++) // {{{ iterations start for cox-aalen model
  {
     if (it>0) {
       vec_zeros(U); mat_zeros(S1);   
       mat_zeros(A); mat_zeros(ZPZ); mat_zeros(ZPX); mat_zeros(ZX); 
       mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); 

       for (j=0;j<antstrat;j++)  {
         S0strata[j]=0; mat_zeros(ZPXs[j]); mat_zeros(ZPZs[j]); 
       } 
     }
     sumscore=0; S0=0; ci=0; 
     
     R_CheckUserInterrupt();

     for (s=1;s<*Ntimes;s++)  // {{{ going through time
     {
         time=times[s]; //  vec_zeros(lamt);
         if (it==(*Nit-1)) { time=times[s]; cu[s]=times[s]; vcu[s]=times[s]; 
                             Rvcu[timegroup[s]]=times[s]; cug[timegroup[s]]=times[s]; }
    
//   if (*detail==1) Rprintf("Starting Data reading, time %d \n",s); 
    // {{{ reading design and computing matrix products
	  if (s==1) { // {{{
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	  {
	   if (( (start[c]<time) && (stop[c]>=time)) ) {
                for(j=0;j<pmax;j++) {
                   if (j<*pg) {   
			   ME(Z,id[c],j)=designG[j*(*ng)+c]; 
			   VE(zi,j)=designG[j*(*ng)+c]; 
		   }
	           if (j<*px) {ME(X,id[c],j)=designX[j*(*nx)+c]; 
			       VE(xi,j)=designX[j*(*nx)+c]; 
		   }
		}
		VE(Gbeta,id[c])=vec_prod(zi,beta)+offs[c]; 
		RR=exp(VE(Gbeta,id[c]));
		S0+=RR*weights[c]; 
		S0strata[stratum[c+2]]+=RR*weights[c]; 
                for(j=0;j<pmax;j++) {
	        if (j<*px) {ME(WX,id[c],j)=RR*designX[j*(*nx)+c];}
	        if (j<*pg) {ME(WZ,id[c],j)=weights[c]*designG[j*(*ng)+c];} 
		}
		if (time==stop[c] && status[c]==1) {pers=id[c];stratpers=stratum[c+2];} 
		if (*mof==1) VE(offset,id[c])=offs[c];  
		if (*mw==1) VE(weight,id[c])=weights[c]; 
	    for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=VE(xi,k)*VE(xi,j)*RR*weights[c]; 
              if ((j<*pg) & (k<*px)) ME(ZX,j,k)+=VE(zi,j)*VE(xi,k)*RR*weights[c]; 
//              if ((j<*px) & (k<*px)) ME(As[stratum[c+2]],j,k)+=VE(xi,k)*VE(xi,j)*RR*weights[c]; 
//              if ((j<*pg) & (k<*px)) ME(ZXs[stratum[c+2]],j,k)+=VE(zi,j)*VE(xi,k)*RR*weights[c]; 
                 if ((j<*pg) & (k<*pg)) ME(ZPZ,j,k)+= VE(zi,j)*VE(zi,k)*weights[c]*RR; 
                 if ((j<*pg) & (k<*px)) ME(ZPX,k,j)+= VE(zi,j)*VE(xi,k)*weights[c]*RR;
                 if ((j<*pg) & (k<*pg)) ME(ZPZs[stratum[c+2]],j,k)+= VE(zi,j)*VE(zi,k)*weights[c]*RR; 
                 if ((j<*pg) & (k<*px)) ME(ZPXs[stratum[c+2]],k,j)+= VE(zi,j)*VE(xi,k)*weights[c]*RR;
	   }
           count=count+1; 
         }		 
	 }
           ci=*nx-1; 
           while ((stop[ci]<time)  & (ci>=0) )  ci=ci-1; 
	 } // }}}

    vec_zeros(rowX); vec_zeros(rowZ); 
    if (s>1)  // {{{ modifying design for next time points
    while ((stop[ci]<time)  & (ci>=0) ) {
	    VE(Gbeta,id[ci])=0; // vec_prod(zi,beta); 
            for(j=0;j<*px;j++) VE(xi,j)=designX[j*(*nx)+ci]; 
            for(j=0;j<*pg;j++) { VE(zi,j)=designG[j*(*nx)+ci]; 
		VE(Gbeta,id[ci])+=VE(zi,j)*VE(beta,j);  
	    }
	    RR=exp(VE(Gbeta,id[ci])+offs[ci]);
	    if (entry[ci]==1)  {
	         replace_row(X,id[ci],xi); 
		 replace_row(Z,id[ci],zi); 
		 scl_vec_mult(RR,xi,tmpv1);
		 replace_row(WX,id[ci],tmpv1);
		 scl_vec_mult(weights[ci],zi,tmpv2);
		 replace_row(WZ,id[ci],tmpv2); 
		 VE(weight,id[ci])=weights[ci]; 
		 if (*mof==1) VE(offset,id[ci])=offs[ci];  
	    } 
	    else { 
		    replace_row(X,id[ci],rowX);
		    replace_row(Z,id[ci],rowZ);
		    replace_row(WX,id[ci],rowX); 
		    replace_row(WZ,id[ci],rowZ);
		    VE(Gbeta,id[ci])=0; 
		    VE(weight,id[ci])=0; 
		    if (*mof==1) VE(offset,id[ci])=0;
	    }
	    S0+=entry[ci]*RR*weights[ci]; 
	    S0strata[stratum[ci+2]]+=entry[ci]*RR*weights[ci]; 
	  for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=entry[ci]*VE(xi,j)*VE(xi,k)*RR*weights[ci]; 
              if ((j<*pg) & (k<*px)) ME(ZX,j,k)+=entry[ci]*VE(zi,j)*VE(xi,k)*RR*weights[ci]; 
//              if ((j<*px) & (k<*px)) ME(As[stratum[ci+2]],j,k)+=entry[ci]*VE(xi,j)*VE(xi,k)*RR*weights[ci]; 
//              if ((j<*pg) & (k<*px)) ME(ZXs[stratum[ci+2]],j,k)+=entry[ci]*VE(zi,j)*VE(xi,k)*RR*weights[ci]; 
                 if ((j<*pg) & (k<*pg)) ME(ZPZ,j,k)+=entry[ci]*VE(zi,j)*VE(zi,k)*weights[ci]*RR;
                 if ((j<*pg) & (k<*px)) ME(ZPX,k,j)+=entry[ci]*VE(zi,j)*VE(xi,k)*weights[ci]*RR;
                 if ((j<*pg) & (k<*pg)) ME(ZPZs[stratum[ci+2]],j,k)+= entry[ci]*VE(zi,j)*VE(zi,k)*weights[ci]*RR; 
                 if ((j<*pg) & (k<*px)) ME(ZPXs[stratum[ci+2]],k,j)+= entry[ci]*VE(zi,j)*VE(xi,k)*weights[ci]*RR;
	  }
	  ci=ci-1; 
	  pers=id[ci]; 
	  stratpers=stratum[ci+2]; 
    }
    // }}}
   ipers[s]=pers;

//   if (s<3) { head_matrix(X); head_matrix(WX); head_matrix(Z); head_matrix(WZ); }
//   if (*detail==1) Rprintf(" Data read, time %d \n",s); 

   if (*detail==2) {
           Rprintf("___________ s=%d jump.time=%lf jump.person=%d \n",s,time,pers); 
	   Rprintf("Z matrix\n"); 
	   print_mat(Z); 
	   Rprintf("X matrix (at risk)\n"); 
	   print_mat(X); 
   }

   scl_mat_mult(1/S0,ZPZ,ZPZo); scl_mat_mult(1/S0,ZPX,ZPXo);
   if (stratum[0]==1 && *px>1) {
        scl_mat_mult(1/S0strata[stratpers],ZPZs[stratpers],ZPZo); 
        scl_mat_mult(1/S0strata[stratpers],ZPXs[stratpers],ZPXo);
   }


   // }}}

  if (s<0) { Rprintf("======================================================= %d \n",s);
	     print_mat(A); print_mat(ZPX); print_mat(ZX); 
	     print_mat(A);  print_mat(ZPZ); 
           }

   if (stratum[0]==0)  invertS(A,AI,*silent); 
   if (ME(AI,0,0)==0 && stratum[0]==0 && *silent==0) {
      Rprintf("additive design X'X not invertible at time (number, value): %d %lf \n",s,time); print_mat(A);
   }
   if (ME(AI,0,0)==0 && stratum[0]==0 && *silent==2) {
      Rprintf("additive design X'X not invertible at time (number, value) : %d %lf \n",s,time); print_mat(A);
      Rprintf("print only first time with non-invertible design X'X\n"); 
      silent[0]=0; 
   }

   if (stratum[0]==1)  {
    for (k=0;k<*px;k++) 
    if (fabs(ME(A,k,k))<0.000001)  ME(AI,k,k)=0; else ME(AI,k,k)=1/ME(A,k,k);
   }

    scale=VE(weight,pers); 
    extract_row(X,pers,xi); 
    scl_vec_mult(scale,xi,xi); 
    Mv(AI,xi,dA); 
    MxA(ZX,AI,ZXAI); 

//  if (*detail==3) {print_vec(xi); print_mat(A); print_mat(AI); }

   if (it==(*Nit-1)) {
	replace_row(dAt,s,dA); 
        for (j=0;j<*pg;j++) for (i=0;i<*px;i++) 
		ME(ZXAIn,j,(s-1)*(*px)+i)=ME(ZXAI,j,i); 
   }

  if (s<0) {print_mat(A); print_mat(AI); print_mat(ZX); }

  /* First derivative U and Second derivative S  */ 
  extract_row(Z,pers,zi);  
  Mv(ZX, dA, zav); 

  scl_vec_mult(scale,zi,zi); 
  vec_subtr(zi,zav,difzzav); 

  vec_add(difzzav,U,U); 

  if (it==((*Nit)-1)) 
  if (*detail==3) {Rprintf(" time %d %lf %lf Dl contribution \n",s,scale,times[s]); print_vec(difzzav);
  }

    if (s<0) {  // {{{ 
       Rprintf(" %d %d %lf %lf \n",pers,s,time,scale); 
       print_vec(xi); print_vec(dA); print_vec(zi); print_vec(zav); 
       print_vec(difzzav); print_vec(U); print_mat(A); print_mat(AI); 
  } // }}}

  if (*betafixed==0)  // {{{
  if (stratum[0]==0) 
  if ( (((*exactderiv==1) && (it==(*Nit-1)) && (*px>1))) || ((*exactderiv==2) && (*px>1)) ) 
  {
  if (*detail==3) Rprintf("Computation of second derivative \n"); 
  mat_zeros(ZPZ1); mat_zeros(ZPX1); 
  for (i=0;i<*antpers;i++)
  {
    extract_row(WX,i,xi);  // er det weight her, nej da !!
    VE(lamt,i)=vec_prod(xi,dA); 
    extract_row(Z,i,zi); 
    scl_vec_mult(VE(lamt,i),zi,rowZ); 
    replace_row(ZP,i,rowZ);
    extract_row(X,i,xi); 
     for(j=0;j<pmax;j++)  for(k=0;k<*pg;k++)   {
         if ((j<*px)) ME(ZPX1,j,k)+=VE(weight,i)*VE(lamt,i)*VE(xi,j)*VE(zi,k); 
         if ((j<*pg)) ME(ZPZ1,j,k)+=VE(weight,i)*VE(lamt,i)*VE(zi,k)*VE(zi,j); 
     }
  } 
  scl_mat_mult(1,ZPZ1,ZPZo); scl_mat_mult(1,ZPX1,ZPXo);
  } // }}}
//  printf(" gamle version \n"); print_mat(ZPZo); print_mat(ZPXo);

if (*betafixed==0)  {
//	printf("=====================\n"); print_mat(ZXAI); print_mat(ZPXo); 
   MxA(ZXAI,ZPXo,tmp2); mat_subtr(ZPZo,tmp2, dS); 
   if (*mw==1) {scale=VE(weight,pers); scl_mat_mult(scale,dS,dS); }
   mat_add(dS,S1,S1);  
   if (it==((*Nit)-1)) 
   if (*detail==4) {
	   Rprintf(" time %d %d %lf D2l contribution \n",s,stratpers,times[s]); 
	   print_mat(ZPZo); print_mat(ZPXo); print_mat(ZXAI); 
	   Rprintf("============================================ \n"); 
	   print_mat(tmp2); print_mat(dS);
	   print_mat(S1); 
//   for (j=0;j<antstrat;j++) { 
//   printf(" j %d %lf \n",j,S0strata[j]); 
//   MxA(ZXAI,ZPXs[j],tmp2); mat_subtr(ZPZs[j],tmp2, dS); 
//   if (*mw==1) {scale=VE(weight,pers); scl_mat_mult(scale,dS,dS); }
//   scale=VE(weight,pers)/S0strata[stratpers]; scl_mat_mult(scale,dS,dS); 
//   print_mat(dS); 
//     }
   }

   scl_mat_mult(1,S1,Stg[timegroup[s]]);
 }


  /* varians beregninger */ 
  if (it==((*Nit)-1)) { // {{{
    replace_row(Utt,timegroup[s],U);
    vscore[0*(*Ntimes)+s]=times[s]; 

    for (i=0;i<*px;i++) for (j=0;j<*pg;j++) ME(dM1M2,j,i)=VE(dA,i)*VE(difzzav,j);
    for (i=0;i<*pg;i++) { 
      if (*betafixed==1) 
	vscore[(i+1)*(*Ntimes)+s]=vscore[(i+1)*(*Ntimes)+s-1]+VE(difzzav,i)*VE(difzzav,i); 
      for (j=0;j<*pg;j++) ME(VU,i,j)=ME(VU,i,j)+VE(difzzav,i)*VE(difzzav,j); 
    }

    MxA(AI,ZPXo,dYI); mat_subtr(Ct,dYI,Ct); 
    scl_mat_mult(1,Ct,Cg[timegroup[s]]); 

    vec_star(dA,dA,VdA); mat_add(dM1M2,M1M2t,M1M2t); 
    
    if (basesim>=0) 
    for (j=0;j<*pg;j++) for (i=0;i<*px;i++) {
	    ME(M1M2n,j,(s-1)*(*px)+i)=ME(M1M2t,j,i); 
	    ME(Cn,(s-1)*(*px)+i,j)=ME(Ct,i,j); 
    }

    for (k=1;k<=*px;k++) {
      cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dA,k-1); 
      cug[k*(*maxtimepoint)+timegroup[s]]=cu[k*(*Ntimes)+s]; 
      vcu[k*(*Ntimes)+s]=VE(VdA,k-1)+vcu[k*(*Ntimes)+s-1];
    }
    if (*robust==1 && basesim>=0) { 
         for (j=0;j<*px;j++) for (i=0;i<*px;i++) ME(AIn,(s-1)*(*px)+j,i)=ME(AI,j,i); 
    }
  } // }}}

   } // }}} /* Ntimes */ 

  if (timing==1) { // {{{
  c1=clock();
  Rprintf ("\telapsed CPU time:  going through times %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}

/* for (k=0;k<*pg;k++) ME(S1,k,k)=ME(S1,k,k)+*ridge;  */
  if (*betafixed==0 )  {
      invertS(S1,SI,*silent); Mv(SI,U,delta); MxA(SI,VU,S2); MxA(S2,SI,VU); 
  }

      if (*detail==1) { 
        Rprintf("=============Iteration %d =============== \n",it);
	Rprintf("Estimate beta \n"); print_vec(beta); 
	Rprintf("delta beta \n"); print_vec(delta); 
	Rprintf("Score D l\n"); print_vec(U); 
	Rprintf("Information -D^2 l\n"); print_mat(SI); 
      };

      if (*betafixed==0 && (*Nit>0)) vec_add(beta,delta,beta); 

      for (k=0;k<*pg;k++) sumscore=sumscore+fabs(VE(U,k)); 
      if ((sumscore<0.0000001) & (it<(*Nit)-2)) {  it=*Nit-2; }
 } /* it */ // }}}

  if (*detail==2) Rprintf("Fitting done \n"); 

  if (timing==2) { // {{{
  c1=clock();
  Rprintf ("\telapsed CPU time: fitting done %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}

  R_CheckUserInterrupt();

  lle=0; llo=0; ci=0; 
  for (k=0;k<*pg;k++) score[k]=VE(U,k); 

  mat_zeros(A); mat_zeros(ZPZ); mat_zeros(ZPX); mat_zeros(ZX); 
  mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); 
  vec_zeros(zav); 

  if (*detail==2) Rprintf("robust==%d \n",*robust); 

  if (*robust==1)  // {{{ 
  {  
  for (s=1;s<*Ntimes;s++) // {{{ terms for robust variances 
  {
    time=times[s]; cu[s]=times[s]; vcu[s]=times[s]; 
    if (*robust==1) {
    Rvcu[timegroup[s]]=times[s]; cug[timegroup[s]]=times[s]; 
    timesg[timegroup[s]]=times[s]; Ut[timegroup[s]]=times[s]; 
    }

    R_CheckUserInterrupt();

    sumscore=0; S0=0; for (j=0;j<antstrat;j++)  S0strata[j]=0; 

    // {{{ reading design and computing matrix products
	  if (s==1) { // {{{
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	  {
	   if (( (start[c]<time) && (stop[c]>=time)) ) {
                for(j=0;j<pmax;j++) {
                   if (j<*pg) {   
			   ME(Z,id[c],j)=designG[j*(*ng)+c]; 
			   VE(zi,j)=designG[j*(*ng)+c]; 
		   }
	           if (j<*px) {ME(X,id[c],j)=designX[j*(*nx)+c]; 
			       VE(xi,j)=designX[j*(*nx)+c]; 
		   }
		}
		VE(Gbeta,id[c])=vec_prod(zi,beta)+offs[c]; 
		RR=exp(VE(Gbeta,id[c]));
		S0+=RR*weights[c]; 
	        S0strata[stratum[c+2]]+=RR*weights[c]; 
                for(j=0;j<pmax;j++) {
	        if (j<*px) {ME(WX,id[c],j)=RR*designX[j*(*nx)+c];}
	        if (j<*pg) {ME(WZ,id[c],j)=weights[c]*designG[j*(*ng)+c];} 
		}
		if (time==stop[c] && status[c]==1) {pers=id[c];stratpers=stratum[c+2];} 
		if (*mof==1) VE(offset,id[c])=offs[c];  
		if (*mw==1) VE(weight,id[c])=weights[c]; 
           count=count+1; 
         }		 
	 }
           ci=*nx-1; 
           while ((stop[ci]<time)  & (ci>=0) )  ci=ci-1; 
	 } // }}}

    vec_zeros(rowX); vec_zeros(rowZ); 
    if (s>1)  // {{{ modifying design for next time points
    while ((stop[ci]<time)  & (ci>=0) ) {
	    VE(Gbeta,id[ci])=0; // vec_prod(zi,beta); 
            for(j=0;j<*px;j++) VE(xi,j)=designX[j*(*nx)+ci]; 
            for(j=0;j<*pg;j++) { 
		VE(zi,j)=designG[j*(*nx)+ci]; 
		VE(Gbeta,id[ci])+=VE(zi,j)*VE(beta,j);  
	    }
	    RR=exp(VE(Gbeta,id[ci])+offs[ci]);
	    if (entry[ci]==1)  {
	         replace_row(X,id[ci],xi); 
		 replace_row(Z,id[ci],zi); 
//		 scl_vec_mult(RR*weights[ci],xi,tmpv1);
		 scl_vec_mult(RR,xi,tmpv1);
		 replace_row(WX,id[ci],tmpv1);
		 scl_vec_mult(weights[ci],zi,tmpv2);
		 replace_row(WZ,id[ci],tmpv2); 
		 if (*mw==1) VE(weight,id[ci])=weights[ci]; 
		 if (*mof==1) VE(offset,id[ci])=offs[ci];  
	    } 
	    else { 
		    replace_row(X,id[ci],rowX);
		    replace_row(WX,id[ci],rowX); 
		    replace_row(Z,id[ci],rowZ);
		    replace_row(WZ,id[ci],rowZ);
		    VE(Gbeta,id[ci])=0; 
		    if (*mw==1) VE(weight,id[ci])=0; 
		    if (*mof==1) VE(offset,id[ci])=0;  
	    }
	    S0+=entry[ci]*RR*weights[ci]; 
	    S0strata[stratum[ci+2]]+=entry[ci]*RR*weights[ci]; 
	  ci=ci-1; 
	  pers=id[ci]; 
	  stratpers=stratum[ci+2]; 
    }
    // }}}
   ipers[s]=pers;
    // }}} 
    
//   if (s<3) { head_matrix(X); head_matrix(WX); head_matrix(Z); head_matrix(WZ); }
//    extract_row(WX,pers,xi); 
     extract_row(dAt,s,dA); 
//    hati=vec_prod(xi,dA);    lle=lle+log(hati);

    for (j=0;j<*pg;j++) for (i=0;i<*px;i++) ME(ZXAI,j,i)=ME(ZXAIn,j,(s-1)*(*px)+i);
    if (basesim>=0) {
        for (j=0;j<*px;j++) for (i=0;i<*px;i++)  ME(AI,j,i)=ME(AIn,(s-1)*(*px)+j,i);
    }
//    print_mat(ZXAI); print_vec(dA); 

    if (*ratesim==1 || *retur>=1)
    for (i=0;i<*antpers;i++)   // {{{
    {
      cin=cluster[i]; 
      extract_row(WX,i,rowX);  // RR*xi
      extract_row(Z,i,zi); 
      extract_row(X,i,xi); 
      hati=vec_prod(rowX,dA); 

     if (*ratesim==1) {
      Mv(ZXAI,xi,tmpv2);  
      vec_subtr(zi,tmpv2,tmpv2); 
      if (*mw==1) { scl_vec_mult(VE(weight,i),tmpv2,tmpv2);}  // weight adjustment

      if (i==pers) vec_add(tmpv2,W2[cin],W2[cin]);
      if (*ratesim==1) {scl_vec_mult(hati,tmpv2,rowZ); vec_subtr(W2[cin],rowZ,W2[cin]); }

      if (basesim>=0) {
	      Mv(AI,xi,rowX); 
	      if (*mw==1) scl_vec_mult(VE(weight,i),rowX,rowX); 
	      if (i==pers) {vec_add(rowX,W3[cin],W3[cin]); }
	      llo=llo+hati;
	      if (*ratesim==1) {scl_vec_mult(hati,rowX,rowX); vec_subtr(W3[cin],rowX,W3[cin]);}
      }
    }

      if (*retur==1) dhatMit[i*(*Ntimes)+s]=1*(i==pers)-hati;
      if (*retur==2) dhatMit[i]=dhatMit[i]+1*(i==pers)-hati;
    } // }}}

    if (*ratesim==1) 
    for (j=0;j<*antclust;j++)  
    {
     replace_row(W2t[j],timegroup[s],W2[j]); 
     if (basesim>=0)  replace_row(W3t[j],timegroup[s],W3[j]);  
    } 

    if ((*ratesim==0))  // {{{  compute resampling counting process LWY style version
    {
      cin=cluster[pers]; 
      extract_row(WX,pers,rowX);  // RR*xi
      extract_row(Z,pers,zi); 
      extract_row(X,pers,xi); 

      Mv(ZXAI,xi,tmpv2);  
      vec_subtr(zi,tmpv2,tmpv2); 
      if (*mw==1) scl_vec_mult(VE(weight,pers),tmpv2,tmpv2); 
       vec_add(tmpv2,W2[cin],W2[cin]);

//      if (mjump==1) 
//      for (j=0;j<*pg;j++) for (i=0;i<*pg;i++) 
//         ME(Uicluster[cin],j,i)+=VE(tmpv2,j)*VE(tmpv2,i); 

      if (basesim>=0) {
	      Mv(AI,xi,rowX); 
	      scl_vec_mult(VE(weight,pers),rowX,rowX); 
	      vec_add(rowX,W3[cin],W3[cin]);
      }

      // distrubes the increments to the end for each process with jumps
      for (s1=timegroup[s];s1<*maxtimepoint;s1++) // {{{ 
      { 
//	      printf("W2t %d %d %d  \n",cin,s1,*maxtimepoint); 
//	      print_mat(W2t[cin]); 
         replace_row(W2t[cin],s1,W2[cin]); 
//	      printf("2 %d %d \n",cin,s1); 
      if (basesim>=0)  replace_row(W3t[cin],s1,W3[cin]);  
//         if (mjump==1) {
//	     cholesky(Uicluster[cin],tmp2); 
////             if (s1==timegroup[s])  {
////	        printf(" tmp2 %d \n",cin); 
////		print_vec(tmpv2); 
////		print_mat(Uicluster[cin]); 
////	        print_mat(tmp2); 
////		MtM(tmp2,dS); 
////	        print_mat(dS); 
////	     }
//	     for (j=0;j<*pg;j++) for (i=0;i<*pg;i++) ME(Uiclustert[cin],s1*(*pg)+j,i)=ME(tmp2,i,j); 
//	 }
      } // }}} 

    } // }}} 

    /* MG baseret varians beregning */

    if (basesim>=0)   // {{{ 
    {
    for (j=0;j<*pg;j++) for (i=0;i<*px;i++) {
         ME(M1M2t,j,i)=ME(M1M2n,j,(s-1)*(*px)+i);
         ME(Ct,i,j)= ME(Cn,(s-1)*(*px)+i,j);
    }

    MxA(Ct,VU,tmp3); MAt(tmp3,Ct,CtVUCt);
    MxA(Ct,SI,tmp3); MxA(tmp3,M1M2t,COV); 

    for (k=1;k<=*px;k++) {
      if (*betafixed==0) 
	vcu[k*(*Ntimes)+s]+=ME(CtVUCt,k-1,k-1)+2*ME(COV,k-1,k-1); 
    }
    } // }}} 

    for (k=1;k<=*pg;k++) Ut[k*(*maxtimepoint)+timegroup[s]]=ME(Utt,timegroup[s],k-1);
  }   // }}}
  } // }}}   

  if (timing==2) { // {{{
  c1=clock();
  Rprintf ("\telapsed CPU time: robust variance terms 1 %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}

  if (*detail==2) Rprintf("Robust variances 1 \n"); 

  R_CheckUserInterrupt();

  ll=lle-llo; /* likelihood beregnes */
  if (*detail==2) Rprintf("loglike is  %lf \n",ll);  

// check af score process er ok 
//  int itest=1; 
//  if (itest==1)
// for (s=0;s<*maxtimepoint;s++)   {
//  mat_zeros(S2); mat_zeros(dS); mat_zeros(SI);
//  mat_zeros(VUI);  mat_zeros(ZPZ); 
// for (j=0;j<*antclust;j++) 
// {
//	  extract_row(W2t[j],s,tmpv2); 
//	  for (k=0;k<*pg;k++) for (i=0;i<*pg;i++) ME(S2,k,i)+=VE(tmpv2,k)*VE(tmpv2,i); 
//	  if (mjump==1) {
//	  for (k=0;k<*pg;k++) for (i=0;i<*pg;i++) ME(dS,k,i)=ME(Uiclustert[j],s*(*pg)+k,i);
//	  mat_transp(dS,dS); 
//	  MtM(dS,VUI); 
//	  mat_add(VUI,ZPZ,ZPZ); 
//	  }
//}
//printf("score process variance  %d \n",s); 
//print_mat(S2); 
//print_mat(ZPZ); 
//}


  if ((*robust==1)) // {{{ robust variances 
  {
      for (s=1;s<*maxtimepoint;s++) 
      {
	vec_zeros(VdB); mat_zeros(Vcov);

	for (j=0;j<*antclust;j++)  // {{{ 
	{
		if (s==1 && *detail==4) {
                    Rprintf("========================= %d \n",j); 
			print_mat(W2t[j]); 
			print_vec(W2[j]); 
			print_mat(Stg[s]); 
			print_mat(S1); 
			print_mat(SI); 
		}

	  Mv(SI,W2[j],tmpv2); 
	  if (basesim>=0) {
		  Mv(Cg[s],tmpv2,rowX);
		  extract_row(W3t[j],s,tmpv1); 
		  vec_add(tmpv1,rowX,difX); 
		  if (*betafixed==1) scl_vec_mult(1,tmpv1,difX); 
		  replace_row(W4t[j],s,difX); 
		  vec_star(difX,difX,tmpv1); 
		  vec_add(tmpv1,VdB,VdB);
	  }

          if (s==1) 
          if (*betafixed==0) {
             for (c=0;c<*pg;c++) gamiid[c*(*antclust)+j]=gamiid[c*(*antclust)+j]+VE(tmpv2,c); 
 	  }
	  if (*resample==1 && basesim>=0) {
	    for (c=0;c<*px;c++) {l=j*(*px)+c; 
	      biid[l*(*maxtimepoint)+s]=biid[l*(*maxtimepoint)+s]+VE(difX,c);
	    } 
	  }
	  if (*covariance==1 && basesim>=0) {
	    for (k=0;k<*px;k++) for (c=0;c<*px;c++) 
	      ME(Vcov,k,c)=ME(Vcov,k,c)+VE(difX,k)*VE(difX,c);
	  }

	  Mv(Stg[s],tmpv2,rowZ); 
	  extract_row(W2t[j],s,tmpv2); 
          if (*detail==4) Rprintf("j,s is %d %d \n",j,s);  
	  if (*betafixed==0) {
	    vec_subtr(tmpv2,rowZ,zi); 
	    replace_row(Uti[j],s,zi); 
//	    if (mjump==1 && *ratesim==0) 
//	    {
//               cholesky(Uicluster[j],tmp2); 
//	       mat_transp(tmp2,tmp2); 
//               MxA(SI,tmp2,tmp2); 
//	       MxA(Stg[s],tmp2,dS); 
//	       for (c=0;c<*pg;c++) for (i=0;i<*pg;i++) ME(Uiclustert[j],s*(*pg)+c,i)-=ME(dS,c,i);
//	    }
	  } else replace_row(Uti[j],s,tmpv2);

	  vec_star(zi,zi,tmpv2); 
	  vec_add(tmpv2,varUthat[s],varUthat[s]);
	} // }}} /* j in clusters  */

	if (*betafixed==0) 
	  for (i=0;i<*pg;i++) vscore[(i+1)*(*maxtimepoint)+s]=VE(varUthat[s],i);  

	if (basesim>=0) 
	for (k=1;k<*px+1;k++) { Rvcu[k*(*maxtimepoint)+s]=VE(VdB,k-1);
	  if (*covariance==1) {
	    for (j=0;j<*px;j++)  { l=(k-1)*(*px)+j; 
		    Vcovs[l*(*maxtimepoint)+s]=ME(Vcov,k-1,j); 
	    } 
	  } 
	}
      }  /*  s=1 ..maxtimepoints */ 

    }  // }}} robust variance


//if (mjump==1 && *ratesim==0) 
//for (j=0;j<*antclust;j++) { 
//     printf("observed score ========================== %d \n",j); 
//     print_vec(W2[j]); 
//     print_mat(Uicluster[j]); 
//}


//if (mjump==2 && *ratesim==0) 
//for (j=0;j<*antclust;j++) { 
//     printf("observed score ========================== %d \n",j); 
//     print_mat(Uti[j]); 
//     print_mat(Uiclustert[j]); 
//}

  if (timing==2) { // {{{
  c1=clock();
  Rprintf ("\telapsed CPU time: variance terms 2 %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}

  if (*detail==2) Rprintf("Robust variances 2 \n"); 


    if ((*betafixed==0) && (*robust==1))  {
         for (j=0;j<*antclust;j++) {
	         Mv(SI,W2[j],tmpv2); 
	         for (c=0;c<*pg;c++) for (k=0;k<*pg;k++)
		 ME(RobVbeta,c,k)=ME(RobVbeta,c,k)+VE(W2[j],c)*VE(W2[j],k);
	      for (k=0;k<*pg;k++) gammaiid[j*(*pg)+k]=VE(tmpv2,k); 
	      
         }
    MxA(RobVbeta,SI,ZPZ); MxA(SI,ZPZ,RobVbeta);
    }

  if (timing==2) { // {{{
  c1=clock();
  Rprintf ("\telapsed CPU time: variance terms 3 %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}

  R_CheckUserInterrupt();

  for(j=0;j<*pg;j++) { 
    betaS[j]= VE(beta,j); 
    loglike[0]=lle; 
    loglike[1]=ll;
    for (k=0;k<*pg;k++){ 
      Iinv[k*(*pg)+j]=ME(SI,j,k);
      Vbeta[k*(*pg)+j]=-ME(VU,j,k); 
      RVbeta[k*(*pg)+j]=-ME(RobVbeta,j,k); 
    } 
  } 

//      printf("cholesky"); cholesky(SI,VU); print_mat(SI); print_mat(VU); 
//      MtM(VU,SI); print_mat(SI); 

  if (timing==2) { // {{{
  c1=clock();
  printf ("\telapsed CPU time: variance terms 4 %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}


// check af observed score process er ok , sammenligning af variancer til 
//  int itest1=1; 
//  if (itest1==1 )
// for (s=0;s<*maxtimepoint;s++)   {
//  mat_zeros(S2); mat_zeros(dS); mat_zeros(VUI); mat_zeros(ZPZ);
// for (j=0;j<*antclust;j++) 
// {
//	  extract_row(Uti[j],s,tmpv2); 
//	  for (k=0;k<*pg;k++) for (i=0;i<*pg;i++) ME(S2,k,i)+=VE(tmpv2,k)*VE(tmpv2,i); 
//	  if (mjump==1) {
//	  for (k=0;k<*pg;k++) for (i=0;i<*pg;i++) ME(dS,k,i)=ME(Uiclustert[j],s*(*pg)+k,i);
//	  print_mat(dS); 
////	  mat_transp(dS,dS); 
//	  MtM(dS,ZPZ); 
//	  printf("obs Score %d %d \n",s,j); 
//	  print_mat(ZPZ); 
//	  mat_add(ZPZ,VUI,VUI); 
//	  }
//}
//printf(" %d \n",s); 
//print_mat(S2); 
//print_mat(VUI); 
//}

//  for(j=0;j<*antclust;j++)  print_mat(Uti[j]); 
  if (*detail==2)  Rprintf("Ready for simulations sim=%d\n",sim[0]);

  if (sim[0]>=1) // {{{ score process simulations
  {
    // Rprintf("Simulations start N= %ld \n",(long int) *antsim);

    tau=times[*Ntimes-1]-times[0];
    for (i=1;i<=*px;i++) VE(rowX,i-1)=cug[i*(*maxtimepoint)+(*maxtimepoint-1)];

    for (s=1;s<*maxtimepoint;s++)  // {{{ /* Beregning af OBS teststørrelser */
    {
      time=timesg[s]-times[0];  //  FIX 

      if (basesim==1)  // {{{
      {
      for (i=1;i<=*px;i++) {
	VE(xi,i-1)=fabs(cug[i*(*maxtimepoint)+s])/sqrt(Rvcu[i*(*maxtimepoint)+s]);
	if (VE(xi,i-1)>testOBS[i-1]) testOBS[i-1]=VE(xi,i-1); 
      }

      scl_vec_mult(time/tau,rowX,difX);
      for (i=1;i<=*px;i++) VE(xi,i-1)=cug[i*(*maxtimepoint)+s];
      vec_subtr(xi,difX,difX);
      for (i=0;i<*px;i++) {
	VE(difX,i)=fabs(VE(difX,i)); l=(*px+i);
	if (VE(difX,i)>testOBS[l]) testOBS[l]=VE(difX,i);
      }
      } // }}} 

      if (*wscore>=1) {  /* sup beregnes i R */ 
	if ((s>*wscore) && (s<*maxtimepoint-*wscore)) {extract_row(Utt,s,rowZ);
	  for (i=0;i<*pg;i++) VE(rowZ,i) = VE(rowZ,i)/sqrt(VE(varUthat[s],i));
	  replace_row(Utt,s,rowZ); /* scaled score process */ 
	}
	else {vec_zeros(rowZ); replace_row(Utt,s,rowZ);}
      }
      for (k=1;k<=*pg;k++) Ut[k*(*maxtimepoint)+s]=ME(Utt,s,k-1);
    }  // }}} *s=1..maxtimepoint Beregning af obs teststørrelser 

  if (*detail==2)  Rprintf("Simulations start N= %ld \n",(long int) *antsim);

    for (k=1;k<=*antsim;k++) // {{{  k=1,...,antsim
    {
	 R_CheckUserInterrupt();
	 if (basesim==1) mat_zeros(Delta); 
	  mat_zeros(Delta2); vec_zeros(tmpv1);
     for (i=0;i<*antclust;i++)  // {{{ 
     {
	random=norm_rand();
	if (basesim==1) {
	    scl_mat_mult(random,W4t[i],tmpM1); mat_add(tmpM1,Delta,Delta);
	}
//        if ((mjump==0 && *ratesim==0) || (*ratesim==1)) {
//	   random=norm_rand();
	   scl_mat_mult(random,Uti[i],tmpM2); 
	   mat_add(tmpM2,Delta2,Delta2);
//	} else {
//           for (c=0;c<*pg;c++) VE(ranvec,c)=norm_rand(); 
//	   Mv(Uiclustert[i],ranvec,vectmp); 
//           for (c=0;c<*maxtimepoint;c++) 
//           for (l=0;l<*pg;l++) ME(tmpM2,c,l)=VE(vectmp,c*(*pg)+l); 
//	   mat_add(tmpM2,Delta2,Delta2);
//	}
     } // }}} 


    if (basesim==1) extract_row(Delta,*maxtimepoint-1,tmpv1); 

    for (s=1;s<*maxtimepoint;s++) { 

    time=timesg[s]-times[0];
    if (basesim==1)   // {{{ 
    {
	scl_vec_mult(time/tau,tmpv1,xi); extract_row(Delta,s,rowX);
	vec_subtr(rowX,xi,difX);

	if (*addresamp==1) {
	  if (k<51) { 
	    for (i=0;i<*px;i++) {l=(k-1)*(*px)+i;
	      addproc[l*(*maxtimepoint)+s]=ME(Delta,s,i);}}
	}

	for (i=0;i<*px;i++) {
	  VE(difX,i)=fabs(VE(difX,i));
	  l=(*px+i);
	  if (VE(difX,i)>test[l*(*antsim)+k-1]) test[l*(*antsim)+k-1]=VE(difX,i);
	  VE(xi,i)=fabs(ME(Delta,s,i))/sqrt(Rvcu[(i+1)*(*maxtimepoint)+s]);
	  if (VE(xi,i)>test[i*((*antsim))+k-1]) test[i*((*antsim))+k-1]=VE(xi,i); 
	}
    } // }}} 

	if (*wscore>=1) {
	  extract_row(Delta2,s,zi); 
	  if ((s>*wscore) && (s<*maxtimepoint-*wscore))  {
	    for (i=0;i<*pg;i++) {VE(zi,i)=fabs(ME(Delta2,s,i))/sqrt(VE(varUthat[s],i));
	      if (VE(zi,i)>simUt[i*(*antsim)+k-1]) simUt[i*(*antsim)+k-1]=VE(zi,i); }

	    if (k<50) { 
	      for (i=0;i<*pg;i++) { l=(k-1)*(*pg)+i;
		Uit[l*(*maxtimepoint)+s]=ME(Delta2,s,i)/sqrt(VE(varUthat[s],i));}}
	  } 
	} /* weigted score */
	else {
	  extract_row(Delta2,s,zi); 
	  for (i=0;i<*pg;i++) {
	    if (fabs(VE(zi,i))>simUt[i*(*antsim)+k-1]) simUt[i*(*antsim)+k-1]=fabs(VE(zi,i)); 
	  }

	  if (k<50) { 
	    for (i=0;i<*pg;i++) { l=(k-1)*(*pg)+i;
	      Uit[l*(*maxtimepoint)+s]=ME(Delta2,s,i);}
	  }
	} /* else wscore=0 */ 

      }  /* s=1..Ntims */
    }  // }}} /* k=1..antsim */

  } /* sim==1 */ // }}}

  if (timing==2) { // {{{
  c1=clock();
 printf ("\telapsed CPU time: before freeing %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}
  PutRNGstate();  /* to use R random normals */

  // {{{ freeing 
  if (sim[0]>=1) free_mats(&Delta2,&tmpM2,NULL); 
  if (basesim==1) free_mats(&Delta,&tmpM1,NULL); 

  if (basesim>=0) free_mats(&Cn,&M1M2n,&AIn,NULL); 
  free_mats(&ZXAIn,NULL); 

  free_mats(&dAt,&Utt,&WX,&X,&cdesX,&cdesX2,&cdesX3, &WZ,&ZP,&Z,
     &Vcov,&COV,&A,&AI,&M1,&CtVUCt, &RobVbeta,&ZPZ,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,
     &ZXAI,&ZX,&dM1M2,&M1M2t, &tmp3,&ZPX,&dYI,&Ct, &ZPX1,&ZPZ1, &ZPXo,&ZPZo,NULL); 

  free_vecs(&vectmp,&ranvec,&reszpbeta,&res1dim,&weight,&lamtt,&lamt,&zcol,&Gbeta,&one,&offset,
            &ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,
            &xtilde, &tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,
            &ta,&vrisk,NULL); 

  if (*robust==1) {
    for (j=0;j<*antclust;j++) {
    free_vec(W2[j]);
      if (basesim>=0)  {
         free_mat(W3t[j]); free_mat(W4t[j]); free_vec(W3[j]); 
      }
      free_mat(W2t[j]); free_mat(Uti[j]); free_vec(Ui[j]); 
    }
    for (j=0;j<*maxtimepoint;j++)  free_vec(varUthat[j]);
  }

//  if (*ratesim==0 && mjump==1) {
//     for(j=0;j<*antclust;j++) { free_mat(Uiclustert[j]); free_mat(Uicluster[j]); }
//  }

  for(j=0;j<*maxtimepoint;j++) { free_mat(Cg[j]); free_mat(Stg[j]);}
  free(cluster); free(ipers); free(imin); free(cug); free(timesg); 
  free(S0strata); free(strata); 
  // }}}
 for (j=0;j<antstrat;j++) { free_mat(ZPZs[j]); free_mat(ZPXs[j]); 
//                            free_mat(As[j]); free_mat(ZXs[j]); 
 }
 
  if (timing==2) { // {{{
  c1=clock();
  printf ("\telapsed CPU time: after freeing %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
  c0=clock();
  } // }}}

}
