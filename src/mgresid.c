//#include <stdio.h>
#include <math.h>
#include <R.h>
#include "matrix.h"
#include"R_ext/Random.h"

void mgresid(double *designX,int *nx,int *px,int *antpers,double *start,double *stop,int *status,int *id,
double *mgtimes,int *nmgt,double *dmgresid,int *sim,double *xval,int *ant,
double *univarproc,double *timeproc,double *simunivarproc,double *simtimeproc,
double *unitest,double *unitestOBS,
double *timetest,double *timetestOBS,
double *unitimetest,double *unitimetestOBS,
double *modelmatrix,int *model,int *pm,double *cummgt,double *dNit,double *robvarcum,
double *testOBS,double *test,double *simUt,double *Ut,int *cumresid,int *maxval,int *startdesign,
int *coxaalen,double *dcum,double *beta,double *designG,int *pg,double *Ogammaiid,
int *clusters,int *antclust,double *robvarcumz,double *simcumz,
int *inXZ,int *inXorZ,int *iptot,int *entry,int *stratum,int *silent,double *weights,double *offsets,int *ratesim,double *weightsmg, int *varweighted) 
//double *designG,*dcum,*beta,*designX,*start,*stop,*mgtimes,
//	*dmgresid,*xval,*univarproc,*timeproc,*simunivarproc,
//	*simtimeproc,*unitest,*unitestOBS, *timetest,*timetestOBS,
//	*unitimetest,*unitimetestOBS,*modelmatrix,*Ogammaiid,
//	*cummgt,*robvarcum,*testOBS,*test,*simUt,*Ut,
//	*robvarcumz,*simcumz,*weights,*offsets,*weightsmg,*dNit;
//int *pg,*coxaalen,*nx,*px,*antpers,*nmgt,*sim,*ant,
//    *status,*id,*model,*pm,*cumresid,*maxval,*startdesign,*clusters,*antclust,
//    *inXZ,*inXorZ,*iptot,*entry,*stratum,*silent,*ratesim,*varweighted; 
{ // {{{
// {{{ // memory allocation
  matrix *Delta,*tmpM1,*X,*cummat,*modelMGT[*antclust],*modMGz,*modMGzosdt;
  matrix *A,*AI,*cumX,*cumXAI,*cumZP,*XPZ,*tmp2,*dS,*S,*St[*nmgt]; 
  matrix *Z,*dS1,*S1,*cumX1,*cumXAI1,*cumZP1,*tmp21,*cummat1;
  vector *Deltazsd,*Deltaz,*tmpM1z,*vtmp2,*vtmp1,*cumdB1,*VdB1,*respm1;
  vector *dMGt[*nmgt],*cumdB,*dB,*VdB,*xi,*rowX,*rowcum,*difX,*vtmp,*respm,*gamma;
  vector *risk,*cumA[*antclust],*cum,*vecX;
  vector *Gbeta,*dA,*xtilde,*zi,*gammaiid[*antclust];
  vector *tmpv1,*rowZ,*rvec,*dB1[*antclust];
  vector *dBgam[*antclust];
  vector *weightmg,*weight,*offset; 
  int ci=0,pmax,m,i,j,k,l,s,c=0,s1=0,count,pers=0;
  int ptot,weighted,*cluster=calloc(*antpers,sizeof(int));
  double lamti=1,time,RR=1,vardiv;
  double random,xij,dtime;
//  void smoothB(),comptest(); 
//  void GetRNGstate(),PutRNGstate();

  weighted=*varweighted; ptot=*px+*pg; ptot=*iptot; 

  GetRNGstate();  /* to use R random normals */

  for (s=0;s<*nmgt;s++) malloc_vec(*antpers,dMGt[s]); 
//  for (s=0;s<*nmgt;s++) malloc_vec(*antpers,dNt[s]); 
  for (i=0;i<*antclust;i++) { malloc_mat(*nmgt,*pm,modelMGT[i]);
     malloc_vec(*pg,gammaiid[i]); malloc_vec(*pm,cumA[i]); 
  }
  malloc_mats(*nmgt,*pm,&Delta,&tmpM1,NULL); 
  malloc_mat(*antpers,*px,X); 
  malloc_mat(*antpers,*pg,Z); 
  malloc_mat(*antpers,*pm,cummat); 
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*pm,*px,&cumX,&cumXAI,NULL); 
  malloc_mat(*pm,*pg,cumZP); 
  malloc_mat(*px,*pg,XPZ); 

//  matrix *AIXZ; 
//  malloc_mat(*px,*pg,&AIXZ); 

  malloc_mats(*pm,*pg,&tmp2,&dS,&S,NULL); 
  for (s=0;s<*nmgt;s++) malloc_mat(*pm,*pg,St[s]); 

  malloc_vecs(*pm,&vtmp,&cumdB,&dB,&VdB,&respm,NULL);
  malloc_vecs(*px,&tmpv1,&cum,&dA,&xtilde,&xi,&rowX,&rowcum,&difX,NULL);
  malloc_vecs(*pg,&zi,&gamma,&rowZ,NULL); 
  malloc_vecs(*antpers,&weightmg,&weight,&offset,&vecX,&risk,&Gbeta,NULL); 
  malloc_vecs(*antclust,&rvec,NULL); 
  // }}}

  for (j=0;j<*nx;j++) { 
      m=id[j]; 
      cluster[m]=clusters[j]; 
 }
  if (*coxaalen==1) for(j=0;j<*pg;j++) VE(gamma,j)=beta[j]; 

  if (*coxaalen==1) for (i=0;i<*antclust;i++)  {
     for (j=0;j<*pg;j++) VE(gammaiid[i],j)=Ogammaiid[j*(*antclust)+i]; 
  }
  for(j=0;j<*antpers;j++) {VE(weight,j)=1; VE(offset,j)=1;} 

//     if (*coxaalen==1 && *ratesim==0) 
//     for (i=0;i<*antpers;i++) VE(dNt[s],i)=dNit[i*(*nmgt)+s];

  R_CheckUserInterrupt();
  /*  cumulative martingales Aalen type */ 
  if (*model==1)  // {{{
  {
      pmax=*px; 
      if (*coxaalen==1)  pmax=max(*px,*pg);
      pmax=max(pmax,*pm); 

      for (s=1;s<*nmgt;s++) // {{{ 
      {
       time=mgtimes[s]; dtime=mgtimes[s]-mgtimes[s-1]; 
       R_CheckUserInterrupt();

//     mat_zeros(X);mat_zeros(cummat);vec_zeros(risk);mat_zeros(Z);

       if (s>=2) {
	    mat_zeros(X);mat_zeros(Z); mat_zeros(A); 
	    mat_zeros(cummat);  mat_zeros(cumX); mat_zeros(XPZ); mat_zeros(cumZP); 
	    vec_zeros(weightmg); vec_zeros(weight); vec_zeros(offset); 
       }

   if (*coxaalen==1) for (j=0;j<*px;j++) VE(dA,j)=dcum[j*(*nmgt-1)+s-1]; 

    // {{{ reading design and computing matrix products
	  if (s>=1) { // {{{
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	  {
	   if ((start[c]<time) && (stop[c]>=time)) 
	   {
                for(j=0;j<pmax;j++) {
                   if (*coxaalen==1)
                   if (j<*pg) {   
			   VE(zi,j)=designG[j*(*nx)+c]; 
		   }
	           if (j<*px) {ME(X,id[c],j)=designX[j*(*nx)+c]; 
			       VE(xi,j)=designX[j*(*nx)+c]; 
		   }
	           if (j<*pm) {ME(cummat,id[c],j)=modelmatrix[j*(*nx)+c];
                               VE(vtmp,j)=modelmatrix[j*(*nx)+c]; 
		   }
		}
		if (*coxaalen==1) { VE(Gbeta,id[c])=vec_prod(zi,gamma); 
		                    RR=exp(VE(Gbeta,id[c])+offsets[c]);
	                            lamti=RR*vec_prod(xi,dA);
		 } 
		VE(offset,id[c])=offsets[c];  
		VE(weight,id[c])=weights[c]; 
		VE(weightmg,id[c])=weightsmg[c]; 
		if (time==stop[c] && status[c]==1) pers=id[c];

	    for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=weights[c]*VE(xi,k)*VE(xi,j)*RR; 
              if ((j<*px) & (k<*pm)) ME(cumX,k,j)+=weightsmg[c]*VE(xi,j)*VE(vtmp,k)*RR;
	      if (*coxaalen==1) {
               if ((j<*pg)&(k<*px)) ME(XPZ,k,j)+=lamti*weights[c]*VE(zi,j)*VE(xi,k);
               if ((j<*pm)&(k<*pg)) ME(cumZP,j,k)+=lamti*weightsmg[c]*VE(vtmp,j)*VE(zi,k);
	      }
	   }
           count=count+1; 
         }		 
	 }
           ci=*nx-1; 
           while ((stop[ci]<time)  & (ci>=0) )  ci=ci-1; 
       } // }}}

//Rprintf("%d %d %d %lf %lf %lf \n",s,ci,id[ci],start[ci],stop[ci],time); 
     vec_zeros(rowX); vec_zeros(rowZ); vec_zeros(dB); 
    if (s<1)  // {{{ modifying design for next time points
    while ((stop[ci]<time)  & (ci>=0) ) {
            for(j=0;j<*px;j++) VE(xi,j)=designX[j*(*nx)+ci]; 
	    if (*coxaalen==1) for(j=0;j<*pg;j++) VE(zi,j)=designG[j*(*nx)+ci]; 
            for(j=0;j<*pm;j++) VE(vtmp,j)=modelmatrix[j*(*nx)+ci]; 
	    if (*coxaalen==1) {
		VE(Gbeta,id[ci])=vec_prod(zi,gamma); 
		RR=exp(VE(Gbeta,id[ci])+offsets[ci]);
	        lamti=RR*vec_prod(xi,dA);
	    }
	    if (entry[ci]==1)  {
	         replace_row(X,id[ci],xi); 
		 replace_row(cummat,id[ci],vtmp); 
//		 scl_vec_mult(RR,xi,tmpv1);replace_row(WX,id[ci],tmpv1);
		 VE(weight,id[ci])=weights[ci]; 
		 VE(offset,id[ci])=offsets[ci];  
	    } 
	    else { 
		    replace_row(X,id[ci],rowX);
		    replace_row(cummat,id[ci],dB);
		    VE(weight,id[ci])=0; 
		    VE(offset,id[ci])=0;  
            }	   
	  for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=entry[ci]*weights[ci]*VE(xi,k)*VE(xi,j)*RR; 
//              if ((j<*pg) & (k<*px)) ME(ZX,j,k)+=entry[ci]*VE(zi,j)*VE(xi,k)*RR; 
              if ((j<*px) & (k<*pm)) ME(cumX,k,j)+=entry[ci]*weights[ci]*VE(xi,j)*VE(vtmp,k)*RR; 
	      if (*coxaalen==1) {
               if ((j<*pg)&(k<*px)) ME(XPZ,k,j)+=lamti*entry[ci]*weights[ci]*VE(zi,j)*VE(xi,k);
               if ((j<*pm)&(k<*pg)) ME(cumZP,j,k)+=lamti*entry[ci]*weights[ci]*VE(vtmp,j)*VE(zi,k);
	      }
	  }
	  ci=ci-1; 
	  pers=id[ci]; 
    // }}}
//   ipers[s]=pers;
   } // }}}

   for (i=0;i<*antpers;i++) VE(dMGt[s],i)=dmgresid[i*(*nmgt)+s];  

   cummgt[s]=time; 
   robvarcum[s]=time; 

   if (*stratum==0) invertS(A,AI,*silent); 
   if (ME(AI,0,0)==0 && *stratum==0 && *silent==0) {
       Rprintf("additive design X'X not invertible at time (number, value): %d %lf \n",s,time); print_mat(A);
   }
   if (*stratum==1)  {
       for (k=0;k<*px;k++) 
       if (fabs(ME(A,k,k))<0.000001)  ME(AI,k,k)=0; else ME(AI,k,k)=1/ME(A,k,k);
   }

   MxA(cumX,AI,cumXAI);

    /* observed increment */ 
    vM(cummat,dMGt[s],respm); 
    for (k=1;k<=*pm;k++) cummgt[k*(*nmgt)+s]=cummgt[k*(*nmgt)+s-1]+VE(respm,k-1);

    /* extra terms for cox aalen iid representation */ 
   if (*coxaalen==1) {
      MxA(cumXAI,XPZ,tmp2);
      mat_subtr(cumZP,tmp2,dS);  
      mat_add(dS,S,S); 
      scl_mat_mult(1,S,St[s]); 
//      printf("=========================== %d \n",s); 
//       print_mat(cumZP); print_mat(tmp2); 
//       print_mat(dS); 
//       print_mat(S); 
   } /*coxaalen=1 */ 

   if (*ratesim==1) 
   for (i=0;i<*antpers;i++) 
   {
    j=cluster[i];
    extract_row(cummat,i,respm); 
    extract_row(X,i,xi); 
    scl_vec_mult(VE(weight,i),xi,xi); 
    Mv(cumXAI,xi,vtmp);
    vec_subtr(respm,vtmp,respm); 
    scl_vec_mult(VE(dMGt[s],i),respm,vtmp);
    vec_add(vtmp,cumA[j],cumA[j]); // cumulative  M_i(t)  = + K_z(t) d\hat M_i(t) 
    replace_row(modelMGT[j],s,cumA[j]);  
  }

   if (*ratesim==0) {
      m=cluster[pers]; 
      extract_row(cummat,pers,respm); 
      extract_row(X,pers,xi); 
      scl_vec_mult(VE(weight,pers),xi,xi); 
      Mv(cumXAI,xi,vtmp);
      vec_subtr(respm,vtmp,respm); 
      vec_add(vtmp,cumA[m],cumA[m]); // cumulative  M_k
      if (coxaalen==0)  for (s1=s;s1<*nmgt;s1++) replace_row(modelMGT[m],s1,cumA[m]); 
//      printf(" %d %d \n",s,pers); 
    }

   vec_zeros(VdB);  vec_zeros(respm); 
//   if (*coxaalen==1 || weighted==1)   
//   if (*coxaalen==1 || weighted==1)   
   for (j=0;j<*antclust;j++)  {
//       if (s==1) { printf(" %d \n",j); print_vec(gammaiid[j]); }
       if (*coxaalen==1)  { 
          Mv(dS,gammaiid[j],respm); 
          vec_subtr(cumA[j],respm,cumA[j]);  // - dS gammaiid(j), total St gammaiid(j) 
        }
        replace_row(modelMGT[j],s,cumA[j]); 
//	if (weighted==1) {
           vec_star(cumA[j],cumA[j],vtmp); 
	   vec_add(vtmp,VdB,VdB); 
//	}
   }
//   if (weighted==1) 
   for (k=1;k<*pm+1;k++) robvarcum[k*(*nmgt)+s]=VE(VdB,k-1); 
 
  /* comp observed sup statistics */ 
//  if (weighted==0)
  Ut[s]=time;  
  for (i=1;i<=*pm;i++) {
    if (weighted==1) vardiv=sqrt(robvarcum[i*(*nmgt)+s]); else vardiv=1; 
    xij=cummgt[i*(*nmgt)+s]/vardiv;
    if (fabs(xij)>testOBS[i-1]) testOBS[i-1]=fabs(xij);
    Ut[i*(*nmgt)+s]=xij; 
    c=(*pm)+i-1; testOBS[c]=testOBS[c]+xij*xij*dtime; 
  } 

 } // }}}  s=1..Loop 

//   vec_zeros(VdB);  vec_zeros(respm); 
//   if (*coxaalen==1 || weighted==1) 
//   for (s=1;s<*nmgt;s++) // {{{ 
//   {
//   vec_zeros(VdB);  vec_zeros(respm); 
//   for (j=0;j<*antclust;j++)  {
//       if (*coxaalen==1)  { 
//          Mv(St[s],gammaiid[j],respm); 
//	  extract_row(modelMGT[j],s,cumdB); 
//          vec_subtr(cumdB,respm,vtmp);
//          replace_row(modelMGT[j],s,vtmp); 
//        } else  extract_row(modelMGT[j],s,vtmp); 
//	if (weighted==1) {
//           vec_star(vtmp,vtmp,cumdB); 
//	   vec_add(cumdB,VdB,VdB); 
//	}
//   }
//   if (weighted==1) 
//   for (k=1;k<*pm+1;k++) robvarcum[k*(*nmgt)+s]=VE(VdB,k-1); 
//
///* comp observed sup statistics */ 
//  if (weighted==1)
//  Ut[s]=time;  
//  for (i=1;i<=*pm;i++) {
//    if (weighted==1) vardiv=sqrt(robvarcum[i*(*nmgt)+s]); else vardiv=1; 
//    xij=cummgt[i*(*nmgt)+s]/vardiv;
//    if (fabs(xij)>testOBS[i-1]) testOBS[i-1]=fabs(xij);
//    Ut[i*(*nmgt)+s]=xij; 
//    c=(*pm)+i-1; testOBS[c]=testOBS[c]+xij*xij*dtime; 
//  }
//
//  }
 
//    for (k=0;k<*antclust;k++) {
//        printf(" %d \n",k);  
//        print_vec(gammaiid[k]); 
//        print_mat(modelMGT[k]); 
//    }


  R_CheckUserInterrupt();
      /* simulation of processes under the model */ 
      for (k=0;k<*sim;k++)  // {{{
      {
        R_CheckUserInterrupt();
	mat_zeros(Delta); 
	for (i=0;i<*antclust;i++) { 
	  random=norm_rand();
	  scl_mat_mult(random,modelMGT[i],tmpM1); mat_add(tmpM1,Delta,Delta); 
	}

	for (s=1;s<*nmgt;s++) { 
	  dtime=mgtimes[s]-mgtimes[s-1]; 

	  for (i=1;i<=*pm;i++) {
	    if (weighted==1) vardiv=sqrt(robvarcum[i*(*nmgt)+s]); else vardiv=1; 
	    xij=ME(Delta,s,i-1)/vardiv; 
	    if (fabs(xij)>test[(*sim)*(i-1)+k]) test[(*sim)*(i-1)+k]=fabs(xij);

	    if (k<50) {l=k*(*pm)+i-1; simUt[l*(*nmgt)+s]=xij;}
      
	    c=*pm+i-1; 
	    test[(*sim)*c+k]=test[(*sim)*c+k]+xij*xij*dtime; 
	    c=2*(*pm)+i-1; 
	    xij=xij/sqrt(robvarcum[i*(*nmgt)+s]); 
	    if (fabs(xij)>test[(*sim)*c+k]) test[(*sim)*c+k]=fabs(xij);
	  }
	}

      } // }}}

  } //  modelmatrix loop }}}
  // }}}
  
    mat_zeros(X); mat_zeros(cummat); vec_zeros(risk);
    mat_zeros(Z);
    mat_zeros(XPZ); 

    R_CheckUserInterrupt();

//  /* LWY cumulative residuals versus covariates */ 
  if (*cumresid>0)  // {{{
  {
 int detail=0; 

  for (l=0;l<ptot;l++) if (ant[l]>2) {
    R_CheckUserInterrupt();
// Rprintf(" %d %d %d %d ================ \n",ptot,ant[l],inXorZ[l],inXZ[l]); 
// for (j=0;j<ant[l];j++)  Rprintf(" %lf \n",xval[(*maxval)*l+j]);

//    mat_zeros(X);mat_zeros(Z); mat_zeros(A); 

    R_CheckUserInterrupt();
 // {{{ allokering for each covariate
 malloc_mat(ant[l],*pg,dS1); malloc_mat(ant[l],*pg,S1);
 malloc_mat(ant[l],*px,cumX1);
 malloc_mat(ant[l],*px,cumXAI1); 
 malloc_mat(ant[l],*pg,cumZP1);
 malloc_mat(ant[l],*pg,tmp21);
 malloc_mat(*antpers,ant[l],cummat1); 
 malloc_mat(*antclust,ant[l],modMGz); 
 malloc_mat(*antclust,ant[l],modMGzosdt); 
 malloc_vecs(ant[l],&vtmp2,&vtmp1,&cumdB1,&VdB1,&respm1,NULL);
 malloc_vec(ant[l],Deltaz); 
 malloc_vec(ant[l],Deltazsd);
 malloc_vec(ant[l],tmpM1z);  
 for (k=0;k<*antclust;k++) malloc_vec(ant[l],dB1[k]); 
 for (k=0;k<*antclust;k++) malloc_vec(ant[l],dBgam[k]); 
 // }}}

    pm[0]=ant[l]; pmax=ant[l]; pmax=max(pmax,*px); 
    if (*coxaalen==1)  pmax=max(pmax,*pg);

    for (s=1;s<*nmgt;s++) // {{{ 
    {
	time=mgtimes[s]; dtime=mgtimes[s]-mgtimes[s-1]; 
	cummgt[s]=time; robvarcum[s]=time; 

	if (*coxaalen==1) 
        for (j=0;j<*px;j++) VE(dA,j)=dcum[j*(*nmgt-1)+s-1]; 

    // reads everything at each "s" (so quite slow)
    
    // {{{ reading design and computing matrix products
	  if (s>=1) { // {{{

          mat_zeros(X);mat_zeros(Z); mat_zeros(A); 
          mat_zeros(cummat1); mat_zeros(A); mat_zeros(cumX1); mat_zeros(XPZ); mat_zeros(cumZP1); 
          vec_zeros(weightmg); vec_zeros(weight); vec_zeros(offset); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	  {
	   if ((start[c]<time) && (stop[c]>=time)) 
	   {
                for(j=0;j<pmax;j++) {
                   if (*coxaalen==1) if (j<*pg) { 
		      ME(Z,id[c],j)=designG[j*(*nx)+c]; 
		      VE(zi,j)=designG[j*(*nx)+c]; 
		   }
	           if (j<*px) {
		       ME(X,id[c],j)=designX[j*(*nx)+c]; 
		       VE(xi,j)=designX[j*(*nx)+c]; 
		   }
	           if (j<*pm) {
                   if (inXorZ[l]==0) ME(cummat1,id[c],j)=((designG[inXZ[l]*(*nx)+c])<=xval[(*maxval)*l+j]);
                                else ME(cummat1,id[c],j)=((designX[inXZ[l]*(*nx)+c])<=xval[(*maxval)*l+j]);
		    VE(vtmp1,j)=ME(cummat1,id[c],j); 
		   }
		}
		if (*coxaalen==1) { VE(Gbeta,id[c])=vec_prod(zi,gamma); 
		                    RR=exp(VE(Gbeta,id[c])+offsets[c]);
	                            lamti=RR*vec_prod(xi,dA);
		} 
		if (time==stop[c] && status[c]==1) {pers=id[c];} 
		VE(offset,id[c])=offsets[c];  
		VE(weight,id[c])=weights[c]; 
		VE(weightmg,id[c])=weightsmg[c]; 
//		printf("%d %d %lf \n",s,id[c],lamti); 

	    for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=weights[c]*VE(xi,k)*VE(xi,j)*RR; 
              if ((j<*px) & (k<*pm)) ME(cumX1,k,j)+=weightsmg[c]*VE(xi,j)*VE(vtmp1,k)*RR;
	      if (*coxaalen==1) {
              if ((j<*pg)&(k<*px)) ME(XPZ,k,j)+=lamti*weights[c]*VE(zi,j)*VE(xi,k);
              if ((j<*pm)&(k<*pg)) ME(cumZP1,j,k)+=lamti*weightsmg[c]*VE(vtmp1,j)*VE(zi,k);
	      }
	   }
           count=count+1; 
         }		 
        }
           ci=*nx-1; 
           while ((stop[ci]<time)  & (ci>=0) )  ci=ci-1; 
     } // }}}

//Rprintf("%d %d %d %lf %lf %lf \n",s,ci,id[ci],start[ci],stop[ci],time); 
     vec_zeros(rowX); vec_zeros(rowZ); vec_zeros(vtmp2); 
    if (s<1)  // {{{ modifying design for next time points
    while ((stop[ci]<time)  & (ci>=0) ) {
            for(j=0;j<pmax;j++) {
	       if (j<*px) VE(xi,j)=designX[j*(*nx)+ci]; 
               if (*coxaalen==1) if (j<*pg) VE(zi,j)=designG[j*(*nx)+ci]; 
	    }
            for(j=0;j<*pm;j++) 
            if (inXorZ[l]==0) VE(vtmp1,j)=((designG[inXZ[l]*(*nx)+ci])<=xval[(*maxval)*l+j]); else VE(vtmp1,j)=((designX[inXZ[l]*(*nx)+ci])<=xval[(*maxval)*l+j]);
	    if (*coxaalen==1) {
	       VE(Gbeta,id[ci])=vec_prod(zi,gamma); 
	       RR=exp(VE(Gbeta,id[ci])+offsets[ci]);
	       lamti=RR*vec_prod(xi,dA);
	    }
	    if (entry[ci]==1)  {
	         replace_row(X,id[ci],xi); 
	         replace_row(Z,id[ci],zi); 
		 replace_row(cummat1,id[ci],vtmp1);
		 VE(weight,id[ci])=weights[ci]; 
		 VE(offset,id[ci])=offsets[ci];  
	    } 
	    else { 
		 replace_row(X,id[ci],rowX);
	         replace_row(Z,id[ci],rowZ); 
		 replace_row(cummat1,id[ci],vtmp2);
		 VE(weight,id[ci])=0; 
		 VE(offset,id[ci])=0;  
	    }
	  for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px)&(k<*px)) ME(A,j,k)+=entry[ci]*weights[ci]*VE(xi,k)*VE(xi,j)*RR; 
              if ((j<*px)&(k<*pm)) ME(cumX1,k,j)+=entry[ci]*weights[ci]*VE(xi,j)*VE(vtmp1,k)*RR; 
	      if (*coxaalen==1) {
               if ((j<*pg)&(k<*px)) ME(XPZ,k,j)+=lamti*entry[ci]*weights[ci]*VE(zi,j)*VE(xi,k);
               if ((j<*pm)&(k<*pg)) ME(cumZP1,j,k)+=lamti*entry[ci]*weights[ci]*VE(vtmp1,j)*VE(zi,k);
	      }
	  }
	  ci=ci-1; 
	  pers=id[ci]; 
    // }}}
//   ipers[s]=pers;
   } // }}}

//    random=0; 
//    if (s > -2) {
//	    head_matrix(X);
//	    head_matrix(Z);
//	    print_mat(Z);
//	    print_ma
//	    for (i=0;i<*antpers;i++) random=random+ME(X,i,0); 
//	    Rprintf(" %d %lf \n",s,random);
//    }

    if (*model==0) for (i=0;i<*antpers;i++) VE(dMGt[s],i)=dmgresid[i*(*nmgt)+s];  

     if (*stratum==0) invertS(A,AI,*silent); 
     if (ME(AI,0,0)==0 && *stratum==0 && *silent==0) {
	   Rprintf("additive design X'X not invertible at time (number, value): %d %lf \n",s,time); print_mat(A);
     }
     if (*stratum==1)  {
          for (k=0;k<*px;k++) 
          if (fabs(ME(A,k,k))<0.000001)  ME(AI,k,k)=0; else ME(AI,k,k)=1/ME(A,k,k);
     }

    MxA(cumX1,AI,cumXAI1);
    if (s < 0) {
//    printf(" %d \n",s); 
    head_matrix(A); 
    head_matrix(AI); 
    head_matrix(cumX1); 
    head_matrix(cummat1); 
    }

    /* observed increment */ 
    vM(cummat1,dMGt[s],respm1); 
    for (j=0;j<ant[l];j++) 
    univarproc[(*maxval)*l+j]=univarproc[(*maxval)*l+j]+VE(respm1,j);

//    vM(Z,dMGt[s],zi); 
//    vec_add(rowZ,zi,rowZ); 
//    printf(" D score proc %d \n",s); 
//    print_vec(zi); 

    /* iid representation increment */ 
    /* extra terms for cox aalen iid representation */ 
    if (*coxaalen==1) { // {{{
	   MxA(cumXAI1,XPZ,tmp21);
	   mat_subtr(cumZP1,tmp21,dS1);  
	   mat_add(dS1,S1,S1);
    }  // }}}

    if (*ratesim==1) 
    for (i=0;i<*antpers;i++)  // {{{
    {
      m=cluster[i]; 
      extract_row(cummat1,i,respm1); 
      extract_row(X,i,xi); 
      scl_vec_mult(VE(weight,i),xi,xi); 
      Mv(cumXAI1,xi,vtmp1);
      vec_subtr(respm1,vtmp1,respm1); 
      scl_vec_mult(VE(dMGt[s],i),respm1,vtmp1); 
      vec_add(dB1[m],vtmp1,dB1[m]); 
    } // }}}

    if (*ratesim==0) {
      m=cluster[pers]; 
      extract_row(cummat1,pers,respm1); 
      extract_row(X,pers,xi); 
      scl_vec_mult(VE(weight,pers),xi,xi); 
      Mv(cumXAI1,xi,vtmp1);
      vec_subtr(respm1,vtmp1,respm1); 
      vec_add(dB1[m],respm1,dB1[m]); 
    }

    }  // }}} /* s=1... *Ntimes */ 

    if (detail==1) Rprintf("Left jump times \n"); 

    vec_zeros(VdB1); 
    for (m=0;m<*antclust;m++)  // {{{
    {
       if (*coxaalen==1){
             Mv(S1,gammaiid[m],respm1);
	     vec_subtr(dB1[m],respm1,dB1[m]);
        }
	for (j=0;j<ant[l];j++)  { ME(modMGz,m,j)=VE(dB1[m],j); 
                                  VE(VdB1,j)+=pow(ME(modMGz,m,j),2);
	}
    }

//    for (m=0;m<1;m++) {
//            Rprintf("%d %d \n",m,*antclust); 
//	    print_vec(dB1[m]); 
//    }
//
//    for (s=0;s<2;s++) {
//            Rprintf("%d %d \n",s,s); 
//	    print_vec(dMGt[s]); 
//    }


//for (j=0;j<ant[l];j++)  Rprintf("%d %lf  \n",l,univarproc[(*maxval)*l+j]);

    /* robust variance */

    for (j=0;j<ant[l];j++) robvarcumz[l*(*maxval)+j]=VE(VdB1,j); 
    for (m=0;m<*antclust;m++)  { 
    for (j=0;j<ant[l];j++) 
       ME(modMGzosdt,m,j)=ME(modMGz,m,j)/sqrt(robvarcumz[l*(*maxval)+j]); 
    }

    for (j=0;j<ant[l];j++) 
    {
      xij=univarproc[(*maxval)*l+j]; 
      if (fabs(xij)>unitimetestOBS[l]) unitimetestOBS[l]=fabs(xij); 
    } 

    R_CheckUserInterrupt();
    /* simulation of testprocesses and teststatistics */ // {{{
   //  Rprintf("Simulations start N= %d \n",*sim);
   
   if (detail==1) Rprintf("Simulations start \n"); 

    for (k=0;k<*sim;k++) { // {{{ 
    R_CheckUserInterrupt();
      for (i=0;i<*antclust;i++) VE(rvec,i)=norm_rand(); 
	vM(modMGz,rvec,Deltaz); 
	vM(modMGzosdt,rvec,Deltazsd); 

	for (j=0;j<ant[l];j++) 
	{
	   xij=VE(Deltaz,j); 
	   if (k<50) {c=k*(ptot)+l; simunivarproc[c*(*maxval)+j]=xij;}
	   if (fabs(xij)>unitest[(*sim)*(l)+k]) unitest[(*sim)*(l)+k]=fabs(xij); 
	   xij=VE(Deltazsd,j); 
	   if (fabs(xij)>simcumz[(*sim)*(l)+k]) simcumz[(*sim)*(l)+k]=fabs(xij); 
	} 
    }  // }}} /* k=1..antsim */ 

    // }}}

       R_CheckUserInterrupt();
// {{{ free allokering local allocation LWY style cum res
free_mats(&S1,&dS1,&cumX1,&cumXAI1,&cumZP1,&tmp21,&cummat1,&modMGz,&modMGzosdt,NULL); 
free_vecs(&vtmp2,&vtmp1,&cumdB1,&VdB1,&respm1,&Deltaz,&Deltazsd,&tmpM1z,NULL);
for (k=0;k<*antclust;k++) free_vec(dB1[k]); 
for (k=0;k<*antclust;k++) free_vec(dBgam[k]); 
// }}}
     
   } /* l=0,...,ptot */ 

  } // }}} cumresid=1

  PutRNGstate();  /* to use R random normals */

// {{{ // freeing variables
  for (i=0;i<*nmgt;i++) { 
     free_vec(dMGt[i]); 
//   free_vec(dNt[i]); 
     free_mat(St[i]); 
  }

  for (i=0;i<*antclust;i++) {
    free_mat(modelMGT[i]); free_vec(gammaiid[i]); free_vec(cumA[i]); 
  }
  free_mats(&cumX,&Delta,&tmpM1,&Z,&X,&A,&AI,&cumXAI,&cumZP,&XPZ,&tmp2,&dS,&S,&cummat,NULL); 
//  free_mats(&AIXZ,NULL); 

  free_vecs(&weightmg,&weight,&offset,&vtmp,&cumdB,&dB,&VdB,&respm,&tmpv1,&cum,&dA,&xtilde,&xi,&rowX,
            &rowcum,&difX,&zi,&gamma,&rowZ,&vecX,&risk,&Gbeta,&rvec,NULL); 
free(cluster); 
// }}}
  
} // }}}
