//#include <stdio.h>
#include <math.h>
#include <R.h>
#include "matrix.h"
#include <time.h>
#include <sys/types.h>
                 
void transsurv(double *times,int *Ntimes,double *designX,int *nx,int *px,int *antpers,double *start,double *stop,double *betaS,int *Nit,double *cu,double *vcu,double *Iinv,
double *Vbeta,int *detail,int *sim,int *antsim,int *rani,double *Rvcu,double *RVbeta,double *test,double *testOBS,double *Ut,double *simUt,double *Uit,int *id,int *status,
int *weighted,int *ratesim,double *score,double *dhatMit,double *dhatMitiid,int *retur,double *loglike,int *profile,int *sym,int *baselinevar,int *clusters,int *antclust,double *biid,double *gamiid)
//double *designX,*times,*betaS,*start,*stop,*cu,*Vbeta,*RVbeta,*vcu,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*score,*dhatMit,*dhatMitiid,*loglike,*biid,*gamiid;
//int *nx,*px,*antpers,*Ntimes,*Nit,*detail,*sim,*antsim,*rani,*id,*status,*weighted,*ratesim,*retur,*profile,*sym,*baselinevar,*clusters,*antclust;
{
// {{{  setting up
  matrix *ldesignX,*WX,*ldesignG,*CtVUCt,*A,*AI;
  matrix *dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *tmp1,*tmp2,*tmp3,*dS1,*SI,*dS2,*S2,*S2pl,*dS2pl,*M1,*VU,*ZXAI,*VUI; 
  matrix *d2S0,*RobVbeta,*tmpM1,*Utt,*S1t,*S1start,*tmpM2,*et,*gt,*qt;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *dotwitowit[*antpers],*W3t[*antclust],*W4t[*antclust],*W2t[*antclust],*AIxit[*antpers],*Uti[*antclust],*d2G[*Ntimes],*Delta,*Delta2; 
  vector *Ctt,*lht,*S1,*dS0,*S0t,*S0start,*dA,*VdA,*dN,*MdA,*delta,*zav,*dlamt,*plamt,*dG[*Ntimes],
    *S1star;
  vector *xav,*difxxav,*xi,*xipers,*zi,*U,*Upl,*beta,*xtilde; 
  vector *Gbeta,*zcol,*one,*difzzav,*difZ; 
  vector *offset,*weight,*ZXdA[*Ntimes],*varUthat[*Ntimes],*Uprofile;
  vector *ahatt,*risk,*tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB; 
  vector *W2[*antclust],*W3[*antclust],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  vector *dLamt[*antpers];
  int *pg=calloc(1,sizeof(int)),c,robust=1,pers=0,ci=0,i,j,k,l,s,t,it,count,*ipers=calloc(*Ntimes,sizeof(int));
  int *cluster=calloc(*antpers,sizeof(int));
  double RR,S0star,time,dummy,ll;
  double S0,tau,random,scale,sumscore;
//  double norm_rand();
//  void GetRNGstate(),PutRNGstate();

  pg[0]=1; 

  for (j=0;j<*antpers;j++) { 
    malloc_vec(*Ntimes,dLamt[j]); malloc_mat(*Ntimes,*px,dotwitowit[j]); 
    malloc_mat(*Ntimes,*px,AIxit[j]); 
  }

  for (j=0;j<*antclust;j++) { 
    malloc_mat(*Ntimes,*pg,W3t[j]); malloc_mat(*Ntimes,*pg,W4t[j]); 
    malloc_mat(*Ntimes,*px,W2t[j]); malloc_vec(*px,W2[j]); malloc_vec(*pg,W3[j]); 
    malloc_mat(*Ntimes,*px,Uti[j]); 
  }

 for(j=0;j<*Ntimes;j++) {
    malloc_mat(*px,*pg,dYIt[j]); malloc_vec(*px,dAt[j]); malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]);
    malloc_mat(*pg,*px,ZXAIs[j]); 
    malloc_vec(*pg,ZXdA[j]); malloc_mat(*px,*px,St[j]); 
    malloc_mat(*px,*px,d2G[j]); malloc_vec(*px,dG[j]); malloc_vec(*px,varUthat[j]);
  }

  malloc_mat(*Ntimes,*pg,tmpM1); 
  malloc_mat(*Ntimes,*px,S1t); malloc_mat(*Ntimes,*px,tmpM2); 
  malloc_mat(*Ntimes,*px,S1start); 
  malloc_mat(*Ntimes,*px,et); malloc_mat(*Ntimes,*px,gt); malloc_mat(*Ntimes,*px,qt); 
  malloc_mat(*Ntimes,*px,Utt); 
  malloc_mat(*Ntimes,*pg,Delta); malloc_mat(*Ntimes,*px,Delta2); 

  malloc_mats(*antpers,*px,&WX,&ldesignX,NULL); 
  malloc_mats(*antpers,*pg,&ZP,&ldesignG,NULL); 
  malloc_mats(*px,*px,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*px,*px,&d2S0,&RobVbeta,&tmp1,&tmp2,&dS1,&S2,&dS2,&S2pl,&dS2pl,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*px,*pg,&tmp3,&ZPX,&dYI,&Ct,NULL); 

  malloc_vec(*Ntimes,S0t);  
  malloc_vec(*Ntimes,S0start); 
  malloc_vec(*Ntimes,lht); 
  malloc_vec(1,reszpbeta); 
  malloc_vec(1,res1dim); 
  malloc_vecs(*antpers,&risk,&weight,&plamt,&dlamt,&dN,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&Ctt,&ahatt,&tmpv1,&difX,&rowX,&xi,&xipers,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&S1,&dS0,&S1star,&xtilde,&xav,&difxxav,NULL); 
  malloc_vecs(*px,&U,&Upl,&beta,&delta,&difzzav,&Uprofile,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&difZ,&zav,&VdB,NULL); 

   ll=0; for(j=0;j<*px;j++) VE(beta,j)=betaS[j]; 
  // }}}

  int timing=0; 
  clock_t c0,c1; 
  c0=clock();

  double plamtj,dlamtj; 
//    mat_zeros(ldesignX); 
    for (c=0;c<*nx;c++) for(j=0;j<*px;j++) ME(WX,id[c],j)=designX[j*(*nx)+c];  

  cu[0]=times[0]; 
  for (it=0;it<*Nit;it++)
  {
      vec_zeros(U); vec_zeros(Upl); mat_zeros(S2pl); mat_zeros(S2); ll=0; sumscore=0; 
      mat_zeros(COV); 

   R_CheckUserInterrupt();

   if (timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

//    for (c=0;c<*nx;c++) for(j=0;j<*px;j++) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];  
    Mv(WX,beta,Gbeta); 

      for (s=1;s<*Ntimes;s++)
      {// {{{  
	  time=times[s];  
	  vec_zeros(dS0);  mat_zeros(d2S0);  mat_zeros(dS1);  vec_zeros(S1star);
	  S0star=0; S0=0; 
	  vec_zeros(S1); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) { // {{{ reading data and computing things 
	    if ((start[c]<time) && (stop[c]>=time))  {
	      for(j=0;j<*px;j++) VE(xi,j)=designX[j*(*nx)+c];  
	      j=id[c];
	      if (time==stop[c] && status[c]==1) {pers=id[c]; scl_vec_mult(1,xi,xipers);} 
	      count=count+1; 
	      RR=exp(-VE(Gbeta,j));

	    scale=(RR+cu[1*(*Ntimes)+s-1]); 
	    dummy=1/scale; 
            if (it==((*Nit)-1)) VE(plamt,j)=dummy; 
	    plamtj=dummy; 
	    dlamtj=dummy*dummy; 
//	    VE(dlamt,j)=dlamtj; 
	    S0star=S0star-dlamtj; 
	    S0=S0+plamtj;

//          S0p=S0p+VE(risk,j)/(RR+cu[1*(*Ntimes)+s]);
//	    S0cox=S0cox+exp(VE(Gbeta,j)); 
//
	    scl_vec_mult(-RR,xi,tmpv1); 
	    vec_add(tmpv1,dG[s-1],dA); 
	    scl_vec_mult(-plamtj,dA,dA); 
	    if (it==(*Nit-1)) { 
	      if (*profile==0) replace_row(dotwitowit[j],s,xi); else replace_row(dotwitowit[j],s,dA); 
	    }
	    scl_vec_mult(plamtj,dA,tmpv1); 
	    vec_add(tmpv1,dS0,dS0); 
	    if (s<0 && j<5 ) {Rprintf(" %d %d \n",s,j); print_vec(tmpv1); }

	    // 16-10-2014                                                   -dlamtj
	    if (*profile==0) scl_vec_mult(-dlamtj,xi,dA); else scl_vec_mult(-plamtj,dA,dA); 
	    vec_add(dA,S1star,S1star); 

	    for (i=0;i<*px;i++) for (k=0;k<*px;k++) { 
		    ME(dS1,i,k)=ME(dS1,i,k)+VE(xi,i)*VE(tmpv1,k); 
	            ME(tmp1,i,k)=(VE(xi,i)*VE(xi,k))*RR;
	    }
	    mat_add(tmp1,d2G[s-1],tmp1);
	    scl_mat_mult(-dlamtj,tmp1,tmp1); 

	    for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(d2S0,i,k)=ME(d2S0,i,k)+ME(tmp1,i,k)+2*scale*(VE(tmpv1,i)*VE(tmpv1,k)); 
	    scl_vec_mult(plamtj,xi,xi); vec_add(S1,xi,S1); 
	    } else VE(plamt,id[c])=0;   
	  } // }}}
	  ipers[s]=pers;

   if (s==1 && timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time:  loop Ntimes  %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

	  replace_row(S1t,s,dS0); 
	  VE(S0start,s)=S0star; 
	  replace_row(S1start,s,S1star); 
	  if (it==((*Nit)-1))  { 
	      VE(S0t,s)=S0; 
	      VE(lht,s)=VE(lht,s-1)-S0star/(S0*S0); 
	  }

	  /* Rprintf(" %ld %lf %lf \n",s,VE(lht,s),ME(AI,0,0)); */

	  scl_vec_mult(S0star,dS0,tmpv1); 
	  scl_vec_mult(S0,S1star,dA); 
	  vec_subtr(tmpv1,dA,dA); 
	  scl_vec_mult(1/S0,dA,dA); 
	  if (it==((*Nit)-1))  replace_row(gt,s,dA); 

	  scl_vec_mult(-1/(S0*S0),dS0,tmpv1); 
	  vec_add(dG[s-1],tmpv1,dG[s]); 

	  if (s<0) { Rprintf(" %lf \n",S0); 
	    print_vec(scl_vec_mult(1/S0,dS0,NULL)); 
	    print_vec(tmpv1); 
	    print_vec(dG[s]); 
	  }

	  scl_mat_mult(-1/(S0*S0),d2S0,A); 
	  for (i=0;i<*px;i++) for (j=0;j<*px;j++) ME(A,i,j)=ME(A,i,j)+2*S0*VE(tmpv1,i)*VE(tmpv1,j);
	  mat_add(d2G[s-1],A,d2G[s]); 

	  /* baseline is computed */ 
	  cu[1*(*Ntimes)+s]=cu[1*(*Ntimes)+s-1]+(1/S0); 
	  if (s<0) Rprintf(" %lf \n",cu[1*(*Ntimes)+s]); 

	  /* First derivative of U ========================================  */ 
//	  extract_row(ldesignX,pers,xi); 
	  scl_vec_mult(1,xipers,xi); 

	  scl_vec_mult(1/S0,S1,xav); vec_subtr(xi,xav,difxxav); vec_add(U,difxxav,U); 
	  if (it==((*Nit)-1))  if (*profile==0)  replace_row(et,s,xav); 

	  /* profile version of score */ 
	  dummy=1/(exp(-VE(Gbeta,pers))+cu[1*(*Ntimes)+s-1]); 
	  scl_vec_mult(-exp(-VE(Gbeta,pers)),xi,tmpv1); 
	  vec_add(tmpv1,dG[s-1],tmpv1); scl_vec_mult(-dummy,tmpv1,tmpv1); 

	  scl_vec_mult(1/S0,dS0,dA); 
	  if (it==((*Nit)-1))  if (*profile==1)  replace_row(et,s,dA); 
	  vec_subtr(tmpv1,dA,dA); 
	  vec_add(Upl,dA,Upl);

	  /* Second derivative S =========================================== */ 

	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(dS2pl,i,k)=(VE(xi,i)*VE(xi,k))*exp(-VE(Gbeta,pers)); 
	  mat_add(dS2pl,d2G[s-1],dS2pl);
	  scl_mat_mult(-dummy,dS2pl,dS2pl); 

	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(dS2pl,i,k)=ME(dS2pl,i,k)+(VE(tmpv1,i)*VE(tmpv1,k)); 

	  scl_mat_mult(-1/S0,d2S0,A); 
	  scl_vec_mult(1/S0,dS0,dA); 
	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(A,i,k)=ME(A,i,k)+(VE(dA,i)*VE(dA,k)); 

	  mat_add(A,dS2pl,dS2pl);
	  mat_add(dS2pl,S2pl,S2pl); 

	  if (*profile==1) St[s]=mat_copy(S2pl,St[s]); 

	  /* simple Second derivative S2 ================================== */
	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(dS2,i,k)=VE(dS0,i)*VE(S1,k); 
	  scl_mat_mult(S0,dS1,M1); 

	  /* */
	  if (s<0) { Rprintf("======================== %lf \n",S0); 
	    print_mat(scl_mat_mult(1/(S0*S0),M1,NULL)); 
	    print_mat(scl_mat_mult(1/(S0*S0),dS2,NULL)); 
	  }

	  mat_subtr(M1,dS2,M1); 

	  if (*sym==1) {
	    scl_mat_mult(-1/(S0*S0),M1,M1); mat_transp(M1,dS2); mat_add(M1,dS2,dS2); 
	    scl_mat_mult(0.5,dS2,dS2); 
	  } else scl_mat_mult(-1/(S0*S0),M1,dS2); 

	  if (s<0) print_mat(dS2); 

	  mat_add(dS2,S2,S2); 
	  if (*profile==0) St[s]=mat_copy(S2,St[s]); 

	  /* ============================================ */
	  /* log-likelihood contributions                 */ 
	  ll=ll+log(dummy)-log(S0);

	  /* scl_mat_mult(1/S0,dS1,dS1);  */

	  if (it==((*Nit)-1)) { // {{{
	    Ut[s]=time; 
	    if (*profile==1) {for (i=1;i<*px+1;i++) Ut[i*(*Ntimes)+s]=VE(Upl,i-1);}  else 
	    {for (i=1;i<*px+1;i++) Ut[i*(*Ntimes)+s]=VE(U,i-1); } 
	    for (i=1;i<*px+1;i++) ME(Utt,s,i-1)=Ut[i*(*Ntimes)+s]; 
	    for (j=0;j<*antpers;j++) VE(dLamt[j],s)=VE(plamt,j)/S0;

	    /* // {{{ 
	      for (i=0;i<*px;i++) for (j=0;j<*pg;j++) ME(dM1M2,j,i)=VE(dA,i)*VE(difzzav,j);
	      for (i=0;i<*pg;i++) 
	      for (j=0;j<*pg;j++) ME(VU,i,j)=ME(VU,i,j)+VE(difzzav,i)*VE(difzzav,j); 

	      MxA(AI,ZPX,dYIt[s]); mat_subtr(Ct,dYIt[s],Ct); C[s]=mat_copy(Ct,C[s]); 

	      vec_star(dA,dA,VdA); mat_add(dM1M2,M1M2t,M1M2t); 
	      M1M2[s]=mat_copy(M1M2t,M1M2[s]); 

	      for (k=1;k<=*px;k++) { cu[k*(*Ntimes)+s]=VE(dA,k-1);
	      vcu[k*(*Ntimes)+s]=VE(VdA,k-1)+vcu[k*(*Ntimes)+s-1]; }
	    */ // }}} 

	  } // }}} 

	} // }}} /* s= .... Ntimes */ 

   if (timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time:  loop Ntimes  %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

      if (*profile==1)  scl_mat_mult(-1,S2pl,A); else scl_mat_mult(-1,S2,A);
      invertS(A,SI,1);
      if (*profile==1) Mv(SI,Upl,delta);  else Mv(SI,U,delta); 
      vec_add(beta,delta,beta); 

      if (*detail==1) {  // {{{ 
	Rprintf("====================Iteration %d ==================== \n",it); 
	Rprintf("log-likelihood "); Rprintf(" %lf \n",ll); 
	Rprintf("Estimate beta "); print_vec(beta); 
	if (*profile==1) {
	  Rprintf("modified partial likelihood Score D l"); print_vec(Upl); }
	if (*profile==0) {Rprintf("simple Score D l"); print_vec(U);  }
	Rprintf("Information -D^2 l\n"); print_mat(SI); 
	if (*profile==1) {Rprintf("simple D2 l");  print_mat(S2pl); }
	if (*profile==0) {Rprintf("simple D2 l");  print_mat(S2); }
      }; // }}} 

      for (k=0;k<*px;k++) sumscore= sumscore+
	(*profile==1)*fabs(VE(Upl,k))+(*profile==0)*fabs(VE(U,k));  

      if ((fabs(sumscore)<0.000001) & (it<*Nit-2)) it=*Nit-2; 
   }  /* it */
  loglike[0]=ll; 

   if (timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time:  out of loop  %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

   R_CheckUserInterrupt();

  for (s=1;s<*Ntimes;s++) { // {{{ /* computation of q(t) ===================================== */ 
    vec_zeros(xi); 
    for (t=s;t<*Ntimes;t++) {
      extract_row(gt,t,dA); 
      scl_vec_mult(exp(VE(lht,t))/VE(S0t,t),dA,dA); 
//      if (s<0) {Rprintf("exp %d %d %lf \n",s,t,exp(-VE(lht,t)+VE(lht,s))); print_vec(dA); }
      vec_add(dA,xi,xi); 
    }
    scl_vec_mult(exp(-VE(lht,s))/VE(S0t,s),xi,xi); 
    replace_row(qt,s,xi); 
  } // }}}

   if (timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time: q(t)     %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

   R_CheckUserInterrupt();

  for (c=0;c<*antpers;c++) cluster[id[c]]=clusters[c]; 

  if (robust==1) { // {{{ terms for robust variances ============================ 
    for (s=1;s<*Ntimes;s++) {
      time=times[s]; 
      cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; Ut[s]=times[s]; 

      extract_row(qt,s,tmpv1); 
      extract_row(et,s,xtilde); 

      for (i=0;i<*antpers;i++) {
        ci=cluster[i]; 
	extract_row(dotwitowit[i],s,rowX); 
	vec_subtr(rowX,xtilde,rowX); 
	if (s==0) { print_vec(rowX); print_vec(tmpv1); }
	vec_add(rowX,tmpv1,rowX); 

	if (i==ipers[s]) for (j=0;j<*px;j++) for (k=0;k<*px;k++)
          	ME(VU,j,k)=ME(VU,j,k)+VE(rowX,j)*VE(rowX,k);

	scl_vec_mult(VE(dLamt[i],s),rowX,xi); 

	vec_subtr(W2[ci],xi,W2[ci]); 
	if (i==ipers[s]) vec_add(rowX,W2[ci],W2[ci]);
	if (*ratesim==1) {scl_vec_mult(VE(dLamt[i],s),tmpv2,rowZ); vec_subtr(W2[ci],rowZ,W2[ci]);}
	replace_row(W2t[ci],s,W2[ci]); 

	VE(rowZ,0)=exp(-VE(lht,s))/VE(S0t,s); 

	scl_vec_mult(VE(dLamt[i],s),rowZ,zi); 
	vec_subtr(W3[ci],zi,W3[ci]); 
	if (i==ipers[s]) vec_add(rowZ,W3[ci],W3[ci]);

	if (*ratesim==1) {scl_vec_mult(VE(dLamt[i],s),rowX,rowX); 
		          vec_subtr(W3[ci],rowX,W3[ci]);
	}
	replace_row(W3t[ci],s,W3[ci]);  

	if (*retur==1)  dhatMit[i*(*Ntimes)+s]=1*(i==pers)-VE(dLamt[i],s);
      } /* i=1..antpers */ 

      /* Martingale baseret variance */
      /*
	MxA(C[s],VU,tmp3); MAt(tmp3,C[s],CtVUCt);
	MxA(C[s],SI,tmp3); MxA(tmp3,M1M2[s],COV); 

	for (k=1;k<=*px;k++) {
  	   cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+cu[k*(*Ntimes)+s]; 
	   vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s]+ME(CtVUCt,k-1,k-1) +2*ME(COV,k-1,k-1); 
	}
      */
    } /* s=1 ..Ntimes */ 



   R_CheckUserInterrupt();


    /* ROBUST VARIANCES  Estimation   */
    for (s=1;s<*Ntimes;s++) {
      vec_zeros(VdB);

      extract_row(S1t,s,rowX); 
      scl_vec_mult(-1.0/(VE(S0t,s)*VE(S0t,s)),rowX,xi); 
      vec_add(xi,Ctt,Ctt); replace_col(Ct,0,Ctt); 

      if (s<0) print_vec(Ctt); 

      for (i=0;i<*antclust;i++) {

	Mv(SI,W2[i],tmpv1); vM(Ct,tmpv1,rowZ);
	extract_row(W3t[i],s,zi); 
	VE(zi,0)= exp(VE(lht,s))*VE(zi,0);

	vec_add(rowZ,zi,zi); 
	if (i==-5) print_vec(zi); 
	replace_row(W4t[i],s,zi);
	biid[i*(*Ntimes)+s]=VE(zi,0);
	vec_star(zi,zi,rowZ); vec_add(rowZ,VdB,VdB);

	extract_row(W2t[i],s,xi); Mv(St[s],tmpv1,rowX); 
	vec_add(xi,rowX,tmpv1); replace_row(Uti[i],s,tmpv1);
  
	vec_star(tmpv1,tmpv1,xi); vec_add(xi,varUthat[s],varUthat[s]);

	if (s==1) { 
	  for (j=0;j<*px;j++) 
	  {
	  gamiid[j*(*antclust)+i]=VE(W2[i],j); 
	  for (k=0;k<*px;k++)
	  ME(RobVbeta,j,k)=ME(RobVbeta,j,k)+VE(W2[i],j)*VE(W2[i],k); 
	  }
	}

	if (*retur==1 && j==0) { // {{{ 

	  for (j=0;j<*antpers;j++)
	    { extract_row(WX,j,xi); // extract_row(ldesignG,j,zi);
	      dummy=exp(VE(Gbeta,j)); // *VE(weight,j)*VE(offset,j); 
	      scl_vec_mult(dummy,xi,xtilde); 
	      replace_row(ldesignX,j,xtilde); 
	    }

	  Mv(ldesignX,dAt[s],dlamt);  
	  for (j=0;j<*antpers;j++)
	  {extract_row(ldesignG,j,zi); scl_vec_mult(VE(dlamt,j),zi,zi); replace_row(ZP,j,zi);} 

	  Mv(ZP,W2[ci],reszpbeta);
	  Mv(dYIt[s],W2[ci],xi); Mv(ldesignX,xi,res1dim); 

	  dhatMitiid[i*(*Ntimes)+s]=dhatMit[i*(*Ntimes)+s]-(VE(reszpbeta,0)- VE(res1dim,0)); 
	} // }}} /* retur ==1 */ 

      } /* i =1 ..Antpers */
      for (k=1;k<*pg+1;k++) Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      for (k=1;k<*pg+1;k++) vcu[k*(*Ntimes)+s]=VE(VdB,k-1);
    }  /*  s=1 ..Ntimes */ 

    MxA(RobVbeta,SI,tmp1); MxA(SI,tmp1,RobVbeta);
  } // }}} /* Robust =1 , default */ 

//   for (i=0;i<*antpers;i++) print_vec(W2[i]); 
//   for (i=0;i<1;i++) print_vec(dLamt[i]); 
//   print_vec(dLamt[0]); 

   if (timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time: variance     %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

  MxA(VU,SI,tmp1); MxA(SI,tmp1,VU);

  for(j=0;j<*px;j++) { 
    betaS[j]= VE(beta,j); 
    if (*profile==1) score[j]=VE(Upl,j); else score[j]=VE(U,j); 
    for (k=0;k<*px;k++){ Iinv[k*(*px)+j]=ME(SI,j,k);
      Vbeta[k*(*px)+j]=-ME(VU,j,k); 
      RVbeta[k*(*px)+j]=-ME(RobVbeta,j,k); 
    } 
  } 

   R_CheckUserInterrupt();

  if (*sim==1) { // {{{ // Rprintf("Simulations start N= %d \n",*antsim);
    GetRNGstate();  /* to use R random normals */

    tau=times[*Ntimes-1]-times[0];
    for (i=1;i<=*pg;i++) VE(rowZ,i-1)=cu[i*(*Ntimes)+(*Ntimes-1)];

    /* Beregning af OBS teststorrelser */
    for (s=1;s<*Ntimes;s++) { time=times[s]-times[0]; 

      for (i=1;i<=*pg;i++) {
	VE(zi,i-1)=fabs(cu[i*(*Ntimes)+s])/sqrt(Rvcu[i*(*Ntimes)+s]);
	if (VE(zi,i-1)>testOBS[i-1]) testOBS[i-1]=VE(zi,i-1); }

      scl_vec_mult(time/tau,rowZ,difZ);
      for (i=1;i<=*pg;i++) VE(zi,i-1)=cu[i*(*Ntimes)+s];
      vec_subtr(zi,difZ,difZ);
      for (i=0;i<*pg;i++) {
	VE(difZ,i)=fabs(VE(difZ,i)); l=(*pg+i);
	if (VE(difZ,i)>testOBS[l]) testOBS[l]=VE(difZ,i);}

      if (*weighted>=1) {  /* sup beregnes i R */ 
	if ((s>*weighted) && (s<*Ntimes-*weighted)) {extract_row(Utt,s,rowX);
	  for (i=0;i<*px;i++) VE(rowX,i)=VE(rowX,i)/sqrt(VE(varUthat[s],i));
	  replace_row(Utt,s,rowX); /* scaled score process */ 
	}
	else {vec_zeros(rowX); replace_row(Utt,s,rowX);}
      }
      for (k=1;k<=*px;k++) Ut[k*(*Ntimes)+s]=ME(Utt,s,k-1);
    } /*s=1..Ntimes Beregning af obs teststorrelser */


    for (k=1;k<=*antsim;k++) {
      mat_zeros(Delta); mat_zeros(Delta2);
      for (i=0;i<*antclust;i++) {
	/* random=gasdev(&idum); */ 
	random=norm_rand();
	scl_mat_mult(random,W4t[i],tmpM1); mat_add(tmpM1,Delta,Delta);
	scl_mat_mult(random,Uti[i],tmpM2); mat_add(tmpM2,Delta2,Delta2);
      }

      extract_row(Delta,*Ntimes-1,zav); 

      for (s=1;s<*Ntimes;s++) { time=times[s]-times[0];
	scl_vec_mult(time/tau,zav,zi); extract_row(Delta,s,rowZ);
	vec_subtr(rowZ,zi,difZ);

	for (i=0;i<*pg;i++) {
	  VE(difZ,i)=fabs(VE(difZ,i));
	  l=(*pg+i);
	  if (VE(difZ,i)>test[l*(*antsim)+k-1]) test[l*(*antsim)+k-1]=VE(difZ,i);
	  VE(zi,i)=fabs(ME(Delta,s,i))/sqrt(Rvcu[(i+1)*(*Ntimes)+s]);
	  if (VE(zi,i)>test[i*(*antsim)+k-1]) test[i*(*antsim)+k-1]=VE(zi,i); }

	if (*weighted>=1) {
	  extract_row(Delta2,s,xi); 
	  if ((s>*weighted) && (s<*Ntimes-*weighted))  {
	    for (i=0;i<*px;i++) {VE(xi,i)=fabs(ME(Delta2,s,i))/sqrt(VE(varUthat[s],i));
	      if (VE(xi,i)>simUt[i*(*antsim)+k-1]) simUt[i*(*antsim)+k-1]=VE(xi,i); }

	    if (k<50) { 
	      for (i=0;i<*px;i++) { l=(k-1)*(*px)+i;
		Uit[l*(*Ntimes)+s]=ME(Delta2,s,i)/sqrt(VE(varUthat[s],i));}}
	  } 
	} /* weigted score */
	else {
	  extract_row(Delta2,s,xi); 
	  for (i=0;i<*px;i++) {
	    if (fabs(VE(xi,i))>simUt[i*(*antsim)+k-1]) simUt[i*(*antsim)+k-1]=fabs(VE(xi,i)); }

	  if (k<50) { 
	    for (i=0;i<*px;i++) { l=(k-1)*(*px)+i;
	      Uit[l*(*Ntimes)+s]=ME(Delta2,s,i);}}
	} /* else wscore=0 */ 

      }  /* s=1..Ntims */
    }  /* k=1..antsim */
    PutRNGstate();  /* to use R random normals */
  } // }}} sim==1 

  // {{{ freeing 
  free_mats(&tmpM1,&S1t,&tmpM2,&S1start,&et,&gt,&qt,&Utt,&Delta,&Delta2,&ldesignX,&ZP,&WX,&ldesignG,&COV,&A,&AI,&M1,&CtVUCt
  ,&d2S0,&RobVbeta,&tmp1,&tmp2,&dS1,&S2,&dS2,&S2pl,&dS2pl,&SI,&VU,&VUI, &ZXAI,&ZX,&dM1M2,&M1M2t
  ,&tmp3,&ZPX,&dYI,&Ct,NULL); 

 free_vecs(&S0t  ,&S0start,&lht ,&reszpbeta,&res1dim, &risk,&weight,&plamt,&dlamt,&dN,&zcol,&Gbeta,&one,&offset
  ,&Ctt,&ahatt,&tmpv1,&difX,&rowX,&xi,&xipers,&dA,&VdA,&MdA,&S1,&dS0,&S1star,&xtilde,&xav,&difxxav
  ,&U,&Upl,&beta,&delta,&difzzav,&Uprofile, &tmpv2,&rowZ,&zi,&difZ,&zav,&VdB,NULL); 

  for (j=0;j<*antpers;j++) {
    free_vec(dLamt[j]); free_mat(dotwitowit[j]); free_mat(AIxit[j]); 
  }
  for (j=0;j<*antclust;j++) {
    free_mat(W3t[j]); free_mat(W4t[j]); 
    free_mat(W2t[j]);free_vec(W2[j]); free_vec(W3[j]); 
    free_mat(Uti[j]); 
  }

  for (j=0;j<*Ntimes;j++) {
    free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]);free_mat(M1M2[j]);free_mat(ZXAIs[j]);
    free_vec(ZXdA[j]); free_mat(St[j]); free_mat(d2G[j]); free_vec(dG[j]);  
    free_vec(varUthat[j]);
  } 
    free(ipers); free(pg); free(cluster); 
    // }}}
    
}
