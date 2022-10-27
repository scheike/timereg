#include <stdio.h>
#include <math.h>
#include <R.h>
#include "matrix.h"
#include <time.h>
#include <sys/types.h>
                 
void posubdist2(double *times,int *Ntimes,double *designX,int *nx,int *px,int *antpers,double *start,double *stop,double *betaS,int *Nit,double *cu,double *vcu,double *Iinv,
double *Vbeta,int *detail,int *sim,int *antsim,int *rani,double *Rvcu,double *RVbeta,double *test,double *testOBS,double *Ut,double *simUt,double *Uit,int *id,int *status,
int *weighted,int *ratesim,double *score,double *dhatMit,double *dhatMitiid,int *retur,double *loglike,int *profile,int *sym,
double *KMtimes,double *KMti,double *etime,int *causeS,int *ipers,int *baselinevar,int *clusters,int *antclust,int *ccode,double *biid,double *gamiid,double *wweights)
//double *designX,*times,*betaS,*start,*stop,*cu,*Vbeta,*RVbeta,*vcu,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*score,*dhatMit,*dhatMitiid,*loglike,
//       *KMtimes,*KMti,*etime,*biid,*gamiid,*wweights;
//int *nx,*px,*antpers,*Ntimes,*Nit,*detail,*sim,*antsim,*rani,*id,*status,*weighted,*ratesim,*retur,*profile,*sym,*causeS,*ipers,*baselinevar,*clusters,*antclust,*ccode; 
{
// {{{  setting up
  matrix *ldesignX,*WX,*ldesignG,*CtVUCt,*A,*AI;
  matrix *dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *tmp1,*tmp2,*tmp3,*dS1,*SI,*dS2,*S2,*S2pl,*dS2pl,*M1,*VU,*ZXAI,*VUI; 
  matrix *d2S0,*RobVbeta,*tmpM1,*Utt,*dS0t,*S1start,*tmpM2,*et,*gt,*qt;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *dotwitowit[*antpers], // *W3tmg[*antclust],
	 *W3t[*antclust],*W4t[*antclust],*W2t[*antclust],*AIxit[*antpers],*Uti[*antclust],*d2G[*Ntimes],*Delta,*Delta2; 
  vector *Ctt,*lht,*S1,*dS0,*incS0t,*S0t,*S0start,*dA,*VdA,*dN,*MdA,*delta,*zav,*dlamt,*plamt,*dG[*Ntimes],
    *S1star;
  vector *xav,*difxxav,*xi,*zi,*U,*Upl,*beta,*xtilde; 
  vector *Gbeta,*zcol,*one,*difzzav,*difZ,*neta2[*antclust]; 
  vector *offset,*weight,*ZXdA[*Ntimes],*varUthat[*Ntimes],*Uprofile;
  vector *ahatt,*risk,*tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB,*VdBmg; 
  vector *W2[*antclust],*W3[*antclust],*reszpbeta,*res1dim,*dAt[*Ntimes],*eta2; 
//  vector *W2[*antclust],*W3[*antclust],*W3mg[*antclust],*reszpbeta,*res1dim,*dAt[*Ntimes],*eta2; 
  vector *dLamt[*antpers];
  int *pg=calloc(1,sizeof(int)),c,robust=1,pers=0,ci,i,j,k,l,s,it;
  double weights,risks,RR,S0star,time,alpha,ll;
  double S0,tau,random,scale,sumscore;
//  double norm_rand();
//  void GetRNGstate(),PutRNGstate();

  pg[0]=1; 

  for (j=0;j<*antpers;j++) { 
    malloc_vec(*Ntimes,dLamt[j]); malloc_mat(*Ntimes,*px,dotwitowit[j]); 
    malloc_mat(*Ntimes,*px,AIxit[j]); 
  }

  for (j=0;j<*antclust;j++) { 
    malloc_mat(*Ntimes,*pg,W3t[j]); 
//    malloc_mat(*Ntimes,*pg,W3tmg[j]); 
    malloc_mat(*Ntimes,*pg,W4t[j]); 
    malloc_mat(*Ntimes,*px,W2t[j]); malloc_mat(*Ntimes,*px,Uti[j]); 
    malloc_vec(*px,W2[j]); 
    malloc_vec(*pg,W3[j]); 
//    malloc_vec(*pg,W3mg[j]); 
    malloc_vec(*Ntimes,neta2[j])
  }

  malloc_mat(*Ntimes,*pg,tmpM1); 
  malloc_mat(*Ntimes,*px,dS0t); malloc_mat(*Ntimes,*px,tmpM2); 
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

  malloc_vec(*Ntimes,S0t);  malloc_vec(*Ntimes,incS0t);  malloc_vec(*Ntimes,eta2);  
  malloc_vec(*Ntimes,S0start); malloc_vec(*Ntimes,lht); malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 
  malloc_vecs(*antpers,&risk,&weight,&plamt,&dlamt,&dN,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&Ctt,&ahatt,&tmpv1,&difX,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&S1,&dS0,&S1star,&xtilde,&xav,&difxxav,NULL); 
  malloc_vecs(*px,&U,&Upl,&beta,&delta,&difzzav,&Uprofile,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&difZ,&zav,&VdB,&VdBmg,NULL); 

  for(j=0;j<*Ntimes;j++) {
    malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]);
    malloc_mat(*pg,*px,ZXAIs[j]); malloc_mat(*px,*pg,dYIt[j]); malloc_vec(*px,dAt[j]);
    malloc_vec(*pg,ZXdA[j]); malloc_mat(*px,*px,St[j]); 
    malloc_mat(*px,*px,d2G[j]); malloc_vec(*px,dG[j]); malloc_vec(*px,varUthat[j]);
  }

  ll=0; for(j=0;j<*px;j++) VE(beta,j)=betaS[j]; 
  // }}}

  int timing=0; 
  clock_t c0,c1; 
  c0=clock();

  double dummy,plamtj,dlamtj,weightp=0; 
  // reading design once and for all
  for (c=0;c<*nx;c++) for(j=0;j<*px;j++) ME(WX,id[c],j)=designX[j*(*nx)+c];  

  cu[0]=times[0]; 
  for (it=0;it<*Nit;it++)
  { // {{{ 
      vec_zeros(U); vec_zeros(Upl); mat_zeros(S2pl); mat_zeros(S2); mat_zeros(COV); 
      ll=0; sumscore=0; 

      R_CheckUserInterrupt();

      Mv(WX,beta,Gbeta); 

      for (s=1;s<*Ntimes;s++)
      {// {{{  
	  time=times[s];  
	  pers=ipers[s];  // person with type 1 jump
          // printf(" pers=%d weight=%lf  cause=%d \n",pers,wweights[pers],status[pers]); 
	  vec_zeros(dS0);  mat_zeros(d2S0);  mat_zeros(dS1);  vec_zeros(S1star); vec_zeros(S1); 
	  S0star=0; S0=0; // S0p=0; S0cox=0; 
	  weightp=1; 

          for (j=0;j<*antpers;j++) { // {{{ 
		  int other=((status[j]!=*causeS) &&  (status[j]!=*ccode))*1; 
		  weights=1; 
		  if (etime[j]<time && other==1) weights=KMtimes[s]/KMti[j]; 
                  if (etime[j]<time) {
		     if (other==1) risks=1; else risks=0; 
		  } else risks=1; 
		  weights=weights*risks; // censoring weights
		  weights=weights*wweights[j];
	          if (j==pers) weightp=weights; 

//                  if (isnan(weights))  {
//                    Rprintf("%lf %lf %d %d %d \n",etime[j],time,*ccode,status[j],*causeS);
//		      Rprintf("%lf %lf \n",risks,weights); 
//                  }

		  extract_row(WX,j,xi); 
		  RR=exp(-VE(Gbeta,j));
		  scale=(RR+cu[1*(*Ntimes)+s-1]); 
		  alpha=1/scale; 
		  plamtj=alpha; 
		  dlamtj=alpha*alpha; 
		  S0star=S0star-dlamtj*weights; 
		  S0=S0+plamtj*weights;
		  if (it==((*Nit)-1)) { 
		     VE(plamt,j)=plamtj; 
	             VE(dLamt[j],s)=weights*plamtj; 
		  }

		  scl_vec_mult(-RR,xi,tmpv1); 
		  vec_add(tmpv1,dG[s-1],dA); 
		  scl_vec_mult(-plamtj,dA,dA); 
		  if (it==(*Nit-1)) { 
		      if (*profile==0) replace_row(dotwitowit[j],s,xi); else replace_row(dotwitowit[j],s,dA); 
		  }
		  scl_vec_mult(plamtj,dA,tmpv1); 
		  vec_add_mult(dS0,tmpv1,weights,dS0); 
		  // dA= (DG- xi RR)/alpha, tmpv1= dA/alpha
                  
		  if (*profile==0) scl_vec_mult(-dlamtj,xi,dA); else scl_vec_mult(-dlamtj,dA,dA); 
		  vec_add_mult(S1star,dA,weights,S1star); 

		  for (i=0;i<*px;i++) for (k=0;k<*px;k++) { 
		       ME(dS1,i,k)=ME(dS1,i,k)+VE(xi,i)*VE(tmpv1,k)*weights;
		       ME(tmp1,i,k)=(VE(xi,i)*VE(xi,k))*RR;
		  }
		  mat_add(tmp1,d2G[s-1],tmp1);
		  scl_mat_mult(-dlamtj,tmp1,tmp1); 

		  for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
		  ME(d2S0,i,k)=ME(d2S0,i,k)+ME(tmp1,i,k)*weights+2*scale*(VE(tmpv1,i)*VE(tmpv1,k))*weights; 

		  vec_add_mult(S1,xi,plamtj*weights,S1); 
	    }   // }}} j=1..antpers 

	  replace_row(dS0t,s,dS0); 
	  VE(S0start,s)=S0star; 
	  replace_row(S1start,s,S1star); 
	  if (it==((*Nit)-1))  { 
//		  printf(" %lf %lf \n",weightp,S0); 
	      VE(incS0t,s)=weightp/S0; 
	      VE(S0t,s)=S0; 
	      VE(lht,s)=VE(lht,s-1)+weightp*S0star/(S0*S0); 
	  }
	  // int S0* / S0 dG_0 =  int ( S0* / S0) (sum_i weight_i dN_i) 

	  scl_vec_mult(S0star,dS0,tmpv1); 
	  scl_vec_mult(S0,S1star,dA); 
	  vec_subtr(tmpv1,dA,dA); 
	  scl_vec_mult(1/S0,dA,dA); 
	  if (it==((*Nit)-1))  replace_row(gt,s,dA); 
	  // g(t) = (SO^* dS0 - S1* S0)/S0 

	  scl_vec_mult(-1/(S0*S0),dS0,tmpv1); 
	  vec_add_mult(dG[s-1],tmpv1,weightp,dG[s]); 
	  // dG(s) = int dS0 / S0^2 (sum_i weight_i dN_i)

	  scl_mat_mult(-1/(S0*S0),d2S0,A); 
	  for (i=0;i<*px;i++) for (j=0;j<*px;j++) ME(A,i,j)=ME(A,i,j)+2*S0*VE(tmpv1,i)*VE(tmpv1,j);
	  scl_mat_mult(weightp,A,A); 
	  mat_add(d2G[s-1],A,d2G[s]); 
	  // d2G(s) = 

	  /* baseline is computed */ 
	  cu[1*(*Ntimes)+s]=cu[1*(*Ntimes)+s-1]+(weightp/S0); 
          if (s<0) Rprintf(" %lf %lf %d %lf \n",cu[1*(*Ntimes)+s],weightp,pers,wweights[pers]); 

	  /* First derivative of U ========================================  */ 
	  // {{{ 
	  extract_row(WX,pers,xi); 
	  scl_vec_mult(1/S0,S1,xav); 
	  vec_subtr(xi,xav,difxxav); 
	  scl_vec_mult(weightp,difxxav,difxxav); 
	  vec_add(U,difxxav,U); 
	  if (it==((*Nit)-1))  if (*profile==0)  replace_row(et,s,xav); 

	  /* profile version of score */ 
	  alpha=1/(exp(-VE(Gbeta,pers))+cu[1*(*Ntimes)+s-1]); 
	  scl_vec_mult(-exp(-VE(Gbeta,pers)),xi,tmpv1); 
	  vec_add(tmpv1,dG[s-1],tmpv1); 
	  scl_vec_mult(-alpha,tmpv1,tmpv1); 
//	  printf("----------- %ld %lf %lf %lf  \n",pers,alpha,S0,weightp); 
//	  print_vec(tmpv1); 

	  scl_vec_mult(1/S0,dS0,dA); 
	  if (it==((*Nit)-1))  if (*profile==1)  replace_row(et,s,dA); 
	  vec_subtr(tmpv1,dA,dA); 
	  scl_vec_mult(weightp,dA,dA); 
//	  print_vec(dA); 
//	  print_vec(dS0); 
	  vec_add(Upl,dA,Upl);
	  // }}} 

	  /* Second derivative S =========================================== */ 
	  // {{{ 
	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(dS2pl,i,k)=(VE(xi,i)*VE(xi,k))*exp(-VE(Gbeta,pers)); 
	  mat_add(dS2pl,d2G[s-1],dS2pl);
	  scl_mat_mult(-alpha,dS2pl,dS2pl); 

	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(dS2pl,i,k)=ME(dS2pl,i,k)+(VE(tmpv1,i)*VE(tmpv1,k)); 

	  scl_mat_mult(-1/S0,d2S0,A); 
	  scl_vec_mult(1/S0,dS0,dA); 
	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(A,i,k)=ME(A,i,k)+(VE(dA,i)*VE(dA,k)); 

	  mat_add(A,dS2pl,dS2pl);
	  scl_mat_mult(weightp,dS2pl,dS2pl); 
	  mat_add(dS2pl,S2pl,S2pl); 

	  if (*profile==1) St[s]=mat_copy(S2pl,St[s]); 

	  /* simple Second derivative S2 ================================== */
	  // {{{ 
	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) ME(dS2,i,k)=VE(dS0,i)*VE(S1,k)*weightp; 
	  scl_mat_mult(S0*weightp,dS1,M1); 

	  if (s<0) { Rprintf("======================== %lf \n",S0); 
	    print_mat(scl_mat_mult(1/(S0*S0),M1,NULL)); 
	    print_mat(scl_mat_mult(1/(S0*S0),dS2,NULL)); 
	  }

          mat_subtr(M1,dS2,M1); 
          scl_mat_mult(-1/(S0*S0),M1,M1); 

	  if (*sym==1) {
	    mat_transp(M1,dS2); mat_add(M1,dS2,dS2); 
	    scl_mat_mult(0.5,dS2,dS2); 
	  } else scl_mat_mult(1,M1,dS2); 
	  mat_add(dS2,S2,S2); 
	  if (*profile==0) St[s]=mat_copy(S2,St[s]); 
	  // }}} 
	  
	  // }}} 

	  /* ============================================ */
	  /* log-likelihood contributions   related to profile score     */ 
	  ll=ll+weightp*(log(alpha)-log(S0));

	  /* scl_mat_mult(1/S0,dS1,dS1);  */

	  if (it==((*Nit)-1)) { // {{{
	    Ut[s]=time; 
	    if (*profile==1) {for (i=1;i<*px+1;i++) Ut[i*(*Ntimes)+s]=VE(Upl,i-1);}  else 
	    {for (i=1;i<*px+1;i++) Ut[i*(*Ntimes)+s]=VE(U,i-1); } 
	    for (i=1;i<*px+1;i++) ME(Utt,s,i-1)=Ut[i*(*Ntimes)+s]; 
	    // dLamt(j,s)   alpha dH = alpha (sum_i w_i(s) dN_i(s)) 

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
//      scl_vec_mult(*step,delta,delta); 

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

      // updating beta 
if (it<(*Nit-1))  vec_add(beta,delta,beta); 

      for (k=0;k<*px;k++) sumscore= sumscore+
	(*profile==1)*fabs(VE(Upl,k))+(*profile==0)*fabs(VE(U,k));  

      if (isnan(sumscore)) oops("missing values in score \n"); 

      if ((fabs(sumscore)<0.0001) & (it<*Nit-2)) it=*Nit-2; 
   }  /* it */ // }}} 

  loglike[0]=ll; 

   if (timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time:   out of loop  %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

   R_CheckUserInterrupt();

//  if (*baselinevar==1) 
    vec_zeros(xi);dummy=0;  
  if (robust==1) 
  for (s=*Ntimes-1;s>0;s--) { // {{{ /* computation of q(t) =============   */ 
      extract_row(gt,s,dA); 
      scl_vec_mult(exp(-VE(lht,s))*VE(incS0t,s),dA,dA); 
      vec_add(dA,xi,xi); 
      replace_row(qt,s,xi); 
      dummy=dummy+exp(-VE(lht,s))*VE(incS0t,s)/VE(S0t,s);
      VE(eta2,s)=dummy; 
  } 
  // q(t)       = int_t^infty ((s_o^star(s)  s_1(s) - s_1^star(s) s_0(s) )/ s_0(s)) k(s) dH_0(s) 
  // q(t)       = int_t^infty g(s) * k(s) dH_0(s) 
  // hat q(t)   = int_t^infty g(s) k(s) (sum_i  w_i(s) (1/s_0(s)) dN_i(s)) 
  // hat eta(t)   = int_t^infty 1 /( s_0(s) k(s)) (sum_i  w_i(s) (1/s_0(s)) dN_i(s)) 
  // k(s) = exp(-lht(s))
  // }}}


   if (timing==1) { // {{{ 
	  c1=clock();
	  Rprintf ("\telapsed CPU time: q(t)     %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
	  c0=clock();
   } // }}} 	

  R_CheckUserInterrupt();

  int *cluster=calloc(*antpers,sizeof(int));
  for (c=0;c<*antpers;c++) cluster[c]=clusters[c]; 

  if (robust==1) { // {{{ terms for robust variances ============================ 
    for (s=1;s<*Ntimes;s++) {

      R_CheckUserInterrupt();

      time=times[s]; 
      cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; Ut[s]=times[s]; 
      pers=ipers[s];  // person with type 1 jump

      extract_row(qt,s,xav); 
      extract_row(et,s,xtilde); 

      for (i=0;i<*antpers;i++) {

	  VE(dLamt[i],s)=VE(dLamt[i],s)*VE(incS0t,s);
	  int other=((status[i]!=*causeS) &&  (status[i]!=*ccode))*1; 
	  weights=1; 
	  if (etime[i]<time && other==1) weights=KMtimes[s]/KMti[i]; 
	  if (etime[i]<time) {
	     if (other==1) risks=1; else risks=0; 
	  } else risks=1; 
	  weights=weights*risks*wweights[i]; // censoring weights

        ci=cluster[i]; 
	extract_row(dotwitowit[i],s,rowX); 
	vec_subtr(rowX,xtilde,rowX); 
	if (s==0) { print_vec(rowX); print_vec(tmpv1); }

	scl_vec_mult(exp(VE(lht,s))/VE(S0t,s),xav,tmpv1); 
	
	vec_subtr(rowX,tmpv1,rowX); 
	scl_vec_mult(weights,rowX,difX); 
	if (i==ipers[s]) vec_add(difX,W2[ci],W2[ci]);

	if (i==ipers[s]) for (j=0;j<*px;j++) for (k=0;k<*px;k++)
          	ME(VU,j,k)=ME(VU,j,k)+VE(difX,j)*VE(difX,k);

	if (*ratesim==1) {
		scl_vec_mult(VE(dLamt[i],s),rowX,tmpv1); 
		vec_subtr(W2[ci],tmpv1,W2[ci]);
	}
	replace_row(W2t[ci],s,W2[ci]); 

	VE(rowZ,0)=exp(VE(lht,s))/VE(S0t,s); 
        if (i==ipers[s])  
	{ 
		VE(W3[ci],0)+=VE(rowZ,0)*weights; 
//		vec_add(rowZ,W3[ci],W3[ci]);
//		vec_add(rowZ,W3mg[ci],W3mg[ci]);
//	        replace_row(W3tmg[ci],s,W3[ci]); 
	}


//	scl_vec_mult(VE(dLamt[i],s),rowZ,zi); 
//		vec_subtr(W3[ci],zi,W3[ci]); 
//        vec_subtr(W3[ci],zi,W3[ci]); 

	if (*ratesim==1) {
//	printf("%ld %ld %ld %lf %lf  \n",ci,s,status[i],weights,etime[i]); 
	scl_vec_mult(VE(dLamt[i],s),rowZ,zi); 
        vec_subtr(W3[ci],zi,W3[ci]); 
	} 
	replace_row(W3t[ci],s,W3[ci]); 

	if (*retur==1)  dhatMit[i*(*Ntimes)+s]=(1*(i==pers)-VE(dLamt[i],s))*weights;
      } /* i=1..antpers */ 

      dummy=0; 
//      for (i=0;i<*antpers;i++) dummy+=VE(dLamt[i],s);  
//      printf("%ld  %lf %lf %lf \n",s,dummy,VE(S0t,s),VE(incS0t,s)); 


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
      R_CheckUserInterrupt();
      vec_zeros(VdB);
//      vec_zeros(VdBmg);

      extract_row(dS0t,s,rowX); 
      scl_vec_mult((-1.0/VE(S0t,s))*VE(incS0t,s),rowX,xi); 
      vec_add(xi,Ctt,Ctt); replace_col(Ct,0,Ctt); 
      // C(t) = D_beta H() = int DS_0(t) / S_0^2(t) (sum_i w_i dN_i) 

      if (s<0) print_vec(Ctt); 
      extract_row(qt,s,xav); 

      for (i=0;i<*antclust;i++) {

	Mv(SI,W2[i],tmpv1); vM(Ct,tmpv1,rowZ);
	extract_row(W3t[i],s,zi); 
	VE(zi,0)= exp(-VE(lht,s))*VE(zi,0);
//	extract_row(W3tmg[i],s,tmpv2); 
//	VE(tmpv2,0)= exp(-VE(lht,s))*VE(tmpv2,0);
//	TESTER uden beta var 
 	vec_add(rowZ,zi,zi); 

	if (i==-5) print_vec(zi); 
	replace_row(W4t[i],s,zi);
	biid[i*(*Ntimes)+s]=VE(zi,0);
	vec_star(zi,zi,rowZ); 
	vec_add(rowZ,VdB,VdB); 
//	vec_star(tmpv2,tmpv2,difZ); vec_add(difZ,VdBmg,VdBmg);

	// U_i(s,hat beta)  for simulations for simulations
        // must subtract W3i(t) * q(t) from W2i(t)	
	// n_i(t) = int_0^t w_i(s) ( Z_i - e(s) + w(s,t)/(k(s) s_0(s))) dM^1_i(s) 
	// w(s,t)= q(s) - q(t) = int_s^t g(t) * k(s) dH_0(s) 
	// W2it= int_0^t w_i(s) ( Z_i - e(s) + q(s)/(k(s) s_0(s))) dM^1_i(s) 
	// W3it= int_0^t w_i(s) 1 /(k(s) s_0(s))) dM^1_i(s) 
	extract_row(W2t[i],s,xi); 
	scl_vec_mult(ME(W3t[i],s,0),xav,xtilde); 
	vec_subtr(xi,xtilde,xi); 
	Mv(St[s],tmpv1,rowX); 
	vec_add(xi,rowX,tmpv1); 
	replace_row(Uti[i],s,tmpv1);

	vec_star(tmpv1,tmpv1,xi); vec_add(xi,varUthat[s],varUthat[s]);

	if (s==1) { 
	  for (j=0;j<*px;j++) {
	  gamiid[j*(*antclust)+i]=VE(W2[i],j); 
	  for (k=0;k<*px;k++) ME(RobVbeta,j,k)=ME(RobVbeta,j,k)+VE(W2[i],j)*VE(W2[i],k); 
	  }
	}

	if (*retur==1 && j==0) { // {{{ 

	  for (j=0;j<*antpers;j++)
	  { extract_row(WX,j,xi); // extract_row(ldesignG,j,zi);
	      alpha=exp(VE(Gbeta,j)); // *VE(weight,j)*VE(offset,j); 
	      scl_vec_mult(alpha,xi,xtilde); 
	      replace_row(ldesignX,j,xtilde); 
	  }

	  Mv(ldesignX,dAt[s],dlamt);  
	  for (j=0;j<*antpers;j++)
	  {extract_row(ldesignG,j,zi); scl_vec_mult(VE(dlamt,j),zi,zi); replace_row(ZP,j,zi);} 

	  Mv(ZP,W2[i],reszpbeta);
	  Mv(dYIt[s],W2[i],xi); Mv(ldesignX,xi,res1dim); 

	  dhatMitiid[i*(*Ntimes)+s]=dhatMit[i*(*Ntimes)+s]-(VE(reszpbeta,0)- VE(res1dim,0)); 
	} // }}} /* retur ==1 */ 

      } /* i =1 ..Antclust*/
      for (k=1;k<*pg+1;k++) Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      for (k=1;k<*pg+1;k++) vcu[k*(*Ntimes)+s]=VE(VdB,k-1);
     /*  s=1 ..Ntimes */ 
    }
//   }
//   else {
//       for (i=0;i<*antclust;i++) 
//       for (j=0;j<*px;j++) 
//        {
//	  gamiid[j*(*antclust)+i]=VE(W2[i],j); 
//	  for (k=0;k<*px;k++)
//	  ME(RobVbeta,j,k)=ME(RobVbeta,j,k)+VE(W2[i],j)*VE(W2[i],k); 
//	  }
//   }

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
  free_mats(&tmpM1,&dS0t,&tmpM2,&S1start,&et,&gt,&qt,&Utt,&Delta,&Delta2,&ldesignX,&ZP,&WX,&ldesignG,&COV,&A,&AI,&M1,&CtVUCt
  ,&d2S0,&RobVbeta,&tmp1,&tmp2,&dS1,&S2,&dS2,&S2pl,&dS2pl,&SI,&VU,&VUI, &ZXAI,&ZX,&dM1M2,&M1M2t
  ,&tmp3,&ZPX,&dYI,&Ct,NULL); 

 free_vecs(&eta2,&incS0t,&S0t,&S0start,&lht ,&reszpbeta,&res1dim, &risk,&weight,&plamt,&dlamt,&dN,&zcol,&Gbeta,&one,&offset
  ,&Ctt,&ahatt,&tmpv1,&difX,&rowX,&xi,&dA,&VdA,&MdA,&S1,&dS0,&S1star,&xtilde,&xav,&difxxav
  ,&U,&Upl,&beta,&delta,&difzzav,&Uprofile, &tmpv2,&rowZ,&zi,&difZ,&zav,&VdB,&VdBmg,NULL); 

  for (j=0;j<*antpers;j++) {
    free_vec(dLamt[j]); free_mat(dotwitowit[j]); free_mat(AIxit[j]); 
  }
  for (j=0;j<*antclust;j++) {
    free_mat(W3t[j]); free_mat(W4t[j]); 
    free_mat(W2t[j]); free_mat(Uti[j]); 
    free_vec(W2[j]); free_vec(W3[j]); 
    free_vec(neta2[j]); 
//    free_mat(W3tmg[j]); free_vec(W3mg[j]); 
  }

  for (j=0;j<*Ntimes;j++) {
    free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]);free_mat(M1M2[j]);free_mat(ZXAIs[j]);
    free_vec(ZXdA[j]); free_mat(St[j]); free_mat(d2G[j]); free_vec(dG[j]);  
    free_vec(varUthat[j]);
  } 
//    free(ipers); 
    free(pg); free(cluster); 
    // }}}
    
}
