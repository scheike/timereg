#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
void transsurv(times,Ntimes,designX,nx,px,antpers,start,stop,betaS,Nit,cu,vcu,Iinv,
Vbeta,detail,sim,antsim,rani,Rvcu,RVbeta,test,testOBS,Ut,simUt,Uit,id,status,
weighted,ratesim,score,dhatMit,dhatMitiid,retur,loglike,profile,sym)
double *designX,*times,*betaS,*start,*stop,*cu,*Vbeta,*RVbeta,*vcu,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*score,*dhatMit,*dhatMitiid,*loglike;
int *nx,*px,*antpers,*Ntimes,*Nit,*detail,*sim,*antsim,*rani,*id,*status,*weighted,*ratesim,*retur,*profile,*sym;
{
  matrix *ldesignX,*ldesignG,*CtVUCt,*A,*AI;
  matrix *dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *tmp1,*tmp2,*tmp3,*dS1,*SI,*dS2,*S2,*S2pl,*dS2pl,*M1,*VU,*ZXAI,*VUI; 
  matrix *d2S0,*RobVbeta,*tmpM1,*Utt,*S1t,*S1start,*tmpM2,*et,*gt,*qt;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *dotwitowit[*antpers],*W3t[*antpers],*W4t[*antpers],*W2t[*antpers],*AIxit[*antpers],*Uti[*antpers],*d2G[*Ntimes],*Delta,*Delta2; 
  vector *Ctt,*lht,*S1,*dS0,*S0t,*S0start,*dA,*VdA,*dN,*MdA,*delta,*zav,*dlamt,*plamt,*dG[*Ntimes],
    *S1star;
  vector *xav,*difxxav,*xi,*zi,*U,*Upl,*beta,*xtilde; 
  vector *Gbeta,*zcol,*one,*difzzav,*difZ; 
  vector *offset,*weight,*ZXdA[*Ntimes],*varUthat[*Ntimes],*Uprofile;
  vector *ahatt,*risk,*tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB; 
  vector *W2[*antpers],*W3[*antpers],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  vector *dLamt[*antpers];
  int *pg=calloc(1,sizeof(int)),c,robust=1,pers=0,i,j,k,l,s,t,it,count,sing;
  double S0p,S0star,dtime,time,dummy,ll;
  double S0cox,S0,tau,random,scale,sumscore;
  int idum,*ipers=calloc(*Ntimes,sizeof(int)),nap;
  double norm_rand();
  void GetRNGstate(),PutRNGstate();

  pg[0]=1; 

  for (j=0;j<*antpers;j++) { 
    malloc_vec(*Ntimes,dLamt[j]); malloc_mat(*Ntimes,*px,dotwitowit[j]); 
    malloc_mat(*Ntimes,*pg,W3t[j]); malloc_mat(*Ntimes,*pg,W4t[j]); 
    malloc_mat(*Ntimes,*px,W2t[j]); malloc_mat(*Ntimes,*px,Uti[j]); 
    malloc_vec(*px,W2[j]); malloc_vec(*pg,W3[j]); malloc_mat(*Ntimes,*px,AIxit[j]); }
  malloc_vec(*Ntimes,S0t);  
  malloc_mat(*Ntimes,*pg,tmpM1); 
  malloc_mat(*Ntimes,*px,S1t); malloc_mat(*Ntimes,*px,tmpM2); 
  malloc_mat(*Ntimes,*px,S1start); 
  malloc_mat(*Ntimes,*px,et); malloc_mat(*Ntimes,*px,gt); malloc_mat(*Ntimes,*px,qt); 
  malloc_mat(*Ntimes,*px,Utt); 
  malloc_mat(*Ntimes,*pg,Delta); malloc_mat(*Ntimes,*px,Delta2); 

  malloc_vec(*Ntimes,S0start); malloc_vec(*Ntimes,lht); 
  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 

  idum=*rani; nap=floor(*antsim/50);

  malloc_mats(*antpers,*px,&ldesignX,NULL); 
  malloc_mats(*antpers,*pg,&ZP,&ldesignG,NULL); 
  malloc_mats(*px,*px,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*px,*px,&d2S0,&RobVbeta,&tmp1,&tmp2,&dS1,&S2,&dS2,&S2pl,&dS2pl,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*px,*pg,&tmp3,&ZPX,&dYI,&Ct,NULL); 

  malloc_vecs(*antpers,&risk,&weight,&plamt,&dlamt,&dN,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&Ctt,&ahatt,&tmpv1,&difX,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&S1,&dS0,&S1star,&xtilde,&xav,&difxxav,NULL); 
  malloc_vecs(*px,&U,&Upl,&beta,&delta,&zav,&difzzav,&Uprofile,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&difZ,&zav,&VdB,NULL); 

  for(j=0;j<*Ntimes;j++) {
    malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]);
    malloc_mat(*pg,*px,ZXAIs[j]); malloc_mat(*px,*pg,dYIt[j]); malloc_vec(*px,dAt[j]);
    malloc_vec(*pg,ZXdA[j]); malloc_mat(*px,*px,St[j]); 
    malloc_mat(*px,*px,d2G[j]); malloc_vec(*px,dG[j]); malloc_vec(*px,varUthat[j]);}

  ll=0; for(j=0;j<*px;j++) VE(beta,j)=betaS[j]; 

  cu[0]=times[0]; 
  for (it=0;it<*Nit;it++)
    {
      vec_zeros(U); vec_zeros(Upl); mat_zeros(S2pl); mat_zeros(S2); ll=0; sumscore=0; 
      mat_zeros(COV); 

      for (s=1;s<*Ntimes;s++)
	{
	  time=times[s]; vec_zeros(dN); sing=0; mat_zeros(ldesignX); vec_zeros(risk); 
	  vec_zeros(dS0);  mat_zeros(d2S0);  mat_zeros(dS1);  vec_zeros(S1star);
	  S0star=0; S0=0; S0p=0; 
	  vec_zeros(S1); S0cox=0; 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	    if ((start[c]<time) && (stop[c]>=time))  {
	      VE(risk,id[c])=1; 
	      for(j=0;j<*px;j++) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];  
	      if (time==stop[c] && status[c]==1) {VE(dN,id[c])=1; pers=id[c];} 
	      count=count+1; } 
	  }
	  ipers[s]=pers;
	  Mv(ldesignX,beta,Gbeta); 

	  for (j=0;j<*antpers;j++) { 
	    extract_row(ldesignX,j,xi); 
	    scale=(exp(-VE(Gbeta,j))+cu[1*(*Ntimes)+s-1]); 
	    dummy=1/scale; 
	    VE(plamt,j)=VE(risk,j)*dummy; 
	    VE(dlamt,j)=VE(risk,j)*dummy*dummy; 
	    S0star=S0star-VE(dlamt,j); 
	    S0=S0+VE(plamt,j);
	    S0p=S0p+VE(risk,j)/(exp(-VE(Gbeta,j))+cu[1*(*Ntimes)+s]);
	    S0cox=S0cox+exp(VE(Gbeta,j)); 

	    scl_vec_mult(-exp(-VE(Gbeta,j)),xi,tmpv1); vec_add(tmpv1,dG[s-1],dA); 
	    scl_vec_mult(-VE(plamt,j),dA,dA); 
	    if (*profile==0) 
	      replace_row(dotwitowit[j],s,xi); else replace_row(dotwitowit[j],s,dA); 
  
	    scl_vec_mult(VE(plamt,j),dA,tmpv1); vec_add(tmpv1,dS0,dS0); 
	    if (s<0 && j<5 ) {printf(" %d %d \n",s,j); print_vec(tmpv1); }

	    if (*profile==0) scl_vec_mult(-VE(dlamt,j),xi,dA); 
	    else scl_vec_mult(-VE(dlamt,j),dA,dA); 
	    vec_add(dA,S1star,S1star); 

	    for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
	      ME(dS1,i,k)=ME(dS1,i,k)+VE(xi,i)*VE(tmpv1,k); 

	    for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
	      ME(tmp1,i,k)=(VE(xi,i)*VE(xi,k))*exp(-VE(Gbeta,j));
	    mat_add(tmp1,d2G[s-1],tmp1);
	    scl_mat_mult(-VE(dlamt,j),tmp1,tmp1); 

	    for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
	      ME(d2S0,i,k)=ME(d2S0,i,k)+ME(tmp1,i,k)+2*scale*(VE(tmpv1,i)*VE(tmpv1,k)); 
   
	    scl_vec_mult(VE(plamt,j),xi,xi); vec_add(S1,xi,S1); 
	  }
	  VE(S0t,s)=S0; replace_row(S1t,s,dS0); 
	  VE(S0start,s)=S0star; replace_row(S1start,s,S1star); 
	  VE(lht,s)=VE(lht,s-1)-S0star/(S0*S0); 

	  /* printf(" %ld %lf %lf \n",s,VE(lht,s),ME(AI,0,0)); */

	  scl_vec_mult(S0star,dS0,tmpv1); scl_vec_mult(S0,S1star,dA); 
	  vec_subtr(tmpv1,dA,dA); scl_vec_mult(1/S0,dA,dA); 
	  replace_row(gt,s,dA); 


	  scl_vec_mult(-1/(S0*S0),dS0,tmpv1); vec_add(dG[s-1],tmpv1,dG[s]); 
	  if (s<0) { printf(" %lf \n",S0); 
	    print_vec(scl_vec_mult(1/S0,dS0,NULL)); 
	    print_vec(tmpv1); 
	    print_vec(dG[s]); }

	  scl_mat_mult(-1/(S0*S0),d2S0,A); 
	  for (i=0;i<*px;i++) for (j=0;j<*px;j++) 
	    ME(A,i,j)=ME(A,i,j)+2*S0*VE(tmpv1,i)*VE(tmpv1,j);
	  mat_add(d2G[s-1],A,d2G[s]); 

	  /* baseline is computed */ 
	  cu[1*(*Ntimes)+s]=cu[1*(*Ntimes)+s-1]+(1/S0); 
	  if (s<0) printf(" %lf \n",cu[1*(*Ntimes)+s]); 

	  /* First derivative of U ========================================  */ 
	  extract_row(ldesignX,pers,xi); 

	  scl_vec_mult(1/S0,S1,xav); vec_subtr(xi,xav,difxxav); vec_add(U,difxxav,U); 
	  if (*profile==0)  replace_row(et,s,xav); 

	  /* profile version of score */ 
	  dummy=1/(exp(-VE(Gbeta,pers))+cu[1*(*Ntimes)+s-1]); 
	  scl_vec_mult(-exp(-VE(Gbeta,pers)),xi,tmpv1); 
	  vec_add(tmpv1,dG[s-1],tmpv1); scl_vec_mult(-dummy,tmpv1,tmpv1); 

	  scl_vec_mult(1/S0,dS0,dA); 
	  if (*profile==1)  replace_row(et,s,dA); 
	  vec_subtr(tmpv1,dA,dA); 
	  vec_add(Upl,dA,Upl);

	  /* Second derivative S =========================================== */ 

	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
	    ME(dS2pl,i,k)=(VE(xi,i)*VE(xi,k))*exp(-VE(Gbeta,pers)); 
	  mat_add(dS2pl,d2G[s-1],dS2pl);
	  scl_mat_mult(-dummy,dS2pl,dS2pl); 

	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
	    ME(dS2pl,i,k)=ME(dS2pl,i,k)+(VE(tmpv1,i)*VE(tmpv1,k)); 

	  scl_mat_mult(-1/S0,d2S0,A); 
	  scl_vec_mult(1/S0,dS0,dA); 
	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
	    ME(A,i,k)=ME(A,i,k)+(VE(dA,i)*VE(dA,k)); 

	  mat_add(A,dS2pl,dS2pl);
	  mat_add(dS2pl,S2pl,S2pl); 

	  if (*profile==1) St[s]=mat_copy(S2pl,St[s]); 

	  /* simple Second derivative S2 ================================== */
	  for (i=0;i<*px;i++) for (k=0;k<*px;k++) 
	    ME(dS2,i,k)=VE(dS0,i)*VE(S1,k); 
	  scl_mat_mult(S0,dS1,M1); 

	  /* */
	  if (s<0) { printf("======================== %lf \n",S0); 
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

	  if (it==((*Nit)-1)) { 
	    Ut[s]=time; 
	    for (i=1;i<*px+1;i++) 
	      if (*profile==1) Ut[i*(*Ntimes)+s]=VE(Upl,i-1);  
	      else Ut[i*(*Ntimes)+s]=VE(U,i-1);  
	    for (i=1;i<*px+1;i++) ME(Utt,s,i-1)=Ut[i*(*Ntimes)+s]; 
	    for (j=0;j<*antpers;j++) VE(dLamt[j],s)=VE(plamt,j)/S0;
	    /*
	      for (i=0;i<*px;i++) for (j=0;j<*pg;j++) ME(dM1M2,j,i)=VE(dA,i)*VE(difzzav,j);
	      for (i=0;i<*pg;i++) 
	      for (j=0;j<*pg;j++) ME(VU,i,j)=ME(VU,i,j)+VE(difzzav,i)*VE(difzzav,j); 

	      MxA(AI,ZPX,dYIt[s]); mat_subtr(Ct,dYIt[s],Ct); C[s]=mat_copy(Ct,C[s]); 

	      vec_star(dA,dA,VdA); mat_add(dM1M2,M1M2t,M1M2t); 
	      M1M2[s]=mat_copy(M1M2t,M1M2[s]); 

	      for (k=1;k<=*px;k++) { cu[k*(*Ntimes)+s]=VE(dA,k-1);
	      vcu[k*(*Ntimes)+s]=VE(VdA,k-1)+vcu[k*(*Ntimes)+s-1]; }
	    */
	  }
	} /* s= .... Ntimes */ 

      if (*profile==1)  scl_mat_mult(-1,S2pl,A); else scl_mat_mult(-1,S2,A);
      invert(A,SI);
      if (*profile==1) Mv(SI,Upl,delta);  else Mv(SI,U,delta); 
      vec_add(beta,delta,beta); 

      if (*detail==1) { 
	printf("====================Iteration %d ==================== \n",it); 
	printf("log-likelihood "); printf(" %lf \n",ll); 
	printf("Estimate beta "); print_vec(beta); 
	if (*profile==1) {
	  printf("modified partial likelihood Score D l"); print_vec(Upl); }
	if (*profile==0) {printf("simple Score D l"); print_vec(U);  }
	printf("Information -D^2 l\n"); print_mat(SI); 
	if (*profile==1) {printf("simple D2 l");  print_mat(S2pl); }
	if (*profile==0) {printf("simple D2 l");  print_mat(S2); }
      };

      for (k=0;k<*px;k++) sumscore= sumscore+
	(*profile==1)*fabs(VE(Upl,k))+(*profile==0)*fabs(VE(U,k));  

      if ((fabs(sumscore)<0.000001) & (it<*Nit-2)) it=*Nit-2; 
    } /* it */
  loglike[0]=ll; 

   R_CheckUserInterrupt();


  /* computation of q(t) ===================================== */ 
  for (s=1;s<*Ntimes;s++) {
    vec_zeros(xi); 
    for (t=s;t<*Ntimes;t++) {
      extract_row(gt,t,dA); scl_vec_mult(exp(VE(lht,t))/VE(S0t,t),dA,dA); 
      if (s<0) {printf("exp %d %d %lf \n",s,t,exp(-VE(lht,t)+VE(lht,s))); 
	print_vec(dA);}
      vec_add(dA,xi,xi); }
    scl_vec_mult(exp(-VE(lht,s))/VE(S0t,s),xi,xi); 
    replace_row(qt,s,xi); 
  }

   R_CheckUserInterrupt();
  /* terms for robust variances ============================ */
  if (robust==1) {
    for (s=1;s<*Ntimes;s++) {
      time=times[s]; vec_zeros(dN);dtime=time-times[s-1]; 
      cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; Ut[s]=times[s]; 

      extract_row(qt,s,tmpv1); extract_row(et,s,xtilde); 

      for (i=0;i<*antpers;i++) {
	extract_row(dotwitowit[i],s,rowX); vec_subtr(rowX,xtilde,rowX); 
	if (s==0) { print_vec(rowX); print_vec(tmpv1); }
	vec_add(rowX,tmpv1,rowX); 
	if (s==0) {print_vec(rowX);}

	if (i==ipers[s]) for (j=0;j<*px;j++) for (k=0;k<*px;k++)
	  ME(VU,j,k)=ME(VU,j,k)+VE(rowX,j)*VE(rowX,k);

	scl_vec_mult(VE(dLamt[i],s),rowX,xi); 

	vec_subtr(W2[i],xi,W2[i]); 
	if (i==ipers[s]) vec_add(rowX,W2[i],W2[i]);
	if (*ratesim==1) {scl_vec_mult(VE(dLamt[i],s),tmpv2,rowZ); vec_subtr(W2[i],rowZ,W2[i]);}
	replace_row(W2t[i],s,W2[i]); 

	VE(rowZ,0)=exp(-VE(lht,s))/VE(S0t,s); 

	scl_vec_mult(VE(dLamt[i],s),rowZ,zi); 
	vec_subtr(W3[i],zi,W3[i]); 
	if (i==ipers[s]) vec_add(rowZ,W3[i],W3[i]);

	if (*ratesim==1) {scl_vec_mult(VE(dLamt[i],s),rowX,rowX); 
		          vec_subtr(W3[i],rowX,W3[i]);}
	replace_row(W3t[i],s,W3[i]);  

	if (*retur==1)  dhatMit[i*(*Ntimes)+s]=1*(i==pers)-VE(dLamt[i],s);
      } /* i=1..antpers */ 

      /* Martingale baseret variance */
      /*
	MxA(C[s],VU,tmp3); MAt(tmp3,C[s],CtVUCt);
	MxA(C[s],SI,tmp3); MxA(tmp3,M1M2[s],COV); 

	for (k=1;k<=*px;k++) {
	cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+cu[k*(*Ntimes)+s]; 
	vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s]+ME(CtVUCt,k-1,k-1)
	+2*ME(COV,k-1,k-1); }
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

      for (i=0;i<*antpers;i++) {

	Mv(SI,W2[i],tmpv1); vM(Ct,tmpv1,rowZ);

	extract_row(W3t[i],s,zi); 
	VE(zi,0)= exp(VE(lht,s))*VE(zi,0);

	vec_add(rowZ,zi,zi); 
	if (i==-5) print_vec(zi); 
	replace_row(W4t[i],s,zi);
	vec_star(zi,zi,rowZ); vec_add(rowZ,VdB,VdB);

	extract_row(W2t[i],s,xi); Mv(St[s],tmpv1,rowX); 
	vec_add(xi,rowX,tmpv1); replace_row(Uti[i],s,tmpv1);
  
	vec_star(tmpv1,tmpv1,xi); vec_add(xi,varUthat[s],varUthat[s]);

	if (s==1) { 
	  for (j=0;j<*px;j++) for (k=0;k<*px;k++)
	    ME(RobVbeta,j,k)=ME(RobVbeta,j,k)+VE(W2[i],j)*VE(W2[i],k); }

	if (*retur==1) {

	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignX,j,xi); extract_row(ldesignG,j,zi);
	      dummy=exp(VE(Gbeta,j))*VE(weight,j)*VE(offset,j); 
	      scl_vec_mult(dummy,xi,xtilde); replace_row(ldesignX,j,xtilde); }

	  Mv(ldesignX,dAt[s],dlamt);  
	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignG,j,zi); scl_vec_mult(VE(dlamt,j),zi,zi); replace_row(ZP,j,zi);} 

	  Mv(ZP,W2[i],reszpbeta);
	  Mv(dYIt[s],W2[i],xi); Mv(ldesignX,xi,res1dim); 

	  dhatMitiid[i*(*Ntimes)+s]=dhatMit[i*(*Ntimes)+s]-
	    (VE(reszpbeta,0)- VE(res1dim,0)); 
	} /* retur ==1 */ 

      } /* i =1 ..Antpers */
      for (k=1;k<*pg+1;k++) Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      for (k=1;k<*pg+1;k++) vcu[k*(*Ntimes)+s]=VE(VdB,k-1);
    }  /*  s=1 ..Ntimes */ 

    MxA(RobVbeta,SI,tmp1); MxA(SI,tmp1,RobVbeta);
  } /* Robust =1 , default */ 


  MxA(VU,SI,tmp1); MxA(SI,tmp1,VU);

  for(j=0;j<*px;j++) { 
    betaS[j]= VE(beta,j); 
    if (*profile==1) score[j]=VE(Upl,j); else score[j]=VE(U,j); 
    for (k=0;k<*px;k++){ Iinv[k*(*px)+j]=ME(SI,j,k);
      Vbeta[k*(*px)+j]=-ME(VU,j,k); 
      RVbeta[k*(*px)+j]=-ME(RobVbeta,j,k); } } 

   R_CheckUserInterrupt();

  if (*sim==1) {
    // printf("Simulations start N= %d \n",*antsim);
    GetRNGstate();  /* to use R random normals */

    tau=times[*Ntimes-1]-times[0];
    for (i=1;i<=*pg;i++) VE(rowZ,i-1)=cu[i*(*Ntimes)+(*Ntimes-1)];

    /* Beregning af OBS teststørrelser */
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
    } /*s=1..Ntimes Beregning af obs teststørrelser */


    for (k=1;k<=*antsim;k++) {
      mat_zeros(Delta); mat_zeros(Delta2);
      for (i=0;i<*antpers;i++) {
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
	  if (VE(difZ,i)>test[l*(*antsim)+k]) test[l*(*antsim)+k]=VE(difZ,i);
	  VE(zi,i)=fabs(ME(Delta,s,i))/sqrt(Rvcu[(i+1)*(*Ntimes)+s]);
	  if (VE(zi,i)>test[i*(*antsim)+k]) test[i*(*antsim)+k]=VE(zi,i); }

	if (*weighted>=1) {
	  extract_row(Delta2,s,xi); 
	  if ((s>*weighted) && (s<*Ntimes-*weighted))  {
	    for (i=0;i<*px;i++) {VE(xi,i)=fabs(ME(Delta2,s,i))/sqrt(VE(varUthat[s],i));
	      if (VE(xi,i)>simUt[i*(*antsim)+k]) simUt[i*(*antsim)+k]=VE(xi,i); }

	    if (k<50) { 
	      for (i=0;i<*px;i++) { l=(k-1)*(*px)+i;
		Uit[l*(*Ntimes)+s]=ME(Delta2,s,i)/sqrt(VE(varUthat[s],i));}}
	  } 
	} /* weigted score */
	else {
	  extract_row(Delta2,s,xi); 
	  for (i=0;i<*px;i++) {
	    if (fabs(VE(xi,i))>simUt[i*(*antsim)+k]) 
	      simUt[i*(*antsim)+k]=fabs(VE(xi,i)); }

	  if (k<50) { 
	    for (i=0;i<*px;i++) { l=(k-1)*(*px)+i;
	      Uit[l*(*Ntimes)+s]=ME(Delta2,s,i);}}
	} /* else wscore=0 */ 

      }  /* s=1..Ntims */
    }  /* k=1..antsim */
    PutRNGstate();  /* to use R random normals */
  } /* sim==1 */


  /*
    free_mats(&Delta,&Delta2,&et,&gt,&qt,&S1start,&d2S0,&Utt,&tmpM2,&VU,&ZX,&COV,&dM1M2,&AI,&A,&ZXAI,&tmp1,&tmp2,&tmp3,&ldesignX,&ldesignG,&M1,&dS1,&SI,&dS2,&S2,&S2pl,&dS2pl,&VUI,&ZP,&ZPX,&dYI,&Ct,&M1M2t,&RobVbeta,&S1t,&tmpM1,&CtVUCt,NULL); 
    free_vecs(&Ctt,&risk,&difZ,&S1star,&lht,&S0start,&S0t,&S1,&dS0,&xav,&difxxav,&ta,&ahatt,&Uprofile,&plamt,&dlamt,&dN,&one,&xi,&zcol,&Gbeta,&VdA,&dA,&MdA,&xtilde,&zi,&U,&Upl,&beta,&delta,&zav,&difzzav,&weight,&offset,&tmpv1,&tmpv2,&rowX,&rowZ,&difX,&VdB,&reszpbeta,&res1dim,NULL); 
  */

  for (j=0;j<*antpers;j++) {
    free_vec(dLamt[j]); free_mat(dotwitowit[j]); free_mat(W3t[j]); free_mat(W4t[j]); 
    free_mat(W2t[j]);free_vec(W2[j]); free_vec(W3[j]); free_mat(AIxit[j]); free_mat(Uti[j]); }
  for (j=0;j<*Ntimes;j++) {
    free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]);free_mat(M1M2[j]);free_mat(ZXAIs[j]);
    free_vec(ZXdA[j]); free_mat(St[j]); free_mat(d2G[j]); free_vec(dG[j]);  
    free_vec(varUthat[j]);} 
    free(ipers); 
    free(pg); 
}
