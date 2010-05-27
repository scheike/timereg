#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
/* ====================================================== */
void score(times,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,
betaS,Nit,cu,vcu,w,mw,loglike,Iinv,Vbeta,detail,offs,mof,sim,antsim,
rani,Rvcu,RVbeta,test,testOBS,Ut,simUt,Uit,XligZ,aalen,nb,id,status,wscore,ridge,ratesim,score,dhatMit,gammaiid,dmgiid,retur,robust,covariance,Vcovs,addresamp,addproc,
resample,gamiid,biid,clusters,antclust,vscore,betafixed)
double
*designX,*designG,*times,*betaS,*start,*stop,*cu,*w,*loglike,*Vbeta,*RVbeta,*vcu,*offs,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*aalen,*ridge,*score,*dhatMit,*gammaiid,*dmgiid,*Vcovs,*addproc,*gamiid,*biid,*vscore;
int*covariance,*nx,*px,*ng,*pg,*antpers,*Ntimes,*mw,*Nit,*detail,*mof,*sim,*antsim,*rani,*XligZ,*nb,*id,*status,*wscore,*ratesim,*retur,*robust,*addresamp,*resample,*clusters,*antclust,*betafixed;
{ 
// {{{ setting up memory 
  matrix *ldesignX,*ldesignG,*cdesX,*cdesX2,*cdesX3,*CtVUCt,*A,*AI;
  matrix *Vcov,*dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *tmp1,*tmp2,*tmp3,*dS,*S1,*SI,*S2,*M1,*VU,*ZXAI,*VUI; 
  matrix *RobVbeta,*Delta,*tmpM1,*Utt,*Delta2,*tmpM2;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *W3t[*antclust],*W4t[*antclust],*W2t[*antclust],*AIxit[*antpers],*Uti[*antclust]; 
  vector *dA,*VdA,*MdA,*delta,*zav,*lamt,*lamtt;
  vector *xi,*zi,*U,*beta,*xtilde; 
  vector *Gbeta,*zcol,*one,*difzzav;
  vector *offset,*weight,*ZXdA[*Ntimes],*varUthat[*Ntimes],*Uprofile;
  vector *ta,*ahatt,*vrisk; 
  vector *tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB; 
  vector *W2[*antclust],*W3[*antclust],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  int ci,c,pers=0,i,j,k,l,s,it,count,sing,pmax,
      *imin=calloc(1,sizeof(int)),
      *cluster=calloc(*antpers,sizeof(int));
  double dtime,time=0,dummy,ll,lle,llo;
  double tau,hati,random,scale,sumscore;
  int idum,*ipers=calloc(*Ntimes,sizeof(int)),nap;
  double norm_rand();
  void GetRNGstate(),PutRNGstate();

  /* float gasdev(),expdev(),ran1(); */ 

  GetRNGstate();  /* to use R random normals */

  if (*robust==1) {
    for (j=0;j<*antclust;j++) { malloc_mat(*Ntimes,*px,W3t[j]); 
      malloc_mat(*Ntimes,*px,W4t[j]); malloc_mat(*Ntimes,*pg,W2t[j]); 
      malloc_mat(*Ntimes,*pg,Uti[j]); malloc_vec(*pg,W2[j]); malloc_vec(*px,W3[j]); }
    for (j=0;j<*antpers;j++) {malloc_mat(*Ntimes,*px,AIxit[j]);
                              cluster[j]=0; }
    for(j=0;j<*Ntimes;j++) malloc_vec(*pg,varUthat[j]);
  }

  if (*sim==1) {
    malloc_mat(*Ntimes,*px,Delta);  malloc_mat(*Ntimes,*px,tmpM1); 
    malloc_mat(*Ntimes,*pg,Delta2); malloc_mat(*Ntimes,*pg,tmpM2);  }
  malloc_mat(*Ntimes,*pg,Utt); 

  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 

  idum=*rani; nap=floor(*antsim/50);

  malloc_mats(*antpers,*px,&ldesignX,&cdesX,&cdesX2,&cdesX3,NULL); 
  malloc_mats(*px,*antpers,NULL); 
  malloc_mats(*antpers,*pg,&ZP,&ldesignG,NULL); 
  malloc_mats(*px,*px,&Vcov,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*pg,*pg,&RobVbeta,&tmp1,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*px,*pg,&tmp3,&ZPX,&dYI,&Ct,NULL); 

  malloc_vecs(*antpers,&weight,&lamtt,&lamt,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&xtilde,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,NULL); 

  for(j=0;j<*Ntimes;j++) { malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]); 
    malloc_mat(*pg,*px,ZXAIs[j]); malloc_mat(*px,*pg,dYIt[j]); malloc_vec(*px,dAt[j]);
    malloc_vec(*pg,ZXdA[j]); malloc_mat(*pg,*pg,St[j]); } 
  malloc_vec(*nb,ta); 
  malloc_vec(*antpers,vrisk); 

  if (*px>=*pg) pmax=*px; else pmax=*pg; ll=0; 
  for(j=0;j<*pg;j++) VE(beta,j)=betaS[j]; 
  for(j=0;j<*antpers;j++) {VE(one,j)=1; VE(weight,j)=1; VE(offset,j)=1;} 
  // }}}
  
  cu[0]=times[0]; 
  for (it=0;it<*Nit;it++) // {{{ iterations start for cox-aalen model
    {
      vec_zeros(U); mat_zeros(S1);  sumscore=0; 
      for (s=1;s<*Ntimes;s++)
	{
	  time=times[s]; sing=0; mat_zeros(ldesignX); mat_zeros(ldesignG); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	    {
	      if ((start[c]<time) && (stop[c]>=time))  {
		cluster[id[c]]=clusters[c];
		for(j=0;j<pmax;j++) {
		  if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		  if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; } 
		if (time==stop[c] && status[c]==1) {pers=id[c];} 
		if (*mof==1) VE(offset,id[c])=offs[c];  
		if (*mw==1) VE(weight,id[c])=w[c]; 

		count=count+1; }
	    }
	  ipers[s]=pers;

	  if (aalen[0]!=1) {printf("aalen \n"); 
	    for(j=0;j<*nb;j++) VE(ta,j)=fabs(aalen[j]-time); dummy=vec_min(ta,imin);
	    for(j=1;j<=*px;j++) VE(ahatt,j-1)=aalen[j*(*nb)+(*imin)];
	    Mv(ldesignX,ahatt,lamtt);
	    for (j=0;j<*antpers;j++) if (VE(lamtt,j)!=0) VE(weight,j)=1/VE(lamtt,j);}

	  Mv(ldesignG,beta,Gbeta); 

	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignX,j,xi); 
	      dummy=exp(VE(Gbeta,j));           /* W D(exp(z*beta))X */ 
	      scl_vec_mult(VE(weight,j)*dummy,xi,xtilde); replace_row(cdesX,j,xtilde); 
	    }

	  scale=VE(weight,pers); 
  
	  MtA(cdesX,ldesignX,A); invert(A,AI); 
	  if (ME(AI,0,0)==0) printf(" X'X not invertible at time %lf \n",time);
	  if (s==0) {
	    print_mat(cdesX); 
	    print_mat(ldesignX); 
	    print_vec(Gbeta); print_mat(A); print_mat(AI); 
	  }

	  extract_row(ldesignX,pers,xi); scl_vec_mult(scale,xi,xi); 
	  Mv(AI,xi,dA); MtA(ldesignG,cdesX,ZX); 
	  MxA(ZX,AI,ZXAIs[s]); 
	  Mv(ZX, dA, ZXdA[s]);  scl_vec_mult(1,dA,dAt[s]); 

	  if (s<1) {print_mat(A); print_mat(AI); print_mat(ZX); print_mat(ZXAIs[s]); }

	  /* First derivative U and Second derivative S  */ 
  
	  extract_row(ldesignG,pers,zi); scl_vec_mult(scale,zi,zi); 
	  Mv(ZX, dA, zav);  vec_subtr(zi,zav,difzzav); vec_add(difzzav,U,U); 

	  if (s<1) {print_vec(zi); print_vec(zav); print_vec(difzzav);}

	  Mv(cdesX,dA,lamt);  
	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignG,j,zi); 
	      scl_vec_mult(VE(lamt,j),zi,zi); replace_row(ZP,j,zi);} 

	  MtA(ldesignX,ZP,ZPX); MxA(ZXAIs[s],ZPX,tmp2); 

	  MtA(ZP,ldesignG, tmp1); mat_subtr( tmp1,tmp2, dS); 
	  mat_add(dS,S1,S1);  St[s]=mat_copy(S1,St[s]);

	  /* varians beregninger */ 
	  if (it==((*Nit)-1)) { 
	    replace_row(Utt,s,U);
	    vscore[0*(*Ntimes)+s]=times[s]; 

	    for (i=0;i<*px;i++) for (j=0;j<*pg;j++) ME(dM1M2,j,i)=VE(dA,i)*VE(difzzav,j);
	    for (i=0;i<*pg;i++) { 
	      if (*betafixed==1) 
		vscore[(i+1)*(*Ntimes)+s]= vscore[(i+1)*(*Ntimes)+s-1]+VE(difzzav,i)*VE(difzzav,i); 
	      for (j=0;j<*pg;j++) ME(VU,i,j)=ME(VU,i,j)+VE(difzzav,i)*VE(difzzav,j); }

	    MxA(AI,ZPX,dYIt[s]); mat_subtr(Ct,dYIt[s],Ct); C[s]=mat_copy(Ct,C[s]); 

	    vec_star(dA,dA,VdA); mat_add(dM1M2,M1M2t,M1M2t); M1M2[s]=mat_copy(M1M2t,M1M2[s]); 

	    for (k=1;k<=*px;k++) {cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dA,k-1); 
	      vcu[k*(*Ntimes)+s]=VE(VdA,k-1)+vcu[k*(*Ntimes)+s-1];}
	  }
	  if (*robust==1) 
	    for (j=0;j<*antpers;j++)
	      {extract_row(ldesignX,j,xi); Mv(AI,xi,rowX); replace_row(AIxit[j],s,rowX);}
	} /* Ntimes */ 

      /* for (k=0;k<*pg;k++) ME(S1,k,k)=ME(S1,k,k)+*ridge;  */
      invert(S1,SI); 

      Mv(SI,U,delta); if (*betafixed==0) vec_add(beta,delta,beta); 
      MxA(SI,VU,S2); MxA(S2,SI,VU); 

      if (*detail==1) { 
	printf("====================Iteration %ld ==================== \n",(long int) it);
	printf("Estimate beta \n"); print_vec(beta); 
	printf("Score D l\n"); print_vec(U); 
	printf("Information -D^2 l\n"); print_mat(SI); };

      for (k=0;k<*pg;k++) sumscore= sumscore+fabs(VE(U,k)); 

      if ((fabs(sumscore)<0.0000000001) & (it<*Nit-2)) it=*Nit-2; 
    } /* it */ // }}}

  for (k=0;k<*pg;k++) score[k]=VE(U,k); 
  lle=0; llo=0;
  for (s=1;s<*Ntimes;s++) { // {{{ terms for robust variances 
    time=times[s]; dtime=time-times[s-1]; 
    cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; Ut[s]=times[s]; 
    vec_zeros(vrisk); 

    mat_zeros(ldesignX); mat_zeros(ldesignG); 
    for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
      {
	if ((start[c]<time) && (stop[c]>=time))  { 
	  cluster[id[c]]=clusters[c];
	  VE(vrisk,id[c])=1.0; 
	  for(j=0;j<pmax;j++) {
	    if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	    if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; } 
	  if (time==stop[c] && status[c]==1) {pers=id[c];} 
	  if (*mof==1) VE(offset,id[c])=offs[c];  
	  if (*mw==1) VE(weight,id[c])=w[c]; 
	  count=count+1; }
      }
    Mv(ldesignG,beta,Gbeta); 

    for (j=0;j<*antpers;j++)
      {extract_row(ldesignX,j,xi); extract_row(ldesignG,j,zi);
	dummy=exp(VE(Gbeta,j))*VE(weight,j)*VE(offset,j); 
	scl_vec_mult(dummy,xi,xtilde); replace_row(cdesX,j,xtilde); 
	if (j==pers) {vec_star(xtilde,dAt[s],tmpv1); hati=vec_sum(tmpv1); 
	  lle=lle+log(hati);}
      }


    /* terms for robust variance   */ 
    if (*robust==1 || *retur==1 ) {
	for (i=0;i<*antpers;i++)  
	{
         ci=cluster[i]; 
	 extract_row(cdesX,i,rowX); scl_vec_mult(1/VE(weight,i),rowX,rowX); 
	 extract_row(ldesignG,i,zi); extract_row(ldesignX,i,xi); 
	 vec_star(rowX,dAt[s],tmpv1); hati=vec_sum(tmpv1); 

	 if (*robust==1) {
	 Mv(ZXAIs[s],xi,tmpv2);  vec_subtr(zi,tmpv2,tmpv2); 
	 scl_vec_mult(VE(weight,i),tmpv2,tmpv2); 

	 if (i==pers) vec_add(tmpv2,W2[ci],W2[ci]);
	 if (*ratesim==1) {scl_vec_mult(hati,tmpv2,rowZ); vec_subtr(W2[ci],rowZ,W2[ci]); }

	 extract_row(AIxit[i],s,rowX);
	 scl_vec_mult(VE(weight,i),rowX,rowX); 

	 if (i==pers) {vec_add(rowX,W3[ci],W3[ci]); }
	 llo=llo+hati;

	 if (*ratesim==1) {scl_vec_mult(hati,rowX,rowX); vec_subtr(W3[ci],rowX,W3[ci]);}
	 }

	 if (*retur==1)  { 
	    dhatMit[i*(*Ntimes)+s]= dhatMit[i*(*Ntimes)+s]+1*(i==pers)-hati;
	      /*if ((i==pers) | (hati>0)) 
		printf(" %ld %ld %ld %lf \n",i,s,(i==pers),hati); */
	 }

       } /* i 1.. antpers */

       if (*robust==1) {
         for (j=0;j<*antclust;j++) 
  	 {
	  replace_row(W2t[j],s,W2[j]); replace_row(W3t[j],s,W3[j]);  
	 } /* j and i=1..antclust */ 
       }
    }


    /* MG baseret varians beregning */
    MxA(C[s],VU,tmp3); MAt(tmp3,C[s],CtVUCt);
    MxA(C[s],SI,tmp3); MxA(tmp3,M1M2[s],COV); 

    for (k=1;k<=*px;k++) {
      if (*betafixed==0) 
	vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s]+ME(CtVUCt,k-1,k-1)
	  +2*ME(COV,k-1,k-1); else vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s];
    }

    for (k=1;k<=*pg;k++) Ut[k*(*Ntimes)+s]=ME(Utt,s,k-1);
    /* */
  } /* s=1 ..Ntimes */  // }}}


  ll=lle-llo; /* likelihood beregnes */
  if (*detail==1) printf("loglike er %lf \n",ll);  

  if (*robust==1) // {{{ robust variances 
  {
      for (s=1;s<*Ntimes;s++) {
	vec_zeros(VdB); mat_zeros(Vcov);

	for (j=0;j<*antclust;j++) {

	  Mv(SI,W2[j],tmpv2); 
	  if (s==-1) {scl_vec_mult(1,tmpv2,W2[j]);}

	  Mv(C[s],tmpv2,rowX);
	  extract_row(W3t[j],s,tmpv1); vec_add(tmpv1,rowX,difX); 
	  if (*betafixed==1) scl_vec_mult(1,tmpv1,difX); 
	  replace_row(W4t[j],s,difX); vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	  if (*resample==1) {
	    if (s==1) 
	      if (*betafixed==0) {
		for (c=0;c<*pg;c++) gamiid[c*(*antclust)+j]=
		  gamiid[c*(*antclust)+j]+VE(tmpv2,c); }
	    for (c=0;c<*px;c++) {l=j*(*px)+c; 
	      biid[l*(*Ntimes)+s]=biid[l*(*Ntimes)+s]+VE(difX,c);} }

	  if (*covariance==1) {
	    for (k=0;k<*px;k++) for (c=0;c<*px;c++) 
	      ME(Vcov,k,c)=ME(Vcov,k,c)+VE(difX,k)*VE(difX,c);}

	  Mv(St[s],tmpv2,rowZ); extract_row(W2t[j],s,tmpv2); 
	  if (*betafixed==0) {
	    vec_subtr(tmpv2,rowZ,zi); replace_row(Uti[j],s,zi);} else replace_row(Uti[j],s,tmpv2);

	  vec_star(zi,zi,tmpv2); vec_add(tmpv2,varUthat[s],varUthat[s]);

	  if (s==1) { for (c=0;c<*pg;c++) for (k=0;k<*pg;k++)
			ME(RobVbeta,c,k)=ME(RobVbeta,c,k)+VE(W2[j],c)*VE(W2[j],k);}

	} /* j in clusters  */

	if (*betafixed==0) 
	  for (i=0;i<*pg;i++) vscore[(i+1)*(*Ntimes)+s]=VE(varUthat[s],i);  

	if (*retur==-1) {
	  mat_zeros(ldesignX); mat_zeros(ldesignG); 
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	    {
	      if ((start[c]<time) && (stop[c]>=time))  {
		cluster[id[c]]=clusters[c];
		for(j=0;j<pmax;j++) {
		  if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		  if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; } 
		if (time==stop[c] && status[c]==1) {pers=id[c];} 
		if (*mof==1) VE(offset,id[c])=offs[c];  
		if (*mw==1) VE(weight,id[c])=w[c]; 
		count=count+1; }
	    }
	  Mv(ldesignG,beta,Gbeta); 

	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignX,j,xi); extract_row(ldesignG,j,zi);
	      dummy=exp(VE(Gbeta,j))*VE(weight,j)*VE(offset,j); 
	      scl_vec_mult(dummy,xi,xtilde); replace_row(cdesX,j,xtilde); }

	  Mv(cdesX,dAt[s],lamt);  

	  for (i=0;i<*antpers;i++) {
	    extract_row(ldesignG,i,zi); scl_vec_mult(VE(lamt,i),zi,zi); 
	    Mv(SI,W2[cluster[i]],rowZ);  
	    vec_star(zi,rowZ,tmpv2); 
	    VE(reszpbeta,0)=vec_sum(tmpv2); Mv(dYIt[s],rowZ,xi); 
	    extract_row(cdesX,i,xtilde); 
	    vec_star(xtilde,xi,tmpv1); VE(res1dim,0)=vec_sum(tmpv1); 

	    dmgiid[i*(*Ntimes)+s]=dhatMit[i*(*Ntimes)+s]+
	      VE(weight,i)*(VE(reszpbeta,0)- VE(res1dim,0)); 
	  }
	} /* retur == -1 */ 


	for (k=1;k<*px+1;k++) { Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
	  if (*covariance==1) {
	    for (j=0;j<*px;j++)  {
	      l=(k-1)*(*px)+j; 
	      Vcovs[l*(*Ntimes)+s]=ME(Vcov,k-1,j); } } }

      }  /*  s=1 ..Ntimes */ 

      MxA(RobVbeta,SI,tmp1); MxA(SI,tmp1,RobVbeta);

      if (*retur==1) {
	for (i=0;i<*antclust;i++) {
	  Mv(SI,W2[i],rowZ);  
	  for (j=0;j<*pg;j++) gammaiid[i*(*pg)+j]=VE(rowZ,j); }
      }
    } /* if robust==1 */  // }}}

  for(j=0;j<*pg;j++) { betaS[j]= VE(beta,j); 
    loglike[0]=lle; loglike[1]=ll;
    for (k=0;k<*pg;k++){ Iinv[k*(*pg)+j]=ME(SI,j,k);
      Vbeta[k*(*pg)+j]=-ME(VU,j,k); 
      RVbeta[k*(*pg)+j]=-ME(RobVbeta,j,k); } } 

  if (*sim==1) { // {{{ score process simulations
    // printf("Simulations start N= %ld \n",(long int) *antsim);

    tau=times[*Ntimes-1]-times[0];
    for (i=1;i<=*px;i++) VE(rowX,i-1)=cu[i*(*Ntimes)+(*Ntimes-1)];

    /* Beregning af OBS teststørrelser */
    for (s=1;s<*Ntimes;s++) { time=times[s]-times[0]; 

      for (i=1;i<=*px;i++) {
	VE(xi,i-1)=fabs(cu[i*(*Ntimes)+s])/sqrt(Rvcu[i*(*Ntimes)+s]);
	if (VE(xi,i-1)>testOBS[i-1]) testOBS[i-1]=VE(xi,i-1); }

      scl_vec_mult(time/tau,rowX,difX);
      for (i=1;i<=*px;i++) VE(xi,i-1)=cu[i*(*Ntimes)+s];
      vec_subtr(xi,difX,difX);
      for (i=0;i<*px;i++) {
	VE(difX,i)=fabs(VE(difX,i)); l=(*px+i);
	if (VE(difX,i)>testOBS[l]) testOBS[l]=VE(difX,i);}

      if (*wscore>=1) {  /* sup beregnes i R */ 
	if ((s>*wscore) && (s<*Ntimes-*wscore)) {extract_row(Utt,s,rowZ);
	  for (i=0;i<*pg;i++) VE(rowZ,i) = VE(rowZ,i)/sqrt(VE(varUthat[s],i));
	  replace_row(Utt,s,rowZ); /* scaled score process */ 
	}
	else {vec_zeros(rowZ); replace_row(Utt,s,rowZ);}
      }
      for (k=1;k<=*pg;k++) Ut[k*(*Ntimes)+s]=ME(Utt,s,k-1);
    } /*s=1..Ntimes Beregning af obs teststørrelser */

    for (k=1;k<=*antsim;k++) {
      mat_zeros(Delta); mat_zeros(Delta2); vec_zeros(tmpv1);
      for (i=0;i<*antclust;i++) { /* random=gasdev(&idum); */ 
	random=norm_rand();
	scl_mat_mult(random,W4t[i],tmpM1); mat_add(tmpM1,Delta,Delta);
	scl_mat_mult(random,Uti[i],tmpM2); mat_add(tmpM2,Delta2,Delta2);
      }

      extract_row(Delta,*Ntimes-1,tmpv1); 

      for (s=1;s<*Ntimes;s++) { time=times[s]-times[0];
	scl_vec_mult(time/tau,tmpv1,xi); extract_row(Delta,s,rowX);
	vec_subtr(rowX,xi,difX);

	if (*addresamp==1) {
	  if (k<51) { 
	    for (i=0;i<*px;i++) {l=(k-1)*(*px)+i;
	      addproc[l*(*Ntimes)+s]=ME(Delta,s,i);}}
	}

	for (i=0;i<*px;i++) {
	  VE(difX,i)=fabs(VE(difX,i));
	  l=(*px+i);
	  if (VE(difX,i)>test[l*(*antsim)+k]) test[l*(*antsim)+k]=VE(difX,i);
	  VE(xi,i)=fabs(ME(Delta,s,i))/sqrt(Rvcu[(i+1)*(*Ntimes)+s]);
	  if (VE(xi,i)>test[i*(*antsim)+k]) test[i*(*antsim)+k]=VE(xi,i); }

	if (*wscore>=1) {
	  extract_row(Delta2,s,zi); 
	  if ((s>*wscore) && (s<*Ntimes-*wscore))  {
	    for (i=0;i<*pg;i++) {VE(zi,i)=fabs(ME(Delta2,s,i))/sqrt(VE(varUthat[s],i));
	      if (VE(zi,i)>simUt[i*(*antsim)+k]) simUt[i*(*antsim)+k]=VE(zi,i); }

	    if (k<50) { 
	      for (i=0;i<*pg;i++) { l=(k-1)*(*pg)+i;
		Uit[l*(*Ntimes)+s]=ME(Delta2,s,i)/sqrt(VE(varUthat[s],i));}}
	  } 
	} /* weigted score */
	else {
	  extract_row(Delta2,s,zi); 
	  for (i=0;i<*pg;i++) {
	    if (fabs(VE(zi,i))>simUt[i*(*antsim)+k]) 
	      simUt[i*(*antsim)+k]=fabs(VE(zi,i)); }

	  if (k<50) { 
	    for (i=0;i<*pg;i++) { l=(k-1)*(*pg)+i;
	      Uit[l*(*Ntimes)+s]=ME(Delta2,s,i);}}
	} /* else wscore=0 */ 

      }  /* s=1..Ntims */
    }  /* k=1..antsim */

  } /* sim==1 */ // }}}

  PutRNGstate();  /* to use R random normals */

  // {{{ freeing 
  
  if (*sim==1) free_mats(&Delta,&Delta2,&tmpM2,&tmpM1,NULL); 

  free_mats(&Vcov,&Utt,&VU,&ZX,&COV,&dM1M2,&AI,&A,&ZXAI,&tmp1,&tmp2,&tmp3,NULL);
  free_mats(&ldesignX,&cdesX,&cdesX2,&cdesX3,&ldesignG,&M1,&dS,&S1,&SI,&S2,NULL);
  free_mats(&ZP,&ZPX,&dYI,&Ct,&M1M2t,&RobVbeta,&CtVUCt,NULL); 
  free_vecs(&ta,&ahatt,&Uprofile,&lamtt,&lamt,&one,&xi,&zcol,&Gbeta,NULL);
  free_vecs(&VdA,&dA,&MdA,&xtilde,&zi,&U,&beta,&delta,&zav,&difzzav,&weight,NULL);
  free_vecs(&offset,&tmpv1,&tmpv2,&rowX,&rowZ,&difX,&VdB,&reszpbeta,NULL);
  free_vecs(&vrisk,&res1dim,NULL); 

  if (*robust==1) {
    for (j=0;j<*antclust;j++) {
      free_mat(W3t[j]); free_mat(W4t[j]); free_mat(W2t[j]);free_vec(W2[j]); free_vec(W3[j]);
      free_mat(Uti[j]); }
    for (j=0;j<*antpers;j++) { free_mat(AIxit[j]); }
    for(j=0;j<*Ntimes;j++) free_vec(varUthat[j]); 
  }

  for (j=0;j<*Ntimes;j++) {
    free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]);free_mat(M1M2[j]);free_mat(ZXAIs[j]);
    free_vec(ZXdA[j]);free_mat(St[j]); } 
    free(cluster); free(ipers); free(imin); 
    // }}}
}
