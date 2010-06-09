#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
/* ====================================================== */
void twostagereg(times,Ntimes,designX,nx,px,designG,ng,pg,
antpers,start,stop, betaS,Nit,cu,
vcu,Iinv,Vbeta, detail,Rvcu,RVbeta,
id,status, ratesim, score, robust, clusters,
antclust,betafixed, theta, vartheta,thetascore, inverse,
clustsize,desthetaI, ptheta,SthetaI,step,idiclust)
double *designX,*designG,*times,*betaS,*start,*stop,*cu,*Vbeta,*RVbeta,*vcu,*Rvcu,*Iinv,*score,*theta,*vartheta,*thetascore,*desthetaI,*SthetaI,*step;
int *nx,*px,*ng,*pg,*antpers,*Ntimes,*Nit,*detail,*id,*status,*ratesim,*robust,
*clusters,*antclust,*betafixed,*inverse,*clustsize,*ptheta,*idiclust;
{
// {{{ defining variables
  matrix *ldesignX,*ldesG0,*ldesignG,*cdesX,*cdesX2,*cdesX3,*CtVUCt,*A,*AI;
  matrix *Vcov,*dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *tmp1,*tmp2,*tmp3,*dS,*S1,*SI,*S2,*M1,*VU,*ZXAI,*VUI,*RobVbeta;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *W3t[*antclust],*W4t[*antclust],*W2t[*antclust],*AIxit[*antpers],*destheta,*d2UItheta,*d2Utheta,*varthetascore,*Sthetaiid[*antclust],*Stheta; 
  matrix *Ftilde,*Gtilde;
  vector *dA,*VdA,*dN,*MdA,*delta,*zav,*lamt,*lamtt,*xi,*zi,*U,*beta,*xtilde; 
  vector *Gbeta,*zcol,*one,*difzzav; 
  vector *offset,*weight,*ZXdA[*Ntimes],*Uprofile;
  vector *ahatt,*phit[*Ntimes],*dbetaNH[*antclust],*dANH[*antclust]; 
  vector *tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB,*atrisk[*antpers]; 
  vector *W2[*antclust],*W3[*antclust],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  vector *vthetascore,*vtheta1,*vtheta2,*dtheta,*thetaiid[*antclust]; 
  vector *dAiid[*antclust]; 
  int c,pers=0,i,j,k,l,s,it,count,sing,pmax,v,
      *cluster=calloc(*antpers,sizeof(int));
  double *Nt=calloc(*antclust,sizeof(double)),dtime,time,dummy,ll,lle,llo;
  double tau,hati,scale,sumscore=999,*Nti=calloc(*antpers,sizeof(double)), theta0=0;
  double *thetaiidscale=calloc(*antclust,sizeof(double)),
         *NH=calloc(*antclust,sizeof(double)),
	 *HeH=calloc(*antclust,sizeof(double)),
	 *H2eH=calloc(*antclust,sizeof(double)),
	 *Rtheta=calloc(*antclust,sizeof(double)),
         *Hik=calloc(*antpers,sizeof(double)),
	 Dthetanu=1,DDthetanu=1; 
  int *ipers=calloc(*Ntimes,sizeof(int));

  for (j=0;j<*antclust;j++) { Nt[j]=0; NH[j]=0; 
    malloc_vec(*pg,dbetaNH[j]); malloc_vec(*px,dANH[j]);  malloc_vec(*ptheta,dAiid[j]); 
    malloc_vec(*ptheta,thetaiid[j]); malloc_mat(*ptheta,*ptheta,Sthetaiid[j]); 
  }
  for (i=0;i<*antpers;i++) { Nti[i]=0; malloc_vec(*Ntimes,atrisk[i]); }
  malloc_mat(*antpers,*ptheta,destheta); 
  malloc_mat(*ptheta,*ptheta,d2Utheta); malloc_mat(*ptheta,*ptheta,d2UItheta); 
  malloc_mat(*ptheta,*ptheta,Stheta); 
  malloc_mat(*ptheta,*ptheta,varthetascore); 
  malloc_vecs(*ptheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,NULL);

  malloc_mat(*ptheta,*px,Ftilde); malloc_mat(*ptheta,*pg,Gtilde); 
  if (*robust==1) {
    for (j=0;j<*antclust;j++) { malloc_mat(*Ntimes,*px,W3t[j]); 
      malloc_mat(*Ntimes,*px,W4t[j]); malloc_mat(*Ntimes,*pg,W2t[j]); 
      malloc_vec(*pg,W2[j]); malloc_vec(*px,W3[j]); }
    for (j=0;j<*antpers;j++) {malloc_mat(*Ntimes,*px,AIxit[j]);}
  }
  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 

  malloc_mats(*antpers,*px,&ldesignX,&cdesX,&cdesX2,&cdesX3,NULL); 
  malloc_mats(*antpers,*pg,&ZP,&ldesignG,&ldesG0,NULL); 
  malloc_mats(*px,*px,&Vcov,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*pg,*pg,&RobVbeta,&tmp1,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*px,*pg,&tmp3,&ZPX,&dYI,&Ct,NULL); 

  malloc_vecs(*antpers,&weight,&lamtt,&lamt,&dN,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&xtilde,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,NULL); 

  for(j=0;j<*Ntimes;j++) { malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]); 
    malloc_mat(*pg,*px,ZXAIs[j]); malloc_mat(*px,*pg,dYIt[j]); malloc_vec(*px,dAt[j]);
    malloc_vec(*pg,ZXdA[j]); malloc_mat(*pg,*pg,St[j]); malloc_vec(*px,phit[j]);  } 

  if (*px>=*pg) pmax=*px; else pmax=*pg; ll=0; 
  for(j=0;j<*pg;j++) VE(beta,j)=betaS[j]; 
  for(j=0;j<*antpers;j++) {Hik[j]=0; VE(one,j)=1; VE(weight,j)=1; VE(offset,j)=1;} 
  // }}}
  

  for (it=0;it<*Nit;it++) // {{{ cox aalen it
    {
      vec_zeros(U); mat_zeros(S1);  sumscore=0;   
      for (s=1;s<*Ntimes;s++)
	{
	  time=times[s]; vec_zeros(dN); sing=0; mat_zeros(ldesignX); mat_zeros(ldesignG); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	    {
	      if ((start[c]<time) && (stop[c]>=time))  {
		cluster[id[c]]=clusters[c];
		for(j=0;j<pmax;j++) {
		  if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		  if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; } 
		if (s==1 && it==0) for(j=0;j<*ptheta;j++)
		  ME(destheta,id[c],j)=desthetaI[j*(*ng)+c];
		if (time==stop[c] && status[c]==1) {VE(dN,id[c])=1; pers=id[c];} 
		if (it==*Nit-1) VE(atrisk[id[c]],s)=1; 
		count=count+1; }
	    }
	  ipers[s]=pers;

	  Mv(ldesignG,beta,Gbeta); 

	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignX,j,xi); 
	      dummy=exp(VE(Gbeta,j));           /* W D(exp(z*beta))X */ 
	      scl_vec_mult(VE(weight,j)*dummy,xi,xtilde); replace_row(cdesX,j,xtilde); 
	    }
	  if (it==(*Nit-1) && s==1) { ldesG0=mat_copy(ldesignG,ldesG0); 
	    cdesX2=mat_copy(cdesX,cdesX2); }

	  scale=VE(weight,pers); 

	  MtA(cdesX,ldesignX,A); invert(A,AI); 
	  if (ME(AI,0,0)==0) printf(" X'X not invertible at time %lf \n",time);

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

	    for (i=0;i<*px;i++) for (j=0;j<*pg;j++) ME(dM1M2,j,i)=VE(dA,i)*VE(difzzav,j);
	    for (i=0;i<*pg;i++) { 
	      for (j=0;j<*pg;j++) ME(VU,i,j)=ME(VU,i,j)+VE(difzzav,i)*VE(difzzav,j); }

	    MxA(AI,ZPX,dYIt[s]); mat_subtr(Ct,dYIt[s],Ct); C[s]=mat_copy(Ct,C[s]); 

	    vec_star(dA,dA,VdA); mat_add(dM1M2,M1M2t,M1M2t); M1M2[s]=mat_copy(M1M2t,M1M2[s]); 

	    for (k=1;k<=*px;k++) {cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dA,k-1); 
	      vcu[k*(*Ntimes)+s]=VE(VdA,k-1)+vcu[k*(*Ntimes)+s-1];}
	    if (*robust==1) 
	      for (j=0;j<*antpers;j++)
		{extract_row(ldesignX,j,xi); Mv(AI,xi,rowX); replace_row(AIxit[j],s,rowX);}

	    for (i=0;i<*antpers;i++) {
	      extract_row(cdesX,i,rowX); scl_vec_mult(1/VE(weight,i),rowX,rowX); 
	      vec_star(rowX,dAt[s],tmpv1); hati=vec_sum(tmpv1); 
	      /*  printf(" %ld %ld %lf %lf \n",i,s,Hik[i],hati); */
	      Hik[i]=Hik[i]+hati; 
	      extract_row(ldesignG,i,rowZ); 
	    }
	    Nti[pers]=Nti[pers]+1; 

	    for (i=0;i<*antpers;i++) {
	      j=cluster[i]; Nt[j]=Nt[j]+(i==pers); }

	  } /* it==*Nit-1 */ 
	} /* Ntimes */ 


      /* for (k=0;k<*pg;k++) ME(S1,k,k)=ME(S1,k,k)+*ridge;  */
      if (*betafixed==0) {
      invert(S1,SI); 

      Mv(SI,U,delta); if (*betafixed==0) vec_add(beta,delta,beta); 
      MxA(SI,VU,S2); MxA(S2,SI,VU); 

      for (k=0;k<*pg;k++) sumscore= sumscore+VE(U,k); 

      if ((fabs(sumscore)<0.0000000001) & (it<*Nit-1)) it=*Nit-2; 
      }
    } /* it */ // }}}

  for (k=0;k<*pg;k++) score[k]=VE(U,k); 


  for (i=0;i<*antpers;i++) 
  {
      j=cluster[i]; 
      // for (i=0;i<*antpers;i++) if (cluster[i]==j)  { }
      NH[j]=NH[j]+Nti[i]*Hik[i]; 
      extract_row(ldesG0,i,zi); extract_row(cdesX2,i,xi);
      vec_add_mult(dbetaNH[j],zi,Nti[i]*Hik[i],dbetaNH[j]);  
      vec_add_mult(dANH[j],xi,Nti[i]*Hik[i],dANH[j]);  
      if (*betafixed==1) vec_zeros(dbetaNH[j]);  
  }

  lle=0; llo=0;
  for (s=1;s<*Ntimes;s++)  // /* terms for robust variance   */  // {{{ 
    {
      time=times[s]; vec_zeros(dN);dtime=time-times[s-1]; 
      cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; 

      mat_zeros(ldesignX); mat_zeros(ldesignG); 
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	{
	  if ((start[c]<time) && (stop[c]>=time))  {
	    cluster[id[c]]=clusters[c];
	    for(j=0;j<pmax;j++) {
	      if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	      if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; } 
	    if (time==stop[c] && status[c]==1) {pers=id[c];} 
	    count=count+1; }
	}
      Mv(ldesignG,beta,Gbeta); 

      for (j=0;j<*antpers;j++)
	{extract_row(ldesignX,j,xi); extract_row(ldesignG,j,zi);
	  dummy=exp(VE(Gbeta,j))*VE(weight,j); 
	  scl_vec_mult(dummy,xi,xtilde); replace_row(cdesX,j,xtilde); }

      if (*robust==1) {
	for (i=0;i<*antpers;i++) 
	{
	  j=cluster[i]; 

	      extract_row(cdesX,i,rowX); scl_vec_mult(1/VE(weight,i),rowX,rowX); 
	      extract_row(ldesignG,i,zi); extract_row(ldesignX,i,xi); 
	      vec_star(rowX,dAt[s],tmpv1); hati=vec_sum(tmpv1); 

	      Mv(ZXAIs[s],xi,tmpv2);  vec_subtr(zi,tmpv2,tmpv2); 
	      scl_vec_mult(VE(weight,i),tmpv2,tmpv2); 

	      if (i==pers) vec_add(tmpv2,W2[j],W2[j]);
	      if (*ratesim==1) {scl_vec_mult(hati,tmpv2,rowZ); vec_subtr(W2[j],rowZ,W2[j]); }

	      extract_row(AIxit[i],s,rowX);
	      scl_vec_mult(VE(weight,i),rowX,rowX); 

	      if (i==pers) {vec_add(rowX,W3[j],W3[j]); if (hati>0) lle=lle+log(hati);}
	      llo=llo+hati;

	      if (*ratesim==1) {scl_vec_mult(hati,rowX,rowX); vec_subtr(W3[j],rowX,W3[j]);}

	    replace_row(W2t[j],s,W2[j]); 
	    replace_row(W3t[j],s,W3[j]);  

	  } /* j and i=1..antclust */ 
      }

      /* MG baseret varians beregning */
      MxA(C[s],VU,tmp3); MAt(tmp3,C[s],CtVUCt);
      MxA(C[s],SI,tmp3); MxA(tmp3,M1M2[s],COV); 

      for (k=1;k<=*px;k++) {
	vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s]+ME(CtVUCt,k-1,k-1)
	  +2*ME(COV,k-1,k-1); }

      /* */
    } /* s=1 ..Ntimes */  // }}}

  ll=lle-llo; /* likelihood beregnes */

  if (*robust==1) { // {{{ /* ROBUST VARIANCES   */
      for (s=1;s<*Ntimes;s++) {
	vec_zeros(VdB); mat_zeros(Vcov);vec_zeros(zi); 
	for (j=0;j<*antclust;j++) {

	  Mv(SI,W2[j],tmpv2); if (s==1) W2[j]=vec_copy(tmpv2,W2[j]); 
	  Mv(C[s],tmpv2,rowX);
	  if (*betafixed==1) vec_zeros(rowX); 
	  extract_row(W3t[j],s,tmpv1); 
	  vec_add(tmpv1,rowX,difX); 
	  replace_row(W4t[j],s,difX);
	  vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	  Mv(St[s],tmpv2,rowZ); extract_row(W2t[j],s,tmpv2); 

	  if (s==1) { for (c=0;c<*pg;c++) for (k=0;k<*pg;k++)
			ME(RobVbeta,c,k)=ME(RobVbeta,c,k)+VE(W2[j],c)*VE(W2[j],k);}
	} /* j in clusters  */


	for (k=1;k<*px+1;k++) { Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1); }

      }  /*  s=1 ..Ntimes */ 
      /* MxA(RobVbeta,SI,tmp1); MxA(SI,tmp1,RobVbeta); */ 
    } /* if robust==1 */  // }}}

  for(j=0;j<*pg;j++) { betaS[j]= VE(beta,j); 
    for (k=0;k<*pg;k++){ Iinv[k*(*pg)+j]=ME(SI,j,k);
      Vbeta[k*(*pg)+j]=-ME(VU,j,k); 
      RVbeta[k*(*pg)+j]=-ME(RobVbeta,j,k); } } 

if (*detail==1) {
printf(" Cox-Aalen fit \n"); 
printf(" beta : \n"); print_vec(beta); 
printf(" var(beta) : \n"); print_mat(RobVbeta); 
}

  /*===================Estimates theta, two stage approach of glidden ==== */
  for (i=0;i<*ptheta;i++) VE(vtheta1,i)=theta[i]; 

  for (it=0;it<*Nit;it++) // {{{ frailty parameter Newton-Raphson
  {
      for (j=0;j<*antclust;j++) {Rtheta[j]=1; HeH[j]=0;H2eH[j]=0;}

      vec_zeros(vthetascore); mat_zeros(d2Utheta); 
      Mv(destheta,vtheta1,lamtt); 

      for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      for (k=0;k<clustsize[j];k++) {
	   i= idiclust[k*(*antclust)+j]; 
	    theta0=VE(lamtt,i); 
            if (*inverse==1) theta0=exp(theta0); 
	    Rtheta[j]=Rtheta[j]+exp(theta0*Hik[i])-1; 
	    HeH[j]=HeH[j]+Hik[i]*exp(theta0*Hik[i]); 
	    H2eH[j]=H2eH[j]+pow(Hik[i],2)*exp(theta0*Hik[i]); } }
	

      for (j=0;j<*antclust;j++)  if (clustsize[j]>=2) {
         for (k=0;k<clustsize[j];k++) {
	      i= idiclust[k*(*antclust)+j]; 
              extract_row(destheta,i,vtheta2); theta0=VE(lamtt,i); 
              if (*inverse==1){theta0=exp(VE(lamtt,i));Dthetanu=theta0;
		               DDthetanu=pow(theta0,1);}
         }

	 sumscore=0;  ll=0; 
	 if (Nt[j]>=2) 
	 for (k=2;k<=Nt[j];k++) {
	     tau=(Nt[j]-1)/(1+theta0*(Nt[j]-1));
	     lle=-pow((Nt[j]-1),2)/pow((1+theta0*(Nt[j]-1)),2);
	     sumscore=sumscore+tau; ll=ll+lle;}

	 thetaiidscale[j]=sumscore+log(Rtheta[j])/(theta0*theta0)
	    -(1/theta0+Nt[j])*HeH[j]/Rtheta[j]+NH[j]; 

	  scl_vec_mult(thetaiidscale[j]*Dthetanu,vtheta2,thetaiid[j]); 

	  vec_add(vthetascore,thetaiid[j],vthetascore); 

	  for (i=0;i<*ptheta;i++)
	    for (k=0;k<*ptheta;k++)
	      ME(Sthetaiid[j],i,k)=VE(vtheta2,i)*VE(vtheta2,k)*pow(Dthetanu,2)*
		(ll+(2/pow(theta0,2))*HeH[j]/Rtheta[j]
		 -(2/pow(theta0,3))*log(Rtheta[j])-(1/theta0+Nt[j])*
		 (H2eH[j]*Rtheta[j]-HeH[j]*HeH[j])/pow(Rtheta[j],2));

	  mat_add(d2Utheta,Sthetaiid[j],d2Utheta); 
	}
      if (*inverse==-10) { 
	for (i=0;i<*ptheta;i++) for (k=0;k<*ptheta;k++)
	    ME(Stheta,i,k)=VE(vthetascore,i)*DDthetanu/Dthetanu; 
	mat_add(d2Utheta,Stheta,d2Utheta); }
     
      LevenbergMarquardt(d2Utheta,d2UItheta,vthetascore,dtheta,step,step);
// invert(d2Utheta,d2UItheta); 

      if (*detail==1) { 
	printf("====================Iteration %d ==================== \n",it);
	printf("Estimate theta \n"); print_vec(vtheta1); 
	printf("Score D l\n");  print_vec(vthetascore); 
	printf("Information D^2 l\n"); print_mat(d2UItheta); 
      }

//    Mv(d2UItheta,vthetascore,dtheta); scl_vec_mult(*step,delta,delta); 
      vec_subtr(vtheta1,dtheta,vtheta1); 

      sumscore=0; 
      for (k=0;k<*pg;k++) sumscore= sumscore+fabs(VE(vthetascore,k));

      if ((sumscore<0.0000000001) & (it<*Nit-1)) it=*Nit-1; 
    } /* it theta Newton-Raphson */  // }}}


  for (i=0;i<*ptheta;i++) { theta[i]=VE(vtheta1,i);
    thetascore[i]=VE(vthetascore,i); }

  /* terms for robust variances ============================ */
  if (*robust==1) {
    mat_zeros(Gtilde);   
    for (k=0;k<*antclust;k++) if (clustsize[k]>=2) {
      for (j=0;j<clustsize[k];j++) {
	   i= idiclust[j*(*antclust)+k]; 
	theta0=VE(lamtt,i); 
	dummy=(1/(theta0*Rtheta[k]))*exp(theta0*Hik[i])-
	  (1/theta0+Nt[k])*(1+theta0*Hik[i])*exp(theta0*Hik[i])
	  /Rtheta[k]+Nti[i]+
	  (1+theta0*Nt[k])*exp(theta0*Hik[i])*HeH[k]/pow(Rtheta[k],2);
        extract_row(ldesG0,i,zi); extract_row(destheta,i,vtheta1); 
        for (c=0;c<*ptheta;c++) for (l=0;l<*pg;l++) 
	  ME(Gtilde,c,l)=ME(Gtilde,c,l)+ 
	    VE(zi,l)*VE(vtheta1,c)*dummy*Hik[i];  }

      for (s=1;s<*Ntimes;s++) {
	if (k==0) {
	  mat_zeros(Ftilde); 
	  for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {

      for (v=0;v<clustsize[j];v++) {
	   i= idiclust[v*(*antclust)+j]; 
           if (VE(atrisk[i],s)!=0) {
//	    for (i=0;i<*antpers;i++) if (cluster[i]==j && VE(atrisk[i],s)!=0) {
	      theta0=VE(lamtt,i); 
	      dummy=(1/(theta0*Rtheta[j]))*exp(theta0*Hik[i])-
		(1/theta0+Nt[j])*(1+theta0*Hik[i])*exp(theta0*Hik[i])
		/Rtheta[j]+Nti[i]+
		(1+theta0*Nt[j])*exp(theta0*Hik[i])*HeH[j]/pow(Rtheta[j],2);

	      extract_row(destheta,i,vtheta1); 
	      extract_row(cdesX2,i,xi); scl_vec_mult(VE(atrisk[i],s),xi,xi); 
	      for (c=0;c<*ptheta;c++) for (l=0;l<*px;l++) 
		ME(Ftilde,c,l)=ME(Ftilde,c,l)+ 
		  VE(xi,l)*VE(vtheta1,c)*dummy;  
	    }
	  }
	}  
	}
	extract_row(W4t[k],s,tmpv1); extract_row(W4t[k],s-1,xi); 
	vec_subtr(tmpv1,xi,xi); Mv(Ftilde,xi,vtheta2); 
	vec_add(dAiid[k],vtheta2,dAiid[k]); 
      } /* s=1..Ntimes */ 
    } /* k=1..antclust */ 

    /* printf(" g er Gtilde \n"); print_mat(Gtilde);  */

    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      if (*betafixed==1) vec_zeros(W2[j]); 

      Mv(Gtilde,W2[j],vtheta2); 
/* printf(" %ld \n",j);print_vec(W2[j]);print_vec(vtheta2); print_vec(dAiid[j]);*/

      vec_add_mult(thetaiid[j],vtheta2,Dthetanu,thetaiid[j]);
      vec_add_mult(thetaiid[j],dAiid[j],Dthetanu,thetaiid[j]);

      for (i=0;i<*ptheta;i++) for (k=0;k<*ptheta;k++)
	ME(varthetascore,i,k) = ME(varthetascore,i,k) +
	  VE(thetaiid[j],i)*VE(thetaiid[j],k);
    }
  }

  MxA(d2UItheta,varthetascore,d2Utheta); 
  MxA(d2Utheta,d2UItheta,varthetascore);  

  for (j=0;j<*ptheta;j++) for (k=0;k<*ptheta;k++)
    {
      SthetaI[k*(*ptheta)+j]=ME(d2UItheta,j,k); 
      vartheta[k*(*ptheta)+j]=ME(varthetascore,j,k); 
    }

  // {{{ freeing everything
  
  if (*robust==1) {
    for (j=0;j<*antclust;j++) {
      free_mat(W3t[j]); free_mat(W4t[j]); free_mat(W2t[j]);free_vec(W2[j]); free_vec(W3[j]);}
    for (j=0;j<*antpers;j++) { free_mat(AIxit[j]); }
  }
  for (j=0;j<*antpers;j++) { free_vec(atrisk[j]); }
  for (j=0;j<*Ntimes;j++) {
    free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]);free_mat(M1M2[j]);free_mat(ZXAIs[j]);
    free_vec(ZXdA[j]);free_mat(St[j]); free_vec(phit[j]);  } 

  for (j=0;j<*antclust;j++) {
    free_vec(dbetaNH[j]); free_vec(dANH[j]);  free_vec(dAiid[j]); 
    free_vec(thetaiid[j]); free_mat(Sthetaiid[j]); }


  free_mats(&destheta, &d2Utheta, &d2UItheta, &Stheta, &Ftilde, &Gtilde,
	      &varthetascore, &ldesignX,&cdesX,&cdesX2,&cdesX3,
	      &ZP,&ldesignG,&ldesG0,
	      &Vcov,&COV,&A,&AI,&M1,&CtVUCt,
	      &RobVbeta,&tmp1,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,
	      &ZXAI,&ZX,&dM1M2,&M1M2t,&tmp3,&ZPX,&dYI,&Ct,NULL); 

  free_vecs(&weight,&lamtt,&lamt,&dN,&zcol,&Gbeta,&one,&offset,
	      &ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,
	      &xtilde, &tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,
	      &vthetascore,&vtheta1,&dtheta,&vtheta2,
	      &reszpbeta, &res1dim,NULL); 

  free(Nt); free(Nti); free(thetaiidscale); free(NH); free(HeH); free(H2eH); 
  free(Rtheta); free(Hik); 
  free(cluster);  free(ipers); 
  // }}}
}
