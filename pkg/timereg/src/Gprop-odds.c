//#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
/* ====================================================== */
void Gtranssurv(times,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,
betaS,Nit,cu,vcu,loglike,Iinv,Vbeta,detail,sim,antsim,
rani,Rvcu,RVbeta,test,testOBS,Ut,simUt,Uit,id,status,wscore,
score,dhatMit,dhatMitiid,retur,exppar,sym,mlestart)
double *designX,*designG,*times,*betaS,*start,*stop,*cu,*loglike,*Vbeta,*RVbeta,
*vcu,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*score,*dhatMit,*dhatMitiid;
int *nx,*px,*ng,*pg,*antpers,*Ntimes,*Nit,*detail,*sim,*antsim,*rani,*id,*status,
*wscore,*retur,*exppar,*sym,*mlestart;
{
  matrix *ldesignX,*cdesG,*ldesignG,*cdesX,*cdesX2,*cdesX3,*cdesX4,*CtVUCt,*A,*AI;
  matrix *dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ddesG,*ZP,*ZPX; 
  matrix *tmp1,*tmp2,*tmp5,*tmp3,*dS,*S1,*SI,*S2,*M1,*VU,*ZXAI,*VUI, *tmp6; // Added tmp6 
  matrix *RobVbeta,*Delta,*tmpM1,*Utt,*Delta2,*tmpM2;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes];
  matrix *dW3t[*antpers],*W3t[*antpers],*W4t[*antpers],*W2t[*antpers],*AIxit[*antpers],*Uti[*antpers],*tmp4,*Fst[(*Ntimes)*(*Ntimes)]; 
  matrix *dG[*Ntimes],*cumdG,*Ft[*Ntimes],*ZcX2AIs[*Ntimes],*ZcX2[*Ntimes],*S0tI[*Ntimes],*Ident,*gt[*Ntimes],*q2t[*Ntimes],*G1mG2t[*Ntimes],*q1t[*antpers]; 
  vector *dLamt[*antpers]; 
  vector *dA,*VdA,*dN,*MdA,*delta,*zav,*lamt,*plamt,*dlamt;
  vector *xi,*zi,*U,*beta,*xtilde; 
  vector *Gbeta,*zcol,*one,*difzzav; 
  vector *offset,*weight,*ZXdA[*Ntimes],*varUthat[*Ntimes],*Uprofile;
  vector *Lplamt,*ta,*ahatt,*risk; 
  vector *tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB,*lht; 
  vector *W2[*antpers],*W3[*antpers],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  int t,c,robust=1,pers=0,i,j,k,l,s,it,count,pmax; 
  double time=0,dummy,ll; 
  double tau,hati=0,random,sumscore; 
  int *ipers=calloc(*Ntimes,sizeof(int)); 
  double norm_rand(); 
  void GetRNGstate(),PutRNGstate();  

  for (j=0;j<*antpers;j++) { 
    malloc_vec(*Ntimes,dLamt[j]); malloc_mat(*Ntimes,*px,W3t[j]);
    malloc_mat(*Ntimes,*px,dW3t[j]); malloc_mat(*Ntimes,*px,W4t[j]);
    malloc_mat(*Ntimes,*pg,W2t[j]); malloc_mat(*Ntimes,*pg,Uti[j]);
    malloc_vec(*pg,W2[j]); malloc_vec(*px,W3[j]);
    malloc_mat(*Ntimes,*pg,q1t[j]); malloc_mat(*Ntimes,*px,AIxit[j]);
  }
  malloc_mat(*Ntimes,*px,Delta); malloc_mat(*Ntimes,*px,tmpM1);
  malloc_mat(*Ntimes,*pg,Delta2); malloc_mat(*Ntimes,*pg,tmpM2); malloc_mat(*Ntimes,*pg,Utt);
  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim);

  malloc_vec(*Ntimes,lht);
  malloc_mats(*antpers,*px,&ldesignX,&cdesX,&cdesX2,&cdesX3,&cdesX4,NULL);
  malloc_mats(*antpers,*pg,&ZP,&cdesG,&ldesignG,&ddesG,NULL); 
  malloc_mats(*px,*px,&tmp4,&Ident,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*pg,*pg,&RobVbeta,&tmp1,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&tmp5,&tmp3,&ZX,&dM1M2,&M1M2t,NULL); // changed dims of tmp3 and tmp5 
  malloc_mats(*px,*pg,&cumdG,&ZPX,&dYI,&Ct,NULL); // changed dims of tmp3 and tmp5
  malloc_mat(*px,*pg,tmp6);

  malloc_vecs(*antpers,&Lplamt,&risk,&weight,&dlamt,&plamt,&lamt,&dN,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&ta,&xtilde,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,NULL); 

  identity_matrix(Ident);

  for(j=0;j<*Ntimes;j++) {
    malloc_mat(*px,*px,Ft[j]); malloc_mat(*pg,*px,ZcX2AIs[j]); malloc_mat(*pg,*px,gt[j]);
    malloc_mat(*pg,*px,G1mG2t[j]); malloc_mat(*pg,*px,q2t[j]); malloc_mat(*pg,*px,ZcX2[j]);
    malloc_mat(*px,*px,S0tI[j]); malloc_mat(*px,*pg,dG[j]); malloc_mat(*px,*pg,C[j]);
    malloc_mat(*pg,*px,M1M2[j]); malloc_mat(*pg,*px,ZXAIs[j]); malloc_mat(*px,*pg,dYIt[j]);
    malloc_vec(*px,dAt[j]); malloc_vec(*pg,ZXdA[j]); malloc_mat(*pg,*pg,St[j]);
    malloc_vec(*pg,varUthat[j]);
    for(i=0;i<=j;i++){ malloc_mat(*px,*px,Fst[j*(*Ntimes)+i]); }
  }

//  if (*px>=*pg){ pmax=*px; } else { pmax=*pg; } 
  pmax=max(*px,*pg); ll=0; vec_ones(one);
  for(j=0;j<*pg;j++){ VE(beta,j)=betaS[j]; }
  vec_ones(difX); cu[0]=times[0]; 

  /* Main procedure ================================== */
  for (it=0;it<*Nit;it++){
    vec_zeros(U); mat_zeros(S1);  sumscore=0; 

    for (s=1;s<*Ntimes;s++){
      time=times[s]; 
      mat_zeros(ldesignX); mat_zeros(ldesignG); vec_zeros(risk); 

      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time))  {
	  VE(risk,id[c])=1.0;
	  for(j=0;j<pmax;j++) {
	    if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	    if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; 
	  } 
	  if (time==stop[c] && status[c]==1) {
	    pers=id[c];
	  } 
	  count=count+1; 
	}
      }


//readXZtsimple(antpers,nx,px,designX,pg,designG,
//		start,stop,status,pers, ldesignX, ldesignG,time, s,id); 

      ipers[s]=pers;

      Mv(ldesignG,beta,Gbeta); 
      for (i=0;i<*px;i++){
	VE(tmpv1,i)=cu[(i+1)*(*Ntimes)+s-1];
      }
      Mv(ldesignX,tmpv1,lamt); 

      if (s < 0) {
	print_vec(tmpv1); 
	print_vec(lamt);
      }

      for (j=0;j<*antpers;j++) {
	extract_row(ldesignX,j,xi); 
	extract_row(ldesignG,j,zi); 
	hati=VE(lamt,j); 
	// VE(plamt,j)=VE(risk,j)/(exp(-VE(Gbeta,j))+hati); 
	VE(plamt,j)=1/(exp(-VE(Gbeta,j))+hati); 
	VE(dlamt,j)=-VE(plamt,j)*VE(plamt,j); 
	scl_vec_mult(VE(plamt,j),xi,xtilde); 
	replace_row(cdesX,j,xtilde); 
	scl_vec_mult(VE(dlamt,j),xi,xtilde); 
	replace_row(cdesX2,j,xtilde); 

	scl_vec_mult(-exp(-VE(Gbeta,j)),zi,tmpv2); 
	vM(dG[s-1],xi,rowZ); 
	vec_add(tmpv2,rowZ,rowZ); 
	scl_vec_mult(VE(dlamt,j),rowZ,rowZ); 
	replace_row(ddesG,j,rowZ); 
      }

      if (s < 0) {
	print_vec(plamt); print_vec(dlamt); print_mat(cdesX); 
      }

      MtA(cdesX,ldesignX,A); 
      invert(A,AI); 
      scl_mat_mult(1.0,AI,S0tI[s]); 
      extract_row(ldesignX,pers,xi); 
      Mv(AI,xi,dA); 
      MtA(ldesignG,cdesX,ZX); 
      MxA(ZX,AI,ZXAIs[s]);
      scl_vec_mult(1.0,dA,dAt[s]); 
      if (s<0) { 
	print_mat(AI); 
	print_vec(dA); 
      } 

      /* baseline computed */ 
      for (k=1;k<=*px;k++) {
	cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dA,k-1);
      }
      if (s<0) Rprintf(" %lf \n",cu[1*(*Ntimes)+s]);

      /* First derivative U and Second derivative S  */ 
      extract_row(ldesignG,pers,zi); 
      Mv(ZX, dA, ZXdA[s]); 
      vec_subtr(zi,ZXdA[s],difzzav); 
      vec_add(difzzav,U,U); 

      Mv(ldesignX,dA,lamt);  
      for (j=0;j<*antpers;j++){
	extract_row(ddesG,j,zi); 
	if (s<0 && j<5 ) {
	  Rprintf(" %ld %ld \n",(long int) s, (long int)j); 
	  print_vec(zi); 
	}
	scl_vec_mult(VE(lamt,j),zi,zi); 
	replace_row(ZP,j,zi);
      }

      MtA(ldesignX,ZP,ZPX); 

      MxA(AI,ZPX,tmp6); // Note the use of tmp6 here, instead of tmp3 
      mat_subtr(dG[s-1],tmp6,dG[s]); // Note the use of tmp6 here, instead of tmp3 
      if (s<0) { 
	Rprintf(" %lf \n",ME(A,0,0)); 
	print_mat(ZPX); 
	print_mat(tmp3); 
	print_mat(dG[s]); 
      }

      MxA(ZXAIs[s],ZPX,SI); 
      mat_transp(SI,tmp2); 
      MtA(ldesignG,ZP,tmp1); 
      mat_subtr( tmp1,tmp2, dS); 

      if (s<0) { 
	Rprintf("=================== %lf \n",ME(A,0,0)); 
	print_mat(tmp1); 
	print_mat(tmp2); 
	print_mat(dS); 
      }

      if (*sym==1) {
	mat_transp(dS,tmp1);
	mat_add(tmp1,dS,dS);
	scl_mat_mult(0.5,dS,dS);
      } 
      /* else {m_transp(dS,tmp1); sm_mlt(1,tmp1,dS); }
       */

      mat_add(dS,S1,S1);  
      scl_mat_mult(1.0,S1,St[s]); 

      /* variance and other things */ 
      if (it==((*Nit)-1)) { 
	replace_row(Utt,s,U);


	for (j=0;j<*px;j++)  {
	  for (i=0;i<*antpers;i++){
	    dummy=ME(ldesignX,i,j); 
	    extract_row(cdesX2,i,xi); 
	    scl_vec_mult(dummy,xi,xi); 
	    replace_row(cdesX3,i,xi); 
	  }
	  MtA(ldesignX,cdesX3,A); 
	  MxA(AI,A,tmp4); 
	  Mv(tmp4,dA,xi); 
	  for (k=0;k<*px;k++){
	    ME(Ft[s],j,k)=VE(xi,k); 
	  }

	  VE(lht,s)=VE(lht,s-1)-ME(A,0,0)*(ME(AI,0,0)*ME(AI,0,0));


	  /* Rprintf(" %ld %lf %lf \n",s,lht->ve[s],AI->me[0][0]); */

	  MtA(ldesignG,cdesX3,ZcX2[s]); 
	  /* m_mlt(ZcX2[s],AI,ZcX2AIs[s]);  */

	  MxA(ZX,tmp4,tmp3); 
	  mat_subtr(tmp3,ZcX2[s],tmp5); 

	  Mv(tmp5,dA,zi); 
	  for (k=0;k<*pg;k++) {
	    ME(G1mG2t[s],k,j)=ME(G1mG2t[s],k,j)+VE(zi,k); 
	  }
	}


	/*
	  for (i=0;i<*px;i++){ for (j=0;j<*pg;j++) dM1M2->me[j][i]=dA->ve[i]*difzzav->ve[j];
	  for (i=0;i<*pg;i++) 
	  for (j=0;j<*pg;j++) VU->me[i][j]=VU->me[i][j]+difzzav->ve[i]*difzzav->ve[j]; 

	  m_mlt(AI,ZPX,dYIt[s]); m_sub(Ct,dYIt[s],Ct); C[s]=m_copy(Ct,C[s]); 

	  v_star(dA,dA,VdA); m_add(dM1M2,M1M2t,M1M2t); 
	  M1M2[s]=m_copy(M1M2t,M1M2[s]); 
	  for (k=1;k<=*px;k++)  vcu[k*(*Ntimes)+s]=VdA->ve[k-1]+vcu[k*(*Ntimes)+s-1]; 
	*/

	for (j=0;j<*antpers;j++){
	  extract_row(ldesignX,j,xi);
	  Mv(S0tI[s],xi,rowX);
	  replace_row(AIxit[j],s,rowX);
	  extract_row(ldesignG,j,zi);  
	  Mv(ZX,rowX,rowZ); 
	  vec_subtr(zi,rowZ,zi); 
	  replace_row(q1t[j],s,zi); 
	  VE(dLamt[j],s)=VE(plamt,j)*vec_sum(vec_star(xi,dA,rowX));
	}
      }/* if (it==((*Nit)-1)) */

    } /* Ntimes */ 

    invert(S1,SI); 
    Mv(SI,U,delta); 
    vec_add(beta,delta,beta); 

    if (*detail>=1) { 
      Rprintf("====================Iteration %ld ==================== \n",(long int) it);
      Rprintf("delta \n"); 
      print_vec(delta); 
      Rprintf("Estimate beta \n"); 
      print_vec(beta); 
      Rprintf("Score D l\n"); 
      print_vec(U); 
      Rprintf("Information -D^2 l\n"); 
      print_mat(SI); 
      Rprintf("simple D2 l\n");  
      print_mat(S1); 
    }

    for (k=0;k<*pg;k++) {
      sumscore += VE(U,k); 
    }

    if ((fabs(sumscore)<0.000001) & (it<*Nit-2)){ 
      it=*Nit-2; 
    }
  } /* it */
  for (k=0;k<*pg;k++) { 
    score[k]=VE(U,k); 
  }


  /* computation of q(t) */
  for (s=1;s<*Ntimes;s++) {
    mat_zeros(M1M2t); 

    for (t=s;t<*Ntimes;t++) { 
      identity_matrix(tmp4); 
      identity_matrix(M1); 
      for (k=s;k<t;k++) {
	if (k>s) {
	  scl_mat_mult(1,M1,tmp4); 
	}
	mat_subtr(Ident,Ft[k],A); 
	MxA(tmp4,A,M1); 
      }
      if (s<0) { 
	Rprintf(" %ld %ld %lf \n",(long int) s,(long int) t,ME(M1,0,0));
	matrix *tempTranspose;
	malloc_mat(ncol_matrix(G1mG2t[t]),
		      nrow_matrix(G1mG2t[t]),tempTranspose);
        print_mat(mat_transp(G1mG2t[t],tempTranspose));
	free_mat(tempTranspose);
      }
      MxA(G1mG2t[t],M1,dM1M2);
      mat_add(dM1M2,M1M2t,M1M2t); 
    }
    scl_mat_mult(1,M1M2t,q2t[s]); 
    /* m_mlt(M1M2t,S0tI[s],q2t[s]);  */
    if (s<0){
	matrix *tempTranspose;
	malloc_mat(ncol_matrix(q2t[s]),
		      nrow_matrix(q2t[s]),tempTranspose);
        print_mat(mat_transp(q2t[s],tempTranspose));
	free_mat(tempTranspose);
    }
  }

  /* terms for robust variances ============================ */
  if (robust==1) {
    for (s=1;s<*Ntimes;s++) {
      time=times[s]; 
      cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; Ut[s]=times[s]; 

      /* terms for robust variance   */ 

      for (i=0;i<*antpers;i++) {
	extract_row(AIxit[i],s,xi); 
	Mv(q2t[s],xi,rowZ); 
	extract_row(q1t[i],s,zi); 
	if (s==0) { 
	  print_vec(rowZ); 
	print_vec(zi); 
	}
	vec_add(zi,rowZ,rowZ); 
	if (s==0) {
	  print_vec(rowZ);
	}

	/* mv_mlt(ZXAIs[s],xi,tmpv2);  v_sub(zi,tmpv2,tmpv2); */
	if (i==ipers[s])  { 
	  for (j=0;j<*pg;j++) { 
	    for (k=0;k<*pg;k++) {
	      ME(VU,j,k) += VE(rowZ,j)*VE(rowZ,k);
	    }
	  }
	}

	scl_vec_mult(VE(dLamt[i],s),rowZ,tmpv2);
	vec_subtr(W2[i],tmpv2,W2[i]);
	if (i==ipers[s]) {
	  vec_add(rowZ,W2[i],W2[i]);
	}
	/* if (*ratesim==1) {sv_mlt(hati,tmpv2,rowZ); v_sub(W2[i],rowZ,W2[i]);}  */
	replace_row(W2t[i],s,W2[i]); 
 
	vec_zeros(W3[i]); 
	for (t=1;t<=s;t++) { 
	  if (i==0) {
	    identity_matrix(tmp4); 
	    identity_matrix(M1); 
	    for (k=t;k<=s;k++) {
	      if (k>t) {
		scl_mat_mult(1.0,M1,tmp4); 
	      }
	      if (k>t || t==s) {
		mat_subtr(Ident,Ft[k],A); 
		MxA(tmp4,A,M1);
	      } 
	    }
	    scl_mat_mult(1,M1,Fst[s*(*Ntimes)+t]); 
	  }
	  /*  Fst[s*(*Ntimes)+t]->me[0][0]=exp(-lht->ve[t]+lht->ve[s]); */
	  extract_row(AIxit[i],t,xi); 
	  vM(Fst[s*(*Ntimes)+t],xi,rowX); 
	  scl_vec_mult(VE(dLamt[i],t),rowX,tmpv1);
	  vec_subtr(W3[i],tmpv1,W3[i]); 
	  if (i==ipers[t]){
	    vec_add(rowX,W3[i],W3[i]); 
	  }
	}
	replace_row(W3t[i],s,W3[i]);  

	/* if (hati>0) lle=lle+log(hati); llo=llo+hati; */
	/* if (*ratesim==1) {sv_mlt(hati,rowX,rowX); v_sub(W3[i],rowX,W3[i]);} */

	if (*retur==1){
	  dhatMit[i*(*Ntimes)+s]=1*(i==pers)-hati;
	}
      } /* i=1..antpers */ 

    } /* s=1 ..Ntimes */ 

    MxA(SI,VU,S2); 
    MxA(S2,SI,VU); 

    /* ROBUST VARIANCES   */
    for (s=1;s<*Ntimes;s++) {

      if (s<0){ 
	print_mat(dG[s]); 
      }
      vec_zeros(VdB);
      for (i=0;i<*antpers;i++) {
	Mv(SI,W2[i],tmpv2);   
	Mv(dG[s],tmpv2,rowX);

	extract_row(W3t[i],s,xi); 
	if (s>*Ntimes-5 && i<0){
	  print_vec(xi);
	}

	vec_add(xi,rowX,difX); 
	replace_row(W4t[i],s,difX);
	if (i==-5){
	print_vec(difX);
	}
	vec_star(difX,difX,tmpv1); 
	vec_add(tmpv1,VdB,VdB);

	Mv(St[s],tmpv2,rowZ); 
	extract_row(W2t[i],s,tmpv2); 
	vec_subtr(tmpv2,rowZ,zi); 
	replace_row(Uti[i],s,zi); 

	vec_star(zi,zi,tmpv2); 
	vec_add(tmpv2,varUthat[s],varUthat[s]);

	if (s==1) { 
	  for (j=0;j<*pg;j++){ 
	    for (k=0;k<*pg;k++){
	      ME(RobVbeta,j,k) += VE(W2[i],j)*VE(W2[i],k);
	    }
	  }
	}

	if (*retur==1) {
	  mat_zeros(ldesignX); mat_zeros(ldesignG); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	    if ((start[c]<time) && (stop[c]>=time))  {
	      // VE(risk,id[c])=1.0;
	      for(j=0;j<pmax;j++) {
		if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; 
	      } 
	      if (time==stop[c] && status[c]==1) { pers=id[c]; } 
	      count=count+1; 
	    }
	  }

//readXZtsimple(antpers,nx,px,designX,pg,designG,
//		start,stop,status,pers, ldesignX, ldesignG,time, s,id); 

	  Mv(ldesignG,beta,Gbeta); 

	  for (j=0;j<*antpers;j++) {
	    extract_row(ldesignX,j,xi); extract_row(ldesignG,j,zi);
	    for (i=0;i<*px;i++){ VE(tmpv1,i)=cu[(i+1)*(*Ntimes)+s-1]; }
	    vec_star(tmpv1,xi,rowX); 
	    hati=vec_sum(rowX);
	    VE(plamt,j)=exp(VE(Gbeta,j)+hati)/(1+exp(VE(Gbeta,j)+hati)); 
	    scl_vec_mult(VE(plamt,j),xi,xtilde); 
	    replace_row(cdesX,j,xtilde);       
	  }

	  Mv(cdesX,dAt[s],lamt);  
	  for (j=0;j<*antpers;j++){
	    extract_row(ldesignG,j,zi); scl_vec_mult(VE(lamt,j),zi,zi); 
	    replace_row(ZP,j,zi);
	  } 

	  Mv(ZP,W2[i],reszpbeta); Mv(dYIt[s],W2[i],xi); Mv(cdesX,xi,res1dim); 

	  dhatMitiid[i*(*Ntimes)+s]=dhatMit[i*(*Ntimes)+s]- (VE(reszpbeta,0)- VE(res1dim,0)); 
	} /* retur ==1 */ 

      } /* i =1 ..Antpers */
      for (k=1;k<*px+1;k++) { 
	Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1); vcu[k*(*Ntimes)+s]=VE(VdB,k-1); 
      }
    }  /*  s=1 ..Ntimes */ 
    MxA(RobVbeta,SI,tmp1); 
    MxA(SI,tmp1,RobVbeta);
  }

  for(j=0;j<*pg;j++) { 
    betaS[j]= VE(beta,j); loglike[0]=ll;
    for (k=0;k<*pg;k++){ 
      Iinv[k*(*pg)+j]=ME(SI,j,k); Vbeta[k*(*pg)+j]=-ME(VU,j,k); 
      RVbeta[k*(*pg)+j]=-ME(RobVbeta,j,k); 
    } 
  } 

  if (*sim==1) {
    Rprintf("Simulations start N= %ld \n",(long int) *antsim);
    GetRNGstate();  /* to use R random normals */

    tau=times[*Ntimes-1]-times[0];
    for (i=1;i<=*px;i++){
      VE(rowX,i-1)=cu[i*(*Ntimes)+(*Ntimes-1)];
    }

    /* Beregning af OBS teststørrelser */
    for (s=1;s<*Ntimes;s++) { 
      time=times[s]-times[0]; 

      for (i=1;i<=*px;i++) {
	VE(xi,i-1)=fabs(cu[i*(*Ntimes)+s])/sqrt(Rvcu[i*(*Ntimes)+s]);
	if (VE(xi,i-1)>testOBS[i-1]) testOBS[i-1]=VE(xi,i-1); 
      }

      scl_vec_mult(time/tau,rowX,difX);
      for (i=1;i<=*px;i++) {
	VE(xi,i-1)=cu[i*(*Ntimes)+s];
      }
      vec_subtr(xi,difX,difX);
      for (i=0;i<*px;i++) {
	VE(difX,i)=fabs(VE(difX,i)); 
	l=(*px+i);
	if (VE(difX,i)>testOBS[l]) testOBS[l]=VE(difX,i);
      }

      if (*wscore>=1) {  /* sup beregnes i R */ 
	if ((s>*wscore) && (s<*Ntimes-*wscore)) {
	  extract_row(Utt,s,rowZ);
	  for (i=0;i<*pg;i++) {
	    VE(rowZ,i)=VE(rowZ,i)/sqrt(VE(varUthat[s],i));
	  }
	  replace_row(Utt,s,rowZ); /* scaled score process */ 
	} else {
	  vec_zeros(rowZ); 
	  replace_row(Utt,s,rowZ);
	}
      }
      for (k=1;k<=*pg;k++){
	Ut[k*(*Ntimes)+s]=ME(Utt,s,k-1);
      }
    } /*s=1..Ntimes Beregning af obs teststørrelser */

    for (k=1;k<*antsim;k++) {
      mat_zeros(Delta); 
      mat_zeros(Delta2); 
      vec_zeros(tmpv1);
      for (i=0;i<*antpers;i++) {
	/* random=gasdev(&idum); */ 
	random=norm_rand();
	scl_mat_mult(random,W4t[i],tmpM1); 
	mat_add(tmpM1,Delta,Delta);
	scl_mat_mult(random,Uti[i],tmpM2); 
	mat_add(tmpM2,Delta2,Delta2);
      }

      extract_row(Delta,*Ntimes-1,tmpv1); 

      for (s=1;s<*Ntimes;s++) { 
	time=times[s]-times[0];
	scl_vec_mult(time/tau,tmpv1,xi); 
	extract_row(Delta,s,rowX);
	vec_subtr(rowX,xi,difX);

	for (i=0;i<*px;i++) {
	  VE(difX,i)=fabs(VE(difX,i));
	  l=(*px+i);
	  if (VE(difX,i)>test[l*(*antsim)+k]) test[l*(*antsim)+k]=VE(difX,i);
	  VE(xi,i)=fabs(ME(Delta,s,i))/sqrt(Rvcu[(i+1)*(*Ntimes)+s]);
	  if (VE(xi,i)>test[i*(*antsim)+k]) test[i*(*antsim)+k]=VE(xi,i); 
	}

	if (*wscore>=1) {
	  extract_row(Delta2,s,zi); 
	  if ((s>*wscore) && (s<*Ntimes-*wscore))  {
	    for (i=0;i<*pg;i++) {
	      VE(zi,i)=fabs(ME(Delta2,s,i))/sqrt(VE(varUthat[s],i));
	      if (VE(zi,i)>simUt[i*(*antsim)+k]) simUt[i*(*antsim)+k]=VE(zi,i); 
	    }

	    if (k<50) { 
	      for (i=0;i<*pg;i++) { 
		l=(k-1)*(*pg)+i;
		Uit[l*(*Ntimes)+s]=ME(Delta2,s,i)/sqrt(VE(varUthat[s],i));
	      }
	    }
	  } 
	} else { /* weigted score */
	  extract_row(Delta2,s,zi); 
	  for (i=0;i<*pg;i++) {
	    if (fabs(VE(zi,i))>simUt[i*(*antsim)+k]) simUt[i*(*antsim)+k]=fabs(VE(zi,i)); 
	  }

	  if (k<50) { 
	    for (i=0;i<*pg;i++) { 
	      l=(k-1)*(*pg)+i;
	      Uit[l*(*Ntimes)+s]=ME(Delta2,s,i);
	    }
	  }
	} /* else wscore=0 */ 

      }  /* s=1..Ntims */
    }  /* k=1..antsim */
    PutRNGstate();  /* to use R random normals */
  } /* sim==1 */
  
  free_mats(&cumdG,&tmp4,&Ident,&ddesG,&Utt,&tmpM2,&VUI,&ZX,&COV,
		&dM1M2,&AI,&A,&ZXAI,&tmp1,&tmp2,&tmp3,&ldesignX,&cdesX,
		&cdesX2,&cdesX4,&cdesX3,&cdesG,&ldesignG,&M1,&dS,&S1,&SI,NULL);
  free_mats(&S2,&VU,&ZP,&ZPX,&dYI,&Ct,&M1M2t,&RobVbeta,&Delta,&Delta2,
		&tmpM1,&CtVUCt,NULL); 

  free_vecs(&risk,&ta,&ahatt,&Uprofile,&dlamt,&plamt,&lamt,&one,&xi,&zcol,&Gbeta,&VdA,&dA,&MdA,&xtilde,&zi,&U,&beta,&delta,&zav,&difzzav,&weight,&offset,&tmpv1,&tmpv2,&rowX,&rowZ,&difX,&VdB,&reszpbeta,&res1dim,NULL); 


  free_mat(tmp6);

  for (j=0;j<*antpers;j++) {
    free_mat(W3t[j]); free_mat(dW3t[j]); 
    free_mat(W4t[j]); free_mat(W2t[j]); free_vec(W2[j]); free_vec(W3[j]);
    free_mat(AIxit[j]); free_mat(Uti[j]); free_vec(dLamt[j]); free_mat(q1t[j]); 
  }

  for (j=0;j<*Ntimes;j++) {
    free_mat(gt[j]); free_mat(G1mG2t[j]); free_mat(q2t[j]);
    free_mat(Ft[j]); free_mat(ZcX2AIs[j]); free_mat(ZcX2[j]); free_mat(S0tI[j]); 
    free_mat(dG[j]); free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]); free_mat(M1M2[j]);
    free_mat(ZXAIs[j]); free_vec(ZXdA[j]); free_mat(St[j]); free_vec(varUthat[j]); 
    for(i=0;i<=j;i++) free_mat(Fst[j*(*Ntimes)+i]); 
  } 

  free(ipers); 

}
