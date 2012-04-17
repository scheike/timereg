//#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
void score(times,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,
betaS,Nit,cu,vcu,w,mw,loglike,Iinv,Vbeta,detail,offs,mof,sim,antsim,
rani,Rvcu,RVbeta,
test,testOBS,Ut,simUt,Uit,XligZ,aalen,nb,id,status,wscore,ridge,ratesim,score,dhatMit,gammaiid,dmgiid,
retur,robust,covariance,Vcovs,addresamp,addproc,
resample,gamiid,biid,clusters,antclust,vscore,betafixed,weights,entry,exactderiv,
timegroup,maxtimepoint,stratum)
double
*designX,*designG,*times,*betaS,*start,*stop,*cu,*w,*loglike,*Vbeta,*RVbeta,*vcu,*offs,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*aalen,*ridge,*score,*dhatMit,*gammaiid,*dmgiid,*Vcovs,*addproc,*gamiid,*biid,*vscore,*weights;
int*covariance,*nx,*px,*ng,*pg,*antpers,*Ntimes,*mw,*Nit,*detail,*mof,*sim,*antsim,*rani,*XligZ,*nb,*id,*status,*wscore,*ratesim,*retur,*robust,*addresamp,*resample,*clusters,*antclust,*betafixed,*entry,*exactderiv,*timegroup,*maxtimepoint,*stratum;
{ 
// {{{ setting up memory 
  matrix *X,*Z,*WX,*WZ,*cdesX,*cdesX2,*cdesX3,*CtVUCt,*A,*AI;
  matrix *Vcov,*dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*ZPX; 
  matrix *ZPZ,*tmp2,*tmp3,*dS,*S1,*SI,*S2,*M1,*VU,*ZXAI,*VUI; 
  matrix *RobVbeta,*Delta,*tmpM1,*Utt,*Delta2,*tmpM2;
//  matrix *St[*maxtimepoint],*M1M2[*Ntimes],*C[*maxtimepoint],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes]; 
  matrix *Stg[*maxtimepoint],*Cg[*maxtimepoint]; 
  matrix *W3t[*antclust],*W4t[*antclust],*W2t[*antclust],*AIs[*Ntimes],*Uti[*antclust]; 
  matrix *ZPX1,*ZPZ1,*ZPXo,*ZPZo; 
  vector *dA,*VdA,*MdA,*delta,*zav,*lamt,*lamtt;
  vector *xi,*zi,*U,*beta,*xtilde,*Gbeta,*zcol,*one,*difzzav;
  vector *offset,*weight,*ZXdA[*Ntimes],*varUthat[*maxtimepoint],*Uprofile;
  vector *ta,*ahatt,*vrisk,*tmpv1,*tmpv2,*rowX,*rowZ,*difX,*VdB; 
  vector *W2[*antclust],*W3[*antclust],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  int cin=0,ci=0,c,pers=0,i=0,j,k,l,s,it,count,pmax,
      *imin=calloc(1,sizeof(int)),
      *cluster=calloc(*antpers,sizeof(int)),
      *ipers=calloc(*Ntimes,sizeof(int)); 
  double S0,RR=1,time=0,ll,lle,llo;
  double tau,hati,random,scale,sumscore;
  double *cug=calloc((*maxtimepoint)*(*px+1),sizeof(double)),
         *timesg=calloc((*maxtimepoint),sizeof(double));
  double norm_rand();
  void GetRNGstate(),PutRNGstate();

  /* float gasdev(),expdev(),ran1(); */ 

  GetRNGstate();  /* to use R random normals */

  if (*robust==1) {
    for (j=0;j<*antclust;j++) { malloc_mat(*maxtimepoint,*px,W3t[j]); 
      malloc_mat(*maxtimepoint,*px,W4t[j]); malloc_mat(*maxtimepoint,*pg,W2t[j]); 
      malloc_mat(*maxtimepoint,*pg,Uti[j]);  malloc_vec(*px,W3[j]); 
    }
    for (j=0;j<*Ntimes;j++) {malloc_mat(*px,*px,AIs[j]); }
    for(j=0;j<*maxtimepoint;j++) malloc_vec(*pg,varUthat[j]);
  }
  for (j=0;j<*antclust;j++) malloc_vec(*pg,W2[j]);
  for (c=0;c<*nx;c++) cluster[id[c]]=clusters[c]; 

  if (*sim==1) {
    malloc_mat(*maxtimepoint,*px,Delta);  malloc_mat(*maxtimepoint,*px,tmpM1); 
    malloc_mat(*maxtimepoint,*pg,Delta2); malloc_mat(*maxtimepoint,*pg,tmpM2);  
  }

  malloc_mat(*maxtimepoint,*pg,Utt); 
  malloc_mats(*antpers,*px,&WX,&X,&cdesX,&cdesX2,&cdesX3,NULL); 
//  malloc_mats(*px,*antpers,NULL); 
  malloc_mats(*antpers,*pg,&WZ,&ZP,&Z,NULL); 
  malloc_mats(*px,*px,&Vcov,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mats(*pg,*pg,&RobVbeta,&ZPZ,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,NULL); 
  malloc_mats(*pg,*px,&ZXAI,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*px,*pg,&tmp3,&ZPX,&dYI,&Ct,NULL); 
  malloc_mats(*px,*pg,&ZPX1,NULL); malloc_mats(*pg,*pg,&ZPZ1,NULL); 
  malloc_mats(*px,*pg,&ZPXo,NULL); malloc_mats(*pg,*pg,&ZPZo,NULL); 

  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 
  malloc_vecs(*antpers,&weight,&lamtt,&lamt,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*px,&ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*px,&xtilde,NULL); 
  malloc_vecs(*pg,&tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,NULL); 
  malloc_vec(*nb,ta); 
  malloc_vec(*antpers,vrisk); 

  for(j=0;j<*maxtimepoint;j++) { malloc_mat(*px,*pg,Cg[j]); malloc_mat(*pg,*pg,Stg[j]);}
  for(j=0;j<*Ntimes;j++) { 
    malloc_mat(*px,*pg,C[j]); malloc_mat(*pg,*px,M1M2[j]); 
    malloc_mat(*pg,*px,ZXAIs[j]); malloc_mat(*px,*pg,dYIt[j]); 
    malloc_vec(*px,dAt[j]); malloc_vec(*pg,ZXdA[j]); malloc_mat(*pg,*pg,St[j]); 
  } 

  pmax=max(*px,*pg); ll=0; 
  for(j=0;j<*pg;j++) VE(beta,j)=betaS[j]; 
  for(j=0;j<*antpers;j++) {VE(one,j)=1; VE(weight,j)=1; VE(offset,j)=1;} 
  // }}}
  
  R_CheckUserInterrupt();

//  for (j=0;j<*antclust;j++) { Rprintf("0 w2w2 %d \n",j); print_vec(W2[j]); }

  cu[0]=times[0]; 
  for (it=0;it<*Nit;it++) // {{{ iterations start for cox-aalen model
  {
     if (it>0) {
      vec_zeros(U); mat_zeros(S1);   
      mat_zeros(A); mat_zeros(ZPZ); mat_zeros(ZPX); mat_zeros(ZX); 
      mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); 
      vec_zeros(Gbeta); 
      }
     sumscore=0; S0=0; ci=0; 
     R_CheckUserInterrupt();

      for (s=1;s<*Ntimes;s++)  // {{{ going through time
      {
         time=times[s]; vec_zeros(lamt);

    // Rprintf(" %d it %d s \n",it,s); 
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
		VE(Gbeta,id[c])=vec_prod(zi,beta); 
		RR=exp(VE(Gbeta,id[c]));
		S0+=RR*weights[c]; 
                for(j=0;j<pmax;j++) {
	        if (j<*px) {ME(WX,id[c],j)=weights[c]*RR*designX[j*(*nx)+c];}
	        if (j<*pg) {ME(WZ,id[c],j)=weights[c]*designG[j*(*ng)+c];} 
		}
		if (time==stop[c] && status[c]==1) {pers=id[c];} 
		if (*mof==1) VE(offset,id[c])=offs[c];  
		if (*mw==1) VE(weight,id[c])=weights[c]; 

	    for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=VE(xi,k)*VE(xi,j)*RR*weights[c]; 
              if ((j<*pg) & (k<*px)) ME(ZX,j,k)+=VE(zi,j)*VE(xi,k)*RR*weights[c]; 
	      if ((*exactderiv<2) || (*px==1)) {
                 if ((j<*pg) & (k<*pg)) ME(ZPZ,j,k)+= VE(zi,j)*VE(zi,k)*weights[c]*RR; 
                 if ((j<*pg) & (k<*px)) ME(ZPX,k,j)+= VE(zi,j)*VE(xi,k)*weights[c]*RR;
	      }
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
            for(j=0;j<*px;j++) VE(xi,j)=designX[j*(*nx)+ci]; 
            for(j=0;j<*pg;j++) VE(zi,j)=designG[j*(*nx)+ci]; 
		VE(Gbeta,id[ci])=vec_prod(zi,beta); 
		RR=exp(VE(Gbeta,id[ci]));
	    if (entry[ci]==1)  {
	         replace_row(X,id[ci],xi); replace_row(Z,id[ci],zi); 
		 scl_vec_mult(RR*weights[ci],xi,tmpv1);replace_row(WX,id[ci],tmpv1);
		 scl_vec_mult(weights[ci],zi,tmpv2);replace_row(WZ,id[ci],tmpv2); 
		 VE(weight,id[ci])=weights[ci]; 
		 if (*mof==1) VE(offset,id[ci])=offs[ci];  
	    } 
	    else { 
		    replace_row(X,id[ci],rowX);replace_row(Z,id[ci],rowZ);
		    replace_row(WX,id[ci],rowX); replace_row(WZ,id[ci],rowZ);
		    VE(Gbeta,id[ci])=0; VE(weight,id[ci])=0; 
		    if (*mof==1) VE(offset,id[ci])=offs[ci];  
	    }
	    S0+=entry[ci]*RR*weights[ci]; 
	  for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=entry[ci]*VE(xi,k)*VE(xi,j)*RR*weights[ci]; 
              if ((j<*pg) & (k<*px)) ME(ZX,j,k)+=entry[ci]*VE(zi,j)*VE(xi,k)*RR*weights[ci]; 
	      if ((*exactderiv<2) || (*px==1)) {
                 if ((j<*pg) & (k<*pg)) ME(ZPZ,j,k)+=entry[ci]*VE(zi,j)*VE(zi,k)*weights[ci]*RR;
                 if ((j<*pg) & (k<*px)) ME(ZPX,k,j)+=entry[ci]*VE(zi,j)*VE(xi,k)*weights[ci]*RR;
	      }
	  }
	  ci=ci-1; 
	  pers=id[ci]; 
    }
    // }}}
   ipers[s]=pers;

// Rprintf("___________ %d %d  %lf %d \n",it,s,time,pers); print_mat(Z); print_mat(X); 

   scl_mat_mult(1/S0,ZPZ,ZPZo); scl_mat_mult(1/S0,ZPX,ZPXo);
   // }}}
   
  if (s<0) { Rprintf("======================================================= %d \n",s);
	     print_mat(A); print_mat(ZPX); print_mat(ZX); 
	     print_mat(A);  print_mat(ZPZ); 
           }

   invertS(A,AI,1); 
   if (ME(AI,0,0)==0 && *stratum==0) {Rprintf(" X'X not invertible at time %d %lf \n",s,time); print_mat(A);}
   if (*stratum==1)  {
    for (k=0;k<*px;k++) 
    if (fabs(ME(A,k,k))<0.000001)  ME(AI,k,k)=0; else ME(AI,k,k)=1/ME(A,k,k);
   }

    scale=VE(weight,pers); 
    extract_row(X,pers,xi); scl_vec_mult(scale,xi,xi); 
    Mv(AI,xi,dA); scl_vec_mult(1,dA,dAt[s]); 
    MxA(ZX,AI,ZXAIs[s]); Mv(ZX, dA, ZXdA[s]);  

  if (s<0) {print_mat(A); print_mat(AI); print_mat(ZX); print_mat(ZXAIs[s]); }

  /* First derivative U and Second derivative S  */ 
 
  extract_row(Z,pers,zi);  scl_vec_mult(scale,zi,zi); 
  Mv(ZX, dA, zav); vec_subtr(zi,zav,difzzav); vec_add(difzzav,U,U); 

  if (s<0) { Rprintf(" %d \n",pers); 
	     print_vec(xi); print_vec(zi); print_vec(zav); print_vec(difzzav);
  }


  if (*betafixed==0) 
  if ( (((*exactderiv==1) && (it==(*Nit-1)) && (*px>1))) || ((*exactderiv==2) && (*px>1)) ) {
  mat_zeros(ZPZ1); mat_zeros(ZPX1); 
  for (i=0;i<*antpers;i++)
  {
    extract_row(WX,i,xi); VE(lamt,i)=vec_prod(xi,dA); 
    extract_row(Z,i,zi); scl_vec_mult(VE(lamt,i),zi,rowZ); replace_row(ZP,i,rowZ);
    extract_row(X,i,xi); 
     for(j=0;j<pmax;j++)  for(k=0;k<*pg;k++)   {
         if ((j<*px)) ME(ZPX1,j,k)+=VE(weight,i)*VE(lamt,i)*VE(xi,j)*VE(zi,k); 
         if ((j<*pg)) ME(ZPZ1,j,k)+=VE(weight,i)*VE(lamt,i)*VE(zi,k)*VE(zi,j); 
     }
  } 
  scl_mat_mult(1,ZPZ1,ZPZo); scl_mat_mult(1,ZPX1,ZPXo);
  }

  if (*betafixed==0)  {
   MxA(ZXAIs[s],ZPXo,tmp2); mat_subtr(ZPZo,tmp2, dS); 
   mat_add(dS,S1,S1);  
   scl_mat_mult(1,S1,St[s]); scl_mat_mult(1,S1,Stg[timegroup[s]]);
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

    MxA(AI,ZPXo,dYIt[s]); mat_subtr(Ct,dYIt[s],Ct); 
    scl_mat_mult(1,Ct,C[s]); scl_mat_mult(1,Ct,Cg[timegroup[s]]); 

    vec_star(dA,dA,VdA); mat_add(dM1M2,M1M2t,M1M2t); M1M2[s]=mat_copy(M1M2t,M1M2[s]); 

    for (k=1;k<=*px;k++) {
      cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dA,k-1); 
      cug[k*(*maxtimepoint)+timegroup[s]]=cu[k*(*Ntimes)+s]; 
      vcu[k*(*Ntimes)+s]=VE(VdA,k-1)+vcu[k*(*Ntimes)+s-1];
    }
    if (*robust==1) AIs[s]=mat_copy(AI,AIs[s]); 
  } // }}}

  } // }}} /* Ntimes */ 

/* for (k=0;k<*pg;k++) ME(S1,k,k)=ME(S1,k,k)+*ridge;  */
  if (*betafixed==0)  {
      invert(S1,SI); Mv(SI,U,delta); MxA(SI,VU,S2); MxA(S2,SI,VU); 
  }

      if (*detail==1) { 
        Rprintf("=============Iteration %d =============== \n",it);
	Rprintf("Estimate beta \n"); print_vec(beta); 
	Rprintf("delta beta \n"); print_vec(delta); 
	Rprintf("Score D l\n"); print_vec(U); 
	Rprintf("Information -D^2 l\n"); print_mat(SI); 
      };

      if (*betafixed==0) vec_add(beta,delta,beta); 

      for (k=0;k<*pg;k++) sumscore=sumscore+fabs(VE(U,k)); 
      if ((sumscore<0.0000001) & (it<(*Nit)-2)) {  it=*Nit-2; }
 } /* it */ // }}}

  if (*detail>=2) Rprintf("Fitting done \n"); 

  R_CheckUserInterrupt();

  mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); vec_zeros(Gbeta); 
  lle=0; llo=0; ci=0; 
  for (k=0;k<*pg;k++) score[k]=VE(U,k); 

  mat_zeros(A); mat_zeros(ZPZ); mat_zeros(ZPX); mat_zeros(ZX); 
  mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); 
  vec_zeros(Gbeta); 

  for (s=1;s<*Ntimes;s++) { // {{{ terms for robust variances 
    time=times[s]; 
    cu[s]=times[s]; vcu[s]=times[s]; Rvcu[timegroup[s]]=times[s]; 
    cug[timegroup[s]]=times[s]; timesg[timegroup[s]]=times[s]; 
    Ut[timegroup[s]]=times[s]; 

    if (*robust==1) {
    // {{{ reading design and computing matrix products
	  if (s==1) { // {{{
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	  {
	      if ((start[c]<time) && (stop[c]>=time)) {
                for(j=0;j<*pg;j++) { 
			ME(Z,id[c],j)=designG[j*(*ng)+c]; VE(zi,j)=designG[j*(*ng)+c]; 
		}
		VE(Gbeta,id[c])=vec_prod(zi,beta); 
		RR=exp(VE(Gbeta,id[c]));
                for(j=0;j<pmax;j++) {
	           if (j<*px) {ME(X,id[c],j)=designX[j*(*nx)+c]; }
	           if (j<*px) {ME(WX,id[c],j) =RR*designX[j*(*nx)+c];}
	           if (j<*pg) {ME(WZ,id[c],j)=weights[c]*designG[j*(*ng)+c];} 
//	        if (j<*pg) {ME(Z,id[c],j)=designG[j*(*ng)+c];} 
		}
		if (time==stop[c] && status[c]==1) {pers=id[c];} 
		if (*mof==1) VE(offset,id[c])=offs[c];  
		if (*mw==1) VE(weight,id[c])=weights[c]; 
		} 
		count=count+1; 
	   }
           ci=*nx-1; 
           while ((stop[ci]<time)  & (ci>=0) )  ci=ci-1; 
	  } // }}}

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
		VE(Gbeta,id[c])=vec_prod(zi,beta); 
		RR=exp(VE(Gbeta,id[c]));
		S0+=RR*weights[c]; 
                for(j=0;j<pmax;j++) {
	        if (j<*px) {ME(WX,id[c],j) =weights[c]*RR*designX[j*(*nx)+c];}
	        if (j<*pg) {ME(WZ,id[c],j)=weights[c]*designG[j*(*ng)+c];} 
		}
		if (time==stop[c] && status[c]==1) {pers=id[c];} 
		if (*mof==1) VE(offset,id[c])=offs[c];  
		if (*mw==1) VE(weight,id[c])=weights[c]; 

	    for(j=0;j<pmax;j++) for(k=0;k<pmax;k++)  {
              if ((j<*px) & (k<*px)) ME(A,j,k)+=VE(xi,k)*VE(xi,j)*RR*weights[c]; 
              if ((j<*pg) & (k<*px)) ME(ZX,j,k)+=VE(zi,j)*VE(xi,k)*RR*weights[c]; 
	      if ((*exactderiv<2) || (*px==1)) {
                 if ((j<*pg) & (k<*pg)) ME(ZPZ,j,k)+= VE(zi,j)*VE(zi,k)*weights[c]*RR; 
                 if ((j<*pg) & (k<*px)) ME(ZPX,k,j)+= VE(zi,j)*VE(xi,k)*weights[c]*RR;
	      }
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
            for(j=0;j<*px;j++) VE(xi,j)=designX[j*(*nx)+ci]; 
            for(j=0;j<*pg;j++) VE(zi,j)=designG[j*(*nx)+ci]; 
		VE(Gbeta,id[ci])=vec_prod(zi,beta); 
		RR=exp(VE(Gbeta,id[ci]));
	    if (entry[ci]==1)  {
	         replace_row(X,id[ci],xi); replace_row(Z,id[ci],zi); 
	         scl_vec_mult(RR,xi,tmpv1);replace_row(WX,id[ci],tmpv1);
		 scl_vec_mult(weights[ci],zi,tmpv2);replace_row(WZ,id[ci],tmpv2); 
	    } 
	    else { replace_row(X,id[ci],rowX);replace_row(Z,id[ci],rowZ);
		   replace_row(WX,id[ci],rowX); replace_row(WZ,id[ci],rowZ);
		   VE(Gbeta,id[ci])=0; 
	    }
	  ci=ci-1; 
	  pers=id[ci]; 
    }
    // }}}
	  ipers[s]=pers;
   // }}}
    }
   
//  print_mat(WX); print_mat(X); print_mat(Z); 
// Rprintf("=============== %d  %lf %d \n",s,time,pers); print_mat(Z); print_mat(X); 

    extract_row(WX,pers,xi); hati=vec_prod(xi,dAt[s]); 
    lle=lle+log(hati);

    /* terms for robust variance   */ 
    if (*robust==1) {
    for (i=0;i<*antpers;i++)   // {{{
    {
      cin=cluster[i]; 
      extract_row(WX,i,rowX); extract_row(Z,i,zi); extract_row(X,i,xi); 
      hati=vec_prod(rowX,dAt[s]); 

//      Rprintf("%d %d %d  %d %lf \n",s,i,ipers[s],pers,hati);  
      Mv(ZXAIs[s],xi,tmpv2);  vec_subtr(zi,tmpv2,tmpv2); 
      scl_vec_mult(VE(weight,i),tmpv2,tmpv2); 

      if (i==pers) vec_add(tmpv2,W2[cin],W2[cin]);
      if (*ratesim==1) {scl_vec_mult(hati,tmpv2,rowZ); vec_subtr(W2[cin],rowZ,W2[cin]); }

      if (*robust==1) 
      {
	 Mv(AIs[s],xi,rowX); scl_vec_mult(VE(weight,i),rowX,rowX); 
	 if (i==pers) {vec_add(rowX,W3[cin],W3[cin]); }
	 llo=llo+hati;
	 if (*ratesim==1) {scl_vec_mult(hati,rowX,rowX); vec_subtr(W3[cin],rowX,W3[cin]);}
      }

      if (*retur==1) dhatMit[i*(*Ntimes)+s]=1*(i==pers)-hati;
      if (*retur==2) dhatMit[i]=dhatMit[i]+1*(i==pers)-hati;

    } /* i 1.. antpers */ // }}}
    }

    if (*robust==1) 
    for (j=0;j<*antclust;j++) 
    {
      replace_row(W2t[j],timegroup[s],W2[j]); replace_row(W3t[j],timegroup[s],W3[j]);  
    } 

    /* MG baseret varians beregning */
    MxA(C[s],VU,tmp3); MAt(tmp3,C[s],CtVUCt);
    MxA(C[s],SI,tmp3); MxA(tmp3,M1M2[s],COV); 

    for (k=1;k<=*px;k++) {
      if (*betafixed==0) 
	vcu[k*(*Ntimes)+s]+=ME(CtVUCt,k-1,k-1) +2*ME(COV,k-1,k-1); 
//      else vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s];
    }
    for (k=1;k<=*pg;k++) Ut[k*(*maxtimepoint)+timegroup[s]]=ME(Utt,timegroup[s],k-1);
  }   // }}}

  if (*detail>=2) Rprintf("Robust variances \n"); 

  R_CheckUserInterrupt();

  ll=lle-llo; /* likelihood beregnes */
  if (*detail==1) Rprintf("loglike is  %lf \n",ll);  

  if ((*robust==1)) // {{{ robust variances 
  {
      for (s=1;s<*maxtimepoint;s++) {
	vec_zeros(VdB); mat_zeros(Vcov);

	for (j=0;j<*antclust;j++) {
	  Mv(SI,W2[j],tmpv2); 
	  Mv(Cg[s],tmpv2,rowX);
	  extract_row(W3t[j],s,tmpv1); vec_add(tmpv1,rowX,difX); 
	  if (*betafixed==1) scl_vec_mult(1,tmpv1,difX); 
	  replace_row(W4t[j],s,difX); 
	  vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	  if (*resample==1) {
	    if (s==1) 
	      if (*betafixed==0) {
		for (c=0;c<*pg;c++) gamiid[c*(*antclust)+j]=gamiid[c*(*antclust)+j]+VE(tmpv2,c); 
	      }
	    for (c=0;c<*px;c++) {l=j*(*px)+c; 
	      biid[l*(*maxtimepoint)+s]=biid[l*(*maxtimepoint)+s]+VE(difX,c);
	    } 
	  }

	  if (*covariance==1) {
	    for (k=0;k<*px;k++) for (c=0;c<*px;c++) 
	      ME(Vcov,k,c)=ME(Vcov,k,c)+VE(difX,k)*VE(difX,c);
	  }

	  Mv(Stg[s],tmpv2,rowZ); extract_row(W2t[j],s,tmpv2); 
	  if (*betafixed==0) {
	    vec_subtr(tmpv2,rowZ,zi); replace_row(Uti[j],s,zi); 
	  } else replace_row(Uti[j],s,tmpv2);

	  vec_star(zi,zi,tmpv2); vec_add(tmpv2,varUthat[s],varUthat[s]);

	} /* j in clusters  */

	if (*betafixed==0) 
	  for (i=0;i<*pg;i++) vscore[(i+1)*(*maxtimepoint)+s]=VE(varUthat[s],i);  

	for (k=1;k<*px+1;k++) { Rvcu[k*(*maxtimepoint)+s]=VE(VdB,k-1);
	  if (*covariance==1) {
	    for (j=0;j<*px;j++)  { l=(k-1)*(*px)+j; 
		    Vcovs[l*(*maxtimepoint)+s]=ME(Vcov,k-1,j); 
	    } 
	  } 
	}
      }  /*  s=1 ..Ntimes */ 

    } /* if robust==1 */  // }}}

  if (*detail>=2) Rprintf("Robust variances \n"); 

    // always compute robust se and return iid of gamma iid
    if (*betafixed==0) 
    for (j=0;j<*antclust;j++) {
      Mv(SI,W2[j],tmpv2); 
      for (c=0;c<*pg;c++) for (k=0;k<*pg;k++)
		ME(RobVbeta,c,k)=ME(RobVbeta,c,k)+VE(W2[j],c)*VE(W2[j],k);
      for (k=0;k<*pg;k++) gammaiid[j*(*pg)+k]=VE(tmpv2,k); 
    }
     MxA(RobVbeta,SI,ZPZ); MxA(SI,ZPZ,RobVbeta);

  R_CheckUserInterrupt();

  for(j=0;j<*pg;j++) { betaS[j]= VE(beta,j); 
    loglike[0]=lle; loglike[1]=ll;
    for (k=0;k<*pg;k++){ Iinv[k*(*pg)+j]=ME(SI,j,k);
      Vbeta[k*(*pg)+j]=-ME(VU,j,k); 
      RVbeta[k*(*pg)+j]=-ME(RobVbeta,j,k); 
    } 
  } 

  R_CheckUserInterrupt();

  if (*sim==1) { // {{{ score process simulations
    // Rprintf("Simulations start N= %ld \n",(long int) *antsim);

    tau=times[*Ntimes-1]-times[0];
    for (i=1;i<=*px;i++) VE(rowX,i-1)=cug[i*(*maxtimepoint)+(*maxtimepoint-1)];

    /* Beregning af OBS teststørrelser */
    for (s=1;s<*maxtimepoint;s++) { 
      time=timesg[s]-times[0];  //  FIX 

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

      if (*wscore>=1) {  /* sup beregnes i R */ 
	if ((s>*wscore) && (s<*maxtimepoint-*wscore)) {extract_row(Utt,s,rowZ);
	  for (i=0;i<*pg;i++) VE(rowZ,i) = VE(rowZ,i)/sqrt(VE(varUthat[s],i));
	  replace_row(Utt,s,rowZ); /* scaled score process */ 
	}
	else {vec_zeros(rowZ); replace_row(Utt,s,rowZ);}
      }
      for (k=1;k<=*pg;k++) Ut[k*(*maxtimepoint)+s]=ME(Utt,s,k-1);
    } /*s=1..maxtimepoint Beregning af obs teststørrelser */

    for (k=1;k<=*antsim;k++) {
    R_CheckUserInterrupt();
      mat_zeros(Delta); mat_zeros(Delta2); vec_zeros(tmpv1);
      for (i=0;i<*antclust;i++) { /* random=gasdev(&idum); */ 
	random=norm_rand();
	scl_mat_mult(random,W4t[i],tmpM1); mat_add(tmpM1,Delta,Delta);
	scl_mat_mult(random,Uti[i],tmpM2); mat_add(tmpM2,Delta2,Delta2);
      }

      extract_row(Delta,*maxtimepoint-1,tmpv1); 

      for (s=1;s<*maxtimepoint;s++) { time=timesg[s]-times[0];
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
    }  /* k=1..antsim */

  } /* sim==1 */ // }}}

  PutRNGstate();  /* to use R random normals */

  // {{{ freeing 
  if (*sim==1) free_mats(&Delta,&Delta2,&tmpM2,&tmpM1,NULL); 

  free_mats(&Utt,&WX,&X,&cdesX,&cdesX2,&cdesX3, &WZ,&ZP,&Z,
     &Vcov,&COV,&A,&AI,&M1,&CtVUCt, &RobVbeta,&ZPZ,&tmp2,&dS,&S1,&S2,&SI,&VU,&VUI,
     &ZXAI,&ZX,&dM1M2,&M1M2t, &tmp3,&ZPX,&dYI,&Ct, &ZPX1,&ZPZ1, &ZPXo,&ZPZo,NULL); 

  free_vecs(&reszpbeta,&res1dim,&weight,&lamtt,&lamt,&zcol,&Gbeta,&one,&offset,
            &ahatt,&tmpv1,&difX,&VdB,&rowX,&xi,&dA,&VdA,&MdA,
            &xtilde, &tmpv2,&rowZ,&zi,&U,&beta,&delta,&zav,&difzzav,&Uprofile,
            &ta,&vrisk,NULL); 

  if (*robust==1) {
    for (j=0;j<*antclust;j++) {
      free_mat(W3t[j]); free_mat(W4t[j]); free_mat(W2t[j]); 
      free_vec(W3[j]); free_mat(Uti[j]); 
    }
    for (j=0;j<*Ntimes;j++)   free_mat(AIs[j]); 
    for (j=0;j<*maxtimepoint;j++)  free_vec(varUthat[j]);
  }
  for (j=0;j<*antclust;j++) free_vec(W2[j]);

  for (j=0;j<*Ntimes;j++) {
    free_mat(dYIt[j]); free_vec(dAt[j]); free_mat(C[j]);free_mat(M1M2[j]);
    free_mat(ZXAIs[j]); free_vec(ZXdA[j]);free_mat(St[j]); 
  } 
  for(j=0;j<*maxtimepoint;j++) { free_mat(Cg[j]); free_mat(Stg[j]);}
  free(cluster); free(ipers); free(imin); free(cug); free(timesg); 
  // }}}
}
