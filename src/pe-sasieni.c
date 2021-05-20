//#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include"R_ext/Random.h"

void pes(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,cu,vcu,gamma,Vgamma,status,Ut,intZHZ,intZHdN,mof,offset,mw,weight,Nit,detail,rani,nsim,test)
double *designX,*alltimes,*start,*stop,*cu,*vcu,*designG,*gamma,*Vgamma,*Ut,*intZHZ,*intZHdN,*offset,*weight,*test;
int *detail,*nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*mof,*mw,*Nit,*rani,*nsim;
{
  matrix *S2,*S2I,*Vcov,*X,*WX,*A,*AI,*AIXW,*Z,*WZ;
  matrix *dCGam,*CGam,*Ct,*ICGam,*VarKorG,*dC,*XWZ,*ZWZ,*XWZAI; 
  matrix *Acorb[*Nalltimes],*Vargam,*dVargam,*M1M2[*Ntimes],*GCdM1M2; 
  matrix *C[*Nalltimes],*dM1M2,*M1M2t,*tmpM2,*tmpM3,*tmpM4;
  matrix *Delta,*tmpM1,*dUt[*Ntimes],*Uti[*antpers],*Utiid[*antpers]; 
  vector *VdB,*difX,*xi,*tmpv1,*tmpv2,*gamoff; 
  vector *dA,*rowX,*dN,*AIXWdN,*bhatt,*pbhat,*plamt;
  vector *S1,*korG,*pghat,*rowZ,*gam,*dgam,*ZHdN,*VZHdN,*IZHdN,*zi,*offsets;
  int it,i,j,k,l,c,s,count,pers=0,pmax;
  int stat, *ls=calloc(*Ntimes,sizeof(int)); 
  double S0,sumscore,time=0,dummy,dtime,random,fabs(),sqrt();
  double *weights=calloc(*antpers,sizeof(double)),
         *times=calloc(*Ntimes,sizeof(double)),
	 *cumoff=calloc((*Nalltimes)*(*px+1),sizeof(double));
  double norm_rand(); 
  void GetRNGstate(),PutRNGstate();  

  malloc_mats(*antpers,*px,&X,&WX,NULL);
  malloc_mats(*antpers,*pg,&Z,&WZ,NULL); 
  malloc_mats(*px,*px,&Vcov,&A,&AI,&GCdM1M2,&VarKorG,NULL); 
  malloc_mats(*pg,*pg,&S2,&S2I,&tmpM2,&ZWZ,&Vargam,&dVargam,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*antpers,&AIXW,NULL);
  malloc_mats(*px,*pg,&tmpM4,&tmpM3,&Ct,&dC,&XWZ,&XWZAI,&dM1M2,&M1M2t,NULL);

  for (j=0;j<*antpers;j++) { 
    malloc_mat(*Ntimes,*pg,Uti[j]); malloc_mat(*Ntimes,*pg,Utiid[j]); 
  }
  malloc_mat(*Ntimes,*pg,tmpM1); 
  malloc_mat(*Ntimes,*pg,Delta); 

  for (j=0;j<*Nalltimes;j++) {
    malloc_mat(*px,*pg,Acorb[j]);malloc_mat(*px,*pg,C[j]);}
  for (j=0;j<*Ntimes;j++) {malloc_mat(*px,*pg,M1M2[j]);
    malloc_mat(*pg,*pg,dUt[j]); }

  malloc_vecs(*px,&dA,&VdB,&difX,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&bhatt,NULL);
  malloc_vecs(*pg,&S1,&gamoff,&zi,&tmpv2,&rowZ,&gam,&dgam,&ZHdN,&IZHdN,&VZHdN,NULL);
  malloc_vecs(*antpers,&offsets,&dN,&pbhat,&pghat,&plamt,NULL);

  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  times[0]=alltimes[0]; 
  for (s=0;s<*pg;s++) VE(gam,s)=gamma[s]; 


  cu[0]=times[0]; 
  for (it=0;it<*Nit;it++)
    {
      mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN); l=0; sumscore=0; S0=0; 
      vec_zeros(gamoff); vec_zeros(offsets); 

      for (s=1;s<*Nalltimes;s++)
	{
	  time=alltimes[s]; dtime=time-alltimes[s-1]; 
	  mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); stat=0; 
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	    {
	      if ((start[c]<time) && (stop[c]>=time)) {
		if (*mof==1) VE(offsets,count)=offset[c]; 
		if (*mw==1)  weights[count]=weight[c]; else weights[count]=1; 
		for(j=0;j<pmax;j++) {
		  if (j<*px) ME(X,count,j)=designX[j*(*nx)+c];
		  if (j<*pg) ME(Z,count,j)=designG[j*(*ng)+c]; }
		if (time==stop[c] && status[c]==1) {pers=count;stat=1;l=l+1;
		  ls[l]=s;}
		count=count+1;}
	    }
	  Mv(Z,gam,pghat); 

//	  print_mat(X); print_mat(Z); 

	  vec_zeros(S1); mat_zeros(S2); S0=0; 
	  for (j=0;j<count;j++)
	    {dummy=exp(VE(pghat,j))*ME(X,j,0); S0=S0+dummy*weights[j]; 
	      extract_row(X,j,xi); 
	      scl_vec_mult(weights[j]*dummy,xi,rowX); replace_row(WX,j,rowX); 
	      extract_row(Z,j,zi); 
	      scl_vec_mult(weights[j]*dummy,zi,rowZ); replace_row(WZ,j,rowZ); 
	      for (k=0;k<*pg;k++) for (i=0;i<*pg;i++) 
		ME(S2,k,i)= ME(S2,k,i)+VE(zi,i)*VE(zi,k)*dummy*weights[j];
	      vec_add_mult(S1,zi,dummy,S1); 
	    }
	  scl_vec_mult(1/S0,S1,S1); 
	  scl_mat_mult(1/S0,S2,S2); 

	  MtA(X,WX,A); invertS(A,AI,1); 
	  MtA(Z,WZ,ZWZ);MtA(WX,Z,XWZ);MxA(AI,XWZ,XWZAI);
	  MtA(XWZAI,XWZ,tmpM2); mat_subtr(ZWZ,tmpM2,dCGam); 
	  scl_mat_mult(1/S0,dCGam,dCGam); 

	  if (stat==1) mat_add(CGam,dCGam,CGam); 

	  for (k=0;k<*pg;k++) 
	    for (j=0;j<*pg;j++) ME(S2,k,j)=ME(S2,k,j)-VE(S1,k)*VE(S1,j); 

	  /* correction from offsets calculated here for gamma komponent */
	  if (*mof==1)  {vM(X,offsets,rowX); Mv(AI,rowX,tmpv1); 
	    vM(Z,offsets,rowZ);  
	    vM(XWZAI,rowX,dgam); 
	    for (k=0;k<*pg;k++)  ME(dC,0,k)=dtime*VE(dgam,k)/S0; 
	    vec_subtr(rowZ,dgam,dgam); 
	    vec_add_mult(gamoff,dgam,dtime,gamoff); 
            for (k=1;k<=*px;k++) cumoff[k*(*Nalltimes)+s]=VE(tmpv1,k-1)*dtime; 
            scl_mat_mult(dtime*VE(rowX,0),dCGam,dCGam); 
	    mat_subtr(CGam,dCGam,CGam); 
	  }

	  if (stat==1) {
	    extract_row(X,pers,tmpv1); Mv(AI,tmpv1,AIXWdN); 
	    extract_row(Z,pers,zi); vM(XWZ,AIXWdN,tmpv2); 
	    vec_subtr(zi,tmpv2,ZHdN);
//	    scl_vec_mult(weights[pers],ZHdN,ZHdN); 
	    vec_add(ZHdN,IZHdN,IZHdN); 

	    scl_mat_mult(1/S0,XWZAI,XWZAI); 
	    mat_subtr(XWZAI,dC,dC); 
//	    printf(" %lf \n",S0); print_vec(tmpv1); print_vec(AIXWdN); 
	  }
	  mat_add(Ct,dC,Ct); C[s]=mat_copy(Ct,C[s]); 

	  if (it==*Nit-1) {

	    if (stat==1) {
	      vcu[l]=time; cu[l]=time; times[l]=time; 

	      for (k=0;k<*pg;k++) 
		{ for (j=0;j<*pg;j++) 
		    { ME(dVargam,k,j)= VE(ZHdN,j)*VE(ZHdN,k);
		      ME(dUt[l],k,j)=ME(CGam,k,j);   
		    }
		  replace_row(Uti[pers],l,ZHdN); 

		  for (j=0;j<*px;j++) ME(dM1M2,j,k)=VE(ZHdN,k)*VE(AIXWdN,j); }
	      mat_add(dVargam,Vargam,Vargam); mat_add(dM1M2,M1M2t,M1M2t); 
	      M1M2[l]=mat_copy(M1M2t,M1M2[l]); 

	      for (k=1;k<=*px;k++) {
		cu[k*(*Ntimes)+l]=VE(AIXWdN,k-1); 
		vcu[k*(*Ntimes)+l]=vcu[k*(*Ntimes)+l-1]+VE(AIXWdN,k-1)*VE(AIXWdN,k-1);}
	    }

	    for (k=0;k<*pg;k++) Ut[(k+1)*(*Nalltimes)+s]=VE(IZHdN,k)-VE(gamoff,k); 
	    Ut[s]=time; 
	  } /* it==*Nit-1 */ 
	} /* s =1...Ntimes */ 

      invertS(CGam,ICGam,1); 
      MxA(Vargam,ICGam,tmpM2); MxA(ICGam,tmpM2,Vargam); 
      vec_subtr(IZHdN,gamoff,IZHdN);  
      Mv(ICGam,IZHdN,rowZ); vec_add(gam,rowZ,gam);

      if (*detail==1) {
	Rprintf("====================Iteration %d ==================== \n",it);
	Rprintf("Estimate beta \n"); print_vec(gam);
	Rprintf("Score D l\n"); print_vec(IZHdN);
	Rprintf("Information -D^2 l\n");  print_mat(ICGam); };

      for (k=0;k<*pg;k++) sumscore=sumscore+fabs(VE(IZHdN,k));

      if ((fabs(sumscore)<0.0000000001) & (it<*Nit-2)) it=*Nit-2;
    } /* it<*Nit */


  l=0; 
  for (s=1;s<*Nalltimes;s++) {
    time=alltimes[s]; dtime=time-alltimes[s-1]; 
//    mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); 
    stat=0;  
    for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
    {
	if ((start[c]<time) && (stop[c]>=time)) {
	  if (*mof==1) VE(offsets,count)=offset[c]; 
	  if (*mw==1) weights[count]=weight[c]; else weights[count]=1; 
	  if (time==stop[c] && status[c]==1) {pers=count;stat=1;l=l+1;} 
	  count=count+1; 
	}
    }

    if (stat==1) 
    for (k=1;k<=*px;k++) VE(dA,k-1)=cu[k*(*Ntimes)+l]; else vec_zeros(dA); 
    if (*mof==1) 
    for (k=1;k<=*px;k++) VE(dA,k-1)=VE(dA,k-1)-cumoff[k*(*Nalltimes)+s];

    if (*mof==1) 
      for (k=1;k<=*px;k++) cumoff[k*(*Ntimes)+s]=cumoff[k*(*Ntimes)+s-1]+
	cumoff[k*(*Nalltimes)+s];

    if (stat==1) {
     for (k=1;k<=*px;k++){cu[k*(*Ntimes)+l]=
	     cu[k*(*Ntimes)+l-1]+cu[k*(*Ntimes)+l];
     if (*mof==1) cu[k*(*Ntimes)+l]=cu[k*(*Ntimes)+l]-
	                            cumoff[k*(*Nalltimes)+s];}

      MxA(C[ls[l]],Vargam,tmpM4); MAt(tmpM4,C[ls[l]],VarKorG);
      MxA(M1M2[l],ICGam,tmpM4); MAt(C[ls[l]],tmpM4,GCdM1M2);
      /* MxA(C[ls[l]],ICGam,tmpM4); MxA(tmpM4,M1M2[l],GCdM1M2);*/

      for (k=1;k<=*px;k++) {vcu[k*(*Ntimes)+l]= 
	  vcu[k*(*Ntimes)+l]+ME(VarKorG,k-1,k-1)-2*ME(GCdM1M2,k-1,k-1); } 
    }
  } /* s=1 ..Ntimes */ 

  for (j=0;j<*pg;j++) {gamma[j]=VE(gam,j); intZHdN[j]=VE(IZHdN,j); 
    for (k=0;k<*pg;k++) {Vgamma[k*(*pg)+j]=ME(Vargam,j,k);
      intZHZ[k*(*pg)+j]=ME(ICGam,j,k); }}
  cu[0]=times[0]; vcu[0]=times[0];

  if (*nsim>0) { 
//    Rprintf(" simulation starts, no resampling = %d \n",*nsim); 
    GetRNGstate();  /* to use R random normals */

    for (s=1;s<*Ntimes;s++)  
      for (i=0;i<*antpers;i++) {
	extract_row(Uti[i],s-1,gam); extract_row(Uti[i],s,ZHdN); 
	vec_add(gam,ZHdN,IZHdN); replace_row(Uti[i],s,IZHdN); 
	extract_row(Uti[i],s,ZHdN);  Mv(ICGam,IZHdN,ZHdN); 
	Mv(dUt[s],ZHdN,IZHdN); replace_row(Utiid[i],s,IZHdN); 
      }

    for (k=1;k<=*nsim;k++) {
      mat_zeros(Delta); 
      for (i=0;i<*antpers;i++) {
	random=norm_rand();
	scl_mat_mult(random,Utiid[i],tmpM1); mat_add(tmpM1,Delta,Delta); }

      for (s=1;s<*Ntimes;s++) { extract_row(Delta,s,zi); 
	for (i=0;i<*pg;i++) {VE(zi,i)=fabs(VE(zi,i)); 
	  if (VE(zi,i)>test[i*(*nsim)+k-1]) test[i*(*nsim)+k-1]=VE(zi,i); 
	} 
      }
    } 
    PutRNGstate();  /* to use R random normals */
  } /* if nsim >0 */ 

  free_mats(&X,&WX,&Z,&WZ,&Vcov,&A,&AI,&GCdM1M2,&S2,&S2I,&tmpM2,&ZWZ,&VarKorG,&Vargam,&dVargam,&ICGam,&CGam,&dCGam,&AIXW,&tmpM4,&tmpM3,&Ct,&dC,&XWZ,&XWZAI,&dM1M2,&M1M2t,&tmpM1,&Delta,NULL);

  for (j=0;j<*Nalltimes;j++) { free_mat(Acorb[j]); free_mat(C[j]);}
  for (j=0;j<*Ntimes;j++) {free_mat(M1M2[j]); free_mat(dUt[j]); }
  for (j=0;j<*antpers;j++) { free_mat(Uti[j]);  free_mat(Utiid[j]); }

  free_vecs(&dA,&VdB,&difX,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&bhatt,&S1,&gamoff,&zi,&tmpv2,&rowZ,&gam,&dgam,&ZHdN,&IZHdN,&VZHdN,&offsets,&dN,&pbhat,&pghat,&plamt,NULL);
  free(ls); free(times); free(cumoff); free(weights); 
}
