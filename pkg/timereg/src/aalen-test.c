//#include <stdio.h>
#include <math.h>
#include "matrix.h"

void robaalentest(times,Ntimes,designX,nx,p,antpers,start,stop,cu,vcu,robvcu,sim,antsim,retur,cumAit,test,rani,testOBS,status,Ut,simUt,id,weighted,robust,covariance,covs,resample,Biid,clusters,antclust,loglike,mof,offset,mw,weight,silent)
double *designX,*times,*start,*stop,*cu,*vcu,*robvcu,*cumAit,*test,*testOBS,*Ut,*simUt,*covs,*Biid,*loglike,*offset,*weight; 
int *nx,*p,*antpers,*Ntimes,*sim,*retur,*rani,*antsim,*status,*id,*covariance,
    *weighted,*robust,*resample,*clusters,*antclust,*mw,*mof,*silent;
{
  matrix *WX,*ldesignX,*A,*AI,*Vcov,*cumAt[*antclust];
  vector *diag,*dB,*dN,*VdB,*xi,*rowX,*rowcum,*difX,*vtmp,*cum,*offsets; 
  vector *vrisk,*cumhatA[*antclust],*cumA[*antclust];
  int i,j,k,l,s,c,count,pers=0;
  int stat,*cluster=calloc(*antpers,sizeof(int));
  double time,ahati,dtime;
  double *vcudif=calloc((*Ntimes)*(*p+1),sizeof(double)),
	 *weights=calloc(*antpers,sizeof(double)); 
  double fabs(),sqrt();

  if (*robust==1) {
    for (i=0;i<*antclust;i++) {
      malloc_vec(*p,cumhatA[i]); malloc_vec(*p,cumA[i]); malloc_mat(*Ntimes,*p,cumAt[i]);} }

  malloc_mat(*antpers,*p,ldesignX); malloc_mat(*antpers,*p,WX); 
  malloc_mats(*p,*p,&Vcov,&A,&AI,NULL); 
  malloc_vecs(*antpers,&vrisk,&dN,&offsets,NULL); 
  malloc_vecs(*p,&cum,&diag,&dB,&VdB,&xi,&rowX,&rowcum,&difX,&vtmp,NULL);
  for (j=0;j<*antpers;j++) cluster[j]=0;

  for (s=1;s<*Ntimes;s++)
  {
      time=times[s]; mat_zeros(ldesignX); dtime=time-times[s-1]; 
      mat_zeros(WX); stat=0; vec_zeros(vrisk); 

      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++)
      {
	  if ((start[c]<time) && (stop[c]>=time)) {
	    if (*mof==1) VE(offsets,id[c])=offset[c]; 
	    if (*mw==1)  weights[id[c]]=weight[c]; else weights[id[c]]=1; 
	    cluster[id[c]]=clusters[c]; 
	    VE(vrisk,id[c])=1.0; 
	    for(j=0;j<*p;j++) {ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	      ME(WX,id[c],j)=weights[id[c]]*designX[j*(*nx)+c];}
	    if (time==stop[c] && status[c]==1) {pers=id[c];stat=1;}
	    count=count+1; } 
      }

      if (*mof==1 || stat==1) {MtA(ldesignX,WX,A); invertS(A,AI,silent[0]);
      if (ME(AI,0,0)==0 && *silent==0) Rprintf("X'X not invertible at time %lf \n",time);
      if (s<-1) {print_mat(AI); print_mat(A);} }

      if (stat==1) {extract_row(WX,pers,xi); Mv(AI,xi,dB);} 
      else vec_zeros(dB);

      vec_star(dB,dB,VdB); 

     if (*mof==1) {vM(WX,offsets,rowX); Mv(AI,rowX,diag);
         scl_vec_mult(dtime,diag,diag); vec_subtr(dB,diag,dB); }

      vec_star(xi,dB,vtmp); ahati=vec_sum(vtmp); 
      if (*mof==1) ahati=ahati+dtime*VE(offsets,pers); 
      loglike[0]=loglike[0]+log(ahati); 
      loglike[1]=loglike[1]-ahati; 

      for (k=1;k<*p+1;k++) {
	cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dB,k-1); 
	vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s-1]+VE(VdB,k-1);
	VE(cum,k-1)=cu[k*(*Ntimes)+s];}

/*
      Rprintf(" %lf %lf  \n",time,dtime); 
      Rprintf(" %lf \n",vec_sum(offsets)); 
      print_mat(AI); print_vec(rowX); print_vec(xi); print_vec(dB);  
      print_vec(diag); 
      for (k=1;k<*p+1;k++) { Rprintf(" %lf  ",cu[k*(*Ntimes)+s]);  } 
      Rprintf(" \n "); 
      */

     robvcu[s]=time; cu[s]=time; vcu[s]=time;  

    if (*robust==1 || *retur==1) {
      vec_zeros(VdB); mat_zeros(Vcov);

	for (i=0;i<*antpers;i++)   // {{{
	{
           j=cluster[i]; 
	   extract_row(ldesignX,i,xi); 
	   Mv(AI,xi,rowX); if (*mw==1) scl_vec_mult(weights[i],rowX,rowX); 
	   vec_star(xi,dB,vtmp); ahati=vec_sum(vtmp); 
	   if (*mof==1) ahati=ahati+dtime*VE(offsets,i); 
	   /* loglike[0]=loglike[0]-ahati;    */
	   if (*robust==1) {
	   if (i==pers) vec_add(rowX,cumhatA[j],cumhatA[j]);
	    scl_vec_mult(ahati,rowX,rowX); vec_add(rowX,cumA[j],cumA[j]);

	   }

	   if (*retur==1) cumAit[i*(*Ntimes)+s]=cumAit[i*(*Ntimes)+s]+
	      weights[i]*(1*(i==pers)-ahati); 
	  
	} // i=0 ... antpers  // }}}

       if (*robust==1) { // {{{
       for (j=0;j<*antclust;j++) {
	  vec_subtr(cumhatA[j],cumA[j],difX); 
	  replace_row(cumAt[j],s,difX); 
	  vec_star(difX,difX,vtmp); vec_add(vtmp,VdB,VdB);

         if (*resample==1) { 
	    for (k=0;k<*p;k++) { l=j*(*p)+k; Biid[l*(*Ntimes)+s]=VE(difX,k); }}

	  if (*covariance==1) {
	    for (k=0;k<*p;k++) for (c=0;c<*p;c++)
	      ME(Vcov,k,c)=ME(Vcov,k,c)+VE(difX,k)*VE(difX,c);}

	} /* j in cluster */ 

	for (k=1;k<*p+1;k++) {robvcu[k*(*Ntimes)+s]=VE(VdB,k-1); 
	  if (*covariance==1) {
	    for (c=0;c<*p;c++)  {
	      l=(k-1)*(*p)+c; covs[l*(*Ntimes)+s]=ME(Vcov,k-1,c);}} }
      } /* if robust==1 */ // }}}
    }

    } /* s = 1..Ntimes */ 

  if (*sim==1) {
    comptest(times,Ntimes,p,cu,robvcu,vcudif,antsim,test,testOBS,Ut,simUt,cumAt,weighted,antclust);
  }

  cu[0]=times[0]; vcu[0]=times[0]; robvcu[0]=times[0]; 
  free_vecs(&offsets,&xi,&rowX,&diag,&dB,&VdB,&rowcum,&cum,&vtmp,NULL); 
  free_mats(&Vcov,&ldesignX,NULL);
  if (*robust==1)
    for (i=0;i<*antclust;i++) {free_vec(cumA[i]);free_vec(cumhatA[i]);free_mat(cumAt[i]);}
  free(weights); free(vcudif); free(cluster); 
} /* robaalentest  main */ 


void semiaalentest(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,nb,bhat,cu,vcu,Robvcu,gamma,Vgamma,RobVgamma,sim,antsim,test,rani,testOBS,robust,status,Ut,simUt,id,weighted,cumAit,retur,covariance,covs,resample,gammaiid,Biid,clusters,antclust,loglike,intZHZ,intZHdN,deltaweight,mof,offset,mw,weight,gamfix,
pseudoscore,pscoret,pscoretiid,intZHZt,ptest,intHdN,intHZ,varintHdN,silent)
double *designX,*alltimes,*start,*stop,*cu,*vcu,*bhat,*designG,*gamma,*Vgamma,*RobVgamma,*Robvcu,*test,*testOBS,*Ut,*simUt,*cumAit,*covs,*Biid,*gammaiid,*loglike,*intZHZ,*intZHdN,*offset,*weight,*pscoret,*pscoretiid,*intZHZt,*ptest,
*intHdN,*intHZ,*varintHdN;
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*nb,*ng,*pg,*sim,*antsim,*rani,*robust,*status,*id,*weighted,*retur,*covariance,*resample,*clusters,*antclust,*deltaweight,*mof,*mw,*gamfix,*pseudoscore,*silent;
{
  matrix *Vcov,*X,*WX,*A,*AI,*AIXW,*Z,*WZ,*tmpM1,*Delta,*iZHZt[*Ntimes];
  matrix *dCGam,*CGam,*Ct,*ICGam,*VarKorG,*dC,*ZH,*XWZ,*ZWZ,*XWZAI,*HZ; 
  matrix *Acorb[*Nalltimes],*Vargam,*dVargam,*M1M2[*Ntimes],*GCdM1M2; 
  matrix *C[*Nalltimes],*dM1M2,*M1M2t,*RobVargam,*tmpM2,*tmpM3,*tmpM4;
  matrix *W3t[*antclust],*W4t[*antclust],*AIxit[*antpers],*Scoret[*antclust]; 
  vector *W2[*antclust],*W3[*antclust]; 
  vector *dB,*VdB,*difX,*xi,*tmpv1,*tmpv2,*gamoff,*vrisk;
  vector *dAoff,*dA,*rowX,*dN,*AIXWdN,*bhatt,*pbhat,*plamt; 
  vector *korG,*pghat,*rowZ,*gam,*dgam,*ZHdN,*VZHdN,*IZHdN,*zi,*offsets;
  int m,i,j,k,l,c,s,count,pers=0,pmax,stat,
        *cluster=calloc(*antpers,sizeof(int)),
        *ipers=calloc(*Ntimes,sizeof(int)),
        *ls=calloc(*Ntimes,sizeof(int));
  double dptime,ddtime,time,dtime,random,fabs(),sqrt();
  double ahati,ghati,hati;
  double *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double)),
         *times=calloc((*Ntimes),sizeof(double)),
         *cumoff=calloc((*Nalltimes)*(*px+1),sizeof(double)),
         *cuL=calloc((*Nalltimes)*(*px+1),sizeof(double)),
         *weights=calloc((*antpers),sizeof(double));
  double norm_rand(); 
  dptime=alltimes[0]; 

  for (j=0;j<=*px;j++) cuL[j*(*Nalltimes)+0]=0.0; 
  if (*pseudoscore>=1) {
    malloc_mat(*Ntimes,*pg,tmpM1); malloc_mat(*Ntimes,*pg,Delta); 
    for (j=0;j<*Ntimes;j++) malloc_mat(*pg,*pg,iZHZt[j]);
    for (j=0;j<*antclust;j++) malloc_mat(*Ntimes,*pg,Scoret[j]);}

  malloc_mats(*antpers,*px,&X,&WX,NULL);
  malloc_mats(*antpers,*pg,&Z,&WZ,&HZ,NULL); 
  malloc_mats(*px,*px,&Vcov,&A,&AI,&GCdM1M2,&VarKorG,NULL); 
  malloc_mats(*pg,*pg,&tmpM2,&ZWZ,&RobVargam,&Vargam,&dVargam,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*antpers,&AIXW,NULL);
  malloc_mats(*pg,*antpers,&ZH,NULL);
  malloc_mats(*px,*pg,&tmpM4,&tmpM3,&Ct,&dC,&XWZ,&XWZAI,&dM1M2,&M1M2t,NULL);

  for (j=0;j<*Nalltimes;j++) {malloc_mat(*px,*pg,Acorb[j]);malloc_mat(*px,*pg,C[j]);}
  for (j=0;j<*Ntimes;j++) {malloc_mat(*px,*pg,M1M2[j]);}

  malloc_vecs(*px,&dAoff,&dA,&dB,&VdB,&difX,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&bhatt,NULL);
  malloc_vecs(*pg,&gamoff,&zi,&tmpv2,&rowZ,&gam,&dgam,&ZHdN,&IZHdN,&VZHdN,NULL);
  malloc_vecs(*antpers,&vrisk,&offsets,&dN,&pbhat,&pghat,&plamt,NULL);

  if (*robust==1) {
    for (j=0;j<*antclust;j++) {malloc_mat(*Ntimes,*px,W3t[j]); 
      malloc_mat(*Ntimes,*px,W4t[j]); malloc_vec(*pg,W2[j]);
      malloc_vec(*px,W3[j]); }
    for (j=0;j<*antpers;j++) malloc_mat(*Nalltimes,*px,AIxit[j]);  } 

  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN); times[0]=alltimes[0]; l=0; 
  for (s=0;s<*pg;s++) VE(gam,s)=gamma[s]; 
  for (j=0;j<*antpers;j++) cluster[j]=0;

  for (s=1;s<*Nalltimes;s++)
    {
      time=alltimes[s]; dtime=time-alltimes[s-1]; 
      mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); stat=0;  
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	{
	  if ((start[c]<time) && (stop[c]>=time)) {
	    if (*mof==1) VE(offsets,id[c])=offset[c]; 
	    if (*mw==1)  weights[id[c]]=weight[c]; else weights[id[c]]=1; 
	    for(j=0;j<pmax;j++) {
	      if (j<*px) ME(X,id[c],j)=designX[j*(*nx)+c];
	      if (j<*px) ME(WX,id[c],j)=weights[id[c]]*designX[j*(*nx)+c];
	      if (j<*pg) ME(Z,id[c],j)=designG[j*(*ng)+c]; 
	      if (j<*pg) ME(WZ,id[c],j)=weights[id[c]]*designG[j*(*ng)+c]; }
	    if (time==stop[c] && status[c]==1) {pers=id[c];stat=1;l=l+1;
	      ipers[l]=pers; ls[l]=s;}
	    count=count+1;}
	}

      MtA(X,WX,A); invertS(A,AI,silent[0]); 
      MtA(Z,WZ,ZWZ);MtA(X,WZ,XWZ);MxA(AI,XWZ,XWZAI);
      MtA(XWZAI,XWZ,tmpM2); mat_subtr(ZWZ,tmpM2,dCGam); 
      scl_mat_mult(dtime,dCGam,dCGam); 
      if (*deltaweight==0) scl_mat_mult(dtime,dCGam,dCGam); 
      mat_add(CGam,dCGam,CGam); 

      MxA(WX,XWZAI,HZ); 
      for (c=0;c<*antpers;c++) for (m=0;m<*pg;m++) 
      intHZ[m*(*antpers)+c]=intHZ[m*(*antpers)+c]+dtime*(ME(WZ,c,m)-ME(HZ,c,m));

      if (stat==1) { 
	extract_row(WX,pers,tmpv1); Mv(AI,tmpv1,AIXWdN); 
        extract_row(WZ,pers,zi); vM(XWZ,AIXWdN,tmpv2); 
        vec_subtr(zi,tmpv2,ZHdN);
	if (*deltaweight==0){scl_vec_mult(dtime,ZHdN,ZHdN);scl_vec_mult(dtime,VZHdN,VZHdN);}
        vec_add(ZHdN,IZHdN,IZHdN); 

        Mv(WX,AIXWdN,dN); 
	for (c=0;c<*antpers;c++) {intHdN[c]=intHdN[c]+
	    (pers==c)-VE(dN,c); 
	  for (m=0;m<*antpers;m++) 
	    varintHdN[c*(*antpers)+m]=varintHdN[c*(*antpers)+m]+
	      (pers==c)*(pers==m)*(1-2*VE(dN,pers))+VE(dN,c)*VE(dN,m); 
	} 

        if (*pseudoscore>=1) {
	  pscoret[l]=time; scl_mat_mult(1,CGam,iZHZt[l]); 
	  for (k=0;k<*pg;k++) {
	    pscoret[(k+1)*(*Ntimes)+l]=VE(IZHdN,k); 
	    for (c=0;c<*pg;c++)  {
	      m=k*(*pg)+c; intZHZt[m*(*Ntimes)+l]=ME(CGam,k,c); } }
	}
      } else vec_zeros(AIXWdN); 

      /* correction from offsets calculated here */
      if (*mof==1)  {vM(WX,offsets,rowX); Mv(AI,rowX,tmpv1); 
        scl_vec_mult(dtime,tmpv1,tmpv1); vec_subtr(AIXWdN,tmpv1,dB); 
	vM(WZ,offsets,rowZ);  vM(XWZAI,rowX,dgam); vec_subtr(rowZ,dgam,dgam); 
	vec_add_mult(gamoff,dgam,dtime,gamoff); 
	for (k=1;k<=*px;k++) cumoff[k*(*Nalltimes)+s]=VE(tmpv1,k-1); 
      } 
      
/*
      Rprintf(" %lf %lf \n",time,dtime); 
      Rprintf(" %lf \n",vec_sum(offsets)); 
      print_mat(AI); print_vec(rowX); 
      extract_row(WX,pers,rowX); 
      print_vec(rowX); 
      print_vec(dB);  
      print_vec(tmpv1); 
      for (k=1;k<*px+1;k++) {cuL[k*(*Nalltimes)+s]=cuL[k*(*Nalltimes)+s-1]+VE(dB,k-1); 
      Rprintf(" %lf  ",cuL[k*(*Nalltimes)+s]);  }
      Rprintf(" \n "); 
      */

      Acorb[s]=mat_copy(XWZAI,Acorb[s]); 

      scl_mat_mult(dtime,XWZAI,tmpM4);mat_add(tmpM4,Ct,Ct); 
      C[s]=mat_copy(Ct,C[s]); 

      if (stat==1) {
	vcu[l]=time; cu[l]=time; times[l]=time; 

	for (k=0;k<*pg;k++) 
	  { for (j=0;j<*pg;j++) ME(dVargam,k,j)= VE(ZHdN,j)*VE(ZHdN,k);
	    for (j=0;j<*px;j++) ME(dM1M2,j,k)=VE(ZHdN,k)*VE(AIXWdN,j); }
	mat_add(dVargam,Vargam,Vargam); mat_add(dM1M2,M1M2t,M1M2t); 
	M1M2[l]=mat_copy(M1M2t,M1M2[l]); 

	for (k=1;k<=*px;k++) {
	  cu[k*(*Ntimes)+l]=VE(AIXWdN,k-1); 
	  vcu[k*(*Ntimes)+l]=vcu[k*(*Ntimes)+l-1]+VE(AIXWdN,k-1)*VE(AIXWdN,k-1);}
      }

      if (*robust==1) for (i=0;i<*antpers;i++) {
	extract_row(WX,i,xi); Mv(AI,xi,rowX); replace_row(AIxit[i],s,rowX); }
    } /* s =1...Ntimes */ 

  invertS(CGam,ICGam,silent[0]); 
  if (*gamfix==0) Mv(ICGam,IZHdN,gam); 
  if ((*mof==1) & (*gamfix==0)) {Mv(ICGam,gamoff,dgam); vec_subtr(gam,dgam,gam);}
  MxA(Vargam,ICGam,tmpM2); MxA(ICGam,tmpM2,Vargam); 
  // if (*gamfix==1) mat_zeros(Vargam); 

  l=0; vec_zeros(dAoff); 
  for (s=1;s<*Nalltimes;s++) {
    time=alltimes[s]; vec_zeros(dN);dtime=time-alltimes[s-1]; 
    mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ); stat=0; 
    for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
      {
	if ((start[c]<time) && (stop[c]>=time)) {
	  cluster[id[c]]=clusters[c]; 
	  if (*mof==1) VE(offsets,id[c])=offset[c]; 
	  if (*mw==1)  weights[id[c]]=weight[c]; else weights[id[c]]=1; 
	  for(j=0;j<pmax;j++) {
	    if (j<*px) ME(X,id[c],j)=designX[j*(*nx)+c];
	    if (j<*pg) ME(Z,id[c],j)=designG[j*(*ng)+c];  }
	  if (time==stop[c] && status[c]==1) {pers=id[c];stat=1;l=l+1;
		  ghati=0;} 
	  count=count+1; }
      }
    vec_zeros(dA); 

    if (stat==1) {for (k=1;k<=*px;k++) VE(dA,k-1)=cu[k*(*Ntimes)+l];}
    if (*mof==1) {for (k=1;k<=*px;k++) { VE(dA,k-1)=VE(dA,k-1)-cumoff[k*(*Nalltimes)+s];
                   VE(dAoff,k-1)=VE(dAoff,k-1)-cumoff[k*(*Nalltimes)+s];
                   cumoff[k*(*Nalltimes)+s]=cumoff[k*(*Nalltimes)+s-1]+cumoff[k*(*Nalltimes)+s];
		  }}

    /* terms for robust variance   */ 
    if (*robust==1 || *retur==1) 
    {
	mat_subtr(C[s],C[s-1],tmpM4); Mv(tmpM4,gam,korG); 

	  for (i=0;i<*antpers;i++) 
	  {
	    j=cluster[i]; 
	    extract_row(X,i,xi); extract_row(Z,i,zi); 
	    vec_star(xi,dA,tmpv1); ahati=vec_sum(tmpv1);
	    vec_star(zi,gam,rowZ); ghati=dtime*vec_sum(rowZ); 
	    vec_star(xi,korG,tmpv1); 
	    hati=ahati+ghati-vec_sum(tmpv1);
	    if (*mof==1) hati=hati+dtime*VE(offsets,i); 
	      /* loglike[0]=loglike[0]-hati;  */ 
            if (*robust==1) {
	      vM(Acorb[s],xi,tmpv2); vec_subtr(zi,tmpv2,tmpv2); 
	      if (*mw==1) scl_vec_mult(weights[i],tmpv2,tmpv2); 

	      if (i==pers && stat==1) vec_add(tmpv2,W2[j],W2[j]);
	      scl_vec_mult(hati,tmpv2,rowZ); vec_subtr(W2[j],rowZ,W2[j]); 
	      if (*pseudoscore>=1 && stat==1) {
		for (k=0;k<*pg;k++) {
		  m=j*(*pg)+k; pscoretiid[m*(*Ntimes)+l]=VE(W2[j],k);} 
		replace_row(Scoret[j],l,W2[j]); 
	      }

	      extract_row(AIxit[i],s,rowX); 
	      if (i==pers && stat==1) vec_add(rowX,W3[j],W3[j]); 
	      scl_vec_mult(hati,rowX,rowX); vec_subtr(W3[j],rowX,W3[j]); 
	    }

	    if (*retur==1) {
	      if (stat==0) cumAit[i*(*Ntimes)+l+1]=
		           cumAit[i*(*Ntimes)+l+1]+1*(i==pers)*stat-hati;
               else cumAit[i*(*Ntimes)+l]=
		           cumAit[i*(*Ntimes)+l]+1*(i==pers)*stat-hati;
	      }
	    // if (*gamfix==1) vec_zeros(W2[j]); 
	  if (stat==1) replace_row(W3t[j],l,W3[j]);  
	} /* j=1..antclust */ 
      } /* robust ==1 */


    if (stat==1) {
      ddtime=time-dptime; 
      extract_row(X,pers,xi); vec_star(xi,dA,rowX); ahati=vec_sum(rowX); 
      extract_row(Z,pers,zi); vec_star(zi,gam,rowZ); 
      ahati=ahati+vec_sum(rowZ)*ddtime; 
      if (*mof==1) ahati=ahati+VE(offsets,pers)*ddtime; 
      loglike[0]=loglike[0]+log(ahati); 
      if (*deltaweight==1) loglike[1]=loglike[1]-ahati/dtime; 
      if (*deltaweight==0) loglike[1]=loglike[1]-ahati; 
      dptime=time; 

     for (k=1;k<=*px;k++)  cu[k*(*Ntimes)+l]=cu[k*(*Ntimes)+l-1]+
                              cu[k*(*Ntimes)+l]+VE(dAoff,k-1); 
     vec_zeros(dAoff); 

      MxA(C[ls[l]],Vargam,tmpM4); MAt(tmpM4,C[ls[l]],VarKorG);
      MxA(M1M2[l],ICGam,tmpM4); MAt(C[ls[l]],tmpM4,GCdM1M2);

      for (k=1;k<=*px;k++) {vcu[k*(*Ntimes)+l]= 
	  vcu[k*(*Ntimes)+l]+ME(VarKorG,k-1,k-1)-2*ME(GCdM1M2,k-1,k-1); } 
    }
  } /* s=1 ..Ntimes */ 

  vec_star(IZHdN,gam,rowZ);    

  loglike[1]=loglike[1]-vec_sum(rowZ); 

  for (s=1;s<*Ntimes;s++) {

    Mv(C[ls[s]],gam,korG); 
    for (k=1;k<=*px;k++) cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s]-VE(korG,k-1); 

    /* ROBUST VARIANCES   */ 
    if (*robust==1) {
      vec_zeros(VdB); mat_zeros(Vcov); 

      for (j=0;j<*antclust;j++) {
        Mv(ICGam,W2[j],tmpv2);  
        if (*gamfix==0) vec_zeros(rowX);  else Mv(C[ls[s]],tmpv2,rowX); 
	extract_row(W3t[j],s,tmpv1); 
	vec_subtr(tmpv1,rowX,difX); 
	replace_row(W4t[j],s,difX); 
	vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	if (*resample==1) {
	  if (s==1) for (k=0;k<*pg;k++) gammaiid[k*(*antclust)+j]=VE(tmpv2,k);
	  for (k=0;k<*px;k++) {
	    l=j*(*px)+k; Biid[l*(*Ntimes)+s]=VE(difX,k);} }

	if (*covariance==1) {
	  for (k=0;k<*px;k++) for (c=0;c<*px;c++)
	    ME(Vcov,k,c)=ME(Vcov,k,c)+VE(difX,k)*VE(difX,c);}

	if (s==1) { for (c=0;c<*pg;c++) for (k=0;k<*pg;k++) 
		      ME(RobVargam,c,k)=ME(RobVargam,c,k)+VE(W2[j],c)*VE(W2[j],k);}
      } /* j =1 ..Antclust */ 
    } /* robust==1 */

    for (k=1;k<*px+1;k++) { Robvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      if (*covariance==1) {
	for (j=0;j<*px;j++)  {
	  l=(k-1)*(*px)+j; covs[l*(*Ntimes)+s]=ME(Vcov,k-1,j); }}}
  }  /*  s=1 ..Ntimes */ 

  if (*robust==1) {
    MxA(RobVargam,ICGam,tmpM2); MxA(ICGam,tmpM2,RobVargam); 
  }

  for (j=0;j<*pg;j++) {gamma[j]=VE(gam,j); intZHdN[j]=VE(IZHdN,j); 
    for (k=0;k<*pg;k++) {Vgamma[k*(*pg)+j]=ME(Vargam,j,k);
      RobVgamma[k*(*pg)+j]=ME(RobVargam,j,k);
      intZHZ[k*(*pg)+j]=ME(CGam,j,k); }}

  cu[0]=times[0]; vcu[0]=times[0];

  if (*sim==1) {
    comptest(times,Ntimes,px,cu,Robvcu,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antclust); 
  }

  if (*pseudoscore>=1) {
    GetRNGstate();  /* to use R random normals */

    for (s=1;s<*Ntimes;s++) {
      for (k=0;k<*pg;k++) VE(tmpv2,k)=pscoret[(k+1)*(*Ntimes)+s];
      Mv(iZHZt[s],gam,zi); 
      vec_subtr(tmpv2,zi,tmpv2); 
      for (k=0;k<*pg;k++) pscoret[(k+1)*(*Ntimes)+s]=VE(tmpv2,k);
    }

    for (j=0;j<*antclust;j++) {
      for (k=0;k<*pg;k++) VE(tmpv2,k)=gammaiid[k*(*antclust)+j];
      for (s=1;s<*Ntimes;s++) {
	extract_row(Scoret[j],s,rowZ); 
	Mv(iZHZt[s],tmpv2,zi); 
	vec_subtr(rowZ,zi,rowZ); replace_row(Scoret[j],s,rowZ); } }

    for (k=1;k<=*pseudoscore;k++) {
      mat_zeros(Delta); 
      for (i=0;i<*antclust;i++) {
	/* random=gasdev(&idum); */ 
	random=norm_rand();
	scl_mat_mult(random,Scoret[i],tmpM1); mat_add(tmpM1,Delta,Delta);
      }

      for (s=1;s<*Ntimes;s++) { 
	extract_row(Delta,s,rowZ);
       
	for (i=0;i<*pg;i++) {
	  VE(rowZ,i)=fabs(VE(rowZ,i));
	  if (VE(rowZ,i)>ptest[i*(*pseudoscore)+k]) 
	    ptest[i*(*pseudoscore)+k]=VE(rowZ,i);
	} }

    }
    PutRNGstate();  /* to use R random normals */
  }

  free_mats(&HZ,&X,&WX,&Z,&WZ, &Vcov,&A,&AI,&GCdM1M2,
	      &tmpM2,&ZWZ,&VarKorG,&RobVargam,&Vargam,&dVargam,&ICGam,&CGam,&dCGam,
	      &AIXW, &ZH, &tmpM4,&tmpM3,&Ct,&dC,&XWZ,&XWZAI,&dM1M2,&M1M2t,NULL);

  for (j=0;j<*Nalltimes;j++) {free_mat(Acorb[j]); free_mat(C[j]); }
  for (j=0;j<*Ntimes;j++) {free_mat(M1M2[j]);}

  free_vecs(&vrisk,&dAoff,&dB,&dA,&VdB,&difX,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&bhatt,
	      &gamoff,&zi,&tmpv2,&rowZ,&gam,&dgam,&ZHdN,&IZHdN,&VZHdN,
	      &offsets,&dN,&pbhat,&pghat,&plamt,NULL);

  if (*robust==1) {
    for (j=0;j<*antclust;j++) {free_mat(W3t[j]); 
      free_mat(W4t[j]); free_vec(W2[j]); free_vec(W3[j]); }
    for (j=0;j<*antpers;j++) free_mat(AIxit[j]); } 
  if (*pseudoscore>=1) {
    free_mats(&tmpM1,&Delta,NULL); 
    for (j=0;j<*Ntimes;j++) free_mat(iZHZt[j]); 
    for (j=0;j<*antclust;j++) free_mat(Scoret[j]); }
  free(vcudif); free(times); free(cumoff); free(cuL); free(ipers);
  free(cluster); free(weights); 
}
