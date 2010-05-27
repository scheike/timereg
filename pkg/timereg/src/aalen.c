#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

void aalen(times,Ntimes,designX,nx,p,antpers,start,stop,cu,vcu,status)
double *designX,*times,*start,*stop,*cu,*vcu;
int *nx,*p,*antpers,*Ntimes,*status;
{
  matrix *ldesignX, *A, *AI;
  vector *dB, *VdB, *tmpv, *xi;
  int j,k,s,c,count,pers=0;
  double time; 

  malloc_mat(*antpers,*p,ldesignX);
  malloc_mat(*p,*p,A);
  malloc_mat(*p,*p,AI);
  
  malloc_vec(*p,xi);
  malloc_vec(*p,dB);
  malloc_vec(*p,VdB);
  malloc_vec(*p,tmpv);

  for (s=1;s<*Ntimes;s++){
    time=times[s]; mat_zeros(ldesignX);

    for (c=0,count=0;((c<*nx) && (count!=*antpers));c++)
    {
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<*p;j++){ ME(ldesignX,count,j) = designX[j*(*nx)+c]; }
	  if (time==stop[c] && status[c]==1) {
	    pers=count; for(j=0;j<*p;j++) { VE(xi,j)=designX[j*(*nx)+c]; }
	  }
	  count=count+1; 
	} 
    }
    //readXt2(antpers,nx,p,designX,start,stop,status,pers,ldesignX,time); 

    extract_row(ldesignX,pers,xi); 
    MtM(ldesignX,A); invert(A,AI); 
      
    Mv(AI,xi,dB); vec_star(dB,dB,VdB); 
      
    if (vec_sum(dB)==0.0){
      printf("Aalen:Singular matrix for time=%lf \n",time); 
    }
   
    cu[s]=time; 
    vcu[s]=time; 
    for (k=1;k<*p+1;k++) {
      cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dB,k-1);
      vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s-1]+VE(VdB,k-1);
    }
  }
  cu[0]=times[0]; vcu[0]=times[0]; 

  free_vec(dB); free_vec(VdB); free_mat(ldesignX);
  free_mat(A); free_mat(AI); free_vec(xi); free_vec(tmpv);
}

void robaalen(times,Ntimes,designX,nx,p,antpers,start,stop,cu,vcu,
	      robvcu,sim,antsim,retur,cumAit,test,rani,testOBS,status,
	      Ut,simUt,id,weighted,robust,covariance,covs,resample,
	      Biid,clusters,antclust,loglike,silent) 
double *designX,*times,*start,*stop,*cu,*vcu,*robvcu,*cumAit,*test,*testOBS,*Ut,*simUt,*covs,*Biid,*loglike; 
int *nx,*p,*antpers,*Ntimes,*sim,*retur,*rani,*antsim,*status,*id,*covariance,
    *weighted,*robust,*resample,*clusters,*antclust,*silent;
{

  matrix *ldesignX, *QR, *R, *A, *AI, *Vcov;
  matrix *cumAt[*antclust];
  vector  *diag,*dB,*dN,*VdB,*xi,*rowX,*rowcum,*difX,*vtmp;
  vector *cumhatA[*antclust],*cumA[*antclust],*cum;
  int ci,i,j,k,l,s,c,count,pers=0,*cluster=calloc(*antpers,sizeof(int));
  int var1,var2;
  double time,ahati,dt,tau,*vcudif=calloc((*Ntimes)*(*p+1),sizeof(double));
  double fabs(),sqrt();
  long idum; 
  
  *rani = -8021; idum=*rani; 
  dt=times[*Ntimes-1]-times[0]; var1=1; var2=0; 

  if (*robust==1) {
    for (i=0;i<*antclust;i++) { malloc_vec(*p,cumhatA[i]); 
    malloc_vec(*p,cumA[i]); 
    if (*sim==1) malloc_mat(*Ntimes,*p,cumAt[i]); } }

/*   print_clock(&debugTime, 0); */

  malloc_vec(*p,cum); malloc_mat(*antpers,*p,ldesignX); malloc_mat(*p,*p,QR);
  malloc_mat(*p,*p,Vcov); malloc_mat(*p,*p,A); malloc_mat(*p,*p,AI);
  malloc_vec(*antpers,dN); malloc_mat(*antpers,*p,R);

  malloc_vec(*p,diag); malloc_vec(*p,dB); malloc_vec(*p,VdB); malloc_vec(*p,xi);
  malloc_vec(*p,rowX); malloc_vec(*p,rowcum); malloc_vec(*p,difX); malloc_vec(*p,vtmp);

  for (j=0;j<*antpers;j++) cluster[j]=0;

/*   print_clock(&debugTime, 1); */

  R_CheckUserInterrupt();

  for (s=1;s<*Ntimes;s++){
    time=times[s]; mat_zeros(ldesignX); 

    for (c=0,count=0;((c<*nx) && (count!=*antpers));c++){
      if ((start[c]<time) && (stop[c]>=time)) {
	for(j=0;j<*p;j++) {
	  ME(ldesignX,id[c],j) = designX[j*(*nx)+c]; }
	  cluster[id[c]]=clusters[c]; 
	if (time==stop[c] && status[c]==1) { pers=id[c]; }
	count=count+1; } 
    }
    
// readXt(antpers,nx,p,designX,start,stop,status,pers,ldesignX,time,clusters,cluster,id);

    MtM(ldesignX,A); 
    invertS(A,AI,silent[0]); 

    if (ME(AI,0,0)==0.0 && *silent==0){ 
       printf(" X'X not invertible at time %lf \n",time); }
       if (s < -1) { print_mat(AI); print_mat(A);	}

    extract_row(ldesignX,pers,xi);
      
    Mv(AI,xi,dB); vec_star(dB,dB,VdB); 

    vec_star(xi,dB,vtmp); 
    ahati = vec_sum(vtmp);
    loglike[0]=loglike[0]-ahati/(time-times[s-1]); 

    for (k=1;k<*p+1;k++) {
      cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dB,k-1);
      vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s-1]+VE(VdB,k-1);
      VE(cum,k-1)=cu[k*(*Ntimes)+s];
    }
    cu[s]=time; vcu[s]=time; robvcu[s]=time; 

    if (*robust==1 || *retur==1) {
      vec_zeros(VdB); mat_zeros(Vcov);

	for (i=0;i<*antpers;i++)
	{
          ci=cluster[i]; 
	  extract_row(ldesignX,i,xi); ahati=vec_prod(xi,dB);
	  Mv(AI,xi,rowX); 
	  if (*robust==1) {
	  if (i==pers) { vec_add(rowX,cumhatA[ci],cumhatA[ci]); }
	    scl_vec_mult(ahati,rowX,rowX);
	    vec_add(rowX,cumA[ci],cumA[ci]);
	  }

	  if (*retur==1){ 
	     cumAit[i*(*Ntimes)+s]= cumAit[i*(*Ntimes)+s]+1*(i==pers)-ahati;  
	  }
       }

       if (*robust==1) {
       for (i=0;i<*antclust;i++) {

	   vec_subtr(cumhatA[i],cumA[i],difX);
	   if (*sim==1) replace_row(cumAt[i],s,difX);
	   vec_star(difX,difX,vtmp); vec_add(vtmp,VdB,VdB);

	   if (*resample==1) {
	     for (k=0;k<*p;k++) {l=i*(*p)+k; Biid[l*(*Ntimes)+s]=VE(difX,k);}
	   }

	   if (*covariance==1) {
	     for (k=0;k<*p;k++) for (c=0;c<*p;c++)
	      ME(Vcov,k,c) = ME(Vcov,k,c) + VE(difX,k)*VE(difX,c); }

       }
           for (k=1;k<*p+1;k++) { 
	       robvcu[k*(*Ntimes)+s]=VE(VdB,k-1); 
	   if (*covariance==1) {
	       for (c=0;c<*p;c++)  {
	       l=(k-1)*(*p)+c; 
	       covs[l*(*Ntimes)+s]=ME(Vcov,k-1,c);
	       }
	   }
           }
     }
    } /* if robust==1  || retur==1*/ 


    R_CheckUserInterrupt();

  } /* s = 1..Ntimes */ 

  tau=time; 

  if (*sim==1) {
    comptest(times,Ntimes,p,cu,robvcu,vcudif,antsim,test,testOBS,Ut,simUt,cumAt,weighted,antclust);
  }

  cu[0]=times[0]; vcu[0]=times[0]; robvcu[0]=times[0]; 
  free_vec(xi); free_vec(rowX); free_vec(diag); free_vec(dB); free_vec(VdB);
  free_vec(rowcum); free_vec(cum); free_vec(vtmp); free_mat(Vcov);
  free_mat(ldesignX); free_mat(QR);

  if (*robust==1){
    for (i=0;i<*antclust;i++) {
      free_vec(cumA[i]); free_vec(cumhatA[i]); if (*sim==1) free_mat(cumAt[i]); } }
  free(cluster); 
  free(vcudif); 
}

void semiaalen(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,nb,bhat,cu,vcu,Robvcu,gamma,Vgamma,RobVgamma,sim,antsim,test,rani,testOBS,robust,status,Ut,simUt,id,weighted,cumAit,retur,covariance,covs,resample,gammaiid,Biid,clusters,antclust,loglike,intZHZ,intZHdN,deltaweight,silent)
double *designX,*alltimes,*start,*stop,*cu,*vcu,*bhat,*designG,*gamma,*Vgamma,*RobVgamma,*Robvcu,*test,*testOBS,*Ut,*simUt,*cumAit,*covs,*Biid,*gammaiid,*loglike,*intZHZ,*intZHdN; 
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*nb,*ng,*pg,*sim,*antsim,*rani,*robust,*status,*id,*weighted,*retur,*covariance,*resample,*clusters,*antclust,*deltaweight,*silent;
{
  matrix *Vcov,*X,*WX,*A,*AI,*AIXW,*Z,*WZ;
  matrix *dCGam,*CGam,*Ct,*ICGam,*VarKorG,*dC,*ZH,*XWZ,*ZWZ,*XWZAI;
  matrix *Acorb[*Nalltimes],*Vargam,*dVargam,*M1M2[*Ntimes],*GCdM1M2;
  matrix *C[*Nalltimes],*dM1M2,*M1M2t,*RobVargam; 
  matrix *tmpM2,*tmpM3,*tmpM4;
  matrix *W3t[*antclust],*W4t[*antclust];
  vector *W2[*antclust],*W3[*antclust];
  matrix *AIxit[*antpers]; 
  vector *VdB,*difX,*xi,*tmpv1,*tmpv2; 
  vector *dA,*rowX,*dN,*AIXWdN,*bhatt,*pbhat,*plamt;
  vector *korG,*pghat,*rowZ,*gam,*dgam,*ZHdN,*IZHdN,*zi;
  int ci,i,j,k,l,c,s,count,pers=0,pmax,stat,maxtime,
      *cluster=calloc(*antpers,sizeof(int)),
      *ls=calloc(*Ntimes,sizeof(int)); 
  double time,dtime,dtime1,fabs(),sqrt();
  double ahati,ghati,hati,tau,dMi;
  double *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double)),
	 *times=calloc(*Ntimes,sizeof(double));
  long idum; 
  idum=*rani; 

  malloc_mat(*antpers,*px,X); malloc_mat(*antpers,*px,WX); malloc_mat(*antpers,*pg,Z);
  malloc_mat(*antpers,*pg,WZ); malloc_mat(*px,*px,Vcov); malloc_mat(*px,*px,A);
  malloc_mat(*px,*px,AI); malloc_mat(*px,*px,GCdM1M2); malloc_mat(*pg,*pg,tmpM2);
  malloc_mat(*pg,*pg,ZWZ); malloc_mat(*px,*px,VarKorG); malloc_mat(*pg,*pg,RobVargam);
  malloc_mat(*pg,*pg,Vargam); malloc_mat(*pg,*pg,dVargam); malloc_mat(*pg,*pg,ICGam);
  malloc_mat(*pg,*pg,CGam); malloc_mat(*pg,*pg,dCGam); malloc_mat(*px,*antpers,AIXW);
  malloc_mat(*pg,*antpers,ZH);
  malloc_mats(*px,*pg,&tmpM3,&Ct,&dC,&XWZ,&XWZAI,&dM1M2,&M1M2t,NULL);
  malloc_mat(*px,*pg,tmpM4);
  for (j=0;j<*Nalltimes;j++) {
	  malloc_mat(*px,*pg,Acorb[j]);
	  malloc_mat(*px,*pg,C[j]);
  }
  for (j=0;j<*Ntimes;j++) {
	  malloc_mat(*px,*pg,M1M2[j]);
  }

  malloc_vec(*px,dA); malloc_vec(*px,VdB); malloc_vec(*px,difX);
  malloc_vec(*px,xi); malloc_vec(*px,tmpv1); malloc_vec(*px,korG);
  malloc_vec(*px,rowX); malloc_vec(*px,AIXWdN); malloc_vec(*px,bhatt);
  malloc_vec(*pg,zi); malloc_vec(*pg,tmpv2); malloc_vec(*pg,rowZ);
  malloc_vec(*pg,gam); malloc_vec(*pg,dgam); malloc_vec(*pg,ZHdN);
  malloc_vec(*pg,IZHdN); malloc_vec(*antpers,dN); malloc_vec(*antpers,pbhat);
  malloc_vec(*antpers,pghat); malloc_vec(*antpers,plamt); 

  if (*robust==1) {
	  for (j=0;j<*antclust;j++) { malloc_mat(*Ntimes,*px,W3t[j]);
		  malloc_mat(*Ntimes,*px,W4t[j]); malloc_vec(*pg,W2[j]); 
		  malloc_vec(*px,W3[j]); }
    for (j=0;j<*antpers;j++){ 
      malloc_mat(*Nalltimes,*px,AIxit[j]);
    }; 
  }
  for (j=0;j<*antpers;j++) cluster[j]=0;

  pmax=max(*pg,*px); 
  mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN);
  times[0]=alltimes[0]; l=0; 
  maxtime=alltimes[*Nalltimes]; 

  for (s=1;s<*Nalltimes;s++){
      time=alltimes[s]; dtime=time-alltimes[s-1]; 
      mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ);
      stat=0;  

      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<pmax;j++) {
	    if (j<*px) { ME(X,id[c],j)=designX[j*(*nx)+c]; }
	    if (j<*px) { ME(WX,id[c],j) =designX[j*(*nx)+c]; }
	    if (j<*pg) { ME(Z,id[c],j)=designG[j*(*ng)+c]; }
	    if (j<*pg) { ME(WZ,id[c],j)=designG[j*(*ng)+c]; } 
	  }
	  if (time==stop[c] && status[c]==1) {
	    pers=id[c];stat=1;l=l+1; ls[l]=s;
	  }
	  count=count+1; 
	}
      }
      

//readXZt(antpers,nx,px,designX,pg,designG,start,stop,status,pers,X,WX,Z,WZ,time,
//	clusters,cluster,ls,stat,l,id,s,1);

      MtA(X, WX,A); invertS(A,AI,silent[0]);
      if (ME(AI,0,0)==0.0 && *silent==0){ 
	printf(" X'X not invertible at time %lf \n",time);
      }
      MtA(Z, WZ,ZWZ); MtA(X,WZ,XWZ); MxA(AI,XWZ,XWZAI);
      MtA(XWZAI,XWZ,tmpM2);
      mat_subtr(ZWZ,tmpM2,dCGam);
      scl_mat_mult(dtime,dCGam,dCGam);
      if (*deltaweight==0) { scl_mat_mult(dtime,dCGam,dCGam); }
      mat_add(CGam,dCGam,CGam);

      if (stat==1) {
	extract_row(WX,pers,tmpv1); Mv(AI,tmpv1,AIXWdN);
	extract_row(WZ,pers,zi);
	vM(XWZ,AIXWdN,tmpv2); vec_subtr(zi,tmpv2,ZHdN);
	if (*deltaweight==0){ scl_vec_mult(dtime,ZHdN,ZHdN); }
	vec_add(ZHdN,IZHdN,IZHdN);
      }
      mat_copy(XWZAI,Acorb[s]);

      scl_mat_mult(dtime,XWZAI,tmpM4);
      mat_add(tmpM4,Ct,Ct); mat_copy(Ct,C[s]);

      if (stat==1) {
	vcu[l]=time; cu[l]=time; times[l]=time; 

	for (k=0;k<*pg;k++){ 
	  for (j=0;j<*pg;j++){ 
	    ME(dVargam,k,j)= VE(ZHdN,j)*VE(ZHdN,k);
	  }
	  for (j=0;j<*px;j++){ 
	    ME(dM1M2,j,k)=VE(ZHdN,k)*VE(AIXWdN,j);
	  } 
	}
	mat_add(dVargam,Vargam,Vargam);
	mat_add(dM1M2,M1M2t,M1M2t);
	mat_copy(M1M2t,M1M2[l]);

	for (k=1;k<=*px;k++) {
	  cu[k*(*Ntimes)+l]=VE(AIXWdN,k-1); 
	  vcu[k*(*Ntimes)+l]=vcu[k*(*Ntimes)+l-1]+VE(AIXWdN,k-1)*VE(AIXWdN,k-1);
	}
      }

      if (*robust==1) {
	for (i=0;i<*antpers;i++) {
	  extract_row(X,i,xi); Mv(AI,xi,rowX); replace_row(AIxit[i],s,rowX); 
	}
      }
    } /* s =1...Ntimes */ 


  invertS(CGam,ICGam,silent[0]);
  if (ME(ICGam,0,0)==0 && *silent==0) printf(" intZHZ  singular\n"); 
  Mv(ICGam,IZHdN,gam); 
  MxA(Vargam,ICGam,tmpM2); 
  MxA(ICGam,tmpM2,Vargam); 

  l=0; 
  for (s=1;s<*Nalltimes;s++) {
    time=alltimes[s]; vec_zeros(dN); dtime=time-alltimes[s-1]; 
    mat_zeros(X); mat_zeros(Z); mat_zeros(WX); mat_zeros(WZ);
    stat=0; dMi=0; 

//    readXZt(antpers,nx,px,designX,pg,designG,start,stop,status,pers,X,WX,Z,WZ,time,
//	clusters,cluster,ls,stat,l,id,s,0);

    for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time)) {
	  cluster[id[c]]=clusters[c]; 
	  for(j=0;j<pmax;j++) {
	    if (j<*px){ ME(X,id[c],j)=designX[j*(*nx)+c]; }
	    if (j<*pg){ ME(Z,id[c],j)=designG[j*(*ng)+c]; }  
	  }
	  if (time==stop[c] && status[c]==1) {pers=id[c];stat=1;l=l+1;} 
	  count=count+1; }
      }

    if (stat==1) {
      for (k=1;k<=*px;k++) { VE(dA,k-1)=cu[k*(*Ntimes)+l]; }
    } else { vec_zeros(dA); }


    /* terms for robust variance   */ 
    if (*robust==1 || *retur==1) 
    {
	mat_subtr(C[s],C[s-1],tmpM4); 
	Mv(tmpM4,gam,korG); 

	for (i=0;i<*antpers;i++) 
	{
          ci=cluster[i]; 
	  extract_row(X,i,xi); extract_row(Z,i,zi); 
	  if (stat==1) {
	     vec_star(xi,dA,tmpv1); 
	     ahati=vec_sum(tmpv1);
	  } else  { ahati=0.0; }
	  vec_star(zi,gam,rowZ); 
	  ghati=dtime*vec_sum(rowZ); 
	  vec_star(xi,korG,tmpv1); 
	  hati=ahati+ghati-vec_sum(tmpv1);

	  if (*robust==1) {
	  vM(Acorb[s],xi,tmpv2); 
	  vec_subtr(zi,tmpv2,tmpv2); 

	  if (i==pers && stat==1){ 
	     vec_add(tmpv2,W2[ci],W2[ci]);
	  }
	  scl_vec_mult(hati,tmpv2,rowZ); 
	  vec_subtr(W2[ci],rowZ,W2[ci]); 

	  extract_row(AIxit[i],s,rowX); 
	  if (i==pers && stat==1) { vec_add(rowX,W3[ci],W3[ci]); }
	  scl_vec_mult(hati,rowX,rowX); 
	  vec_subtr(W3[ci],rowX,W3[ci]); 
	  }

	  if (*retur==1) { 
		if (stat==0){
		  cumAit[i*(*Ntimes)+l+1]=cumAit[i*(*Ntimes)+l+1]-hati;
		} else {
		  cumAit[i*(*Ntimes)+l]=cumAit[i*(*Ntimes)+l]+1*(i==pers)-hati;
		}
	  }
	  if (stat==1)  replace_row(W3t[ci],l,W3[ci]);   
	} /* i=1..antpers */ 
    } /* robust ==1 */


    if (stat==1) {
      if (*deltaweight==0) dtime1=1; else dtime1=dtime; 
      extract_row(X,pers,xi); 
      vec_star(xi,dA,rowX); 
      ahati=vec_sum(rowX); 
      loglike[0]=loglike[0]-ahati/dtime1; 

      for (k=1;k<=*px;k++) cu[k*(*Ntimes)+l]=cu[k*(*Ntimes)+l-1]+cu[k*(*Ntimes)+l];

      MxA(C[ls[l]],Vargam,tmpM4); 
      MAt(tmpM4,C[ls[l]],VarKorG);
      MxA(M1M2[l],ICGam,tmpM4); 
      MAt(C[ls[l]],tmpM4,GCdM1M2);

      for (k=1;k<=*px;k++) {
	vcu[k*(*Ntimes)+l]= vcu[k*(*Ntimes)+l]+ME(VarKorG,k-1,k-1)-2.0*ME(GCdM1M2,k-1,k-1); 
      } 
    }
  } /* s=1 ..Ntimes */ 

  vec_star(IZHdN,gam,rowZ);    

  loglike[0]=loglike[0]-vec_sum(rowZ); 

  tau=time;
  for (s=1;s<*Ntimes;s++) {

    Mv(C[ls[s]],gam,korG); 
    for (k=1;k<=*px;k++) {
      cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s]-VE(korG,k-1); 
    }

    /* ROBUST VARIANCES   */ 
    if (*robust==1) {
      vec_zeros(VdB); 
      mat_zeros(Vcov); 

      for (j=0;j<*antclust;j++) {
	Mv(ICGam,W2[j],tmpv2); 
	Mv(C[ls[s]],tmpv2,rowX); 
	extract_row(W3t[j],s,tmpv1); 
	vec_subtr(tmpv1,rowX,difX); 

	replace_row(W4t[j],s,difX); 
	vec_star(difX,difX,tmpv1); 
	vec_add(tmpv1,VdB,VdB);

	if (*resample==1) {
	  if (s==1){ 
	    for (k=0;k<*pg;k++){ 
	      gammaiid[k*(*antclust)+j]=VE(tmpv2,k);
	    }
	  }
	  for (k=0;k<*px;k++) {
	    l=j*(*px)+k; 
	    Biid[l*(*Ntimes)+s]=VE(difX,k);
	  } 
	}

	if (*covariance==1) {
	  for (k=0;k<*px;k++) {
	    for (c=0;c<*px;c++){
	      ME(Vcov,k,c)=ME(Vcov,k,c)+VE(difX,k)*VE(difX,c);
	    }
	  }
	}

	if (s==1) { 
	  for (c=0;c<*pg;c++){ 
	    for (k=0;k<*pg;k++){
	      ME(RobVargam,c,k)=ME(RobVargam,c,k)+VE(W2[j],c)*VE(W2[j],k);
	    }
	  }
	}
      } /* j =1 ..Antclust */ 
    } /* robust==1 */

    for (k=1;k<*px+1;k++) { 
      Robvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      if (*covariance==1) {
	for (j=0;j<*px;j++)  {
	  l=(k-1)*(*px)+j; 
	  covs[l*(*Ntimes)+s]=ME(Vcov,k-1,j); 
	}
      }
    }
  }  /*  s=1 ..Ntimes */ 

  if (*robust==1) {
    MxA(RobVargam,ICGam,tmpM2); 
    MxA(ICGam,tmpM2,RobVargam); 
  }

  for (j=0;j<*pg;j++) {
    gamma[j]=VE(gam,j); 
    intZHdN[j]=VE(IZHdN,j); 
    for (k=0;k<*pg;k++) {
      Vgamma[k*(*pg)+j]=ME(Vargam,j,k);
      RobVgamma[k*(*pg)+j]=ME(RobVargam,j,k);
      intZHZ[k*(*pg)+j]=ME(CGam,j,k); 
    }
  }

  cu[0]=times[0]; vcu[0]=times[0];

  if (*sim==1) {
    comptest(times,Ntimes,px,cu,Robvcu,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antclust); 
  }

  free_mat(Vcov); free_mat(X); free_mat(WX); free_mat(Z); free_mat(WZ); free_mat(A); 
  free_mat(AI); free_mat(ZWZ); free_mat(VarKorG);free_mat(Vargam);free_mat(dVargam); 
  free_mat(ICGam); free_mat(CGam); free_mat(dCGam); free_mat(M1M2t); 
  free_mat(dM1M2); free_mat(Ct); free_mat(AIXW); free_mat(ZH); 
  free_mat(dC); free_mat(XWZ); free_mat(XWZAI); free_mat(GCdM1M2); 
  free_mat(tmpM2); free_mat(tmpM3); free_mat(tmpM4);
  free_vec(dA); free_vec(tmpv1); free_vec(tmpv2); free_vec(difX); 
  free_vec(korG); free_vec(rowX); free_vec(AIXWdN); free_vec(bhatt); free_vec(zi); 
  free_vec(rowZ); free_vec(gam); free_vec(dgam); free_vec(ZHdN); free_vec(IZHdN); 
  free_vec(dN); free_vec(pbhat); free_vec(pghat); free_vec(plamt); 
  for (j=0;j<*Ntimes;j++) { free_mat(M1M2[j]); }
  for (j=0;j<*Nalltimes;j++) { free_mat(Acorb[j]); free_mat(C[j]); }
  if (*robust==1) { 
    for (j=0;j<*antpers;j++) { free_mat(AIxit[j]); }
    for (j=0;j<*antclust;j++) { 
      free_mat(W3t[j]); free_mat(W4t[j]); free_vec(W2[j]); free_vec(W3[j]); } 
  }
  free(vcudif); free(times); free(cluster); free(ls); 
}
