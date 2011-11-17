//#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
void OStimecox(times,Ntimes,designX,nx,p,antpers,start,stop,nb,bhat,cu,vcu,it,b,degree,id,status,sim,antsim,cumAit,test,rani,testOBS,Ut,simUt,robvcu,retur,weighted,cumAiid,robust,covariance,covs)
double *designX,*times,*start,*stop,*cu,*vcu,*bhat,*b,
*cumAit,*test,*testOBS,*Ut,*simUt,*robvcu,*cumAiid,*covs;
int *nx,*p,*antpers,*Ntimes,*nb,*it,*degree,*id,*status,*sim,*antsim,*rani,*retur,*weighted,*robust,*covariance;
{
  matrix *Vcov,*ldesignX,*A,*AI,*cdesignX;
  matrix *cumAt[*antpers];
  vector *difX,*rowX,*xi,*vtmp,*vtmp1,*diag,*dB,*dN,*dM,*VdB,*AIXdM,*AIXdN,*AIXlamt,*ta,*bhatt,*pbhat,*plamt,*lrisk,*score;
  vector *cumhatA[*antpers],*cumA[*antpers],*cum,*dN1;
  int i,j,k,l,s,c,count,pers=0,itt,
        *coef=calloc(1,sizeof(int)),
        *ps=calloc(1,sizeof(int)),*imin=calloc(1,sizeof(int));;
  double ahati,time,dummy,dtime;
  double *vcudif=calloc((*Ntimes)*(*p+1),sizeof(double)),
	 *sbhat=calloc((*Ntimes)*(*p+1),sizeof(double)),
         *sscore=calloc((*Ntimes)*(*p+1),sizeof(double)); 

  malloc_mats(*antpers,*p,&ldesignX,&cdesignX,NULL);
  malloc_mats(*p,*p,&Vcov,&A,&AI,NULL);
  malloc_vecs(*p,&vtmp,&VdB,&rowX,&difX,&cum,&xi,&vtmp1,&diag,&dB,&VdB,&AIXdM,&AIXdN,&AIXlamt,&bhatt,&score,NULL);
  malloc_vecs(*antpers,&dN1,&dN,&pbhat,&plamt,&lrisk,&dM,NULL);
  malloc_vecs(*nb,&ta,NULL);

  if (*robust==1) 
    for (i=0;i<*antpers;i++) {malloc_vec(*p,cumhatA[i]); malloc_vec(*p,cumA[i]); 
      malloc_mat(*Ntimes,*p,cumAt[i]);}

  coef[0]=1; ps[0]=*p+1; 
  for (itt=0;itt<*it;itt++)
    {
      vec_zeros(score); 

      for (s=1;s<*Ntimes;s++)
	{
	  time=times[s]; vec_zeros(dN);dtime=time-times[s-1]; 
	  mat_zeros(ldesignX); vec_zeros(lrisk); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++)
	    {
	      if ((start[c]<time) && (stop[c]>=time)) {
		VE(lrisk,id[c])=1; 
		for(j=0;j<*p;j++) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		if (time==stop[c] && status[c]==1) {VE(dN,id[c])=1; pers=id[c];}
		count=count+1; } 
	    }

	  for(j=0;j<*nb;j++) VE(ta,j)=fabs(bhat[j]-time); 
	  dummy=vec_min(ta,imin);
	  for(j=1;j<=*p;j++) VE(bhatt,j-1)=bhat[j*(*nb)+(*imin)];
	  Mv(ldesignX,bhatt,pbhat);

	  for (j=0;j<*antpers;j++) { 
	    VE(plamt,j)=VE(lrisk,j)*exp(VE(pbhat,j)); 
	    scl_vec_mult(VE(plamt,j),extract_row(ldesignX,j,dB),dB);
	    replace_row(cdesignX,j,dB); /* sampling corrected design */ }
				
	  MtA(cdesignX,ldesignX,A); invert(A,AI); 
	  if (ME(AI,0,0)==0)
	    Rprintf(" X'X not invertible at time %lf \n",time);

	  extract_row(ldesignX,pers,xi); Mv(AI,xi,AIXdN); 
	  vM(ldesignX,plamt,vtmp1); Mv(AI,vtmp1,AIXlamt); 
	  scl_vec_mult(dtime,AIXlamt,AIXlamt); vec_subtr(AIXdN,AIXlamt,AIXdM); 

	  extract_row(ldesignX,pers,xi); scl_vec_mult(dtime,vtmp1,vtmp1);  
	  vec_subtr(xi,vtmp1,vtmp1); vec_add(vtmp1,score,score); 

	  for (k=1;k<=*p;k++) sscore[k*(*Ntimes)+s]=VE(score,k-1); 
	   sscore[s]=time; 

	  cu[s]=time; vcu[s]=time; robvcu[s]=time; sbhat[s]=time;
	  if (itt==0) 
	    for (k=1;k<=*p;k++) sbhat[k*(*Ntimes)+s-1]=bhat[k*(*Ntimes)+s-1];

	  for (k=1;k<=*p;k++) {
	    cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+(s-1)]+dtime*VE(bhatt,k-1)+VE(AIXdM,k-1);
	    sbhat[k*(*Ntimes)+s-1]=bhat[k*(*Ntimes)+s-1]+VE(AIXdM,k-1); 
	    /*
	     ssscore[s]=time;
	      Rprintf(" %lf %lf ",sbhat[k*(*Ntimes)+s-1]-bhat[k*(*Ntimes)+s-1],
	      ssscore[k*(*Ntimes)+s-1]); Rprintf("\n"); 
	    */
	    vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+(s-1)]+dtime*ME(AI,k-1,k-1);}

	  if (itt==(*it-1)) {
	    if (*robust==1) 
	      {
		vec_zeros(VdB);
		for (i=0;i<*antpers;i++) { 
		  extract_row(ldesignX,i,xi); Mv(AI,xi,rowX); 
		  ahati=VE(plamt,i)*dtime;
		  vec_star(xi,AIXdM,vtmp1); dummy=vec_sum(vtmp1); 
		  if (i==pers) vec_add(rowX,cumhatA[i],cumhatA[i]);
		  scl_vec_mult(ahati,rowX,rowX);
		  vec_add(rowX,cumA[i],cumA[i]);
		  vec_subtr(cumhatA[i],cumA[i],difX); replace_row(cumAt[i],s,difX);
		  vec_star(difX,difX,vtmp); vec_add(vtmp,VdB,VdB);
		  if (*covariance==1) {
		    for (k=0;k<*p;k++) for (j=0;j<*p;j++)
		      ME(Vcov,k,j)=ME(Vcov,k,j)+VE(difX,k)*VE(difX,j);}

		  if (*retur==1) {
		    cumAit[i*(*Ntimes)+s]=1*(i==pers)-ahati;
		    cumAiid[i*(*Ntimes)+s]=cumAit[i*(*Ntimes)+s]-VE(plamt,i)*dummy;
		    /*
		      VE(dN,i)=cumAiid[i*(*Ntimes)+s];
		      VE(dM,i)=cumAit[i*(*Ntimes)+s];
		    */
		    cumAit[i*(*Ntimes)+s]=cumAiid[i*(*Ntimes)+s]; 
		  }
		}
		for (k=1;k<*p+1;k++){robvcu[k*(*Ntimes)+s]=VE(VdB,k-1); 
		  if (*covariance==1) {
		    for (j=0;j<*p;j++)  {
		      l=(k-1)*(*p)+j; covs[l*(*Ntimes)+s]=ME(Vcov,k-1,j);
		    }}}

	      } /*robust=1 */
	  }
  
	} /* s =1,...,*Ntimes */ 
      /* print_vec(score); */

      smoothB(cu,Ntimes,ps,bhat,nb,b,degree,coef); 

    } /* itterations løkke */ 
  cu[0]=times[0]; vcu[0]=times[0];

  if (*sim==1) {
    comptest(times,Ntimes,p,cu,robvcu,vcudif,antsim,test,testOBS,Ut,simUt,cumAt,weighted,antpers);
  }

  if (*robust==1)
    for (i=0;i<*antpers;i++) {free_mat(cumAt[i]);free_vec(cumA[i]);free_vec(cumhatA[i]);}

  free_mats(&Vcov,&ldesignX,&A,&AI,&cdesignX,NULL);
  free_vecs(&dM,&difX,&rowX,&xi,&vtmp,&vtmp1,&diag,&dB,&dN,&VdB,&AIXdN,&AIXlamt,&ta,&bhatt,&pbhat,&plamt,&lrisk,&cum,NULL);
  free(vcudif); free(sbhat); free(sscore); 
  free(coef); free(ps); free(imin);
}

void OSsemicox(times,Ntimes,designX,nx,px,designG,ng,pg,
antpers,start,stop,nb,bhat,cu,vcu,gamma,Vgamma,b,degree,it,
RobVgamma,robvcu,sim,antsim,retur,cumAit,test,rani,testOBS,status,Ut,simUt,id,weighted,robust,covariance,covs)
double *designX,*times,*start,*stop,*cu,*vcu,*bhat,*b,*designG,
*gamma,*Vgamma,*RobVgamma,*cumAit,*test,*testOBS,*Ut,*simUt,*robvcu,*covs; 
int
*nx,*px,*antpers,*Ntimes,*nb,*ng,*pg,*it,*degree,*id,*status,*sim,*antsim,*retur,*rani,*weighted,*robust,*covariance;
{
  matrix *ldesignX,*A,*AI,*cdesignX,*ldesignG,*cdesignG;
  matrix *S,*dCGam,*CGam,*ICGam,*VarKorG,*dC,*XZ,*ZZ,*ZZI,*XZAI; 
  matrix *Ct,*C[*Ntimes],*Acorb[*Ntimes],*Vcov; 
  matrix *tmpM2,*tmpM3,*tmpM4; 
  matrix *Vargam,*dVargam,*M1M2[*Ntimes];
  matrix *dM1M2,*M1M2t,*RobVargam;
  matrix *W3t[*antpers],*W4t[*antpers],*AIxit[*antpers];
  vector *W2[*antpers],*W3[*antpers];
  vector *diag,*dB,*dN,*VdB,*AIXdN,*AIXlamt,*ta,*bhatt,*pbhat,*plamt;
  vector *korG,*pghat,*rowG,*gam,*dgam,*ZGdN,*IZGdN,*ZGlamt,*IZGlamt;
  vector *xi,*rowX,*rowZ,*difX,*zi,*z1,*tmpv1,*tmpv2,*lrisk;
  int itt,i,j,k,l,s,c,count,pers=0,pmax,
      *coef=calloc(1,sizeof(int)),*ps=calloc(1,sizeof(int)),
      *imin=calloc(1,sizeof(int));;
  double time,dtime,hati;
  double *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double));
  int *ipers=calloc(*Ntimes,sizeof(int));

  if (*robust==1) 
    for (j=0;j<*antpers;j++) { malloc_mat(*Ntimes,*px,W3t[j]);
      malloc_mat(*Ntimes,*px,W4t[j]); malloc_vec(*pg,W2[j]); malloc_vec(*px,W3[j]);
      malloc_mat(*Ntimes,*px,AIxit[j]); }

  malloc_mats(*antpers,*px,&ldesignX,&cdesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,&cdesignG,NULL); 
  malloc_mats(*px,*px,&VarKorG,&Vcov,&A,&AI,NULL);
  malloc_mats(*pg,*pg,&dVargam,&Vargam,&RobVargam,&tmpM2,&ZZ,&ICGam,&CGam,&dCGam,&S,&ZZI,NULL); 
  malloc_mats(*px,*pg,&XZAI,&tmpM3,&Ct,&dC,&XZ,&dM1M2,&M1M2t,NULL);
  malloc_mat(*px,*pg,tmpM4); 
  for (j=0;j<*Ntimes;j++) { malloc_mat(*pg,*px,Acorb[j]); 
    malloc_mat(*px,*pg,C[j]); malloc_mat(*px,*pg,M1M2[j]);}

  malloc_vecs(*px,&xi,&rowX,&difX,&tmpv1,&korG,&diag,&dB,&VdB,&AIXdN,&AIXlamt,&bhatt,NULL);
  malloc_vecs(*pg,&zi,&rowZ,&tmpv2,&zi,&z1,&rowG,&gam,&dgam,&ZGdN,&IZGdN,&ZGlamt,&IZGlamt,NULL);
  malloc_vecs(*antpers,&lrisk,&dN,&pbhat,&pghat,&plamt,NULL);
  malloc_vecs(*nb,&ta,NULL);

  coef[0]=1; ps[0]=*px+1; 
  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 

  for (itt=0;itt<*it;itt++)
    {
      mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZGdN); vec_zeros(IZGlamt); 
      /*
	Rprintf("Itteration, it loop, Number Jumps %ld %ld %ld \n",*it,itt,*Ntimes); 
      */

      for (s=1;s<*Ntimes;s++)
	{
	  time=times[s]; dtime=time-times[s-1]; 
	  mat_zeros(ldesignX);mat_zeros(ldesignG);vec_zeros(lrisk);vec_zeros(dN); 

	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++)
	    {
	      if ((start[c]<time) && (stop[c]>=time)) {
		VE(lrisk,id[c])=1; 
		for(j=0;j<pmax;j++)  {
		  if (j<*px) ME(ldesignX,id[c],j)= designX[j*(*nx)+c];
		  if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; }
		if (time==stop[c] && status[c]==1) {VE(dN,id[c])=1; pers=id[c];}
		count=count+1; } 
	    }
	  ipers[s]=pers;


	  for(j=0;j<*nb;j++) VE(ta,j)=fabs(bhat[j]-time);
	  for(j=1;j<=*px;j++) VE(bhatt,j-1)=bhat[j*(*nb)+(*imin)];
	  Mv(ldesignX,bhatt,pbhat); Mv(ldesignG,gam,pghat);

	  for (j=0;j<*antpers;j++) { 
	    VE(plamt,j)=VE(lrisk,j)*exp(VE(pbhat,j)+VE(pghat,j)); 
	    scl_vec_mult(VE(plamt,j),extract_row(ldesignX,j,dB),dB);
	    replace_row(cdesignX,j,dB); /* sampling corrected design */ 
	    scl_vec_mult(VE(plamt,j),extract_row(ldesignG,j,rowG),rowG);
	    replace_row(cdesignG,j,rowG); /* sampling corrected design */ 
	  }
				
	  MtA(cdesignX,ldesignX,A); invert(A,AI); 
	  if (ME(AI,0,0)==0)
	    Rprintf(" X'X not invertible at time %lf \n",time);

	  extract_row(ldesignX,pers,xi); Mv(AI,xi,AIXdN);

	  vM(ldesignX,plamt,xi); Mv(AI,xi,AIXlamt); 

	  MtA(ldesignG,cdesignG,ZZ); MtA(ldesignX,cdesignG,XZ);
	  MxA(AI,XZ,XZAI); MtA(XZAI,XZ,tmpM2); 
	  mat_subtr(ZZ,tmpM2,dCGam); 
	  scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

	  extract_row(ldesignG,pers,zi); 
	  vM(XZ,AIXdN,tmpv2); vec_subtr(zi,tmpv2,ZGdN); vec_add(ZGdN,IZGdN,IZGdN); 
	  Acorb[s]=mat_transp(XZAI,Acorb[s]); 

	  scl_mat_mult(dtime,XZAI,tmpM4);mat_add(tmpM4,Ct,Ct);C[s]=mat_copy(Ct,C[s]); 

	  vM(ldesignG,plamt,z1); vM(XZ,AIXlamt,tmpv2);
	  vec_subtr(z1,tmpv2,ZGlamt);
	  vec_add(scl_vec_mult(dtime,ZGlamt,tmpv2),IZGlamt,IZGlamt); 

	  for (k=1;k<=*px;k++) {
	    cu[k*(*Ntimes)+s]=
	      cu[k*(*Ntimes)+s-1]+
	      (dtime*VE(bhatt,k-1))+(VE(AIXdN,k-1)-dtime*VE(AIXlamt,k-1)); 
	    vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s-1]+dtime*ME(AI,k-1,k-1);}

	  if (itt==(*it-1)) { 
	    for (k=0;k<*pg;k++) {
	      for (j=0;j<*pg;j++) ME(dVargam,k,j)=VE(ZGdN,j)*VE(ZGdN,k);
	      for (j=0;j<*px;j++) ME(dM1M2,j,k)=VE(ZGdN,k)*VE(AIXdN,j); }

	    mat_add(dVargam,Vargam,Vargam);
	    mat_add(dM1M2,M1M2t,M1M2t); M1M2[s]=mat_copy(M1M2t,M1M2[s]);  

	    if (*robust==1) 
	      for (i=0;i<*antpers;i++) {extract_row(ldesignX,i,xi);
		Mv(AI,xi,rowX);replace_row(AIxit[i],s,rowX);} } 
	} /* s=1,...Ntimes */

      invert(CGam,ICGam); vec_subtr(IZGdN,IZGlamt,tmpv2); 
      Mv(ICGam,tmpv2,dgam); vec_add(gam,dgam,gam); 

      for (s=1;s<*Ntimes;s++) 
	{
	  time=times[s]; vec_zeros(dN);dtime=time-times[s-1]; 

	  Mv(C[s],dgam,korG);

	  cu[s]=times[s]; vcu[s]=times[s]; robvcu[s]=times[s]; 
	  for (k=1;k<=*px;k++) cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s]-VE(korG,k-1);


	} /* s=1,...Ntimes */
      smoothB(cu,Ntimes,ps,bhat,nb,b,degree,coef); 
    } /*itt løkke */ 

  invert(Vargam,dVargam); 

  for (s=1;s<*Ntimes;s++) {
    if (*robust==1)
      {
	time=times[s]; vec_zeros(dN);dtime=time-times[s-1]; 
	mat_zeros(ldesignX); mat_zeros(ldesignG); vec_zeros(lrisk);
	for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) 
	  {
	    if ((start[c]<time) && (stop[c]>=time)) {
	      VE(lrisk,id[c])=1; 
	      for(j=0;j<pmax;j++) {
		if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; }
	      if (time==stop[c] && status[c]==1) pers=id[c];
	      count=count+1;} 
	  }

	for(j=0;j<*nb;j++) VE(ta,j)=fabs(bhat[j]-time);
	for(j=1;j<=*px;j++) VE(bhatt,j-1)=bhat[j*(*nb)+(*imin)];
	Mv(ldesignX,bhatt,pbhat); Mv(ldesignG,gam,pghat);

	for (i=0;i<*antpers;i++) {
	  VE(plamt,i)=VE(lrisk,i)*exp(VE(pbhat,i)+VE(pghat,i)); 
	  extract_row(ldesignX,i,xi); extract_row(ldesignG,i,zi); 

	  hati=dtime*VE(plamt,i);  

	  Mv(Acorb[s],xi,tmpv2); vec_subtr(zi,tmpv2,tmpv2); 

	  if (i==ipers[s]) vec_add(tmpv2,W2[i],W2[i]);
	  scl_vec_mult(hati,tmpv2,rowZ); vec_subtr(W2[i],rowZ,W2[i]); 

	  extract_row(AIxit[i],s,rowX); 
	  if (i==ipers[s]) vec_add(rowX,W3[i],W3[i]); 

	  scl_vec_mult(hati,rowX,rowX); 
	  vec_subtr(W3[i],rowX,W3[i]); 
	  replace_row(W3t[i],s,W3[i]);  

	  if (*retur==1)
	    cumAit[i*(*Ntimes)+s]=1*(i==ipers[s])-hati;
	} /* i=1..antpers */ 
      }

    MxA(C[s],dVargam,tmpM4); MAt(tmpM4,C[s],VarKorG);
    for (k=1;k<=*px;k++) 
      vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s]+ME(VarKorG,k-1,k-1);   
    /* 
       MAt(ICGam, M1M2[s], tmpM4); MxA(C[s],tmpM4,GCdM1M2);
       for (k=1;k<=*px;k++) {vcu[k*(*Ntimes)+s]= 
       vcu[k*(*Ntimes)+s]+ME(VarKorG,k-1,k-1)-2*ME(GCdM1M2,k-1,k-1);} 
    */
  } /* s=1 ..Ntimes */ 


  /* ROBUST VARIANCES   */ 
  if (*robust==1) {
    for (s=1;s<*Ntimes;s++) {
      vec_zeros(VdB); 
      for (i=0;i<*antpers;i++) {
	Mv(ICGam,W2[i],tmpv2); vM(Acorb[s],tmpv2,rowX); 
	extract_row(W3t[i],s,tmpv1); vec_subtr(tmpv1,rowX,difX); 
	replace_row(W4t[i],s,difX); 
	vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	if (*covariance==1) {
	  for (k=0;k<*px;k++) for (j=0;j<*px;j++)
	    ME(Vcov,k,j)=ME(Vcov,k,j)+VE(difX,k)*VE(difX,j);}


	if (s==1) { for (j=0;j<*pg;j++) for (k=0;k<*pg;k++) 
		      ME(RobVargam,j,k)=ME(RobVargam,j,k)+VE(W2[i],j)*VE(W2[i],k);}
      } /* i =1 ..Antpers */ 
      for (k=1;k<*px+1;k++) {robvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
	if (*covariance==1) {
	  for (j=0;j<*px;j++)  {
	    l=(k-1)*(*px)+j; covs[l*(*Ntimes)+s]=ME(Vcov,k-1,j); }}}
    }  /*  s=1 ..Ntimes */ 

    MxA(RobVargam,ICGam,tmpM2); MxA(ICGam,tmpM2,RobVargam); 
  }

  for (j=0;j<*pg;j++) {gamma[j]=VE(gam,j);
    for (k=0;k<*pg;k++) {Vgamma[k*(*pg)+j]=ME(dVargam,j,k);
      RobVgamma[k*(*pg)+j]=ME(RobVargam,j,k);}}
  cu[0]=times[0]; vcu[0]=times[0];

  if (*sim==1) {
    comptest(times,Ntimes,px,cu,robvcu,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antpers);
  }

  free_mats(&ldesignX,&A,&AI,&cdesignX,&ldesignG,&cdesignG,
	      &S,&dCGam,&CGam,&ICGam,&VarKorG,&dC,&XZ,&ZZ,&ZZI,&XZAI, 
	      &Ct,&tmpM2,&tmpM3,&tmpM4,&Vargam,&dVargam,
	      &dM1M2,&M1M2t,&RobVargam,NULL); 

  free_vecs(&diag,&dB,&dN,&VdB,&AIXdN,&AIXlamt,&ta,&bhatt,&pbhat,
	      &plamt,&korG,&pghat,&rowG,&gam,&dgam,&ZGdN,&IZGdN,&ZGlamt,&IZGlamt,
	      &xi,&rowX,&rowZ,&difX,&zi,&z1,&tmpv1,&tmpv2,&lrisk,NULL); 

  for (j=0;j<*Ntimes;j++) {free_mat(Acorb[j]);free_mat(C[j]);free_mat(M1M2[j]);}
  if (*robust==1) 
    for (j=0;j<*antpers;j++) {free_mat(W3t[j]); free_mat(W4t[j]);
      free_mat(AIxit[j]); free_vec(W2[j]); free_vec(W3[j]); }
      free(vcudif); free(ipers); 
      free(coef); free(ps); free(imin);
}
