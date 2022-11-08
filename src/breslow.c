//#include <stdio.h>
#include <math.h>
#include <R.h>
#include "matrix.h"
                 
void OSbreslow(double *times,int *Ntimes,double *designX,int *nx,int *p,int *antpers,double *start,double *stop,int *nb,double *bhat,double *cu,double *vcu,int *it,double *b,int *degree,double *schoen,int *sim,int *antsim,double *test,int *rani,double *testOBS,double *rvcu,double *cumlam,int *nullresid,int *status,int *id,int *sim2,double *Ut,double *simUt,int *weighted,int *robust)
{
  matrix *ldesignX,*A,*AI,*AIX,*cdesignX,*XmavX,*cXmavX,*Aav;
  vector *diag,*dB,*dN,*VdB,*AIXdN,*AIXlamt,*ta,*bhatt,*pbhat,*plamt,*avx,*lrisk;
  vector *ssrow2,*ssrow,*vtmp,*xi,*rowX,*difX,*Btau,*score; 
//  matrix *Delta2,*Delta,*tmpM1,*tmpM2,*varBL;
  int supsup=0,itt,i,j,k,s,c,count,pers=0,
      *imin=calloc(1,sizeof(int)), *coef=calloc(1,sizeof(int)),*ps=calloc(1,sizeof(int));
  double time2,rr,time=0,time1,dummy,dtime,S0,lam0t,sdBt,tau,random;
  double *Basei=calloc(*antpers,sizeof(double)),rvarbase, *vcudif=calloc((*Ntimes)*(*p+2),sizeof(double));
//  double norm_rand(); void GetRNGstate(),PutRNGstate();
  matrix *Delta,*Delta2,*tmpM1,*tmpM2; 
  matrix *varBL,*cumBL[*antpers],*cumB[*antpers],*BLsubbetaLam[*antpers];
  vector *cumi[*antpers],*cumBLi[*antpers],*Base[*antpers]; 

//  if (*sim==1) { }; 
    malloc_mat(*Ntimes,*p,Delta); 
    malloc_mat(*Ntimes,*p,Delta2); 
    malloc_mat(*Ntimes,*p,tmpM1); 
    malloc_mat(*Ntimes,*p,tmpM2);

//  if (*robust==1) { }
    malloc_mat(*Ntimes,*p,varBL); 
    for (j=0;j<*antpers;j++) {
      malloc_mat(*Ntimes,(*p)+1,cumB[j]); 
      malloc_vec(*p,cumi[j]); 
      malloc_mat(*Ntimes,*p,BLsubbetaLam[j]); 
      malloc_mat(*Ntimes,*p,cumBL[j]); 
      malloc_vec(*p,cumBLi[j]); 
      malloc_vec(*Ntimes,Base[j]); 
      Basei[j]=0.0; 
    }

  malloc_mat(*antpers,*p,ldesignX);
  malloc_mat(*antpers,*p,cdesignX);
  malloc_mat(*antpers,*p,XmavX);
  malloc_mat(*antpers,*p,cXmavX);
  malloc_mat(*p,*antpers,AIX);
  malloc_mat(*p,*p,A);
  malloc_mat(*p,*p,AI);
  malloc_mat(*p,*p,Aav);
  malloc_vec(*p,score);
  malloc_vec(*p,ssrow2);
  malloc_vec(*p,ssrow);
  malloc_vec(*p,Btau);
  malloc_vec(*p,vtmp);
  malloc_vec(*p,difX);
  malloc_vec(*p,xi);
  malloc_vec(*p,rowX);
  malloc_vec(*p,avx);
  malloc_vec(*p,diag);
  malloc_vec(*p,dB);
  malloc_vec(*p,VdB);
  malloc_vec(*p,AIXdN);
  malloc_vec(*p,AIXlamt);
  malloc_vec(*p,bhatt);
  malloc_vec(*p,dN);
  malloc_vec(*p,pbhat);
  malloc_vec(*p,plamt);
  malloc_vec(*p,lrisk);
  malloc_vec(*nb,ta);

  coef[0]=1; ps[0]=(*p)+2; tau=times[*Ntimes-1]; 

  for (itt=0;itt<*it;itt++){
    vec_zeros(score); 
    for (s=1;s<*Ntimes;s++){
      time=times[s]; 
      dtime=time-times[s-1]; 
      vec_zeros(lrisk); 
      mat_zeros(ldesignX); 

      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++){
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<*p;j++) {
	    ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	  }
	  VE(lrisk,id[c])=1;
	  if (time==stop[c] && status[c]==1){ 
	    pers=id[c];
	  }
	  count=count+1; 
	} 
      }

      for(j=0;j<*nb;j++){ 
	VE(ta,j)=fabs(bhat[j]-time);
      }
      dummy=vec_min(ta,imin); 
      lam0t=bhat[1*(*nb)+(*imin)]; 
      for(j=2;j<=(*p)+1;j++){ 
	VE(bhatt,j-2)=bhat[j*(*nb)+(*imin)];
      }

      Mv(ldesignX,bhatt,pbhat);

      for (j=0;j<*antpers;j++) { 
	VE(plamt,j)=VE(lrisk,j)*exp(VE(pbhat,j)); 
	scl_vec_mult(VE(plamt,j),extract_row(ldesignX,j,dB),dB);
	replace_row(cdesignX,j,dB); /* sampling corrected design */ 
      }

      S0=vec_sum(plamt); 
      vM(ldesignX,plamt,avx); 
      scl_vec_mult(1/S0,avx,avx); 

      for (j=0;j<*p;j++){
	for (i=0;i<*p;i++) {
	  ME(Aav,j,i)=VE(avx,i)*VE(avx,j)*S0;
	}
      }

      MtA(cdesignX,ldesignX,A); 
      mat_subtr(A,Aav,A); 
      invert(A,AI); 
      extract_row(ldesignX,pers,AIXdN); 
      vec_subtr(AIXdN,avx,AIXdN); 
      Mv(AI,AIXdN,dB);    

      vec_add(dB,score,score); 

      schoen[s]=time; cu[s]=time; vcu[s]=time; rvcu[s]=time; 
      cu[1*(*Ntimes)+s]=cu[1*(*Ntimes)+s-1]+(1/S0); vcu[1*(*Ntimes)+s]=0; 

      for (k=2;k<=(*p)+1;k++){
	cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+dtime*VE(bhatt,k-2)+VE(dB,k-2)/lam0t;   
	cumlam[k*(*Ntimes)+s]=cumlam[k*(*Ntimes)+s-1]+(dtime*VE(bhatt,k-2)*lam0t)+VE(dB,k-2);  
	if (itt==(*it-1)) {
	  schoen[(k-1)*(*Ntimes)+s]=VE(dB,k-2)*S0;
	  vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s-1]+(dtime/lam0t)*ME(AI,k-2,k-2); 
	}
      }
      /* Rprintf(" \n");   */

      if (itt==(*it-1)) {
	if (*robust==1) {
	  vec_zeros(VdB); 
	  rvarbase=0; 
	  for (j=0;j<*antpers;j++) { 
	    extract_row(ldesignX,j,rowX); 
	    vec_subtr(rowX,avx,vtmp); 
	    Mv(AI,vtmp,xi); 

	    if (*nullresid>=0) {
	      k=*nullresid;  
	      rr=VE(pbhat,j)+VE(rowX,k)*(cu[(k+2)*(*Ntimes)+(*Ntimes-1)]/tau-VE(bhatt,k)); 
	      rr=VE(lrisk,j)*exp(rr);
	    } else { 
	      rr=VE(plamt,j); 
	    }

	    scl_vec_mult(rr/S0,xi,rowX);  
	    vec_subtr(cumBLi[j],rowX,cumBLi[j]); 
	    if (j==pers){ 
	      vec_add(xi,cumBLi[j],cumBLi[j]); 
	    }
	    replace_row(cumBL[j],s,cumBLi[j]);  /* BLAM(t) sum iid */

	    vec_star(avx,xi,rowX); 
	    dummy=vec_sum(rowX); 
	    dummy=(1/S0)-dummy; 
	    Basei[j]=Basei[j]-dummy*rr/S0;  
	    if (j==pers){ 
	      Basei[j]=dummy+Basei[j]; 
	    }
	    VE(Base[j],s)=Basei[j];  /* Baseline sum iid */
	    rvarbase=rvarbase+Basei[j]*Basei[j]; 
	    ME(cumB[j],s,0)  =Basei[j];

	    scl_vec_mult(1/lam0t,xi,xi); 
	    scl_vec_mult(rr/S0,xi,rowX);  
 
	    vec_subtr(cumi[j],rowX,cumi[j]); 
	    if (j==pers) {
	      vec_add(xi,cumi[j],cumi[j]); 
	    }
	    /* set_row(cumB[j],s,cumi[j]); */   /* B(t) sum iid */
	    for (k=1;k<(*p)+1;k++) {
	      ME(cumB[j],s,k)=VE(cumi[j],k-1); 
	    }

	    vec_star(cumi[j],cumi[j],difX); 
	    vec_add(difX,VdB,VdB); 
	  } 
	  rvcu[1*(*Ntimes)+s]=rvarbase; 
	  for (k=2;k<(*p)+2;k++) {
	    rvcu[k*(*Ntimes)+s]=VE(VdB,k-2); 
	  }
	}
      } 
    } /* s */

    smoothB(cu,Ntimes,ps,bhat,nb,b,degree,coef);

  } /* itterations lokke */ 
  for (i=2;i<(*p)+2;i++) {
    VE(Btau,i-2)=cu[i*(*Ntimes)+(*Ntimes-1)];
  }

  cu[0]=times[0]; 
  vcu[0]=times[0]; 
  tau=time; 
  rvcu[0]=times[0]; 


  /* Beregning af iid bidrag til BLam(t) - beta Lam(t) */
  if (*robust==1){
      for (s=1;s<*Ntimes;s++) {
	vec_zeros(VdB); 
	for (j=0;j<*antpers;j++) { 
	  scl_vec_mult(1/tau,cumi[j],rowX);
	  scl_vec_mult(cu[1*(*Ntimes)+s],rowX,difX);
	  scl_vec_mult(VE(Base[j],s)/tau,Btau,xi); 
	  extract_row(cumBL[j],s,vtmp); 
	  vec_add(vtmp,xi,xi); 
	  vec_subtr(xi,difX,xi); 
	  replace_row(BLsubbetaLam[j],s,xi); 

	  vec_star(xi,xi,difX); 
	  vec_add(difX,VdB,VdB); 

	  replace_row(varBL,s,VdB); 
	}
      } /*  s */
    }

  /* korrektion for lam0t i \int beta(s) lam0t(s) ds */ 
  /*
    for (s=1;s<*Ntimes;s++) {
    time=times[s]; v_zero(dN);dtime=time-times[s-1]; 
    for(j=0;j<*nb;j++) ta->ve[j]=fabs(bhat[j]-time);
    dummy=v_min(ta,imin); lam0t=bhat[1*(*nb)+(*imin)]; 
    for(j=2;j<=(*p)+1;j++) {bhatt->ve[j-2]=bhat[j*(*nb)+(*imin)]/lam0t;
    cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+dtime*bhatt->ve[k-2];  }; } 
  */

  if (*sim==1) {
    if (*sim2!=1) {
      ps[0]=*p+1; 
      comptest(times,Ntimes,ps,cu,rvcu,vcudif,antsim,test,testOBS,Ut,simUt,
	       cumB,weighted,antpers);
      /*
	for (s=0;s<*Ntimes;s++) {
	for(j=0;j<(*p)+2;j++)  Rprintf(" %lf ",cu[j*(*Ntimes)+s]); 
	Rprintf(" \n"); }
      */
    } else {
      for (i=2;i<=*p+1;i++){ 
	VE(bhatt,i-2)=cu[i*(*Ntimes)+(*Ntimes-1)];
      }
      for (s=1;s<*Ntimes-1;s++){ /* Beregning af obs teststorrelser */ 
	time=times[s]; 
	dtime=time-times[s-1]; 
	for (i=2;i<=*p+1;i++) {
	  VE(xi,i-2)=fabs(cu[i*(*Ntimes)+s])/sqrt(rvcu[i*(*Ntimes)+s]);
	  if (VE(xi,i-2)>testOBS[i-2]){ 
	    testOBS[i-2]=VE(xi,i-2); 
	  }
	}

	scl_vec_mult(time/tau,bhatt,difX);
	for (i=2;i<=*p+1;i++){ 
	  VE(xi,i-2)=cu[i*(*Ntimes)+s];
	}
	vec_subtr(xi,difX,difX); 
	vec_star(difX,difX,ssrow); 
  	for (i=0;i<*p;i++) {c=(*p+i); 
	  VE(difX,i)=fabs(VE(difX,i));/*sqrt(rvcu[(i+2)*(*Ntimes)+s]);*/
	  if (VE(difX,i)>testOBS[c]){ testOBS[c]=VE(difX,i);}
	  c=2*(*p)+i; 
	  testOBS[c]=testOBS[c]+VE(ssrow,i)*dtime;
	} 

	scl_vec_mult(1/tau,bhatt,rowX); 
	scl_vec_mult(cu[1*(*Ntimes)+s],rowX,rowX);
	for (i=2;i<=*p+1;i++) {
	  VE(xi,i-2)=cumlam[i*(*Ntimes)+s];
	}   
	vec_subtr(xi,rowX,difX); 
	vec_star(difX,difX,ssrow); 
	for (i=0;i<*p;i++) { 
	  c=3*(*p)+i; 
	  VE(difX,i)=fabs(VE(difX,i)); 
	  if (VE(difX,i)>testOBS[c]) testOBS[c]=VE(difX,i);
	  c=4*(*p)+i; 
	  testOBS[c]=testOBS[c]+VE(ssrow,i)*dtime;
	} 
	/* sup| BLAM(t)-beta*LAM(t)| */

	/* Beregning af sup_a,s | B(a+t)-B(t) - gam t | */ 
	if (supsup==1) {
	  for (j=s+1;j<*Ntimes;j++){ 
	    time1=times[j];  
	    time2=time1-time; 
	    scl_vec_mult(time2/tau,bhatt,difX);
	    for (i=2;i<=*p+1;i++) {
	      VE(xi,i-2)=cu[i*(*Ntimes)+s]; 
	      VE(vtmp,i-2)=cu[i*(*Ntimes)+j];
	    }
	    vec_subtr(vtmp,xi,rowX); 
	    vec_subtr(rowX,difX,xi); 
	    c=5*(*p)+i; 
	    for (i=2;i<=*p+1;i++) {
	      if (fabs(VE(xi,i-2))>testOBS[c]) testOBS[c]=fabs(VE(xi,i-2)); 
	    }
	  } 
	} /* supsup==1 */ 

      } /*s=1..Ntimes Beregning af obs teststorrelser */ 


      Rprintf(" Simulations start N= %ld \n",(long int) *antsim); 
      GetRNGstate();  /* to use R random normals */

      for (k=1;k<*antsim;k++) {
	if (k%50==0)  Rprintf(" %ld Simulations \n",(long int) k); 
	mat_zeros(Delta); 
	mat_zeros(Delta2); 
	vec_zeros(vtmp); 
	for (i=0;i<*antpers;i++) {
	  random=norm_rand();
	  scl_mat_mult(random,cumB[i],tmpM1); 
	  mat_add(tmpM1,Delta,Delta); 
	  if (*sim2==1) {
	    scl_mat_mult(random,BLsubbetaLam[i],tmpM2); 
	    mat_add(tmpM2,Delta2,Delta2);
	  }

	  extract_row(cumB[i],*Ntimes-1,rowX);
	  scl_vec_mult(random,rowX,rowX); 
	  vec_add(rowX,vtmp,vtmp); 
	}

	for (s=1;s<*Ntimes;s++){
	  time=times[s]; dtime=time-times[s-1]; 
	  scl_vec_mult(time/tau,vtmp,xi); 
	  extract_row(Delta,s,rowX); 
	  vec_subtr(rowX,xi,difX); 
	  vec_star(difX,difX,ssrow); 

	  if (*sim2==1) {
	    extract_row(Delta2,s,dB); 
	    vec_star(dB,dB,ssrow2);  
	  }

	  for (i=0;i<*p;i++) {
	    VE(difX,i)=fabs(VE(difX,i)); 
	    VE(dB,i)=fabs(VE(dB,i)); 

	    sdBt=sqrt(rvcu[(i+2)*(*Ntimes)+s]+0.03);
	    VE(xi,i)=fabs(ME(Delta,s,i))/sdBt;
	    if (VE(xi,i)>test[i*(*antsim)+k]) test[i*(*antsim)+k]=VE(xi,i); 

	    c=(*p+i); 
	    if (VE(difX,i)>test[c*(*antsim)+k]) test[c*(*antsim)+k]=VE(difX,i);
	    c=2*(*p)+i; 
	    test[c*(*antsim)+k]=test[c*(*antsim)+k]+VE(ssrow,i)*dtime; 
	    if (*sim==1) {
	      c=(3*(*p)+i); 
	      if ((VE(dB,i))>test[c*(*antsim)+k]) test[c*(*antsim)+k]=VE(dB,i); 
	      c=4*(*p)+i; 
	      test[c*(*antsim)+k]=test[c*(*antsim)+k]+VE(ssrow2,i)*dtime; 
	    }
	  } 

	  if (supsup==1) {
	    for (j=s+1;j<*Ntimes;j++){ /* Beregning af sup_a,s | B(a+t)-B(t) - gam t | */ 
	      time1=times[j];  time2=time1-time; 
	      scl_vec_mult(time2/tau,vtmp,difX);
	      extract_row(Delta,j,xi); 
	      vec_subtr(xi,rowX,dB); 
	      vec_subtr(difX,dB,dB); 
	      for (i=0;i<=*p+1;i++) {
		c=5*(*p)+i; 
		if (fabs(VE(dB,i))>test[c*(*antsim)+k]) test[c*(*antsim)+k]=fabs(VE(dB,i)); 
	      }
	    } 
	  } /* supsup==1 */ 

	} 
      } 
    }  /* s=1..Ntimes k=1.. antsim, sim2==1*/  
    PutRNGstate();  /* to use R random normals */
  } /* sim==1 */

//  if (*sim==1) { }
    free_mat(Delta);
    free_mat(Delta2);
    free_mat(tmpM1);
    free_mat(tmpM2);

//  if (*robust==1) { }
    free_mat(varBL);
    for (j=0;j<*antpers;j++) { 
      free_mat(cumB[j]);
      free_vec(cumi[j]);
      free_vec(Base[j]);
      free_mat(cumBL[j]);
      free_vec(cumBLi[j]);
      free_mat(BLsubbetaLam[j]);
    } 

  free_vec(diag); free_vec(dB); free_vec(dN); free_vec(VdB); free_vec(AIXdN); free_vec(AIXlamt); 
  free_vec(ta); free_vec(bhatt); free_vec(pbhat); free_vec(plamt); free_vec(avx); free_vec(lrisk); 
  free_vec(ssrow2); free_vec(ssrow); free_vec(vtmp); free_vec(xi); free_vec(rowX); free_vec(difX); 
  free_vec(Btau); 
  free_mat(ldesignX); free_mat(A); free_mat(AI); free_mat(AIX); free_mat(cdesignX); 
  free_mat(XmavX); free_mat(cXmavX); free_mat(Aav); 
  free(coef); free(ps); free(imin);
  free(vcudif); free(Basei); 
}

void semibreslow(double *times,int *Ntimes,double *designX,int *nx,int *px,
		double *designG, int *ng,int *pg,int *antpers,double *start,
		double *stop,int *nb,double *bhat,double *cu,double *vcu,
		double *rvcu,double *gamma,double *Vgamma,double *robVgamma,double *b,
		int *degree,int *it,int *sim,int *antsim,double *test,
		int *rani,double *testOBS,int *status,int *id,double *schoen,
                double *simUt,double *Ut,int *weighted,int *robust)
{
  matrix *ldesignX, *A,*AI,*cdesignX,*ldesignG,*cdesignG;
  matrix *XmavX,*ZmavZ,*E2x,*E2z,*E2xz,*XX;
  matrix *S,*dCGam,*CGam,*ICGam,*VarKorG,*dC,*XZ,*ZZ,*ZZI,*XZAI; 
  matrix *Ct,*C[*Ntimes],*Acorb[*Ntimes],*ZXAI,*tmpM4; 
  matrix *RobVargam,*tmpM3;
//  matrix *W3t[*antpers],*W4t[*antpers],*AIxit[*antpers];
  vector *dB,*dN,*VdB,*AIXdN,*AIXlamt,*ta,*bhatt,*pbhat,*plamt;
  vector *difX,*korG,*pghat,*gam,*dgam,*ZGdN,*IZGdN,*ZGlamt,*IZGlamt;
  vector *zi,*z1,*lrisk,*avx,*avz,*rowG,*xi,*rowX,*rowZ,*tmpv2;
  int itt,i,j,k,s,c,count,pers=0,pmax,
        *imin=calloc(1,sizeof(int)), *coef=calloc(1,sizeof(int)),*ps=calloc(1,sizeof(int));
  double time,dummy,dtime,lam0t,S0,
	 *Basei=calloc((*antpers),sizeof(double)),
         *vcudif=calloc((*Ntimes)*(*px+2),sizeof(double)),dum2,rvarbase; 
  vector *cumi[*antpers],*W2[*antpers],*W3[*antpers],*Base[*antpers]; 
  matrix *W3t[*antpers],*W4t[*antpers],*AIxit[*antpers],*cumB[*antpers]; 

//  if (*robust==1){ }
    for (j=0;j<*antpers;j++) {
      malloc_mat(*Ntimes,*px,cumB[j]); 
      malloc_vec(*px,cumi[j]);
      malloc_mat(*Ntimes,*px,W3t[j]);
      malloc_mat(*Ntimes,*px+1,W4t[j]);
      malloc_mat(*Ntimes,*px,AIxit[j]); 
      malloc_vec(*Ntimes,Base[j]); 
      Basei[j]=0.0; 
      malloc_vec(*pg,W2[j]); 
      malloc_vec(*px,W3[j]);
    }

  malloc_mat(*antpers,*px,XmavX);
  malloc_mat(*antpers,*px,ldesignX);
  malloc_mat(*antpers,*px,cdesignX);
  malloc_mat(*antpers,*pg,ZmavZ);
  malloc_mat(*antpers,*pg,ldesignG);
  malloc_mat(*antpers,*pg,cdesignG);
  malloc_mat(*px,*px,XX);
  malloc_mat(*px,*px,E2x);
  malloc_mat(*px,*px,A);
  malloc_mat(*px,*px,AI);
  malloc_mat(*pg,*pg,tmpM3);
  malloc_mat(*pg,*pg,RobVargam);
  malloc_mat(*pg,*pg,E2z);
  malloc_mat(*pg,*pg,ZZ);
  malloc_mat(*pg,*pg,VarKorG);
  malloc_mat(*pg,*pg,ICGam);
  malloc_mat(*pg,*pg,CGam);
  malloc_mat(*pg,*pg,dCGam);
  malloc_mat(*pg,*pg,S);
  malloc_mat(*pg,*pg,ZZI);
  malloc_mat(*px,*pg,E2xz);
  malloc_mat(*px,*pg,XZ);
  malloc_mat(*px,*pg,XZAI);
  malloc_mat(*pg,*px,Ct);
  malloc_mat(*pg,*px,dC);
  malloc_mat(*pg,*px,ZXAI);
  malloc_mat(*pg,*px,tmpM4);
  for (j=0;j<*Ntimes;j++) {
    malloc_mat(*pg,*px,Acorb[j]);
    malloc_mat(*pg,*px,C[j]);
  }

  malloc_vec(*px,difX);
  malloc_vec(*px,xi);
  malloc_vec(*px,rowX);
  malloc_vec(*px,avx);
  malloc_vec(*px,korG);
  malloc_vec(*px,dB);
  malloc_vec(*px,VdB);
  malloc_vec(*px,AIXdN);
  malloc_vec(*px,AIXlamt);
  malloc_vec(*px,bhatt);
  malloc_vec(*pg,tmpv2);
  malloc_vec(*pg,rowZ);
  malloc_vec(*pg,avz);
  malloc_vec(*pg,rowG);
  malloc_vec(*pg,zi);
  malloc_vec(*pg,z1);
  malloc_vec(*pg,gam);
  malloc_vec(*pg,dgam);
  malloc_vec(*pg,ZGdN);
  malloc_vec(*pg,IZGdN);
  malloc_vec(*pg,ZGlamt);
  malloc_vec(*pg,IZGlamt);
  malloc_vec(*antpers,lrisk);
  malloc_vec(*antpers,dN);
  malloc_vec(*antpers,pbhat);
  malloc_vec(*antpers,pghat);
  malloc_vec(*antpers,plamt);
  malloc_vec(*nb,ta);

  coef[0]=1; ps[0]=(*px)+2; 
  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  for (j=0;j<*pg;j++){ 
    VE(gam,j)=gamma[j]; 
  }

  for (itt=0;itt<*it;itt++){
    mat_zeros(Ct);
    mat_zeros(CGam);
    vec_zeros(IZGdN);
    vec_zeros(IZGlamt);

    for (s=1;s<*Ntimes;s++){
      time=times[s];
      dtime=time-times[s-1]; 
      vec_zeros(lrisk); 
      mat_zeros(ldesignX);
      mat_zeros(ldesignG);

      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++){
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<pmax;j++) {
	    if (j<*px){ 
	      ME(ldesignX,id[c],j)= designX[j*(*nx)+c];
	    }
	    if (j<*pg){ 
	      ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; 
	    }
	  }
	  VE(lrisk,id[c])=1.0; 
	  if (time==stop[c] && status[c]==1) pers=id[c]; 
	  count=count+1; 
	}
      }

      for(j=0;j<*nb;j++){ 
	VE(ta,j)=fabs(bhat[j]-time);
      }
      dummy=vec_min(ta,imin); 
      lam0t=bhat[1*(*nb)+(*imin)];
      for(j=2;j<=*px+1;j++){ 
	VE(bhatt,j-2)=bhat[j*(*nb)+(*imin)];
      }

      Mv(ldesignX,bhatt,pbhat); 
      Mv(ldesignG,gam,pghat);

      for (j=0;j<*antpers;j++) {
	VE(plamt,j)=VE(lrisk,j)*exp(VE(pbhat,j)+VE(pghat,j)); 
      }

      S0=vec_sum(plamt); 
      vM(ldesignX,plamt,avx); 
      scl_vec_mult((1/S0),avx,avx);
      vM(ldesignG,plamt,avz); 
      scl_vec_mult((1/S0),avz,avz);

      for (j=0;j<pmax;j++) {
	for (i=0;i<pmax;i++) {
	  if ((j<*px) && (i<*px)) ME(E2x,j,i)=VE(avx,i)*VE(avx,j)*S0;
	  if ((j<*px) && (i<*pg)) ME(E2xz,j,i)=VE(avx,j)*VE(avz,i)*S0;
	  if ((j<*pg) && (i<*pg)) ME(E2z,j,i)=VE(avz,i)*VE(avz,j)*S0; 
	}
      } 

      for (j=0;j<*antpers;j++) { 
	extract_row(ldesignX,j,dB); 
	scl_vec_mult(VE(plamt,j),dB,dB); 
	replace_row(cdesignX,j,dB);
	extract_row(ldesignG,j,rowG); 
	scl_vec_mult(VE(plamt,j),rowG,rowG); 
	replace_row(cdesignG,j,rowG); 
      }
				
      MtA(cdesignX,ldesignX,XX); 
      mat_subtr(XX,E2x,A); 
      invert(A,AI); 
      extract_row(ldesignX,pers,dB); 
      vec_subtr(dB,avx,dB); 
      Mv(AI,dB,AIXdN); 

      MtA(ldesignG, cdesignG, ZZ); 
      mat_subtr(ZZ,E2z,ZZ); 
      MtA(ldesignX, cdesignG, XZ); 
      mat_subtr(XZ,E2xz,XZ); 

      MtA(XZ,AI,ZXAI); 
      scl_mat_mult(1/(S0*lam0t),ZXAI,tmpM4); 
      /* sm_mlt(dtime,tmpM4,tmpM4); */ 
      mat_add(tmpM4,Ct,Ct); 
      mat_copy(Ct,C[s]); 

      mat_copy(ZXAI,Acorb[s]); 

      MxA(ZXAI,XZ,E2z); 
      mat_subtr(ZZ,E2z,dCGam); 
      scl_mat_mult(1/S0,dCGam,dCGam); 
      mat_add(CGam,dCGam,CGam); 

      extract_row(ldesignG,pers,zi); 
      vec_subtr(zi,avz,zi); 
      vM(XZ,AIXdN,rowG); 
      vec_subtr(zi,rowG,ZGdN); 
      vec_add(ZGdN,IZGdN,IZGdN); 

      cu[1*(*Ntimes)+s]=cu[1*(*Ntimes)+s-1]+(1/S0);
      for (k=2;k<=*px+1;k++) {
	cu[k*(*Ntimes)+s]= cu[k*(*Ntimes)+s-1]+ (VE(bhatt,k-2)*dtime)+(VE(AIXdN,k-2)/lam0t); 
	vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s-1]+ dtime*ME(AI,k-2,k-2)/lam0t;
      }

      if (itt==(*it-1)) { 
	for (i=0;i<*antpers;i++) {
	  extract_row(ldesignX,i,xi); 
	  vec_subtr(xi,avx,xi); 
	  Mv(AI,xi,rowX); 
	  replace_row(AIxit[i],s,rowX);  
	}
      } 
    } /* s=1,...,Ntimes */ 

    invert(CGam,ICGam); 
    Mv(ICGam,IZGdN,dgam); 
    vec_add(dgam,gam,gam); 
    /*
      v_output(gam); 
    */ 

    for (s=1;s<*Ntimes;s++) {
      time=times[s]; dtime=time-times[s-1]; 

      vM(C[s],dgam,korG);

      cu[s]=times[s]; vcu[s]=times[s]; rvcu[s]=times[s];
      for (k=2;k<=*px+1;k++) {
	cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s]-VE(korG,k-2);
      }

      if (itt==(*it-1)) { 
	MxA(ICGam,C[s],tmpM4);
	MtA(tmpM4,C[s],VarKorG); 
	for (k=2;k<=*px+1;k++){
	  vcu[k*(*Ntimes)+s]=vcu[k*(*Ntimes)+s]+ME(VarKorG,k-2,k-2);
	}
      } 
    }
    smoothB(cu,Ntimes,ps,bhat,nb,b,degree,coef);
  } /*itt lokke */ 

  /* ==================ROBUST terms ================================= */ 
  if (*robust==1){
    for (s=1;s<*Ntimes;s++){
      time=times[s]; dtime=time-times[s-1]; 
      vec_zeros(lrisk); 
      mat_zeros(ldesignX); 
      mat_zeros(ldesignG); 

      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++){
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<pmax;j++) {
	    if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	    if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; 
	  }
	  VE(lrisk,id[c])=1.0;  
	  if (time==stop[c]) pers=id[c]; 
	  count=count+1; 
	}
      }

      for(j=0;j<*nb;j++){
	VE(ta,j)=fabs(bhat[j]-time);
      }
      dummy=vec_min(ta,imin); 
      lam0t=bhat[1*(*nb)+(*imin)];
      for(j=2;j<=*px+1;j++) {
	VE(bhatt,j-2)=bhat[j*(*nb)+(*imin)];
      }

      Mv(ldesignX,bhatt,pbhat); 
      Mv(ldesignG,gam,pghat);

      for (j=0;j<*antpers;j++) {
	VE(plamt,j)=VE(lrisk,j)*exp(VE(pbhat,j)+VE(pghat,j)); 
      }

      S0=vec_sum(plamt); 
      vM(ldesignX,plamt,avx); 
      scl_vec_mult((1/S0),avx,avx);
      vM(ldesignG,plamt,avz); 
      scl_vec_mult((1/S0),avz,avz);

      rvarbase=0; 
      for (j=0;j<*antpers;j++) {
	extract_row(ldesignX,j,rowX); 
	vec_subtr(rowX,avx,rowX);
	extract_row(ldesignG,j,rowZ); 
	vec_subtr(rowZ,avz,rowZ);

	Mv(Acorb[s],rowX,tmpv2); 
	vec_subtr(rowZ,tmpv2,tmpv2); 

	if (j==pers) vec_add(tmpv2,W2[j],W2[j]);
	dummy=VE(plamt,j)/S0; 
	scl_vec_mult(dummy,tmpv2,rowZ); 
	vec_subtr(W2[j],rowZ,W2[j]); 

	extract_row(AIxit[j],s,rowX); 
	scl_vec_mult(1/lam0t,rowX,rowX); 
	if (j==pers) vec_add(rowX,W3[j],W3[j]); 

	scl_vec_mult(dummy,rowX,rowX); 
	vec_subtr(W3[j],rowX,W3[j]); 
	replace_row(W3t[j],s,W3[j]);  

	extract_row(AIxit[j],s,rowX);  /* Baseline sum iid */
	vec_star(avx,rowX,xi); 
	dum2=vec_sum(xi); 
	dum2=(1/S0)-dum2; 
	if (j==pers) Basei[j]=dum2+Basei[j]; 
	Basei[j]=Basei[j]-dummy*dum2;   
	VE(Base[j],s)=Basei[j];  
	ME(W4t[j],s,0)=Basei[j];
	rvarbase=rvarbase+Basei[j]*Basei[j]; 
      } /* j=1..antpers */
      rvcu[1*(*Ntimes)+s]=rvarbase; 
    } /* s loekke */ 

    for (s=1;s<*Ntimes;s++){
      vec_zeros(VdB); 
      for (i=0;i<*antpers;i++){
	/* B(t) iid rep */ 
	Mv(ICGam,W2[i],tmpv2); 
	vM(C[s],tmpv2,rowX);
	extract_row(W3t[i],s,xi); 
	vec_subtr(xi,rowX,difX); 

	/* set_row(W4t[i],s,difX); */ 
	for (k=1;k<*px+1;k++) {
	  ME(W4t[i],s,k)=VE(difX,k-1); 
	}
	vec_star(difX,difX,xi); 
	vec_add(xi,VdB,VdB);

	if (s==1) { 
	  for (j=0;j<*pg;j++) {
	    for (k=0;k<*pg;k++){
	      ME(RobVargam,j,k)=ME(RobVargam,j,k)+VE(W2[i],j)*VE(W2[i],k);
	    }
	    /* Rprintf(" %ld %lf %lf \n",i,W2[i]->ve[0],RobVargam->me[0][0]); */ 
	  }
	}
      } /* i =1 ..Antpers */
      for (k=2;k<=*px+1;k++) {
	rvcu[k*(*Ntimes)+s]=VE(VdB,k-2);
      }
    }  /*  s=1 ..Ntimes */ 

    MxA(RobVargam,ICGam,tmpM3); 
    MxA(ICGam,tmpM3,RobVargam);
  }

  for (j=0;j<*pg;j++) {
    gamma[j]=VE(gam,j);
    for (k=0;k<*pg;k++) {
      Vgamma[k*(*pg)+j]=ME(ICGam,j,k);
      robVgamma[k*(*pg)+j]=ME(RobVargam,j,k);
    }
  }

  if (*sim==1) {
    ps[0]=(*px)+1; 
    comptest(times,Ntimes,ps,cu,rvcu,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antpers);
  }

  cu[0]=times[0]; vcu[0]=times[0];

  /* korrektion for lam0t i \int beta(s) lam0t(s) ds */ 
  /* 
     for (s=1;s<*Ntimes;s++)
     {
     time=times[s]; dtime=time-times[s-1]; 
     for(j=0;j<*nb;j++) ta->ve[j]=fabs(bhat[j]-time);
     dummy=v_min(ta,imin); lam0t=bhat[1*(*nb)+(*imin)]; 
     for(j=2;j<=(*px)+1;j++) bhatt->ve[j-2]=bhat[j*(*nb)+(*imin)];
	     
     for (k=2;k<=(*px)+1;k++) {
     cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+dtime*bhatt->ve[k-2]/lam0t;  }; 
     } 
  */ 

  free_mat(ldesignX); free_mat(A); free_mat(AI); free_mat(cdesignX); free_mat(ldesignG); 
  free_mat(XmavX); free_mat(ZmavZ); free_mat(E2x); free_mat(E2xz); free_mat(XX); 
  free_mat(S); free_mat(dCGam); free_mat(CGam); free_mat(ICGam); free_mat(VarKorG); 
  free_mat(dC); free_mat(XZ); free_mat(ZZ); free_mat(ZZI); free_mat(XZAI);
  free_mat(Ct); free_mat(ZXAI); free_mat(tmpM4); free_mat(RobVargam); free_mat(tmpM3); 

  free_vec(dB); free_vec(VdB); free_vec(AIXdN); free_vec(AIXlamt); free_vec(ta); 
  free_vec(bhatt); free_vec(pbhat); free_vec(plamt); free_vec(difX); free_vec(korG); 
  free_vec(pghat); free_vec(gam); free_vec(dgam); free_vec(ZGdN); free_vec(IZGdN); 
  free_vec(IZGlamt); free_vec(zi); free_vec(z1); free_vec(lrisk); free_vec(avx);
  free_vec(avz); free_vec(rowG); free_vec(xi); free_vec(rowX); free_vec(rowZ);
  free_vec(tmpv2); 

  for (j=0;j<*Ntimes;j++) {
    free_mat(Acorb[j]);
    free_mat(C[j]);
  }

//  if (*robust==1){ }
//
    for (j=0;j<*antpers;j++) {
      free_vec(Base[j]);
      free_vec(cumi[j]);
      free_mat(cumB[j]);
      free_vec(W2[j]);
      free_vec(W3[j]);
      free_mat(W3t[j]);
      free_mat(W4t[j]);
      free_mat(AIxit[j]);
    }
  free(vcudif); free(Basei); 
  free(coef); free(ps); free(imin);
}
