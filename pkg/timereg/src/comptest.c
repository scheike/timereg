#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

void comptest(times,Ntimes,px,cu,vcu,vcudif,antsim,test,testOBS,Ut,
	      simUt,W4t,weighted,antpers)
     double *times,*cu,*vcu,*vcudif,*test,*testOBS,*Ut,*simUt;
     int    *px,*Ntimes,*antsim,*weighted,*antpers;
     matrix  **W4t;
{
  matrix *Delta,*tmpM1;
  vector *tmpv1,*rowX,*xi,*difX,*ssrow,*VdB;
  int i,k,l,s,c;
  double xij,vardif,tau,time,dtime,random,fabs(),sqrt();
  double norm_rand();
  void GetRNGstate(),PutRNGstate();

  /* float gasdev(),expdev(),ran1(); */ 

  malloc_vec(*px,tmpv1);
  malloc_vec(*px,rowX);
  malloc_vec(*px,xi);
  malloc_vec(*px,difX);
  malloc_vec(*px,ssrow);
  malloc_vec(*px,VdB); 
  malloc_mat(*Ntimes,*px,Delta);
  malloc_mat(*Ntimes,*px,tmpM1); 

// printf("Simulations start N= %ld \n",(long int) *antsim); 

  GetRNGstate();  /* to use R random normals */ 

  tau=times[(*Ntimes-1)]-times[0];
  Ut[0]=times[0]; 

  if (*weighted>=1) {
    for (s=1;s<*Ntimes;s++) {
      vec_zeros(VdB);
      for (i=0;i<*antpers;i++) {
	extract_row(W4t[i],s,tmpv1);  
	extract_row(W4t[i],*Ntimes-1,rowX);
	scl_vec_mult((times[s]-times[0])/tau,rowX,rowX);
	vec_subtr(tmpv1,rowX,difX);
	vec_star(difX,difX,rowX);
	vec_add(rowX,VdB,VdB); 
      }
      for (k=1;k<=*px;k++) {
	vcudif[k*(*Ntimes)+s]=VE(VdB,k-1); 
      }
    }
  } /* weighted==1 */ 

  for (i=1;i<=*px;i++){ 
    VE(rowX,i-1)=cu[i*(*Ntimes)+(*Ntimes-1)];
  }

  /* Computation of observed teststatistics */ 
  for (s=1;s<*Ntimes;s++){
    time=times[s];dtime=times[s]-times[s-1];
    scl_vec_mult((time-times[0])/tau,rowX,difX);
   
    for (i=1;i<=*px;i++) {
      xij=fabs(cu[i*(*Ntimes)+s])/sqrt(vcu[i*(*Ntimes)+s]);
      /*
	printf(" %lf %lf %ld \n",xij,testOBS[i-1],i);  
	printf(" %lf %lf \n",cu[i*(*Ntimes)+s],vcu[i*(*Ntimes)+s]);  
      */

      if (xij>testOBS[i-1]) {
	testOBS[i-1]=xij;
      }
    } 

    for (i=1;i<=*px;i++){ 
      VE(xi,i-1)=cu[i*(*Ntimes)+s];
    }
    vec_subtr(xi,difX,difX); 
    vec_star(difX,difX,ssrow); 

    Ut[s]=time; 


    for (i=0;i<*px;i++) { 
      if (*weighted>=1) {
	vardif=vcudif[(i+1)*(*Ntimes)+s]; 
      } else {
	vardif=1;
      }
      if (*weighted>=1)  {

	if ((s>*weighted) && (s<*Ntimes-*weighted)){  
	  VE(difX,i)=VE(difX,i)/sqrt(vardif);
	} else {
	  VE(difX,i)=0.0;
	}
      } else {
	VE(difX,i)=VE(difX,i);
      }

      Ut[(i+1)*(*Ntimes)+s]=VE(difX,i);

      c=(*px)+i;
      if (fabs(VE(difX,i))>testOBS[c]) {
	testOBS[c]=fabs(VE(difX,i));
      }
      c=2*(*px)+i;
      if ((s>*weighted) && (s<*Ntimes-*weighted)){
	testOBS[c]=testOBS[c]+VE(ssrow,i)*dtime/vardif;
      } 
    }
  } 

  /* for (i=0;i<3*(*px);i++) printf(" %lf \n",testOBS[i]);  */



  /* simulation of testprocesses and teststatistics */ 
  for (k=1;k<=*antsim;k++) {

    mat_zeros(Delta); 
    vec_zeros(tmpv1); 

    for (i=0;i<*antpers;i++) {

      /* random=gasdev(&idum);  */ 
      random=norm_rand(); 

      scl_mat_mult(random,W4t[i],tmpM1);
      mat_add(tmpM1,Delta,Delta); 

    }

    extract_row(Delta,*Ntimes-1,tmpv1);     

    for (s=1;s<*Ntimes;s++) { 

      time=times[s]-times[0]; 
      dtime=times[s]-times[s-1]; 
      scl_vec_mult(time/tau,tmpv1,xi);

      extract_row(Delta,s,rowX); 
      vec_subtr(rowX,xi,difX); 
      vec_star(difX,difX,ssrow); 

      for (i=0;i<*px;i++) { 
	VE(xi,i)=fabs(ME(Delta,s,i))/sqrt(vcu[(i+1)*(*Ntimes)+s]);

	if (VE(xi,i)>test[i*(*antsim)+k]){ 
	  test[i*(*antsim)+k]=VE(xi,i);
	}

	if (*weighted>=1) {
	  vardif=vcudif[(i+1)*(*Ntimes)+s];
	}  else {
	  vardif=1; 	
	}

	if (*weighted>=1)  {
	  if ((s>*weighted) && (s<*Ntimes-*weighted)){
	    VE(difX,i)=VE(difX,i)/sqrt(vardif);
	  } else {
	    VE(difX,i)=0.0;
	  }
	} else {
	  VE(difX,i)=VE(difX,i);
	}

	if (k<51) {
	  l=(k-1)*(*px)+i; 
	  simUt[l*(*Ntimes)+s]=VE(difX,i);
	}

	c=(*px+i);   
	VE(difX,i)=fabs(VE(difX,i));
	if (VE(difX,i)>test[c*(*antsim)+k]) {
	  test[c*(*antsim)+k]=VE(difX,i);
	}
	c=2*(*px)+i; 
	if ((s>*weighted) && (s<*Ntimes-*weighted)) {

	  test[c*(*antsim)+k]=test[c*(*antsim)+k]+VE(ssrow,i)*dtime/vardif; 
	}

      }
    }  /* s=1..Ntimes */ 
  }  /* k=1..antsim */ 

  PutRNGstate();  /* to use R random normals */

  free_mat(Delta);
  free_mat(tmpM1);
  free_vec(VdB);
  free_vec(rowX);
  free_vec(difX);
  free_vec(xi);
  free_vec(tmpv1);
  free_vec(ssrow); 

}


void comptestM(times,Ntimes,px,cu,vcu,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antpers,cu0,argmax)
double *times,*cu,*vcu,*vcudif,*test,*testOBS,*Ut,*simUt,*cu0,*argmax;
int    *px,*Ntimes,*antsim,*weighted,*antpers;
matrix  *W4t[];
{
  matrix *Delta,*tmpM1;
  vector *tmpv1,*rowX,*xi,*difX,*ssrow,*VdB;
  int i,k,l,s,c,u,t;
  double xij,vardif,tau,time,dtime,random,fabs(),sqrt();
  double ixij,mu,ms,mt,tu,ts,tt,uhat,dmus,dmts,icxij; 
  double norm_rand();
  void GetRNGstate(),PutRNGstate();
  /* float gasdev(),expdev(),ran1(); */

  malloc_vec(*px,tmpv1);
  malloc_vec(*px,rowX);
  malloc_vec(*px,xi);
  malloc_vec(*px,difX);
  malloc_vec(*px,ssrow);
  malloc_vec(*px,VdB); 
  malloc_mat(*Ntimes,*px,Delta); 
  malloc_mat(*Ntimes,*px,tmpM1);

  printf("Simulations start N= %ld \n",(long int) *antsim); 

  GetRNGstate();  /* to use R random normals */

  tau=times[(*Ntimes-1)]-times[0];

  if (*weighted>=1) {
    for (s=1;s<*Ntimes;s++) {
      vec_zeros(VdB);
      for (i=0;i<*antpers;i++) {
        extract_row(W4t[i],s,tmpv1);  
        extract_row(W4t[i],*Ntimes-1,rowX);
        scl_vec_mult((times[s]-times[0])/tau,rowX,rowX);
        vec_subtr(tmpv1,rowX,difX);
        vec_star(difX,difX,rowX);
	vec_add(rowX,VdB,VdB); 
      }
      for (k=1;k<=*px;k++) {
	vcudif[k*(*Ntimes)+s]=VE(VdB,k-1); 
      }
    }
  } /* weighted==1 */ 

  for (i=1;i<=*px;i++) {
    VE(rowX,i-1)=cu[i*(*Ntimes)+(*Ntimes-1)];
  }

  uhat= VE(rowX,0)/tau; 
  Ut[0]=times[0]; 

  /* Computation of observed teststatistics */ 
  for (s=1;s<*Ntimes;s++){
    time=times[s]-times[0]; 
    dtime=times[s]-times[s-1]; 
    scl_vec_mult(time/tau,rowX,difX);
   
    for (i=1;i<=*px;i++) {
      xij=fabs(cu[i*(*Ntimes)+s])/sqrt(vcu[i*(*Ntimes)+s]);
      if (xij>testOBS[i-1]) {
	testOBS[i-1]=xij;
      }

      c=3*(*px);
      testOBS[c]=testOBS[c]+cu[i*(*Ntimes)+s]*cu[i*(*Ntimes)+s]*dtime; 
      /* printf(" %lf \n",testOBS[c]);  */
    } 

    for (i=1;i<=*px;i++){ 
      VE(xi,i-1)=cu[i*(*Ntimes)+s];
    }
    vec_subtr(xi,difX,difX); 
    vec_star(difX,difX,ssrow); 

    Ut[s]=times[s]; 

    for (i=0;i<*px;i++) { 
      if (*weighted>=1){ 
	vardif=vcudif[(i+1)*(*Ntimes)+s];  
      }else{ 
	vardif=1; 
      }
      if (*weighted>=1)  {
	if ((s>*weighted) && (s<*Ntimes-*weighted)) {
	  VE(difX,i)=VE(difX,i)/sqrt(vardif);
	} else {
	  VE(difX,i)=0;
	}
      } else {
	VE(difX,i)=VE(difX,i); 
      }

      Ut[(i+1)*(*Ntimes)+s]=VE(difX,i);

      c=(*px);
      if (fabs(VE(difX,i))>testOBS[c]) {
	testOBS[c]=fabs(VE(difX,i));
      }
      c=2*(*px);
      if ((s>*weighted) && (s<*Ntimes-*weighted)) {
	testOBS[c]=testOBS[c]+VE(ssrow,i)*dtime/vardif; 
      }
    }

    /* konveksitet */ 
    if (s > *Ntimes){
      ts=times[s]; 
      ms=cu[1*(*Ntimes)+s];
      for (t=s+1;t<*Ntimes;t++)
	{
	  tt=times[t]; 
	  mt=cu[1*(*Ntimes)+t];
	  ixij=0; 
	  icxij=0; 
	  for (u=s;u<t;u++){
	    tu=times[u]; 
	    mu=cu[1*(*Ntimes)+u];
	    dtime=times[u]-times[u-1]; 

	    xij=(mu-ms)-uhat*(tu-ts); 
	    c=3*(*px);
	    if (fabs(xij)>testOBS[c])  { 
	      testOBS[c]=fabs(xij); 
	      /* printf(" %lf %lf %lf %lf \n",ts,tt,tu,xij);   */ 
	    }
	    ixij=ixij+dtime*xij*xij; 

	    xij=(mu-ms)-(mt-ms)*(tu-ts)/(tt-ts);  
	    c=5*(*px);
	    if (xij>testOBS[c]) { 
	      testOBS[c]=xij; 
	    }
	    icxij=icxij+dtime*xij; 
	  } 
	  c=4*(*px);
	  if (ixij>testOBS[c]){  
	    testOBS[c]=ixij; 
	  }
	  c=6*(*px);
	  if (icxij>testOBS[c]){  
	    testOBS[c]=icxij;
	  }
	} 
    }
  } 

  /* simulation of testprocesses and teststatistics */ 
  for (k=1;k<=*antsim;k++) {
    mat_zeros(Delta); 
    vec_zeros(tmpv1); 
    for (i=0;i<*antpers;i++) {
      /* random=gasdev(&idum);  */
      random=norm_rand();
      scl_mat_mult(random,W4t[i],tmpM1); 
      mat_add(tmpM1,Delta,Delta); 
    }

    extract_row(Delta,*Ntimes-1,tmpv1); 

    uhat=VE(tmpv1,0)/tau;  

    for (s=1;s<*Ntimes;s++) { 

      time=times[s]-times[0]; 
      dtime=times[s]-times[s-1]; 
      scl_vec_mult(time/tau,tmpv1,xi);
      extract_row(Delta,s,rowX); 
      vec_subtr(rowX,xi,difX); 
      vec_star(difX,difX,ssrow); 

      for (i=0;i<*px;i++) { 
	VE(xi,i)=fabs(ME(Delta,s,i))/sqrt(vcu[(i+1)*(*Ntimes)+s]);
	if (VE(xi,i)>test[i*(*antsim)+k]){ 
	  test[i*(*antsim)+k]=VE(xi,i); 
	}
	c=3*(*px);
	test[c*(*antsim)+k]=test[c*(*antsim)+k]+ME(Delta,s,i)*ME(Delta,s,i)*dtime; 

	if (*weighted>=1){ 
	  vardif=vcudif[(i+1)*(*Ntimes)+s];
	}  else {
	  vardif=1; 
	}
	if (*weighted>=1)  {
	  if ((s>*weighted) && (s<*Ntimes-*weighted)){
	    VE(difX,i)=VE(difX,i)/sqrt(vardif);
	  } else{ 
	    VE(difX,i)=0.0;
	  }
	} else { 
	  VE(difX,i)=VE(difX,i); 
	}

	if (k<51) {
	  l=(k-1)*(*px)+i; 
	  simUt[l*(*Ntimes)+s]=VE(difX,i);
	}

	c=(*px+i);   
	VE(difX,i)=fabs(VE(difX,i));
	if (VE(difX,i)>test[c*(*antsim)+k]){ 
	  test[c*(*antsim)+k]=VE(difX,i);
	}
	c=2*(*px)+i; 
	if ((s>*weighted) && (s<*Ntimes-*weighted)){
	  test[c*(*antsim)+k]=test[c*(*antsim)+k]+VE(ssrow,i)*dtime/vardif; 
	}
      }

      if (s>*Ntimes) {
	ts=times[s]; 
	ms=ME(Delta,0,s); 
	for (t=s+1;t<*Ntimes;t++){ 
	  tt=times[t]; 
	  mt=ME(Delta,0,t);
	  ixij=0;  
	  icxij=0; 
	  for (u=s;u<t;u++){
	    tu=times[u]; 
	    mu=ME(Delta,0,u);
	    dtime=times[u]-times[u-1]; 

	    xij=(mu-ms)-uhat*(tu-ts); 
	    c=3*(*px);
	    if (fabs(xij)>test[c*(*antsim)+k]) { 
	      test[c*(*antsim)+k]=fabs(xij); 
	      /* printf("local %lf %lf %lf %lf \n",ts,tt,tu,xij); */  
	    }
	    ixij=ixij+dtime*xij*xij; 
	    dmus=cu0[i*(*Ntimes)+u]- cu0[i*(*Ntimes)+s];
	    dmts=cu0[i*(*Ntimes)+t]- cu0[i*(*Ntimes)+s];
	    xij=(mu-ms)-(mt-ms)*(tu-ts)/(tt-ts);  
	    xij=dmus+(mu-ms)-(dmts+mt-ms)*(tu-ts)/(tt-ts);  
	    c=5*(*px);
	    if (xij>test[c*(*antsim)+k]) { 
	      test[c*(*antsim)+k]=xij; 
	      /* printf("conveks %lf %lf %lf %lf \n",ts,tt,tu,xij); */  
	    }
	    icxij=icxij+dtime*xij; 
	  } 
	  c=4*(*px);
	  if (ixij>test[c*(*antsim)+k]){  
	    test[c*(*antsim)+k]=ixij; 
	  }
	  c=6*(*px);
	  if (icxij>test[c*(*antsim)+k]){  
	    test[c*(*antsim)+k]=icxij; 
	  }
	} 
      }
    }  /* s=1..Ntimes */ 
  }  /* k=1..antsim */ 

  PutRNGstate();  /* to use R random normals */

  free_mat(Delta);
  free_mat(tmpM1); 
  free_vec(VdB);
  free_vec(rowX);
  free_vec(difX);
  free_vec(xi);
  free_vec(tmpv1);
  free_vec(ssrow); 
}

