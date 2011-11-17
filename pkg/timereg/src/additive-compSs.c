//#include <stdio.h>
#include <math.h>
#include "matrix.h"

void compSs(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,id,status,deltaweight,intZHZ,intZHdN,silent)
double *designX,*alltimes,*start,*stop,*intZHZ,*intZHdN,*designG;
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*deltaweight,*id,*silent;
{ // {{{
  matrix *X,*A,*AI,*AIXW,*dCGam,*CGam,*Ct,*ICGam,*XWZ,*ZWZ,*XWZAI,*tmpM4,*tmpM2;
  vector *xi,*tmpv2,*tmpv1,*PLScomp,*Xi,*dA,*rowX,*AIXWdN,*korG,*rowZ,*gam,*ZHdN,
    *IZHdN,*zi;
  int j,k,l,c,s,count,pers=0,pmax,*ipers=calloc(*Ntimes,sizeof(int)); 
  int stat,*ls=calloc(*Ntimes,sizeof(int)); 
  double time,dtime,fabs(),sqrt();

  malloc_mats(*antpers,*px,&X,NULL);
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*px,*antpers,&AIXW,NULL);
  // malloc_mats(*antpers,*pg,&Z,NULL); 
  malloc_mats(*pg,*pg,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*pg,&Ct,&XWZ,&XWZAI,NULL);
  malloc_mat(*px,*pg,tmpM4); 

  malloc_vecs(*px,&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,NULL);
  malloc_vecs(*pg,&zi,&tmpv2,&rowZ,&gam,&ZHdN,&IZHdN,NULL);
  malloc_vecs(*antpers,&PLScomp,&Xi,NULL);

  if (*px>=*pg) pmax=*px; else pmax=*pg; 

    mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN);  

    // Rprintf(" test \n"); 
    for (s=1;s<*Nalltimes;s++){
    // Rprintf(" test %d %d %d  \n",s,*antpers,*nx); 
      time=alltimes[s]; 
      dtime=time-alltimes[s-1]; mat_zeros(A); stat=0;  
      mat_zeros(ZWZ); mat_zeros(XWZ); 

      l=0; stat=0; 
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
//		Rprintf("times  %lf  %lf %lf \n",time,start[c],stop[c]); 
	if ((start[c]<time) && (stop[c]>=time)) {
//		Rprintf("under risk %d %d %d \n",c,id[c],count); 
	    for(j=0;j<pmax;j++)  {
	    if (j<*px) { 
	    for(k=0;k<*px;k++) 
		    ME(A,j,k)+= designX[j*(*nx)+c]*designX[k*(*nx)+c];
	    for(k=0;k<*pg;k++) 
		    ME(XWZ,j,k)+= designX[j*(*ng)+c]*designG[k*(*ng)+c]; 
	    }
	    if (j<*pg) 
	    for(k=0;k<*pg;k++) 
		    ME(ZWZ,j,k)+= designG[k*(*ng)+c]*designG[j*(*ng)+c]; 
	  }

	  if (time==stop[c] && status[c]==1)
	    {pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;
	    for(k=0;k<*pg;k++) VE(zi,k)=designG[k*(*ng)+c]; 
	    for(k=0;k<*px;k++) VE(xi,k)=designX[k*(*ng)+c]; 
	    }
	  count=count+1; 
      } }

      //Rprintf(" s er %d \n",s); 
      //print_mat(A); print_mat(ZWZ); print_mat(XWZ); 

     // MtA(X,X,A); 
      invertS(A,AI,silent[0]); 
      if (ME(AI,0,0)==0 && *silent==0) Rprintf("time %lf X'X singular \n",time); 
     // MtA(Z,Z,ZWZ);MtA(X,Z,XWZ);
      MxA(AI,XWZ,XWZAI);

      MtA(XWZAI,XWZ,tmpM2); mat_subtr(ZWZ,tmpM2,dCGam); 
      scl_mat_mult(dtime,dCGam,dCGam); 
      if (*deltaweight==0) scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

      if (stat==1) {
   // 	extract_row(X,pers,tmpv1); 
        Mv(AI,xi,AIXWdN); 
   //	extract_row(Z,pers,zi); 
	vM(XWZ,AIXWdN,tmpv2);
	vec_subtr(zi,tmpv2,ZHdN);
	if (*deltaweight==0) scl_vec_mult(dtime,ZHdN,ZHdN); 
	vec_add(ZHdN,IZHdN,IZHdN); 
      }

      // scl_mat_mult(dtime,XWZAI,tmpM4);mat_add(tmpM4,Ct,Ct); 
    } /* s =1...Ntimes */ 

    // invertS(CGam,ICGam,silent[0]); Mv(ICGam,IZHdN,gam); 
    //if (ME(ICGam,0,0)==0 && *silent==0) Rprintf(" intZHZ  singular\n"); 
    //  print_mat(CGam); print_vec(IZHdN); 

    for(k=0;k<*pg;k++)  {
       intZHdN[k]=VE(IZHdN,k); 
       for(j=0;j<*pg;j++) intZHZ[k*(*pg)+j]=ME(CGam,k,j); 
    }



  free_mats(&X,&A,&AI,&AIXW,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam, 
		  &Ct,&XWZ,&XWZAI, &tmpM4,NULL); 
  free_vecs(&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&zi,&tmpv2,&rowZ,&gam,
		  &ZHdN,&IZHdN,&PLScomp,&Xi,NULL);
  free(ipers); free(ls); 
} // }}}

void compSsrev(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,id,status,deltaweight,intZHZ,intZHdN,silent)
double *designX,*alltimes,*start,*stop,*intZHZ,*intZHdN,*designG;
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*deltaweight,*id,*silent;
{ // {{{
// {{{

  matrix *X,*A,*AI,*AIXW,*dCGam,*CGam,*Ct,*ICGam,*XWZ,*ZWZ,*XWZAI,*tmpM4,*tmpM2;
  vector *xi,*tmpv2,*tmpv1,*PLScomp,*Xi,*dA,*rowX,*AIXWdN,*korG,*rowZ,*gam,*ZHdN,
    *IZHdN,*zi;
  int sstop,j,k,l,c,s,count,pers=0,pmax,*ipers=calloc(*Ntimes,sizeof(int)); 
  int stat,*ls=calloc(*Ntimes,sizeof(int)); 
  double time,dtime,fabs(),sqrt();

  malloc_mats(*antpers,*px,&X,NULL);
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*px,*antpers,&AIXW,NULL);
  // malloc_mats(*antpers,*pg,&Z,NULL); 
  malloc_mats(*pg,*pg,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*pg,&Ct,&XWZ,&XWZAI,NULL);
  malloc_mat(*px,*pg,tmpM4); 

  malloc_vecs(*px,&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,NULL);
  malloc_vecs(*pg,&zi,&tmpv2,&rowZ,&gam,&ZHdN,&IZHdN,NULL);
  malloc_vecs(*antpers,&PLScomp,&Xi,NULL);
  // }}}

  if (*px>=*pg) pmax=*px; else pmax=*pg; 

    //mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN);  
    //mat_zeros(A); mat_zeros(ZWZ); mat_zeros(XWZ); 

    count=nx[0]-1; 
    for (s=(*Nalltimes)-1;s>0;s=s-1){
      sstop=0; 
    // Rprintf(" test %d %d %d  \n",s,*antpers,*nx); 
      time=alltimes[s]; dtime=time-alltimes[s-1]; stat=0;  

      l=0; stat=0;  
      if (1==0) {
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) { // {{{
//		Rprintf("times  %lf  %lf %lf \n",time,start[c],stop[c]); 
	if ((start[c]<time) && (stop[c]>=time)) {
//		Rprintf("under risk %d %d %d \n",c,id[c],count); 
	    for(j=0;j<pmax;j++)  {
	    if (j<*px) { 
	    for(k=0;k<*px;k++) 
		    ME(A,j,k)+= designX[j*(*nx)+c]*designX[k*(*nx)+c];
	    for(k=0;k<*pg;k++) 
		    ME(XWZ,j,k)+= designX[j*(*ng)+c]*designG[k*(*ng)+c]; 
	    }
	    if (j<*pg) 
	    for(k=0;k<*pg;k++) 
		    ME(ZWZ,j,k)+= designG[k*(*ng)+c]*designG[j*(*ng)+c]; 
	  }

	  if (time==stop[c] && status[c]==1)
	    {pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;
	    for(k=0;k<*pg;k++) VE(zi,k)=designG[k*(*ng)+c]; 
	    for(k=0;k<*px;k++) VE(xi,k)=designX[k*(*ng)+c]; 
	    }
	  count=count+1; 
      } } // }}}
      } else {
//	      Rprintf("=============== %d \n",count); 
      for (c=count;sstop==0 && c>=0;c=c-1) { // {{{
// Rprintf("times %d %lf  %lf %lf %d %d %d \n",s,time,start[c],stop[c],c,sstop,count); 
	if ((start[c]<time) && (stop[c]>=time)) {
	// Rprintf("under risk %d %d %d \n",c,id[c],count); 
	    for(j=0;j<pmax;j++)  {
	    if (j<*px) { 
	    for(k=0;k<*px;k++) 
		    ME(A,j,k)+= designX[j*(*nx)+c]*designX[k*(*nx)+c];
	    for(k=0;k<*pg;k++) 
		    ME(XWZ,j,k)+= designX[j*(*ng)+c]*designG[k*(*ng)+c]; 
	    }
	    if (j<*pg) 
	    for(k=0;k<*pg;k++) 
		    ME(ZWZ,j,k)+= designG[k*(*ng)+c]*designG[j*(*ng)+c]; 
	  }

	  if (time==stop[c] && status[c]==1)
	    {pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;
	    for(k=0;k<*pg;k++) VE(zi,k)=designG[k*(*ng)+c]; 
	    for(k=0;k<*px;k++) VE(xi,k)=designX[k*(*ng)+c]; 
	    }
	   // count=count-1; 

      } else {sstop=1; count=c; } } // }}}
      
      }


      //Rprintf(" s er %d \n",s); 
      //print_mat(A); print_mat(ZWZ); print_mat(XWZ); 

     // MtA(X,X,A); 
      invertS(A,AI,silent[0]); 
      if (ME(AI,0,0)==0 && *silent==0) Rprintf("time %lf X'X singular \n",time); 
     // MtA(Z,Z,ZWZ);MtA(X,Z,XWZ);
      MxA(AI,XWZ,XWZAI);

      MtA(XWZAI,XWZ,tmpM2); mat_subtr(ZWZ,tmpM2,dCGam); 
      scl_mat_mult(dtime,dCGam,dCGam); 
      if (*deltaweight==0) scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

      if (stat==1) {
   // 	extract_row(X,pers,tmpv1); 
        Mv(AI,xi,AIXWdN); 
   //	extract_row(Z,pers,zi); 
	vM(XWZ,AIXWdN,tmpv2);
	vec_subtr(zi,tmpv2,ZHdN);
	if (*deltaweight==0) scl_vec_mult(dtime,ZHdN,ZHdN); 
	vec_add(ZHdN,IZHdN,IZHdN); 
      }

      // scl_mat_mult(dtime,XWZAI,tmpM4);mat_add(tmpM4,Ct,Ct); 
    } /* s =1...Ntimes */ 

    // invertS(CGam,ICGam,silent[0]); Mv(ICGam,IZHdN,gam); 
    //if (ME(ICGam,0,0)==0 && *silent==0) Rprintf(" intZHZ  singular\n"); 
    //  print_mat(CGam); print_vec(IZHdN); 

    for(k=0;k<*pg;k++)  {
       intZHdN[k]=VE(IZHdN,k); 
       for(j=0;j<*pg;j++) intZHZ[k*(*pg)+j]=ME(CGam,k,j); 
    }


  free_mats(&X,&A,&AI,&AIXW,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam, 
		  &Ct,&XWZ,&XWZAI, &tmpM4,NULL); 
  free_vecs(&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&zi,&tmpv2,&rowZ,&gam,
		  &ZHdN,&IZHdN,&PLScomp,&Xi,NULL);
  free(ipers); free(ls); 
} // }}}

void compSsfix(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,id,status,deltaweight,intZHZ,intZHdN,silent)
double *designX,*alltimes,*start,*stop,*intZHZ,*intZHdN,*designG;
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*deltaweight,*id,*silent;
{ // {{{
  matrix *X,*A,*AI,*AIXW,*dCGam,*CGam,*Ct,*ICGam,*XWZ,*ZWZ,*XWZAI,*tmpM4,*tmpM2;
  vector *xi,*tmpv2,*tmpv1,*PLScomp,*Xi,*dA,*rowX,*AIXWdN,*korG,*rowZ,*gam,*ZHdN,
    *IZHdN,*zi;
  int sstop,j,k,l,c,s,count,pers=0,pmax,*ipers=calloc(*Ntimes,sizeof(int)); 
  int stat,*ls=calloc(*Ntimes,sizeof(int)); 
  double time,dtime,fabs(),sqrt();

  malloc_mats(*antpers,*px,&X,NULL);
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*px,*antpers,&AIXW,NULL);
  // malloc_mats(*antpers,*pg,&Z,NULL); 
  malloc_mats(*pg,*pg,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*pg,&Ct,&XWZ,&XWZAI,NULL);
  malloc_mat(*px,*pg,tmpM4); 

  malloc_vecs(*px,&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,NULL);
  malloc_vecs(*pg,&zi,&tmpv2,&rowZ,&gam,&ZHdN,&IZHdN,NULL);
  malloc_vecs(*antpers,&PLScomp,&Xi,NULL);

  if (*px>=*pg) pmax=*px; else pmax=*pg; 

    //mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN);  
    //mat_zeros(A); mat_zeros(ZWZ); mat_zeros(XWZ); 

    count=nx[0]-1; 
    for (s=(*Nalltimes)-1;s>0;s=s-1){
      sstop=0; 
    // Rprintf(" test %d %d %d  \n",s,*antpers,*nx); 
      time=alltimes[s]; 
      dtime=time-alltimes[s-1]; 
      stat=0;  

      l=0; stat=0;  
      if (1==0) {
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) { // {{{
//		Rprintf("times  %lf  %lf %lf \n",time,start[c],stop[c]); 
	if ((start[c]<time) && (stop[c]>=time)) {
//		Rprintf("under risk %d %d %d \n",c,id[c],count); 
	    for(j=0;j<pmax;j++)  {
	    if (j<*px) { 
	    for(k=0;k<*px;k++) 
		    ME(A,j,k)+= designX[j*(*nx)+c]*designX[k*(*nx)+c];
	    for(k=0;k<*pg;k++) 
		    ME(XWZ,j,k)+= designX[j*(*ng)+c]*designG[k*(*ng)+c]; 
	    }
	    if (j<*pg) 
	    for(k=0;k<*pg;k++) 
		    ME(ZWZ,j,k)+= designG[k*(*ng)+c]*designG[j*(*ng)+c]; 
	  }

	  if (time==stop[c] && status[c]==1)
	    {pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;
	    for(k=0;k<*pg;k++) VE(zi,k)=designG[k*(*ng)+c]; 
	    for(k=0;k<*px;k++) VE(xi,k)=designX[k*(*ng)+c]; 
	    }
	  count=count+1; 
      } } // }}}
      } else {
//	      Rprintf("=============== %d \n",count); 
      for (c=count;sstop==0 && c>=0;c=c-1) { // {{{
// Rprintf("times %d %lf  %lf %lf %d %d %d \n",s,time,start[c],stop[c],c,sstop,count); 
	if ((start[c]<time) && (stop[c]>=time)) {
	// Rprintf("under risk %d %d %d \n",c,id[c],count); 
	    for(j=0;j<pmax;j++)  {
	    if (j<*px) { 
	    for(k=0;k<*px;k++) 
		    ME(A,j,k)+= designX[j*(*nx)+c]*designX[k*(*nx)+c];
	    for(k=0;k<*pg;k++) 
		    ME(XWZ,j,k)+= designX[j*(*ng)+c]*designG[k*(*ng)+c]; 
	    }
	    if (j<*pg) 
	    for(k=0;k<*pg;k++) 
		    ME(ZWZ,j,k)+= designG[k*(*ng)+c]*designG[j*(*ng)+c]; 
	  }

	  if (time==stop[c] && status[c]==1)
	    {pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;
	    for(k=0;k<*pg;k++) VE(zi,k)=designG[k*(*ng)+c]; 
	    for(k=0;k<*px;k++) VE(xi,k)=designX[k*(*ng)+c]; 
	    }
	   // count=count-1; 

      } else {sstop=1; count=c; } } // }}}
      
      }


      //Rprintf(" s er %d \n",s); 
      //print_mat(A); print_mat(ZWZ); print_mat(XWZ); 

     // MtA(X,X,A); 
      invertS(A,AI,silent[0]); 
      if (ME(AI,0,0)==0 && *silent==0) Rprintf("time %lf X'X singular \n",time); 
     // MtA(Z,Z,ZWZ);MtA(X,Z,XWZ);
      MxA(AI,XWZ,XWZAI);

      MtA(XWZAI,XWZ,tmpM2); mat_subtr(ZWZ,tmpM2,dCGam); 
      scl_mat_mult(dtime,dCGam,dCGam); 
      if (*deltaweight==0) scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

      if (stat==1) {
   // 	extract_row(X,pers,tmpv1); 
        Mv(AI,xi,AIXWdN); 
   //	extract_row(Z,pers,zi); 
	vM(XWZ,AIXWdN,tmpv2);
	vec_subtr(zi,tmpv2,ZHdN);
	if (*deltaweight==0) scl_vec_mult(dtime,ZHdN,ZHdN); 
	vec_add(ZHdN,IZHdN,IZHdN); 
      }

      // scl_mat_mult(dtime,XWZAI,tmpM4);mat_add(tmpM4,Ct,Ct); 
    } /* s =1...Ntimes */ 

    // invertS(CGam,ICGam,silent[0]); Mv(ICGam,IZHdN,gam); 
    //if (ME(ICGam,0,0)==0 && *silent==0) Rprintf(" intZHZ  singular\n"); 
    //  print_mat(CGam); print_vec(IZHdN); 

    for(k=0;k<*pg;k++)  {
       intZHdN[k]=VE(IZHdN,k); 
       for(j=0;j<*pg;j++) intZHZ[k*(*pg)+j]=ME(CGam,k,j); 
    }


 free_mats(&X,&A,&AI,&AIXW,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam, 
		  &Ct,&XWZ,&XWZAI, &tmpM4,NULL); 
  free_vecs(&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&zi,&tmpv2,&rowZ,&gam,
		  &ZHdN,&IZHdN,&PLScomp,&Xi,NULL);
  free(ipers); free(ls); 
} // }}}


void compSsforward(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,id,status,deltaweight,intZHZ,intZHdN,silent)
double *designX,*alltimes,*start,*stop,*intZHZ,*intZHdN,*designG;
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*deltaweight,*id,*silent;
{ // {{{
 // {{{ allocating
  matrix *X,*A,*AI,*AIXW,*dCGam,*CGam,*Ct,*ICGam,*XWZ,*ZWZ,*XWZAI,*tmpM4,*tmpM2;
  vector *xi,*tmpv2,*tmpv1,*PLScomp,*Xi,*dA,*rowX,*AIXWdN,*korG,*rowZ,*gam,*ZHdN,
    *IZHdN,*zi;
  int sstop,j,k,l,c,s,count,pers=0,pmax,*ipers=calloc(*Ntimes,sizeof(int)); 
  int stat,*ls=calloc(*Ntimes,sizeof(int)); 
  double time,dtime,fabs(),sqrt();

  malloc_mats(*antpers,*px,&X,NULL);
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*px,*antpers,&AIXW,NULL);
  // malloc_mats(*antpers,*pg,&Z,NULL); 
  malloc_mats(*pg,*pg,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*pg,&Ct,&XWZ,&XWZAI,NULL);
  malloc_mat(*px,*pg,tmpM4); 

  malloc_vecs(*px,&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,NULL);
  malloc_vecs(*pg,&zi,&tmpv2,&rowZ,&gam,&ZHdN,&IZHdN,NULL);
  malloc_vecs(*antpers,&PLScomp,&Xi,NULL);
  // }}}

  if (*px>=*pg) pmax=*px; else pmax=*pg; 

    //mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN);  
    //mat_zeros(A); mat_zeros(ZWZ); mat_zeros(XWZ); 

    count=0; sstop=0; 
    for (s=1;s<*Nalltimes;s++){
    // Rprintf(" test %d %d %d  \n",s,*antpers,*nx); 
      time=alltimes[s]; dtime=time-alltimes[s-1]; stat=0;  

      l=0; stat=0;   sstop=0; 
      if (s==1) {
      for (c=0;c<*nx;c++) { // {{{
	if ((start[c]<time) && (stop[c]>=time)) {
	    for(j=0;j<pmax;j++)  {
	    if (j<*px) { 
	    for(k=0;k<*px;k++) ME(A,j,k)+= designX[j*(*nx)+c]*designX[k*(*nx)+c];
	    for(k=0;k<*pg;k++) ME(XWZ,j,k)+= designX[j*(*ng)+c]*designG[k*(*ng)+c]; 
	    }
	    if (j<*pg) 
	    for(k=0;k<*pg;k++) ME(ZWZ,j,k)+= designG[k*(*ng)+c]*designG[j*(*ng)+c]; 
	  }

	  if (time==stop[c] && status[c]==1)
	    {pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;
	    for(k=0;k<*pg;k++) VE(zi,k)=designG[k*(*ng)+c]; 
	    for(k=0;k<*px;k++) VE(xi,k)=designX[k*(*ng)+c]; 
	    }
      } } // }}}
      } else {
//	      Rprintf("=============== %d \n",count); 
      for (c=count;sstop==0 && c<*nx;c++) { // {{{
// Rprintf("times %d %lf  %lf %lf %d %d %d \n",s,time,start[c],stop[c],c,sstop,count); 
	if ((start[c]>time) || (stop[c]<time)) {
//	 Rprintf("not risk %d %d %d \n",c,id[c],count); 
	    for(j=0;j<pmax;j++)  {
	    if (j<*px) { 
	    for(k=0;k<*px;k++) 
		    ME(A,j,k)-= designX[j*(*nx)+c]*designX[k*(*nx)+c];
	    for(k=0;k<*pg;k++) 
		    ME(XWZ,j,k)-= designX[j*(*ng)+c]*designG[k*(*ng)+c]; 
	    }
	    if (j<*pg) 
	    for(k=0;k<*pg;k++) 
		    ME(ZWZ,j,k)-= designG[k*(*ng)+c]*designG[j*(*ng)+c]; 
	  }

      } else {sstop=1; count=c; } } // }}}
      sstop=0; 

      for (c=count;sstop==0 && c<*nx;c++) { // {{{

      if (time==stop[c] && status[c]==1) 
      {     pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;
	    for(k=0;k<*pg;k++) VE(zi,k)=designG[k*(*ng)+c]; 
	    for(k=0;k<*px;k++) VE(xi,k)=designX[k*(*ng)+c]; 
      } else {sstop=1;}
      if (stop[c]>time) sstop=1; 
      }

      }

      // Rprintf(" s er %d \n",s); 
      //print_mat(A); print_mat(ZWZ); print_mat(XWZ); 
     // MtA(X,X,A); 
     
      invertS(A,AI,silent[0]); 
      if (ME(AI,0,0)==0 && *silent==0) Rprintf("time %lf X'X singular \n",time); 
     // MtA(Z,Z,ZWZ);MtA(X,Z,XWZ);
      MxA(AI,XWZ,XWZAI);

      MtA(XWZAI,XWZ,tmpM2); mat_subtr(ZWZ,tmpM2,dCGam); 
      scl_mat_mult(dtime,dCGam,dCGam); 
      if (*deltaweight==0) scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

      if (stat==1) {
   // 	extract_row(X,pers,tmpv1); 
        Mv(AI,xi,AIXWdN); 
   //	extract_row(Z,pers,zi); 
	vM(XWZ,AIXWdN,tmpv2);
	vec_subtr(zi,tmpv2,ZHdN);
	if (*deltaweight==0) scl_vec_mult(dtime,ZHdN,ZHdN); 
	vec_add(ZHdN,IZHdN,IZHdN); 
      }

      // scl_mat_mult(dtime,XWZAI,tmpM4);mat_add(tmpM4,Ct,Ct); 
    } /* s =1...Ntimes */ 

    // invertS(CGam,ICGam,silent[0]); Mv(ICGam,IZHdN,gam); 
    //if (ME(ICGam,0,0)==0 && *silent==0) Rprintf(" intZHZ  singular\n"); 
    //  print_mat(CGam); print_vec(IZHdN); 

    for(k=0;k<*pg;k++)  {
       intZHdN[k]=VE(IZHdN,k); 
       for(j=0;j<*pg;j++) intZHZ[k*(*pg)+j]=ME(CGam,k,j); 
    }

  free_mats(&X,&A,&AI,&AIXW,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam, 
		  &Ct,&XWZ,&XWZAI, &tmpM4,NULL); 
  free_vecs(&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,&zi,&tmpv2,&rowZ,&gam,
		  &ZHdN,&IZHdN,&PLScomp,&Xi,NULL);
  free(ipers); free(ls); 
} // }}}


