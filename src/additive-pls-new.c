//#include <stdio.h>
#include <math.h>
#include "matrix.h"

void plssemiaddN(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,id,cu,status,deltaweight,
plscov,dimplscov,betapls,plscomp,semipls,weighted,silent,weights,entry)
double *designX,*alltimes,*start,*stop,*cu,*designG,*plscov,*betapls,*plscomp,*weights; 
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*weighted,*deltaweight,*dimplscov,*semipls,*id,*silent,*entry;
{
// {{{ setup and allocating
  matrix *WX,*WZ,*X,*A,*AI,*AIXW,*Z,*dCGam,*CGam,*Ct,*ICGam,*XWZ,*ZWZ,*XWZAI,*C[*Nalltimes],*tmpM4,*tmpM2;
  vector *xi,*tmpv2,*tmpv1,*PLScomp,*Xi,*dA,*rowX,*AIXWdN,*korG,*rowZ,*gam,*ZHdN,
    *IZHdN,*zi;
  int ci=0,i,j,k,l,c,s,count,pers=0,pmax,*ipers=calloc(*Ntimes,sizeof(int)); 
  int stat,*ls=calloc(*Ntimes,sizeof(int)),pls; 
  double time,dtime,fabs(),sqrt();

  if (*semipls==0) px[0]=px[0]+1; 
  malloc_mats(*antpers,*px,&X,&WX,NULL);
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*px,*antpers,&AIXW,NULL);
  malloc_mats(*antpers,*pg,&Z,&WZ,NULL); 
  malloc_mats(*pg,*pg,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*pg,&Ct,&XWZ,&XWZAI,NULL);
  malloc_mat(*px,*pg,tmpM4); 
  for (j=0;j<*Nalltimes;j++) malloc_mat(*px,*pg,C[j]);

  malloc_vecs(*px,&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,NULL);
  malloc_vecs(*pg,&zi,&tmpv2,&rowZ,&gam,&ZHdN,&IZHdN,NULL);
  malloc_vecs(*antpers,&PLScomp,&Xi,NULL);

  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  // }}}

  for (pls=0;pls<*dimplscov;pls++){
    mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN);  
    mat_zeros(X); mat_zeros(WX); mat_zeros(Z); mat_zeros(WZ); 
    mat_zeros(A); mat_zeros(XWZ); mat_zeros(ZWZ); 
    l=0; 

    for (s=1;s<*Nalltimes;s++){
      time=alltimes[s]; dtime=time-alltimes[s-1]; stat=0;  

      vec_zeros(rowX); vec_zeros(rowZ); stat=0;   

  // {{{ reading design and making matrix products
   if (s==1)  { // {{{ reading start design 
     for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<pmax;j++) {
	  if (j<*px) { ME(X,id[c],j)=designX[j*(*nx)+c]; }
	  if (j<*px) { ME(WX,id[c],j) =designX[j*(*nx)+c]; }
	  if (j<*pg) { ME(Z,id[c],j)=designG[j*(*ng)+c]; }
	  if (j<*pg) { ME(WZ,id[c],j)=designG[j*(*ng)+c]; } 
	   VE(Xi,id[c])=plscov[pls*(*ng)+c];
	  if (*semipls==1) {ME(Z,id[c],(*pg)-1)=plscov[pls*(*ng)+c];}
	  if (*semipls==1) {ME(WZ,id[c],(*pg)-1)=plscov[pls*(*ng)+c];}
	  if (*semipls==0) ME(X,id[c],(*px)-1)=plscov[pls*(*ng)+c];
	  if (*semipls==0) ME(WX,id[c],(*px)-1)=plscov[pls*(*ng)+c];
	  }
	  if (time==stop[c] && status[c]==1) {
	    pers=id[c];stat=1;l=l+1; ls[l]=s;
	  }
	  count=count+1; 
	}
      }
    MtA(X,WX,A); MtA(Z, WZ,ZWZ); MtA(X,WZ,XWZ); 
    ci=*nx-1; 
    while ((stop[ci]<time)  & (ci>=0) )  ci=ci-1; 
    } // }}}

     vec_zeros(rowX); vec_zeros(rowZ); 
    if (s>1)  // {{{ modifying design for next time points
    while ((stop[ci]<time)  & (ci>=0) ) {
// Rprintf("ww %d %d  %lf %lf %d \n",ci,id[ci],stop[ci],time,entry[ci]); 
            for(j=0;j<*px;j++) VE(xi,j)=designX[j*(*nx)+ci]; 
            for(j=0;j<*pg;j++) VE(zi,j)=designG[j*(*nx)+ci]; 
	  if (*semipls==1) VE(zi,(*pg)-1)=plscov[pls*(*ng)+ci];
	  if (*semipls==0) VE(xi,(*px)-1)=plscov[pls*(*ng)+ci];
//            print_vec(xi); print_vec(zi); 
	    if (entry[ci]==1)  {
	           replace_row(X,id[ci],xi); replace_row(Z,id[ci],zi); 
	    } 
	    else {replace_row(X,id[ci],rowX);replace_row(Z,id[ci],rowZ);}
//	    Rprintf(" hej \n"); 
	  for(j=0;j<pmax;j++) 
	  for(k=0;k<pmax;k++)  {
//		  Rprintf(" %d %d %d \n",j,k,pmax); 
              if ((j<*px) & (k<*px)) ME(A,j,k)=
		      ME(A,j,k)+entry[ci]*VE(xi,k)*VE(xi,j)*weights[ci]; 
              if ((j<*px) & (k<*pg)) ME(XWZ,j,k)=
		      ME(XWZ,j,k)+entry[ci]*VE(zi,k)*VE(xi,j)*weights[ci]; 
              if ((j<*pg) & (k<*pg)) ME(ZWZ,j,k)=
		      ME(ZWZ,j,k)+entry[ci]*VE(zi,k)*VE(zi,j)*weights[ci]; 
	  }
	  ci=ci-1; 
    }

// Rprintf("ci ci ci  %lf %lf %ld \n",time,stop[ci],status[ci]); 
    if ((s>1) & (time==stop[ci]) & (status[ci]==1)) {
         pers=id[ci]; stat=1;l=l+1; ls[l]=s;
    }

    // }}}
 // }}}

      if (s==0) { print_vec(Xi); 
	for (k=0;k<500;k++) {
	  for (c=0;c<*pg;c++) Rprintf(" %lf ",ME(Z,k,c));  Rprintf("\n"); }
	if (*pg==0) for (k=0;k<0;k++) Rprintf(" %lf \n",ME(Z,k,0)); }

      invertS(A,AI,silent[0]); 
      if (ME(AI,0,0)==0 && *silent==0) Rprintf("time %lf X'X singular \n",time); 
      MxA(AI,XWZ,XWZAI);

      MtA(XWZAI,XWZ,tmpM2); mat_subtr(ZWZ,tmpM2,dCGam); 
      scl_mat_mult(dtime,dCGam,dCGam); 
      if (*deltaweight==0) scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

      if (stat==1) {
	extract_row(X,pers,tmpv1); Mv(AI,tmpv1,AIXWdN); 
	extract_row(Z,pers,zi); vM(XWZ,AIXWdN,tmpv2);
	vec_subtr(zi,tmpv2,ZHdN);
	if (*deltaweight==0) scl_vec_mult(dtime,ZHdN,ZHdN); 
	vec_add(ZHdN,IZHdN,IZHdN); 
      }

      scl_mat_mult(dtime,XWZAI,tmpM4);mat_add(tmpM4,Ct,Ct); 
      // mat_zeros(XWZ); 
       mat_copy(Ct,C[s]); 
    } /* s =1...Ntimes */ 

    invertS(CGam,ICGam,silent[0]); Mv(ICGam,IZHdN,gam); 
    if (ME(ICGam,0,0)==0 && *silent==0) Rprintf(" intZHZ  singular\n"); 

    if (*semipls==1) {
      if (*weighted==0) betapls[pls]=VE(gam,*pg-1); else betapls[pls]=VE(IZHdN,*pg-1);
      if (*weighted==0) vec_add_mult(PLScomp,Xi,betapls[pls],PLScomp);
      else vec_add_mult(PLScomp,Xi,VE(IZHdN,*pg-1),PLScomp); }

    if (*semipls==0) {
      for (s=1;s<*Ntimes;s++) {
	Mv(C[ls[s]],gam,korG);
	for (k=1;k<=*px;k++) cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s]-VE(korG,k-1);

	betapls[pls*(*Ntimes)+s]=cu[(*px)*(*Ntimes)+s];    

	for (i=0;i<*antpers;i++) plscomp[i*(*Ntimes)+s]=plscomp[i*(*Ntimes)+s]+ 
	  VE(Xi,i)*betapls[pls*(*Ntimes)+s]; }
    } /* semipls==0 */ 
  } /* dimcovpls comps */ 

  if (*semipls==1) for (i=0;i<*antpers;i++) plscomp[i]=VE(PLScomp,i); 
  /*====================================================== */

  free_mats(&WX,&WZ,&tmpM2,&X,&Z,&A,&AI,&ZWZ,&ICGam,&CGam,&dCGam,&Ct,
		&AIXW,&XWZ,&XWZAI,&tmpM4,NULL); // removed &C from the list
  free_vecs(&PLScomp,&Xi,&dA,&tmpv1,&tmpv2,&korG,&rowX,&AIXWdN,&zi,&rowZ,&gam,&ZHdN,&IZHdN,NULL);

  for (j=0;j<*Nalltimes;j++) free_mat(C[j]);

  free(ipers); free(ls); 
}
