#include <stdio.h>
#include <math.h>
#include "matrix.h"

void plssemiadd(alltimes,Nalltimes,Ntimes,designX,nx,px,designG,ng,pg,antpers,start,stop,id,cu,status,deltaweight,
plscov,dimplscov,betapls,plscomp,semipls,weighted,silent)
double *designX,*alltimes,*start,*stop,*cu,*designG,*plscov,*betapls,*plscomp; 
int *nx,*px,*antpers,*Nalltimes,*Ntimes,*ng,*pg,*status,*weighted,*deltaweight,*dimplscov,*semipls,*id,*silent;
{
  matrix *X,*A,*AI,*AIXW,*Z,*dCGam,*CGam,*Ct,*ICGam,*XWZ,*ZWZ,*XWZAI,*C[*Nalltimes],*tmpM4,*tmpM2;
  vector *xi,*tmpv2,*tmpv1,*PLScomp,*Xi,*dA,*rowX,*AIXWdN,*korG,*rowZ,*gam,*ZHdN,
    *IZHdN,*zi;
  int i,j,k,l,c,s,count,pers=0,pmax,*ipers=calloc(*Ntimes,sizeof(int)); 
  int stat,maxtime,*ls=calloc(*Ntimes,sizeof(int)),pls; 
  double time,dtime,fabs(),sqrt();

  if (*semipls==0) px[0]=px[0]+1; 
  malloc_mats(*antpers,*px,&X,NULL);
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*px,*antpers,&AIXW,NULL);
  malloc_mats(*antpers,*pg,&Z,NULL); 
  malloc_mats(*pg,*pg,&tmpM2,&ZWZ,&ICGam,&CGam,&dCGam,NULL); 
  malloc_mats(*px,*pg,&Ct,&XWZ,&XWZAI,NULL);
  malloc_mat(*px,*pg,tmpM4); 
  for (j=0;j<*Nalltimes;j++) malloc_mat(*px,*pg,C[j]);

  malloc_vecs(*px,&dA,&xi,&tmpv1,&korG,&rowX,&AIXWdN,NULL);
  malloc_vecs(*pg,&zi,&tmpv2,&rowZ,&gam,&ZHdN,&IZHdN,NULL);
  malloc_vecs(*antpers,&PLScomp,&Xi,NULL);

  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  maxtime=alltimes[*Nalltimes]; 

  for (pls=0;pls<*dimplscov;pls++){
    mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZHdN);  

    for (s=1;s<*Nalltimes;s++){
      time=alltimes[s]; dtime=time-alltimes[s-1]; mat_zeros(X); mat_zeros(Z); stat=0;  

      l=0; stat=0; 
      for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<pmax;j++) {
	    if (j<*px) ME(X,id[c],j)=designX[j*(*nx)+c];
	    if (j<(*pg)-1) ME(Z,id[c],j)= designG[j*(*ng)+c]; }
	  if (s==1) VE(Xi,id[c])=plscov[pls*(*ng)+c];
	  if (*semipls==1) {ME(Z,id[c],(*pg)-1)=plscov[pls*(*ng)+c];}
	  if (*semipls==0) ME(X,id[c],(*px)-1)=plscov[pls*(*ng)+c];
	  if (time==stop[c] && status[c]==1)
	    {pers=id[c];stat=1;l=l+1;ipers[l]=pers; ls[l]=s;}
	  count=count+1; }
      }

      if (s==0) {
	print_vec(Xi); 
	for (k=0;k<500;k++) {
	  for (c=0;c<*pg;c++) printf(" %lf ",ME(Z,k,c));  printf("\n"); }
	if (*pg==0) for (k=0;k<0;k++) printf(" %lf \n",ME(Z,k,0)); }

      MtA(X,X,A); invertS(A,AI,silent[0]); 
      if (ME(AI,0,0)==0 && *silent==0) printf("time %lf X'X singular \n",time); 
      MtA(Z,Z,ZWZ);MtA(X,Z,XWZ);MxA(AI,XWZ,XWZAI);

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
      mat_zeros(XWZ); mat_add(Ct,XWZ,C[s]); 
    } /* s =1...Ntimes */ 

    invertS(CGam,ICGam,silent[0]); Mv(ICGam,IZHdN,gam); 
    if (ME(ICGam,0,0)==0 && *silent==0) printf(" intZHZ  singular\n"); 

    if (*semipls==1) {
      if (*weighted==0) betapls[pls]=VE(gam,*pg-1);
      else betapls[pls]=VE(IZHdN,*pg-1);
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

  free_mats(&tmpM2,&X,&Z,&A,&AI,&ZWZ,&ICGam,&CGam,&dCGam,&Ct,
		&AIXW,&XWZ,&XWZAI,&tmpM4,NULL); // removed &C from the list
  free_vecs(&PLScomp,&Xi,&dA,&tmpv1,&tmpv2,&korG,&rowX,&AIXWdN,&zi,&rowZ,&gam,&ZHdN,&IZHdN,NULL);

  for (j=0;j<*Nalltimes;j++) free_mat(C[j]);

  free(ipers); free(ls); 
}
