#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"

void robaalenC(times,Ntimes,designX,nx,p,antpers,start,stop,cu,vcu,
	      robvcu,sim,antsim,retur,cumAit,test,rani,testOBS,status,
	      Ut,simUt,id,weighted,robust,covariance,covs,resample,
	      Biid,clusters,antclust,loglike,silent) 
double *designX,*times,*start,*stop,*cu,*vcu,*robvcu,*cumAit,*test,*testOBS,*Ut,*simUt,*covs,*Biid,*loglike; 
int *nx,*p,*antpers,*Ntimes,*sim,*retur,*rani,*antsim,*status,*id,*covariance,
    *weighted,*robust,*resample,*clusters,*antclust,*silent;
{ // {{{
  matrix *ldesignX, *QR, *R, *A, *AI, *Vcov;
  matrix *cumAt[*antclust];
  vector  *diag,*dB,*dN,*VdB,*xi,*rowX,*rowcum,*difX,*vtmp;
  vector *cumhatA[*antclust],*cumA[*antclust],*cum;
  int ci,i,j,k,l,s,c,count,pers=0,*cluster=calloc(*antpers,sizeof(int));
  double time,ahati,*vcudif=calloc((*Ntimes)*(*p+1),sizeof(double));
  double fabs(),sqrt();
  

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

  R_CheckUserInterrupt();

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
} // }}}

