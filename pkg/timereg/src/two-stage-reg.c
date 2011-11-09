#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
/* ====================================================== */
void twostagereg(times,Ntimes,designX,nx,px,designG,ng,pg,
antpers,start,stop, Nit,
detail, id,status, ratesim, robust, clusters,
antclust,betafixed, theta, vartheta,thetascore, inverse,
clustsize,desthetaI, ptheta,SthetaI,step,idiclust,notaylor,gamiid,biid,semi,cumhaz,
cumhazleft,lefttrunk,rr,maxtimesim,timegroup)
double *designX,*designG,*times,*start,*stop,*theta,*vartheta,*thetascore,*desthetaI,*SthetaI,*step,*gamiid,*biid,*cumhaz,*cumhazleft,*rr;
int *nx,*px,*ng,*pg,*antpers,*Ntimes,*Nit,*detail,*id,*status,*ratesim,*robust,
*clusters,*antclust,*betafixed,*inverse,*clustsize,*ptheta,*idiclust,*notaylor,*semi,*lefttrunk,*maxtimesim,*timegroup;
{
// {{{ defining variables
  matrix *ldesG0,*cdesX2;
  matrix *destheta,*d2UItheta,*d2Utheta,*varthetascore,*Sthetaiid[*antclust],*Stheta; 
  matrix *Ftilde,*Gtilde;
  matrix *Biid[*antclust]; 
  vector *lamt,*lamtt,*offset,*weight,*one;
  vector *tmpv1,*xi,*zi,*reszpbeta,*res1dim; 
  vector *vthetascore,*vtheta1,*vtheta2,*dtheta,*thetaiid[*antclust]; 
  vector *dbetaNH[*antclust],*dANH[*antclust]; 
  vector *dAiid[*antclust],*gammaiid[*antclust]; 
  int c,i,j,k,l,s,it,pmax,v,
      *cluster=calloc(*antpers,sizeof(int));
  double *Nt=calloc(*antclust,sizeof(double)),dummy,ll,lle;
  double tau,sumscore=999,*Nti=calloc(*antpers,sizeof(double)), theta0=0;
  double *thetaiidscale=calloc(*antclust,sizeof(double)),
   *NH=calloc(*antclust,sizeof(double)), 
   *HeH=calloc(*antclust,sizeof(double)),
   *HeHleft=calloc(*antclust,sizeof(double)),
   *H2eH=calloc(*antclust,sizeof(double)), *H2eHleft=calloc(*antclust,sizeof(double)), 
   *Rtheta=calloc(*antclust,sizeof(double)), *Rthetaleft=calloc(*antclust,sizeof(double)),
   *Hik=calloc(*antpers,sizeof(double)), Dthetanu=1; 
  int *ipers=calloc(*Ntimes,sizeof(int));

  for (j=0;j<*antclust;j++) { Nt[j]=0; NH[j]=0; 
   if (*notaylor==0) { malloc_vec(*pg,gammaiid[j]); 
	   malloc_mat(*maxtimesim,*px,Biid[j]);}
    malloc_vec(*pg,dbetaNH[j]); malloc_vec(*px,dANH[j]);  malloc_vec(*ptheta,dAiid[j]); 
    malloc_vec(*ptheta,thetaiid[j]); malloc_mat(*ptheta,*ptheta,Sthetaiid[j]); 
  }

  malloc_mat(*ptheta,*px,Ftilde); malloc_mat(*ptheta,*pg,Gtilde); 
  malloc_mats(*antpers,*px,&cdesX2,NULL); malloc_mats(*antpers,*pg,&ldesG0,NULL); 
  malloc_mat(*antpers,*ptheta,destheta); 
  malloc_mat(*ptheta,*ptheta,d2Utheta); malloc_mat(*ptheta,*ptheta,d2UItheta); 
  malloc_mat(*ptheta,*ptheta,Stheta); 
  malloc_mat(*ptheta,*ptheta,varthetascore); 
  malloc_vecs(*ptheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,NULL);
  malloc_vecs(*px,&tmpv1,&xi,NULL);
  malloc_vecs(*antpers,&weight,&lamtt,&lamt,&one,&offset,NULL); 
  malloc_vecs(*pg,&zi,NULL); 
  malloc_vec(1,reszpbeta); malloc_vec(1,res1dim); 

  if (*px>=*pg) pmax=*px; else pmax=*pg; ll=0; 
  
  for (i=0;i<*antpers;i++) { Hik[i]=0; Nti[i]=0; 
                             VE(one,i)=1; VE(weight,i)=1; VE(offset,i)=1;} 

  for (c=0;c<*antpers;c++) cluster[id[c]]=clusters[c];

  for (j=0;j<*antpers;j++)  { 
	  Hik[j]=cumhaz[j]; Nti[j]=Nti[j]+status[j]; Nt[cluster[j]]= Nt[cluster[j]]+status[j]; 
  }

//  for (j=0;j<*antpers;j++)  Rprintf("%d %lf %lf %lf  \n",j,cumhaz[j],Nti[j],cumhazleft[j]); 

  if (*notaylor==0) {
  for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++) VE(gammaiid[i],j)=gamiid[j*(*antclust)+i]; 

   for (i=0;i<*antclust;i++) for (s=0;s<*maxtimesim;s++) 
   for (c=0;c<*px;c++) {l=i*(*px)+c; ME(Biid[i],s,c)=biid[l*(*maxtimesim)+s]; }
    }
 

  // assuming that destheta is on this form: antpers x ptheta
  for(i=0;i<*antpers;i++) 
  for(j=0;j<*ptheta;j++) ME(destheta,i,j)=desthetaI[j*(*antpers)+i];

for (c=0;((c<*antpers));c++) 
for(j=0;j<pmax;j++) {
   if (j<*px) ME(cdesX2,c,j)=designX[j*(*ng)+c];  
   if (j<*pg) ME(ldesG0,c,j)=designG[j*(*ng)+c]; 
} 
  // }}}
 
   R_CheckUserInterrupt();

  for (i=0;i<*antpers;i++) 
  {
      j=cluster[i]; 
      NH[j]=NH[j]+Nti[i]*Hik[i]; 
      extract_row(ldesG0,i,zi); extract_row(cdesX2,i,xi);
      vec_add_mult(dbetaNH[j],zi,Nti[i]*Hik[i],dbetaNH[j]);  
      vec_add_mult(dANH[j],xi,Nti[i]*Hik[i],dANH[j]);  
      if (*betafixed==1) vec_zeros(dbetaNH[j]);  
  }

   R_CheckUserInterrupt();

//  for (j=0;j<*antpers;j++)  Rprintf("%d %lf %lf %lf  \n",j,Hik[j],Nti[j],cumhazleft[i]); 
//  Rprintf(" test =============\n");  
//  for (j=0;j<*antclust;j++)  Rprintf("%d %lf %lf \n",j,Nt[j],NH[j]); 
//  Rprintf(" test =============\n");  
//
/*===================Estimates theta, two stage approach of glidden ==== */


  for (i=0;i<*ptheta;i++) VE(vtheta1,i)=theta[i];  // starting values

  for (it=0;it<*Nit;it++) // {{{ frailty parameter Newton-Raphson
  {
   R_CheckUserInterrupt();
   for (j=0;j<*antclust;j++) {
	   Rthetaleft[j]=1; Rtheta[j]=1; HeHleft[j]=0; HeH[j]=0;H2eH[j]=0; H2eHleft[j]=0; }

      vec_zeros(vthetascore); mat_zeros(d2Utheta); 
      Mv(destheta,vtheta1,lamtt); 

      for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      for (k=0;k<clustsize[j];k++) {
	   i=idiclust[k*(*antclust)+j]; 
	    theta0=VE(lamtt,i); 
//	   Rprintf(" %d %d %d  %lf %lf %lf \n",j,k,i,Hik[i],theta0,Nt[j]); 
            if (*inverse==1) theta0=exp(theta0); 
//	    if (theta0<=0.00) Rprintf(" %lf %d %d %d \n",theta0,j,k,i); 
	    Rtheta[j]=Rtheta[j]+exp(theta0*Hik[i])-1; 
	    if (*lefttrunk==1) Rthetaleft[j]=Rthetaleft[j]+exp(theta0*cumhazleft[i])-1; 
	    HeH[j]=HeH[j]+Hik[i]*exp(theta0*Hik[i]); 
	    if (*lefttrunk==1) HeHleft[j]=HeHleft[j]+cumhazleft[i]*exp(theta0*cumhazleft[i]); 
	    H2eH[j]=H2eH[j]+pow(Hik[i],2)*exp(theta0*Hik[i]); 
	    if (*lefttrunk==1) H2eHleft[j]=H2eHleft[j]+pow(cumhazleft[i],2)*exp(theta0*cumhazleft[i]); 
      } 
      }

      for (j=0;j<*antclust;j++)  if (clustsize[j]>=2) {
      for (k=0;k<clustsize[j];k++) {
	      i=idiclust[k*(*antclust)+j]; 
              extract_row(destheta,i,vtheta2); theta0=VE(lamtt,i); 
              if (*inverse==1){theta0=exp(VE(lamtt,i));Dthetanu=theta0; 
	      }
//	      if (theta0<=0.00){  Rprintf("==== %lf %d %d %d \n",theta0,j,k,i); print_vec(vtheta2); }
     }

//      if (it==0) Rprintf(" %d %d %d %lf %lf %lf \n",j,k,i,Hik[i],theta0,Nt[j]); 

	 sumscore=0;  ll=0; 
	 if (Nt[j]>=2) 
	 for (k=2;k<=Nt[j];k++) {
	     tau=(k-1)/(1+theta0*(k-1));
	     lle=-pow((k-1),2)/pow((1+theta0*(k-1)),2);
	     sumscore=sumscore+tau; ll=ll+lle;
	 }

	 thetaiidscale[j]= sumscore+
		 log(Rtheta[j])/(theta0*theta0)-(1/theta0+Nt[j])*HeH[j]/Rtheta[j]+NH[j]; 
	 if (*lefttrunk==1) thetaiidscale[j]= thetaiidscale[j]-
		 log(Rthetaleft[j])/(theta0*theta0)+(1/theta0)*HeHleft[j]/Rthetaleft[j]; 

  scl_vec_mult(thetaiidscale[j]*Dthetanu,vtheta2,thetaiid[j]); 

  if (isnan(thetaiidscale[j])) {
  if (theta0<0) Rprintf("negative value of random effect variances causes problems, try step.size=0.1\n"); 
  Rprintf("nan i score subject=%d %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",j,theta0,Nt[j],NH[j],HeH[j],Rtheta[j],Rthetaleft[j],HeHleft[j],Dthetanu,thetaiidscale[j]); 
  print_vec(vtheta2); 
  oops("missing varlue\n"); 
  }
  if (isnan(thetaiidscale[j])) vec_zeros(thetaiid[j]); 

	  vec_add(vthetascore,thetaiid[j],vthetascore); 

	  tau=ll+(2/pow(theta0,2))*HeH[j]/Rtheta[j]
		-(2/pow(theta0,3))*log(Rtheta[j])
		-(1/theta0+Nt[j])*(H2eH[j]*Rtheta[j]-HeH[j]*HeH[j])/pow(Rtheta[j],2);
	 if (*lefttrunk==1) {
		 tau=tau-((2/pow(theta0,2))*HeHleft[j]/Rthetaleft[j]-
	        (2/pow(theta0,3))*log(Rthetaleft[j])-
	         (1/theta0)*(H2eHleft[j]*Rthetaleft[j]-HeHleft[j]*HeHleft[j])/pow(Rthetaleft[j],2));
	 }

	  for (i=0;i<*ptheta;i++) for (k=0;k<*ptheta;k++) {
	      ME(Sthetaiid[j],i,k)=VE(vtheta2,i)*VE(vtheta2,k)*tau*pow(Dthetanu,2);
	  }

	  mat_add(d2Utheta,Sthetaiid[j],d2Utheta); 
	}

      LevenbergMarquardt(d2Utheta,d2UItheta,vthetascore,dtheta,step,step);
// invert(d2Utheta,d2UItheta); 

      if (*detail==1) { 
	Rprintf("====================Iteration %d ==================== \n",it);
	Rprintf("Estimate theta \n"); print_vec(vtheta1); 
	Rprintf("Score D l\n");  print_vec(vthetascore); 
	Rprintf("Information D^2 l\n"); print_mat(d2UItheta); 
      }

//    Mv(d2UItheta,vthetascore,dtheta); scl_vec_mult(*step,dtheta,dthetaa); 
      vec_subtr(vtheta1,dtheta,vtheta1); 

      sumscore=0; 
      for (k=0;k<*pg;k++) sumscore= sumscore+fabs(VE(vthetascore,k));

      if ((sumscore<0.0000000001) & (it<*Nit-1)) it=*Nit-1; 
    } /* it theta Newton-Raphson */  // }}}

  for (i=0;i<*ptheta;i++) { theta[i]=VE(vtheta1,i);
    thetascore[i]=VE(vthetascore,i); }

   R_CheckUserInterrupt();
  /* terms for robust variances ============================ */
  if (*robust==1) { // {{{
    mat_zeros(Gtilde);   
    double dummyleft=0;
    for (k=0;k<*antclust;k++) if (clustsize[k]>=2) {

      if ((*notaylor==0) && (*semi==1)) {
      for (j=0;j<clustsize[k];j++) {
	   i= idiclust[j*(*antclust)+k]; 
           theta0=VE(lamtt,i); 
        if (*inverse==1) theta0=exp(theta0); 
	dummy= Hik[i]*(
		(1/(theta0*Rtheta[k]))*exp(theta0*Hik[i])-
	  (1/theta0+Nt[k])*(1+theta0*Hik[i])*exp(theta0*Hik[i])/Rtheta[k]+
	  Nti[i]+
	  (1+theta0*Nt[k])*exp(theta0*Hik[i])*HeH[k]/pow(Rtheta[k],2));
	if (*lefttrunk==1) {
		dummyleft=-cumhazleft[i]*(
		(1/(theta0*Rthetaleft[k]))*exp(theta0*cumhazleft[i])-
	  (1/theta0)*(1+theta0*cumhazleft[i])*exp(theta0*cumhazleft[i])/Rthetaleft[k]+
	  (1+theta0*Nt[k])*exp(theta0*cumhazleft[i])*HeHleft[k]/pow(Rthetaleft[k],2));
	}
        extract_row(ldesG0,i,zi); extract_row(destheta,i,vtheta1); 
        for (c=0;c<*ptheta;c++) for (l=0;l<*pg;l++) 
	  ME(Gtilde,c,l)=ME(Gtilde,c,l)+VE(zi,l)*VE(vtheta1,c)*(dummy+dummyleft);  
	}
      }

     if (*notaylor==0) 
      for (s=1;s<*Ntimes;s++) { // {{{
	if (k==0) {
	  mat_zeros(Ftilde); 
	  for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      for (v=0;v<clustsize[j];v++) {
	   i= idiclust[v*(*antclust)+j]; 
           if (stop[i]>=times[s]) {
	      theta0=VE(lamtt,i); 
              if (*inverse==1) theta0=exp(theta0); 
	      dummy=rr[i]* ((1/(theta0*Rtheta[j]))*exp(theta0*Hik[i])
		     -(1/theta0+Nt[j])*(1+theta0*Hik[i])*exp(theta0*Hik[i])/Rtheta[j]
		     +Nti[i]
		     +(1+theta0*Nt[j])*exp(theta0*Hik[i])*HeH[j]/pow(Rtheta[j],2));
    	      if (*lefttrunk==1) {
	       	dummyleft=-rr[i]*
	  ((1/(theta0*Rthetaleft[j]))*exp(theta0*cumhazleft[i])
	  -(1/theta0+Nt[j])*(1+theta0*cumhazleft[i])*exp(theta0*cumhazleft[i])/Rtheta[j]
	  +(1+theta0*Nt[j])*exp(theta0*cumhazleft[i])*HeHleft[j]/pow(Rthetaleft[j],2));
	      }
	      extract_row(destheta,i,vtheta1); extract_row(cdesX2,i,xi); 
            
	      for (c=0;c<*ptheta;c++) for (l=0;l<*px;l++) 
		ME(Ftilde,c,l)=ME(Ftilde,c,l)+ VE(xi,l)*VE(vtheta1,c)*(dummy+dummyleft);  
	    }
	  }
	}  
	}
	extract_row(Biid[k],timegroup[s],tmpv1); 
	extract_row(Biid[k],timegroup[s]-1,xi); 
	vec_subtr(tmpv1,xi,xi); Mv(Ftilde,xi,vtheta2); 
	vec_add(dAiid[k],vtheta2,dAiid[k]); 
      }  // }}}

    } /* k=1..antclust */ 


    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      if (*notaylor==0) {
      if (*semi==1) {
         Mv(Gtilde,gammaiid[j],vtheta2); 
         vec_add_mult(thetaiid[j],vtheta2,Dthetanu,thetaiid[j]);
      }
      vec_add_mult(thetaiid[j],dAiid[j],Dthetanu,thetaiid[j]);
      }

      for (i=0;i<*ptheta;i++) for (k=0;k<*ptheta;k++)
	ME(varthetascore,i,k) = ME(varthetascore,i,k) +
	  VE(thetaiid[j],i)*VE(thetaiid[j],k);
    }
  } // }}} robust==1

  MxA(d2UItheta,varthetascore,d2Utheta); 
  MxA(d2Utheta,d2UItheta,varthetascore);  

  for (j=0;j<*ptheta;j++) for (k=0;k<*ptheta;k++)
    {
      SthetaI[k*(*ptheta)+j]=ME(d2UItheta,j,k); 
      vartheta[k*(*ptheta)+j]=ME(varthetascore,j,k); 
    }

  // {{{ freeing everything
  
  for (j=0;j<*antclust;j++) {
  if (notaylor==0) { free_vec(gammaiid[j]); free_mat(Biid[j]); }
    free_vec(dbetaNH[j]); free_vec(dANH[j]);  free_vec(dAiid[j]); 
    free_vec(thetaiid[j]); free_mat(Sthetaiid[j]); }

  free_mats(&destheta, &d2Utheta, &d2UItheta, &Stheta, &Ftilde, &Gtilde,
	      &varthetascore, &cdesX2, &ldesG0,NULL); 

  free_vecs(&weight,&lamtt,&lamt, &one,&offset,&xi,&zi,&tmpv1,
	      &vthetascore,&vtheta1,&dtheta,&vtheta2,&reszpbeta, &res1dim,NULL); 

  free(Nt); free(Nti); free(thetaiidscale); free(NH); 
  free(HeH); free(H2eH); 
  free(HeHleft); free(H2eHleft); 
  free(Rtheta); free(Rthetaleft); free(Hik); free(cluster);  free(ipers); 
  // }}}
}
