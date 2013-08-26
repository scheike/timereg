//#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
/* ====================================================== */
void twostagereg(times,Ntimes,designX,nx,px,designG,ng,pg,
antpers,start,stop, Nit,
detail, id,status, ratesim, robust, clusters,
antclust,betafixed, theta, vartheta,thetascore, inverse,
clustsize,desthetaI, ptheta,SthetaI,step,idiclust,notaylor,gamiid,biid,semi,cumhaz,
cumhazleft,lefttrunk,rr,maxtimesim,timegroup,secluster,antsecluster,thetiid,timereso)
double *designX,*designG,*times,*start,*stop,*theta,*vartheta,*thetascore,*desthetaI,*SthetaI,*step,*gamiid,*biid,*cumhaz,*cumhazleft,*rr,*thetiid,*timereso;
int *nx,*px,*ng,*pg,*antpers,*Ntimes,*Nit,*detail,*id,*status,*ratesim,*robust,
*clusters,*antclust,*betafixed,*inverse,*clustsize,*ptheta,*idiclust,*notaylor,*semi,*lefttrunk,*maxtimesim,*timegroup,*antsecluster,*secluster;
{
// {{{ defining variables
  matrix *ldesG0,*cdesX2,*Ftilde,*Gtilde;
  matrix *destheta,*d2UItheta,*d2Utheta,*varthetascore,*Stheta,*mattheta; 
  matrix *Biid[*antsecluster],*Sthetaiid[*antsecluster];
  vector *dAiid[*antsecluster],*gammaiid[*antsecluster],*thetaiid[*antsecluster]; 
  vector *lamt,*lamtt,*offset,*weight,*one;
  vector *tmpv1,*xi,*zi,*reszpbeta,*res1dim; 
  vector *vthetascore,*vtheta3,*vtheta1,*vtheta2,*dtheta;
  int    cc,c,i,j,k,l,s,it,pmax,v;
  double dummy,ll,lle;
  double tau,sumscore=999, theta0=0,Dthetanu=1,logl=0; 
  double *thetaiidscale=calloc(*antclust,sizeof(double)),
	 *Nt=calloc(*antclust,sizeof(double)),
         *Nti=calloc(*antpers,sizeof(double)), 
         *NH=calloc(*antclust,sizeof(double)), 
         *HeH=calloc(*antclust,sizeof(double)),
         *HeHleft=calloc(*antclust,sizeof(double)),
         *H2eH=calloc(*antclust,sizeof(double)), 
	 *H2eHleft=calloc(*antclust,sizeof(double)), 
         *Rtheta=calloc(*antclust,sizeof(double)), 
	 *Rthetaleft=calloc(*antclust,sizeof(double)),
         *Hik=calloc(*antpers,sizeof(double)),
         *insecluster=calloc(*antsecluster,sizeof(double));
  int *ipers=calloc(*Ntimes,sizeof(int)),
      *cluster=calloc(*antpers,sizeof(int)),
      *multitrunc=calloc(*antclust,sizeof(int));

  for (j=0;j<*antclust;j++) { 
      Nt[j]=0; NH[j]=0; multitrunc[j]=0; 
   if ((j<*antsecluster)) {
	   insecluster[j]=0; 
      if ((*notaylor==0)) {
	malloc_vec(*pg,gammaiid[j]); 
	malloc_mat(*maxtimesim,*px,Biid[j]);
       }
       malloc_vec(*ptheta,dAiid[j]); 
       malloc_vec(*ptheta,thetaiid[j]); 
       malloc_mat(*ptheta,*ptheta,Sthetaiid[j]); 
   }
  }

  malloc_mat(*ptheta,*px,Ftilde); 
  malloc_mat(*ptheta,*pg,Gtilde); 
  malloc_mats(*antpers,*px,&cdesX2,NULL); 
  malloc_mats(*antpers,*pg,&ldesG0,NULL); 
  malloc_mat(*antpers,*ptheta,destheta); 
  malloc_mat(*ptheta,*ptheta,d2Utheta); 
  malloc_mat(*ptheta,*ptheta,d2UItheta); 
  malloc_mat(*ptheta,*ptheta,Stheta); 
  malloc_mat(*ptheta,*ptheta,mattheta); 
  malloc_mat(*ptheta,*ptheta,varthetascore); 

  malloc_vecs(*ptheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,&vtheta3,NULL);
  malloc_vecs(*px,&tmpv1,&xi,NULL);
  malloc_vecs(*antpers,&weight,&lamtt,&lamt,&one,&offset,NULL); 
  malloc_vecs(*pg,&zi,NULL); 
  malloc_vec(1,reszpbeta); 
  malloc_vec(1,res1dim); 

  if (*px>=*pg) pmax=*px; else pmax=*pg; ll=0; 
  
  for (i=0;i<*antpers;i++) { Hik[i]=0; Nti[i]=0; 
                             VE(one,i)=1; VE(weight,i)=1; VE(offset,i)=1;} 

  for (c=0;c<*antpers;c++) cluster[id[c]]=clusters[c];

  for (j=0;j<*antpers;j++)  { 
	  Hik[j]=cumhaz[j]; Nti[j]=Nti[j]+status[j]; Nt[cluster[j]]= Nt[cluster[j]]+status[j]; 
  }

// for (j=0;j<*antpers;j++) Rprintf("%d %d %d %lf %lf %lf  \n",j,id[j],cluster[id[j]],cumhaz[j],Nti[j],cumhazleft[j]); 

  if (*notaylor==0) {
if (*semi==1) for (i=0;i<*antsecluster;i++) 
	      for (j=0;j<*pg;j++) VE(gammaiid[i],j)=gamiid[j*(*antclust)+i]; 
//     if (*semi==1)  for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++) VE(gammaiid[i],j)=0;
//     for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++) 
//     printf("%d %d %d %lf \n",*antclust,j,i,gamiid[j*(*antclust)+i]); 
      for (i=0;i<*antsecluster;i++) for (s=0;s<*maxtimesim;s++) 
      for (c=0;c<*px;c++) {
	      l=i*(*px)+c; 
	      ME(Biid[i],s,c)=biid[l*(*maxtimesim)+s]; 
//	      ME(Biid[i],s,c)=0;
//	      printf("%d %d %d %d %d %d \n",*maxtimesim,l,s,j,i,l*(*maxtimesim)+s); 
      }
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
      if (start[i]>0) multitrunc[j]+=1; 
  }

   R_CheckUserInterrupt();
/*===================Estimates theta, two stage approach of glidden ==== */
  for (i=0;i<*ptheta;i++) VE(vtheta1,i)=theta[i];  // starting values

  for (it=0;it<*Nit;it++) // {{{ frailty parameter Newton-Raphson
  {
   R_CheckUserInterrupt();
   for (j=0;j<*antclust;j++) {
	   Rthetaleft[j]=1; Rtheta[j]=1; HeHleft[j]=0; HeH[j]=0;H2eH[j]=0; H2eHleft[j]=0; 
   }

  vec_zeros(vthetascore); mat_zeros(d2Utheta); 
  logl=0; 
  Mv(destheta,vtheta1,lamtt); 

      for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      for (k=0;k<clustsize[j];k++) {
	   i=idiclust[k*(*antclust)+j]; 
           theta0=VE(lamtt,i); 
//	   Rprintf(" %d %d %d  %lf %lf %lf \n",j,k,i,Hik[i],theta0,Nt[j]); 
//	   Rprintf(" %lf %lf %lf %lf \n",start[i],stop[i],cumhazleft[i],Hik[i]); 
//	   if (theta0<=0.00) Rprintf(" %lf %d %d %d \n",theta0,j,k,i); 
           if (*inverse==1) theta0=exp(theta0); 
	    Rtheta[j]=Rtheta[j]+exp(theta0*Hik[i])-1; 
	    HeH[j]=HeH[j]+Hik[i]*exp(theta0*Hik[i]); 
	    H2eH[j]=H2eH[j]+pow(Hik[i],2)*exp(theta0*Hik[i]); 
	    if (*lefttrunk==1 && (multitrunc[j]>=2) && (start[i]>0)) {
  	       Rthetaleft[j]=Rthetaleft[j]+exp(theta0*cumhazleft[i])-1; 
	       HeHleft[j]=HeHleft[j]+cumhazleft[i]*exp(theta0*cumhazleft[i]); 
	       H2eHleft[j]=H2eHleft[j]+pow(cumhazleft[i],2)*exp(theta0*cumhazleft[i]); 
	    }
      }
      }

  for (j=0;j<*antclust;j++)  if (clustsize[j]>=2) {
     if (it==0)  {
         cc=idiclust[0*(*antclust)+j]; // takes one index from cluster 
         cc=secluster[cc];             // secluster id related to this cluster
         insecluster[cc]=1;            // something going on in secluster
     }

     //  takes design and parameter to this cluster
      i=idiclust[0*(*antclust)+j];  // index from this cluster 
      extract_row(destheta,i,vtheta2); 
      theta0=VE(lamtt,i); 
      if (*inverse==1){theta0=exp(VE(lamtt,i));Dthetanu=theta0; }

 sumscore=0;  ll=0; 
 if (Nt[j]>=2) 
 for (k=2;k<=Nt[j];k++) {
     tau=(k-1)/(1+theta0*(k-1));
     lle=-pow((k-1),2)/pow((1+theta0*(k-1)),2);
     sumscore=sumscore+tau; 
     ll=ll+lle;
     logl+=log((1+theta0*(k-1)));
 }
 logl=logl+theta0*NH[j]+(1/theta0+Nt[j])*log(Rtheta[j]); 
 if (*lefttrunk==1 && multitrunc[j]>=2) logl=logl-(1/theta0)*log(Rthetaleft[j]); 

  if (it<0  && multitrunc[j]>=2) {
    Rprintf(" %d %d %lf %lf %lf \n",j,clustsize[j],Rtheta[j],HeH[j],H2eH[j]); 
    Rprintf(" %d %d %lf %lf %lf \n",j,multitrunc[j],Rthetaleft[j],HeHleft[j],H2eHleft[j]); 
  }

 thetaiidscale[j]=sumscore+
	          log(Rtheta[j])/(theta0*theta0)-(1/theta0+Nt[j])*HeH[j]/Rtheta[j]+NH[j]; 
// printf(" %d %d %lf \n",j,clustsize[j],thetaiidscale[j]); 
 if (*lefttrunk==1 && multitrunc[j]>=2) 
 thetaiidscale[j]=thetaiidscale[j]-
	          log(Rthetaleft[j])/(theta0*theta0)+(1/theta0)*HeHleft[j]/Rthetaleft[j]; 
// printf("%d %d %lf \n",j,multitrunc[j],thetaiidscale[j]); 

    scl_vec_mult(thetaiidscale[j]*Dthetanu,vtheta2,vtheta3); 

  if (isnan(thetaiidscale[j])) {
  if (theta0<0) Rprintf("negative value of random effect variances causes problems, try step.size=0.1\n"); 
  Rprintf("nan i score subject=%d %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",j,theta0,Nt[j],NH[j],HeH[j],Rtheta[j],Rthetaleft[j],HeHleft[j],Dthetanu,thetaiidscale[j]); 
  print_vec(vtheta3); 
  oops("missing value\n"); 
  }
  if (isnan(thetaiidscale[j])) vec_zeros(vtheta3); 

  if (it==(*Nit-1)) { 
	  vec_add(vtheta3,thetaiid[secluster[i]],thetaiid[secluster[i]]); 
  }
  vec_add(vthetascore,vtheta3,vthetascore); 

   tau=ll+(2/pow(theta0,2))*HeH[j]/Rtheta[j]-(2/pow(theta0,3))*log(Rtheta[j])
 	 -(1/theta0+Nt[j])*(H2eH[j]*Rtheta[j]-HeH[j]*HeH[j])/pow(Rtheta[j],2);
     if (*lefttrunk==1 && multitrunc[j]>=2) {
	 tau=tau-((2/pow(theta0,2))*HeHleft[j]/Rthetaleft[j]-(2/pow(theta0,3))*log(Rthetaleft[j])
  	        -(1/theta0)*(H2eHleft[j]*Rthetaleft[j]-HeHleft[j]*HeHleft[j])/pow(Rthetaleft[j],2));
     }

    for (c=0;c<*ptheta;c++) for (k=0;k<*ptheta;k++) {
//	      ME(Sthetaiid[j],c,k)=VE(vtheta2,c)*VE(vtheta2,k)*tau*pow(Dthetanu,2);
      if (*inverse==0) ME(mattheta,c,k)=VE(vtheta2,c)*VE(vtheta2,k)*tau*pow(Dthetanu,2);
      if (*inverse==1) 
      ME(mattheta,c,k)=VE(vtheta2,c)*VE(vtheta2,k)*(tau*pow(Dthetanu,2)+thetaiidscale[j]*Dthetanu); 
    }

    mat_add(d2Utheta,mattheta,d2Utheta); 
  } // j=1...antclust

//   LevenbergMarquardt(d2Utheta,d2UItheta,vthetascore,dtheta,step,step);
    invertS(d2Utheta,d2UItheta,1); 

      if (*detail==1) { 
	Rprintf("====================Iteration %d ==================== \n",it);
	Rprintf("Log-likelihood  %lf \n",logl); 
	Rprintf("Estimate theta \n"); print_vec(vtheta1); 
	Rprintf("Score D l\n");  print_vec(vthetascore); 
	Rprintf("Information D^2 l\n"); print_mat(d2UItheta); 
      }

    Mv(d2UItheta,vthetascore,dtheta); 
    scl_vec_mult(*step,dtheta,dtheta); 
    vec_subtr(vtheta1,dtheta,vtheta1); 

    sumscore=0; 
    for (k=0;k<*pg;k++) sumscore= sumscore+fabs(VE(vthetascore,k));

    if ((sumscore<0.0000001) & (it<*Nit-2)) it=*Nit-2; 
  } /* it theta Newton-Raphson */  // }}}

  if (*detail==1) Rprintf("Newton-Raphson ok \n"); 
  for (i=0;i<*ptheta;i++) { theta[i]=VE(vtheta1,i); thetascore[i]=VE(vthetascore,i); }

   R_CheckUserInterrupt();
  /* terms for robust variances ============================ */
  if (*robust==1) { // {{{  
    mat_zeros(Gtilde);   
    double dummyleft=0;
    int engang=0; 

      if ((*notaylor==0) && (*semi==1)) { // {{{  derivative D_betaU =  Gtilde
      for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      for (v=0;v<clustsize[j];v++) {
	   i=idiclust[v*(*antclust)+j]; 
	   theta0=VE(lamtt,i); 
	if (*inverse==1) theta0=exp(theta0); 
	dummy= Hik[i]*((1/(theta0*Rtheta[j]))*exp(theta0*Hik[i])
			-(1/theta0+Nt[j])*(1+theta0*Hik[i])*exp(theta0*Hik[i])/Rtheta[j]
			+Nti[i]+(1+theta0*Nt[j])*exp(theta0*Hik[i])*HeH[j]/pow(Rtheta[j],2));
	if (*lefttrunk==1 && multitrunc[j]>=2) {
		dummyleft=-cumhazleft[i]*(
		(1/(theta0*Rthetaleft[j]))*exp(theta0*cumhazleft[i])
		-(1/theta0)*(1+theta0*cumhazleft[i])*exp(theta0*cumhazleft[i])/Rthetaleft[j]
	        +exp(theta0*cumhazleft[i])*HeHleft[j]/pow(Rthetaleft[j],2)
		);
	} else dummyleft=0; 
	if (*notaylor==0) {
	   extract_row(ldesG0,i,zi); 
	   extract_row(destheta,i,vtheta1); 
	   for (c=0;c<*ptheta;c++) for (l=0;l<*pg;l++) 
	      ME(Gtilde,c,l)=ME(Gtilde,c,l)+VE(zi,l)*VE(vtheta1,c)*Dthetanu*(dummy+dummyleft);  
	}
	}
      }
      } // }}} 

    if (*notaylor==0) 
      for (s=1;s<*maxtimesim;s++) { // {{{
         engang=0; 
    for (k=0;k<*antsecluster;k++) if (insecluster[k]>0) { 
//    for (k=0;k<*antclust;k++) if (clustsize[k]>=2) {
//      cc=idiclust[0*(*antclust)+k]; // takes one index from cluster 
//      cc=secluster[cc];             // secluster id related to this cluster
//      insecluster[cc]=1; 

    if (engang==0) { // {{{ derivatve D_baseline(t)
        engang=1; 
        mat_zeros(Ftilde); 
      for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
      for (v=0;v<clustsize[j];v++) {
	   i=idiclust[v*(*antclust)+j]; 
           if ((stop[i]>=timereso[s])) {
	      theta0=VE(lamtt,i); 
              if (*inverse==1) theta0=exp(theta0); 
	      dummy=rr[i]*((1/(theta0*Rtheta[j]))*exp(theta0*Hik[i])
		     -(1/theta0+Nt[j])*(1+theta0*Hik[i])*exp(theta0*Hik[i])/Rtheta[j]+Nti[i]
		     +(1+theta0*Nt[j])*exp(theta0*Hik[i])*HeH[j]/pow(Rtheta[j],2));
	   } else dummy=0; 
//    	   if (*lefttrunk==1 && multitrunc[j]>=2 ) {
    	   if (*lefttrunk==1 && multitrunc[j]>=2 && (start[i]>=timereso[s])) {
	     dummyleft=-rr[i]*
	     ((1/(theta0*Rthetaleft[j]))*exp(theta0*cumhazleft[i])
	     -(1/theta0)*(1+theta0*cumhazleft[i])*exp(theta0*cumhazleft[i])/Rthetaleft[j] 
	      +exp(theta0*cumhazleft[i])*HeHleft[j]/pow(Rthetaleft[j],2));
	   } else dummyleft=0; 
	   extract_row(destheta,i,vtheta1); 
	   extract_row(cdesX2,i,xi); 
	   for (c=0;c<*ptheta;c++) for (l=0;l<*px;l++) 
	   ME(Ftilde,c,l)=ME(Ftilde,c,l)+VE(xi,l)*VE(vtheta1,c)*Dthetanu*(dummy+dummyleft);  
	   } 
       }
      } // }}} 
	extract_row(Biid[k],s,tmpv1); extract_row(Biid[k],s-1,xi); 
	vec_subtr(tmpv1,xi,xi); 
	Mv(Ftilde,xi,vtheta2); 
	vec_add(dAiid[k],vtheta2,dAiid[k]); 
      }  /* k=0..antseclust */ // }}}
    } // s=1, maxtimesim 

    for (j=0;j<*antsecluster;j++) { // {{{ iid 
//	    printf("-------------------- %d \n",j); 
      if (*notaylor==0) {
      if (insecluster[j]>0) 
      if (*semi==1) {
         Mv(Gtilde,gammaiid[j],vtheta2); 
         vec_add(thetaiid[j],vtheta2,thetaiid[j]);
//         print_vec(vtheta2); 
      }
       vec_add(thetaiid[j],dAiid[j],thetaiid[j]);
//      print_vec(dAiid[j]); 
      }
      for (i=0;i<*ptheta;i++) 
      {
          thetiid[j*(*ptheta)+i]=VE(thetaiid[j],i); 
          for (k=0;k<*ptheta;k++)  
	   ME(varthetascore,i,k)=ME(varthetascore,i,k)+VE(thetaiid[j],i)*VE(thetaiid[j],k);
      }
    } // }}} 

  } // }}} robust==1

  MxA(d2UItheta,varthetascore,d2Utheta); 
  MxA(d2Utheta,d2UItheta,varthetascore);  

  for (j=0;j<*ptheta;j++) for (k=0;k<*ptheta;k++)
  {
      SthetaI[k*(*ptheta)+j]=ME(d2UItheta,j,k); 
      vartheta[k*(*ptheta)+j]=ME(varthetascore,j,k); 
  }

  // {{{ freeing everything
  
  for (j=0;j<*antsecluster;j++) {
  if (*notaylor==0) { free_vec(gammaiid[j]); free_mat(Biid[j]); }
    free_vec(dAiid[j]); 
    free_vec(thetaiid[j]); 
    free_mat(Sthetaiid[j]); 
  }

  free_mats(&mattheta,&destheta, &d2Utheta, &d2UItheta, &Stheta, &Ftilde, &Gtilde,
	      &varthetascore, &cdesX2, &ldesG0,NULL); 

  free_vecs(&weight,&lamtt,&lamt, &one,&offset,&xi,&zi,&tmpv1,
	      &vthetascore,&vtheta3,&vtheta1,&dtheta,&vtheta2,&reszpbeta, &res1dim,NULL); 

  free(Nt); free(Nti); free(thetaiidscale); free(NH); free(HeH); free(H2eH); 
  free(HeHleft); free(H2eHleft); free(Rtheta); free(Rthetaleft); 
  free(Hik); free(cluster);  free(ipers); free(insecluster); free(multitrunc);  
  // }}}
}
