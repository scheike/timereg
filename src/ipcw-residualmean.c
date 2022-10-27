#include <stdio.h>
#include <math.h>
#include "matrix.h"
	 
void resmean(double *times,int *Ntimes,double *x,int *delta,int *cause,double *KMc,double *z,int *n,int *px,int *Nit,double *betaS,
double *score,double *hess,double *est,double *var,int *sim,int *antsim,int *rani,double *test,double *testOBS,double *Ut,double *simUt,int *weighted,
double *gamma,double *vargamma,int *semi,double *zsem,int *pg,int *trans,double *gamma2,int *CA,int *line,int *detail,double *biid,double *gamiid,int *resample,
double *timepow,int *clusters,int *antclust,double *timepowtest,int *silent,double *convc,double *tau,int *estimator,int *causeS,double *weights,
double *KMtimes,int *ordertime,int *conservative,int *censcode)
//double *times,*betaS,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
//       *Ut,*simUt,*gamma,*zsem,*gamma2,*biid,*gamiid,*vargamma,*timepow,
//       *weights,*KMtimes,*timepowtest,*convc,*tau;
//int *n,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*rani,*weighted,
//    *semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,*silent,
//    *estimator,*causeS,*ordertime,*conservative,*censcode;
{ // {{{
  // {{{ allocation and reading of data from R
  matrix *X,*cX,*A,*AI,*cumAt[*antclust],*VAR,*Z,*censX;
  vector *VdB,*risk,*SCORE,*W,*Y,*Gc,*DELTA,*CAUSE,*bhat,*pbhat,*beta,*xi,*censXv,
    *rr,*rowX,*difbeta,*qs,*bhatub,*betaub,*dcovs,*pcovs,*zi,*rowZ,*zgam; 
  vector *cumhatA[*antclust],*cumA[*antclust],*bet1,*gam,*dp,*dp1,*dp2; 
  int osilent,convt,ps,sing,c,i,j,k,l,s,it,convproblems=0,clusterj,nrisk; 
  double skm,rit,time,sumscore,totrisk,*vcudif=calloc((*Ntimes)*(*px+1),sizeof(double));
//  float gasdev(),expdev(),ran1();
//  void resmeansemi();
  ps=(*px); 

//  printf(" %d %d %d %d %d %d \n",*px,*semi,*Ntimes,*trans,*antclust,*n); 
//  printf(" %d \n",ps); 

  if (*semi==0) { 
    osilent=silent[0]; silent[0]=0; 
    malloc_mat(*n,*px,X); malloc_mat(*n,*px,cX); 
    if (*trans==2) {malloc_mat(*n,*pg,Z);malloc_vecs(*pg,&zgam,&gam,&zi,&rowZ,NULL);}
    malloc_mats(ps,ps,&A,&AI,&VAR,NULL); 

    malloc_vecs(*n,&rr,&bhatub,&risk,&W,&Y,&Gc,&DELTA,&CAUSE,&bhat,&pbhat,NULL); 
    malloc_vecs(*px,&bet1,&xi,&rowX,&censXv,NULL); 
    malloc_vecs(ps,&dp,&dp1,&dp2,&dcovs,&pcovs,&betaub,&VdB,&qs,&SCORE,&beta,&difbeta,NULL); 

    malloc_mats(*n,*px,&censX,NULL);

    for (i=0;i<*antclust;i++) {
      malloc_vec(ps,cumhatA[i]); malloc_vec(ps,cumA[i]); 
      malloc_mat(*Ntimes,ps,cumAt[i]);
    }

    for (c=0;c<ps;c++) VE(beta,c)=betaS[c]; 
    for (c=0;c<ps;c++) VE(bet1,c)=betaS[c]; 
    for (c=0;c<*n;c++) {VE(Gc,c)=KMc[c]; VE(DELTA,c)=delta[c]; 
      VE(CAUSE,c)=cause[c]; 
      for(j=0;j<*px;j++)  ME(X,c,j)=z[j*(*n)+c]; 
//    printf(" %d %d %d  \n",*causeS,cause[c],cause[c]==*causeS); 
    }

    // }}}
    
//head_matrix(X); 
	 
for (s=0;s<*Ntimes;s++)
{
   time=times[s]; est[s]=time; score[s]=time; var[s]=time;
   convt=1;  

  for (it=0;it<*Nit;it++)
  {
    totrisk=0; 

    // could start later than j=0  if data is sorted after timing
    // and pass this to c as j=startn[s]
    for (j=0;j<*n;j++) { // {{{ computation of P1 and DP1  and observed response 

// {{{ setting up at risk in different situations

//    skm=sqrt(weights[j]*KMtimes[s]/KMc[j]); 
    skm=sqrt(weights[j]*KMtimes[s]/KMc[j]); 

    if (*estimator==1) { // standard conditional residual 
       VE(risk,j)=(x[j]>=time); 
       rit= (x[j]>=time); 
       rit=rit*delta[j]*skm; 
    } else if (*estimator==2) { // cause specific YL to cause  given event
       VE(risk,j)=(x[j]<= (*tau))*(cause[j]==*causeS); 
       rit=(x[j]<= (*tau))*(cause[j]==*causeS);
       rit=rit*delta[j]*skm; 
    } else  if (*estimator==3) { // PKA years lost 
       VE(risk,j)=(x[j]>= time); // *(cause[j]==*causeS); 
       rit=(x[j]>= time); // *(cause[j]==*causeS);
       rit=rit*delta[j]*skm; 
    } else  if (*estimator==4) { // inside weight
       rit= 1; 
       VE(risk,j)=rit; 
       skm=sqrt(weights[j]*KMtimes[s]); 
       rit=weights[j]; 
    }
       // }}}

    totrisk=totrisk+VE(risk,j); 
    extract_row(X,j,xi); 
    VE(bhat,j)=vec_prod(xi,bet1); 

    if (*trans==1) { VE(pbhat,j)=VE(bhat,j); scl_vec_mult(rit,xi,dp); }
    if (*trans==2) { VE(pbhat,j)=exp(VE(bhat,j)); scl_vec_mult(rit*VE(pbhat,j),xi,dp); }
    if ((*trans==1 ) || (*trans==2)) { replace_row(cX,j,dp); }


    if (*estimator==1) VE(Y,j)=(((x[j]-time))-VE(pbhat,j))*rit; 
    else  if (*estimator==2) VE(Y,j)=((*tau-x[j])-VE(pbhat,j))*rit; 
    else  if (*estimator==3) VE(Y,j)=((*tau-x[j])*(cause[j]==*causeS)-VE(pbhat,j))*rit; 
    else  if (*estimator==4) VE(Y,j)=((x[j]-time)*(cause[j]==*causeS)*delta[j]/KMc[j]-VE(pbhat,j))*rit; 

   if (it==(*Nit-1) && (*conservative==0)) { // {{{ for censoring distrubution
	   scl_vec_mult(VE(Y,j),dp,dp1); 
           vec_add(censXv,dp1,censXv); 
           replace_row(censX,j,dp1);
      } // }}}


//    if (it==*Nit-1) {
////	if (KMc[j]<0.00001) vec_zeros(dp); else scl_vec_mult(1/KMc[j],dp,dp); 
//	scl_vec_mult(VE(Y,j),dp,dp); vec_add(dp,qs,qs); 
//    }

    } // }}}

//    head_matrix(cX); 
//    printf("================ %lf \n",vec_sum(Y)); 
//    printf("================ %lf \n",vec_sum(risk)); 
//    printf("================ %lf \n",vec_sum(pbhat)); 
    
    totrisk=vec_sum(risk); 
    MtA(cX,cX,A); 
    invertS(A,AI,osilent); 
//    print_mat(A); 
    sing=0; 

    if (fabs(ME(AI,0,0))<.0000001) {
      convproblems=1; convt=0; silent[s]=1;
      for (c=0;c<ps;c++) VE(beta,c)=0; 
      for (c=0;c<*px;c++) VE(bet1,c)=betaS[c]; 
      sing=1;
      if (osilent==0) Rprintf("Non-invertible design at time %lf \n",time); 
      it=*Nit-1;  
    }
    if (sing==0) {
      /* print_vec(Y); print_vec(SCORE); print_vec(difbeta); */ 
      vM(cX,Y,SCORE); 
//      print_vec(SCORE); 
      Mv(AI,SCORE,difbeta); 
      vec_add(beta,difbeta,beta); 
      for (i=0;i<*px;i++) VE(bet1,i)=VE(beta,i); 

      sumscore=0; 
      for (k=0;k<*px;k++) sumscore=sumscore+fabs(VE(difbeta,k)); 
      if ((sumscore<*convc) & (it<*Nit-2)) it=*Nit-2;

      if (isnan(vec_sum(SCORE))) {
	Rprintf("missing values in SCORE %ld \n",(long int) s); 
	convproblems=1; convt=0; silent[s]=2;
	it=*Nit-1; 
	for (c=0;c<ps;c++) { VE(beta,c)=0; VE(SCORE,c)=99; }
	for (c=0;c<*px;c++) VE(bet1,c)=betaS[c]; 
      }
    }

    if (*detail==1) { 
      Rprintf(" s er %ld, Estimate beta \n",(long int) s); print_vec(beta); 
      Rprintf("Score D l\n"); print_vec(difbeta); 
      Rprintf("Information -D^2 l\n"); print_mat(AI); };

    if (it==*Nit-1) scl_vec_mult(1/totrisk,qs,qs); 
  } /* it */

vec_zeros(VdB); mat_zeros(VAR); 

if (convt==1 ) {
//   for (j=0;j<*antclust;j++) {vec_zeros(cumA[j]);vec_zeros(cumhatA[j]);}
   for (i=0;i<*n;i++) { 
      j=clusters[i]; 
      if (s<-1) Rprintf("%d  %d %d \n",s,i,j);
      extract_row(cX,i,dp); 
      scl_vec_mult(VE(Y,i),dp,dp); 
      vec_add(dp,cumA[j],cumA[j]); 
      if ((*conservative==0)) { // {{{ censoring terms for variance 
 	k=ordertime[i]; nrisk=(*n)-i; clusterj=clusters[k]; 
//	printf(" %d %d %lf %lf %lf %d \n",i,k,nrisk,time,x[k],cause[k]); 
	if (cause[k]==(*censcode)) { 
           for(l=0;l<ps;l++) VE(cumA[clusterj],l)+=VE(censXv,l)/nrisk; 
           for (j=i;j<*n;j++) {
             clusterj=clusters[ordertime[j]]; 	
             for(l=0;l<ps;l++) VE(cumA[clusterj],l)-=VE(censXv,l)/pow(nrisk,2); 
	   }
	}
        // fewer where I(s <= T_i) , because s is increasing
        extract_row(censX,k,xi); vec_subtr(censXv,xi,censXv);  
    } // }}}
//      if ((time==x[i])&(delta[i]==0))vec_add(qs,cumhatA[j],cumhatA[j]);  
      if (s<-1) print_vec(dp2); 
   }

   for (j=0;j<*antclust;j++) { 
//      vec_add(cumhatA[j],cumA[j],dp1); 
      Mv(AI,cumA[j],dp2); 
      replace_row(cumAt[j],s,dp2);  

      for(k=0;k<ps;k++) 
      for(c=0;c<ps;c++) ME(VAR,k,c)=ME(VAR,k,c)+VE(dp2,k)*VE(dp2,c); 

      if (*resample==1) {
      for (c=0;c<*px;c++) {l=j*(*px)+c; biid[l*(*Ntimes)+s]=VE(dp2,c);}
      }
   }
}

   for (i=1;i<ps+1;i++) {
      var[i*(*Ntimes)+s]=ME(VAR,i-1,i-1); 
      est[i*(*Ntimes)+s]=VE(beta,i-1); 
      score[i*(*Ntimes)+s]=VE(SCORE,i-1); 
   }

} /* s=1 ... *Ntimes */ 


    if (*sim==1)
      comptestfunc(times,Ntimes,px,est,var,vcudif,antsim,test,testOBS,Ut,
		   simUt,cumAt,weighted,antclust,gamma2,line,timepowtest); 
  } else {
    resmeansemi(times,Ntimes,x,delta,cause,KMc,z,n,px,Nit,
	      score,hess,est,var,sim,antsim,rani,test,testOBS,Ut,simUt,weighted,
	      gamma,vargamma,semi,zsem,pg,trans,gamma2,CA,line,detail,biid,
	      gamiid,resample,timepow,clusters,antclust,timepowtest,silent,convc,tau,estimator,causeS,
              weights,KMtimes); 
  }
 
  if (convproblems>0) convc[0]=1; 
  if (*semi==0) { 
    free_mats(&censX,&VAR,&X,&cX,&A,&AI,NULL); 
    if (*trans==2) {free_mats(&Z,NULL); free_vecs(&zgam,&gam,&zi,&rowZ,NULL);}

    free_vecs(&censXv,&rr,&bhatub,&risk,&W,&Y,&Gc,&DELTA,&CAUSE,&bhat,&pbhat,NULL); 
    free_vecs(&bet1,&xi,&rowX,NULL); 
    free_vecs(&dp,&dp1,&dp2,&dcovs,&pcovs,&betaub,&VdB,&qs,&SCORE,&beta,&difbeta,NULL); 

    for (i=0;i<*antclust;i++) {free_vec(cumhatA[i]); free_vec(cumA[i]); 
    free_mat(cumAt[i]);}

  }
free(vcudif); 
} // }}}

//double *times,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
//*Ut,*simUt,*gamma,*zsem,*vargamma,*gamma2,*biid,*gamiid,*timepow,*timepowtest,
//               *weights,*KMtimes,
//	*convc,*tau;
//int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*rani,*weighted,
//*semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,*silent,*estimator,*causeS;

void resmeansemi(double *times,int  *Ntimes,double *x,int *delta,int *cause,
	       double *KMc,double *z,int *antpers,int *px,int *Nit,
	       double *score,double *hess,double *est,double *var,int *sim,
	       int *antsim,int *rani,double *test,double *testOBS,double *Ut,
	       double *simUt,int *weighted,double *gamma,double *vargamma,int *semi,
	       double *zsem,int *pg,int *trans,double *gamma2,int *CA,
	       int *line,int *detail,double *biid,double *gamiid,int *resample,
	       double *timepow,int *clusters,int *antclust,double *timepowtest,int *silent,double *convc,double *tau,int *estimator,int *causeS, double *weights,double *KMtimes) { 
  // {{{ allocation and reading of data from R
  matrix *ldesignX,*A,*AI,*cdesignX,*ldesignG,*cdesignG;
  matrix *S,*dCGam,*CGam,*ICGam,*VarKorG,*dC,*XZ,*ZZ,*ZZI,*XZAI; 
  matrix *Ct,*C[*Ntimes],*Acorb[*Ntimes],*tmpM1,*tmpM2,*tmpM3,*tmpM4; 
  matrix *Vargam,*dVargam,*M1M2[*Ntimes],*Delta,*dM1M2,*M1M2t,*RobVargam;
  matrix *W3t[*antclust],*W4t[*antclust];
  vector *W2[*antclust],*W3[*antclust];
  vector *diag,*dB,*dN,*VdB,*AIXdN,*AIXlamt,*bhatt,*pbhat,*plamt;
  vector *korG,*pghat,*rowG,*gam,*dgam,*ZGdN,*IZGdN,*ZGlamt,*IZGlamt;
  vector *covsx,*covsz,*qs,*Y,*rr,*bhatub,*xi,*rowX,*rowZ,*difX,*zi,*z1,*tmpv1,*tmpv2,*lrisk;
  int sing,itt,i,j,k,l,s,c,pmax,totrisk,convproblems=0,fixedcov,osilent, 
      *n= calloc(1,sizeof(int)), *nx= calloc(1,sizeof(int)),
      *robust= calloc(1,sizeof(int));
  double lrr,skm,rit,dtau,time,dummy,dtime;
  double *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double)),
	 *inc=calloc((*Ntimes)*(*px+1),sizeof(double));
  osilent=silent[0]; silent[0]=0; 
  robust[0]=1; fixedcov=1; 
  n[0]=antpers[0]; nx[0]=antpers[0];

//if (*trans==1) for (j=0;j<*pg;j++) if (fabs(timepow[j]-1)>0.0001) {timem=1;break;}
//if (*trans==2) for (j=0;j<*pg;j++) if (fabs(timepow[j])>0.0001) {timem=1;break;}

 for (j=0;j<*antclust;j++) { 
    malloc_mat(*Ntimes,*px,W3t[j]);
    malloc_mat(*Ntimes,*px,W4t[j]); malloc_vec(*pg,W2[j]); malloc_vec(*px,W3[j]);
 }

  malloc_mats(*antpers,*px,&ldesignX,&cdesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,&cdesignG,NULL); 
  malloc_mats(*px,*px,&tmpM1,&A,&AI,NULL);
  malloc_mats(*pg,*pg,&dVargam,&Vargam,&RobVargam,&tmpM2,&ZZ,&VarKorG,&ICGam,&CGam,&dCGam,&S,&ZZI,NULL); 
  malloc_mats(*px,*pg,&XZAI,&tmpM3,&Ct,&dC,&XZ,&dM1M2,&M1M2t,NULL);
  malloc_mat(*px,*pg,tmpM4); 
  for (j=0;j<*Ntimes;j++) { malloc_mat(*pg,*px,Acorb[j]); 
    malloc_mat(*px,*pg,C[j]); malloc_mat(*px,*pg,M1M2[j]);
  }
  malloc_mat(*Ntimes,*px,Delta); malloc_mat(*Ntimes,*px,tmpM1);

  malloc_vecs(*px,&covsx,&xi,&rowX,&difX,&tmpv1,&korG,&diag,&dB,&VdB,&AIXdN,&AIXlamt,&bhatt,NULL);
  malloc_vecs(*pg,&covsz,&zi,&rowZ,&tmpv2,&zi,&z1,&rowG,&gam,&dgam,&ZGdN,&IZGdN,&ZGlamt,&IZGlamt,NULL);
  malloc_vecs(*antpers,&Y,&bhatub,&rr,&lrisk,&dN,&pbhat,&pghat,&plamt,NULL);
  malloc_vec((*px)+(*pg),qs); 

  if (*px>=*pg) pmax=*px; else pmax=*pg; 
  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 
  // }}}

  if (fixedcov==1) {
    for (c=0;c<*antpers;c++) {
      for(j=0;j<pmax;j++)  {
	if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c];
	if (j<*pg) ME(ldesignG,c,j)=zsem[j*(*antpers)+c]; } } 
  }

  for (itt=0;itt<*Nit;itt++)
    {
      mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZGdN); vec_zeros(IZGlamt); 

      Mv(ldesignG,gam,pghat);
      for (s=0;s<*Ntimes;s++)
      {
	  time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 
//	  dtime=1; 
		    if (*tau>0) dtau=*tau-time; else dtau=time; 
//		    printf(" %lf %lf %lf \n",time,dtau,*tau); 

	  for(j=1;j<=*px;j++) VE(bhatt,j-1)=est[j*(*Ntimes)+s];
	  Mv(ldesignX,bhatt,pbhat); 

	  totrisk=0; 
	  for (j=0;j<*antpers;j++) {  // {{{ 
	  skm=sqrt(weights[j]*KMtimes[s]/KMc[j]); 

	    if (*estimator==2) { // cause specific YL to cause  given event
	       VE(lrisk,j)=(x[j]<= (*tau))*(cause[j]==*causeS); 
	       rit=(x[j]<= (*tau))*(cause[j]==*causeS);
	       rit=rit*delta[j]*skm; 
	    } else  if (*estimator==3) { // PKA years lost 
	       VE(lrisk,j)=(x[j]>= time); // *(cause[j]==*causeS); 
	       rit=(x[j]>= time); // *(cause[j]==*causeS);
	       rit=rit*delta[j]*skm; 
	    } else  if (*estimator==4) { // inside weight conditional residual
	       rit= 1; 
	       VE(lrisk,j)=rit; 
	       skm=sqrt(weights[j]*KMtimes[s]); 
	       rit=1; 
	    } else  if (*estimator==1) { // standard conditional residual 
	       VE(lrisk,j)=(x[j]>=time); 
	       rit= (x[j]>=time); 
	       rit=rit*delta[j]*skm; 
	    }


	    totrisk=totrisk+VE(lrisk,j);
	    extract_row(ldesignX,j,xi); 
	    extract_row(ldesignG,j,zi); 

	    lrr=0; 
	    // {{{ compute P_1 and DP_1 
	    if (*trans==1 ) {
		for (l=0;l<*pg;l++) lrr=lrr+VE(gam,l)*VE(zi,l)*pow(dtau,timepow[l]); 
	      VE(plamt,j)=VE(pbhat,j)+lrr; 
	      scl_vec_mult(rit,xi,xi); scl_vec_mult(rit,zi,zi);  
              for (l=0;l<*pg;l++) VE(zi,l)=pow(dtau,timepow[l])*VE(zi,l); 
	    }
	   if (*trans==2) {
	        for (l=0;l<*pg;l++) lrr=lrr+VE(gam,l)*VE(zi,l)*pow(dtau,timepow[l]); 
	      VE(rr,j)=lrr;  
	      VE(plamt,j)=exp(VE(pbhat,j)+lrr); 
	      scl_vec_mult(rit*VE(plamt,j),xi,xi); 
	      scl_vec_mult(rit*VE(plamt,j),zi,zi); 
	      for (l=0;l<*pg;l++) VE(zi,l)= pow(dtau,timepow[l])*VE(zi,l); 
           }
	   // }}}
          
	    if ((*trans==1 ) || (*trans==2)) {
	       replace_row(cdesignX,j,xi); replace_row(cdesignG,j,zi); 
	    }

	  /*
	   if (itt==*Nit-1) {
	   if (KMc[j]<0.00001) vec_zeros(xi); else scl_vec_mult(1/KMc[j],xi,xi); 
	       scl_vec_mult(VE(lrisk,j),xi,xi); vec_add(xi,qs,qs); 
	   }
	  */
	  if (*estimator==1) VE(Y,j)=(((x[j]-time))-VE(pbhat,j))*rit; 
	  else  if (*estimator==3) VE(Y,j)=((*tau-x[j])*(cause[j]==*causeS)-VE(pbhat,j))*rit; 
	  else  if (*estimator==4) VE(Y,j)=((x[j]-time)*(cause[j]==*causeS)*delta[j]/KMc[j]-VE(pbhat,j))*rit; 
	  else  if (*estimator==2) VE(Y,j)=((*tau-x[j])-VE(pbhat,j))*rit; 

	} // j=1..antpers ## }}} 

	  MtA(cdesignX,cdesignX,A); 
	  invertS(A,AI,osilent); sing=0; 

          if (fabs(ME(AI,0,0))<.0000001) {
             convproblems=1;  silent[s]=1; 
             if (osilent==0) Rprintf("Iteration %d: non-invertible design at time %lf\n",itt,time); 
	     for (k=1;k<=*px;k++) inc[k*(*Ntimes)+s]=0; 
	     sing=1;
          }

	  if (sing==0) { 
	  vM(cdesignX,Y,xi); Mv(AI,xi,AIXdN); 
	  MtA(cdesignG,cdesignG,ZZ); MtA(cdesignX,cdesignG,XZ);
	  MxA(AI,XZ,XZAI); MtA(XZAI,XZ,tmpM2); 
	  mat_subtr(ZZ,tmpM2,dCGam); 
	  scl_mat_mult(dtime,dCGam,dCGam); mat_add(CGam,dCGam,CGam); 

	  vM(cdesignG,Y,zi); vM(XZ,AIXdN,tmpv2); 
	  vec_subtr(zi,tmpv2,ZGdN); scl_vec_mult(dtime,ZGdN,ZGdN); 
	  vec_add(ZGdN,IZGdN,IZGdN); 
	  Acorb[s]=mat_transp(XZAI,Acorb[s]); 
	  C[s]=mat_copy(XZ,C[s]); 

	  /* scl_mat_mult(dtime,XZAI,tmpM4);mat_add(tmpM4,Ct,Ct); */
	  for (k=1;k<=*px;k++) inc[k*(*Ntimes)+s]=VE(AIXdN,k-1); 
	  }

	  if (itt==*Nit-1) {
	    for (i=0;i<*antpers;i++) 
            { // vec_zeros(tmpv1); vec_zeros(z1); 
              j=clusters[i]; 	
	      extract_row(cdesignX,i,xi); scl_vec_mult(VE(Y,i),xi,xi); 
	      Mv(AI,xi,rowX);
	      extract_row(cdesignG,i,zi); scl_vec_mult(VE(Y,i),zi,zi); 
	      vM(C[s],rowX,tmpv2); vec_subtr(zi,tmpv2,rowZ); 
	      scl_vec_mult(dtime,rowZ,rowZ); 
	     // vec_add(rowZ,z1,z1); 
	     // vec_add(rowX,tmpv1,tmpv1); 
	      vec_add(rowZ,W2[j],W2[j]); 
	      for (k=0;k<*px;k++) ME(W3t[j],s,k)= ME(W3t[j],s,k)+VE(rowX,k); 
	    }  
	 }
	} /* s=1,...Ntimes */

      invertS(CGam,ICGam,osilent); Mv(ICGam,IZGdN,dgam); vec_add(gam,dgam,gam); 

      if (isnan(vec_sum(dgam))) {
         if (convproblems==1) convproblems=3;  else convproblems=2; 
         if (osilent==1) Rprintf("missing values in dgam %ld \n",(long int) s);
	 vec_zeros(gam); 
      }

      dummy=0; for (k=0;k<*pg;k++)  dummy=dummy+fabs(VE(dgam,k)); 

      for (s=0;s<*Ntimes;s++) {
	vM(Acorb[s],dgam,korG); 
	est[s]=times[s]; var[s]=times[s]; 
	for (k=1;k<=*px;k++)  { 
            est[k*(*Ntimes)+s]=
            est[k*(*Ntimes)+s]+inc[k*(*Ntimes)+s]-VE(korG,k-1); 
	  dummy=dummy+fabs(inc[k*(*Ntimes)+s]-VE(korG,k-1)); 
	  /* printf(" %lf ",est[k*(*Ntimes)+s]); printf(" \n");*/ }
      } /* s=1,...Ntimes */
      if (dummy<*convc && itt<*Nit-2) itt=*Nit-2; 

      if (*detail==1) { 
	Rprintf(" iteration %d %d \n",itt,*Nit); 
	Rprintf("Total sum of changes %lf \n",dummy); 
	Rprintf("Gamma parameters \n"); print_vec(gam); 
	Rprintf("Change in Gamma \n"); print_vec(dgam); }

    } /*itt lokke */ 

  /* ROBUST VARIANCES   */ 
  if (*robust==1) 
    {
      for (s=0;s<*Ntimes;s++) {
	vec_zeros(VdB); 
	 for (i=0;i<*antclust;i++) {

	  Mv(ICGam,W2[i],tmpv2); vM(Acorb[s],tmpv2,rowX); 
	  extract_row(W3t[i],s,tmpv1); vec_subtr(tmpv1,rowX,difX); 
	  replace_row(W4t[i],s,difX); vec_star(difX,difX,tmpv1); 
	  vec_add(tmpv1,VdB,VdB);

	  if (*resample==1) {
	    if (s==1)
	      for (c=0;c<*pg;c++)
	      gamiid[c*(*antclust)+i]=gamiid[c*(*antclust)+i]+VE(tmpv2,c);
	    for (c=0;c<*px;c++) {l=i*(*px)+c; 
	      biid[l*(*Ntimes)+s]=biid[l*(*Ntimes)+s]+VE(difX,c);} }


	  if (s==0) { for (j=0;j<*pg;j++) for (k=0;k<*pg;k++) 
			ME(RobVargam,j,k)=ME(RobVargam,j,k)+VE(tmpv2,j)*VE(tmpv2,k);} 
	}  /* for (i=0;i<*antclust;i++) */ 
	for (k=1;k<*px+1;k++) var[k*(*Ntimes)+s]=VE(VdB,k-1); 

      } /* s=0..Ntimes*/
    }

  /* MxA(RobVargam,ICGam,tmpM2); MxA(ICGam,tmpM2,RobVargam);*/
  /* print_mat(RobVargam);  */ 

  for (j=0;j<*pg;j++) {gamma[j]=VE(gam,j);
    for (k=0;k<*pg;k++) {vargamma[k*(*pg)+j]=ME(RobVargam,j,k);}}

  if (convproblems>=1) convc[0]=convproblems; 
  if (*sim==1) {
    comptestfunc(times,Ntimes,px,est,var,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antclust,gamma2,line,timepowtest);
  }

  free_mats(&ldesignX,&A,&AI,&cdesignX,&ldesignG,&cdesignG,
	      &S,&dCGam,&CGam,&ICGam,&VarKorG,&dC,&XZ,&ZZ,&ZZI,&XZAI, 
	      &Ct,&tmpM1,&tmpM2,&tmpM3,&tmpM4,&Vargam,&dVargam,
	      &Delta,&dM1M2,&M1M2t,&RobVargam,NULL); 

  free_vecs(&qs,&Y,&rr,&bhatub,&diag,&dB,&dN,&VdB,&AIXdN,&AIXlamt,
	      &bhatt,&pbhat,&plamt,&korG,&pghat,&rowG,&gam,&dgam,&ZGdN,&IZGdN,
	      &ZGlamt,&IZGlamt,&xi,&rowX,&rowZ,&difX,&zi,&z1,&tmpv1,&tmpv2,&lrisk,
	      NULL); 

  for (j=0;j<*Ntimes;j++) {free_mat(Acorb[j]);free_mat(C[j]);free_mat(M1M2[j]);}
  for (j=0;j<*antclust;j++) {free_mat(W3t[j]); free_mat(W4t[j]);
    free_vec(W2[j]); free_vec(W3[j]); }
  free(vcudif); free(inc); 
  free(n); free(nx);  free(robust); 
}

