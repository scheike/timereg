#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include <Rdefines.h>
#include <R.h>
#include "haplosurv.h" 
                 

void haplocif(
// {{{
times,Ntimes,x,delta,cause,KMc,z,n,px,Nit,betaS, 
score,hess,est,var,sim,antsim,rani,test,testOBS,Ut,simUt,weighted,
gamma,vargamma,semi,zsem,pg,trans,gamma2,CA,line,detail,biid,gamiid,resample,
timepow,clusters,antclust,fixhaplo,
haplofreq,alphaiid,hapdim,
nph,oh,nphpp,designfuncX,designfuncZ,rhoR,dimxih,dimzih,haplodes,designtest)
double *times,*betaS,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
*Ut,*simUt,*gamma,*zsem,*gamma2,*biid,*gamiid,*vargamma,*haplofreq,*alphaiid,*haplodes,*timepow;
int *n,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*rani,*weighted,
*semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,
*fixhaplo,*nph,*oh,*nphpp,*hapdim,*dimxih,*dimzih,*designtest;
SEXP designfuncX,designfuncZ,rhoR; 
// }}}
{
 // {{{ variable defs 
  matrix *X,*cX,*A,*AI,*cumAt[*antclust],*VAR,*DnuP1;
  matrix *DUeta[*Ntimes],*DUetaD,*DUetaP,*hapDes;
  vector *VdB,*riskV,*SCORE,*W,*Y,*Gc,*DELTA,*CAUSE,*bhat,*pbhat,*beta,*xi,
    *rr,*rowX,*difbeta,*qs,*bhatub,*betaub;
  vector *cumhatA[*antclust],*cumA[*antclust],*bet1,*dp,*dp1,*dp2; 
  vector *hapiid[*antclust]; 
  vector *dph,*xih1,*xih,*rowXh; 
  vector *DPGnuD,*DPHnuD,*DP1nuD; 
  int c1D,ps,sing,sc,c,i,j,k,l,s,it;
  int  iD,pps[1],haplotypeD[2],risk[*antclust]; 
  double phD,time,sumscore,totrisk,
	 *vcudif=calloc((*Ntimes)*(*dimxih+1), sizeof(double)),ph,pgeno; 
  long idum;
  // void Fhaplodes(),comptestfunc(),semihaplo(),fDUeta(); 
  idum=*rani;
  pps[0]=ps=(*dimxih); 
// }}}

if (*semi==0) { 
// {{{ matrix allocation
    malloc_mat(*n,*px,X); 
    malloc_mat(*n,*dimxih,cX); 
    malloc_mat(*n,*dimxih,DnuP1); 
    malloc_mats(ps,ps,&A,&AI,&VAR,NULL); 

    malloc_vecs(*n,&rr,&bhatub,&riskV,&W,&Y,&Gc,&DELTA,&CAUSE,&bhat,&pbhat,NULL); 
    malloc_vecs(*px,&xi,&rowX,NULL); 
    malloc_vecs(ps,&bet1,&dph,&xih1,&rowXh,&xih,NULL); 
    malloc_vecs(ps,&dp,&dp1,&dp2,&betaub,&VdB,&qs,&SCORE,&beta,&difbeta,NULL); 

    for (i=0;i<*antclust;i++) {
      malloc_vec(*hapdim,hapiid[i]); 
      malloc_vec(ps,cumhatA[i]); malloc_vec(ps,cumA[i]); 
      malloc_mat(*Ntimes,ps,cumAt[i]);}

  for (j=0;j<*Ntimes;j++) { malloc_mat(*dimxih,*hapdim,DUeta[j]); }
  malloc_mats(*dimxih,*hapdim,&DUetaD,&DUetaP,NULL); 
  malloc_mat(*nph-1,*hapdim,hapDes); 
  malloc_vecs(*hapdim,&DPGnuD,&DPHnuD,&DP1nuD,NULL);
 // }}}
   
// {{{ initialisering and reading of starting values 

    if (*fixhaplo==0) {
    for (i=0;i<*antclust;i++) 
    for (j=0;j<*hapdim;j++) VE(hapiid[i],j)=alphaiid[j*(*n)+i];
    // print_vec(hapiid[i]); 
    //
    for (i=0;i<*nph-1;i++)
    for(j=0;j<*hapdim;j++)  ME(hapDes,i,j)=haplodes[j*(*nph-1)+i];
    }

    for (c=0;c<ps;c++) VE(beta,c)=betaS[c]; 
    for (c=0;c<ps;c++) VE(bet1,c)=betaS[c]; //  est[c*(*Ntimes)+0];

    for (c=0;c<*n;c++) {VE(Gc,c)=KMc[c]; VE(DELTA,c)=delta[c]; 
      VE(CAUSE,c)=cause[c]; 
      for(j=0;j<*px;j++)  ME(X,c,j)=z[j*(*n)+c]; 
    }

// }}}

  // reading in basic matrices X and Z with observed covariates  {{{
  // making haplobased design for all subject and the relevant haplotypes
  
  int nallH=0; 
  for (i=0;i<*n;i++) nallH+=nphpp[i]; 
  matrix *XallH;
  double phallH[nallH]; 
  int ih,first,indexpersallH[*n],indexpersnph[*n]; 
  //vector *ZallHbeta; 

  malloc_mat(nallH,*dimxih,XallH); 

  // head_matrix(XallH); head_matrix(ZallH); 
  // print_matrix(XallH); print_matrix(ZallH); 
  
  k=0;c1D=0; 
  for (i=0;i<*n;i++) { 
   indexpersnph[i]=c1D; first=0; 
   extract_row(X,i,xi); 

    for (j=0;j<nphpp[i];j++) {
       haplotypeD[0]=oh[c1D]; haplotypeD[1]=oh[c1D+1]; 
       ph=haplofreq[oh[c1D]]*haplofreq[oh[c1D+1]]; 
       phallH[k]=ph; 

       Fhaplodes(xi,haplotypeD,xih,designfuncX,rhoR);; 

       if (*detail>=3) { // print test {{{
	  Rprintf("======Design =================== \n"); 
          Rprintf(" person %d  \n",i); 
          Rprintf(" haplotypes hap1 hap2 %d %d  \n",oh[c1D],oh[c1D+1]);  
          Rprintf(" design, x \n");  
          print_vec(xi); // print_vec(ziI); 
          Rprintf(" design, x(h) \n");  
	  print_vec(xih); // print_vec(zih);
          ph=haplofreq[oh[c1D]]*haplofreq[oh[c1D+1]];
          Rprintf("probability for haplotype pair Donor %lf %d %d \n",ph,oh[c1D],oh[c1D+1]); 
       } // }}}
     
       
       replace_row(XallH,k,xih); 
       if (first==0) {indexpersallH[i]=k; first=1;}
       c1D+=2; k+=1; 
    }
   // Rprintf(" %d %d %d \n",i,indexpersallH[i],indexpersnph[i]); 
  }
  // }}}

 //  head_matrix(XallH); 

sc=0;
for (s=0;s<*Ntimes;s++)
{
	time=times[s]; est[s]=time; score[s]=time; var[s]=time;
   R_CheckUserInterrupt();

       for (it=0;it<*Nit;it++)
       {
        totrisk=0;  
        c1D=0; 
	
       for (j=0;j<*n;j++) { 
            ih=indexpersallH[j]; 
// Rprintf("================================= \n"); 
// Rprintf(" %lf %lf %d %d %d %d %d \n",x[j],time,cause[j],*CA,*n,j,*antclust); 
	    risk[j]=(x[j]>=time); totrisk=totrisk+risk[j];
	    VE(riskV,j)=risk[j]; 
	    // extract_row(X,j,xi); 

	  // initializing vectors for subject sums // {{{
	    pgeno=0; vec_zeros(dph); VE(pbhat,j)=0; 
	    vec_zeros(DPGnuD); vec_zeros(DP1nuD); 
	   // }}}

	  VE(Y,j)=((x[j]<=time) & (cause[j]==*CA))*1;

    // Rprintf("it er %d %d %d  %d \n",it,j,nphpp[j],c1); 
	  
      for (iD=0;iD<nphpp[j];iD++) { 

        haplotypeD[0]=oh[c1D]; haplotypeD[1]=oh[c1D+1]; 
        phD=haplofreq[oh[c1D]]*haplofreq[oh[c1D+1]];
        pgeno=pgeno+phD; 

        if (it==*Nit-1 && *fixhaplo==0) {// computes derivative of P(H_D=h_D) 
           fDUeta(c1D,oh,nph,hapdim,haplofreq,hapDes,DPHnuD); 
           vec_add(DPHnuD,DPGnuD,DPGnuD);
        }

         ph=phD; 
  //       Fhaplodes(xi,haplotypeD,xih,designfuncX,rhoR);; 
         
         extract_row(XallH,ih,xih);  
         vec_star(xih,bet1,rowXh); VE(bhat,j)=vec_sum(rowXh); 

	 if (j<*designtest && it==0) { // print test {{{
	    Rprintf("timepoint=%d person=%d #haplotype-pairs=%d \n",s,j,nphpp[j]); 
	    Rprintf("Donor haplotype1=%d haplotype2=%d  \n",oh[c1D],oh[c1D+1]);  
	    Rprintf("input x and f(x,h) design depending on haplotype \n"); 
	    print_vec(xi); print_vec(xih);
	 }
	   // print test }}}

         // {{{ computes P1 and D_P1  
         if (*trans==1) {VE(pbhat,j)+=ph*(1-exp(-VE(bhat,j)));
                         vec_add_mult(dph,xih,exp(-VE(bhat,j))*ph,dph);
	 }
	    if (*trans==2) {
	       VE(pbhat,j)+=ph*(1-exp(-exp(VE(bhat,j)))); 
	       vec_add_mult(dph,xih,ph*exp(-exp(VE(bhat,j)))*exp(VE(bhat,j)),dph); }
         // }}}
	
       if (it==*Nit-1 && *fixhaplo==0){// {{{ DP1nu 
          if (*trans==1) 
	    vec_add_mult(DP1nuD,DPHnuD,(1-exp(-VE(bhat,j))),DP1nuD); 
          if (*trans==2) 
	    vec_add_mult(DP1nuD,DPHnuD,(1-exp(-exp(VE(bhat,j)))),DP1nuD); 
       } // }}}
 
	   c1D+=2;  
	   ih+=1; 
       } // iD=0 ... nphpp[j] 

       if (it==*Nit-1 && *fixhaplo==0) {// computes derivative of DUnu 
          scl_vec_mult(1/(pgeno),DP1nuD,DP1nuD); 
        }

       VE(pbhat,j)=VE(pbhat,j)/(pgeno); 
       scl_vec_mult(1/(pgeno),dph,dph); 
       replace_row(cX,j,dph); 

      if (it==*Nit-1) {
	if (KMc[j]<0.00001) vec_zeros(dp); else scl_vec_mult(1/KMc[j],dp,dp); 
	scl_vec_mult(VE(Y,j),dp,dp); vec_add(dp,qs,qs); 

        if (*fixhaplo==0) {
           for (i=0;i<*dimxih;i++)
           for(k=0;k<*hapdim;k++)  ME(DUeta[s],i,k)+=VE(dph,i)*
	      (VE(DP1nuD,k)-VE(pbhat,j)*VE(DPGnuD,k)/pgeno ); 
        }
      
      }

      if (KMc[j]<0.001) VE(Y,j)=(VE(Y,j)/0.001)-VE(pbhat,j); 
	      else VE(Y,j)=(VE(Y,j)/KMc[j])-VE(pbhat,j);

   } // j in 1...antpers

    if (0<*designtest && it==0 && s==0) { // print test {{{
       Rprintf(" Designmatrix first iteration step \n"); 
       print_mat(cX); 
    } // }}}

    totrisk=vec_sum(riskV); MtA(cX,cX,A); invert(A,AI); sing=0; 

    // print_mat(A); 

    for (i=0;i<*px;i++) if (fabs(ME(AI,i,i))<.0000001) {sing=1;}

    if (sing==1) {Rprintf(" non-invertible design time %lf\n",time); 
      it=*Nit-1;  
      for (c=0;c<ps;c++) VE(beta,c)=betaS[c]; 
      for (c=0;c<ps;c++) VE(bet1,c)=betaS[c]; 
    }
    if (sing==0) {
      vM(cX,Y,SCORE); 
      Mv(AI,SCORE,difbeta); vec_add(beta,difbeta,beta); 

//	 print_vec(Y); 
 // print_vec(SCORE); print_vec(difbeta); 

      for (i=0;i<ps;i++) VE(bet1,i)=VE(beta,i); 

      sumscore=0; 
      for (k=0;k<ps;k++) sumscore=sumscore+fabs(VE(difbeta,k)); 
      if ((sumscore<0.000001) & (it<*Nit-2)) it=*Nit-2;

      if (isnan(vec_sum(SCORE))) {
	Rprintf("missing values in SCORE %ld \n",(long int) s); 
	for (i=0;i<ps;i++) VE(beta,i)=-99; sim[0]=0;
	it=*Nit-1; 
	for (c=0;c<ps;c++) VE(beta,c)=betaS[c]; 
	for (c=0;c<ps;c++) VE(bet1,c)=betaS[c]; 
	break; }
    }


  if (*detail==1) { 
    Rprintf(" s er %ld,Estimate beta \n",(long int) s); print_vec(beta); 
    Rprintf("Score D l\n"); print_vec(difbeta); 
    Rprintf("Information -D^2 l\n"); print_mat(AI); };

    if (it==*Nit-1) scl_vec_mult(1/totrisk,qs,qs); 

} /* it */

	vec_zeros(VdB); mat_zeros(VAR); 
      for (j=0;j<*antclust;j++) {vec_zeros(cumA[j]);vec_zeros(cumhatA[j]);}
      for (i=0;i<*n;i++) { 
        j=clusters[i]; 
	extract_row(cX,i,dp); scl_vec_mult(VE(Y,i),dp,dp); 
	vec_add(dp,cumA[j],cumA[j]); 

	if (*fixhaplo==0) { // also correct for uncertainty in haplo-pars
	    // if (j==0) { Rprintf("==============\n"); print_mat(DUeta[s]);}
           Mv(DUeta[s],hapiid[j],xih); 
	   vec_subtr(cumA[j],xih,cumA[j]); }

	if ((time==x[i]) & (delta[i]==0)) vec_add(qs,cumhatA[j],cumhatA[j]);  

	}

        for (j=0;j<*antclust;j++) { 

	  vec_add(cumhatA[j],cumA[j],dp1); 
	  Mv(AI,dp1,dp2); replace_row(cumAt[j],s,dp2);  
	if (s<-1) print_vec(dp2); 
	  for(k=0;k<ps;k++) 
	    for(c=0;c<ps;c++) ME(VAR,k,c)=ME(VAR,k,c)+VE(dp2,k)*VE(dp2,c); 
	  if (*resample==1) 
	    for (c=0;c<ps;c++) {l=j*(ps)+c; biid[l*(*Ntimes)+s]=VE(dp2,c);}
      }
      for (i=1;i<ps+1;i++) {
	  var[i*(*Ntimes)+s]=ME(VAR,i-1,i-1); 
	  est[i*(*Ntimes)+s]=VE(beta,i-1); score[i*(*Ntimes)+s]=VE(SCORE,i-1); }


} /* s=1 ... *Ntimes */ 


    if (*sim==1) 
      comptestfunc(times,Ntimes,pps,est,var,vcudif,antsim,test,testOBS,Ut,
		   simUt,cumAt,weighted,antclust,gamma2,line); 
}   else {
    semihaplo(times,Ntimes,x,delta,cause,KMc,z,n,px,Nit,
      score,hess,est,var,sim,antsim,rani,test,testOBS,Ut,simUt,weighted,
      gamma,vargamma,semi,zsem,pg,trans,gamma2,CA,line,detail,biid,
      gamiid,resample,timepow,clusters,antclust,
      haplofreq,alphaiid,hapdim,
      nph,oh,nphpp,designfuncX,designfuncZ,rhoR,dimxih,dimzih,haplodes,
      fixhaplo,designtest); 
 
   }

if (*semi==0) { // {{{ free matrices 
   free_mats(&X,&cX,&DnuP1,&A,&AI,&VAR,&DUetaD,&DUetaP,&hapDes,NULL); 

   free_vecs(&rr,&bhatub,&riskV,&W,&Y,&Gc,&DELTA,&CAUSE,&bhat,&pbhat
   ,&xi,&rowX,&bet1,&dph,&xih1,&xih,&rowXh,
   &dp,&dp1,&dp2,&betaub,&VdB,&qs,&SCORE,&beta,&difbeta,
   &DPGnuD,&DPHnuD,&DP1nuD,NULL);

    for (i=0;i<*antclust;i++) {
      free_vec(hapiid[i]); 
      free_vec(cumhatA[i]); free_vec(cumA[i]); 
      free_mat(cumAt[i]);}

  for (j=0;j<*Ntimes;j++) { free_mat(DUeta[j]); }
 } // }}}
  free(vcudif); 
}


void semihaplo(
// {{{
times,Ntimes,x,delta,cause,KMc,z,antpers,px,Nit,
score,hess,est,var,sim,antsim,rani,test,testOBS,Ut,simUt,weighted,
gamma,vargamma,semi,zsem,pg,trans,gamma2,CA,line,detail,biid,
gamiid,resample,timepow,clusters,antclust,
haplofreq,alphaiid,hapdim,
nph,oh,nphpp,designfuncX,designfuncZ,rhoR,dimxih,dimzih,haplodes,fixhaplo,
designtest)
double *times,*x,*KMc,*z,*score,*hess,*est,*var,*test,*testOBS,
*Ut,*simUt,*gamma,*zsem,*vargamma,*gamma2,*biid,*gamiid,*timepow,
*haplofreq,*alphaiid,*haplodes;
int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*sim,*antsim,*rani,*weighted,
*semi,*pg,*trans,*CA,*line,*detail,*resample,*clusters,*antclust,
*nph,*oh,*nphpp,*hapdim,*dimxih,*dimzih,*fixhaplo,*designtest;
SEXP designfuncX,designfuncZ,rhoR; 
// }}}

{
 // {{{ variable defs 
  matrix *ldesignX,*A,*AI,*cdesignX,*ldesignG,*cdesignG;
  matrix *S,*dCGam,*CGam,*ICGam,*VarKorG,*dC,*XZ,*ZZ,*ZZI,*XZAI; 
  matrix *Ct,*C[*Ntimes],*Acorb[*Ntimes],*tmpM2,*tmpM3,*tmpM4; 
  matrix *Vargam,*dVargam,*M1M2[*Ntimes],*dM1M2,*M1M2t,*RobVargam;
  matrix *DUeta[*Ntimes],*DUeta1,*DUgamma,*DUgamma1,*hapDes;
  matrix *W3t[*antclust],*W4t[*antclust]; 
  vector *W2[*antclust],*W3[*antclust];
  vector *diag,*dB,*dN,*VdB,*AIXdN,*AIXlamt,*bhatt,*pbhat,*plamt;
  vector *korG,*pghat,*rowG,*gam,*dgam,*ZGdN,*IZGdN,*ZGlamt,*IZGlamt;
  vector *covsx,*covsz,*Y,*rr,*bhatub,*ziI,*xiI,*xi,*rowX,*rowZ,*difX,*zi,*z1,
         *tmpv1,*tmpv2,*lrisk;
  vector *xih,*zih,*hapiid[*antclust]; 
  vector *DPGnuD,*DPHnuD,*DP1nuD;
  int nb[1],itt,i,j,k,l,s,c,pmax,coef[1],c1D,iD,
      totrisk,ps[1],n[1],retur[1],haplotypeD[2]; 
  double ph,pgeno,time,dummy,dtime,timem,phD,
	 *vcudif=calloc((*Ntimes)*(*dimxih+1), sizeof(double)), 
	 *inc=calloc((*Ntimes)*(*dimxih+1), sizeof(double)); 
  double lrr,fabs(),pow(); 
  long robust[1],idum,fixedcov; 
  void FhaplodesMM(),comptestfunc(),fDUeta(); 

  idum=*rani; robust[0]=1; fixedcov=1; 
  n[0]=antpers[0]; retur[0]=0; nb[0]=Ntimes[0]; timem=0; 
  if (*trans==1) for (j=0;j<*dimzih;j++) if (fabs(timepow[j]-1)>0.0001) {timem=1;break;}
  if (*trans==2) for (j=0;j<*dimzih;j++) if (fabs(timepow[j])>0.0001) {timem=1;break;}
  // }}}
  

// {{{ matrix allocation
  for (j=0;j<*antclust;j++) { malloc_mat(*Ntimes,*dimxih,W3t[j]);
    malloc_vec(*hapdim,hapiid[j]); 
    malloc_mat(*Ntimes,*dimxih,W4t[j]); malloc_vec(*dimzih,W2[j]); 
    malloc_vec(*dimxih,W3[j]);
  }

  malloc_mats(*antpers,*px,&ldesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,NULL); 
  malloc_mats(*antpers,*dimxih,&cdesignX,NULL);
  malloc_mats(*antpers,*dimzih,&cdesignG,NULL); 
  malloc_mats(*dimxih,*dimxih,&A,&AI,NULL);
  malloc_mats(*dimzih,*dimzih,&dVargam,&Vargam,&RobVargam,&tmpM2,&ZZ,&VarKorG,
	                      &ICGam,&CGam,&dCGam,&S,&ZZI,NULL); 
  malloc_mats(*dimxih,*dimzih,&XZAI,&tmpM3,&Ct,&dC,&XZ,&dM1M2,&M1M2t,&tmpM4,NULL);
  for (j=0;j<*Ntimes;j++) { malloc_mat(*dimzih,*dimxih,Acorb[j]); 
    malloc_mat(*dimxih,*dimzih,C[j]); malloc_mat(*dimxih,*dimzih,M1M2[j]);
    malloc_mat(*dimxih,*hapdim,DUeta[j]); }
  malloc_mat(*dimzih,*hapdim,DUgamma); 
  malloc_mat(*dimxih,*hapdim,DUeta1); 
  malloc_mats(*dimzih,*hapdim,&DUgamma1,NULL); 

  malloc_vec(*px,xiI); malloc_vec(*pg,ziI); 
  malloc_vecs(*dimxih,&xi,&covsx,&rowX,&difX,&tmpv1,&korG,&diag,&dB,&VdB,&AIXdN,&AIXlamt,&bhatt,&xih,NULL);
  malloc_vecs(*dimzih,&zih,&covsz,&rowZ,&tmpv2,&zi,&z1,&rowG,&gam,&dgam,&ZGdN,&IZGdN,&ZGlamt,&IZGlamt,NULL);
  malloc_vecs(*antpers,&Y,&bhatub,&rr,&lrisk,&dN,&pbhat,&pghat,&plamt,NULL);
  malloc_mat(*nph-1,*hapdim,hapDes); 
  malloc_vecs(*hapdim,&DPGnuD,&DPHnuD,&DP1nuD,NULL);
// }}}

coef[0]=1; ps[0]=*dimxih; 
if (*px>=*pg) pmax=*px; else pmax=*pg; 
for (j=0;j<*dimzih;j++) VE(gam,j)=gamma[j]; 

if (fixedcov==1) { 
	for (c=0;c<*antpers;c++) { 
	for(j=0;j<pmax;j++)  { 
	   if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c]; 
	   if (j<*pg) ME(ldesignG,c,j)=zsem[j*(*antpers)+c]; } } 
} 

  // reading in basic matrices X and Z with observed covariates  {{{
  // making haplobased design for all subject and the relevant haplotypes
  
  int nallH=0; 
  for (i=0;i<*n;i++) nallH+=nphpp[i]; 
  matrix *XallH, *ZallH;
  double phallH[nallH]; 
  int ih,first,indexpersallH[*n],indexpersnph[*n]; 

  malloc_mat(nallH,*dimxih,XallH); 
  malloc_mat(nallH,*dimzih,ZallH); 

  // head_matrix(XallH); head_matrix(ZallH); 
  // print_matrix(XallH); print_matrix(ZallH); 
  
  k=0;c1D=0; 
  for (i=0;i<*n;i++) { 
   indexpersnph[i]=c1D; first=0; 
   extract_row(ldesignX,i,xiI); 
   extract_row(ldesignG,i,ziI); 

    for (j=0;j<nphpp[i];j++) {
       haplotypeD[0]=oh[c1D]; haplotypeD[1]=oh[c1D+1]; 
       ph=haplofreq[oh[c1D]]*haplofreq[oh[c1D+1]]; 
       phallH[k]=ph; 

       FhaplodesMM(xiI,ziI,haplotypeD,xih,designfuncX,rhoR);; 
       FhaplodesMM(xiI,ziI,haplotypeD,zih,designfuncZ,rhoR); 

       if (*detail>=3) { // print test {{{
	  Rprintf("======Design =================== \n"); 
          Rprintf(" person %d  \n",i); 
          Rprintf(" haplotypes hap1 hap2 %d %d  \n",oh[c1D],oh[c1D+1]);  
          Rprintf(" design, x,z \n");  
          print_vec(xiI);  print_vec(ziI); 
          Rprintf(" design, x(h),z(h) \n");  
	  print_vec(xih); print_vec(zih);
          ph=haplofreq[oh[c1D]]*haplofreq[oh[c1D+1]];
          Rprintf("probability for haplotype pair Donor %lf %d %d \n",ph,oh[c1D],oh[c1D+1]); 
       } // }}}
       
       replace_row(XallH,k,xih); 
       replace_row(ZallH,k,zih); 
       if (first==0) {indexpersallH[i]=k; first=1;}
       c1D+=2; k+=1; 
    }
   // Rprintf(" %d %d %d \n",i,indexpersallH[i],indexpersnph[i]); 
  }
  // }}}
  
 if (2<*designtest) { // print test {{{
  Rprintf(" Original designmatrices (x,z) \n"); 
  print_mat(ldesignX); print_mat(ldesignG); 
 } // }}}

// {{{ initialisering and reading of starting values 
   if (*fixhaplo==0) {
    for (i=0;i<*antclust;i++) {
    for (j=0;j<*hapdim;j++) VE(hapiid[i],j)=alphaiid[j*(*n)+i]; }
  
  for (i=0;i<*nph-1;i++)
  for(j=0;j<*hapdim;j++)  ME(hapDes,i,j)=haplodes[j*(*nph-1)+i];
   }
 // }}}

  for (itt=0;itt<*Nit;itt++)
    {
   R_CheckUserInterrupt();
      mat_zeros(Ct); mat_zeros(CGam); vec_zeros(IZGdN); vec_zeros(IZGlamt); 

      for (s=0;s<*Ntimes;s++)
      {
	  time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 

	  for(j=1;j<=*dimxih;j++) VE(bhatt,j-1)=est[j*(*Ntimes)+s];

	  totrisk=0; c1D=0; 
	  for (j=0;j<*antpers;j++) { 
            ih=indexpersallH[j]; 
	    VE(lrisk,j)=(x[j]>=time); totrisk=totrisk+VE(lrisk,j);

	  // initializing vectors for subject sums // {{{
	    pgeno=0; 
	    VE(plamt,j)=0; VE(pbhat,j)=0; VE(pghat,j)=0; 
	    vec_zeros(xi); vec_zeros(zi); 
	    vec_zeros(DPGnuD); vec_zeros(DP1nuD); 
	   // }}}
	  

       for (iD=0;iD<nphpp[j];iD++) {
          haplotypeD[0]=oh[c1D]; haplotypeD[1]=oh[c1D+1]; 
          phD=haplofreq[oh[c1D]]*haplofreq[oh[c1D+1]];
          pgeno=pgeno+phD; 

          if (itt==*Nit-1 && *fixhaplo==0) {// computes derivative of P(H_D=h_D) 
             fDUeta(c1D,oh,nph,hapdim,haplofreq,hapDes,DPHnuD); 
             vec_add(DPHnuD,DPGnuD,DPGnuD);
          }

           ph=phD; 
           extract_row(XallH,ih,xih);  
           extract_row(ZallH,ih,zih);  
           vec_star(xih,bhatt,rowX); VE(pbhat,j)=vec_sum(rowX); 
           vec_star(zih,gam,rowZ);   VE(pghat,j)=vec_sum(rowZ); 

	   if (1==*designtest && itt==0 && s==0) { // print test {{{
	      Rprintf("timepoint=%d person=%d #haplotype-pairs=%d \n",s,j,nphpp[j]); 
	      Rprintf("Donor haplotype1=%d haplotype2=%d  \n",oh[c1D],oh[c1D+1]);  
	      Rprintf("input x and f(x,h) design depending on haplotype \n"); 
	      print_vec(xiI); print_vec(xih);
	      Rprintf("input z and f(z,h) design depending on haplotype \n"); 
	      print_vec(ziI); print_vec(zih);
	   }
	   // print test }}}

	    lrr=0; 
	    if (*trans==1) { // {{{ computes P1 and D_P1
	      if (timem>0) for (l=0;l<*dimzih;l++)
		  lrr=lrr+VE(gam,l)*VE(zih,l)*pow(time,timepow[l]); 
	      else lrr=time*VE(pghat,j);   
	      VE(plamt,j)+=phD*(1-exp(-VE(pbhat,j)-lrr)); 
              vec_add_mult(xi,xih,ph*exp(-VE(pbhat,j)-lrr),xi); 
              scl_vec_mult(ph*exp(-VE(pbhat,j)-lrr),zih,zih); 
	      if (timem>0) for(l=0;l<*dimzih;l++) VE(zih,l)=
		  pow(time,timepow[l])*VE(zih,l); else scl_vec_mult(time,zih,zih); 
              vec_add(zi,zih,zi); 
	    }
	    if (*trans==2) { 
	       if (timem>0) for (l=0;l<*dimzih;l++) 
		       lrr=lrr+VE(gam,l)*VE(zih,l)*pow(time,timepow[l]); 
	      else lrr=VE(pghat,j);  
	      VE(rr,j)=exp(lrr);  
	      VE(plamt,j)+=phD*(1-exp(-exp(VE(pbhat,j))*VE(rr,j))); 
              vec_add_mult(xi,xih,phD*exp(-exp(VE(pbhat,j))*VE(rr,j))*
			          exp(VE(pbhat,j))*VE(rr,j),xi); 
              scl_vec_mult(phD*exp(-exp(VE(pbhat,j))*VE(rr,j))*
			           exp(VE(pbhat,j))*VE(rr,j),zih,zih); 
	      if (timem>0) for (l=0;l<*dimzih;l++) VE(zih,l)=
		               pow(time,timepow[l])*VE(zih,l); 
              vec_add(zi,zih,zi); 
	    } // }}}

       if (itt==*Nit-1 && *fixhaplo==0){// {{{computes derivative of P(H_D=h_D) 
          if (*trans==1) 
	    vec_add_mult(DP1nuD,DPHnuD,(1-exp(-VE(pbhat,j)-lrr)),DP1nuD); 
          if (*trans==2) 
	    vec_add_mult(DP1nuD,DPHnuD,(1-exp(-exp(VE(pbhat,j))*VE(rr,j))),DP1nuD); 
       } // }}}

	   c1D+=2;  
	   ih+=1; 
	  } // jD=0 ... nphpp[j] 

       if (itt==*Nit-1 && *fixhaplo==0) {// computes derivative of DUnu 
          scl_vec_mult(1/(pgeno),DP1nuD,DP1nuD); 
        }

    VE(plamt,j)= VE(plamt,j)/(pgeno); 
    scl_vec_mult(1/(pgeno),xi,xi); 
    scl_vec_mult(1/(pgeno),zi,zi); 

    replace_row(cdesignX,j,xi); replace_row(cdesignG,j,zi); 

     VE(Y,j)=((x[j]<=time) & (cause[j]==*CA))*1;
    if (KMc[j]<0.001) VE(Y,j)=(VE(Y,j)/0.001)-VE(plamt,j); 
    else VE(Y,j)=(VE(Y,j)/KMc[j])-VE(plamt,j);

    if (itt==*Nit-1 && *fixhaplo==0) {// DUeta[s] and DUgamma 
       for (i=0;i<*dimxih;i++)
       for(k=0;k<*hapdim;k++) ME(DUeta[s],i,k)+=VE(xi,i)*
       (VE(DP1nuD,k)-VE(plamt,j)*VE(DPGnuD,k)/pgeno ); 

       for (i=0;i<*dimzih;i++)
       for(k=0;k<*hapdim;k++) ME(DUgamma,i,k)+=dtime*VE(zi,i)*
             (VE(DP1nuD,k)-VE(plamt,j)*VE(DPGnuD,k)/pgeno); 
    }
   } // j in 1..antpers 

	 if (0<*designtest && itt==0 && s==0) { // print test {{{
          Rprintf(" Designmatrices first iteration step \n"); 
	  print_mat(cdesignX); print_mat(cdesignG); 
	 } // }}}

	  MtA(cdesignX,cdesignX,A); invert(A,AI); 
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

	  if (itt==*Nit-1 && *fixhaplo==0) {
	     // print_mat(DUgamma); 
             MtA(XZAI,DUeta[s],DUgamma1); 
	     scl_mat_mult(dtime,DUgamma1,DUgamma1); 
	     mat_subtr(DUgamma,DUgamma1,DUgamma); 
	     //print_mat(DUgamma); 
          
	     MxA(AI,DUeta[s],DUeta1); 
	     // MxA(AI,DUeta[s],DUeta[s]); 
	     //print_mat(DUeta1); 
	     mat_copy(DUeta1,DUeta[s]); 
	     //print_mat(DUeta[s]); 
	  }
	  /* scl_mat_mult(dtime,XZAI,tmpM4);mat_add(tmpM4,Ct,Ct); */

	  for (k=1;k<=*dimxih;k++) inc[k*(*Ntimes)+s]=VE(AIXdN,k-1); 

	  if (itt==*Nit-1) {
	    for (i=0;i<*antpers;i++) 
            { 
              j=clusters[i]; 	
	      extract_row(cdesignX,i,xi);
	      scl_vec_mult(VE(Y,i),xi,xi); Mv(AI,xi,rowX);
	      extract_row(cdesignG,i,zi); scl_vec_mult(VE(Y,i),zi,zi); 
	      vM(C[s],rowX,tmpv2); vec_subtr(zi,tmpv2,rowZ); 
	      scl_vec_mult(dtime,rowZ,rowZ); 
	      vec_add(rowZ,W2[j],W2[j]); 
	      for (k=0;k<*dimxih;k++) ME(W3t[j],s,k)= ME(W3t[j],s,k)+VE(rowX,k); 
	    }
          }
	} /* s=1,...Ntimes */

      invert(CGam,ICGam); Mv(ICGam,IZGdN,dgam); vec_add(gam,dgam,gam); 

      if (isnan(vec_sum(dgam))) {Rprintf("missing values in dgam %d \n",s);
	vec_zeros(gam); }

      dummy=0; for (k=0;k<*dimzih;k++)  dummy=dummy+fabs(VE(dgam,k)); 


      for (s=0;s<*Ntimes;s++) {
	vM(Acorb[s],dgam,korG); 
	est[s]=times[s]; var[s]=times[s]; 
	for (k=1;k<=*dimxih;k++)  { est[k*(*Ntimes)+s]=
            est[k*(*Ntimes)+s]+inc[k*(*Ntimes)+s]-VE(korG,k-1); 
	  dummy=dummy+fabs(inc[k*(*Ntimes)+s]-VE(korG,k-1)); 
	  /* Rprintf(" %lf ",est[k*(*Ntimes)+s]); printf(" \n");*/ }
      } /* s=1,...Ntimes */
      if ((dummy<0.000001) & (itt<*Nit-2)) itt=*Nit-2; 

      if (*detail==1) { 
	Rprintf(" iteration %d %d \n",itt,*Nit); 
	Rprintf("Total score %lf \n",dummy); 
	Rprintf("gamma change \n"); print_vec(dgam); }


    } /*itt løkke */ 

    MxA(ICGam,DUgamma,DUgamma1); 

    for (i=0;i<*antclust;i++) {Mv(ICGam,W2[i],zi); scl_vec_mult(1,zi,W2[i]); 
//	    Rprintf("==================================== \n"); 
//	    Rprintf(" %d \n",i); print_vec(W2[i]);

   if (*fixhaplo==0) { // also correct for uncertainty in haplo-pars
           Mv(DUgamma1,hapiid[i],zi); // print_vec(zi); 
	   vec_subtr(W2[i],zi,W2[i]); }
    }



  /* ROBUST VARIANCES   */ 
  if (*robust==1) 
  {
      for (s=0;s<*Ntimes;s++) {
	vec_zeros(VdB); 
	 for (i=0;i<*antclust;i++) {

	vM(Acorb[s],W2[i],rowX); extract_row(W3t[i],s,tmpv1); 
  // scl_vec_mult(0,rowX,rowX);   
	vec_subtr(tmpv1,rowX,difX); 

	if (*fixhaplo==0) { // also correct for uncertainty in haplo-pars
           Mv(DUeta[s],hapiid[i],xi); vec_subtr(difX,xi,difX); }
	replace_row(W4t[i],s,difX); 
	vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	if (*resample==1) {
	  if (s==1)
	    for (c=0;c<*dimzih;c++)
	    gamiid[c*(*antclust)+i]=VE(W2[i],c);
	    for (c=0;c<*dimxih;c++) {l=i*(*dimxih)+c; 
	    biid[l*(*Ntimes)+s]=biid[l*(*Ntimes)+s]+VE(difX,c);} }


	if (s==0) {for (j=0;j<*dimzih;j++) for (k=0;k<*dimzih;k++) 
	           ME(RobVargam,j,k)=ME(RobVargam,j,k)+VE(W2[i],j)*VE(W2[i],k);} 
	}  /* for (i=0;i<*antclust;i++) */ 
	for (k=1;k<*dimxih+1;k++) var[k*(*Ntimes)+s]=VE(VdB,k-1); 

      } /* s=0..Ntimes*/
    }

  /* MxA(RobVargam,ICGam,tmpM2); MxA(ICGam,tmpM2,RobVargam);*/

  for (j=0;j<*dimzih;j++) {gamma[j]=VE(gam,j);
    for (k=0;k<*dimzih;k++) {vargamma[k*(*dimzih)+j]=ME(RobVargam,j,k);}}

  if (*sim==1) {
  comptestfunc(times,Ntimes,ps,est,var,vcudif,antsim,test,testOBS,Ut,
	       simUt,W4t,weighted,antclust,gamma2,line);
  }

// {{{ freeing matrices 
  for (j=0;j<*antclust;j++) { free_mat(W3t[j]);
    free_vec(hapiid[j]); free_mat(W4t[j]); free_vec(W2[j]); free_vec(W3[j]); }

  for (j=0;j<*Ntimes;j++) { free_mat(Acorb[j]); 
    free_mat(C[j]); free_mat(M1M2[j]); free_mat(DUeta[j]); }

  free_mats(&ldesignX, &cdesignX, &ldesignG,&cdesignG, &A,&AI,
  &dVargam,&Vargam,&RobVargam,&tmpM2,&ZZ,&VarKorG,&ICGam,&CGam,&dCGam,&S,&ZZI,
  &XZAI,&tmpM3,&Ct,&dC,&XZ,&dM1M2,&M1M2t,&tmpM4, &hapDes, &DUgamma,
  &DUeta1,&DUgamma1,&XallH,&ZallH,NULL); 

  free_vecs(&xiI,&ziI,
  &xi,&covsx,&rowX,&difX,&tmpv1,&korG,&diag,&dB,&VdB,&AIXdN,&AIXlamt,&bhatt,&xih,
  &zih,&covsz,&rowZ,&tmpv2,&zi,&z1,&rowG,&gam,&dgam,&ZGdN,&IZGdN,&ZGlamt,
  &IZGlamt,&Y,&bhatub,&rr,&lrisk,&dN,&pbhat,&pghat,&plamt,
  &DPGnuD,&DPHnuD,&DP1nuD,NULL);
free(vcudif); free(inc); 
// }}}

}

void fDUeta(int c2,int *oh,int *nph,int *hapdim,
	    double *haplofreq,matrix *hapDes,vector *etascore)
// {{{
{
vector *tempScore; 
int k,loopThroughAll,nphm1=*nph-1; 
double tempSum,tempArray[nphm1]; 

 for(k = 0; k < nphm1; k++) tempArray[k]=0.0; 

  malloc_vecs(nphm1,&tempScore,NULL); 

 // calculates score for haplotype h(j) for subject i  {{{
        vec_zeros(tempScore); 
       if(oh[c2] != nphm1 && oh[c2+1] != nphm1){
          loopThroughAll = 0;
          tempArray[oh[c2]] += haplofreq[oh[c2+1]];
          tempArray[oh[c2+1]] += haplofreq[oh[c2]];
       } else if (oh[c2] != nphm1 && oh[c2+1] == nphm1){
          loopThroughAll = 1;
          tempArray[oh[c2]] += haplofreq[oh[c2+1]];
          for(k = 0; k < nphm1; k++){
             tempArray[k] -= haplofreq[oh[c2]];
          }
       } else if (oh[c2] == nphm1 && oh[c2+1] != nphm1){
          loopThroughAll = 1;
          tempArray[oh[c2+1]] += haplofreq[oh[c2]];
          for(k = 0; k < nphm1; k++){                   
             tempArray[k] -= haplofreq[oh[c2+1]];         
          }
       } else {
         loopThroughAll = 1;
         for(k = 0; k < nphm1; k++){
          tempArray[k] -= 2.0*haplofreq[nphm1];        
         }
       }
       if(loopThroughAll == 1){
         tempSum = 0;
         for(k = 0; k < nphm1; k++){
          tempSum += tempArray[k]*haplofreq[k];
         }
       } else {
       tempSum = tempArray[oh[c2]] * haplofreq[oh[c2]] + tempArray[oh[c2+1]] * haplofreq[oh[c2+1]];
       }

       for(k = 0; k < nphm1; k++){
         VE(tempScore,k) = haplofreq[k]*(tempArray[k] - tempSum);
       }

//for(k = 0; k < nphm1; k++) Rprintf(" %lf %lf -----\n",tempArray[k],haplofreq[k]); 

     tempSum=vec_sum(tempScore);  
     vM(hapDes,tempScore,etascore); 
     // }}}

     free_vecs(&tempScore,NULL); 
}
// }}}

/* defined in haplo-surv
SEXP mkansv(double *x,int dim) 
// {{{
{ 
   int i;
   SEXP ans;
   PROTECT(ans = allocVector(REALSXP, dim));
   for(i=0;i<dim;i++) {REAL(ans)[i] = x[i]; }
   UNPROTECT(1);
   return ans;
} 
// }}}

void designeval(double *x, double *h, SEXP f, double *res, SEXP rho, int dimx,int dimres)
// {{{
{ 
   SEXP ans;
   int i;
   defineVar(install("x"), mkansv(x,dimx),rho);
   defineVar(install("h"), mkansv(h,2),rho);
   ans=eval(f,rho);
   PROTECT(ans);
   for(i=0;i<dimres;i++) { res[i]=REAL(ans)[i]; }
   UNPROTECT(1);
} 
// }}}
*/
void Fdesigneval(double *x, double *z,double *h, SEXP f, double *res, SEXP rho, int dimx,int dimz, int dimres)
// {{{
{ 
   SEXP ans;
   SEXP mkansv(); 
   int i;
   defineVar(install("z"), mkansv(z,dimz),rho);
   defineVar(install("x"), mkansv(x,dimx),rho);
   defineVar(install("h"), mkansv(h,2),rho);
   ans=eval(f,rho);
   PROTECT(ans);
   for(i=0;i<dimres;i++) { res[i]=REAL(ans)[i]; }
   UNPROTECT(1);
} 
// }}}

void Fhaplodes(vector *xi,int *haplotype,vector *xih,SEXP f,SEXP rho) 
 // {{{
{
   int i,j,dimx,dimxih; 
   dimx=length_vector(xi); 
   dimxih=length_vector(xih); 
   double x[dimx],h[2]; 
   double res[dimxih]; 
   void designeval(); 
   for (j=0;j<2;j++) h[j]=haplotype[j]*1.0;
   for (j=0;j<dimx;j++) x[j]=VE(xi,j);
   for (j=0;j<dimxih;j++) res[j]=0; 
   designeval(x, h, f , res,rho,dimx, dimxih);
   for(i=0;i<dimxih;i++) VE(xih,i)=res[i]; 
} 
// }}}


void FhaplodesMM(vector *xi,vector *zi,int *haplotype,vector *xih,SEXP f,SEXP rho)
 // {{{
{
   int i,j,dimx,dimxih,dimz; 
   dimx=length_vector(xi); 
   dimz=length_vector(zi); 
   dimxih=length_vector(xih); 
   double x[dimx],h[2],z[dimz]; 
   double res[dimxih]; 
   void Fdesigneval(); 
   for (j=0;j<2;j++) h[j]=haplotype[j]*1.0;
   for (j=0;j<dimx;j++) x[j]=VE(xi,j);
   for (j=0;j<dimz;j++) z[j]=VE(zi,j);
   for (j=0;j<dimxih;j++) res[j]=0; 
   Fdesigneval(x, z, h, f , res,rho,dimx, dimz, dimxih);
   for(i=0;i<dimxih;i++) VE(xih,i)=res[i]; 
} 
// }}}
