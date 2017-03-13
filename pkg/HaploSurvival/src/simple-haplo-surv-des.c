#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include <Rdefines.h>
#include <R.h>
#include "haplosurv.h" 
                 
void simplehaplosurvdes(
// {{{
times,Ntimes,designX, 
nx,px,designG,
ng,pg,antpers,
start,stop,betaS,
Nit,cu,vcu,
loglike,Iinv,Vbeta,
detail,sim,antsim,
Rvcu,RVbeta,test,
testOBS,Ut,simUt,
Uit,wscore,idtimes,
status,score,dhatMit,
dhatMitiid,retur,sym,
nph,oh,nphpp,
haplopars,fixbeta,fixhaplofreq,
dimzih,dimxih,rho,
scoregeno,d2lgeno, survhaploscore,
twostage,varpar, d2score,
step,params,lm,minlm,designfuncX,designfuncZ,rhoR,haplodes,alpha,dimhap,
haplofreq,biid,gamiid,resample)
SEXP designfuncX,designfuncZ,rhoR; 
double *designX,*designG,*times,*betaS,*start,*stop,*cu,*loglike,*Vbeta,*RVbeta,
       *vcu,*Rvcu,*Iinv,*test,*testOBS,*Ut,*simUt,*Uit,*scoregeno,*dhatMit,
       *dhatMitiid,*haplopars,
       *survhaploscore,*rho,*d2lgeno,*score,*varpar,*d2score,*step,*params,*lm,*minlm,
       *haplodes,*alpha,*haplofreq,*biid,*gamiid; 
int *nx,*px,*ng,*pg,*antpers,*Ntimes,*Nit,*detail,*sim,*antsim,*status,
*retur,*sym,*fixbeta,*fixhaplofreq,
*nph,*oh,*dimzih,*dimxih,*nphpp,*twostage,*wscore,*dimhap,*resample,
	*idtimes;
// }}}
{

// {{{ memory allocation 1
  matrix *ldesignX,*Z,*ldesignG,*cdesX,*cdesX2,*CtVUCt,*A,*AI,*AA;
  matrix *dYI,*Ct,*dM1M2,*M1M2t,*COV,*ZX,*ZP,*hapDes; 
  matrix *tmp1,*tmp2,*tmp5,*tmp3,*dS,*S1hap,*S1f,*S1,*SI,*S2,*M1,*VU,*ZXAI,*VUI, *tmp6; 
  matrix *RobVbeta,*Delta,*tmpM1,*Utt,*Delta2,*tmpM2;
  matrix *St[*Ntimes],*M1M2[*Ntimes],*C[*Ntimes],*ZXAIs[*Ntimes],*dYIt[*Ntimes],*YIt[*Ntimes];
  matrix *dW3t[*antpers],*W3t[*antpers],*W4t[*antpers],*W2t[*antpers],*AIxit[*antpers],*Uti[*antpers],*tmp4,*Fst[(*Ntimes)*(*Ntimes)]; 
  matrix *dGt[*Ntimes],*dFt[*Ntimes],*S0tI[*Ntimes],*Ident,*gt[*Ntimes],*q2t[*Ntimes],*dG[*Ntimes],*q1t[*antpers]; 
  vector *dLamt[*antpers],*zcol,*Gbeta,*one,*difzzav; 
  vector *dA,*VdA,*MdA,*delta,*zav,*lamt,*plamt,*dlamt;
  vector *xi,*zi,*xih2,*xih1,*xih,*zih1,*zih2,*zih,*Uhap,*Vhaplopars,*Uf,*U,*beta; 
  vector *offset,*ZXdA[*Ntimes],*varUthat[*Ntimes],*Uprofile;
  vector *Lplamt,*ahatt; 
  vector *tmpv1,*tmpv2,*rowXh,*rowZbeta1,*rowZbeta2,*rowZbeta1dN,*rowZbeta2dN,*rowZh,*difX,*VdB,*lht; 
  vector *W2[*antpers],*W3[*antpers],*reszpbeta,*res1dim,*dAt[*Ntimes]; 
  vector *pars,*vtheta1,*vtheta2,*vtheta3; 
  vector *vhaplo1,*vhaplo2,*XtempScore,*tempScore,
	 *scorei[(*antpers)],*Xscorei[(*antpers)]; 
  int nphm1=nph[0]-1,pp=0,t,c,robust=1,i,j,k,l,s=0,it,pmax,
  // int pers=0,count; 
      *risk=calloc(*antpers,sizeof(int)), 
      *ipers=calloc(*Ntimes,sizeof(int)),nap; 
  int amount[1],haplotype[2],dimpar,fixmax,c1,c2; 
  double *weight=calloc(*antpers,sizeof(double)),dtime,time,ph,sph,RR,dummy,
	 ll,lle,llo,tael,naevn;
  double S0,tau,random,sumscore=0,LogLikegeno[1],ddsurvbeta,ddhazbeta; 
  double norm_rand(); 
  double surv,shaz,dsurvbeta,dsurvA,dhazbeta,dhazA; 
  double xihAt,maxdelt=10,taelpers=0,naevnpers=0; 
  // double zihbeta; 
  double *scoregenoC=calloc(nphm1,sizeof(double)),
	 *d2lgenoC=calloc(nphm1*nphm1,sizeof(double)); 
   void GetRNGstate(),PutRNGstate(),FhaplodesMM(),genoLogLikeHp(); 
  int nallH,loopThroughAll;
  double tempSum;
// }}}
  // set up dimension of score 
// {{{
  dimpar=0; nphm1=*nph-1; 
  if (*fixbeta==0) dimpar=*dimzih; 
  if (*fixhaplofreq==0) dimpar=dimpar+*dimhap;
  // if (*fixbeta==1 && *fixhaplofreq==1) dimpar=1; 

  //Rprintf(" %d %d %d %d \n",dimpar,*pg,*dimzih,*dimhap); 

  vector *rowZnu1,*rowZnu2,*rowZnu,*DnurowXh,*DbetarowXh; 
  vector *DbetarowXh1[*dimzih],*DbetarowXh2[*dimzih]; 
  vector *DbetaZ1[dimpar],*DbetaZ2[dimpar]; 
  vector *DnurowXh1[*dimhap],*DnurowXh2[*dimhap],*DArowXh,*DArowXh1[*dimxih],*DArowXh2; 
  matrix *DAX[*dimxih],*E[*dimxih],*DAS0[*dimxih],*ZDAX[*dimxih],*Smat[*dimxih]; 
  matrix *DthetaS0[dimpar],*DthetaZ[dimpar],*DthetaZY[dimpar]; 
  vector *DthetaZdN[dimpar]; 
  vector *rowZtheta1, *rowZtheta2, *rowZtheta1dN, *rowZtheta2dN; 
  for (j=0;j<*dimxih;j++) { malloc_mat(*antpers,*dimxih,DAX[j]); malloc_mat(*dimxih,*dimxih,E[j]); 
     malloc_mat(*dimxih,*dimxih,DAS0[j]); malloc_mat(dimpar,*dimxih,ZDAX[j]); 
     malloc_mat(dimpar,*dimxih,Smat[j]); 
  }
  for (j=0;j<dimpar;j++) { 
     malloc_mat(*dimxih,*dimxih,DthetaS0[j]); malloc_mat(*antpers,*dimzih,DthetaZ[j]);      
     malloc_mat(*dimzih,*dimxih,DthetaZY[j]); malloc_vec(*dimzih,DthetaZdN[j]);      
  }
  //}}}
  // matrix allocation 1  
  // {{{
  double tempArray[nphm1]; 
  // matrix allocation 2  
  malloc_mat(*dimhap,*nph-1,AA); 
  for (j=0;j<*antpers;j++) { 
    malloc_vec(*Ntimes,dLamt[j]); malloc_mat(*Ntimes,*dimxih,W3t[j]);
    malloc_mat(*Ntimes,*dimxih,dW3t[j]); malloc_mat(*Ntimes,*dimxih,W4t[j]);
    malloc_mat(*Ntimes,dimpar,W2t[j]); malloc_mat(*Ntimes,*dimzih,Uti[j]);
    malloc_vec(dimpar,W2[j]); malloc_vec(*dimxih,W3[j]);
    malloc_mat(*Ntimes,dimpar,q1t[j]); malloc_mat(*Ntimes,*dimxih,AIxit[j]);
  }
  malloc_mat(nphm1,*dimhap,hapDes); 
  malloc_mat(*Ntimes,*dimxih,Delta); malloc_mat(*Ntimes,*dimxih,tmpM1);
  malloc_mat(*Ntimes,*dimzih,Delta2); malloc_mat(*Ntimes,*dimzih,tmpM2);
  malloc_mat(*Ntimes,dimpar,Utt); malloc_vec(*Ntimes,lht); malloc_vec(1,reszpbeta); malloc_vec(1,res1dim);
  nap=floor(*antsim/50);
  malloc_mats(*antpers,*px,&ldesignX,NULL);
  malloc_mats(*antpers,*dimxih,&cdesX2,&cdesX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,NULL); 
  malloc_mats(*antpers,dimpar,&ZP,&Z,NULL); 
  malloc_mats(*dimxih,*dimxih,&tmp4,&Ident,&COV,&A,&AI,&M1,&CtVUCt,NULL); 
  malloc_mat(*nph-1,*nph-1,S1f); 
  malloc_mat(*dimhap,*dimhap,S1hap); 
  malloc_vecs(*nph-1,&Uf,&Vhaplopars,NULL); 
  malloc_vec(*dimhap,Uhap); 
  malloc_mats(dimpar,dimpar,&RobVbeta,&dS,&S1,&S2,&SI,&VU,&VUI,NULL); 
  malloc_mats(dimpar,*dimxih,&ZXAI,&tmp5,&tmp3,&ZX,&dM1M2,&M1M2t,NULL); 
  malloc_mats(*dimxih,dimpar,&dYI,&Ct,NULL); 
  malloc_mat(*dimxih,dimpar,tmp6);
  malloc_mats(dimpar,*dimxih,&tmp1,&tmp2,NULL);
  malloc_vecs(*antpers,&Lplamt,&dlamt,&plamt,&lamt,&zcol,&Gbeta,&one,&offset,NULL); 
  malloc_vecs(*dimxih,&DArowXh,&ahatt,&tmpv1,&difX,&VdB,&rowXh,&xih1,&xih2,&xih,&dA,&VdA,&MdA,NULL); 
  malloc_vecs(*dimzih,&tmpv2, &rowZbeta1,&rowZbeta2, &rowZbeta1dN,&rowZbeta2dN,
		  &rowZh,&zih1,&zih2,&zih,&beta,NULL); 
  malloc_vecs(dimpar,&pars,&U,&vtheta1,&vtheta2,&vtheta3,&delta,&zav,&difzzav,&Uprofile,NULL); 
  malloc_vecs(dimpar,&rowZtheta1dN,&rowZtheta2dN,&rowZtheta1,&rowZtheta2,NULL); 
  malloc_vec(*px,xi); malloc_vec(*pg,zi); 
  identity_matrix(Ident);
  for(j=0;j<*Ntimes;j++) {
    malloc_mat(*dimxih,*dimxih,dFt[j]); 
    malloc_mat(dimpar,*dimxih,gt[j]); malloc_mat(dimpar,*dimxih,dGt[j]);
    malloc_mat(dimpar,*dimxih,q2t[j]); 
    malloc_mat(*dimxih,*dimxih,S0tI[j]); malloc_mat(*dimxih,dimpar,dG[j]);
    malloc_mat(*dimxih,dimpar,C[j]); malloc_mat(dimpar,*dimxih,M1M2[j]);
    malloc_mat(dimpar,*dimxih,ZXAIs[j]); 
    malloc_mat(*dimxih,dimpar,YIt[j]); malloc_mat(*dimxih,dimpar,dYIt[j]);
    malloc_vec(*dimxih,dAt[j]); malloc_vec(dimpar,ZXdA[j]);
    malloc_mat(dimpar,dimpar,St[j]); malloc_vec(dimpar,varUthat[j]);
    for(i=0;i<*Ntimes;i++){ malloc_mat(*dimxih,*dimxih,Fst[j*(*Ntimes)+i]); }
  }
  for(j=0;j<*antpers;j++)  {malloc_vec(*nph-1,scorei[j]); 
	                    malloc_vec(*dimhap,Xscorei[j]);}
  // }}}
  
  // {{{ // matrix allocation 3  
  pmax=max(*px,*pg); ll=0; vec_ones(one); vec_ones(difX); 
  for(j=0;j<*dimzih;j++) VE(beta,j)=betaS[j]; 
  cu[0]=times[0]; 
  matrix *DbetaX[*dimzih],*DnuX[*dimhap],*DthetaX[dimpar];
  for (j=0;j<*dimzih;j++) { malloc_mat(*antpers,*dimxih,DbetaX[j]); }
  for (j=0;j<dimpar;j++) { malloc_mat(*antpers,*dimxih,DthetaX[j]); }
  for (j=0;j<*dimhap;j++) { malloc_mat(*antpers,*dimxih,DnuX[j]); }
  malloc_vecs(*nph-1,&tempScore,NULL); 
  malloc_vecs(*dimhap,&XtempScore,&rowZnu1,&rowZnu2,&rowZnu,&vhaplo1,&vhaplo2,NULL); 
  if (*fixbeta==0) {for (j=0;j<*dimzih;j++) VE(pars,j)=VE(beta,j);}
  if (*fixhaplofreq==0) {for (j=0;j<*dimhap;j++) VE(pars,(*fixbeta==0)*(*dimzih)+j)=alpha[j];;}
  fixmax=max(*fixhaplofreq,*fixbeta); 
  for(k=0;k<*dimzih;k++) { malloc_vec(*dimxih,DbetarowXh1[k]); malloc_vec(*dimxih,DbetarowXh2[k]); }
  for(k=0;k<dimpar;k++) { malloc_vec(*dimzih,DbetaZ1[k]); malloc_vec(*dimzih,DbetaZ2[k]); }
  for(k=0;k<*dimxih;k++) { malloc_vec(*dimxih,DArowXh1[k]); }
  malloc_vec(*dimxih,DArowXh2); 
  for(k=0;k<*dimhap;k++) { malloc_vec(*dimxih,DnurowXh1[k]); malloc_vec(*dimxih,DnurowXh2[k]); }
  malloc_vec(*dimxih,DnurowXh); malloc_vec(*dimxih,DbetarowXh); 
  // }}}

  // setting up haplotype freqs given haplopars */  {{{
  if (*fixhaplofreq==0) { 
  for (j=0;j<*nph-1;j++) { 
      haplofreq[j]=exp(haplopars[j]); 
      haplofreq[*nph-1]=haplofreq[*nph-1]+haplofreq[j]; }
      haplofreq[*nph-1]=1+haplofreq[*nph-1];
      for (j=0;j<*nph-1;j++) { haplofreq[j]=haplofreq[j]/haplofreq[*nph-1];} 
			        haplofreq[*nph-1]=1/haplofreq[*nph-1]; 
  for (i=0;i<nphm1;i++)
  for(j=0;j<*dimhap;j++)  ME(hapDes,i,j)=haplodes[j*(nphm1)+i];
  }


// }}}

  for (j=0;j<*dimzih;j++) { VE(beta,j)=betaS[j]; } 

  // reading in basic matrices X and Z with observed covariates  {{{
  for (c=0;c<*antpers;c++)
	  for(j=0;j<pmax;j++) {
	    if (j<*px) { ME(ldesignX,c,j)=designX[j*(*nx)+c]; }
	    if (j<*pg) { ME(ldesignG,c,j)=designG[j*(*ng)+c]; }
  }
  if (pp==2 && s==1) {print_mat(ldesignX); print_mat(ldesignG);}
       
  // making haplobased design for all subject and the relevant haplotypes
  //
  nallH=0; 
  for (i=0;i<*antpers;i++) nallH+=nphpp[i]; 
  matrix *XallH, *ZallH;
  double RRallH[nallH]; 
  // double Zbeta[nallH]; 
  int ih,first,indexpersallH[*antpers],indexpersnph[*antpers]; 
  vector *ZallHbeta; 

  malloc_mat(nallH,*dimxih,XallH); 
  malloc_mat(nallH,*dimzih,ZallH); 
  malloc_vec(nallH,ZallHbeta); 

  k=0;c1=0; 
  for (i=0;i<*antpers;i++) {
   indexpersnph[i]=c1; first=0; 
   extract_row(ldesignX,i,xi); extract_row(ldesignG,i,zi);  

    for (j=0;j<nphpp[i];j++) {
       haplotype[0]=oh[c1]; haplotype[1]=oh[c1+1]; 
       ph=haplofreq[oh[c1]]*haplofreq[oh[c1+1]]; 
       FhaplodesMM(xi,zi,haplotype,xih,designfuncX,rhoR); 
       FhaplodesMM(xi,zi,haplotype,zih,designfuncZ,rhoR); 

       if (*detail>=4) { // print test {{{
	  Rprintf("=======Design========================= \n"); 
          Rprintf(" %d %d %d %d \n",i,j,k,c1); 
          Rprintf(" hap1 hap2 %d %d  \n",oh[c1],oh[c1+1]);  
          print_vec(zi); print_vec(zih);
          print_vec(xi); print_vec(xih);
          ph=haplofreq[oh[c1]]*haplofreq[oh[c1+1]];
          Rprintf(" %lf %d %d \n",ph,oh[c1],oh[c1+1]); 
	  Rprintf("=============================== \n"); 
       }
       // print test }}}
       VE(ZallHbeta,k)=vec_prod(zih,beta);
       RRallH[k]=exp(VE(ZallHbeta,k)); 
       replace_row(XallH,k,xih); replace_row(ZallH,k,zih);
       if (first==0) {indexpersallH[i]=k; first=1;}
       c1+=2; k+=1; 
    }
    //Rprintf(" %d %d %d \n",i,indexpersallH[i],indexpersnph[i]); 
  }
  // }}}

  // head_matrix(XallH); head_matrix(ZallH); 
    //print_mat(XallH); print_mat(ZallH); 

  // Main procedure ================================== 
  for (it=0;it<*Nit;it++) {

    // initializing {{{ score and haplopars 
    vec_zeros(U); mat_zeros(S1);  lle=0; llo=0;  S0=0; mat_zeros(dS); 

    if (*fixbeta==0) {
	    for (j=0;j<*dimzih;j++) VE(beta,j)=VE(pars,j);
            if (it>=1) Mv(ZallH,beta,ZallHbeta); 
    }

    if (*fixhaplofreq==0) {
       for (j=0;j<*dimhap;j++) VE(Uhap,j)=VE(pars,(*fixbeta==0)*(*dimzih)+j);
       Mv(hapDes,Uhap,Vhaplopars); 
    haplofreq[*nph-1]=0; 
   for (j=0;j<*nph-1;j++) { 
	 haplopars[j]=VE(Vhaplopars,j);
	 haplofreq[j]=exp(haplopars[j]); 
	 haplofreq[*nph-1]=haplofreq[*nph-1]+haplofreq[j]; }
	 haplofreq[*nph-1]=1+haplofreq[*nph-1];
   for (j=0;j<*nph-1;j++) { 
	 haplofreq[j]=haplofreq[j]/haplofreq[*nph-1];} 
	 haplofreq[*nph-1]=1/haplofreq[*nph-1]; 
    }
    // }}}

  // genologlike  computation {{{
  if (*fixhaplofreq==0)  {
	  for (j=0;j<*nph-1;j++) { scoregenoC[j]=0; 
	    for (l=0;l<*nph-1;l++) d2lgenoC[j*(*nph-1)+l]=0; } 
    //genoLogLike(oh,nphpp,nph,antpers,haplofreq,rho,LogLikegeno,scoregeno,d2lgeno); 
    amount[0]=3; 
 //    Rprintf(" genolog \n"); 
    genoLogLikeHp(oh,nphpp,nph,amount,antpers,haplopars,rho,LogLikegeno,scoregenoC,scorei,d2lgenoC); 

   for (j=0;j<*nph-1;j++) {VE(Uf,j)=scoregenoC[j];}
       for (j=0;j<*nph-1;j++) for (l=0;l<*nph-1;l++) ME(S1f,j,l)=-d2lgenoC[j*(*nph-1)+l]; 

    vM(hapDes,Uf,Uhap); MtA(hapDes,S1f,AA); MxA(AA,hapDes,S1hap);

    if (*detail>=5) {  
     print_vec(Uf); print_vec(Uhap); print_mat(hapDes); 
     print_mat(S1f); 
     print_mat(S1hap); 
    }

   for (j=0;j<*dimhap;j++) {VE(U,*dimzih*(*fixbeta==0)+j)=VE(Uhap,j); 
   for (k=0;k<*dimhap;k++) ME(S1,*dimzih*(*fixbeta==0)+j,*dimzih*(*fixbeta==0)+k)=
	   ME(S1hap,j,k);  }
    } 

    if (*detail>=3) {  
       Rprintf(" Relative risk parameters beta \n"); print_vec(beta); 
       Rprintf("haplo-parameters \n"); 
       for (j=0;j<*dimhap;j++) Rprintf(" %lf ",VE(pars,(*fixbeta==0)*(*dimzih)+j)); 
       Rprintf(" \n"); 
       Rprintf("haplo-frequencies \n"); 
       for(k=0;k<*nph;k++)  Rprintf(" %lf ",haplofreq[k]); printf("\n"); 
       Rprintf("haplo-parameters \n"); 
       for(k=0;k<*nph-1;k++)  Rprintf(" %lf ",haplopars[k]); printf("\n"); 
       Rprintf("score geno \n"); 
       for (j=0;j<*dimhap;j++)  Rprintf(" %lf ",VE(U,j));  printf("\n"); 
       Rprintf("D2loglike geno  \n"); print_mat(S1); 
    }
 // }}}
 
  for (s=1;s<*Ntimes;s++){ // starts with 1 because times[0] is 0 (start.time)
       time=times[s]; 
       mat_zeros(cdesX); mat_zeros(Z); mat_zeros(cdesX2); 

   for (i=0;i<*dimxih;i++) VE(tmpv1,i)=cu[(i+1)*(*Ntimes)+s-1]; 

   // initializing  {{{
   for (k=0;k<*dimxih;k++) {mat_zeros(DAS0[k]);}
   for (k=0;k<dimpar;k++) {mat_zeros(DthetaS0[k]); mat_zeros(DthetaZY[k]);
     mat_zeros(DthetaZ[k]); vec_zeros(DthetaZdN[k]);      
   }
   vec_zeros(rowZtheta1dN); vec_zeros(rowZtheta2dN); 
   // }}}
   
    //Rprintf(" %d \n",idtimes[s-1]); 
   
    for (i=idtimes[s-1];i<*antpers;i++) {
         //    Rprintf(" %d \n",idtimes[s-1]); 
        ih=indexpersallH[i]; c1=indexpersnph[i]; 

      // Rprintf(" %d %d %d %d %d \n",s,i,j,c1,ih); 

   // initializing  {{{
       naevn=0; tael=0;  
       vec_zeros(rowXh); vec_zeros(rowZh); vec_zeros(rowZbeta1); vec_zeros(rowZbeta2); 
       vec_zeros(rowZtheta1); vec_zeros(rowZtheta2); 
       vec_zeros(rowZnu1); vec_zeros(rowZnu2); vec_zeros(rowZnu); 
       vec_zeros(DArowXh); vec_zeros(DArowXh2); vec_zeros(DnurowXh); vec_zeros(DbetarowXh); 
       for (k=0;k<dimpar;k++) { vec_zeros(DbetaZ1[k]); vec_zeros(DbetaZ2[k]); }
       for (k=0;k<*dimzih;k++) { vec_zeros(DbetarowXh1[k]); vec_zeros(DbetarowXh2[k]); }
       for (k=0;k<*dimhap;k++) { vec_zeros(DnurowXh1[k]); vec_zeros(DnurowXh2[k]); tempArray[k]=0; }
       for (k=0;k<*dimxih;k++) { vec_zeros(DArowXh1[k]); }
       ph=0; sph=0; 
       // }}}

     // extract_row(ldesignX,i,xi); extract_row(ldesignG,i,zi); 
     // Rprintf(" hej 00\n"); 

       for (j=0;j<nphpp[i];j++) {
           haplotype[0]=oh[c1]; haplotype[1]=oh[c1+1]; 
           c2=c1; 

     // calculates score for haplotype h(j) for subject i  {{{
     if (*fixhaplofreq==0 && *twostage==0) { 
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
     vM(hapDes,tempScore,XtempScore); 
     } 
     // }}}
     
   //
   extract_row(XallH,ih,xih);  
   extract_row(ZallH,ih,zih);  

   if ( s<=0) { // print test {{{
      Rprintf("============================== \n"); 
      Rprintf(" %d %d %d %d %d \n",s,i,j,c1,ih); 
      Rprintf(" hap1 hap2 %d %d  \n",oh[c1],oh[c1+1]);  
      print_vec(zih);
      print_vec(xih);
      ph=haplofreq[oh[c1]]*haplofreq[oh[c1+1]];
      Rprintf(" %lf %d %d \n",ph,oh[c1],oh[c1+1]); 
   }
   // print test }}}

   // basic quantities for intensity {{{
   RR=exp(VE(ZallHbeta,ih));  
   vec_star(xih,tmpv1,xih2); xihAt=vec_sum(xih2);  xihAt=max(0,xihAt); 
   ph=haplofreq[oh[c1]]*haplofreq[oh[c1+1]];
   sph=sph+ph; surv=exp(-xihAt*RR)*ph;
   if (surv<0) Rprintf(" it s i surv ph %d %d %d %lf %lf \n",it,s,i,surv,ph); 
   shaz=exp(-xihAt*RR)*RR*ph; naevn=naevn+surv; tael=tael+shaz; 
   vec_add_mult(rowXh,xih,shaz,rowXh); 
   if (pp==1) Rprintf("shaz surv ===== %d %d %d %lf %lf %lf %lf %lf  \n",s,i,j,shaz,surv,naevn,tael,xihAt); 
   if (pp==1) { Rprintf("%lf %lf \n",ph,RR);  print_vec(xih); print_vec(rowXh); }
   // }}}

   // matrices are set up {{{ 
   // part of the tilde Z matrix  with beta scores  
   dsurvbeta=(-1)*(xihAt*RR)*surv; 
   dhazbeta=RR*dsurvbeta+surv*RR;  
   ddsurvbeta=(-1)*(xihAt*RR)*surv+(-xihAt*RR)*dsurvbeta; 
   ddhazbeta=RR*dsurvbeta+RR*ddsurvbeta+dsurvbeta*RR+RR*surv;  
   if (*fixbeta==0 ) {
      vec_add_mult(rowZbeta1,zih,dhazbeta,rowZbeta1);   // D_beta T
      vec_add_mult(rowZbeta2,zih,dsurvbeta,rowZbeta2);  // D_beta N
      // for beta derivative of Z 
      for (k=0;k<*dimzih;k++) {
         vec_add_mult(DbetaZ1[k],zih,VE(zih,k)*ddhazbeta,DbetaZ1[k]); 
          vec_add_mult(DbetaZ2[k],zih,VE(zih,k)*ddsurvbeta,DbetaZ2[k]); 
      } 
   }

   // part of the tilde Z matrix  with nu scores 
   if (*fixhaplofreq==0 && *twostage==0 ) {
   if (pp==1) Rprintf("i j shaz surv  %d %d %lf %lf \n",i,j,shaz,surv); 
   if (pp==1) print_vec(XtempScore); 
        vec_add_mult(rowZnu1,XtempScore,shaz,rowZnu1);     // D_nu T 
        vec_add_mult(rowZnu2,XtempScore,surv,rowZnu2);     // D_nu N 
   } 

  // for second derivatives, for D_beta X and so forth    
  if (*fixbeta==0 ) {
  for (k=0;k<*dimzih;k++) {
     vec_add_mult(DbetarowXh1[k],xih,VE(zih,k)*dhazbeta,DbetarowXh1[k]);  
     vec_add_mult(DbetarowXh2[k],xih,VE(zih,k)*dsurvbeta,DbetarowXh2[k]); } 
   }
   dsurvA=-RR*surv; dhazA=RR*dsurvA; 
           for (k=0;k<*dimxih;k++) {
               vec_add_mult(DArowXh1[k],xih,VE(xih,k)*dhazA,DArowXh1[k]);  // D_A_k Y 
	   } 
             vec_add_mult(DArowXh2,xih,dhazA,DArowXh2); // D_A T
	     vec_add_mult(DArowXh,xih,dsurvA,DArowXh);  // D_A N
           if (*fixhaplofreq==0 &&  *twostage==0  ) {
             for (k=0;k<*dimhap;k++) {
                vec_add_mult(DnurowXh1[k],xih,VE(XtempScore,k)*shaz,DnurowXh1[k]);  
                vec_add_mult(DnurowXh2[k],xih,VE(XtempScore,k)*surv,DnurowXh2[k]); } 
	   } 

	 // }}}
	 ih+=1; c1+=2;  
         }  // for (j=0;j<nphpp[i];j++)  

       // set up  X and Z matrices {{{
      weight[i]=naevn/tael;  
      if (i==idtimes[s-1]) {naevnpers=naevn; taelpers=tael;}
        // if (naevn<=0) naevn=1; if (tael<=0) tael=1; 
        scl_vec_mult(1/naevn,rowXh,rowXh); replace_row(cdesX,i,rowXh); 
        scl_vec_mult(naevn/tael,rowXh,xih); replace_row(cdesX2,i,xih); 


      if (*fixbeta==0 ) { 
         scl_vec_mult(1/tael,rowZbeta1,zih1); scl_vec_mult(1/naevn,rowZbeta2,zih2); 
         vec_subtr(zih1,zih2,rowZh); 
         for (j=0;j<*dimzih;j++) ME(Z,i,j)=VE(rowZh,j); 
	 extract_row(cdesX,i,xih1); 
         for (k=0;k<*dimzih;k++) {
            for (c=0;c<*dimzih;c++)   {
	      VE(DbetaZ1[k],c)=(tael*VE(DbetaZ1[k],c)-VE(rowZbeta1,c)*VE(rowZbeta1,k))/(tael*tael); 
	      VE(DbetaZ2[k],c)=(naevn*VE(DbetaZ2[k],c)-VE(rowZbeta2,c)*VE(rowZbeta2,k))/(naevn*naevn); 
	    }
            vec_subtr(DbetaZ1[k],DbetaZ2[k],rowZh); // NY version med W 
            for (c=0;c<*dimzih;c++) for (l=0;l<*dimxih;l++)  
	       ME(DthetaZY[k],c,l)+=VE(rowZh,c)*VE(xih1,l);
	    if (i==idtimes[s-1]) for (j=0;j<*dimzih;j++) VE(DthetaZdN[k],j)=VE(rowZh,j); 
            VE(rowZtheta1,k)=VE(rowZbeta1,k); VE(rowZtheta2,k)=VE(rowZbeta2,k); 
	    if (i==idtimes[s-1]) {VE(rowZtheta1dN,k)=VE(rowZbeta1,k);VE(rowZtheta2dN,k)=VE(rowZbeta2,k);}
         }
	 }

      if (*fixhaplofreq==0 && *twostage==0  ) {
         scl_vec_mult(1/tael,rowZnu1,vhaplo1); 
	 scl_vec_mult(1/naevn,rowZnu2,vhaplo2); 
         vec_subtr(vhaplo1,vhaplo2,vhaplo1); 
         for (j=0;j<*dimhap;j++) {ME(Z,i,*dimzih*(*fixbeta==0)+j)=VE(vhaplo1,j); 
                                 VE(rowZtheta1,*dimzih*(*fixbeta==0)+j)=VE(rowZnu1,j); 
		                 VE(rowZtheta2,*dimzih*(*fixbeta==0)+j)=VE(rowZnu2,j); 
	    if (i==idtimes[s-1]) {VE(rowZtheta1dN,*dimzih*(*fixbeta==0)+j)=VE(rowZnu1,j); 
		          VE(rowZtheta2dN,*dimzih*(*fixbeta==0)+j)=VE(rowZnu2,j); }
	 }
      }

      // set up matrices with joint derivatives 
       if (*fixbeta==0) {
          for (k=0;k<*dimzih;k++) {
             scl_vec_mult(1/naevn,DbetarowXh1[k],DbetarowXh1[k]);
             scl_vec_mult(VE(rowZbeta2,k)/(naevn),rowXh,DbetarowXh2[k]); 
             vec_subtr(DbetarowXh1[k],DbetarowXh2[k],DbetarowXh1[k]); 
             for (c=0;c<*dimxih;c++)  ME(DthetaX[k],i,c)=VE(DbetarowXh1[k],c); 
	  }
       }

       if (*fixhaplofreq==0 && *twostage==0 )
       for (k=0;k<*dimhap;k++) {
            scl_vec_mult(1/naevn,DnurowXh1[k],DnurowXh1[k]);
            scl_vec_mult(VE(rowZnu2,k)/(naevn),rowXh,DnurowXh2[k]);   // check this 
            vec_subtr(DnurowXh1[k],DnurowXh2[k],DnurowXh1[k]); 
            for (c=0;c<*dimxih;c++) ME(DthetaX[*dimzih*(*fixbeta==0)+k],i,c)=VE(DnurowXh1[k],c); 
	 }
//      Rprintf(" mig io 4 \n");

       for (k=0;k<*dimxih;k++) {
          scl_vec_mult(1/naevn,DArowXh1[k],xih1);
          scl_vec_mult(VE(DArowXh,k)/naevn,rowXh,xih2);   // check this 
          vec_subtr(xih1,xih2,xih2); 
          replace_row(DAX[k],i,xih2);  
       }

	  extract_row(cdesX,i,xih1); 
          for (k=0;k<*dimxih;k++)  {
	     extract_row(DAX[k],i,xih2); 
             for (c=0;c<*dimxih;c++)  for (l=0;l<*dimxih;l++)  {
	         ME(DAS0[k],l,c)+= (VE(xih1,l)*VE(xih2,c)+VE(xih2,l)*VE(xih1,c))*(naevn/tael) +
	                  VE(xih1,l)*VE(xih1,c)*
		          (VE(DArowXh,k)*tael-VE(DArowXh2,k)*naevn)/(tael*tael); //check NY VERSION  med W
	     } 
	  }


          if (*fixhaplofreq==0 || *fixbeta==0) 
          for (k=0;k<dimpar;k++)  {
	     extract_row(DthetaX[k],i,xih2); 
             for (c=0;c<*dimxih;c++)  for (l=0;l<*dimxih;l++)  {
	        ME(DthetaS0[k],l,c)+= (VE(xih1,l)*VE(xih2,c)+VE(xih2,l)*VE(xih1,c))*(naevn/tael)+
	         VE(xih1,l)*VE(xih1,c)*
		 (VE(rowZtheta2,k)*tael-VE(rowZtheta1,k)*naevn)/(tael*tael); // NY VERSION  med W
	     }
           }
	  // Rprintf("00 lhhhhhhhhhhhhhhhh\n"); 
       // }}}

      } /* i = 0 ... *anpers */

      
      // baseline computed  {{{
      MtA(cdesX,cdesX2,A); invert(A,AI); scl_mat_mult(1.0,AI,S0tI[s]);          // NY version med W 
      extract_row(cdesX2,idtimes[s-1],xih);  Mv(AI,xih,dA); 
      scl_vec_mult(1.0,dA,dAt[s]); 
     for (k=1;k<=*dimxih;k++) {cu[k*(*Ntimes)+s]=cu[k*(*Ntimes)+s-1]+VE(dA,k-1);}
      // }}}

     if (s<0) { // {{{ test 
        Rprintf(" sss er %d %d \n",s,idtimes[s-1]); 
        Rprintf("X matrix er  \n"); head_matrix(cdesX); 
        Rprintf("X matrix er  \n"); head_matrix(cdesX2); 
        Rprintf("Z matrix er  \n"); head_matrix(Z);
         extract_row(Z,idtimes[s-1],vtheta1); 
	 print_vec(vtheta1); 
        print_mat(AI); 
        print_mat(A); 
     } // }}}

      // First derivative U   {{{
      if (( *fixbeta==0 || *fixhaplofreq==0 ) ) {
	 MtA(Z,cdesX,ZX); MxA(ZX,AI,ZXAIs[s]);
         extract_row(Z,idtimes[s-1],vtheta1); 
         Mv(ZX,dA,ZXdA[s]); vec_subtr(vtheta1,ZXdA[s],difzzav); vec_add(difzzav,U,U); 
      }
      // }}}
      
      // derivatives with respect to A_k for taylor series {{{
      for (k=0;k<*dimxih;k++) { 
         MtA(cdesX2,DAX[k],E[k]);  
	 if (*fixbeta==0 || *fixhaplofreq==0) {
	    MtA(Z,DAX[k],ZDAX[k]);  MxA(ZXAIs[s],E[k],tmp3); 
            mat_subtr(ZDAX[k],tmp3,Smat[k]); scl_mat_mult(-1,Smat[k],Smat[k]); 
	 }
         MxA(AI,E[k],E[k]);  
      }
      // }}}

     // Second derivative S   {{{
      for (k=0;k<dimpar;k++) {
          MtA(Z,DthetaX[k],tmp1); Mv(tmp1,dA,vtheta2); 
          MxA(ZXAIs[s],DthetaS0[k],tmp2); Mv(tmp2,dA,vtheta1);  
          extract_row(DthetaX[k],idtimes[s-1],rowXh); 
          scl_vec_mult(weight[idtimes[s-1]],rowXh,xih1); // NY version med W 
          dummy=(VE(rowZtheta2dN,k)*taelpers-VE(rowZtheta1dN,k)*naevnpers)/(taelpers*taelpers); // NY 
          extract_row(cdesX,idtimes[s-1],rowXh);  
	  scl_vec_mult(dummy,rowXh,rowXh);  // NY verssion med W
	  vec_add(xih1,rowXh,rowXh); 
	  Mv(ZXAIs[s],rowXh,vtheta3);  // NY verssion med W
          vec_add(vtheta1,vtheta3,vtheta3); 
	  // Rprintf(" D2 %d %d \n",s,k); print_vec(vtheta3); 
          vec_subtr(vtheta2,vtheta3,vtheta2); 
	  if (*fixbeta==-1 && k<*dimzih) { // kun anden afledet mht beta 
             Mv(DthetaZY[k],dA,zih1); vec_subtr(DthetaZdN[k],zih1,zih2); 
	     // Rprintf(" s k %d %d \n",s,k); print_vec(zih1); print_vec(zih2); 
	     for (j=0;j<*dimzih;j++)  VE(vtheta2,j)+=VE(zih2,j); 
	  }
	  scl_vec_mult(1,vtheta2,vtheta2); 
          replace_row(dS,k,vtheta2); 
       }
       if (*sym==1) {mat_transp(dS,S2); mat_add(S2,dS,dS); scl_mat_mult(0.5,dS,dS);} 
       mat_add(dS,S1,S1);  scl_mat_mult(1.0,S1,St[s]); 
    // }}}

   // variance and other things  {{{
    if (it==((*Nit)-1)) { 
    if (*fixhaplofreq==0 && s==1)  {
	  for (j=0;j<*nph-1;j++) { scoregenoC[j]=0; 
	    for (l=0;l<*nph-1;l++) d2lgenoC[j*(*nph-1)+l]=0; } 
      amount[0]=2; 
      genoLogLikeHp(oh,nphpp,nph,amount,antpers,haplopars,rho,LogLikegeno,scoregenoC,scorei,d2lgenoC); 

     for (j=0;j<*antpers;j++) vM(hapDes,scorei[j],Xscorei[j]);  
    }
        // if (*fixbeta==0 ) { for (j=0;j<*pg;j++)  ME(Utt,s,j+1)=VE(U,j);  }
        replace_row(Utt,s,U);
        // computation of dF and dG increments 
	for (j=0;j<*dimxih;j++)  {
	  Mv(E[j],dA,xih); replace_col(dFt[s],j,xih); 
	  if (*fixbeta==0 || *fixhaplofreq==0) {
	  Mv(Smat[j],dA,vtheta1); 
	  replace_col(dGt[s],j,vtheta1); 
	  }
	}
        // computation of YIt  
	if (*fixbeta==0 || *fixhaplofreq==0) {
	   for (j=0;j<dimpar;j++)  {
	   //   Rprintf(" YIt %d \n",j); print_mat(S0tI[s]); print_mat(DthetaS0[j]); 
              MxA(S0tI[s],DthetaS0[j],tmp4); Mv(tmp4,dA,xih2);  
              extract_row(DthetaX[j],idtimes[s-1],xih1); 
	      scl_vec_mult(weight[idtimes[s-1]],xih1,xih1); // NY version med W 
              extract_row(cdesX,idtimes[s-1],rowXh);  // NY version med W
	      dummy=(VE(rowZtheta2dN,j)*taelpers-VE(rowZtheta1dN,j)*naevnpers)/(taelpers*taelpers); 
	      scl_vec_mult(dummy,rowXh,rowXh);  // NY verssion med W
	      vec_add(rowXh,xih1,xih1); 
	      Mv(S0tI[s],xih1,xih1); 
	      vec_subtr(xih1,xih2,xih1); 
	      replace_col(dYIt[s],j,xih1); extract_col(YIt[s-1],j,xih2); 
	      vec_add(xih2,xih1,xih1); replace_col(YIt[s],j,xih1); 
	  //     Rprintf(" %d %d \n",s,j); Mv(ZX,xih1,vtheta1); print_vec(vtheta1); 
	   }
	   // Rprintf(" %d \n",s); print_mat(YIt[s]); 
        } 
        // computation of q1(t) and dLam_i(t) (cumulative increments of compensator) 
	for (j=0;j<*antpers;j++){
	   extract_row(cdesX2,j,xih);
	   Mv(S0tI[s],xih,rowXh); replace_row(AIxit[j],s,rowXh);
           if (*fixbeta==0 || *fixhaplofreq==0) {
	      extract_row(Z,j,vtheta1);  Mv(ZX,rowXh,vtheta2); vec_subtr(vtheta1,vtheta2,vtheta1); 
	      replace_row(q1t[j],s,vtheta1); 
	   }
	   extract_row(cdesX,j,xih);
	   VE(dLamt[j],s)=vec_sum(vec_star(xih,dA,rowXh)); // define plamt exp(zih 

        }
      }
    // }}}

    } // Ntimes  

    // Newton-Raphson step  {{{
     if (*fixbeta==0 || *fixhaplofreq==0) {
     if (*lm>0) { // Levenberg-Marquards algorithm 
        if (*detail>=1)  Rprintf(" Levenberg-Marquardt steps, sumscore, maxdelt %d %lf %lf \n",
			it,sumscore,maxdelt);
         mat_transp(S1,S2); MxA(S1,S1,S2);
         for (k=0;k<dimpar;k++) ME(S2,k,k)=ME(S2,k,k)+lm[0];
         invert(S2,SI); MxA(SI,S1,S2); Mv(S2,U,delta);
         if (maxdelt<step[0]) {lm[0]=minlm[0]; step[0]=1;}
     }
     else { // Newton-Raphson step 
      if (*detail>=1)  Rprintf(" Newton-Raphson steps, sumscore %d %lf %lf \n",it,sumscore,maxdelt);
       invert(S1,SI); Mv(SI,U,delta);
     }

   maxdelt=0; 
   for (k=0;k<dimpar;k++) { maxdelt+=fabs(VE(delta,k));}

   if (*lm==0)  scl_vec_mult(step[0],delta,delta);

    if (*detail>=2) { 
        Rprintf("====================Iteration %d ==================== \n",it);
        Rprintf("Estimate beta \n"); print_vec(pars); 
        Rprintf("Score D l\n"); print_vec(U); 
        Rprintf("Information D^2 l\n"); print_mat(SI); 
        Rprintf("simple D2 l\n");  print_mat(S1); 
        Rprintf("delta \n"); print_vec(delta); 
    }
    vec_add(pars,delta,pars); 
    sumscore=0; 
    for (k=0;k<dimpar;k++) {sumscore +=fabs(VE(U,k));}
     if ((fabs(sumscore)<0.000000001) & (it<*Nit-2)){ it=*Nit-2; }
    }
  // }}}

  } // it 

     // computation of q2(t) {{{
     if (*fixbeta==0 || (*fixhaplofreq==0 && twostage==0)) {
     for (s=1;s<*Ntimes;s++) {
       mat_zeros(M1M2t); 
       for (t=s;t<*Ntimes;t++) { 
         identity_matrix(tmp4); identity_matrix(M1); 
         for (k=s;k<t;k++) {
	   if (k>s) {scl_mat_mult(1,M1,tmp4); }
	   mat_subtr(Ident,dFt[k],A); MxA(tmp4,A,M1); 
         }
         MxA(dGt[t],M1,dM1M2);
         mat_add(dM1M2,M1M2t,M1M2t); 
       }
       scl_mat_mult(1,M1M2t,q2t[s]); 
     }
  } 
  // }}}

  lle=0; llo=0;
// terms for robust variances  {{{
if (robust==1) {
    for (s=1;s<*Ntimes;s++) {
      time=times[s]; dtime=time-times[s-1]; 
      cu[s]=times[s]; vcu[s]=times[s]; Rvcu[s]=times[s]; Ut[s]=times[s]; 
      // terms for robust variance   
      vec_zeros(delta); 

      for (i=0;i<*antpers;i++) {
        extract_row(AIxit[i],s,xih); 

        if (*fixbeta==0 || (*fixhaplofreq==0 && twostage==0)) {
           Mv(q2t[s],xih,vtheta1);  
	   extract_row(q1t[i],s,vtheta2); 
//	    Rprintf(" %ld %ld %lf \n",s,i,VE(dLamt[i],s)); 
//	    print_vec(vtheta1); print_vec(vtheta2); 
	   vec_add(vtheta2,vtheta1,vtheta1); 

	   if (i==idtimes[s-1])  { 
	     for (j=0;j<dimpar;j++)  for (k=0;k<dimpar;k++) ME(VU,j,k) += VE(vtheta1,j)*VE(vtheta1,k); 
	   }

	   if (i==idtimes[s-1]) {vec_add(vtheta1,W2[i],W2[i]);}
	   scl_vec_mult(VE(dLamt[i],s),vtheta1,vtheta3);
	   vec_subtr(W2[i],vtheta3,W2[i]);
	   // if (*ratesim==1) {sv_mlt(hati,tmpv2,rowZ); v_sub(W2[i],rowZ,W2[i]);}  
	   replace_row(W2t[i],s,W2[i]); 
	}

	if  (*fixhaplofreq==0) {
           if (s==1)  { 
	      for (j=0;j<*dimhap;j++) 
	      for (k=0;k<*dimhap;k++) 
	         ME(VU,(*dimzih)*(*fixbeta==0)+j,(*dimzih)*(*fixbeta==0)+k)+=VE(Xscorei[i],j)*VE(Xscorei[i],k); 
	   } 
	}

	vec_zeros(W3[i]); 
	for (t=1;t<=s;t++) { 
	   if (i==0) {
              identity_matrix(tmp4); identity_matrix(M1); 
	      for (k=t;k<=s;k++) {
		if (k>t) { scl_mat_mult(1.0,M1,tmp4); }
		if (k>t || t==s) {
	            mat_subtr(Ident,dFt[k],A); 
		    MxA(tmp4,A,M1); } 
	       }
	       scl_mat_mult(1,M1,Fst[s*(*Ntimes)+t]); 
           }
	   // print_mat(Fst[s*(*Ntimes)+t]); 
	   //  Fst[s*(*Ntimes)+t]->me[0][0]=exp(-lht->ve[t]+lht->ve[s]); 
	   extract_row(AIxit[i],t,xih); 
	   vM(Fst[s*(*Ntimes)+t],xih,rowXh); 
	   scl_vec_mult(VE(dLamt[i],t),rowXh,tmpv1);
	   vec_subtr(W3[i],tmpv1,W3[i]); 
	   if (i==idtimes[t-1]){vec_add(rowXh,W3[i],W3[i]);}
	}

	replace_row(W3t[i],s,W3[i]);  
//	Rprintf("W3  %ld %ld \n",s,i); 
////	print_vec(W3[i]); 
//	extract_row(AIxit[i],s,xih); 
//	print_vec(xih); 
	// if (hati>0) lle=lle+log(hati); llo=llo+hati; 
	// if (*ratesim==1) {sv_mlt(hati,rowX,rowX); v_sub(W3[i],rowX,W3[i]);} 
	if (*retur==1){ dhatMit[i*(*Ntimes)+s]=1*(i==idtimes[s-1])-VE(dLamt[i],s); }
        if (*fixhaplofreq==0) {
           for (j=0;j<*dimhap;j++) {
	      if (s==*Ntimes-1) VE(W2[i],(*dimzih)*(*fixbeta==0)+j)+=VE(Xscorei[i],j);
	      ME(W2t[i],s,(*dimzih)*(*fixbeta==0)+j)+=VE(Xscorei[i],j);
	   } 
        }
      } // i=1..antpers 
    } // s=1 ..Ntimes 
    if (*fixbeta==0 || *fixhaplofreq==0) { MxA(SI,VU,S2); MxA(S2,SI,VU); }
}
// }}}

// print test {{{ 
/*
    for (i=0;i<*antpers;i++) {
      Rprintf(" %ld \n",i); 
      print_vec(W2[i]); 
    }
    vec_zeros(delta); 	
    for (i=0;i<*antpers;i++) {
     // Rprintf(" %ld \n",i); 
     // print_vec(Xscorei[i]); 
     vec_add(Xscorei[i],delta,delta); }
Rprintf(" individual scores added \n"); 
print_vec(delta); 
print_vec(U); 
*/
// }}} 

    // ROBUST VARIANCES computation {{{ 
    vec_zeros(rowXh); 
    for (s=1;s<*Ntimes;s++) {
      vec_zeros(VdB);
      for (i=0;i<*antpers;i++) {
        if (*fixbeta==0 || *fixhaplofreq==0) {
          Mv(SI,W2[i],vtheta1);   
          if (s==1) {   
	     for (j=0;j<dimpar;j++){ 
	     for (k=0;k<dimpar;k++) ME(RobVbeta,j,k) += VE(vtheta1,j)*VE(vtheta1,k); } 
	  }
	  Mv(YIt[s],vtheta1,rowXh);
        }
	else if (s==1) {vec_zeros(rowXh); vec_zeros(vtheta1);}
	extract_row(W3t[i],s,xih); 
	vec_add(xih,rowXh,difX); 
	replace_row(W4t[i],s,difX);
	vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

        if (*resample==1) {
        if (s==1) 
	      if (*fixbeta==0) {
		for (c=0;c<*dimzih;c++) gamiid[c*(*antpers)+i]=VE(vtheta1,c); }
	    for (c=0;c<*dimxih;c++) {l=i*(*dimxih)+c; 
		                     biid[l*(*Ntimes)+s]=VE(difX,c);} 
	} 

        if (*fixbeta==0 || *fixhaplofreq==0) {
           // observed score process asymptotics 
	   Mv(St[s],vtheta1,vtheta2); 
	   extract_row(W2t[i],s,vtheta3); 
	   vec_subtr(vtheta3,vtheta2,vtheta3); 
	   if (*fixbeta==0) for(j=0;j<*dimzih;j++) { ME(Uti[i],s,j)=VE(vtheta3,j); }
	   // replace_row(Uti[i],s,vtheta3); 
	   if (i==-10)   {
	    print_mat(St[s]); print_vec(vtheta2); print_vec(vtheta3); print_vec(W2[i]); 
	   }
	   vec_star(vtheta3,vtheta3,vtheta2); 
	   vec_add(vtheta2,varUthat[s],varUthat[s]);
	}

/* 	if (*retur==1) { // {{{ does not work
	  mat_zeros(ldesignX); mat_zeros(ldesignG); 
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	    if ((start[c]<time) && (stop[c]>=time))  {
	      risk[id[c]]=1;
	      for(j=0;j<pmax;j++) {
		if (j<*px) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		if (j<*pg) ME(ldesignG,id[c],j)=designG[j*(*ng)+c]; 
	      } 
	      if (time==stop[c] && status[c]==1) { pers=id[c]; } 
	      count=count+1; } }
	  Mv(ldesignG,beta,Gbeta); 
	  for (j=0;j<*antpers;j++) {
	    extract_row(ldesignX,j,xi); 
	    extract_row(ldesignG,j,zi);
	    for (i=0;i<*dimxih;i++){
	      VE(tmpv1,i)=cu[(i+1)*(*Ntimes)+s-1];
	    }
	    vec_star(tmpv1,xi,rowXh); 
	    hati=vec_sum(rowXh);
	    VE(plamt,j)=exp(VE(Gbeta,j)+hati)/(1+exp(VE(Gbeta,j)+hati)); 
	    scl_vec_mult(VE(plamt,j),xi,xtilde); 
	    replace_row(cdesX,j,xtilde);       
	  }
	  Mv(cdesX,dAt[s],lamt);  
	  for (j=0;j<*antpers;j++){
	    extract_row(ldesignG,j,zi); 
	    scl_vec_mult(VE(lamt,j),zi,zi); 
	    replace_row(ZP,j,zi);
	  } 
	  Mv(ZP,W2[i],reszpbeta);
	  Mv(dYIt[s],W2[i],xi); 
	  Mv(cdesX,xi,res1dim); 
	  dhatMitiid[i*(*Ntimes)+s]=dhatMit[i*(*Ntimes)+s]- (VE(reszpbeta,0)- VE(res1dim,0)); 
	} // retur ==1  // }}} */ 

      } // i =1 ..Antpers 

      for (k=1;k<*dimxih+1;k++) { 
	Rvcu[k*(*Ntimes)+s]=VE(VdB,k-1); vcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      }
    }  //  s=1 ..Ntimes  
// }}}

// returning some arguments to R  {{{
  if (*fixbeta==0 || *fixhaplofreq==0) {
     for (k=0;k<dimpar;k++) score[k]=VE(U,k);
     if (*fixbeta==0) {for (j=0;j<*dimzih;j++) VE(beta,j)=VE(pars,j);}
     if (*fixhaplofreq==0) {
            for (j=0;j<*dimhap;j++) VE(Uhap,j)=VE(pars,(*fixbeta==0)*(*dimzih)+j);
	    Mv(hapDes,Uhap,Vhaplopars); 
	     for (j=0;j<*dimhap;j++)  alpha[j]=VE(Uhap,j);
	     // for (j=0;j<*nph-1;j++)  haplopars[j]=VE(Vhaplopars,j);
     }
   for(j=0;j<dimpar;j++) { 
     params[j]= VE(pars,j); loglike[0]=ll;
     for (k=0;k<dimpar;k++){ 
       Iinv[k*(dimpar)+j]=ME(SI,j,k); Vbeta[k*(dimpar)+j]=ME(VU,j,k); 
       RVbeta[k*(dimpar)+j]=ME(RobVbeta,j,k); 
     } 
   } 
  }
// }}}

  if (*sim==1) { // {{{
    // Rprintf("Simulations start N= %d \n",*antsim);
    GetRNGstate();  // to use R random normals 
    tau=times[*Ntimes-1]-times[0];
    for (i=1;i<=*dimxih;i++) VE(rowXh,i-1)=cu[i*(*Ntimes)+(*Ntimes-1)];

    // Beregning af OBS teststørrelser // {{{
    for (s=1;s<*Ntimes;s++) {  
      time=times[s]-times[0]; 
      for (i=1;i<=*dimxih;i++) {
	VE(xih,i-1)=fabs(cu[i*(*Ntimes)+s])/sqrt(Rvcu[i*(*Ntimes)+s]);
	if (VE(xih,i-1)>testOBS[i-1]) testOBS[i-1]=VE(xih,i-1); 
      }
      scl_vec_mult(time/tau,rowXh,difX);
      for (i=1;i<=*dimxih;i++) {
	VE(xih,i-1)=cu[i*(*Ntimes)+s];
      }
      vec_subtr(xih,difX,difX);
      for (i=0;i<*dimxih;i++) {
	VE(difX,i)=fabs(VE(difX,i)); 
	l=(*dimxih+i);
	if (VE(difX,i)>testOBS[l]) testOBS[l]=VE(difX,i);
      }
      if (*fixbeta==0 || *fixhaplofreq==0) {
         if (*wscore>=1) {  // sup beregnes i R 
	   if ((s>*wscore) && (s<*Ntimes-*wscore)) {
	     extract_row(Utt,s,vtheta1);
	     for (i=0;i<dimpar;i++) {
	       VE(vtheta1,i)=VE(vtheta1,i)/sqrt(VE(varUthat[s],i));
	     }
	     replace_row(Utt,s,vtheta1); // scaled score process 
	   } else {
	     vec_zeros(vtheta1); 
	     replace_row(Utt,s,vtheta1);
	   }
         }
         for (k=1;k<=dimpar;k++){ Ut[k*(*Ntimes)+s]=ME(Utt,s,k-1); }
      }
    } // s=1..Ntimes Beregning af obs teststørrelser  }}}

    for (k=1;k<*antsim;k++) { // {{{
      mat_zeros(Delta); mat_zeros(Delta2); vec_zeros(tmpv1);
      for (i=0;i<*antpers;i++) {
	// random=gasdev(&idum);  
	random=norm_rand();
	scl_mat_mult(random,W4t[i],tmpM1); 
	mat_add(tmpM1,Delta,Delta);
        if (*fixbeta==0 || *fixhaplofreq==0) {
           scl_mat_mult(random,Uti[i],tmpM2); 
	   mat_add(tmpM2,Delta2,Delta2);
	}
      }
      extract_row(Delta,*Ntimes-1,tmpv1); 
       for (s=1;s<*Ntimes;s++) { 
	time=times[s]-times[0];
	scl_vec_mult(time/tau,tmpv1,xih); 
	extract_row(Delta,s,rowXh);
	vec_subtr(rowXh,xih,difX);
	for (i=0;i<*dimxih;i++) {
	  VE(difX,i)=fabs(VE(difX,i));
	  l=(*dimxih+i);
	  if (VE(difX,i)>test[l*(*antsim)+k]) test[l*(*antsim)+k]=VE(difX,i);
	  VE(xih,i)=fabs(ME(Delta,s,i))/sqrt(Rvcu[(i+1)*(*Ntimes)+s]);
	  if (VE(xih,i)>test[i*(*antsim)+k]) test[i*(*antsim)+k]=VE(xih,i); 
	}
        if (*fixbeta==0) { // only for zih part of score
	  extract_row(Delta2,s,zih); 
	  if (*wscore>=1) {
	  if ((s>*wscore) && (s<*Ntimes-*wscore))  {
	  for (i=0;i<*dimzih;i++) {
	    VE(zih,i)=fabs(ME(Delta2,s,i))/sqrt(VE(varUthat[s],i));
	    if (VE(zih,i)>simUt[i*(*antsim)+k]) simUt[i*(*antsim)+k]=VE(zih,i); 
	 }
	 if (k<50) { 
	         for (i=0;i<*dimzih;i++) { l=(k-1)*(*dimzih)+i;
		   Uit[l*(*Ntimes)+s]=ME(Delta2,s,i)/sqrt(VE(varUthat[s],i)); } }
	     } 
	   } else { // weigted score 
	     for (i=0;i<*dimzih;i++) {
	       if (fabs(VE(zih,i))>simUt[i*(*antsim)+k]) simUt[i*(*antsim)+k]=
		                                            fabs(VE(zih,i)); 
	     }
	     if (k<50) { 
	       for (i=0;i<*dimzih;i++) {l=(k-1)*(*dimzih)+i;  
		        Uit[l*(*Ntimes)+s]=ME(Delta2,s,i); } }
	   } // else wscore=0 
	}
      }  // s=1..Ntims
    }  // k=1..antsim  }}}

    PutRNGstate();  // to use R random normals 
  } // sim==1  }}}

 // freeing  all matrices and vectors  {{{
  for (j=0;j<*dimxih;j++) { free_mat(DAX[j]); free_mat(E[j]); 
     free_mat(DAS0[j]); free_mat(ZDAX[j]); free_mat(Smat[j]); 
  }
  for (j=0;j<dimpar;j++) {
     free_mat(DthetaS0[j]); free_mat(DthetaZ[j]);  free_mat(DthetaZY[j]); free_vec(DthetaZdN[j]);      
  }
  for (j=0;j<*antpers;j++) { 
    free_vec(dLamt[j]); free_mat(W3t[j]); free_mat(dW3t[j]); free_mat(W4t[j]);
    free_mat(W2t[j]); free_mat(Uti[j]); free_vec(W2[j]); free_vec(W3[j]);
    free_mat(q1t[j]); free_mat(AIxit[j]);
  }
  free_mat(Delta); free_mat(tmpM1);
  free_mat(Delta2); free_mat(tmpM2);
  free_mat(Utt); free_vec(lht); free_vec(reszpbeta); free_vec(res1dim);
  free_mats(&ldesignX,&cdesX2,&cdesX,NULL);
  free_mats(&ldesignG,NULL); 
  free_mats(&ZP,&Z,NULL); 
  free_mats(&tmp4,&Ident,&COV,&A,&AA,&AI,&M1,&CtVUCt,NULL); 
  free_mats(&RobVbeta,&dS,&S1f,&S1,&S2,&SI,&VU,&VUI,NULL); 
  free_mats(&ZXAI,&tmp5,&tmp3,&ZX,&dM1M2,&M1M2t,NULL); 
  free_mats(&dYI,&Ct,NULL); 
  free_mat(tmp6);
  free_mat(hapDes);
  free_mats(&tmp1,&tmp2,NULL);
  free_vecs(&Lplamt,&dlamt,&plamt,&lamt,&zcol,&Gbeta,&one,&offset,NULL); 
  free_vecs(&DArowXh,&ahatt,&tmpv1,&difX,&VdB,&rowXh,&xih1,&xih2,&xih,&dA,&VdA,&MdA,NULL); 
  free_vecs(&tmpv2, &rowZbeta1,&rowZbeta2, &rowZbeta1dN,&rowZbeta2dN,
		  &rowZh,&zih1,&zih2,&zih,&beta,NULL); 
  free_vecs(&pars,&Uhap,&Vhaplopars,&Uf,&U,&vtheta1,&vtheta2,&vtheta3,&delta,&zav,&difzzav,&Uprofile,NULL); 
  free_vecs(&rowZtheta1dN,&rowZtheta2dN, &rowZtheta1,&rowZtheta2,NULL); 
  free_vec(xi); free_vec(zi); 
  for(j=0;j<*Ntimes;j++) {
    free_mat(dFt[j]); free_mat(gt[j]); free_mat(dGt[j]);
    free_mat(q2t[j]); free_mat(S0tI[j]); free_mat(dG[j]);
    free_mat(C[j]); free_mat(M1M2[j]); free_mat(ZXAIs[j]); 
    free_mat(YIt[j]); free_mat(dYIt[j]); free_vec(dAt[j]); free_vec(ZXdA[j]);
    free_mat(St[j]); free_vec(varUthat[j]);
    for(i=0;i<=j;i++){ free_mat(Fst[j*(*Ntimes)+i]); }
  }
  for(j=0;j<*antpers;j++)  { free_vec(scorei[j]); free_vec(Xscorei[j]); }
  for (j=0;j<*dimzih;j++) { free_mat(DbetaX[j]); }
  for (j=0;j<dimpar;j++) { free_mat(DthetaX[j]); }
  for (j=0;j<*dimhap;j++) { free_mat(DnuX[j]); }
  free_vecs(&XtempScore,&tempScore,&rowZnu1,&rowZnu2,&rowZnu,&vhaplo1,&vhaplo2,NULL); 
  for(k=0;k<*dimzih;k++) { free_vec(DbetarowXh1[k]); free_vec(DbetarowXh2[k]); }
  for(k=0;k<dimpar;k++) { free_vec(DbetaZ1[k]); free_vec(DbetaZ2[k]); }
  for(k=0;k<*dimxih;k++) { free_vec(DArowXh1[k]); }
  free_vec(DArowXh2); 
  for(k=0;k<*dimhap;k++) { free_vec(DnurowXh1[k]); free_vec(DnurowXh2[k]); }
  free_vec(DnurowXh); free_vec(DbetarowXh); 

   free(scoregenoC); free(d2lgenoC); free(weight);  free(ipers); free(risk); 
  // }}}

}

/*
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

void haplodesign(vector *xi,int *haplotype,vector *xih,SEXP f,SEXP rho) 
 // {{{
{
   int i,j,dimx,dimxih; 
   dimx=length_vector(xi); 
   dimxih=length_vector(xih); 
   double x[dimx],h[2]; 
   double res[dimxih]; 
   for (j=0;j<2;j++) h[j]=haplotype[j]*1.0;
   for (j=0;j<dimx;j++) x[j]=VE(xi,j);
   for (j=0;j<dimxih;j++) res[j]=0; 
   designeval(x, h, f , res,rho,dimx, dimxih);
   for(i=0;i<dimxih;i++) VE(xih,i)=res[i]; 
} 
// }}}
*/
