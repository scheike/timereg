//#include <stdio.h>
#include <math.h>
#include "matrix.h"
                 
void dynadd(times,y,Ntimes,designX,nx,px,designA,na,pa,ahat,bhat,bhatny,nxval,antpers,
start,stop,cu0,cuf,cuMS,vcu0,vcuf,robvcu,w,mw,rani,sim,antsim,cumBit,test,
testOBS,status,Ut,simUt,b,cumly,retur,id,smoothXX,weighted,vculy,clusters,antclust)
double *bhatny,*bhat,*ahat,*designX,*designA,*times,*y,*start,*stop,*cu0,*cuf,*cuMS,
*vcu0,*vcuf,*w,*robvcu,*cumBit,*test,*testOBS,*Ut,*simUt,*b,*cumly,*vculy;
int *sim,*antsim,*retur,*nxval,*nx,*px,*na,*pa,*antpers,*Ntimes,*mw,*rani,*status,*id,*smoothXX,*weighted,*clusters,*antclust;
{
  matrix *ldesignX,*ldesignA,*cdesignX,*cdesignA,*Aa,*AaI,*A,*AI; 
  matrix *XbXa,*XWX;   
  vector *korf,*dB,*dA,*dR,*ahatt,*xt,*pdA,*diag,*xai,*sumx,*vone,*itot;
  vector *VdBly,*VdB,*VdB0,*VdBf,*fkor,*dkorB,*tmpv,*tmpv1,*tmpv2,*tmpv3,*tmpv4; 
  vector *pahat,*pbhat,*bhatt,*pbahat,*pdbahat,*dBly;
  vector *dAt[*Ntimes];
  matrix *cumBt[*antpers];
  vector *cumhatB[*antpers],*cumB[*antpers],*cum;
  int pers=0,i,j,k,s,c,count,pmax,nmax,risk;
  int *coef=calloc(1,sizeof(int)),*imin=calloc(1,sizeof(int)),
      *ps=calloc(1,sizeof(int)),*degree=calloc(1,sizeof(int));
  double time,zpers=0,dif,dtime,YoneN,kia;
  double *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double));

  for (i=0;i<*antpers;i++) { malloc_vec(*px,cumhatB[i]); malloc_vec(*px,cumB[i]);
    malloc_mat(*Ntimes,*px,cumBt[i]); }
  for (i=0;i<*Ntimes;i++) malloc_vec(*pa,dAt[i]);

  malloc_vec(*px,cum); malloc_mat(*px,*pa,XbXa); 
  malloc_mat(*antpers,*px,ldesignX); malloc_mat(*antpers,*px,cdesignX); 
  malloc_mat(*antpers,*pa,ldesignA);  malloc_mat(*antpers,*pa,cdesignA);  
  malloc_mat(*pa,*pa,Aa); malloc_mat(*pa,*pa,AaI);
  malloc_mats(*px,*px,&XWX,&A,&AI,NULL);

  malloc_vecs(*px,&dB,&diag,&sumx,&itot,&korf,&dBly,&bhatt,&fkor,&VdBly,&VdB,&VdB0,&VdBf,&dkorB,&tmpv,&tmpv1,&tmpv2,&tmpv3,&tmpv4,NULL);
  malloc_vecs(*antpers,&pbhat,&pbahat,&pahat,&pdbahat,&vone,&dR,&pdA,NULL); 
 
  malloc_vec(*nxval,xt);vone=vec_ones(vone);malloc_vec(*pa,ahatt);malloc_vec(*pa,xai);malloc_vec(*pa,dA); 

  if (*px>=*pa) pmax=*px; else pmax=*pa; 
  if (*nx>=*na) nmax=*nx; else nmax=*na; 
   
  R_CheckUserInterrupt();

  for (s=1;s<*Ntimes;s++)
    {
      time=times[s]; risk=0; dtime=time-times[s-1]; 
      mat_zeros(ldesignX); mat_zeros(ldesignA); vec_zeros(dR); 

      for (c=0,count=0;((c<nmax));c++) 
	{
	  if (status[c]==1) {
	    kia=tukey(stop[c]-time,b[0]);
	    for(j=0;j<*px;j++) 
	      for(i=0;i<*px;i++) ME(XWX,j,i)=ME(XWX,j,i)+
		kia*designX[j*(*nx)+c]*designX[i*(*nx)+c]; }

	  if ((start[c]<time) && (stop[c]>=time)) {
	    for(j=0;j<pmax;j++) {
	      if (j<*px) {ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		if (*mw==1) ME(ldesignX,id[c],j)=
		  ME(ldesignX,id[c],j); } 
	      if (j<*pa) ME(ldesignA,id[c],j)=designA[j*(*na)+c]; }
	    VE(dR,id[c])=y[c];  
	    if (status[c]==1) risk=risk+1; 
	    if (time==stop[c] && status[c]==1) {pers=id[c]; zpers=y[c];
	      /* Rprintf(" %ld %lf %ld \n",s,zpers,pers); */  }
	    count=count+1;} 
	}

      MtA(ldesignA,ldesignA,Aa); invert(Aa,AaI);
      if (ME(AaI,0,0)==0.0) 
	Rprintf("Dynadd: Aalen design not invertible at time %lf",time); 
      extract_row(ldesignA,pers,xai); 
      Mv(AaI,xai,dA); Mv(ldesignA,dA,pdA);

      for (k=0;k<*pa;k++) VE(dAt[s],k)=VE(dA,k); 

      for(j=0;j<*nxval;j++) VE(xt,j)=fabs(ahat[j]-time);
      for(j=1;j<=*pa;j++) VE(ahatt,j-1)=ahat[j*(*nxval)+(*imin)];
      Mv(ldesignA,ahatt,pahat);

      for(j=1;j<=*px;j++) VE(bhatt,j-1)=bhat[j*(*nxval)+(*imin)];
      Mv(ldesignX,bhatt,pbhat);

      for (j=0;j<*antpers;j++) /* sampling corrected design */ 
	{scl_vec_mult(VE(pahat,j),extract_row(ldesignX,j,dB),dB);
	  replace_row(cdesignX,j,dB); VE(pbahat,j)=VE(pbhat,j)*VE(pdA,j); } 

      MtA(cdesignX,ldesignX,A); invert(A,AI);
      if (ME(AI,0,0)==0.0) 
	Rprintf("Dynadd: Regression design not invertible at time %lf",time); 
      extract_row(ldesignX,pers,tmpv); Mv(AI,tmpv,tmpv1); 
      scl_vec_mult(zpers,tmpv1,dB); 

      /* korrektions faktor ================================= */ 
      dif=(zpers-VE(pbhat,pers)); 

      /* Rprintf(" %lf %lf %lf \n",dif,zpers,VE(pbhat,pers));  */ 

      scl_vec_mult(dtime,bhatt,tmpv4); scl_vec_mult(dif,tmpv1,tmpv2); 
      vec_add(tmpv2,tmpv4,dkorB); 

      vM(ldesignX,pbahat,tmpv3); Mv(AI,tmpv3,tmpv4); 
      vec_add(tmpv2,tmpv4,korf); 

      /* LY estimate */ 
      YoneN=vec_sum(dR); 
      vM(ldesignX,pdA,tmpv3); Mv(AI,tmpv3,tmpv); 
      scl_vec_mult(YoneN/risk,tmpv,tmpv); 
      scl_vec_mult(zpers-YoneN/risk,tmpv1,dBly); 
      vec_add(dBly,tmpv,dBly); 

      for (k=1;k<*px+1;k++) {
	cu0[k*(*Ntimes)+s]=VE(dB,k-1)+cu0[k*(*Ntimes)+s-1];
	cuf[k*(*Ntimes)+s]=VE(dkorB,k-1)+cuf[k*(*Ntimes)+s-1];
	cuMS[k*(*Ntimes)+s]=VE(korf,k-1)+cuMS[k*(*Ntimes)+s-1];
	cumly[k*(*Ntimes)+s]=VE(dBly,k-1)+cumly[k*(*Ntimes)+s-1]; }
      cuf[s]=time; cu0[s]=time; cuMS[s]=time; cumly[s]=time; 

      vec_subtr(dB,tmpv4,tmpv3); vec_star(tmpv3,tmpv3,VdB0);  
      vec_star(tmpv2,tmpv2,VdBf); 
      vec_subtr(tmpv4,tmpv,tmpv4); vec_subtr(dBly,tmpv4,tmpv4); vec_star(tmpv4,tmpv4,VdBly); 

      vec_zeros(VdB); 
      for (i=0;i<*antpers;i++) {
	if (i==pers) vec_add(tmpv2,cumhatB[i],cumhatB[i]);
	/* sv_mlt(ahati,tmpv1,tmpv1); v_add(tmpv1,cumA[i],cumA[i]); 
	   v_sub(cumhatA[i],cumA[i],tmpv2); */ 
	vec_star(cumhatB[i],cumhatB[i],tmpv); vec_add(tmpv,VdB,VdB);

	if (*retur==1) 
	  for (k=0;k<*px;k++) {c=i*(*px)+k; 
	    cumBit[c*(*Ntimes)+s]=VE(cumhatB[i],k);} 
	replace_row(cumBt[i],s,cumhatB[i]); 
      } 

      for (k=1;k<*px+1;k++) {
	vcu0[k*(*Ntimes)+s]=VE(VdB0,k-1)+vcu0[k*(*Ntimes)+s-1]; 
	vcuf[k*(*Ntimes)+s]=VE(VdBf,k-1)+vcuf[k*(*Ntimes)+s-1]; 
	vculy[k*(*Ntimes)+s]=VE(VdBly,k-1)+vculy[k*(*Ntimes)+s-1]; 
	robvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      }  
      vcu0[s]=time; vcuf[s]=time; robvcu[s]=time; vculy[s]=time; 
    }

  /* variance comp done based on preliminary estimate of bhat */
//  sing=0; 
//  if (sing==1) { degree[0]=1; coef[0]=1; ps[0]=(*px)+1; 
//    if (vec_sum(pbhat)==0) 
//      smoothB(cumly,Ntimes,ps,bhatny,nxval,b,degree,coef); 
//    else smoothB(cuf,Ntimes,ps,bhatny,nxval,b,degree,coef); 
//  } /* sing==1 */ 

 R_CheckUserInterrupt();

  /* ====================SIMULATIONS ============================= */
  if (*sim==1) {
   comptest(times,Ntimes,px,cuf,robvcu,vcudif,antsim,test,testOBS,Ut,simUt,cumBt,weighted,antpers); 
  } /* sim==1 */ 

  cu0[0]=times[0];vcu0[0]=times[0];robvcu[0]=times[0];cuf[0]=times[0];
  vcuf[0]=times[0]; vculy[0]=times[0]; cumly[0]=times[0]; 

  for (i=0;i<*antpers;i++) {
    free_vec(cumB[i]); free_vec(cumhatB[i]); free_mat(cumBt[i]); }
  for (i=0;i<*Ntimes;i++) free_vec(dAt[i]); 

  free_mats(&XWX,&cdesignA,&ldesignX,&cdesignX,&ldesignA,&Aa,&AaI,&A,&AI,&XbXa,NULL); 
  free_vecs(&VdBly,&korf,&dB,&dA,&dR,&ahatt,&xt,&pdA,&diag,&xai,&sumx,&vone,&itot,&VdB,&VdB0,&VdBf,&fkor,&dkorB,&tmpv,&tmpv1,&tmpv2,&tmpv3,&tmpv4,&pahat,&pbhat,&bhatt,&pbahat,&pdbahat,&dBly,&cum,NULL);
  free(vcudif); 
free(coef); free(ps); free(degree); free(imin); 

}


void semidynadd(times,y,Ntimes,designX,nx,px,designG,ng,pg,designA,na,pa,
ahat,naval,bhat,nxval,antpers,start,stop,cu0,cu,cums,robvcu,robvcue,
gamma,gamma2,gamLY,gamKOR,gameffi,gameffims,
Vgamma,Vkorgam,Vgamef,robvargam,robvargame,w,mw,rani,sim,antsim,
mgresid,test,testOBS,Ut,simUt,b,id,status,weighted,vargamLY,clusters,antclust,resample,
gammaiid,Biid)
double *b,*bhat,*ahat,*designX,*designA,*designG,*times,*y,*start,*stop,*cu,*w,
*gamLY,*gamKOR,*gamma,*gamma2,*gameffi,*gameffims,
*Vgamma,*Vkorgam,*Vgamef,*robvargam,*robvargame,*vargamLY,
*cums,*cu0,*robvcu,*robvcue,*Ut,*simUt,*test,*testOBS,*mgresid,*gammaiid,*Biid; 
int *naval,*nxval,*nx,*px,*na,*pa,*ng,*pg,*antpers,*Ntimes,*mw,
*sim,*antsim,*rani,*id,*status,*weighted,*clusters,*antclust,*resample;
{
  matrix *ldesignX,*ldesignA,*cdesignX,*ldesignG,*cdesignG,*A,*AI;
  matrix *dC,*Cg,*VarG,*dVarG,*VarGly,*dVarGly,*CI,*Vargam,*dCdt,*CGam,*dCGam,*Cgdt;
  matrix *korVarG,*dkorVarG,*CIdt,*C2[*Ntimes],*C2dA[*Ntimes],*dC2dA; 
  matrix *XX,*XXI,*ZX,*cZX,*ZZ; 
  matrix *XZhat,*ZZhat,*Zhat; 
  matrix *Robvar,*RobvarE;
  matrix *Delta,*tmpM1,*tmpM2,*tmpM3,*tmpM4,*tmpsim;
  matrix *W3t[*antpers],*W4t[*antpers],*W4te[*antpers],*AIxit[*antpers];
  vector *W2[*antpers],*W3[*antpers],*W2e[*antpers];
  vector *dB,*dB0,*dN,*dR,*VdB,*VdBe,*ahatt,*bhatt;
  vector *xtb,*xta,*pahat,*pbghat,*pbhat,*pghat,*pbahat,*dA,*pdA;
  vector *gam,*gam2,*gamly,*gamkor,*gamef,*gamefms,*corBg,*dVgam,*pbghatpdA;
  vector *tmpv,*tmpv1,*tmpv2,*tmpv3,*tmpv4,*zi,*xi,*ZHdp,*IZHdp,*IZHdN,*ZHdN;
  vector *korgamly,*korgam,*gamstart;
  vector *rowX,*rowZ,*difX,*dgamef,*gammsd,*ai; 
  int l,i,j,k,s,c,count,sing,pmax,nmax,pers=0;
  int robust=1,*ipers=calloc(*Ntimes,sizeof(int)), 
      *imin=calloc(1,sizeof(int)); 
  double time,dtime,zpers,risk,YoneN,dif,dif2,ctime;
  double *C=calloc((*pg)*(*pg),sizeof(double));
  double *vcudif=calloc((*Ntimes)*(*px+1),sizeof(double));
  void comptest(); 
  ctime=0;  
 
  malloc_mats(*antpers,*px,&ldesignX,&cdesignX,NULL);
  malloc_mats(*antpers,*pg,&Zhat,&ldesignG,&cdesignG,NULL);
  malloc_mats(*antpers,*pa,&ldesignA,NULL);
  malloc_mats(*px,*px,&XX,&XXI,NULL);
  malloc_mats(*pg,*pg,&Robvar,&RobvarE,&ZZhat,&tmpM2,&tmpM3,&ZZ,&dCGam,&CGam,&Cg,&Cgdt,&dC,&CI,&CIdt,&dCdt,&Vargam,&korVarG,&VarGly,&dVarGly,&VarG,&dVarG,&dkorVarG,NULL);

  malloc_mat(*pa,*pa,A); malloc_mat(*pa,*pa,AI); 
  malloc_mat(*pg,*px,tmpM4); malloc_mat(*Ntimes,*px,Delta); 
  malloc_mat(*Ntimes,*px,tmpsim);
  malloc_mat(*pg,*px,ZX); malloc_mat(*pg,*px,cZX); malloc_mat(*pg,*px,tmpM1); 
  malloc_mat(*px,*pg,XZhat); malloc_mat(*px,*pg,dC2dA); 
  for (j=0;j<*Ntimes;j++) {malloc_mat(*px,*pg,C2[j]); malloc_mat(*px,*pg,C2dA[j]);}

  malloc_vecs(*antpers,&pbahat,&pbghatpdA,&pdA,&dN,&dR,&pahat,&pbghat,&pbhat,&pghat,NULL); 
  malloc_vecs(*px,&rowX,&difX,&xi,&tmpv1,&tmpv2,&dB0,&dB,&VdB,
		 &VdBe,&bhatt,&corBg,NULL); // changed length of tmpv
  malloc_vecs(*pa,&ahatt,&dA,&ai,NULL);
  malloc_vecs(*pg,&gammsd,&dgamef,&rowZ,&gamef,&gamefms,&gamly,&gamkor,&korgam,
		 &korgamly,&ZHdN,&IZHdN,&ZHdp,&IZHdp,&zi,&tmpv3,&tmpv4,&dVgam,
		 &gam,&gamstart,&gam2,&tmpv,NULL); // changed length of tmpv
  malloc_vec(*nxval,xtb); malloc_vec(*naval,xta);

  for (j=0;j<*antpers;j++) { malloc_mat(*Ntimes,*px,W3t[j]); 
    malloc_mat(*Ntimes,*px,W4t[j]); malloc_vec(*pg,W2[j]); 
    malloc_mat(*Ntimes,*px,W4te[j]); malloc_vec(*pg,W2e[j]); 
    malloc_vec(*px,W3[j]); malloc_mat(*Ntimes,*px,AIxit[j]); }

  if (*px>=*pa) pmax=*px; else pmax=*pa; if (*pg>=pmax) pmax=*pg; 
  if (*nx>=*na) nmax=*nx; else nmax=*na; 
	 
  /* Prelim. est. of gamma for var. est. loaded from  (B(t)/t */
  for(j=0;j<*pg;j++) {VE(gam,j)=gamma[j];VE(gamstart,j)=gamma[j];}


  R_CheckUserInterrupt();

  for (s=1;s<*Ntimes;s++)
    {
      vec_zeros(dR); sing=0; zpers=0; risk=0;
      time=times[s]; dtime=time-times[s-1];
      mat_zeros(ldesignX); mat_zeros(ldesignG); mat_zeros(ldesignA); 
      ctime=dtime+ctime; 

      for (c=0,count=0;((c<nmax) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<pmax;j++)  {
	    if (j<*px) {ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
	      if (*mw==1) ME(ldesignX,id[c],j)= 
		ME(ldesignX,id[c],j)*sqrt(w[c]);} 
	    if (j<*pg) {ME(ldesignG,id[c],j)=designG[j*(*ng)+c];
	      if (*mw==1) ME(ldesignG,id[c],j)=
		ME(ldesignG,id[c],j)*sqrt(w[c]);} 
	    if (j<*pa)  ME(ldesignA,id[c],j)=designA[j*(*na)+c]; } 
	  VE(dR,id[c])=y[c]; 
	  if (y[c]!=0) risk=risk+1;
	  if (time==stop[c] && status[c]==1) {pers=id[c]; zpers=y[c];} 
	  count=count+1;} }
      ipers[s]=pers; YoneN=vec_sum(dR);  YoneN=YoneN/risk; /* LY korrektion */

      for(j=0;j<*naval;j++) VE(xta,j)=fabs(ahat[j]-time);
      for(j=1;j<=*pa;j++) VE(ahatt,j-1)=ahat[j*(*naval)+(*imin)];
      Mv(ldesignA,ahatt,pahat);

      for(j=0;j<*nxval;j++) VE(xtb,j)=fabs(bhat[j]-time);
      for(j=1;j<=*px;j++) VE(bhatt,j-1)=bhat[j*(*nxval)+(*imin)];
      Mv(ldesignX,bhatt,pbhat); 
      Mv(ldesignG,gam,pghat);
      vec_add(pbhat,pghat,pbghat); 

      MtA(ldesignA,ldesignA,A); invert(A,AI);
      if (ME(AI,0,0)==0.0) 
	Rprintf("Dynadd: Aalen design not invertible at time %lf\n",time); 
      extract_row(ldesignA,pers,ai); 
      Mv(AI,ai,dA);
      Mv(ldesignA,dA,pdA);

      for (j=0;j<*antpers;j++) {
	extract_row(ldesignX,j,dB); scl_vec_mult(VE(pahat,j),dB,dB);
	replace_row(cdesignX,j,dB); 
	extract_row(ldesignG,j,zi); scl_vec_mult(VE(pahat,j),zi,zi);
	replace_row(cdesignG,j,zi); 
	extract_row(ldesignG,j,zi); scl_vec_mult(VE(pdA,j),zi,zi);
	replace_row(Zhat,j,zi);
	VE(pbahat,j)=VE(pbghat,j)*VE(pdA,j);}/*sampling cor. design*/

      MtA(ldesignX,cdesignX,XX); invert(XX,XXI); 
      if (ME(XXI,0,0)==0.0) 
	Rprintf("Dynadd: Regression design not invertible at time %lf\n",time); 
      MtA(ldesignG,ldesignX,ZX); MtA(cdesignG,ldesignX,cZX); 
      MtA(ldesignG,cdesignG,ZZ); 
      MxA(cZX,XXI,tmpM1); MAt(tmpM1,cZX,tmpM2);

      mat_subtr(ZZ,tmpM2,dCGam);scl_mat_mult(dtime,dCGam,dCGam);
      mat_add(dCGam,CGam,CGam); 

      MtA(ldesignG,Zhat,ZZhat); MtA(ldesignX,Zhat,XZhat); 
      MxA(tmpM1,XZhat,tmpM3); mat_subtr(ZZhat,tmpM3,dC); mat_add(dC,Cg,Cg); 

      extract_row(ldesignG,pers,zi); extract_row(ldesignX,pers,xi); 

      Mv(tmpM1,xi,tmpv3); vec_subtr(zi,tmpv3,ZHdN); 
      scl_vec_mult(zpers,ZHdN,ZHdp); vec_add(ZHdp,IZHdp,IZHdp); 
      if (s<0) { Rprintf(" %d %d ±n",pers,s);print_vec(zi);print_vec(xi);} 

      /* LY korrektion */ 
      scl_vec_mult(zpers-YoneN,ZHdN,korgamly); 

      vM(ldesignG,pdA,rowZ); vM(ldesignX,pdA,rowX); 
      Mv(tmpM1,rowX,zi); 
      vec_subtr(rowZ,zi,zi); scl_vec_mult(YoneN,zi,tmpv); 

      vec_add(korgamly,tmpv,korgamly); vec_add(korgamly,gamly,gamly); 

      /* MS korrektion */ 
      dif=zpers-VE(pbghat,pers); dif2=dif*dif;
      scl_vec_mult(dif,ZHdN,gamkor); vec_add(gamkor,korgam,korgam); 

      vM(ldesignX,pbahat,tmpv1); vM(ldesignG,pbahat,tmpv4); 
      Mv(tmpM1,tmpv1,rowZ); vec_subtr(tmpv4,rowZ,tmpv4); 

      vec_add(tmpv4,gamkor,rowZ); vec_add(rowZ,gam2,gam2); 

      /* Efficient estimates when sigma^2(s) depends on s */
      /*
	scl_mat_mult(1/dtime,dCGam,dCGam); 
	if (problems==0) { 
	invert(dCGam,tmpM3); 
	if (ME(tmpM3,0,0)==0.0) {problems=1; vec_zeros(gamef);
	Rprintf("Semi-Dynadd: gamma.ef contains non-invertible\n");} 
	Mv(tmpM3,gamkor,dgamef); vec_add(dgamef,gamef,gamef); 
	}
	vec_add(gamkor,zi,zi);  Mv(tmpM3,zi,rowZ); 
	vec_add(rowZ,gamefms,gamefms); 
      */

      /*  preparing nonparametric estimates  */ 
      extract_row(ldesignX,pers,xi); Mv(XXI,xi,dB0); 
      scl_vec_mult(dif,dB0,dB); scl_vec_mult(dtime,bhatt,tmpv2); vec_add(dB,tmpv2,tmpv2);

      Mv(XXI,tmpv1,rowX); 
      vec_add(dB,rowX,tmpv1); 

      for (k=0;k<*px;k++) {
	cu0[(k+1)*(*Ntimes)+s]=zpers*VE(dB0,k)+cu0[(k+1)*(*Ntimes)+s-1];
	cu[(k+1)*(*Ntimes)+s]= VE(tmpv2,k)+ cu[(k+1)*(*Ntimes)+s-1];
	cums[(k+1)*(*Ntimes)+s]=VE(tmpv1,k)+ cums[(k+1)*(*Ntimes)+s-1];}
      cu0[s]=time; cu[s]=time; cums[s]=time; 

      MAt(XXI,cZX,dC2dA); scl_mat_mult(dtime,dC2dA,dC2dA); 
      mat_add(dC2dA,C2[s-1],C2[s]);
      MxA(XXI,XZhat,dC2dA);  mat_add(dC2dA,C2dA[s-1],C2dA[s]); 

      /*mindividual MG resid gamma term */ 
      /* terms for robust variance   */ 
      for (i=0;i<*antpers;i++) {
	if (i==ipers[s]) { 
	  vec_add(gamkor,W2[i],W2[i]);
	  /* Mv(tmpM3,gamkor,rowZ); vec_add(rowZ,W2e[i],W2e[i]); */
	  vec_add(dB,W3[i],W3[i]); }
	replace_row(W3t[i],s,W3[i]);  
      } /* i=1..antpers */ 

      /* variance calculations based on non-parametric model (bhat) */
      vec_subtr(ZHdp,tmpv4,tmpv4); /* variance of gamma0 */
      scl_vec_mult(zpers-YoneN,ZHdN,zi); 
      vec_subtr(zi,tmpv,tmpv); /* variance of gammaLY */

      for (k=0;k<*pg;k++) 
	for (j=0;j<*pg;j++) {
	  ME(dVarG,k,j)=VE(tmpv4,k)*VE(tmpv4,j);
	  ME(dVarGly,k,j)=VE(tmpv,k)*VE(tmpv,j);
	  ME(dkorVarG,k,j)=dif2*VE(ZHdN,k)*VE(ZHdN,j); 
	  Vgamef[k*(*pg)+j]=Vgamef[k*(*pg)+j]+VE(dgamef,k)*VE(dgamef,j);
	} 
      mat_add(VarG, dVarG, VarG); mat_add(VarGly, dVarGly, VarGly); 
      mat_add(korVarG, dkorVarG, korVarG); 
    }   /* for s in Ntimes -------------------------------------------   */
  /* if (*gamdt==1) {invert(CGam,CIdt); Mv(CIdt,IZHdp,gam2);}*/ 

  R_CheckUserInterrupt();
  
  invert(Cg,CI); 
  Mv(CI,IZHdp,gam);
  Mv(CI,gamly,korgamly);
  Mv(CI,gam2,gammsd); 

  Mv(CI,korgam,tmpv3);     vec_add(gamstart,tmpv3,gamkor); 
  scl_vec_mult(1/ctime,gamef,gamef); vec_add(gamstart,gamef,gamef); 

  MxA(CI,VarG,tmpM3); MAt(tmpM3,CI,VarG);  
  MxA(CI,VarGly,tmpM3); MAt(tmpM3,CI,VarGly);  
  MxA(CI,korVarG,tmpM3); MAt(tmpM3,CI,korVarG);  

  for (k=0;k<*pg;k++) {gamma[k]=VE(gam,k); 
    gamLY[k]=VE(korgamly,k); 
    gamKOR[k]=VE(gamkor,k); 
    gamma2[k]=VE(gammsd,k); 
    gameffi[k]=VE(gamef,k); 
    gameffims[k]=VE(gamefms,k)/ctime; 
    for (j=0;j<*pg;j++) { C[k*(*pg)+j]=ME(CI,j,k);
      Vgamma[k*(*pg)+j]= ME(VarG,j,k); 
      vargamLY[k*(*pg)+j]= ME(VarGly,j,k); 
      Vkorgam[k*(*pg)+j]= ME(korVarG,j,k); 
      Vgamef[k*(*pg)+j]=Vgamef[k*(*pg)+j]/(ctime*ctime); 
    } } 

  /* =========================================================== */
  /* Robust variances and Estimates of cumulative reg. functions */ 

  R_CheckUserInterrupt();

  if (robust==1) {
    for (s=1;s<*Ntimes;s++) {

      time=times[s]; 
      /* vec_subtr(gamstart,gamkor,tmpv3); */
      Mv(C2dA[s],gamef,corBg);  
      for (k=0;k<*px;k++) {
	cu[(k+1)*(*Ntimes)+s]=cu[(k+1)*(*Ntimes)+s]-VE(corBg,k);
	cu0[(k+1)*(*Ntimes)+s]=cu0[(k+1)*(*Ntimes)+s]-VE(corBg,k);
	cums[(k+1)*(*Ntimes)+s]=cums[(k+1)*(*Ntimes)+s]-VE(corBg,k);}
      cu[s]=time; cums[s]=time; 

      vec_zeros(VdB); 
      for (i=0;i<*antpers;i++) {
        Mv(CI,W2[i],tmpv3); // Note that tmpv2 was changed to tmpv3
	Mv(C2[s],tmpv3,rowX); // Note that tmpv2 was changed to tmpv3 

	extract_row(W3t[i],s,tmpv1); vec_subtr(tmpv1,rowX,difX); 
	replace_row(W4t[i],s,difX); 
	vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdB,VdB);

	if (*resample==1) {
	  if (s==1){ 
	    for (k=0;k<*pg;k++) gammaiid[k*(*antpers)+i]=VE(tmpv3,k); 
	  }
	  for (k=0;k<*px;k++) {
	    l=i*(*px)+k; 
	    Biid[l*(*Ntimes)+s]=VE(difX,k);
	  } 
	}


	/*
	  Mv(C2dA[s],W2e[i],rowX); extract_row(W3t[i],s,tmpv1); 
	  vec_subtr(tmpv1,rowX,difX); replace_row(W4te[i],s,difX); 
	  vec_star(difX,difX,tmpv1); vec_add(tmpv1,VdBe,VdBe);
	*/

	if (s==1) { for (j=0;j<*pg;j++) for (k=0;k<*pg;k++)  {
		      ME(RobvarE,j,k)=ME(RobvarE,j,k)+VE(W2e[i],j)*VE(W2e[i],k);
		      ME(Robvar,j,k)=ME(Robvar,j,k)+VE(W2[i],j)*VE(W2[i],k);} }

      } /* i =1 ..Antpers */ 

      /*
	MAt(Vargam,C[s],tmpM4); MxA(C[s],tmpM4,VarKorG);
	MAt(CI, M1M2[s], tmpM4); MxA(C[s], tmpM4, GCdM1M2);
	for (k=1;k<=*px;k++) vcu[k*(*Ntimes)+s]= 
	vcu[k*(*Ntimes)+s]+ME(VarKorG,k-1,k-1)-2*ME(GCdM1M2,k-1,k-1);
      */

      for (k=1;k<*px+1;k++) robvcu[k*(*Ntimes)+s]=VE(VdB,k-1);
      for (k=1;k<*px+1;k++) robvcue[k*(*Ntimes)+s]=VE(VdBe,k-1);
      robvcu[s]=times[s]; robvcue[s]=times[s];
    } /* s=1 ..Ntimes */ 

    MxA(Robvar,CI,tmpM2); MxA(CI,tmpM2,Robvar); 

    for (j=0;j<*pg;j++) for (k=0;k<*pg;k++) {
      robvargam[k*(*pg)+j]=ME(Robvar,j,k);
      robvargame[k*(*pg)+j]=ME(RobvarE,j,k)/(ctime*ctime);}

  } /* robust==1 */ 

  R_CheckUserInterrupt();

  if (*sim==1) {
    comptest(times,Ntimes,px,cu,robvcu,vcudif,antsim,test,testOBS,Ut,simUt,W4t,weighted,antpers);
  } /* sim==1 */ 

  /* Free memory ---------------------------------------------------- */
  for (j=0;j<*antpers;j++) {
    free_mat(W3t[j]); free_mat(W4t[j]); free_mat(W4te[j]); free_mat(AIxit[j]); 
    free_vec(W2[j]); free_vec(W3[j]); free_vec(W2e[j]); }
  for (j=0;j<*Ntimes;j++) {free_mat(C2[j]); free_mat(C2dA[j]); }

  free_mats(&ldesignX,&ldesignA,&cdesignX,&A,&AI,&ldesignG,&cdesignG,
	      &dC,&Cg,&VarG,&dVarG,&VarGly,&dVarGly,&CI,&Vargam,&dCdt,&CGam,&dCGam,&Cgdt,
	      &korVarG,&dkorVarG,&CIdt,&dC2dA, 
	      &XX,&XXI,&ZX,&cZX,&ZZ,&XZhat,&ZZhat,&Zhat,&Robvar,&RobvarE,
	      &Delta,&tmpsim,&tmpM1,&tmpM2,&tmpM3,&tmpM4,NULL);

  free_vecs(&dB0,&ai,&gammsd,&dgamef,&dB,&dN,&dR,&VdB,&VdBe,
	      &ahatt,&bhatt,&xtb,&xta,&pahat,&pbghat,&pbhat,&pghat,&pbahat,&dA,&pdA,
	      &gam,&gam2,&gamly,&gamkor,&gamef,&gamefms,&corBg,&dVgam,&pbghatpdA,
	      &tmpv,&tmpv1,&tmpv2,&tmpv3,&tmpv4,&zi,&xi,&ZHdp,&IZHdp,&IZHdN,&ZHdN,
	      &korgamly,&korgam,&gamstart,
	      &rowX,&rowZ,&difX,NULL); 

  free(C); free(vcudif); free(ipers); 
  free(imin); 
}
