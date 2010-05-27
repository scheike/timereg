#include <stdio.h>
#include <math.h>
#include "matrix.h"

void l1boost(D,px,d,lambda,Nit,beta,detail)
double *D,*d,*beta,*lambda; 
int *px,*detail,*Nit;
{
  matrix *iD;
  vector *Db,*betav,*dL,*tmpv1;
  int index,it,i,j;
  double k,lg,bDb,bd,dummy,fabs(),sqrt(),L0,L1,Lk,fdL; 

  malloc_mats(*px,*px,&iD,NULL);
  malloc_vecs(*px,&tmpv1,&betav,&dL,&Db,&bd,NULL);

  for (i=0;i<*px;i++) { 
    VE(betav,i)=beta[i]; 
    for(j=0;j<*px;j++) {
      ME(iD,i,j)=D[j*(*px)+i]; 
    }
  }

  for (it=0;it<*Nit;it++){
    Mv(iD,betav,Db); 
    bd=0; bDb=0; 
    for (i=0;i<*px;i++) { 
      bd=bd+VE(betav,i)*d[i]; 
      bDb=bDb+VE(Db,i)*VE(betav,i); 
    }

    index=0; dummy=0; 
    for (i=0;i<*px;i++) {
      VE(dL,i)=d[i]-VE(Db,i); 
      fdL=fabs(VE(dL,i)); 
      if (fdL>dummy)  {dummy=fdL; index=i; }
    }
    if (*detail==1) printf(" %ld \n",(long int) index); 

    if (VE(dL,index)<0) lg=-1*(*lambda); else lg=*lambda; 

    if (fabs(VE(dL,index))<.00000000001) {it=*Nit-1; break;}

    k=(VE(Db,index)*lg-bDb+bd-lg*d[index]) /(-bDb-lg*lg*ME(iD,index,index)+2*lg*VE(Db,index)); 
    if (*detail==1) printf(" %lf %lf \n",k,lg); 
    if (*detail==1) printf(" %lf %lf \n",bDb,bd); 
    if (*detail==1) printf(" %lf %lf %lf %lf %lf %lf %lf %lf \n",
			   lg,ME(iD,index,index),d[index],
			   0.5*lg*lg*ME(iD,index,index)-lg*d[index],
			   0.5*lg*lg*ME(iD,index,index),
			   0.5*lg*lg,ME(iD,index,index),
			   -lg*d[index]); 

    L1=0.5*lg*lg*ME(iD,index,index)-lg*d[index]; 
    L0=0.5*bDb-bd;
    Lk=0.5*((1-k)*(1-k)*bDb+k*k*lg*lg*ME(iD,index,index)+2*k*(1-k)*lg*VE(Db,index))-(1-k)*bd-k*lg*d[index];

    if (L0<=Lk && L0<=L1) k=0; if (L1<=Lk && L1<=L0) k=1; 

    if (*detail==1) printf(" %lf %lf %lf %lf \n",L0,L1,Lk,k); 

    for (i=0;i<*px;i++) VE(betav,i)*=(1-k); 
    VE(betav,index)+=k*lg;
  }

  for (i=0;i<*px;i++) beta[i]=VE(betav,i); 

  free_mats(&iD,NULL);
  free_vecs(&tmpv1,&betav,&dL,&Db,NULL);
}

void addiboost(D,px,d,Nit,beta,detail,index,step,var,met)
double *D,*d,*beta,*step,*var; 
int *met,*px,*detail,*Nit,*index;
{
  int k,itta,pls; 
  double fabs(),sqrt(),L,Lp=0,pval,pvalp=0;
  double gam,cgam,kor,izhdn; 

  for (itta=0;itta<*Nit-1;itta++)
    {
      for (pls=0;pls<*px;pls++)
	{
	  cgam=D[pls*(*px)+pls]; izhdn=d[pls]; 

	  kor=0; 
	  for (k=0;k<itta;k++){kor=kor+D[index[k]*(*px)+pls]*beta[k];}
	  izhdn=izhdn-step[0]*kor; 
	  gam=izhdn/cgam; 

	  L=0.5*gam*gam*cgam-gam*izhdn; 
	  pval=gam*gam/var[pls];  

	  if (pls==0) {Lp=L+1;pvalp=pval-1; }
	  if (*met==0)
	    if (L< Lp){beta[itta]=gam; index[itta]=pls; Lp=L;}
	  if (*met==1)
	    if (pval> pvalp){beta[itta]=gam;index[itta]=pls;pvalp=pval;}
	}
    } /* dimcovpls comps */ 
}

/*
void Lc(D,px,d,Dd,db,beta,bDb,L)
MAT *D; 
VEC *Dd,*beta; 
int *px; 
double *db,*bDb,*d,*L; 
{
mv_mlt(D,betav,Db); 
bDb[0]=0; db[0]=0; 
for (i=0;i<*px;i++) { 
bd[0]=bd[0]+beta->ve[i]*d[i]; bDb[0]=bDb[0]+Db->ve[i]*betav->ve[i]; }

L[0]=0.5*bDb[0]-bd[0]; 
}
*/
