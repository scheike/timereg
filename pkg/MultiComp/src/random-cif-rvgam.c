#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "multicomp.h"
                 
void rcifdes(times,Ntimes,x,delta,cause,CA1,KMc,z,antpers,px,Nit,score,hess,est,
gamma,semi,zsem,pg,detail,biid,gamiid,timepow,theta,vartheta,thetades,ptheta,antclust,
cluster,clustsize,clusterindex,maxclust,step,inverse,dscore,rvdes,prv,notaylor,samecens,
trunkp, entryage,cif1lin 
)
double *theta,*times,*x,*KMc,*z,*score,*hess,*est,*gamma,*zsem,*vartheta,*biid,*gamiid,*timepow,*thetades,*step,*rvdes, 
       *trunkp, *entryage,*cif1lin ; 
int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*semi,*pg,*CA1,*detail,*ptheta,
*antclust,*cluster,*clustsize,*clusterindex,*maxclust,*inverse,*dscore,*prv,*notaylor,*samecens;
{ // {{{
// {{{
 matrix *ldesignX,*A,*AI,*cdesignX,*ldesignG,*cdesignG;
 matrix *S,*dCGam,*CGam,*ICGam,*VarKorG,*XZ,*ZZ,*ZZI,*XZAI;
 matrix *Vargam,*dVargam,*dM1M2,*M1M2t,*RobVargam;
 matrix *DUeta[*Ntimes],*DUgamma,*RVdes;
 matrix *Biid[*antclust]; 
 matrix *pardes[*antpers]; 
 matrix *destheta,*d2UItheta,*d2Utheta,*varthetascore,*Stheta; // *Sthetaiid[*antclust]; 
 vector *W2[*antclust];
 vector *gammaiid[*antclust], *W3[*antclust], *W4[*antclust]; 
 vector *diag,*dB,**VdB,*AIXdN,*AIXlamt,*bhatt,*pbhat,*lamtt;
 vector *korG,*pghat0,*pghat,*gam;
 vector *xk,*xi,*rowX,*rowZ,*difX,*zk,*zi,*z1;
 vector *Utheta,*vthetascore,*vtheta1,*vtheta2,*vtheta3,*dtheta; //*thetaiid[*antclust]; 
 vector *alphai,*alphaj,*alphatot,*rvvec,*rvvec1,
	*rvvec2, *rvvec2vt, *rvvec2tv,
	*rvvec2vv;

 int naprint=0,pmax,v,itt,i,j,k,l,l1,s,c;
 double Li,Lk,ithetak=0,thetak=0,response,time,dtime;
 double edd,test,fabs(),Dinverse,DDinverse;
 double sumscore,sdj,pow(),diff,
        *ckij=calloc(1,sizeof(double)),
        *ckijvv=calloc(1,sizeof(double)), *ckijtv=calloc(1,sizeof(double)),
        *ckijvt=calloc(1,sizeof(double)), *dckij=calloc(1,sizeof(double));
 void LevenbergMarquardt(),ckrvdes2(),DUetagammarv(); 

 test=1; 
	
//if (*trans==1) for (j=0;j<*pg;j++) if (timepow[j]!= 1) {timem=1;break;}
//if (*trans==2) for (j=0;j<*pg;j++) if (timepow[j]!= 0) {timem=1;break;}

  malloc_mat(*antpers,*prv,RVdes); 
  malloc_vecs(*prv,&rvvec,&rvvec1,&rvvec2,&rvvec2vv,&rvvec2vt,&rvvec2tv,&alphai,&alphaj,&alphatot,NULL); 

//Rprintf("m1 hej med mig \n"); 

  if ((*notaylor==1)) 
  for (j=0;j<*Ntimes;j++) { malloc_mat(*px,*ptheta,DUeta[j]); }
  for (j=0;j<*antclust;j++) { 

   if ((*notaylor==1)) {
      malloc_vec(*pg,gammaiid[j]); 
      malloc_mat(*Ntimes,*px,Biid[j]); 
      malloc_vec(*ptheta,W3[j]); malloc_vec(*ptheta,W4[j]); 
   }
   malloc_vec(*ptheta,W2[j]); 
//   malloc_vec(*ptheta,thetaiid[j]); malloc_mat(*ptheta,*ptheta,Sthetaiid[j]); 
   }

//Rprintf("m2 hej med mig \n"); 

  malloc_mat(*antpers,*ptheta,destheta); malloc_mat(*ptheta,*ptheta,Stheta); 
  malloc_mats(*ptheta,*ptheta,&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  malloc_mats(*antpers,*px,&ldesignX,&cdesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,&cdesignG,NULL); 
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*pg,*pg,&dVargam,&Vargam,&RobVargam,&ZZ,&VarKorG,&ICGam,&CGam,&dCGam,
		      &S,&ZZI,NULL); 
  malloc_mats(*px,*pg,&XZAI,&XZ,&dM1M2,&M1M2t,NULL);
  malloc_mats(*pg,*ptheta,&DUgamma,NULL);
  for (c=0;c<*antpers;c++) malloc_mat(*prv,*ptheta,pardes[c]); 

  malloc_vecs(*ptheta,&Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta3,&vtheta2,NULL);
  malloc_vecs(*px,&xk,&xi,&rowX,&difX,&korG,&diag,&dB,&VdB,&AIXdN,
		  &AIXlamt,&bhatt,NULL);
  malloc_vecs(*pg,&zk,&zi,&rowZ,&z1,&gam,NULL);
  malloc_vecs(*antpers,&pbhat,&pghat0,&pghat,&lamtt,NULL);
  
 
  pmax=max(*px,*pg); pmax=max(pmax,*ptheta); 
  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 
  for (j=0;j<*ptheta;j++) VE(vtheta1,j)=theta[j]; 

  if ((*notaylor==1)) {
  if (*semi==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++)
        VE(gammaiid[i],j)=gamiid[j*(*antclust)+i]; 
  if (*semi==2) for (i=0;i<*antclust;i++) print_vec(gammaiid[i]); 

   for (i=0;i<*antclust;i++) {
   for (s=0;s<*Ntimes;s++) 
   for (c=0;c<*px;c++) {l=i*(*px)+c; ME(Biid[i],s,c)=biid[l*(*Ntimes)+s]; }
   }
   }

    for (c=0;c<*antpers;c++) {
//	if (j<*ptheta) ME(destheta,c,j)= thetades[j*(*antpers)+c];
      for(j=0;j<*prv;j++)   {
         ME(RVdes,c,j)= rvdes[j*(*antpers)+c];
         for(i=0;i<*ptheta;i++)  { l=i*(*prv)+j; 
	      ME(pardes[c],j,i)=thetades[l*(*antpers)+c];
         }
      }
      for(j=0;j<pmax;j++)  {
	if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c];
        if (*semi==1) if (j<*pg) { ME(ldesignG,c,j)=zsem[j*(*antpers)+c];}
      }
   } 
// }}}

if (*semi==1) Mv(ldesignG,gam,pghat0);

if (test<0) {
    Rprintf(" theta \n"); print_vec(vtheta1); 
    head_matrix(RVdes); 
    for (c=0;c<10;c++) print_mat(pardes[c]); 
}

  for (itt=0;itt<*Nit;itt++) // {{{
  { 
      R_CheckUserInterrupt();
      mat_zeros(d2Utheta); vec_zeros(Utheta); 
      sumscore=0; Dinverse=1; DDinverse=1; 

      if (*inverse==1)  {
            for (l1=0;l1<*ptheta;l1++) VE(vtheta2,l1)= exp(VE(vtheta1,l1)); 
//            for (l1=0;l1<*ptheta;l1++) VE(vtheta2,l1)= pow(VE(vtheta1,l1),2); 
	    } else scl_vec_mult(1,vtheta1,vtheta2); 

      for (s=0;s<*Ntimes;s++)
      {
     // Rprintf("times  s %d %d %d \n",s,*Ntimes,*antclust); 
	  time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 
	  for(j=1;j<=*px;j++) {VE(bhatt,j-1)=est[j*(*Ntimes)+s];}
	  Mv(ldesignX,bhatt,pbhat); 
	  if (*semi==1) {scl_vec_mult(time,pghat0,pghat);vec_add(pbhat,pghat,pbhat);}
	
    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
          diff=0;sdj=0; vec_zeros(rowX);vec_zeros(rowZ); 

          for (c=0;c<clustsize[j];c++) for (v=0;v<clustsize[j];v++) 
	  if (v!=c) { 
	    vec_zeros(rvvec2); 
	    vec_zeros(rvvec2vv); vec_zeros(rvvec2tv); vec_zeros(rvvec2vt); 
	    i=clusterindex[c*(*antclust)+j]; k=clusterindex[v*(*antclust)+j];

 if ((entryage[i] < time) && (entryage[k]< time)) { 
            Mv(pardes[i],vtheta2,alphai); Mv(pardes[k],vtheta2,alphaj); 

       if ((test<0) & (j==1) & (s==1)) {
	  Rprintf(" pars =========================\n"); 
          print_vec(vtheta1); print_vec(vtheta2); 
	  print_mat(pardes[i]); 
          print_vec(alphai); print_vec(alphaj); 
       }

     response=(((x[i]<=time) & (cause[i]==*CA1))*1*((x[k]<=time) & (cause[k]==*CA1))*1); 

    if (*samecens==1) response=response/min(KMc[i],KMc[k]); 
                 else response=response/(KMc[i]*KMc[k]);

       Li=VE(pbhat,i); Lk=VE(pbhat,k); 
       extract_row(RVdes,i,rvvec); extract_row(RVdes,k,rvvec1); 


if (trunkp[i]<1) {
       ckrvdes2(alphai,alphaj,1.0,Li,Lk,ckij,rvvec2,rvvec,rvvec1); 
       ckrvdes2(alphai,alphaj,1.0,cif1lin[i],cif1lin[k],ckijvv,rvvec2vv,rvvec,rvvec1); 
       ckrvdes2(alphai,alphaj,1.0,Li,cif1lin[k],ckijtv,rvvec2tv,rvvec,rvvec1); 
       ckrvdes2(alphai,alphaj,1.0,cif1lin[i],Lk,ckijvt,rvvec2vt,rvvec,rvvec1); 
       vec_zeros(rvvec); 
//	  ddd=(dckij[0]+dckijvv[0]-dckijtv[0]-dckijvt[0])/trunkp[i]; 
	  edd=(ckij[0]+ckijvv[0]-ckijtv[0]-ckijvt[0])/trunkp[i]; 
	  vec_add(rvvec2vv,rvvec2,rvvec2); 
	  vec_subtr(rvvec2,rvvec2tv,rvvec2); 
	  vec_subtr(rvvec2,rvvec2vt,rvvec2); 
          diff=(response-edd); 
         vM(pardes[i],rvvec2,vthetascore); 
        } else {
       ckrvdes2(alphai,alphaj,1.0,Li,Lk,ckij,rvvec2,rvvec,rvvec1); 
       diff=(response-ckij[0]); 
       vM(pardes[i],rvvec2,vthetascore); 
       }



//       ckrvdes2(alphai,alphaj,1.0,Li,Lk,ckij,rvvec2,rvvec,rvvec1); 
//       diff=(response-ckij[0]); 
//       vM(pardes[i],rvvec2,vthetascore); 

	if (isnan(response))   { // removes these from score equations
	   diff=0; sdj=0; 
	}


       if (*inverse==1)  
//  for (l1=0;l1<*ptheta;l1++) VE(vthetascore,l1)=VE(vthetascore,l1)*2*VE(vtheta1,l1); 
for (l1=0;l1<*ptheta;l1++) VE(vthetascore,l1)=VE(vthetascore,l1)*VE(vtheta2,l1); 

       if (test<0) {Rprintf("thomas lomas lomas \n"); 
       print_vec(rvvec2); print_mat(pardes[i]); print_vec(vthetascore); 
       }

	if ((isnan(response)) && (naprint==0))   { // {{{ print diverse na information
	  	 Rprintf(" Missing values, probably censoring dist. problems\n"); 
	         Rprintf(" %d %d %d %d \n",j,clustsize[j],i,k); 
                 Rprintf(" %lf %lf %lf \n",x[i],x[k],time); 
		 Rprintf(" resp %lf  %lf %lf  \n",response,KMc[i],KMc[k]); 
                 Rprintf(" %d %d  \n",i,k); 
                 Rprintf(" %lf %lf %lf \n",thetak,Li,Lk); 
                 Rprintf(" %lf %lf \n",dckij[0],ckij[0]); 
		 print_vec(rvvec2); 
		 naprint=1; 
		 Rprintf("============================== \n"); 
	 } // }}}

	    
       	if ((itt==*Nit-1) & (*notaylor==1)) {
	   extract_row(ldesignX,i,xi);extract_row(ldesignX,k,xk); 
           DUetagammarv(ithetak,Li,Lk,xi,xk); 
	   if (*dscore==1) scl_vec_mult(dckij[0],xi,xi); 
	   vec_add(xi,rowX,rowX); 
           if (*semi==1)  
	   {
               extract_row(ldesignG,i,zi);extract_row(ldesignG,k,zk); 
               DUetagammarv(ithetak,Li,Lk,zi,zk); 
	       if (*dscore==1) scl_vec_mult(dckij[0],zi,zi); 
	       scl_vec_mult(time,zi,zi); vec_add(zi,rowZ,rowZ); 
	   } 
	}

      scl_vec_mult(diff,vthetascore,vtheta3); 
      vec_add(vtheta3,Utheta,Utheta); 

      if (itt==*Nit-1) {
	  vec_add(vtheta3,W2[j],W2[j]);
       	if ((*notaylor==1)) {
         for (l=0;l<*px;l++) for (l1=0;l1<*ptheta;l1++) 
	   ME(DUeta[s],l,l1)+=VE(rowX,l)*VE(vthetascore,l1); 
	   if (*semi==1) for (l=0;l<*pg;l++) for (l1=0;l1<*ptheta;l1++) 
	     ME(DUgamma,l,l1)+=VE(rowZ,l)*VE(vthetascore,l1);
	}
      }

      for (l=0;l<*ptheta;l++) 
      for (l1=0;l1<*ptheta;l1++) {
//      if (itt==*Nit-1) ME(Sthetaiid[j],l,l1)+=VE(vthetascore,k)*VE(vthetascore,l1);  
	 ME(d2Utheta,l,l1)+= VE(vthetascore,l)*VE(vthetascore,l1);
      }
 } // entryage
  } /* for (c=0....... */  
 } /* j in antclust */ 

    if ((itt==*Nit-1) & (*notaylor==1)) {
        for (j=0;j<*antclust;j++) {
       extract_row(Biid[j],s,rowX); 
       //Rprintf(" %ld %ld  \n",s,j); print_vec(rowX); 
       vM(DUeta[s],rowX,dtheta); vec_add(dtheta,W3[j],W3[j]);}
    }


 } /* s=1,...Ntimes */

//   Rprintf("Score D l\n"); print_vec(Utheta);
ckij[0]=1000; 
   LevenbergMarquardt(d2Utheta,d2UItheta,Utheta,dtheta,step,step);

  if (ME(d2UItheta,0,0)==0.0 ){ 
       Rprintf(" second derivative not invertible at iteration %d \n",itt); 
       itt=*Nit; 
  }


//  invert(d2Utheta,d2UItheta); Mv(d2UItheta,Utheta,dtheta);
//  scl_vec_mult(step[0],dtheta,dtheta); 

 if (*detail==1) {
    Rprintf("===============Iteration %d ==================== \n",itt);
     Rprintf("Estimate theta \n"); print_vec(vtheta1);
     Rprintf("Score D l\n"); print_vec(Utheta);
     Rprintf("D Score \n"); print_mat(d2Utheta);
     Rprintf("Information D^2 l\n"); print_mat(d2UItheta); }

    for (k=0;k<*ptheta;k++) sumscore= sumscore+fabs(VE(Utheta,k));

    if ((sumscore<0.000001) & (itt<*Nit-2)) itt=*Nit-2;
    if (isnan(vec_sum(Utheta)))  itt=*Nit-1;

    vec_add(vtheta1,dtheta,vtheta1); // adding dtheta

} // }}} /*itt løkke */ 

   R_CheckUserInterrupt();

   vec_zeros(dtheta); 
   for (j=0;j<*antclust;j++)  // {{{
   {
	   
   if ((*notaylor==1)) {
  // Rprintf("W2  ========= \n"); print_vec(W2[j]); 
  // Rprintf("W3  ========= \n"); print_vec(W3[j]); 
   // vec_star(W3[j],dtheta); 
//    vec_subtr(W2[j],W3[j],W2[j]); 
if (*semi==1) { //Rprintf(" =W4======== \n"); print_vec(W4[j]); 
               vM(DUgamma,gammaiid[j],W4[j]); 
//               vec_subtr(W2[j],W4[j],W2[j]);
              }
    }
    for (k=0;k<*ptheta;k++) 
    for (c=0;c<*ptheta;c++) 
    ME(varthetascore,k,c)=ME(varthetascore,k,c)+VE(W2[j],c)*VE(W2[j],k); 
   } // }}}

   MxA(varthetascore,d2UItheta,d2Utheta); 
   MxA(d2UItheta,d2Utheta,varthetascore);
   // print_mat(varthetascore);

  for (j=0;j<*ptheta;j++) {theta[j]=VE(vtheta1,j); score[j]=VE(Utheta,j);
    for (k=0;k<*ptheta;k++) {vartheta[k*(*ptheta)+j]=ME(varthetascore,j,k);
                             hess[k*(*ptheta)+j]=ME(d2UItheta,j,k);}}

// {{{ freeing 
			     
   if ((*notaylor==1)) for (j=0;j<*Ntimes;j++) { free_mat(DUeta[j]); }

  for (j=0;j<*antclust;j++) {
  if ((*notaylor==1)) {
     free_vec(gammaiid[j]); free_mat(Biid[j]); 
     free_vec(W3[j]); free_vec(W4[j]); 
  }
  free_vec(W2[j]); 
  }
  for (c=0;c<*antpers;c++) free_mat(pardes[c]);

free_mats(&RVdes,&destheta,&Stheta,
  &varthetascore,&d2Utheta,&d2UItheta,
  &ldesignX,&cdesignX, &ldesignG,&cdesignG, &A,&AI,
  &dVargam,&Vargam,&RobVargam,&ZZ,&VarKorG,&ICGam,&CGam,&dCGam,&S,&ZZI,
  &XZAI,&XZ,&dM1M2,&M1M2t,&DUgamma,NULL);

free_vecs(&xk,&xi,&rowX,&difX,&korG,&diag,&dB,&VdB,&AIXdN,&AIXlamt,&bhatt,
  &zk,&zi,&rowZ,&z1,&gam,
  &pbhat,&pghat0,&pghat,&lamtt,
  &Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,&vtheta3,NULL);
  free_vecs(&rvvec2tv,&rvvec2vt,&rvvec2vv,&rvvec2,&rvvec,&rvvec1,&alphai,&alphaj,&alphatot,NULL); 

  free(ckij); free(ckijvv); free(ckijtv); free(ckijvt); 
  free(dckij); // }}}
} // }}}

void ckrvdes(vector *alphai,vector *alphak, // {{{
		double beta, double x,double y,double *ckij,
		vector *dckij,vector *rvi,vector *rvk)
{
double val,val1,val2,val3,alphi,alphk,alph;
double test=1,lapgam(),ilapgam(),Dtlapgam(),Dalphalapgam(),Dilapgam();
int prv,k; 
vector *Dphi,*Dphk;

if (test<1) {
Rprintf("ckr \n"); 
print_vec(dckij); print_vec(rvk); print_vec(rvi); print_vec(alphai); print_vec(alphak); 
}

alphi=vec_prod(rvi,alphai); alphk=vec_prod(rvk,alphak); 

if (test<1) Rprintf("=============================ckr %lf %lf \n",alphi,alphk); 

prv=length_vector(rvi); 
malloc_vec(prv,Dphi); malloc_vec(prv,Dphk); 

val=1; 
for (k=0;k<prv;k++) if (VE(rvi,k)+VE(rvk,k)>0) 
{
val1=VE(rvi,k)*ilapgam(alphi,beta,exp(-x))+
     VE(rvk,k)*ilapgam(alphk,beta,exp(-y)); 
if (VE(rvi,k)>0) alph=VE(alphai,k); else alph=VE(alphak,k); 
val1=lapgam(alph,beta,val1); 
val=val*val1; 
}
ckij[0]=1-exp(-x)-exp(-y)+val; 

if (test<1) Rprintf(" %lf ckij \n",ckij[0]); 

val1=0;
for (k=0;k<prv;k++) if (VE(rvi,k)+VE(rvk,k)>0) 
{
if (VE(rvi,k)>0) alph=VE(alphai,k); else alph=VE(alphak,k); 

val2=VE(rvi,k)*ilapgam(alphi,beta,exp(-x))+VE(rvk,k)*ilapgam(alphk,beta,exp(-y)); 
val1= Dtlapgam(alph,beta,val2);
val3= lapgam(alph,beta,val2);

VE(dckij,k)=VE(dckij,k)+Dalphalapgam(alph,beta,val2)/val3;

scl_vec_mult(val1*VE(rvi,k)*Dilapgam(alphi,beta,exp(-x))/val3,rvi,Dphi); 
scl_vec_mult(val1*VE(rvk,k)*Dilapgam(alphk,beta,exp(-y))/val3,rvk,Dphk); 

vec_add(Dphi,dckij,dckij); 
vec_add(Dphk,dckij,dckij); 
};
scl_vec_mult(val,dckij,dckij); 

//val2=rvi[k]*ilapgam(alphi,beta,exp(-x))+rvk[k]*ilapgam(alphk,beta,exp(-y)); 
//val1= Dtlapgam(alph,beta,val2);
//val3= lapgam(alph,beta,val2);
//
//dckij[k]=dckij[k]+Dalphalapgam(alph,beta,val2)/val3;
//
//Dphi=val1*rvi*rvi[k]*Dilapgam(alphi,beta,exp(-x))/val3; 
//Dphk=val1*rvk*rvk[k]*Dilapgam(alphk,beta,exp(-y))/val3; 
//
//dckij=(dckij+Dphi+Dphk);
//dckij=val*dckij; 

if (test<1) print_vec(dckij); 
free_vecs(&Dphi,&Dphk,NULL); 
if (test<1) Rprintf("=============================================== ude af cvrks \n");
} // }}}

void funkdes2(vector *alphai,vector *alphak, // {{{
		double beta, double x,double y,double *ckij,
		vector *dckij,vector *rvi,vector *rvk)
{
double val,val1,alphi,alphk,alph,betai,betak;
double test=1,lapgam(),ilapgam(),Dtlapgam(),
       Dalphalapgam(),Dilapgam(),Dbetalapgam(),Dbetailapgam();
int prv,k; 
vector *Dphi,*Dphk;

if (test<1) {
Rprintf("ckr \n"); 
print_vec(dckij); print_vec(rvk); print_vec(rvi); print_vec(alphai); print_vec(alphak); 
}

alphi=vec_prod(rvi,alphai); alphk=vec_prod(rvk,alphak); 
betai=alphi; betak=alphk;

if (test<1) Rprintf("=============================ckr %lf %lf \n",alphi,alphk); 

prv=length_vector(rvi); 
malloc_vec(prv,Dphi); malloc_vec(prv,Dphk); 

val=1; 
for (k=0;k<prv;k++) if (VE(rvi,k)+VE(rvk,k)>0) 
{
val1=VE(rvi,k)*ilapgam(alphi,betai,exp(-x))+
     VE(rvk,k)*ilapgam(alphk,betak,exp(-y)); 
if (VE(rvi,k)>0) alph=VE(alphai,k); else alph=VE(alphak,k); 
val1=lapgam(alph,betai,val1); 
val=val*val1; 
}
ckij[0]=1-exp(-x)-exp(-y)+val; 
free_vec(Dphi); free_vec(Dphk); 
} // }}}

void ckrvdes2(vector *alphai,vector *alphak, // {{{
		double beta, double x,double y,double *ckij,
		vector *dckij,vector *rvi,vector *rvk)
{
double val,val1,val2,val3,alphi,alphk,alph,betai,betak;
double *d1=calloc(1,sizeof(double)),*d2=calloc(1,sizeof(double));  
double test=1,lapgam(),ilapgam(),Dtlapgam(),
       Dalphalapgam(),Dilapgam(),Dbetalapgam(),Dbetailapgam();
int prv,k; 
vector *Dphi,*Dphk;
void funkdes2(); 

if (test<1) {
Rprintf("ckr \n"); 
print_vec(dckij); print_vec(rvk); print_vec(rvi); print_vec(alphai); print_vec(alphak); 
}

alphi=vec_prod(rvi,alphai); alphk=vec_prod(rvk,alphak); 
betai=alphi; betak=alphk;

prv=length_vector(rvi); 
malloc_vec(prv,Dphi); malloc_vec(prv,Dphk); 

val=1; 
for (k=0;k<prv;k++) if (VE(rvi,k)+VE(rvk,k)>0) 
{
val1=VE(rvi,k)*ilapgam(alphi,betai,exp(-x))+
     VE(rvk,k)*ilapgam(alphk,betak,exp(-y)); 
if (VE(rvi,k)>0) alph=VE(alphai,k); else alph=VE(alphak,k); 
val1=lapgam(alph,betai,val1); 
val=val*val1; 
}
ckij[0]=1-exp(-x)-exp(-y)+val; 

for (k=0;k<prv;k++) if (VE(rvi,k)+VE(rvk,k)>0) 
{
if (VE(rvi,k)>0) alph=VE(alphai,k); else alph=VE(alphak,k); 
val2=VE(rvi,k)*ilapgam(alphi,betai,exp(-x))+
     VE(rvk,k)*ilapgam(alphk,betak,exp(-y)); 
val1= Dtlapgam(alph,betai,val2);
val3= lapgam(alph,betai,val2);

VE(dckij,k)=VE(dckij,k)+Dalphalapgam(alph,betai,val2)/val3;
scl_vec_mult(val1*VE(rvi,k)*(Dilapgam(alphi,betai,exp(-x))+
	Dbetailapgam(alphi,betai,exp(-x)))/val3,rvi,Dphi); 
scl_vec_mult(val1*VE(rvk,k)*(Dilapgam(alphk,betak,exp(-y))+
        Dbetailapgam(alphk,betai,exp(-y)))/val3,rvk,Dphk); 
vec_add(Dphi,dckij,dckij); vec_add(Dphk,dckij,dckij); 

scl_vec_mult(Dbetalapgam(alph,betai,val2)/val3,rvi,Dphi);
vec_add(Dphi,dckij,dckij); 
};
scl_vec_mult(val,dckij,dckij); 

free_vecs(&Dphi,&Dphk,NULL); 
free(d1); free(d2); 
} // }}}

// {{{ laplace and derivatives
double lapgam(alpha,beta,t)
double alpha,beta,t;
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
return(val); 
}

double ilapgam(alpha,beta,y)
double alpha,beta,y;
{
double val; 
val=beta*(exp(-log(y)/alpha)-1); 
return(val); 
}

double Dilapgam(alpha,beta,y)
double alpha,beta,y;
{
double val; 
val=beta*exp(-log(y)/alpha)*(log(y)/(alpha*alpha)); 
return(val); 
}

double Dbetailapgam(alpha,beta,y)
double alpha,beta,y;
{
double val; 
val=(exp(-log(y)/alpha)-1);
return(val); 
}

double Dalphalapgam(alpha,beta,t)
double alpha,beta,t;
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
val=(log(beta)-log(beta+t))*val; 
return(val); 
}

double Dbetalapgam(alpha,beta,t)
double alpha,beta,t;
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
val=alpha*(1/beta-1/(beta+t))*val; 
return(val); 
}


double Dtlapgam(alpha,beta,t)
double alpha,beta,t;
{
double val; 
val=exp(alpha*(log(beta)-log(beta+t))); 
val=-alpha*(1/(beta+t))*val; 
return(val); 
}
// }}}

void DUetagammarv(double t, double x,double y,vector *xi,vector *xk) 
{

vec_add(xi,xk,xi); 
vec_zeros(xi); 
}
