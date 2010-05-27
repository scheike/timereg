#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "multicomp.h"
                 
/* ====================================================== */
void randomcif(times,Ntimes,x,delta,cause,CA1,KMc,z,antpers,px,Nit,score,hess,est,
gamma,semi,zsem,pg,detail,biid,gamiid,timepow,theta,vartheta,thetades,ptheta,antclust,
cluster,clustsize,clusterindex,maxclust,step,inverse,dscore)
double *theta,*times,*x,*KMc,*z,*score,*hess,*est,*gamma,*zsem,*vartheta,*biid,*gamiid,*timepow,*thetades,*step;
int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*semi,*pg,*CA1,*detail,*ptheta,
*antclust,*cluster,*clustsize,*clusterindex,*maxclust,*inverse,*dscore;
{
 matrix *ldesignX,*A,*AI,*cdesignX,*ldesignG,*cdesignG;
 matrix *S,*dCGam,*CGam,*ICGam,*VarKorG,*XZ,*ZZ,*ZZI,*XZAI;
 matrix *Vargam,*dVargam,*dM1M2,*M1M2t,*RobVargam;
 matrix *DUeta[*Ntimes],*DUgamma;
 matrix *Biid[*antclust]; 
 matrix *destheta,*d2UItheta,*d2Utheta,*varthetascore,*Sthetaiid[*antclust],*Stheta; 
 vector *gammaiid[*antclust], *W2[*antclust], *W3[*antclust], *W4[*antclust];
 vector *diag,*dB,**VdB,*AIXdN,*AIXlamt,*bhatt,*pbhat,*lamtt;
 vector *korG,*pghat0,*pghat,*gam;
 vector *xk,*xi,*rowX,*rowZ,*difX,*zk,*zi,*z1;
 vector *Utheta,*vthetascore,*vtheta1,*vtheta2,*dtheta,*thetaiid[*antclust]; 

 int pmax;
 int v,nb[1],itt,i,j,k,l,s,c,coef[1],ps[1],n[1],nx[1],retur[1];
 double Li,Lk,ithetak,thetak,response,time,dtime,timem;
 double fabs(),Dinverse,DDinverse;
 double sumscore,sdj,pow(),ckij[1],dckij[1],diff;
 long robust[1],fixedcov; 
 float gasdev(),expdev(),ran1();
 void ck(),DUetagamma(); 

 robust[0]=2; fixedcov=1;
 n[0]=antpers[0]; retur[0]=0; nx[0]=antpers[0]; nb[0]=Ntimes[0]; timem=0;
	
//if (*trans==1) for (j=0;j<*pg;j++) if (timepow[j]!= 1) {timem=1;break;}
//if (*trans==2) for (j=0;j<*pg;j++) if (timepow[j]!= 0) {timem=1;break;}

  for (j=0;j<*Ntimes;j++) { malloc_mat(*px,*ptheta,DUeta[j]); }

  for (j=0;j<*antclust;j++) { malloc_vec(*pg,gammaiid[j]); 
  malloc_mat(*Ntimes,*px,Biid[j]); 
  malloc_vec(*ptheta,W2[j]); malloc_vec(*ptheta,W3[j]); malloc_vec(*ptheta,W4[j]); 
  malloc_vec(*ptheta,thetaiid[j]); malloc_mat(*ptheta,*ptheta,Sthetaiid[j]); }

  malloc_mat(*antpers,*ptheta,destheta); malloc_mat(*ptheta,*ptheta,Stheta); 
  malloc_mats(*ptheta,*ptheta,&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  malloc_mats(*antpers,*px,&ldesignX,&cdesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,&cdesignG,NULL); 
  malloc_mats(*px,*px,&A,&AI,NULL);
  malloc_mats(*pg,*pg,&dVargam,&Vargam,&RobVargam,&ZZ,&VarKorG,&ICGam,&CGam,&dCGam,
		      &S,&ZZI,NULL); 
  malloc_mats(*px,*pg,&XZAI,&XZ,&dM1M2,&M1M2t,NULL);
  malloc_mats(*pg,*ptheta,&DUgamma,NULL);

  malloc_vecs(*ptheta,&Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,NULL);
  malloc_vecs(*px,&xk,&xi,&rowX,&difX,&korG,&diag,&dB,&VdB,&AIXdN,
		  &AIXlamt,&bhatt,NULL);
  malloc_vecs(*pg,&zk,&zi,&rowZ,&z1,&gam,NULL);
  malloc_vecs(*antpers,&pbhat,&pghat0,&pghat,&lamtt,NULL);

  coef[0]=1; ps[0]=*px+1; 
  pmax=max(*px,*pg); pmax=max(pmax,*ptheta); 

  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 
  for (j=0;j<*ptheta;j++) VE(vtheta1,j)=theta[j]; 

  if (*semi==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++)
        VE(gammaiid[i],j)=gamiid[j*(*antclust)+i]; 
  if (*semi==2) for (i=0;i<*antclust;i++) print_vec(gammaiid[i]); 

   for (i=0;i<*antclust;i++) {
   for (s=0;s<*Ntimes;s++) 
   for (c=0;c<*px;c++) {l=i*(*px)+c; ME(Biid[i],s,c)=biid[l*(*Ntimes)+s]; }
  // printf(" %ld \n",i); print_mat(Biid[i]); 
   }

    for (c=0;c<*antpers;c++) {
      for(j=0;j<pmax;j++)  {
	if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c];
	if (j<*ptheta) ME(destheta,c,j)= thetades[j*(*antpers)+c];
        if (*semi==1) if (j<*pg) { ME(ldesignG,c,j)=zsem[j*(*antpers)+c];}}} 

// print_mat(ldesignX); print_mat(ldesignG); print_mat(destheta); 

if (*semi==1) Mv(ldesignG,gam,pghat0);

  for (itt=0;itt<*Nit;itt++)
    {
      mat_zeros(d2Utheta); vec_zeros(Utheta); Mv(destheta,vtheta1,lamtt);
      sumscore=0; 

      Dinverse=1; DDinverse=1; 
      for (s=0;s<*Ntimes;s++)
      {
     // printf("times  s %d %d %d \n",s,*Ntimes,*antclust); 

	  time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 
	  for(j=1;j<=*px;j++) {VE(bhatt,j-1)=est[j*(*Ntimes)+s];}
	  Mv(ldesignX,bhatt,pbhat); 
	  if (*semi==1) {scl_vec_mult(time,pghat0,pghat);vec_add(pbhat,pghat,pbhat);}

    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {
          vec_zeros(vtheta2);diff=0;sdj=0;
	  vec_zeros(rowX);vec_zeros(rowZ); 

          for (c=0;c<clustsize[j];c++) for (v=0;v<clustsize[j];v++) 
	  if (v!=c) { 
	    i=clusterindex[c*(*antclust)+j]; k=clusterindex[v*(*antclust)+j];
            response=(((x[i]<=time) & (cause[i]==*CA1))*1*
	              ((x[k]<=time) & (cause[k]==*CA1))*1); 

	    if (KMc[i]<0.001) response=response/0.001; else 
		    response=response/( KMc[i] );
	    if (KMc[k]<0.001) response=response/0.001; else 
	            response=response/( KMc[k] );

       thetak=VE(lamtt,i); Li=VE(pbhat,i); Lk=VE(pbhat,k); 
       ithetak=thetak;  

       ck(ithetak,Li,Lk,ckij,dckij); 
       if (*dscore==1) response=dckij[0]*Dinverse*(response-ckij[0]); 
       else  response=Dinverse*(response-ckij[0]); 
       diff=diff+response; 
       if (*dscore==1) sdj=sdj+DDinverse*dckij[0]*dckij[0]; 
       else  sdj=sdj+DDinverse*dckij[0]; 


        if (isnan(response))   {
	         printf(" %d %d %d \n",clustsize[j],i,k); 
                 printf(" %lf %lf %lf \n",x[i],x[k],time); 
		 printf(" resp %lf  %lf %lf  \n",response,KMc[i],KMc[k]); 
                 printf(" %d %d  \n",i,k); 
                 printf(" %lf %lf %lf \n",thetak,Li,Lk); 
                 printf(" %lf %lf \n",dckij[0],ckij[0]); 
		 printf("============================== \n"); 
	 }

	    
       	if (itt==*Nit-1) {
	   extract_row(ldesignX,i,xi);extract_row(ldesignX,k,xk); 
           DUetagamma(ithetak,Li,Lk,xi,xk); 
	   // printf(" ===================== \n"); print_vec(xi); 
	   if (*dscore==1) scl_vec_mult(dckij[0],xi,xi); 
	   vec_add(xi,rowX,rowX); 
           if (*semi==1)  
	   {
               extract_row(ldesignG,i,zi);extract_row(ldesignG,k,zk); 
               DUetagamma(ithetak,Li,Lk,zi,zk); 
	       if (*dscore==1) scl_vec_mult(dckij[0],zi,zi); 
	       scl_vec_mult(time,zi,zi); vec_add(zi,rowZ,rowZ); 
	   } 
	}
        } /* for (c=0....... */  

	 extract_row(destheta,clusterindex[j],vthetascore); 

   if (itt==*Nit-1) {
      for (k=0;k<*px;k++) for (c=0;c<*ptheta;c++) 
         ME(DUeta[s],k,c)= ME(DUeta[s],k,c)+VE(rowX,k)*VE(vthetascore,c);  
     if (*semi==1) for (k=0;k<*pg;k++) for (c=0;c<*ptheta;c++) 
	 ME(DUgamma,k,c)= ME(DUgamma,k,c)+VE(rowZ,k)*VE(vthetascore,c);  }


         for (k=0;k<*ptheta;k++) 
	 for (c=0;c<*ptheta;c++) {
	 if (itt==*Nit-1) ME(Sthetaiid[j],k,c)=ME(Sthetaiid[j],k,c)+
	               sdj*VE(vthetascore,k)*VE(vthetascore,c);  
	 ME(d2Utheta,k,c)= ME(d2Utheta,k,c)+ 
	               sdj*VE(vthetascore,k)*VE(vthetascore,c);}

         scl_vec_mult(diff,vthetascore,vthetascore);  
         vec_add(vthetascore,Utheta,Utheta); 

	if (itt==*Nit-1) {vec_add(vthetascore,W2[j],W2[j]);}

        } /* j in antclust */ 

   if (itt==(Nit[0]-1)) for (j=0;j<*antclust;j++) {
       extract_row(Biid[j],s,rowX); 
       //printf(" %ld %ld  \n",s,j); print_vec(rowX); 
       vM(DUeta[s],rowX,dtheta); vec_add(dtheta,W3[j],W3[j]);}

      } /* s=1,...Ntimes */

//   printf(" %ld \n",*inverse); 
// print_mat(d2Utheta); // print_vec(Utheta); 

  invert(d2Utheta,d2UItheta); Mv(d2UItheta,Utheta,dtheta);
  scl_vec_mult(step[0],dtheta,dtheta); 

 if (*detail==1) {
    printf("===============Iteration %d ==================== \n",itt);
     printf("Estimate theta \n"); print_vec(vtheta1);
     printf("Score D l\n"); print_vec(Utheta);
     printf("Information D^2 l\n"); print_mat(d2UItheta); }

    for (k=0;k<*ptheta;k++) sumscore= sumscore+fabs(VE(Utheta,k));

    if ((sumscore<0.000001) & (itt<*Nit-2)) itt=*Nit-2;
    if (isnan(vec_sum(Utheta)))  itt=*Nit-1;

    vec_add(vtheta1,dtheta,vtheta1); // adding dtheta

} /*itt løkke */ 


//print_mat(DUgamma); 
//printf("===============lomse lomse ==================== \n",itt);

vec_zeros(dtheta); 
   for (j=0;j<*antclust;j++) 
   {
  // printf("W2  ========= \n"); print_vec(W2[j]); 
  // printf("W3  ========= \n"); print_vec(W3[j]); 
   // vec_star(W3[j],dtheta); 
    vec_subtr(W2[j],W3[j],W2[j]); 
if (*semi==1) { //printf(" =W4======== \n"); print_vec(W4[j]); 
               vM(DUgamma,gammaiid[j],W4[j]); 
               vec_subtr(W2[j],W4[j],W2[j]);
              }
    
    for (k=0;k<*ptheta;k++) 
    for (c=0;c<*ptheta;c++) 
    ME(varthetascore,k,c)=ME(varthetascore,k,c)+VE(W2[j],c)*VE(W2[j],k); 
   }


   MxA(varthetascore,d2UItheta,d2Utheta); 
   MxA(d2UItheta,d2Utheta,varthetascore);
   // print_mat(varthetascore);

  for (j=0;j<*ptheta;j++) {theta[j]=VE(vtheta1,j); score[j]=VE(Utheta,j);
    for (k=0;k<*ptheta;k++) {vartheta[k*(*ptheta)+j]=ME(varthetascore,j,k);
                             hess[k*(*ptheta)+j]=ME(d2UItheta,j,k);}}


// {{{ freeing 
			     
  for (j=0;j<*Ntimes;j++) { free_mat(DUeta[j]); }

  for (j=0;j<*antclust;j++) { free_vec(gammaiid[j]); 
  free_mat(Biid[j]); free_vec(W2[j]); free_vec(W3[j]); free_vec(W4[j]); 
  free_vec(thetaiid[j]); free_mat(Sthetaiid[j]); }

free_mats(&destheta,&Stheta,
  &varthetascore,&d2Utheta,&d2UItheta,
  &ldesignX,&cdesignX, &ldesignG,&cdesignG, &A,&AI,
  &dVargam,&Vargam,&RobVargam,&ZZ,&VarKorG,&ICGam,&CGam,&dCGam,&S,&ZZI,
  &XZAI,&XZ,&dM1M2,&M1M2t,&DUgamma,NULL);

free_vecs(&xk,&xi,&rowX,&difX,&korG,&diag,&dB,&VdB,&AIXdN,&AIXlamt,&bhatt,
  &zk,&zi,&rowZ,&z1,&gam,
  &pbhat,&pghat0,&pghat,&lamtt,
  &Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,NULL);

// }}}


}

double mypow(double x,double p)
{
double val,log(),exp();

val=exp(log(x)*p);
/* printf(" i mypow %lf %lf %lf  \n",x,p,val); */
return(val);
}

void ck(double t,double x,double y,double *ckij,double *dckij)
{
double val,val2,val3,val4,log(),exp(); 
double laplace(),ilaplace(),Dilaplace(),Dlaplace(),D2laplace(); 
double t0;

if (x<0) x=0.0001; if (y<0) y=0.0001; 

val=ilaplace(t,exp(-x))+ilaplace(t,exp(-y)); val2=laplace(t,val); 
ckij[0]=1-exp(-x)-exp(-y)+val2; 

val3=exp(x*t)+exp(y*t)-1; 
val4=val3*log(val3)+exp(x*t)*(-x*t)+exp(y*t)*(-y*t); 

// t0 =exp(-log(t)*2)*exp(log(val3)*(-1/t-1))*val4; 
t0 =pow(1/t,2)*exp(log(val3)*(-1/t-1))*val4; 

dckij[0]=t0; 
}

double laplace(t,x)
double t,x;
{
double val,val1; 
val=(1+x*t); 
if (val<0) oops(" val < 0 \n"); 

//if (fabs(t)< 0.000000000000001) val1=0; else val1=exp(-log(val)*(1/t)); 
val1=exp(-log(val)*(1/t)); 

// printf("laplace %lf %lf  \n",val,val3); 

return(val1); 
}

double ilaplace(t,y)
double t,y;
{
double val,laplace(); 

val=exp(-log(y)*t); val= (val-1)/t;  
// printf("ilaplace y^(1/t)  %lf %lf  \n",exp(log(y)/t),pow(y,1.0/t)); 
// printf("ilaplace  %lf %lf  \n",val,val1); 
return(val); 
}

double Dilaplace(theta,y)
double theta,y;
{
double val4,val2,val,val1,log(),exp(); 
val=exp(log(y)/theta); 
val2=-log(y)*val/(theta*theta); 
val4=(1-val)+theta*val2; 
val1=(val4+log(y)*(1-val)/theta)/val;  
return(val1); 
}

double Dlaplace(double theta,double t)
{
double val,val1,log(),exp(),laplace(); 

val=1+t/theta; val1=theta*val-log(val); 
val=val1*laplace(theta,t); 
return(val); 
}

double D2laplace(double theta,double t)
{
double val,val1,val2,val3,log(),exp(),laplace(); 

val=1+t/theta; val1=theta*val-log(val); 
val3=(t/(theta*theta))/val+(val+t/theta)/(val*val); 
val2=Dlaplace(theta,t)*val1+laplace(theta,t)*val3; 
return(val2); 
}

void DUetagamma(double t, double x,double y,vector *xi,vector *xk) 
{
double y1,y2,t1,val3,val4;

y1=exp(-x); y2=exp(-y); 
val3=exp(x*t)+exp(y*t)-1; 

val4 =exp(log(val3)*(-1/t)); 

t1=val4/(val3); 
if (isnan(t1)) {
printf(" missing values in DUetagamma \n"); 
printf(" t x y val3=exp(x*t)+exp(y*t)-1 %lf %lf %lf %lf  \n",t,x,y,val3); 
print_vec(xi); 
print_vec(xk); 
}; 

scl_vec_mult(y1-t1*exp(t*x),xi,xi); 
scl_vec_mult(y2-t1*exp(t*y),xk,xk); 

vec_add(xi,xk,xi); 
}
