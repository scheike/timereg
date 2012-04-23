#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "multicomp.h"
                 
/* ====================================================== */
void randomcif2cause(times,Ntimes,x, delta,cause,CA1,
		KMc,z,antpers, px,Nit,score,
		hess,est,gamma, semi,zsem,pg,
		detail,biid,gamiid,timepow,theta,vartheta,
		thetades,ptheta,antclust, cluster,clustsize,clusterindex,
		maxclust,step,inverse,CA2,x2,px2,
		semi2,z2,pg2,est2,gamma2,b2iid,
		gam2iid,dscore,squarepar,c1fc2,samecens,cifmodel
)
double *theta,*times,*x,*KMc,*z,*score,*hess,*est,*gamma,*zsem,*vartheta,*biid,*gamiid,*timepow,*thetades,*step,*x2,*z2,*est2,*gamma2,*b2iid,*gam2iid;
int *antpers,*px,*Ntimes,*Nit,*cause,*delta,*semi,*pg,*CA1,*CA2,*detail,*ptheta,
*antclust,*cluster,*clustsize,*clusterindex,*maxclust,*inverse,
	*pg2,*px2,*semi2,*dscore,*squarepar,*c1fc2,*samecens,*cifmodel;
{ // {{{
// {{{ allocating and setting up memory
 matrix *ldesignX,*ldesignG,*X2,*Z2;
 matrix *DUeta[*Ntimes],*DUeta2[*Ntimes],*DUgamma2,*DUgamma;
 matrix *Biid[*antclust],*B2iid[*antclust]; 
 matrix *destheta,*d2UItheta,*d2Utheta,*varthetascore;//*Sthetaiid[*antclust]; 
 vector *gamma2iid[*antclust],*gammaiid[*antclust],
	*W2[*antclust],*W3[*antclust],*W4[*antclust];
 vector *dB,*VdB,*bhatt,*pbhat,*plamt,*lamtt,*pghat0,*pghat,*gam;
 vector *xk,*xi,*rowX,*rowZ,*difX,*zk,*zi;
 vector *Utheta,*vthetascore,*vtheta1,*vtheta2,*dtheta;  
 vector *gam2,*bhatt2,*pbhat2,*pghat2,*pghat02,*rowX2,*xi2,*xk2,*zi2,*zk2,
	*rowZ2; 

 int naprint=0,l2,pmax,v,itt,i,j,k,l,s,c;
 double Li2,Lk2,Li,Lk,ithetak,thetak,response,time,dtime;
 double fabs(),Dinverse,DDinverse,count24=0;
 double sumscore,sdj,pow(), diff, 
        *ckij=calloc(1,sizeof(double)), *dckij=calloc(1,sizeof(double)),
        *ckij2=calloc(1,sizeof(double)), *dckij2=calloc(1,sizeof(double));
 long indc1fc2ik=0,indc1fc2ki=0;
 float gasdev(),expdev(),ran1();
 void ck(),DUetagamma(); 
	
//if (*trans==1) for (j=0;j<*pg;j++) if (timepow[j]!= 1) {timem=1;break;}
//if (*trans==2) for (j=0;j<*pg;j++) if (timepow[j]!= 0) {timem=1;break;}

  for (j=0;j<*Ntimes;j++) { 
    malloc_mat(*px,*ptheta,DUeta[j]); malloc_mat(*px2,*ptheta,DUeta2[j]); 
  }
  for (j=0;j<*antclust;j++) { 
     malloc_vec(*pg,gammaiid[j]); malloc_mat(*Ntimes,*px,Biid[j]); 
     malloc_vec(*pg2,gamma2iid[j]); malloc_mat(*Ntimes,*px2,B2iid[j]); 
     malloc_vec(*ptheta,W2[j]); malloc_vec(*ptheta,W3[j]); 
     malloc_vec(*ptheta,W4[j]); 
//      malloc_mat(*ptheta,*ptheta,Sthetaiid[j]); 
  }
  malloc_mats(*antpers,*px,&ldesignX,NULL);
  malloc_mats(*antpers,*pg,&ldesignG,NULL); 
  malloc_mats(*antpers,*pg2,&Z2,NULL);
  malloc_mats(*antpers,*px2,&X2,NULL);
  malloc_mats(*pg,*ptheta,&DUgamma,NULL);
  malloc_mats(*pg2,*ptheta,&DUgamma2,NULL);
  malloc_mats(*ptheta,*ptheta,&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  malloc_mat(*antpers,*ptheta,destheta); 

  malloc_vecs(*ptheta,&Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,NULL);
  malloc_vecs(*antpers,&pbhat2,&pghat02,&pghat2,NULL);
  malloc_vecs(*px2,&bhatt2,&rowX2,&xi2,&xk2,NULL);
  malloc_vecs(*pg2,&gam2,&rowZ2,&zi2,&zk2,NULL);
  malloc_vecs(*px,&xk,&xi,&rowX,&difX,&dB,&VdB,&bhatt,NULL);
  malloc_vecs(*pg,&zk,&zi,&rowZ,&gam,NULL);
  malloc_vecs(*antpers,&pbhat,&pghat0,&pghat,&plamt,&lamtt,NULL);

  pmax=max(*px,*pg); pmax=max(pmax,*ptheta); 
  int pmax2=max(*px2,*pg2); 

  for (j=0;j<*pg;j++) VE(gam,j)=gamma[j]; 
  for (j=0;j<*ptheta;j++) VE(vtheta1,j)=theta[j]; 

  if (*semi==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg;j++)
        VE(gammaiid[i],j)=gamiid[j*(*antclust)+i]; 
  if (*semi2==1) for (i=0;i<*antclust;i++) for (j=0;j<*pg2;j++)
        VE(gamma2iid[i],j)=gam2iid[j*(*antclust)+i]; 

   for (i=0;i<*antclust;i++) 
   for (s=0;s<*Ntimes;s++) {
   for (c=0;c<*px;c++) {l=i*(*px)+c; ME(Biid[i],s,c)=biid[l*(*Ntimes)+s]; }
   for (c=0;c<*px2;c++) {l2=i*(*px2)+c; ME(B2iid[i],s,c)=b2iid[l2*(*Ntimes)+s]; }
   }

    for (c=0;c<*antpers;c++) {
      for(j=0;j<pmax;j++)  {
	if (j<*px) ME(ldesignX,c,j)= z[j*(*antpers)+c];
	if (j<*ptheta) ME(destheta,c,j)= thetades[j*(*antpers)+c];
        if (*semi==1) if (j<*pg) { ME(ldesignG,c,j)=zsem[j*(*antpers)+c];}} 
        if (*CA1!=*CA2) {
           for(j=0;j<pmax2;j++)  {
              if (j<*px2) ME(X2,c,j)=x2[j*(*antpers)+c];
              if (*semi2==1) if (j<*pg2) { ME(Z2,c,j)=z2[j*(*antpers)+c];}}
              for (j=0;j<*pg2;j++) VE(gam2,j)=gamma2[j]; 
        } 
    }

    // }}}


if (*semi==1) Mv(ldesignG,gam,pghat0);
if (*CA1!=*CA2 && *semi2==1) Mv(Z2,gam2,pghat02);

  for (itt=0;itt<*Nit;itt++) // {{{
  { 
     R_CheckUserInterrupt();
     sumscore=0; 
     mat_zeros(d2Utheta); vec_zeros(Utheta); Mv(destheta,vtheta1,lamtt);

      for (s=0;s<*Ntimes;s++)
      {
	  time=times[s]; if (s==0) dtime=0; else dtime=time-times[s-1]; 
	  for(j=1;j<=*px;j++) VE(bhatt,j-1)=est[j*(*Ntimes)+s];
	  Mv(ldesignX,bhatt,pbhat); 
	  if ((*semi==1) & (*cifmodel==1)) {scl_vec_mult(time,pghat0,pghat);vec_add(pbhat,pghat,pbhat);}
	  if ((*semi==1) & (*cifmodel==2)) for (c=0;c<*antpers;c++)  VE(pbhat,c)=VE(pbhat,c)*exp(VE(pghat0,c)); 
          if (*CA1!=*CA2) {
	     for(j=1;j<=*px2;j++) {VE(bhatt2,j-1)=est2[j*(*Ntimes)+s];}
	     Mv(X2,bhatt2,pbhat2); 
	  if ((*semi2==1) & (*cifmodel==1)) {scl_vec_mult(time,pghat02,pghat2);vec_add(pbhat2,pghat2,pbhat2);}
	  if ((*semi2==1) & (*cifmodel==2)) for(c=0;c<*antpers;c++)  VE(pbhat2,c)=VE(pbhat2,c)*exp(VE(pghat02,c)); 

	  }

    for (j=0;j<*antclust;j++) if (clustsize[j]>=2) {

          vec_zeros(vtheta2);diff=0;sdj=0;
	  vec_zeros(rowX);vec_zeros(rowZ); 
	  vec_zeros(rowX2);vec_zeros(rowZ2); 

          for (c=0;c<clustsize[j];c++) for (v=0;v<clustsize[j];v++) 
	  if (c!=v) { 
	    i=clusterindex[c*(*antclust)+j]; k=clusterindex[v*(*antclust)+j];
	    if (*c1fc2==1) indc1fc2ik=(x[i]<x[k]); else  indc1fc2ik=1; 
	    if (*c1fc2==1) indc1fc2ki=(x[k]<x[i]); else  indc1fc2ki=1; 
            response=(
            ((x[i]<=time) && (x[k]<=time) && ((cause[i]==*CA1) && (cause[k]==*CA2)))*
	    indc1fc2ik
            + 
            ((x[i]<=time) && (x[k]<=time) && ((cause[k]==*CA1) && (cause[i]==*CA2)))* 
	    indc1fc2ki); 
 
//    Rprintf(" s j %d  %d\n",s,j); 
     if (*samecens==1) response=response/min(KMc[i],KMc[k]); else response=response/(KMc[i]*KMc[k]);

       thetak=VE(lamtt,i); 
       Li=VE(pbhat,i); Lk=VE(pbhat,k); Li2=VE(pbhat2,i); Lk2=VE(pbhat2,k); 

      if (*squarepar==0) {
         ithetak=thetak;  
         Dinverse=1; DDinverse=1; 
      } else { ithetak=pow(thetak,2); Dinverse=2*VE(vtheta1,0); DDinverse=2; }


       ck(ithetak,Li,Lk2,ckij,dckij); 
       ck(ithetak,Li2,Lk,ckij2,dckij2); 

//   Rprintf(" %d %d %d %d %d \n",s,j,clustsize[j],i,k); 
//   Rprintf(" %lf %lf %lf \n",thetak,Li,Lk); 
//   Rprintf(" %lf %lf %lf \n",thetak,Li2,Lk2); 
//  Rprintf(" %lf %lf %lf %lf \n",ckij[0],ckij2[0],dckij[0],dckij2[0]); 


       if (*dscore==1) response=(indc1fc2ik*dckij[0]+indc1fc2ki*dckij2[0])*
	               Dinverse*(response-indc1fc2ik*ckij[0]-indc1fc2ki*ckij2[0]); 
       else  response=Dinverse*(response-indc1fc2ik*ckij[0]-indc1fc2ki*ckij2[0]); 
       diff=diff+response; 
       if (*dscore==1) sdj=sdj+DDinverse*pow((indc1fc2ik*dckij[0]+indc1fc2ki*dckij2[0]),2); 
       else  sdj=sdj+DDinverse*(indc1fc2ik*dckij[0]+indc1fc2ki*dckij2[0]); 

	if (isnan(response))   { // removes these from score equations
	   diff=0; sdj=0; 
	}

	if ((isnan(response)) && (naprint==0))   { // {{{ print diverse na information
	         Rprintf(" %d %d %d \n",clustsize[j],i,k); 
                 Rprintf(" %lf %lf %lf \n",x[i],x[k],time); 
		 Rprintf(" resp, cens1, cens2  %lf  %lf %lf  \n",response,KMc[i],KMc[k]); 
                 Rprintf(" %lf %lf %lf \n",thetak,Li,Lk); 
                 Rprintf(" %lf %lf %lf \n",thetak,Li2,Lk2); 
                 Rprintf(" %lf %lf %lf %lf \n",
				 ckij[0],ckij2[0],dckij[0],dckij2[0]); 
		 Rprintf("============================== \n"); 
		 naprint=1; 
	 } // }}}


        if (itt==*Nit-1) { // {{{
	   if (*CA1!=*CA2 ) {
	      extract_row(ldesignX,i,xi);extract_row(ldesignX,k,xk); 
	      extract_row(X2,i,xi2);extract_row(X2,k,xk2); 
              DUetagamma(ithetak,Li,Lk,xi,xk); 
	      if (*dscore==1) scl_vec_mult(dckij[0],xi,xi); 
	      scl_vec_mult(indc1fc2ik,xi,xi); 
	      vec_add(xi,rowX,rowX); 
              DUetagamma(ithetak,Li2,Lk2,xi2,xk2); 
	      if (*dscore==1) scl_vec_mult(dckij2[0],xi2,xi2); 
	      scl_vec_mult(indc1fc2ki,xi2,xi2); 
	      vec_add(xi2,rowX2,rowX2); 
	   } 

           if (*semi==1)  {
	      if (*CA1!=*CA2 ) {
                  extract_row(ldesignG,i,zi);extract_row(ldesignG,k,zk); 
                  extract_row(Z2,i,zi2);extract_row(Z2,k,zk2); 
                  DUetagamma(ithetak,Li,Lk,zi,zk); 
	          if (*dscore==1) scl_vec_mult(dckij[0],zi,zi); 
	          scl_vec_mult(indc1fc2ik,zi,zi); 
                  DUetagamma(ithetak,Li2,Lk2,zi2,zk2); 
	          if (*dscore==1) scl_vec_mult(dckij2[0],zi2,zi2); 
	          scl_vec_mult(indc1fc2ki,zi2,zi2); 
	          scl_vec_mult(time,zi,zi); vec_add(zi,rowZ,rowZ); 
	          scl_vec_mult(time,zi2,zi2); vec_add(zi2,rowZ2,rowZ2); 
	      } 
	      }
	   } // }}}

        } /* for (c=0....... */  

	 extract_row(destheta,clusterindex[j],vthetascore); 

   if (itt==*Nit-1) {
      if (*CA1!=*CA2) {
         for (k=0;k<*px;k++) for (c=0;c<*ptheta;c++) 
            ME(DUeta[s],k,c)= ME(DUeta[s],k,c)+indc1fc2ik*VE(rowX,k)*VE(vthetascore,c);  
         for (k=0;k<*px2;k++) for (c=0;c<*ptheta;c++) 
            ME(DUeta2[s],k,c)= ME(DUeta2[s],k,c)+indc1fc2ki*VE(rowX2,k)*VE(vthetascore,c);  
        if (*semi==1) for (k=0;k<*pg;k++) for (c=0;c<*ptheta;c++) 
    	 ME(DUgamma,k,c)= ME(DUgamma,k,c)+indc1fc2ik*VE(rowZ,k)*VE(vthetascore,c);  
        if (*semi2==1) for (k=0;k<*pg2;k++) for (c=0;c<*ptheta;c++) 
    	 ME(DUgamma2,k,c)= ME(DUgamma2,k,c)+indc1fc2ki*VE(rowZ2,k)*VE(vthetascore,c);  
      }
      else  {
      }
   }

   for (k=0;k<*ptheta;k++) 
   for (c=0;c<*ptheta;c++) {
//     if (itt==*Nit-1) ME(Sthetaiid[j],k,c)=ME(Sthetaiid[j],k,c)+
//	               sdj*VE(vthetascore,k)*VE(vthetascore,c);  
	 ME(d2Utheta,k,c)= ME(d2Utheta,k,c)+ 
	               sdj*VE(vthetascore,k)*VE(vthetascore,c);}
         scl_vec_mult(diff,vthetascore,vthetascore);  
         vec_add(vthetascore,Utheta,Utheta); 

	if (itt==*Nit-1) {vec_add(vthetascore,W2[j],W2[j]);}

   } /* j in antclust */ 


   if (itt==*Nit-1 && *detail==2) { Rprintf(" s er %d \n",s); print_mat(DUeta[s]);  }

   if (itt==*Nit) for (j=0;j<*antclust;j++) {
       extract_row(Biid[j],s,rowX); 
       vM(DUeta[s],rowX,dtheta); vec_add(dtheta,W3[j],W3[j]);
       extract_row(B2iid[j],s,rowX2); 
       vM(DUeta2[s],rowX2,dtheta); vec_add(dtheta,W3[j],W3[j]);
   }

      } /* s=1,...Ntimes */

  //   Rprintf(" %ld \n",*inverse); 
  // print_mat(d2Utheta); // print_vec(Utheta); 


  invert(d2Utheta,d2UItheta); Mv(d2UItheta,Utheta,dtheta);
  scl_vec_mult(step[0],dtheta,dtheta); 

 if (*detail==1) {
    Rprintf("===============Iteration %d ==================== \n",itt);
    Rprintf(" %lf \n",count24); 
     Rprintf("Estimate theta \n"); print_vec(vtheta1);
     Rprintf("Score D l\n"); print_vec(Utheta);
     Rprintf("Information D^2 l\n"); print_mat(d2UItheta); }

     for (k=0;k<*ptheta;k++) sumscore= sumscore+fabs(VE(Utheta,k));

    if ((sumscore<0.000001) & (itt<*Nit-2)) itt=*Nit-2;
    if (isnan(vec_sum(Utheta)))  itt=*Nit-1;

     vec_add(vtheta1,dtheta,vtheta1);

} // }}} /*itt løkke */ 

   vec_zeros(dtheta); 
   for (j=0;j<*antclust;j++) 
   {
      vec_subtr(W2[j],W3[j],W2[j]); 
      if (*semi==1) { //Rprintf(" =W4======== \n"); print_vec(W4[j]); 
                   vM(DUgamma,gammaiid[j],dtheta); 
		   vec_subtr(W2[j],dtheta,W2[j]);
      }
      if (*semi2==1) { //Rprintf(" =W4======== \n"); print_vec(W4[j]); 
                   vM(DUgamma2,gamma2iid[j],dtheta); 
		   vec_subtr(W2[j],dtheta,W2[j]);
      }
    
       for (k=0;k<*ptheta;k++) 
       for (c=0;c<*ptheta;c++) 
       ME(varthetascore,k,c)=ME(varthetascore,k,c)+VE(W2[j],c)*VE(W2[j],k); 
   }

   MxA(varthetascore,d2UItheta,d2Utheta); 
   MxA(d2UItheta,d2Utheta,varthetascore);

  for (j=0;j<*ptheta;j++) {theta[j]=VE(vtheta1,j); score[j]=VE(Utheta,j);
    for (k=0;k<*ptheta;k++) {vartheta[k*(*ptheta)+j]=ME(varthetascore,j,k);
                             hess[k*(*ptheta)+j]=ME(d2UItheta,j,k);}}

   // {{{ freeing 
  free(ckij); free(dckij); free(ckij2); free(dckij2); 
   
  for (j=0;j<*Ntimes;j++) { free_mat(DUeta[j]); free_mat(DUeta2[j]); }

  for (j=0;j<*antclust;j++) { 
     free_vec(gammaiid[j]); free_mat(Biid[j]); 
     free_vec(gamma2iid[j]); free_mat(B2iid[j]); 
     free_vec(W2[j]); free_vec(W3[j]); free_vec(W4[j]); 
//     free_mat(Sthetaiid[j]); 
  }

  free_mats(&ldesignX,&ldesignG,&Z2,&X2,&DUgamma,&DUgamma2,NULL);
  free_mats(&varthetascore,&d2Utheta,&d2UItheta,NULL); 
  free_mats(&destheta,NULL); 

  free_vecs(&pbhat2,&pghat02,&pghat2,
  &Utheta,&vthetascore,&vtheta1,&dtheta,&vtheta2,
  &bhatt2,&rowX2,&xi2,&xk2,&gam2,&rowZ2,&zi2,&zk2,NULL); 

  free_vecs(&xk,&xi,&rowX,&difX,&dB,&VdB,&bhatt,&zk,&zi,&rowZ,&gam,
            &pbhat,&pghat0,&pghat,&plamt,&lamtt,NULL);
// }}}
  // }}}
}

