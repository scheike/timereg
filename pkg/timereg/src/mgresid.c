#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include"R_ext/Random.h"

void mgresid(designX,nx,p,antpers,start,stop,status,id,
mgtimes,nmgt,mgresid,sim,rani,xval,ant,
univarproc,timeproc,simunivarproc,simtimeproc,
unitest,unitestOBS,
timetest,timetestOBS,
unitimetest,unitimetestOBS,
modelmatrix,model,pm,cummgt,mgresidiid,robvarcum,
testOBS,test,simUt,Ut,cumresid,maxval,startdesign,
coxaalen,dcum,beta,designG,pg,Ogammaiid,
clusters,antclust,robvarcumz,simcumz,
inXZ,inXorZ,iptot) 
double *designG,*dcum,*beta,*designX,*start,*stop,*mgtimes,
*mgresid,*xval,*univarproc,*timeproc,*simunivarproc,
*simtimeproc,*unitest,*unitestOBS, *timetest,*timetestOBS,
*unitimetest,*unitimetestOBS,*modelmatrix,*Ogammaiid,
*cummgt,*mgresidiid,*robvarcum,*testOBS,*test,*simUt,*Ut,
	*robvarcumz,*simcumz;
int *pg,*coxaalen,*nx,*p,*antpers,*nmgt,*sim,*rani,*ant,
*status,*id,*model,*pm,*cumresid,*maxval,*startdesign,*clusters,*antclust,
*inXZ,*inXorZ,*iptot;
{
  matrix *Delta,*tmpM1,*ldesignX,*cummat,*modelMGT[*antpers],
	 *modMGz,*modMGzosdt;
  matrix *cdesX,*ldesignG,*ZP,*A,*AI,*cumX,*cumXAI,*cumZP,*XPZ,*tmp2,*dS,*St[*nmgt]; 
  matrix *dS1,*cumX1,*cumXAI1,*cumZP1,*tmp21,*cummat1;
  vector *Deltazsd,*Deltaz,*tmpM1z;
  vector *vtmp1,*cumdB1,*VdB1,*respm1;
  vector *dMGt[*antpers],*dMGtiid[*antpers],
    *cumdB,*diag,*dB,*VdB,*xi,*rowX,*rowcum,*difX,*vtmp,*respm,*gamma;
  vector *risk,*cumA[*antclust],*cum,*vecX,*MGt,*MGtiid;
  vector *lamt,*Gbeta,*dA,*xtilde,*zi,*gammaiid[*antpers];
  vector *rvec,*covlesszi,*dB1[*antclust];
  int m,i,j,k,l,s,c,count,pers;
  int ptot,weighted,*cluster=calloc(*antpers,sizeof(int));
  double time,dummy;
  double vardiv;
  double random,fabs(),sqrt(),xij,dtime; 
  void smoothB(),comptest(); 
  long idum; idum=*rani; 
  double norm_rand();
  void GetRNGstate(),PutRNGstate();

  weighted=0; 
  ptot=*p+*pg; 
  ptot=*iptot; 

  GetRNGstate();  /* to use R random normals */
  for (j=0;j<*antpers;j++) cluster[j]=0;

  for (i=0;i<*nmgt;i++){malloc_vec(*antpers,dMGt[i]); 
    malloc_vec(*antpers,dMGtiid[i]); malloc_mat(*pm,*pg,St[i]);}

  for (i=0;i<*antclust;i++) {malloc_mat(*nmgt,*pm,modelMGT[i]);
    malloc_vec(*pg,gammaiid[i]); malloc_vec(*pm,cumA[i]); }
  malloc_mat(*nmgt,*pm,Delta); malloc_mat(*nmgt,*pm,tmpM1); 
  malloc_vec(*p,cum); malloc_mat(*antpers,*p,ldesignX); malloc_vec(*antpers,vecX); 

  malloc_mat(*antpers,*p,cdesX); malloc_mat(*antpers,*pg,ldesignG); 
  malloc_mat(*antpers,*pg,ZP); malloc_mat(*p,*p,A); malloc_mat(*p,*p,AI); 
  malloc_mat(*pm,*p,cumX); malloc_mat(*pm,*p,cumXAI); malloc_mat(*pm,*pg,cumZP); 
  malloc_mat(*p,*pg,XPZ); malloc_mat(*pm,*pg,tmp2); malloc_mat(*pm,*pg,dS); malloc_vec(*p,dA); 
  malloc_vec(*pg,zi); malloc_vec(*pg,gamma); malloc_vec(*antpers,Gbeta); malloc_vec(*antpers,lamt); 
  malloc_vec(*antpers,covlesszi); 
  malloc_vec(*antclust,rvec); 


  malloc_mat(*antpers,*pm,cummat); 
  malloc_vecs(*pm,&vtmp,&cumdB,&dB,&VdB,&respm,NULL);
  malloc_vec(*antpers,risk); 
  malloc_vecs(*p,&diag,&xtilde,&xi,&rowX,&rowcum,&difX,NULL);
  malloc_vec(*antpers,MGt); 
  malloc_vec(*antpers,MGtiid); 

  for(j=0;j<*pg;j++) VE(gamma,j)=beta[j]; 

  if (*coxaalen==1) {
    for (i=0;i<*antclust;i++) { 
      for (j=0;j<*pg;j++) 
	VE(gammaiid[i],j)=Ogammaiid[i*(*pg)+j]; }
  }

  /*  cumulative martingales Aalen type */ 
  if (*model==1) { // {{{
      for (s=1;s<*nmgt;s++)
	{
	  time=mgtimes[s]; dtime=mgtimes[s]-mgtimes[s-1]; 

	 /* if ((*startdesign==1 & s==1) | (*startdesign==0)) {}*/

	 mat_zeros(ldesignX);mat_zeros(cummat);vec_zeros(risk);mat_zeros(ldesignG);
	  for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) // {{{
	    {
	      if ((start[c]<time) && (stop[c]>=time)) {
		cluster[id[c]]=clusters[c];
		if (*coxaalen==1)
		  for(j=0;j<*pg;j++) ME(ldesignG,id[c],j)=designG[j*(*nx)+c];
		for(j=0;j<*p;j++) ME(ldesignX,id[c],j)=designX[j*(*nx)+c];
		for(j=0;j<*pm;j++) ME(cummat,id[c],j)=modelmatrix[j*(*nx)+c];
		VE(risk,id[c])=1;   
		if (time==stop[c] && status[c]==1) {pers=id[c];}
		count=count+1; } 
	    } // }}}

	  for (i=0;i<*antpers;i++) {
	    VE(dMGt[s],i)=mgresid[i*(*nmgt)+s];
	    if (*coxaalen==-1) 
	      VE(dMGtiid[s],i)=mgresidiid[s*(*nmgt)+i];}

	  cummgt[s]=time; robvarcum[s]=time; 

	  vM(cummat,dMGt[s],respm); 
	  if (s<=-2) { 
	    print_mat(cummat); printf(" %lf \n",vec_sum(dMGt[s])); 
	    print_vec(dMGt[s]); }
	  if (s==-2) { print_vec(respm); }

	  for (k=1;k<=*pm;k++) cummgt[k*(*nmgt)+s]=cummgt[k*(*nmgt)+s-1]+VE(respm,k-1);

	  if (*coxaalen==1) Mv(ldesignG,gamma,Gbeta);
	  for (j=0;j<*antpers;j++)
	    {extract_row(ldesignX,j,xi); dummy=exp(VE(Gbeta,j));
	      scl_vec_mult(dummy,xi,xtilde); replace_row(cdesX,j,xtilde); } 

	  MtA(cdesX,ldesignX,A); invert(A,AI);
	  if (ME(AI,0,0)==0)
	    printf(" X'X not invertible at time %lf \n",time);

	  MtA(cummat,cdesX,cumX); MxA(cumX,AI,cumXAI);

	  /* extra terms for cox aalen iid representation */ 
	  if (*coxaalen==1) {
	    for (j=0;j<*p;j++) VE(dA,j)=dcum[j*(*nmgt-1)+s-1]; 

	    Mv(cdesX,dA,lamt);
	    for (j=0;j<*antpers;j++)
	      {extract_row(ldesignG,j,zi); scl_vec_mult(VE(lamt,j),zi,zi);
		replace_row(ZP,j,zi);}

	    MtA(cummat,ZP,cumZP); 
	    MtA(ldesignX,ZP,XPZ); 
	    MxA(cumXAI,XPZ,tmp2);
	    mat_subtr(cumZP,tmp2,dS); /* mat_add(dS,St[s-1],St[s]);  */
	  } /*coxaalen=1 */ 


	  vec_zeros(VdB); 
	  for (i=0;i<*antpers;i++) 
	  {
	   j=cluster[i];
	   extract_row(cummat,i,respm); 
	   extract_row(ldesignX,i,xi); Mv(cumXAI,xi,vtmp);
	   vec_subtr(respm,vtmp,respm); 
	   scl_vec_mult(VE(dMGt[s],i),respm,vtmp);
	   vec_add(vtmp,cumA[j],cumA[j]); 
	  }

          for (j=0;j<*antclust;j++)  {
          if (*coxaalen==1)  {
             Mv(dS,gammaiid[j],respm);vec_subtr(cumA[j],respm,cumA[j]);}

	     replace_row(modelMGT[j],s,cumA[j]); 
	     vec_star(cumA[j],cumA[j],vtmp); vec_add(vtmp,VdB,VdB); 
	  }
	  for (k=1;k<*pm+1;k++) robvarcum[k*(*nmgt)+s]=VE(VdB,k-1); 
 
	  /* comp observed sup statistics */ 
	  Ut[s]=time; 
	  for (i=1;i<=*pm;i++) {
	    if (weighted==1) vardiv=sqrt(robvarcum[i*(*nmgt)+s]); else vardiv=1; 
	    xij=cummgt[i*(*nmgt)+s]/vardiv;
	    if (fabs(xij)>testOBS[i-1]) testOBS[i-1]=fabs(xij);
	    Ut[i*(*nmgt)+s]=xij; 
	    c=(*pm)+i-1; 
	    testOBS[c]=testOBS[c]+xij*xij*dtime; 
	  }
	}


      /* simulation of processes under the model */ 
      for (k=0;k<*sim;k++) {
	mat_zeros(Delta); 
	for (i=0;i<*antclust;i++) { 
	  /*  random=gasdev(&idum);  */ 
	  random=norm_rand();
	  scl_mat_mult(random,modelMGT[i],tmpM1); mat_add(tmpM1,Delta,Delta); }

	for (s=1;s<*nmgt;s++) { dtime=mgtimes[s]-mgtimes[s-1]; 

	  for (i=1;i<=*pm;i++) {
	    if (weighted==1) vardiv=sqrt(robvarcum[i*(*nmgt)+s]); else vardiv=1; 
	    xij=ME(Delta,s,i-1)/vardiv; 
	    if (fabs(xij)>test[(*sim)*(i-1)+k]) test[(*sim)*(i-1)+k]=fabs(xij);

	    if (k<50) {l=k*(*pm)+i-1; simUt[l*(*nmgt)+s]=xij;}
      
	    c=*pm+i-1; 
	    test[(*sim)*c+k]=test[(*sim)*c+k]+xij*xij*dtime; 
	    c=2*(*pm)+i-1; 
	    xij=xij/sqrt(robvarcum[i*(*nmgt)+s]); 
	    if (fabs(xij)>test[(*sim)*c+k]) test[(*sim)*c+k]=fabs(xij);
	  }
	}
      }
    } // }}}

  /* =================================================   */


  /* LWY cumulative residuals versus covariates */ 
  if (*cumresid>0) { // {{{

  for (l=0;l<ptot;l++)  {
     // printf(" %d %d ================ \n",ptot,ant[l]); 

  if (ant[l]>2) {
 // {{{ allokering
 malloc_mat(ant[l],*pg,dS1);malloc_mat(ant[l],*p,cumX1);
 malloc_mat(ant[l],*p,cumXAI1); 
 malloc_mat(ant[l],*pg,cumZP1);
 malloc_mat(ant[l],*pg,tmp21);malloc_mat(*antpers,ant[l],cummat1); 
 malloc_vecs(ant[l],&vtmp1,&cumdB1,&VdB1,&respm1,NULL);
 for (k=0;k<*antclust;k++) malloc_vec(ant[l],dB1[k]); 
 malloc_mat(*antclust,ant[l],modMGz); 
 malloc_mat(*antclust,ant[l],modMGzosdt); 
    malloc_vec(ant[l],Deltaz); malloc_vec(ant[l],Deltazsd);
    malloc_vec(ant[l],tmpM1z);  
    for (i=0;i<*antclust;i++) vec_zeros(dB1[i]);
    // }}}

    for (s=1;s<*nmgt;s++)
    {
	time=mgtimes[s]; dtime=mgtimes[s]-mgtimes[s-1]; 
	cummgt[s]=time; robvarcum[s]=time; 
	mat_zeros(ldesignX);mat_zeros(ldesignG); 
	for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) // {{{
	{
	    if ((start[c]<time) && (stop[c]>=time)) {
	      cluster[id[c]]=clusters[c];
	      if (*coxaalen==1)
		for(k=0;k<*pg;k++) ME(ldesignG,id[c],k)=designG[k*(*nx)+c];
	      for(k=0;k<*p;k++) ME(ldesignX,id[c],k)=designX[k*(*nx)+c];
	      if (time==stop[c] && status[c]==1) {pers=id[c];}
	      count=count+1; } 
	} // }}}

 if ((*startdesign==1 && s==1) || (*startdesign==0)) {// {{{making Kz matrix 
   for (j=0;j<ant[l];j++) { 
     for (i=0;i<*antpers;i++) 
        if (inXorZ[l]==0) 
	ME(cummat1,i,j)=(ME(ldesignG,i,inXZ[l])<=xval[(*maxval)*l+j]);
        else ME(cummat1,i,j)=(ME(ldesignX,i,inXZ[l])<=xval[(*maxval)*l+j]); 
   }// }}}
 }

	if (*model==0)  { 
	  for (i=0;i<*antpers;i++) {
	    VE(dMGt[s],i)=mgresid[i*(*nmgt)+s];
	    if (*coxaalen==-1) VE(dMGtiid[s],i)= mgresidiid[s*(*nmgt)+i] ;} }

	if (*coxaalen==1)  Mv(ldesignG,gamma,Gbeta);

	for (k=0;k<*antpers;k++)
	{extract_row(ldesignX,k,xi); dummy=exp(VE(Gbeta,k));
         scl_vec_mult(dummy,xi,xtilde); replace_row(cdesX,k,xtilde); 
	} 

	MtA(cdesX,ldesignX,A); invert(A,AI);
	if (ME(AI,0,0)==0) printf(" X'X not invertible at time %lf \n",time);

	/* observed increment */ 
	vM(cummat1,dMGt[s],respm1); 
	for (j=0;j<ant[l];j++) 
	univarproc[(*maxval)*l+j]=univarproc[(*maxval)*l+j]+VE(respm1,j);

	/* iid representation increment */ 
	MtA(cummat1,cdesX,cumX1); MxA(cumX1,AI,cumXAI1);

        if (s<-2) {
	    //print_mat(cummat1); 
	    //print_mat(cdesX); 
            printf(" %d %d %d \n",s,j,l); 
	    print_mat(cumX1); 
	    print_vec(dMGt[s]); 
	}
	/* extra terms for cox aalen iid representation */ 
	if (*coxaalen==1) { // {{{
	   for (k=0;k<*p;k++) VE(dA,k)=dcum[k*(*nmgt-1)+s-1]; 

	   Mv(cdesX,dA,lamt);
	   for (k=0;k<*antpers;k++)
	   {  
	     extract_row(ldesignG,k,zi); scl_vec_mult(VE(lamt,k),zi,zi); 
	     replace_row(ZP,k,zi);
	    }

	   MtA(cummat1,ZP,cumZP1); 
	   MtA(ldesignX,ZP,XPZ); 
	   MxA(cumXAI1,XPZ,tmp21);
	   mat_subtr(cumZP1,tmp21,dS1);/*mat_add(dS1,St1[s-1],St1[s]);*/
          }  // }}}

      for (i=0;i<*antclust;i++) vec_zeros(dB1[i]); 

       for (i=0;i<*antpers;i++)  // {{{
       {
	m=cluster[i]; 
	extract_row(cummat1,i,respm1); 
	extract_row(ldesignX,i,xi); Mv(cumXAI1,xi,vtmp1);
	vec_subtr(respm1,vtmp1,respm1); 
	scl_vec_mult(VE(dMGt[s],i),respm1,vtmp1); 
	vec_add(dB1[m],vtmp1,dB1[m]); 
       } // }}}

       for (m=0;m<*antclust;m++)  // {{{
       {
          if (*coxaalen==1) { Mv(dS1,gammaiid[m],respm1);
	     vec_subtr(dB1[m],respm1,dB1[m]);}
	  for (j=0;j<ant[l];j++) 
	  ME(modMGz,m,j)=ME(modMGz,m,j)+VE(dB1[m],j); 
      } // }}}

    } /* s=1... *Ntimes */ 


    //for (j=0;j<ant[l];j++)  printf("%d %lf  \n",l,univarproc[(*maxval)*l+j]);

    /* robust variance not computed in this case */
    for (m=0;m<*antclust;m++)  { // {{{
      for (j=0;j<ant[l];j++) VE(VdB1,j)+=pow(ME(modMGz,m,j),2);
      //vec_star(Deltaz,Deltaz,respm1); vec_add(respm1,VdB1,VdB1); 
    }
    for (j=0;j<ant[l];j++) robvarcumz[l*(*maxval)+j]=VE(VdB1,j); 
    for (m=0;m<*antclust;m++)  { 
    for (j=0;j<ant[l];j++) 
    ME(modMGzosdt,m,j)=ME(modMGz,m,j)/sqrt(robvarcumz[l*(*maxval)+j]); 
    }

    for (j=0;j<ant[l];j++) 
    {
      xij=univarproc[(*maxval)*l+j]; 
      if (fabs(xij)>unitimetestOBS[l]) unitimetestOBS[l]=fabs(xij); 
    }

    /* simulation of testprocesses and teststatistics */ 
   //  printf("Simulations start N= %d \n",*sim);

    for (k=0;k<*sim;k++) {
      for (i=0;i<*antclust;i++) VE(rvec,i)=norm_rand(); 
	vM(modMGz,rvec,Deltaz); vM(modMGzosdt,rvec,Deltazsd); 

	for (j=0;j<ant[l];j++) 
	{
	   xij=VE(Deltaz,j); 
	   if (k<50) {c=k*(ptot)+l; simunivarproc[c*(*maxval)+j]=xij;}
	   if (fabs(xij)>unitest[(*sim)*(l)+k]) unitest[(*sim)*(l)+k]=fabs(xij); 
	   xij=VE(Deltazsd,j); 
	   if (fabs(xij)>simcumz[(*sim)*(l)+k]) simcumz[(*sim)*(l)+k]=fabs(xij); 
	} 
    }  /* k=1..antsim */ 


// {{{ free allokering
free_mats(&dS1,&cumX1, &cumXAI1, &cumZP1,&tmp21,&cummat1,
		 &modMGz,&modMGzosdt,NULL); 
free_vecs(&vtmp1,&cumdB1,&VdB1,&respm1,&Deltaz,&Deltazsd,&tmpM1z,NULL);
for (k=0;k<*antclust;k++) free_vec(dB1[k]); 
// }}}
    
    } 
   } /* l=0,...,ptot */ 

  } // }}} cumresid=1


  PutRNGstate();  /* to use R random normals */

  /*
  free_vecs(&MGt,&MGtiid,&risk,&xi,&rowX,&diag,&dB,&VdB,&rowcum,&cum,&vtmp,&lamt,&dA,&Gbeta,&gamma,&vtmp1,&cumdB1,&VdB1,&respm1,NULL);
  free_mats(&Delta,&tmpM1,&ldesignX,&cdesX,&ldesignG,&ZP,&A,&AI,&XPZ,
	      &cumX, &cumXAI, &cumZP, &tmp2, &dS,&cummat,NULL); 
//	      &cumX1,&cumXAI1,&cumZP1,&tmp21,&dS1,&cummat1,NULL); 
  for (i=0;i<*antclust;i++) {free_vec(cumA[i]); free_mat(modelMGT[i]); }
  for (i=0;i<*nmgt;i++) {free_vec(dMGt[i]);free_mat(St[i]);free_vec(dMGtiid[i]);} 
  */
  free(cluster); 
}
