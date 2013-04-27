/*gcc -c -munaligned-doubles -I/opt/local/lib/splus-3.3/include -I/pack/meschach/dist/meschach/ -O fil.navn*/

/* ld -r -o allfunctions.o funktioner.o meschach.a -lm */
/* ld -r -o allfunctions.o funktioner.o randomC.o meschach.a -lm */

/* -L/pack/meschach */
/* gcc -c source.c */
/* #include <matrix.h> */
/* #include <matrix2.h> */

/*   TS-help 30/5-2005 
gcc -O -c  addmult.c
gcc -shared -Wl -o  addmult.so addmult.o -L/coll/local/lib -lm meschach.a
*/


//#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include"R_ext/Random.h"

/*#include <S.h>*/


/* #########################################################*/


void addmult(time,status,Xinp,Xtilinp,Zinp,Uinp,dUinp,optinp,
             excess,phi,stid,beta,n,p,q,k,tol,alpha,Psiinp,
             CoVarPsiinp,VarPsiinp,rani,testinp,testinpHW,testinpCM,testinpGOFCM,
             Scoreinp,antsim,k1)
   double *time,*status,*Xinp,*Zinp,*Xtilinp,*Uinp,*dUinp,*optinp,
          *excess,*phi,*stid,*beta,*tol,*alpha,*Psiinp,*CoVarPsiinp,
          *VarPsiinp,*testinp,*testinpHW,*testinpCM,*testinpGOFCM,*Scoreinp;
   int *n,*p,*q,*k,*rani,*k1, *antsim ;
{ int i,j,l,l1,l2,it,init_it,nummer;
  double alpha_tmp,y[*n], y1[*n], y2[*n],tmp1_sc,del,del_old,betaZ[*n],tau,b,beta_tmp[*q],random,testOBSGOFCM,dum1;

  vector *vtmp1,*vtmp2,*vtmp3,*vtmp4,*vtmp5,*vtmp6,*vtmp7;
  vector *testOBS, *testOBSHW, *testOBSCM, *testtmp,*testtmpHW,*testtmpCM,*testtmp1,
    *testtmp2, *testGOFCM,*testtmp1GOFCM,*Ufunkdim1,*Ufunkdim12, *Delta2tmp, 
    *Delta2tmp1;
  matrix *mtmp1,*mtmp2,*mtmp2t,*mtmp3,*mtmp3m,*mtmp3mm,*mtmp3mmm,*mtmp4,*mtmp5,*mtmp6,
    *tmtmp1,*testOBS1,*testHW,*testCM;
  matrix *Xt,*PhiXt,*tmp10_mat,*DdN_Xt,*dM2m;
  matrix *tXt,*tDdN_Xt,*tdM2m;
  matrix *X;
  matrix *Z,*Phi_Z,*DdN_Phi_Z,*tmp3_mat,*tmp4_mat,*tmp5_mat,*tmp6_mat,*ttmpI;
  matrix *dU_dbeta,*dU_dbeta1,*dU_dbeta_tmp,*dU_dbeta_I,*dU_dbeta_I1,*tmp7_mat,
    *dU_dbeta_tilv,*tmp13_mat,
    *opt_tilv,*opt,*dM1,*dM1tmp;
  matrix *tZ,*tZ_Phi,*tDdN_Phi_Z,*tZ_PhiXt_tXt_Xt_I_tDdN_Xt,*tmpI;
  matrix *tXt_Xt,*tXt_Xt_I,*tXt_Xt_tmp;
  matrix *dN,*Phi_dN,*Q_dN,*Xt_tXt_Xt_I_tXt_dN;
  matrix *Ubeta_tilv,*Ubeta,*tZ_Phi_dN,
    *tZ_PhiXt_tXt_Xt_I_tXt_dN,*beta_m,*beta_m1,*dU_dbeta_I_Ubeta;
  matrix *tZ_PhiXt,*tZ_PhiXt_tXt_Xt_I;
  matrix *tXt_dN,*tXt_Xt_I_tXt_dN;
  matrix *tmp2_mat,*tmp11_mat,*tmp12_mat;
  matrix *Psi,*Ufunk,*C1,*dC1,*tdC1, *tC1,*C1_dU_I;
  matrix *dM1M2,*dM1M2tmp,*dM2,*tdM2;
  matrix *M1M2,*M2, *M1,*M1tmp1,*M1tmp2, *C1M1M2,*dC1M1M2,*tC1M1M2 ,*C1M1, *C1M1tC1,
    *VarPsi,*VarPsi1,*VarPsi2,*VarPsi_out;
  matrix *W1,*W2t[*n],*Ut[*n],*Ui[*k],*Utm[*n],*dU_dbeta_i[*k],*dW2,*Delta,*Delta1,*tmpM1,*Delta2,*tmpM2;
  double norm_rand();
  void GetRNGstate(),PutRNGstate();
  

  malloc_vecs((*p+1),&vtmp1,&vtmp3,&testOBS,&testtmp,&testOBSHW,&testOBSCM,
	     &testtmpHW,&testtmpCM,&testtmp1,&testtmp2,NULL);  
  malloc_vecs(*q,&vtmp2,&vtmp4,&vtmp5,NULL);  
  malloc_vecs((*p+1)*(*p+1),&vtmp7,NULL);  
  malloc_vecs(1,&vtmp6,NULL);  
  malloc_vecs(*k,&Ufunkdim1,&Ufunkdim12,&Delta2tmp,&Delta2tmp1,NULL); 
  malloc_vecs(*antsim,&testGOFCM,&testtmp1GOFCM,NULL);
  malloc_mats(*antsim,(*p+1),&testOBS1,&testHW,&testCM,NULL);
  malloc_mats((*p+1),1,&mtmp1,&mtmp6,NULL);
  malloc_mats(1,(*p+1),&tmtmp1,NULL);
  malloc_mats(1,1,&mtmp4,&mtmp5,NULL);
  malloc_mats(*q,1,&mtmp2,&mtmp3,&mtmp3m,&mtmp3mm,&mtmp3mmm,NULL);
  malloc_mats(1,*q,&mtmp2t,NULL);
  malloc_mats(*n,(*p+1),&Xt,&PhiXt,&tmp10_mat,&DdN_Xt,&dM2m,&dW2,NULL);
  malloc_mats((*p+1),*n,&tXt,&tDdN_Xt,&tdM2m,NULL);
  malloc_mats(*n,*p,&X,NULL);
  malloc_mats(*n,*q,&Z,&Phi_Z,&DdN_Phi_Z,&tmp3_mat,&tmp4_mat,&tmp5_mat,
		  &tmp6_mat,&ttmpI,&W1,&tmp12_mat,NULL); // changed dims of tmp12_mat
  malloc_mats(*q,*q,&dU_dbeta,&dU_dbeta1,&dU_dbeta_tmp,&dU_dbeta_I,&dU_dbeta_I1,&tmp7_mat,&dU_dbeta_tilv,
	     &tmp13_mat,&opt_tilv,&opt,&dM1,&dM1tmp,&M1,&M1tmp1,&M1tmp2,NULL);
  malloc_mats(*q,*n,&tZ,&tZ_Phi,&tDdN_Phi_Z,
	     &tZ_PhiXt_tXt_Xt_I_tDdN_Xt,&tmpI,NULL);
  malloc_mats((*p+1),(*p+1),&tXt_Xt,&tXt_Xt_tmp,&tXt_Xt_I,NULL);
  malloc_mats(*n,1,&dN,&Phi_dN,&Q_dN,&Xt_tXt_Xt_I_tXt_dN,NULL);
  malloc_mats(*q,1,&Ubeta_tilv,&Ubeta,&tZ_Phi_dN,
	     &tZ_PhiXt_tXt_Xt_I_tXt_dN,&beta_m,&beta_m1,&dU_dbeta_I_Ubeta,NULL);
  malloc_mats(*q,(*p+1),&tZ_PhiXt,&tZ_PhiXt_tXt_Xt_I,&dM1M2,&dM1M2tmp,&M1M2,&tdC1,&tC1,NULL);
  malloc_mats((*p+1),1,&tXt_dN,&tXt_Xt_I_tXt_dN,NULL);
  malloc_mats((*p+1),*q,&tmp2_mat,&tmp11_mat,&dC1,&C1,&C1M1,&C1_dU_I,NULL); // changed the dims of tmp12_mat
  malloc_mats((*p+1),*k,&Psi,NULL);
  malloc_mats((*p+1),(*p+1),&dM2,&M2,&tdM2,&C1M1M2,&dC1M1M2,&tC1M1M2,&C1M1tC1,
	     &VarPsi,&VarPsi1,&VarPsi2,NULL);
  malloc_mats(((*p+1)*(*p+1)),*k,&VarPsi_out,NULL); /* OBS OBS */

  for (j=0;j<*n;j++) { 
    malloc_mat(*k,*p+1,W2t[j]);    malloc_mat(*k,*q,Ut[j]);    malloc_mat(*k,*q,Utm[j]);    
  }
  for (j=0;j<*k;j++) { 
    malloc_mat(*q,*q,dU_dbeta_i[j]);  malloc_mat(*q,1,Ui[j]);  
  }
  malloc_mats(*k,*p+1,&Delta,&tmpM1,&Delta1,NULL);
  malloc_mats(*k,*q,&Delta2,&tmpM2,&Ufunk,NULL);




  /*for (l=0;l<*q;++l){ U_beta[l]=0;}*/
  /* for (l=0;l<((*q)*(*q));++l){dU_beta[l]=0;}*/

  /** Newton-iteration **/

  for (l=0;l<*q;++l){beta_tmp[l]=beta[l];} /* beta_0 */
//  eps=0.01;
  tau=0;
  dum1=0; dum1=dum1+0; 
  nummer=50; /* Antal iterationer der højst udføres */
  del=1000;
  del_old=10000;
  it=0;  
//  kny=0; /* max(tau_i) for hvilket at design ej sing. */
  alpha_tmp=(*alpha);
  init_it=2;


  while ( (del>*tol) && (it<nummer)&&(dum1<1)){

    *alpha=(it<init_it) ? 0.3:alpha_tmp;

    /* Rprintf("%2d    \n",it);*/

    /** Beregning af U, dU og [U] **/


    /*** Initialiseringer  ****/


    for (j=0;j<*n;++j){ 
      betaZ[j]=0;
      for (l=0;l<*q;++l){betaZ[j]+=excess[j]*
	  Zinp[l*(*n)+j]*beta[l];}
      phi[j]=excess[j]*exp(betaZ[j]);}
             
    for (l=0;l<*q;++l){ Uinp[l]=0;} /* Score */
    for (l=0;l<*q;++l){for (l1=0;l1<*q;++l1){
	dUinp[l*(*q)+l1]=0;}} /* dU */
    for (l=0;l<*q;++l){for (l1=0;l1<*q;++l1){
	optinp[l*(*q)+l1]=0;}}/* [U] */

    mat_zeros(dU_dbeta); mat_zeros(dU_dbeta1); 
    mat_zeros(Ufunk); mat_zeros(opt); vec_zeros(Ufunkdim1);
    
    for (i=0;i<*k;++i){  /* Start gennemløb tau_1 til tau_k */
   

      for (j=0;j<*n;++j){
	y[j]=(time[j]>=stid[i]) ? 1:0;
	y1[j]=(time[j]>=stid[i]) ? 1:0;          
	y2[j]=(time[j]<=stid[i]) ? 1:0;
	ME(dN,j,0)=y1[j]*y2[j];     
	ME(Phi_dN,j,0)=phi[j]*y1[j]*y2[j];
      }

      /******* Konstruktion design-matricer *******/
      for (l=0;l<*p;++l){ 
	for (j=0;j<*n;++j){
	  ME(X,j,l)=y[j]*Xinp[l*(*n)+j];
	}
      }
      for (l=0;l<(*p+1);++l){ 
	for (j=0;j<*n;++j){
	  ME(Xt,j,l)=(l<(*p)) ? y[j]*Xinp[l*(*n)+j]:y[j]*phi[j]; /** X.tilde **/
	  ME(PhiXt,j,l)=(l<(*p)) ? phi[j]*y[j]*Xinp[l*(*n)+j]:phi[j]*y[j]*phi[j]; /** Phi*X.tilde **/
	  ME(DdN_Xt,j,l)=(l<(*p)) ? y1[j]*y2[j]*Xinp[l*(*n)+j]:y1[j]*y2[j]*phi[j]; /** Diag(dN)*X.tilde **/
	}
      }
      for (l=0;l<*q;++l){ 
	for (j=0;j<*n;++j){
	  ME(Z,j,l)=y[j]*excess[j]*Zinp[l*(*n)+j];
	  ME(Phi_Z,j,l)=phi[j]*y[j]*excess[j]*Zinp[l*(*n)+j];
	  ME(DdN_Phi_Z,j,l)=y1[j]*y2[j]*phi[j]*y[j]*excess[j]*Zinp[l*(*n)+j];
	}
      }

      /******** Slut konstruktion design-matricer *******/


      /******* Beregning af score U_beta ************/
      mat_transp(Z,tZ);mat_transp(Xt,tXt);

      MxA(tZ,Phi_dN,tZ_Phi_dN);
      MxA(tZ,PhiXt,tZ_PhiXt);
      MxA(tXt,Xt,tXt_Xt);
      /* mat_copy(tXt_Xt,tXt_Xt_tmp);
	 QRfactor(tXt_Xt_tmp,vtmp1);
	 a=1;
	 for (l=0;l<(*p+1);++l){a=a*ME(tXt_Xt_tmp,l,l);}
	 a=sqrt(a*a);
	 b=(a>eps) ? 1:0;*/
      b=1;
      /*mat_zeros(tXt_Xt_I);*/
      /*  if (a>eps){*/
      invert(tXt_Xt,tXt_Xt_I);
     
      MxA(tXt,dN,tXt_dN);

      MxA(tXt_Xt_I,tXt_dN,tXt_Xt_I_tXt_dN);
      MxA(tZ_PhiXt,tXt_Xt_I_tXt_dN,tZ_PhiXt_tXt_Xt_I_tXt_dN);

      mat_subtr(tZ_Phi_dN,tZ_PhiXt_tXt_Xt_I_tXt_dN,Ubeta_tilv);
      scl_mat_mult(b,Ubeta_tilv,Ubeta_tilv);

      if (i<1){mat_copy(Ubeta_tilv,Ui[i]);}
      if (i>0){mat_add(Ubeta_tilv,Ui[i-1],Ui[i]);} /* Ui[i] er værd. af scoren til tid tau_i*/

      for (l=0;l<*q;++l){ 
	Uinp[l]+=ME(Ubeta_tilv,l,0);}

      /* Rprintf("%2d %6.4f   \n",i,Uinp[0]);*/

      /******* Slut beregning af score U_beta ************/


      /****** Beregning af dU_dbeta************/
      mat_transp(Phi_Z,tZ_Phi);

      tmp1_sc=ME(tXt_Xt_I_tXt_dN,(*p),0);

      MxA(Xt,tXt_Xt_I_tXt_dN,Xt_tXt_Xt_I_tXt_dN);
      mat_subtr(dN,Xt_tXt_Xt_I_tXt_dN,Q_dN);



      mat_zeros(tmp2_mat);
      for (l=0;l<*q;++l){for (j=0;j<*n;++j){
	  ME(tmp2_mat,(*p),l)=ME(tmp2_mat,(*p),l)+
	    phi[j]*(ME(Z,j,l))*(ME(Q_dN,j,0));}}


      scl_mat_mult(tmp1_sc,Phi_Z,tmp3_mat);

      MxA(Xt,tXt_Xt_I,tmp10_mat);
      MxA(tXt,tmp3_mat,tmp11_mat);
      MxA(tmp10_mat,tmp11_mat,tmp12_mat); 

      mat_subtr(tmp3_mat,tmp12_mat,tmp4_mat);
      MxA(tmp10_mat,tmp2_mat,tmp5_mat);
      /*mat_subtr(tmp4_mat,tmp5_mat,tmp6_mat);*/
      mat_add(tmp4_mat,tmp5_mat,tmp6_mat);

      mat_zeros(tmp7_mat);
      for (l=0;l<*q;++l){for (l1=0;l1<*q;++l1){for (j=0;j<*n;++j){
			     ME(tmp7_mat,l,l1)=ME(tmp7_mat,l,l1)+
			       phi[j]*(ME(Z,j,l))*(ME(Z,j,l1))*(ME(Q_dN,j,0));}}}
      /* print_mat(tmp7_mat);*/

      MxA(tZ_Phi,tmp6_mat,tmp13_mat);  
      /*m_mlt(tZ_Phi,tmp5_mat,tmp13_mat); if (i>(*k-2)){m_output(tmp13_mat);}*/ 

      MxA(tZ_Phi,tmp4_mat,dU_dbeta_tilv);scl_mat_mult(-1,dU_dbeta_tilv,dU_dbeta_tilv);/* Hovedleddet i dU_dbeta*/
      scl_mat_mult(b,dU_dbeta_tilv,dU_dbeta_tilv);
      mat_subtr(tmp7_mat, tmp13_mat,dU_dbeta_tilv);
      scl_mat_mult(b,dU_dbeta_tilv,dU_dbeta_tilv);  
      mat_add(dU_dbeta_tilv,dU_dbeta1,dU_dbeta1); 
      mat_copy(dU_dbeta1,dU_dbeta_i[i]);
      /*  if (i>(*k-2)) {print_mat(dU_dbeta1); }*/
  


      mat_subtr(tmp7_mat, tmp13_mat,dU_dbeta_tilv);
      scl_mat_mult(b,dU_dbeta_tilv,dU_dbeta_tilv);  
      mat_add(dU_dbeta_tilv,dU_dbeta,dU_dbeta);
      /* if (i>(*k-2)) {print_mat(dU_dbeta); }  */


      /****** Slut beregning af dU_dbeta************/

      /***** Beregning af [Ubeta] ************/

      mat_transp(DdN_Xt,tDdN_Xt);  
      mat_transp(DdN_Phi_Z,tDdN_Phi_Z);  

      MxA(tZ_PhiXt,tXt_Xt_I,tZ_PhiXt_tXt_Xt_I);
      MxA(tZ_PhiXt_tXt_Xt_I,tDdN_Xt,tZ_PhiXt_tXt_Xt_I_tDdN_Xt);
      mat_subtr(tDdN_Phi_Z,tZ_PhiXt_tXt_Xt_I_tDdN_Xt,tmpI);
      mat_transp(tmpI,ttmpI);
      MxA(tmpI,ttmpI,opt_tilv);
      scl_mat_mult(b,opt_tilv,opt_tilv);
      for (l=0;l<*q;++l){for (l1=0;l1<*q;++l1){
	  optinp[l*(*q)+l1]+=ME(opt_tilv,l1,l);}}
      mat_add(opt,opt_tilv,opt);
  
  

                                                            
      /*if (i<1){Rprintf("%6.4f   \n",Uinp[0]);print_mat(opt);}*/
  
      /***** Slut beregning af [Ubeta] ************/

      /** Start beregning af ene komponent i robust varians***/
      for (l=0;l<*n;++l){
        extract_row(Xt,l,vtmp3);scl_vec_mult(ME(dN,l,0),vtmp3,vtmp3);
        replace_col(mtmp1,0,vtmp3); 
        MxA(tZ_PhiXt_tXt_Xt_I,mtmp1,mtmp2);

        extract_row(Phi_Z,l,vtmp4);scl_vec_mult(ME(dN,l,0),vtmp4,vtmp4);
        replace_col(mtmp3,0,vtmp4);
        mat_subtr(mtmp3,mtmp2,mtmp3);extract_col(mtmp3,0,vtmp4);


	for (l1=0;l1<*q;++l1){
	  ME(Utm[l],i,l1)=(i>0) ? (VE(vtmp4,l1)+ME(Utm[l],i-1,l1)):VE(vtmp4,l1);
	}
	/* tau_2,...,tau_k */ 
    
	if (i<1){replace_row(W1,l,vtmp4);} /* tau_1 */
        if (i>0){for (l1=0;l1<*q;++l1){
	    ME(W1,l,l1)=ME(W1,l,l1)+VE(vtmp4,l1);}} /* tau_2,...,tau_k */      
      }
      /** Slut beregning af ene komponent i robust varians***/
           
    }
    /* Slut gennemløb tau_1 til tau_k */
    /** Slut beregning af U, dU og [U] **/


    for (l=0;l<*q;++l){ME(beta_m,l,0)=beta[l];}
    /*Rprintf("%6.4f   \n",beta[0]);*/
    for (l=0;l<*q;++l){ME(Ubeta,l,0)=Uinp[l];}
    for (l=0;l<*q;++l){for (l1=0;l1<*q;++l1){
	dUinp[l*(*q)+l1]=ME(dU_dbeta,l1,l);}}/* dU_dbeta */
    tau=(it<init_it)? 1:0; 
    for (l=0;l<*q;++l){
      ME(dU_dbeta,l,l)=ME(dU_dbeta,l,l)-tau;} /* tau subtraheres til diagonalen af dU_dbeta*/
    /* mat_copy(dU_dbeta,dU_dbeta_tmp);   
       QRfactor(dU_dbeta_tmp,vtmp2);
       a=1;
       for (l=0;l<*q;++l){a=a*ME(dU_dbeta_tmp,l,l);}
       c=sqrt(a*a); 
       if (c>eps){*/
    invert(dU_dbeta,dU_dbeta_I);invert(dU_dbeta1,dU_dbeta_I1);
    /*Rprintf("%6.4f   \n",c);*/
    /*print_mat(dU_dbeta_I1);  print_mat(dU_dbeta_I);*/
    for (i=0;i<*k;++i){
      VE(Ufunkdim1,i)=0; VE(Ufunkdim12,i)=0;
      for (l=0;l<*q;++l){
	ME(Ufunk,i,l)=(ME(Ui[i],l,0))*sqrt(fabs(ME(dU_dbeta_I1,l,l)));
	VE(Ufunkdim1,i)=VE(Ufunkdim1,i)+fabs(ME(Ufunk,i,l));
	VE(Ufunkdim12,i) = VE(Ufunkdim12,i)+ME(Ufunk,i,l);
      }/* Normeret score til tid tau_i     */   


      for (l=0;l<*n;++l){
	for (l1=0;l1<*q;++l1){
	  ME(mtmp3m,l1,0)=ME(Utm[l],(*k-1),l1);}
	MxA(dU_dbeta_I1,mtmp3m,mtmp3mm); 
	MxA(dU_dbeta_i[i],mtmp3mm,mtmp3mmm);
	for (l1=0;l1<*q;++l1){
	  ME(Ut[l],i,l1)=(ME(Utm[l],i,l1)-ME(mtmp3mmm,l1,0)); 
	}
      }

    }


    scl_mat_mult(*alpha,dU_dbeta_I,dU_dbeta_I); /* alpha er skridtlængden */
    /* i Newton-iterationen   */   
    mat_copy(opt,M1tmp1);
    MxA(dU_dbeta_I,M1tmp1,M1tmp2);
    MxA(M1tmp2,dU_dbeta_I,M1);/*  i [M_1](tau_k) */
    MxA(dU_dbeta_I,Ubeta,dU_dbeta_I_Ubeta);
    mat_subtr(beta_m,dU_dbeta_I_Ubeta,beta_m1);
    /*print_mat(beta_m1);*/
    if (del<*tol){dum1=2;}
    del=0;
    for (l=0;l<*q;++l){del+=Uinp[l]*Uinp[l];}
    del=sqrt(del);
   
    for (l=0;l<*q;++l){beta[l]=ME(beta_m1,l,0);}
    /*if (it<7){ Rprintf("  it alpha beta del \t");
      Rprintf("%2d %6.4f  %6.4f %6.4f   \n",it,*alpha,beta[0],del);
      }*/
    if (del<*tol){for (l=0;l<*q;++l){beta[l]=ME(beta_m,l,0);}}
    /* Rprintf("  it alpha beta del \t");
       Rprintf("%2d %6.4f  %6.4f %6.4f   \n",it,*alpha,beta[0],del);*/
    if ((del_old<del)&&(it>init_it) ){ 
      alpha_tmp=0.67*(*alpha);
      for (l=0;l<*q;++l){beta[l]=beta_tmp[l];}}
    del_old=del;
    it=it+1;
  }

  /*   for (l=0;l<*q;++l){beta[l]=ME(beta_m,l,0);}*/

  /** Slut Newton-iteration **/

  /*print_vec(Ufunkdim1);*/

  /** Start beregning af hat Psi samt varians **/

  /*** Initialiseringer  ****/

  for (j=0;j<*n;++j){ 
    betaZ[j]=0;
    for (l=0;l<*q;++l){betaZ[j]+=excess[j]*
	Zinp[l*(*n)+j]*beta[l];}
    phi[j]=excess[j]*exp(betaZ[j]);}
             

  for (l=0;l<(*p+1);++l){for (i=0;i<*k;++i){
      ME(Psi,l,i)=0;}
  }

  for (l=0;l<(*p+1);++l){for (i=0;i<*q;++i){
      ME(C1,l,i)=0;}
  }

  /*   for (l=0;l<*q;++l){for (i=0;i<*q;++i){
       ME(M1,l,i)=0;}
       }*/

  for (l=0;l<(*p+1);++l){for (i=0;i<(*p+1);++i){
      ME(M2,l,i)=0;}
  }

  for (l=0;l<*q;++l){for (i=0;i<(*p+1);++i){
      ME(M1M2,l,i)=0;}
  }


  for (i=0;i<*k;++i){  /* Start gennemløb tau_1 til tau_k */

    /*** Initialiseringer  ****/


    /*** Slut initialiseringer  ****/

    
    for (j=0;j<*n;++j){
      y[j]=(time[j]>=stid[i]) ? 1:0;
      y1[j]=(time[j]>=stid[i]) ? 1:0;          
      y2[j]=(time[j]<=stid[i]) ? 1:0;
      ME(dN,j,0)=y1[j]*y2[j];     
      ME(Phi_dN,j,0)=phi[j]*y1[j]*y2[j];
    }


    /******* Konstruktion design-matricer *******/
    for (l=0;l<*p;++l){ 
      for (j=0;j<*n;++j){
	ME(X,j,l)=y[j]*Xinp[l*(*n)+j];
      }
    }
    for (l=0;l<(*p+1);++l){ 
      for (j=0;j<*n;++j){
	ME(Xt,j,l)=(l<(*p)) ? y[j]*Xtilinp[l*(*n)+j]:y[j]*phi[j]; /** X.tilde **/
	ME(PhiXt,j,l)=(l<(*p)) ? phi[j]*y[j]*Xtilinp[l*(*n)+j]:phi[j]*y[j]*phi[j]; /** Phi*X.tilde **/
	ME(DdN_Xt,j,l)=(l<(*p)) ? y1[j]*y2[j]*Xtilinp[l*(*n)+j]:y1[j]*y2[j]*phi[j]; /** Diag(dN)*X.tilde **/
      }
    } 
 

    for (l=0;l<*q;++l){ 
      for (j=0;j<*n;++j){
	ME(Z,j,l)=y[j]*excess[j]*Zinp[l*(*n)+j]; 
	ME(Phi_Z,j,l)=phi[j]*y[j]*excess[j]*Zinp[l*(*n)+j];
	ME(DdN_Phi_Z,j,l)=y1[j]*y2[j]*phi[j]*y[j]*excess[j]*Zinp[l*(*n)+j];
      }
    }

    
    /******** Slut konstruktion design-matricer *******/


    /******* Beregning af score U_beta ************/

    mat_transp(Z,tZ);  
    mat_transp(Xt,tXt);
    

    MxA(tZ,Phi_dN,tZ_Phi_dN);
    MxA(tZ,PhiXt,tZ_PhiXt);
    MxA(tXt,Xt,tXt_Xt);/*print_mat(tXt_Xt);*/
    mat_copy(tXt_Xt,tXt_Xt_tmp);
    /*   QRfactor(tXt_Xt_tmp,vtmp1);
	 a=1;
	 for (l=0;l<(*p+1);++l){a=a*ME(tXt_Xt_tmp,l,l);}
	 a=sqrt(a*a);

	 b=(a>eps) ? 1:0;*/
    b=1;
    /*  kny=(a>eps) ? i:kny;   
	if (a>eps){*/

    invert(tXt_Xt,tXt_Xt_I);

    MxA(tXt,dN,tXt_dN);

    MxA(tXt_Xt_I,tXt_dN,tXt_Xt_I_tXt_dN);
    MxA(tZ_PhiXt,tXt_Xt_I_tXt_dN,tZ_PhiXt_tXt_Xt_I_tXt_dN);

    mat_subtr(tZ_Phi_dN,tZ_PhiXt_tXt_Xt_I_tXt_dN,Ubeta_tilv);

    for (l=0;l<*q;++l){ 
      Uinp[l]+=b*ME(Ubeta_tilv,l,0);}


    /******* Slut beregning af score U_beta ************/


    /****** Beregning af dU_dbeta************/
    mat_transp(Phi_Z,tZ_Phi);

    tmp1_sc=ME(tXt_Xt_I_tXt_dN,(*p),0);

    MxA(Xt,tXt_Xt_I_tXt_dN,Xt_tXt_Xt_I_tXt_dN);
    mat_subtr(dN,Xt_tXt_Xt_I_tXt_dN,Q_dN);

    mat_zeros(tmp2_mat);
    for (l=0;l<*q;++l){for (j=0;j<*n;++j){
	ME(tmp2_mat,(*p),l)+=phi[j]*(ME(Z,j,l))*(ME(Q_dN,j,0));}}

    scl_mat_mult(tmp1_sc,Phi_Z,tmp3_mat);

    MxA(Xt,tXt_Xt_I,tmp10_mat);
    MxA(tXt,tmp3_mat,tmp11_mat);
    MxA(tmp10_mat,tmp11_mat,tmp12_mat);

    mat_subtr(tmp3_mat,tmp12_mat,tmp4_mat);
    MxA(tmp10_mat,tmp2_mat,tmp5_mat);
    mat_add(tmp4_mat,tmp5_mat,tmp6_mat);
   
    mat_zeros(tmp7_mat);
    for (l=0;l<*q;++l){for (l1=0;l1<*q;++l1){for (j=0;j<*n;++j){
			   ME(tmp7_mat,l,l1)+=phi[j]*(ME(Z,j,l))*
			     (ME(Z,j,l1))*(ME(Q_dN,j,0));}}}

    MxA(tZ_Phi,tmp6_mat,tmp13_mat);
    mat_subtr(tmp7_mat, tmp13_mat,dU_dbeta_tilv);
    scl_mat_mult(b,dU_dbeta_tilv,dU_dbeta_tilv); 

    /*m_mlt(tZ_Phi,tmp4_mat,dU_dbeta_tilv);sm_mlt(-1,dU_dbeta_tilv,dU_dbeta_tilv); 
      sm_mlt(b,dU_dbeta_tilv,dU_dbeta_tilv);*/ 
    mat_add(dU_dbeta_tilv,dU_dbeta,dU_dbeta);

    /****** Slut beregning af dU_dbeta************/

    /***** Beregning af [Ubeta] ************/

    mat_transp(DdN_Xt,tDdN_Xt);  
    mat_transp(DdN_Phi_Z,tDdN_Phi_Z);  

    MxA(tZ_PhiXt,tXt_Xt_I,tZ_PhiXt_tXt_Xt_I);
    MxA(tZ_PhiXt_tXt_Xt_I,tDdN_Xt,tZ_PhiXt_tXt_Xt_I_tDdN_Xt);
    mat_subtr(tDdN_Phi_Z,tZ_PhiXt_tXt_Xt_I_tDdN_Xt,tmpI);


    /***** Slut beregning af [Ubeta] ************/
           
    
    /** Beregning af Psi **/
  
    for (l=0;l<(*p+1);++l){
      ME(Psi,l,i)=(i>0) ? ME(Psi,l,i-1)+b*ME(tXt_Xt_I_tXt_dN,l,0):b*ME(tXt_Xt_I_tXt_dN,l,0);
    }

    /** Slut beregning Psi **/

    /** Beregning af Var(\hat Psi) **/

    scl_mat_mult(ME(tXt_Xt_I_tXt_dN,*p,0),tZ_PhiXt_tXt_Xt_I,tdC1);
    mat_transp(tdC1,dC1);
    scl_mat_mult(b,dC1,dC1);
    mat_add(C1,dC1,C1); /* C1(tau_i) */
    mat_transp(C1,tC1);

    MxA(DdN_Xt,tXt_Xt_I,dM2m);mat_transp(dM2m,tdM2m);
    MxA(tdM2m,dM2m,dM2); /* tilvækst i [M_2] */
    scl_mat_mult(b,dM2,dM2);
    mat_add(M2,dM2,M2);     /* [M_2](tau_i) */


    MxA(C1,M1,C1M1);
    MxA(C1M1,tC1,C1M1tC1);  /* C1(tau_i)[M_1](tau_k) C1^T(tau_i) */


    MxA(tmpI,dM2m,dM1M2tmp);  
    MxA(dU_dbeta_I,dM1M2tmp,dM1M2);  /* tilvækst i [M_1,M_2] */   
    scl_mat_mult(b,dM1M2,dM1M2);
    mat_add(M1M2,dM1M2,M1M2);  /*[M_1,M_2](tau_i) */
    MxA(C1,M1M2,C1M1M2);  /* C1(tau_i)[M_1,M_2](tau_i) */
    mat_transp(C1M1M2,tC1M1M2); 
    mat_add(C1M1M2,tC1M1M2,dC1M1M2);


    mat_add(M2,C1M1tC1,VarPsi1);
    mat_subtr(VarPsi1,dC1M1M2,VarPsi);/* Est. Var(hat Phi)(tau_i) */
    /* en (p+1)x(p+1)-matrix    */


    for (l1=0;l1<(*p+1);++l1){
      for (l2=0;l2<(*p+1);++l2){  
	ME(VarPsi_out,l1*(*p+1)+l2,i)=ME(VarPsi,l2,l1);}}
    
    /*  print_mat(VarPsi);*/
 

 
    /** Start beregning af anden komponent i robust varians***/
    vec_zeros(vtmp7);

    for (l=0;l<*n;++l){
      extract_row(Xt,l,vtmp3);
      scl_mat_mult(b,tXt_Xt_I_tXt_dN,tXt_Xt_I_tXt_dN);
      vM(tXt_Xt_I_tXt_dN,vtmp3,vtmp6);
      VE(vtmp6,0)=y1[l]*y2[l]-VE(vtmp6,0);
      scl_vec_mult(VE(vtmp6,0),vtmp3,vtmp3);
      replace_col(mtmp1,0,vtmp3);
      MxA(tXt_Xt_I,mtmp1,mtmp6);
         
      for (l1=0;l1<(*p+1);++l1){
	ME(dW2,l,l1)=(i<1) ? ME(mtmp6,l1,0):ME(dW2,l,l1)+ME(mtmp6,l1,0);}
           
	 
      MxA(C1,dU_dbeta_I,C1_dU_I);   
      extract_row(W1,l,vtmp5);Mv(C1_dU_I,vtmp5,vtmp3);
      extract_row(dW2,l,vtmp1);
      vec_subtr(vtmp1,vtmp3,vtmp1);
      replace_row(W2t[l],i,vtmp1);
   


      for (l1=0;l1<(*p+1);++l1){
	for (l2=0;l2<(*p+1);++l2){  
	  VE(vtmp7,l1*(*p+1)+l2)=VE(vtmp7,l1*(*p+1)+l2)+(VE(vtmp1,l1))*(VE(vtmp1,l2));}}
    }

    /*print_vec(vtmp7);*/

    for (l1=0;l1<(*p+1);++l1){
      for (l2=0;l2<(*p+1);++l2){  
	ME(VarPsi_out,l1*(*p+1)+l2,i)=VE(vtmp7,l1*(*p+1)+l2);}}

    /** Slut beregning af Var(\hat Psi) **/

    /*print_mat(VarPsi);*/

  }

  /* Slut gennemløb tau_1 til tau_k */
     
 
  /* Beregning af robust varians af beta */
  mat_zeros(M1);mat_zeros(M1tmp1);mat_zeros(M1tmp2);

  for (l=0;l<*n;++l){
    extract_row(W1,l,vtmp5);
    replace_col(mtmp2,0,vtmp5);mat_transp(mtmp2,mtmp2t);
    MxA(mtmp2,mtmp2t,tmp7_mat);
    mat_add(M1tmp1,tmp7_mat,M1tmp1);
  }



  mat_transp(dU_dbeta_I,dU_dbeta_tmp);
  MxA(dU_dbeta_I,M1tmp1,M1tmp2);
  MxA(M1tmp2,dU_dbeta_tmp,M1);

  /* print_mat(M1); */
 
  /* Slut beregning af robust varians af beta */


  /* Beregning af obs. teststørrelse */

  for (l=0;l<*p+1;++l){VE(testOBS,l)=0;VE(testOBSHW,l)=0;VE(testOBSCM,l)=0;
    for (i=0;i<*k1;++i){
      VE(testtmp,l)=fabs(ME(Psi,l,i))/sqrt(ME(VarPsi_out,l*(*p+1+1),i));
      VE(testOBS,l)=(VE(testOBS,l)>VE(testtmp,l)) ? VE(testOBS,l):VE(testtmp,l);
      VE(testtmpHW,l)=fabs(ME(Psi,l,i))*sqrt(ME(VarPsi_out,l*(*p+1+1),*k1-1))/
	(ME(VarPsi_out,l*(*p+1+1),i)+ME(VarPsi_out,l*(*p+1+1),*k1-1));
      VE(testtmpCM,l)=(ME(Psi,l,i))*(ME(Psi,l,i))/ME(VarPsi_out,l*(*p+1+1),i);
      VE(testOBSHW,l)=(VE(testOBSHW,l)>VE(testtmpHW,l)) ? VE(testOBSHW,l):VE(testtmpHW,l);
      VE(testOBSCM,l)=VE(testOBSCM,l)+VE(testtmpCM,l);
    }  
  }
 
  testOBSGOFCM=0;
  for (i=0;i<*k;++i){
    testOBSGOFCM=(testOBSGOFCM> sqrt(VE(Ufunkdim1,i)*VE(Ufunkdim1,i))) ? 
      testOBSGOFCM:sqrt(VE(Ufunkdim1,i)*VE(Ufunkdim1,i));
  }

  /*  Rprintf("%6.4f   \n",testOBSGOFCM);*/

  /***  Start simulationer   ***/
  Rprintf("Simulations start N=\t");
  Rprintf("%2d    \n",*antsim);
  GetRNGstate();  /* to use R random normals */


  for (j=0;j<*antsim;++j){ 

    mat_zeros(Delta);mat_zeros(Delta2);
    for (i=0;i<*n;++i){ 
      random=norm_rand();
      scl_mat_mult(random,W2t[i],tmpM1);mat_add(tmpM1,Delta,Delta); 
      scl_mat_mult(random,Ut[i],tmpM2);mat_add(tmpM2,Delta2,Delta2); 
    }/* Delta indeholder Simuleret Psi(t) */
    /* Delta indeholder Simuleret U(t) */
    vec_zeros(Delta2tmp);       vec_zeros(Delta2tmp1);
    for (i=0;i<*k;++i){ 
      for (l=0;l<*q;++l){ 
	VE(Delta2tmp,i)=VE(Delta2tmp,i)+fabs(ME(Delta2,i,l))*
	  sqrt(fabs(ME(dU_dbeta_I1,l,l)));
	VE(Delta2tmp1,i)=VE(Delta2tmp1,i)+ME(Delta2,i,l)*
	  sqrt(fabs(ME(dU_dbeta_I1,l,l)));
      }
    }

    if (j<51) { /*if (j<1){print_vec(Ufunkdim1);} if (j<6){print_vec(Delta2tmp);}*/
      for (l1=0;l1<*q;++l1){  
	for (l=0;l<*k;++l){  
	  Scoreinp[(j*(*q)+l1)*(*k)+l]=(j<1) ? ME(Ui[l],l1,0):ME(Delta2,l,l1);
	}}}

    /*Ui[i] er værd. af scoren til tid tau_i*/ 

    for (l=0;l<*p+1;++l){ME(testOBS1,j,l)=0;ME(testHW,j,l)=0;ME(testCM,j,l)=0;
      for (i=0;i<*k1;++i){
	VE(testtmp,l)=fabs(ME(Delta,i,l))/sqrt(ME(VarPsi_out,l*(*p+1+1),i));
	/* testtmp: Simulations baseret konf.bånd */
	VE(testtmp1,l)=fabs(ME(Delta,i,l))*sqrt(ME(VarPsi_out,l*(*p+1+1),*k1-1))/
	  (ME(VarPsi_out,l*(*p+1+1),i)+ME(VarPsi_out,l*(*p+1+1),*k1-1));
	/* testtmp1: Simulering af Hall-Wellner band */
	VE(testtmp2,l)=(ME(Delta,i,l))*(ME(Delta,i,l))/ME(VarPsi_out,l*(*p+1+1),i);
        ME(testOBS1,j,l)=(ME(testOBS1,j,l)>VE(testtmp,l)) ? ME(testOBS1,j,l):VE(testtmp,l);
        ME(testHW,j,l)=(ME(testHW,j,l)>VE(testtmp1,l)) ? ME(testHW,j,l):VE(testtmp1,l);
        ME(testCM,j,l)=ME(testCM,j,l)+VE(testtmp2,l);
      }
    } 
    
    VE(testGOFCM,j)=0;
    for (i=0;i<*k;++i){
      VE(testGOFCM,j)=(VE(testGOFCM,j)>sqrt(VE(Delta2tmp,i)*VE(Delta2tmp,i))) ?
	VE(testGOFCM,j):sqrt(VE(Delta2tmp,i)*VE(Delta2tmp,i));
    }

    for (l=0;l<(*p+1);++l){
      testinp[j*(*p+1)+l]=(j<1) ? VE(testOBS,l):ME(testOBS1,j-1,l);
      testinpHW[j*(*p+1)+l]=(j<1) ? VE(testOBSHW,l):ME(testHW,j-1,l);
      testinpCM[j*(*p+1)+l]=(j<1) ? VE(testOBSCM,l):ME(testCM,j-1,l);
    }
    testinpGOFCM[j]=(j<1) ? testOBSGOFCM:VE(testGOFCM,j-1);     

 
  }
  /* Slut simulationer   */
  PutRNGstate();  /* to use R random normals */
   
  /*print_vec(testGOFCM);*/
  for (l=0;l<(*p+1);++l){
    testinp[(*antsim)*(*p+1)+l]=ME(testOBS1,*antsim-1,l);
    testinpHW[(*antsim)*(*p+1)+l]=ME(testHW,*antsim-1,l);
    testinpCM[(*antsim)*(*p+1)+l]=ME(testCM,*antsim-1,l);
  }
  testinpGOFCM[*antsim]=VE(testGOFCM,*antsim-1);     




  for (l=0;l<*q;++l){
    for (l1=0;l1<*q;++l1){
      optinp[l*(*q)+l1]=ME(M1,l1,l);
    }}  /* Robust varians for beta returneres */

  for (l=0;l<*k;++l){
    for (l1=0;l1<(*p+1);++l1){
      Psiinp[l*(*p+1)+l1]=ME(Psi,l1,l);
    }}  /* Psi returneres */
  
  for (l=0;l<*k;++l){
    for (l1=0;l1<(*p+1);++l1){
      for (l2=0;l2<(*p+1);++l2){  
	CoVarPsiinp[l*(*p+1)*(*p+1)+l1*(*p+1)+l2]=ME(VarPsi_out,l1*(*p+1)+l2,l);
      }}} /* Robust covar-matrix for Psi returneres */


  for (l=0;l<*k;++l){
    for (l2=0;l2<(*p+1);++l2){   
      VarPsiinp[l*(*p+1)+l2]=ME(VarPsi_out,l2*(*p+1+1),l); 
    }}  /* Robust varians for enkeltkomp. af  Psi returneres */

  *n=it;

  /*  print_mat(VarPsi_out);*/

  free_mats(&mtmp1,&mtmp2,&mtmp2t,&mtmp3,&mtmp3m,&mtmp3mm,&mtmp3mmm,
              &mtmp4,&mtmp5,&mtmp6,&tmtmp1,
              &testOBS1,&testHW,&testCM,
              &Xt,&PhiXt,&tmp10_mat,&DdN_Xt,&dM2m,&tXt,&tDdN_Xt,&tdM2m,&X,
              &Z,&Phi_Z,&DdN_Phi_Z,&tmp3_mat,&tmp4_mat,&tmp5_mat,&tmp6_mat,&ttmpI,
              &dU_dbeta,&dU_dbeta1,&dU_dbeta_tmp,&dU_dbeta_I,&dU_dbeta_I1,&tmp7_mat,&dU_dbeta_tilv,&tmp13_mat,
              &opt_tilv,&opt,&dM1,&dM1tmp,
              &tZ,&tZ_Phi,&tDdN_Phi_Z,&tZ_PhiXt_tXt_Xt_I_tDdN_Xt,&tmpI,
              &tXt_Xt,&tXt_Xt_I,&tXt_Xt_tmp,&dN,&Phi_dN,&Q_dN,&Xt_tXt_Xt_I_tXt_dN,
              &Ubeta_tilv,&Ubeta,&tZ_Phi_dN,&tZ_PhiXt_tXt_Xt_I_tXt_dN,&beta_m,&beta_m1,
              &dU_dbeta_I_Ubeta,&tZ_PhiXt,&tZ_PhiXt_tXt_Xt_I,
              &tXt_dN,&tXt_Xt_I_tXt_dN,&tmp2_mat,&tmp11_mat,&tmp12_mat,
              &Psi,&C1,&dC1,&tdC1, &tC1,&C1_dU_I,&dM1M2,&dM1M2tmp,&dM2,&tdM2,
              &M1M2,&M2, &M1,&M1tmp1,&M1tmp2, &C1M1M2,&dC1M1M2,&tC1M1M2 ,&C1M1, &C1M1tC1,
              &VarPsi,&VarPsi1,&VarPsi2,&VarPsi_out,&W1,&Delta,&Delta1,&tmpM1,&dW2,
              &Delta2,&tmpM2,&Ufunk,NULL);
  for (j=0;j<*n;j++) { free_mat(W2t[j]); free_mat(Ut[j]);free_mat(Utm[j]);}
  for (j=0;j<*k;j++) { free_mat(dU_dbeta_i[j]); free_mat(Ui[j]);}
  free_vecs(&vtmp1,&vtmp2,&vtmp3,&vtmp4,&vtmp5,&vtmp6,&vtmp7,&testtmp,
              &testtmp1,&testtmp2,&testOBS,&testOBSHW,&testOBSCM,
              &testtmpHW,&testtmpCM,&testGOFCM,&testtmp1GOFCM,
              &Delta2tmp,&Delta2tmp1,&Ufunkdim1,&Ufunkdim12,NULL);
}



