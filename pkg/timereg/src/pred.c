//#include <stdio.h>
#include <math.h>


/* organize indeces to different clusters in matrix  of size nclust x maxclust */ 

void nclusters(int  *npers,int  *clusters, int  *nclust, int  *uniqueclust, int  *mclust)
{
  int i,maxclust=0;
  for (i=0;i<*npers;i++){
      if (nclust[clusters[i]]==0) uniqueclust[0]+=1; 
      nclust[clusters[i]]+=1; 
      if (nclust[clusters[i]]>maxclust) maxclust=nclust[clusters[i]]; 
  } 
  mclust[0]=maxclust; 
}

void clusterindex(int *clusters,int *nclust,int *npers,int *idclust,int *clustsize,int *mednum,int *num)
{
  int i;
  if (*mednum==0) {
     for (i=0;i<*npers;i++){
         idclust[(clustsize[clusters[i]])*(*nclust)+clusters[i]]=i; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (i=0;i<*npers;i++){
        idclust[num[i]*(*nclust)+clusters[i]]=i; 
        clustsize[clusters[i]]+=1; 
     } 
  }
}

/* compute the values of a step function, 
   ie how many of the jumps are smaller or
   equal to the eval points  */

void sindex(int *index, double *jump, double *eval, int *njump, int *neval){
  int n,nt,i,t;

  n = *njump;
  
  nt = *neval;

  index[0] = 0;

  i = 0;
  
  for (t=0;t<nt;t++){
    
    while(jump[i]<=eval[t]
	  && i<n)
      i++;

    index[t] = i;
  }
}


void Cpred(cum,nx,px,xval,nxval,pred,tminus)
double *cum,*xval;
int *nxval,*nx,*px,*tminus,*pred;
{ // {{{
int s,c;
double timex,sc1,sc2,smax; 

smax=cum[*nx-1]; 
for (s=0;s<*nxval;s++)
{
   timex=xval[s]; 
//   pred[s]=timex; 
   c=*nx-1; 
   sc1=cum[*nx-1];  sc2=smax+xval[*nxval-1];  

   if (*tminus==0) {
   if (timex< cum[0])  { pred[s]=0; 
//	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=0; 
   }
   else if (timex> cum[*nx-1]) { pred[s]=*nx;
//	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+(*nx-1)]; 
   }
   else {
   while ((!((timex<sc2) && (timex>=sc1))) && (c>=0)) {
   /* Rprintf(" %lf %lf %lf %ld \n",timex,sc2,sc1,c); */ 
   sc1=cum[c-1];sc2=cum[c];c=c-1; } 
   /* Rprintf("før pred  %lf %lf %lf %ld \n",timex,sc2,sc1,c); */
//   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+c]; 
   pred[s]=c+1; 
   }
   } else  { // tminus=TRUE
   if (timex<= cum[0])  { pred[s]=0; 
//	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=0; 
   }
   else if (timex> cum[*nx-1]) { pred[s]=*nx;
//	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+(*nx-1)]; 
   }
   else {
   while ((!((timex<=sc2) && (timex>sc1))) && (c>=0)) {
   /* Rprintf(" %lf %lf %lf %ld \n",timex,sc2,sc1,c); */ 
   sc1=cum[c-1];sc2=cum[c];c=c-1; } 
   /* Rprintf("før pred  %lf %lf %lf %ld \n",timex,sc2,sc1,c); */
   pred[s]=c+1; 
//   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+c]; 
   }
   }

}
} // }}}

