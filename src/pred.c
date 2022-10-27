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

void clusterindex(int *clusters,int *nclust,int *npers,int *idclust,int *clustsize,int *mednum,int *num,
		  int *firstclustid)
{
  int i;
  if (*mednum==0) {
     for (i=0;i<*npers;i++){
         idclust[(clustsize[clusters[i]])*(*nclust)+clusters[i]]=i; 
         clustsize[clusters[i]]+=1; 
	 if (clustsize[clusters[i]]==1) firstclustid[clusters[i]]=i; 
      } 
  } else {
    for (i=0;i<*npers;i++){
        idclust[num[i]*(*nclust)+clusters[i]]=i; 
        clustsize[clusters[i]]+=1; 
        if (clustsize[clusters[i]]==1) firstclustid[clusters[i]]=i; 
     } 
  }
}

void atriskindex(double *start,double *stop,int *id,int *n,double *times,int *ntimes,int *nrisk,int *riskindex)
{
  int i,j;
    for (j=0;j<*ntimes;j++)
    for (i=0;i<*n;i++)
    if ((start[i]<times[j]) && (stop[i]>=times[j]))
    {
       riskindex[(nrisk[j])*(*ntimes)+j]=id[i]; 
       nrisk[j]+=1; 
    } 
}

void clusterindexdata(int *clusters,int *nclust,int *npers,int *idclust,int *clustsize,int *mednum,
		int *num,double *data, int *p,double *nydata)
{
  int i,j;
  if (*mednum==0) {
     for (i=0;i<*npers;i++){
         idclust[(clustsize[clusters[i]])*(*nclust)+clusters[i]]=i; 
     for (j=0;j<*p;j++) nydata[(clustsize[clusters[i]]*(*p)+j)*(*nclust)+clusters[i]]=data[(*npers)*j+i]; 
         clustsize[clusters[i]]+=1; 
      } 
  } else {
    for (i=0;i<*npers;i++){
        idclust[num[i]*(*nclust)+clusters[i]]=i; 
        for (j=0;j<*p;j++) nydata[(num[i]*(*p)+j)*(*nclust)+clusters[i]]=data[(*npers)*j+i]; 
        clustsize[clusters[i]]+=1; 
     } 
  }
}

/* compute the values of a step function, 
   ie how many of the jumps are smaller or
   equal to the eval points from prodlim THomas Gerds */

void sindex(int *index, double *jump, double *eval, int *N, int *NT, int *strict){
int i,t;
index[0] = 0;
i = 0;
if (*strict==0){
for (t=0;t<*NT;t++){
	     while(i<*N && jump[i]<=eval[t]) i++;
		   index[t] = i;
		       }
 }
else{
 for (t=0;t<*NT;t++){
	       while(i<*N && jump[i] < eval[t]) i++;
		     index[t] = i;
			 }
   }
}

void bubble_sort(double *val,int *list,int n)
{
  int c, d, t;
 
  // ini
  for (c = 0 ; c < ( n - 1 ); c++) list[c]=c; 

  for (c = 0 ; c < ( n - 1 ); c++)
  {
    for (d = 0 ; d < n - c - 1; d++)
    {
      if (val[list[d]] > val[list[d+1]])
      {
        /* Swapping */
        t         = list[d];
        list[d]   = list[d+1];
        list[d+1] = t;
      }
    }
  }
}


void Cpred(double *cum,int *nx,int *px,double *xval,int *nxval,double *pred,int *tminus)
{ // {{{
int j,s,c;
double timex,sc1,sc2,smax; 

smax=cum[*nx-1]; 
for (s=0;s<*nxval;s++)
{
   timex=xval[s]; pred[s]=timex; 
   c=*nx-1; 
   sc1=smax;  
   sc2=smax+1; 

   if (*tminus==0) { // {{{ 
   if (timex< cum[0])  {
//	   pred[s]=0; 
	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=0; 
   }
   else if (timex> cum[*nx-1]) { 
//	   pred[s]=*nx;
	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+(*nx-1)]; 
   }
   else {
   while ((!((timex<sc2) && (timex>=sc1))) && (c>=0)) {
   sc1=cum[c-1];sc2=cum[c];c=c-1; } 
   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+c]; 
//   pred[s]=c+1; 
   } // }}} 
   } else  { // tminus=TRUE
   if (timex<= cum[0])  { 
	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=0; 
   }
   else if (timex> smax) { 
	   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+(*nx-1)]; 
   }
   else {
   while ((!((timex<=sc2) && (timex>sc1))) && (c>=0)) {
   /* Rprintf(" %lf %lf %lf %ld \n",timex,sc2,sc1,c); */ 
   sc1=cum[c-1];sc2=cum[c];c=c-1; } 
   for(j=1;j<*px;j++) pred[j*(*nxval)+s]=cum[j*(*nx)+c]; 
   }
   }

}
} // }}}

