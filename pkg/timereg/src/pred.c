//#include <stdio.h>
#include <math.h>
                 
void Cpred(cum,nx,px,xval,nxval,pred,tminus)
double *cum,*xval;
int *nxval,*nx,*px,*tminus,*pred;
{
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
}
