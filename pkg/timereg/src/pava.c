/* pava.c: R extension, PAVA (Pool Adjacent Violators Algorithm)  */ 
/* By Bahjat Qaqish */ 
/************************************************************/ 

#include <R.h> 

typedef  double DBL; 

static void wpool (DBL *y, DBL *w, int i, int j) 
/*  
Pool y[i:j] using weights w[i:j] 
*/ 
{   
  int k;
  DBL s0=0, s1=0;

  for (k=i; k<=j; k++) {s1 += y[k]*w[k]; s0 += w[k];}
  s1 /= s0;
  for (k=i; k<=j; k++) y[k] = s1;
}

/*************************************************/

static void wpava  (DBL *y, DBL *w, int *np)
/*
Apply weighted pava to y[0:n-1] using weights w[0:n-1]
*/
{
  int npools, n = *np;

  if (n <= 1) return;
  n--;

  /* keep passing through the array until pooling is not needed */
  do {
    int i = 0;
    npools = 0;
    while (i < n) {
      int k = i;
      /* starting at y[i], find longest non-increasing sequence y[i:k] */
      while (k < n && y[k] >= y[k+1])  k++;
      if (y[i] != y[k]) {wpool(y, w, i, k); npools++;}
      i = k+1;
    }
  } while (npools > 0);
}

/*************************************************/

static void upool (DBL *y, int i, int j)
/*
Pool y[i:j]
*/
{
  int k;
  DBL s=0;

  for (k=i; k<=j; k++) {s += y[k];}
  s /= (j-i+1);
  for (k=i; k<=j; k++) y[k] = s;
}

/*************************************************/

static void  upava  (DBL *y, int *np)
/*
Apply pava to y[0:n-1]
*/
{
  int npools, n = *np;

  if (n <= 1) return;
  n--;

  /* keep passing through the array until pooling is not needed */
  do {
    int i = 0;
    npools = 0;
    while (i < n) {
      int k = i;
      /* starting at y[i], find longest non-increasing sequence y[i:k] */
      while (k < n && y[k] >= y[k+1])  k++;
      if (y[i] != y[k]) {upool(y, i, k); npools++;}
      i = k+1;
    }
  } while (npools > 0);
} 

/*************************************************/

void  pava  (DBL *y, DBL *w, int *np) 
/*
Apply pava to y[0:n-1] using weights w[0:n-1]
Calls an unweighted version if all weights are equal and != 0
Does nothing if all weights are == 0
Calls a weighted version otherwise
*/
{
  int n = *np, i=1;
  DBL w0;

  if (n <= 1) return;

  w0 = w[0];
  while (i < n && w[i] == w0) i++;
  if (i == n) {
    if (w0 == 0.0) return;  /* all weights are == 0 */
    else upava(y, np);      /* unweighted */
  }
  else wpava(y, w, np);     /* weighted   */
} 
