#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h" 
#include "haplosurv.h" 


void genoLogLike(int *oh, // an integer vector
		 int *nphppp, // an integer vector
		 int *nph, // an integer
		 int *np, // an integer
		 double *hf, // a numeric vector
		 double *rho, // a numeric vector
		 double *logLike, // a numeric
		 double *score, // a numeric vector: this must be zero to start with
		 vector *personScoreLL[*np], // a vector of (numeric) vectors (of length *nph-1) /* Comment out for the example */
		 double *d2l){  // a numeric vector: so must this

  int i,j,k,l;
  int c1,c2, c3; // counters to keep track of the position within the
		 // List of Lists
  double withinLogSum;
  double factor;
  double v1, v2;
  int nphm1 = *nph-1;
  double *v = calloc(*np, sizeof(double));
  // v[k] is the denominator for the score vector for person k

  // u[m + nphm1*k] is the numerator from mth element of the score vector,
  // for person k
  double *u = calloc(nphm1*(*np), sizeof(double));

  c1 = 0; // This allows us to loop through the lists of lists, HPIordered
  c2 = 0; // which is stored as the integer vector orderedHaplos 
  c3 = 0; // (which is now called 'oh')
  for(i = 0; i< *np; i++){
    // First pass: calculate the sum within the log term in the log-likelihood
    withinLogSum = 0.0;

    for(j = 0; j< nphppp[i]; j++){
      if(oh[c1] == oh[c1+1]){
	factor = rho[i] * hf[oh[c1]];
      } else {
	factor = 0.0;
      }

      withinLogSum += (hf[oh[c1]] * hf[oh[c1 + 1]] * (1.0 - rho[i])) + factor;
      c1 += 2;
    }

    // update log-likelihood accordingly
    v[i] = withinLogSum;
    *logLike += log(withinLogSum);

    // second pass: use the results of the first pass to 
    // update the score function
    v1 = (1.0-rho[i]);
    v2 = (1.0-rho[i])/v[i];
    for(j = 0; j< nphppp[i]; j++){
      if(oh[c2] != nphm1 &&
	 oh[c2+1] != nphm1){

	if(oh[c2] == oh[c2+1]){
	  //u[i + nphm1*oh[c2]] += rho[i];
	  u[nphm1*i + oh[c2]] += rho[i];
	}
	u[nphm1*i + oh[c2]] += hf[oh[c2+1]]*v1;
	u[nphm1*i + oh[c2+1]] += hf[oh[c2]]*v1;

      } else if (oh[c2] != nphm1 &&
		 oh[c2+1] == nphm1){

	u[nphm1*i + oh[c2]] += hf[oh[c2+1]]*v1;
	for(k = 0; k < nphm1; k++){
	  u[nphm1*i + k] -= hf[oh[c2]]*v1;
	}

      } else if (oh[c2] == nphm1 &&
		 oh[c2+1] != nphm1){

	u[nphm1*i + oh[c2+1]] += hf[oh[c2]]*v1;
	for(k = 0; k < nphm1; k++){
	  u[nphm1*i + k] -= hf[oh[c2+1]]*v1;
	}

      } else {
	for(k = 0; k < nphm1; k++){
	  u[nphm1*i + k] -= (rho[i] + 2.0*hf[nphm1]*v1);
	}
      }
      c2 += 2;
    }

    // Third pass: calculate part of the score function, using the results
    // from the first and second passes
    for(j = 0; j< nphppp[i]; j++){
      if(oh[c3] != nphm1 &&
	 oh[c3+1] != nphm1){
	  
	d2l[oh[c3]*(nphm1)+oh[c3+1]] += v2;
	d2l[oh[c3+1]*(nphm1)+oh[c3]] += v2;
	    
      } else if (oh[c3] != nphm1 &&
		 oh[c3+1] == nphm1){

	for(k = 0; k < nphm1; k++){
	  d2l[oh[c3]*(nphm1)+k] -= v2;
	  d2l[k*(nphm1)+oh[c3]] -= v2;
	}

      } else if (oh[c3] == nphm1 &&
		 oh[c3+1] != nphm1){

	for(k = 0; k < nphm1; k++){
	  d2l[oh[c3+1]*(nphm1)+k] -= v2;
	  d2l[k*(nphm1)+oh[c3+1]] -= v2;
	}

      } else {
	for(k = 0; k < nphm1; k++){
	  for(l = 0; l < nphm1; l++){
	    d2l[k*(nphm1)+l] += 2.0*v2;
	  }
	}
      }
      c3 += 2;
    }

  }

  // Add the final touch to the D2L matrix

  for(i = 0; i< *np; i++){
    v2 = v[i]*v[i];
    for(k = 0; k < nphm1; k++){
      score[k] += u[nphm1*i + k]/v[i];
      VE(personScoreLL[i],k) += u[nphm1*i + k]/v[i]; /* Comment out for the example */
      for(l = 0; l < nphm1; l++){
	d2l[k*(nphm1)+l] -= u[k + nphm1*i]*u[l + nphm1*i]/v2;
      }
    }
  }

  free(u);
  free(v);

}


void genoLogLikeHp(int *oh, // stands for orderedHaplos, an integer
			  // vector, giving the indices of the
			  // possible haplotype pairs for each
			  // individual, this is in fact a list of
			  // lists convereted to an integer vector;
			  // each element of the outer list
			  // corresponds to a given individual, each
			  // element of the inner list corresponds to
			  // a pair of possible haplotypes for this
			  // individual, and is an integer vector of
			  // length two, giving the indices of this
			  // pair of haplotypes
		   int *nphppp, // stands for number of possible
			      // haplotypes pairs per person, it is an
			      // integer vector giving the lengths of
			      // the inner lists of orderedHaplos; it
			      // is used to keep track of the position
			      // when looping over oh[i]
		   int *nph, // stands for number of possible
			   // haplotypes, an integer giving; it is one
			   // more than the number of haplotype
			   // parameters
		   int *amountToReturn, // if 0, only calculates the
					// log-likelihood. if 1, also
					// calculates the score
					// vector. if 2, also
					// calculates the individual
					// score vectors. if 3, does
					// not calculate the
					// individual score vector,
					// but calculates the score
					// vector and the hessian matrix.
		   int *np, // stands for the number of people, an
			  // integer
		   double *hp, // stands for haplotype parameters, a
			     // numeric vector of length *nph-1, that
			     // is used to compute the haplotype
			     // frequency vector
		   double *rho, // is a vector of inbreeding
			      // coefficients, a numeric vector of
			      // length *np, one entry per person
		   double *logLike, // the log likelihood, a double
		   double *score, // the score vector, with respect to
				// the haplotype parameters, it is a
				// numeric vector of length *nhp-1:
				// this must be zero to start with
		   vector *personScoreLL[*np], // person-specific  /* Comment out for the example */
					     // derivatives of the
					     // score vector with
					     // respect to the
					     // haplotype parameters,
					     // not the hapltype
					     // frequencies: this is a
					     // vector of (numeric)
					     // vectors (of length
					     // *nph-1)
		   double *d2l){  // the matrix of second derivatives of
				// the log-likelihood with respect to
				// the haplotype parameters a numeric
				// vector: this must be zero to start
				// with

  int i,j,k,l;
  int c1,c2, c3; // counters to keep track of the position within the
		 // List of Lists
  double withinLogSum;
  double factor;
  double v1, v2;
  int nphm1 = *nph-1; // stands for the number of possible haplotypes
		      // minus 1
  double *v = calloc(*np, sizeof(double));
  // v[k] is the denominator for the score vector for person k

  // u[m + nphm1*k] is the numerator from mth element of the score vector,
  // for person k
  double *u = calloc(nphm1*(*np), sizeof(double));

  // The score with respect to the _haplotype frequencies_, rather
  // than the haplotype parameters themselves
  double *scoreWRThf = calloc(nphm1, sizeof(double));

  // The hessian with respect to the _haplotype frequencies_, rather
  // than the haplotype parameters themselves
  double *d2lWRThf = calloc(nphm1*nphm1, sizeof(double));

  // used for preliminary calculations of the hessian matrix
  double D2lpipi = 0.0;
  double scorepi = 0.0;

  // Calculate haplotype frequencies based on haplotype parameters
  double *hf = calloc(*nph, sizeof(double));
  factor = 1.0;
  for(i = 0; i< nphm1; i++){
    factor += exp(hp[i]);
  }
  hf[nphm1] = 1.0/factor;
  for(i = 0; i< nphm1; i++){
    hf[i] = exp(hp[i])/factor;
  }

  
  c1 = 0; // This allows us to loop through the lists of lists, HPIordered
  c2 = 0; // which is stored as the integer vector orderedHaplos 
  c3 = 0; // (which is now called 'oh')
  for(i = 0; i< *np; i++){
    // First pass: calculate the sum within the log term in the log-likelihood
    withinLogSum = 0.0;


    for(j = 0; j< nphppp[i]; j++){
      if(oh[c1] == oh[c1+1]){
	factor = rho[i] * hf[oh[c1]];
      } else {
	factor = 0.0;
      }

      withinLogSum += (hf[oh[c1]] * hf[oh[c1 + 1]] * (1.0 - rho[i])) + factor;
      c1 += 2;
    }

    // update log-likelihood accordingly
    v[i] = withinLogSum;
    *logLike += log(withinLogSum);
    
    if(*amountToReturn>0){
      // second pass: use the results of the first pass to 
      // update the score function
      v1 = (1.0-rho[i]);
      for(j = 0; j< nphppp[i]; j++){
	if(oh[c2] != nphm1 &&
	   oh[c2+1] != nphm1){

	  if(oh[c2] == oh[c2+1]){
	    //u[i + nphm1*oh[c2]] += rho[i];
	    u[nphm1*i + oh[c2]] += rho[i];
	  }
	  u[nphm1*i + oh[c2]] += hf[oh[c2+1]]*v1;
	  u[nphm1*i + oh[c2+1]] += hf[oh[c2]]*v1;

	} else if (oh[c2] != nphm1 &&
		   oh[c2+1] == nphm1){

	  u[nphm1*i + oh[c2]] += hf[nphm1]*v1;
	  for(k = 0; k < nphm1; k++){
	    u[nphm1*i + k] -= hf[oh[c2]]*v1;
	  }

	} else if (oh[c2] == nphm1 &&
		   oh[c2+1] != nphm1){

	  u[nphm1*i + oh[c2+1]] += hf[nphm1]*v1;
	  for(k = 0; k < nphm1; k++){
	    u[nphm1*i + k] -= hf[oh[c2+1]]*v1;
	  }

	} else {
	  for(k = 0; k < nphm1; k++){
	    u[nphm1*i + k] -= (rho[i] + 2.0*hf[nphm1]*v1);
	  }
	}
	c2 += 2;
      }
    }
    if(*amountToReturn>1){
      // Third pass: calculate part of the score function, using the results
      // from the first and second passes
      v2 = (1.0-rho[i])/v[i];
      for(j = 0; j< nphppp[i]; j++){
	if(oh[c3] != nphm1 &&
	   oh[c3+1] != nphm1){
	  
	  d2lWRThf[oh[c3]*(nphm1)+oh[c3+1]] += v2;
	  d2lWRThf[oh[c3+1]*(nphm1)+oh[c3]] += v2;
	    
	} else if (oh[c3] != nphm1 &&
		   oh[c3+1] == nphm1){

	  for(k = 0; k < nphm1; k++){
	    d2lWRThf[oh[c3]*(nphm1)+k] -= v2;
	    d2lWRThf[k*(nphm1)+oh[c3]] -= v2;
	  }

	} else if (oh[c3] == nphm1 &&
		   oh[c3+1] != nphm1){

	  for(k = 0; k < nphm1; k++){
	    d2lWRThf[oh[c3+1]*(nphm1)+k] -= v2;
	    d2lWRThf[k*(nphm1)+oh[c3+1]] -= v2;
	  }

	} else {
	  for(k = 0; k < nphm1; k++){
	    for(l = 0; l < nphm1; l++){
	      d2lWRThf[k*(nphm1)+l] += 2.0*v2;
	    }
	  }
	}
	c3 += 2;
      }
    }
  }

  if(*amountToReturn==2){ // compute the score vector and the subject-specific score vectors, but not the hessian
    for(i = 0; i< *np; i++){
      v2 = v[i]*v[i];
      // calculate the person-specific score vector, with respect to
      // the haplotype parameters, rather than the haplotype frequencies.
      VE(personScoreLL[i],0) = -1.0*hf[0]*u[nphm1*i]/v[i]; /* Comment out for the example */
      for(k = 1; k < nphm1; k++){
	VE(personScoreLL[i],0) -= hf[k]*u[nphm1*i + k]/v[i]; /* Comment out for the example */
      }
      for(k = 1; k < nphm1; k++){
	VE(personScoreLL[i],k) = VE(personScoreLL[i],0); /* Comment out for the example */
      }
      for(k = 0; k < nphm1; k++){
	scoreWRThf[k] += u[nphm1*i + k]/v[i];
	VE(personScoreLL[i],k) += u[nphm1*i + k]/v[i]; /* Comment out for the example */
	VE(personScoreLL[i],k) *= hf[k]; /* Comment out for the example */
      }
    }
  } else if(*amountToReturn==1){ // compute only the score vector, but not the hessian or the subject-specific score vectors
    for(i = 0; i< *np; i++){
      for(k = 0; k < nphm1; k++){
	scoreWRThf[k] += u[nphm1*i + k]/v[i];
      }
    }
  } else if(*amountToReturn==3){ // compute the score vector and the hessian matrix, but not the subject-specific score vectors
    for(i = 0; i< *np; i++){
      v2 = v[i]*v[i];
      for(k = 0; k < nphm1; k++){
	scoreWRThf[k] += u[nphm1*i + k]/v[i];
	for(l = 0; l < nphm1; l++){
	  d2lWRThf[k*(nphm1)+l] -= u[k + nphm1*i]*u[l + nphm1*i]/v2;   // Add the final touch to the d2lWRThr matrix
	}
      }
    }
  }

  free(u);
  free(v);

  if(*amountToReturn==0){
    free(scoreWRThf);
    free(d2lWRThf);
    free(hf);
    return;
  }

  // calculate the score with respect to the haplotype parameters,
  // rather than the haplotype frequencies
  for(i = 0; i < nphm1; i++){
    scorepi += scoreWRThf[i]*hf[i];
  }

  for(i = 0; i < nphm1; i++){
    score[i] = hf[i]*(scoreWRThf[i]-scorepi);
  }

  if(*amountToReturn==1 || *amountToReturn==2){
    free(scoreWRThf);
    free(d2lWRThf);
    free(hf);
    return;
  }

  // calculate the product of two different configurations of the
  // hessian and the frequency vector...
  double *D2lpi = calloc(nphm1, sizeof(double));
  for(i = 0; i < nphm1; i++){
    for(j = 0; j < nphm1; j++){
      D2lpipi += d2lWRThf[i*nphm1+j]*hf[i]*hf[j];
      D2lpi[i] += d2lWRThf[i*nphm1+j]*hf[j];
    }
  }  
  // ... then use this to calculate the hessian matrix with respect to the
  // haplotype parameters, rather than with respect to the haplotype
  // frequencies
  for(i = 0; i < nphm1; i++){
    for(j = 0; j < nphm1; j++){
      d2l[j*nphm1+i] = D2lpipi - D2lpi[i] - D2lpi[j] + d2lWRThf[j*nphm1+i]
	+ 2.0*scorepi - scoreWRThf[i]- scoreWRThf[j];
      d2l[j*nphm1+i] *= hf[i]*hf[j];
    }
    d2l[i*nphm1+i] += hf[i]*(scoreWRThf[i] - scorepi);
  }

  free(scoreWRThf);
  free(d2lWRThf);
  free(D2lpi);
  free(hf);

}

void genoLogLikeHpCall(int *oh, // stands for orderedHaplos, an integer
		       int *nphppp, // stands for number of possible
		       int *nph, // stands for number of possible
		       int *amountToReturn, // if 0, only calculates the
		       int *np, // stands for the number of people, an
		       double *hp, // stands for haplotype parameters, a
		       double *rho, // is a vector of inbreeding
		       double *logLike, // the log likelihood, a double
		       double *score, // the score vector, with respect to
		       double *personScoreLLdouble, // person-specific  /* Comment out for the example */
		       double *d2l){  // the matrix of second derivatives of

  vector *personScoreLL[*np];
  int i,j;
  int nphm1 = *nph - 1;
  
  for(i = 0; i < *np; i++){
    malloc_vec(nphm1,personScoreLL[i]);
  }

  genoLogLikeHp(oh, // stands for orderedHaplos, an integer
		nphppp, // stands for number of possible
		nph, // stands for number of possible
		amountToReturn, // if 0, only calculates the
		np, // stands for the number of people, an
		hp, // stands for haplotype parameters, a
		rho, // is a vector of inbreeding
		logLike, // the log likelihood, a double
		score, // the score vector, with respect to
		personScoreLL, // person-specific  /* Comment out for the example */
		d2l);  // the matrix of second derivatives of

  for(i = 0; i < *np; i++){
    for(j = 0; j < nphm1; j++){
      personScoreLLdouble[i*nphm1+j] = VE(personScoreLL[i],j);
    }
  }
  for(i = 0; i < *np; i++){
    free_vec(personScoreLL[i]);
  }

  return;

}

void genoLogLikeCall(int *oh, // stands for orderedHaplos, an integer
		     int *nphppp, // stands for number of possible
		     int *nph, // stands for number of possible
		     int *amountToReturn, // if 0, only calculates the
		     int *np, // stands for the number of people, an
		     double *hf, // stands for haplotype parameters, a
		     double *rho, // is a vector of inbreeding
		     double *logLike, // the log likelihood, a double
		     double *score, // the score vector, with respect to
		     double *personScoreLLdouble, // person-specific  /* Comment out for the example */
		     double *d2l){  // the matrix of second derivatives of

  vector *personScoreLL[*np];
  int i,j;
  int nphm1 = *nph - 1;
  
  for(i = 0; i < *np; i++){
    malloc_vec(nphm1,personScoreLL[i]);
  }

  genoLogLike(oh, // stands for orderedHaplos, an integer
	      nphppp, // stands for number of possible
	      nph, // stands for number of possible
	      np, // stands for the number of people, an
	      hf, // stands for haplotype parameters, a
	      rho, // is a vector of inbreeding
	      logLike, // the log likelihood, a double
	      score, // the score vector, with respect to
	      personScoreLL, // person-specific  /* Comment out for the example */
	      d2l);  // the matrix of second derivatives of

  for(i = 0; i < *np; i++){
    for(j = 0; j < nphm1; j++){
      personScoreLLdouble[i*nphm1+j] = VE(personScoreLL[i],j);
    }
  }
  for(i = 0; i < *np; i++){
    free_vec(personScoreLL[i]);
  }

  return;

}


