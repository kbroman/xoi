/**********************************************************************
 *
 * simStahl.c
 *
 * copyright (c) 2006, Karl W Broman
 * 
 * last modified Dec, 2006
 * first written Nov, 2006
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 * 
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * Part of the R/xoi package
 * Contains: simGamma
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GammaS.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>


void simStahl(int *n_sim, double *nu, double *p, double *L,
	      int *nxo, double *loc, int *max_nxo, 
	      int *n_bins4start)
{
  double **Loc, scale;
  double curloc=0.0, u;
  double *startprob, step;
  int i, j, n_nixo;

  /* re-organize loc as a doubly index array */
  Loc = (double **)R_alloc(*n_sim, sizeof(double *));
  Loc[0] = loc;
  for(i=1; i < *n_sim; i++) 
    Loc[i] = Loc[i-1] + *max_nxo;
  
  GetRNGstate();

  if(fabs(*nu - 1.0) < 1e-8) { /* looks like a Poisson model */
    for(i=0; i< *n_sim; i++) {
      R_CheckUserInterrupt(); /* check for ^C */

      nxo[i] = rpois(*L);
      if(nxo[i] > *max_nxo)
	error("Exceeded maximum number of crossovers.");
      
      for(j=0; j < nxo[i]; j++) 
	Loc[i][j] = runif(0.0, *L);
    }
  }
  else {
    scale = 1.0 / (2.0 * *nu * (1.0 - *p)); 

    /* set up starting distribution */
    startprob = (double *)R_alloc(*n_bins4start, sizeof(double));
    step = *L/(double)*n_bins4start;

    startprob[0] = 2.0*(1.0 - *p)*pgamma(((double)i+0.5)*step, *nu, scale, 0, 0)*step;
    for(i=1; i< *n_bins4start; i++) {
      R_CheckUserInterrupt(); /* check for ^C */

      startprob[i] = startprob[i-1] + 
	2.0*(1.0 - *p)*pgamma(((double)i+0.5)*step, *nu, scale, 0, 0)*step;
    }

    for(i=0; i< *n_sim; i++) { 
      R_CheckUserInterrupt(); /* check for ^C */

      nxo[i] = 0;

      /* locations of chiasmata from the gamma model */
      /* shape = nu, rate = 2*nu*(1-p) [scale = 1/{2*nu*(1-p)}] */

      u = unif_rand();
      if( u > startprob[*n_bins4start-1] )
	curloc = *L+1;
      else {
	for(j=0; j< *n_bins4start; j++) {
	  if(u <= startprob[j]) {
	    curloc = ((double)j+0.5)*step;
	    if(unif_rand() < 0.5) {
	      nxo[i] = 1;
	      Loc[i][0] = curloc;
	    }
	    break;
	  }
	}
      }

      if(curloc < *L) {
	while(curloc < *L) {
	  curloc += rgamma(*nu, scale);
	  if(curloc < *L && unif_rand() < 0.5) {
	    if(nxo[i] > *max_nxo)
	      error("Exceeded maximum number of crossovers.");
	    
	    Loc[i][nxo[i]] = curloc;
	    (nxo[i])++;
	  }
	}
      }
      
      /* locations of crossovers from the no interference mechanism */
      if(*p > 0) {
	n_nixo = rpois(*L * *p);
	if(n_nixo > *max_nxo)
	  error("Exceeded maximum number of crossovers.");
	
	for(j=0; j < n_nixo; j++) 
	  Loc[i][nxo[i]+j] = runif(0.0, *L);
	nxo[i] += n_nixo;
      }
    }
  }
    
  /* sort the results */
  for(i=0; i< *n_sim; i++) 
    R_rsort(Loc[i], nxo[i]);

  PutRNGstate();
}

/* end of simStahl.c */
