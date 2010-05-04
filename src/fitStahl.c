/**********************************************************************
 *
 * fitStahl.c
 *
 * copyright (c) 2009-2010, Karl W Broman
 * 
 * last modified May, 2010
 * first written Jun, 2009
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
 * Contains: R_stahl_loglik, stahl_loglik
 *
 * Functions to calculate and maximize the likelihood for the Stahl
 * model, given data on the crossover locations
 *
 **********************************************************************/

#include <math.h>
#include "GammaS.h"
#include "fitStahl.h"
#include "stahl_util.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>

void R_stahl_loglik(int *n_ind, int *n_xo, double *xoloc, double *chrlen,
		    int *n_nu, double *nu, double *p, double *loglik,
		    int *max_conv, double *intgr_tol, int *max_subd, 
		    int *min_subd)
{
  int i;
  double **XOloc;

  /* reorganize XOloc as 2-d array */
  XOloc = (double **)R_alloc(*n_ind, sizeof(double *));
  XOloc[0] = xoloc;
  for(i=1; i<*n_ind; i++) {
    XOloc[i] = XOloc[i-1] + n_xo[i-1];
  }

  stahl_loglik(*n_ind, n_xo, XOloc, chrlen, *n_nu, nu, p,
	       loglik, *max_conv, *intgr_tol, *max_subd, *min_subd);
}

/**********************************************************************
 * stahl_loglik: calculate log likelihood for the Stahl model
 * 
 * n_ind  = no. individuals
 * n_xo   = vector with no. crossovers for each individual
 * xoloc  = crossover locations xoloc[i][j] = XO j in individual i
 * chrlen = vector of chromosome lengths
 * n_nu   = length of nu and p vectors
 * nu     = interference parameter at which loglik is to be calculated
 * p      = proportion of crossovers from the no-interference pathway
 * loglik = (on output) the log likelihood (of length n_nu)
 * max_conv  = maximum number of convolutions
 * intgr_tol = tolerance for integration
 * max_subd  = maximum number of subdivisions in integration
 * min_subd  = minimum number of subdivisions in integration
 **********************************************************************/

void stahl_loglik(int n_ind, int *n_xo, double **xoloc, double *chrlen,
		  int n_nu, double *nu, double *p, double *loglik,
		  int max_conv, double intgr_tol, int max_subd, int min_subd)
{
  int i, j, k, s, n_pat, max_n_xo;
  double *y, *z, cury, curz, *patll, *indll;
  int n_y, n_z;
  struct gamma_data info;
  
  max_n_xo=1;
  for(i=0; i<n_ind; i++) 
    if(max_n_xo < n_xo[i]) max_n_xo = n_xo[i];
  y = (double *)R_alloc(max_n_xo, sizeof(double));
  z = (double *)R_alloc(max_n_xo, sizeof(double));
  patll = (double *)R_alloc(n_nu, sizeof(double));
  indll = (double *)R_alloc(n_nu, sizeof(double));

  info.max_conv = max_conv;

  setup_integr_par(intgr_tol, max_subd, min_subd, &(info.integr_info));

  for(i=0; i<n_nu; i++) loglik[i] = 0.0;

  for(j=0; j<n_ind; j++) {
    R_CheckUserInterrupt(); /* check for ^C */

    if(n_xo[j] == 0) { /* no crossovers */
      for(i=0; i<n_nu; i++) {
	loglik[i] += log(oneminus_Gstar_stahl(chrlen[j], 1.0, p[i], max_conv, info.integr_info)) +
	  log(oneminus_Gstar_stahl(chrlen[j], nu[i], 1.0-p[i], max_conv, info.integr_info));
      }
    }
    else {
      for(i=0; i<n_nu; i++) indll[i] = 0.0;

      n_pat = 1 << n_xo[j]; /* no. patterns */
      for(k=0; k<n_pat; k++) { /* loop over patterns */

	cury = curz = 0;
	n_y = n_z = 0;
	for(s=0; s<n_xo[j]; s++) {
	  if(k & (1 << s)) { /* is XO from NI pathway? */
	    y[n_y] = xoloc[j][s] - cury;
	    cury = xoloc[j][s];
	    n_y++;
	  }
	  else {
	    z[n_z] = xoloc[j][s] - curz;
	    curz = xoloc[j][s];
	    n_z++;
	  }
	}

	/* no interference pathway */
	if(n_y==0) { /* no crossovers */
	  for(i=0; i<n_nu; i++) 
	    patll[i] = log(oneminus_Gstar_stahl(chrlen[j], 1.0, p[i], 
						max_conv, info.integr_info));
	}
	else if(n_y==1) { /* exactly one crossover */
	  for(i=0; i<n_nu; i++) 
	    patll[i] = log(gstar_stahl(y[0], 1.0, p[i], max_conv)) +
	      log(oneminus_Fstar_stahl(chrlen[j]-y[0], 1.0, p[i], max_conv));
	}
	else { /* more than one crossover */
	  for(i=0; i<n_nu; i++) {
	    patll[i] = log(gstar_stahl(y[0], 1.0, p[i], max_conv));
	    for(s=1; s<n_y; s++)
	      patll[i] += log(fstar_stahl(y[s], 1.0, p[i], max_conv));
	    patll[i] += log(oneminus_Fstar_stahl(chrlen[j]-cury, 1.0, p[i], max_conv));
	  }
	}

	/* interference pathway */
	if(n_z==0) { /* no crossovers */
	  	  for(i=0; i<n_nu; i++) 
	  	    patll[i] += log(oneminus_Gstar_stahl(chrlen[j], nu[i], 1.0-p[i], 
	  						 max_conv, info.integr_info));
	}
	else if(n_z==1) { /* exactly one crossover */
	  for(i=0; i<n_nu; i++) 
	    patll[i] += log(gstar_stahl(z[0], nu[i], 1.0 - p[i], max_conv)) +
	      log(oneminus_Fstar_stahl(chrlen[j]-z[0], nu[i], 1.0 - p[i], max_conv));
	}
	else { /* more than one crossover */
	  for(i=0; i<n_nu; i++) {
	    patll[i] += log(gstar_stahl(z[0], nu[i], 1.0 - p[i], max_conv));
	    for(s=1; s<n_z; s++)
	      patll[i] += log(fstar_stahl(z[s], nu[i], 1.0 - p[i], max_conv));
	    patll[i] += log(oneminus_Fstar_stahl(chrlen[j]-curz, nu[i], 1.0 - p[i], max_conv));
	  }
	}
	
	if(k==0) 
	  for(i=0; i<n_nu; i++) indll[i] = patll[i];
	else
	  for(i=0; i<n_nu; i++) indll[i] = addlog_stahl(indll[i], patll[i]);

      } /* loop over patterns */
      for(i=0; i<n_nu; i++) loglik[i] += indll[i];

    } /* there are crossovers */
  } /* loop over individuals */

}



/* end of fitStahl.c */

