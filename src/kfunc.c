/**********************************************************************
 * 
 * kfunc.c
 *
 * Code to calculate the 1-d version of Ripley's K function
 *
 * Karl W Broman
 * First written 7 April 2006
 * Last modified Nov 2006
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
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "kfunc.h"

/**********************************************************************
 * 
 * ngrp   number of groups of locations
 * 
 * n      vector giving number of locations in each group
 * 
 * loc    vector of locations (assumed sorted); of length sum(n)
 *
 * maxl   maximum length studied in each group
 *
 * n_d    number of points at which to calculate the k function
 *
 * d      values at which to calculate the k function
 *
 * exclude length of region around each point that should be excluded
 *
 * k      on return, the values of the k function
 *
 * area   on return, contains the area covered for each d[i]
 *
 * rate   on return, the estimated rate
 *
 **********************************************************************/

void R_kfunc(int *ngrp, int *n, double *loc, double *maxl, 
	     int *n_d, double *d, double *exclude, double *k, 
	     double *area, double *rate, double *tol)
{
  int i, tot;
  double **Loc;

  /* count total number of points */
  tot=0; 
  for(i=0; i< *ngrp; i++) tot += n[i];

  /* reorganize loc as doubly index array */
  Loc = (double **)R_alloc(*ngrp, sizeof(double *));
  Loc[0] = loc;
  for(i=1; i< *ngrp; i++) Loc[i] = Loc[i-1] + n[i-1];

  kfunc(*ngrp, n, Loc, maxl, *n_d, d, *exclude, k, area, rate, tot, *tol);
}


void kfunc(int ngrp, int *n, double **Loc, double *maxl,
	   int n_d, double *d, double exclude, double *k, 
	   double *area, double *rate, int tot, double tol)
{
  int i, j, k1, k2;
  double temp;

  /* first estimate the rate; exclude points == 0 or == maxl */
  *rate = temp = 0.0;
  for(i=0; i<ngrp; i++) {
    temp += maxl[i];

    for(j=0; j<n[i]; j++) {
      if(fabs(Loc[i][j]) > tol && fabs(Loc[i][j] - maxl[i]) > tol) 
	*rate += 1.0;
    }
  }
  *rate /= temp; /* this is the estimated no points per unit length */

  for(i=0; i<n_d; i++) {
    area[i] = 0.0; /* to contain area covered */
    k[i] = 0.0; /* initially to contain no. points */
    for(j=0; j<ngrp; j++) {

      for(k1=0; k1<n[j]; k1++) {

	if(Loc[j][k1] < d[i]) {
	  if(Loc[j][k1] < exclude) 
	    area[i] += (d[i]-exclude);
	  else 
	    area[i] += (Loc[j][k1]+d[i]-2.0*exclude);
	}

	else if(Loc[j][k1] > maxl[j]-d[i]) {
	  if(Loc[j][k1] > maxl[j]-exclude) 
	    area[i] += (d[i]-exclude);
	  else
	    area[i] += (maxl[j] - Loc[j][k1] + d[i] - 2.0*exclude);
	}

	else  /* entire ball observed */
	  area[i] += 2.0*(d[i]-exclude);

	for(k2=0; k2<n[j]; k2++) {
	  if(k1 != k2 && fabs(Loc[j][k2] - Loc[j][k1]) <= d[i] &&
	     fabs(Loc[j][k2] - Loc[j][k1]) > exclude) {/* within ball */
	    k[i] += 1.0; 
	  }
	}

      }

    } /* end loop over groups */

    /*    Rprintf("%d %d %lf %lf\n", i, (int)k[i], area[i], *rate); */
    k[i] /= (area[i] * *rate);
  }
}
  
/* end of kfunc.c */
