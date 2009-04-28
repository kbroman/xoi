/**********************************************************************
 *
 * GammaDensities.c
 *
 * copyright (c) 1998-2007, Karl W Broman
 * 
 * last modified Apr, 2007
 * first written ~Oct, 1998
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
 * Contains: location_given_one, first_given_two, distance_given_two
 *           xoprob, xoprob_sub, xoprob_subsub, ioden, firstden,
 *           GammaCoincidence, StahlCoincidence
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GammaS.h"

/* dist'n of location of XO given exactly one XO */
void location_given_one(double *nu, double *x, double *y, int *n, double *L, 
			int *max_conv, double *integr_tol, int *maxsubd,
			int *minsubd) 
{
  int i;
  double denom;
  struct gamma_data info;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  info.nu = *nu;
  info.length = L;
  info.max_conv = *max_conv;

  denom = my_integrate( (void (*)(double *, int, void *))lg1_sub, 
  			(void *)(&info), 0, *L, info.integr_info);

  for(i=0; i<*n; i++) y[i] = x[i];

  lg1_sub(y, *n, &info);
  for(i=0; i<*n; i++) y[i] /= denom;
}

void lg1_sub(double *x, int n, struct gamma_data *info)
{
  int i;

  for(i=0; i<n; i++) 
    x[i] = exp(ll((*info).nu, 1, x[i],                     (*info).max_conv, (*info).integr_info) +
	       ll((*info).nu, 1, (*info).length[0] - x[i], (*info).max_conv, (*info).integr_info));
}


/* dist'n of location of first XO given exactly two XOs */
void first_given_two(double *nu, double *L, double *x,
		     double *y, int *n, int *max_conv, 
		     double *integr_tol, int *maxsubd, int *minsubd)
{
  double temp[2], denom;
  struct gamma_data info;
  int i;

  temp[0] = *L;

  info.max_conv = *max_conv;
  info.nu = *nu;
  info.length = temp;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  denom = my_integrate((void (*)(double *, int, void *))xoprob_sub, 
  		       (void *)(&info), 0.0, *L, info.integr_info);

  for(i=0; i < *n; i++) y[i] = x[i];
  xoprob_sub(y, *n, &info);
  for(i=0; i < *n; i++) y[i] /= denom;
}



/* joint dist'n of XO locations given exactly two XOs */
void joint_given_two(double *nu, double *L, double *x, double *y,
		     double *z, int *n, int *max_conv, 
		     double *integr_tol, int *maxsubd, int *minsubd)
{
  double temp[2], denom;
  struct gamma_data info;
  int i;

  temp[0] = *L;

  info.max_conv = *max_conv;
  info.nu = *nu;
  info.length = temp;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  denom = my_integrate((void (*)(double *, int, void *))xoprob_sub, 
  		       (void *)(&info), 0.0, *L, info.integr_info);

  for(i=0; i < *n; i++) 
    z[i] = exp(ll(*nu, 1, x[i], *max_conv, info.integr_info) +
	       ll(*nu, 0, y[i]-x[i], *max_conv, info.integr_info) +
	       ll(*nu, 1, *L - y[i], *max_conv, info.integr_info))/denom;
}





/* dist'n of distance between XOs given two XOs */
void distance_given_two(double *nu, double *L, double *x,
			double *y, int *n, int *max_conv,
			double *integr_tol, int *maxsubd,
			int *minsubd)
{
  double temp[2], denom;
  struct gamma_data info;
  int i;

  temp[0] = *L;

  info.max_conv = *max_conv;
  info.nu = *nu;
  info.length = temp;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  denom = my_integrate((void (*)(double *, int, void *))xoprob_sub_b, 
  		       (void *)(&info), 0.0, *L, info.integr_info);

  for(i=0; i < *n; i++) y[i] = x[i];
  xoprob_sub_b(y, *n, &info);
  for(i=0; i < *n; i++) y[i] /= denom;
}


void distance_given_two_sub(double *x, int n, struct gamma_data *info)
{
  int i;

  for(i=0; i<n; i++) 
    x[i] = exp(ll((*info).nu, 2, x[i], (*info).max_conv, (*info).integr_info) +
	       ll((*info).nu, 1, (*info).length[0]-x[i], (*info).max_conv, (*info).integr_info));
}


/* calculate probabilities of 0, 1, 2 and >2 XOs */
void xoprob(double *nu, double *L, double *pr, int *max_conv,
	    double *integr_tol, int *maxsubd, int *minsubd)
{
  double temp[2];
  struct gamma_data info;

  temp[0] = *L;

  info.max_conv = *max_conv;
  info.nu = *nu;
  info.length = temp;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  /* pr of 0 XOs */
  pr[0] = exp(ll(*nu, 3, *L, *max_conv, info.integr_info));

  /* pr of 1 XO */
  pr[1] = my_integrate((void (*)(double *, int, void *))lg1_sub,
		       (void *)(&info), 0.0, *L, info.integr_info);

  /* pr of 2 XO */
  pr[2] = my_integrate((void (*)(double *, int, void *))xoprob_sub, 
		       (void *)(&info), 0.0, *L, info.integr_info);

  pr[3] = 1.0 - pr[0] - pr[1] - pr[2];

}

/* calculates un-scaled joint density of locations of two XO given that there are two */
void xoprob_subsub(double *x, int n, struct gamma_data *info)
{
  /* (*info).length[0] is chr length */
  /* (*info).length[1] is location of first XO */
  int i;
  double x0, L;

  L = (*info).length[0];
  x0 = (*info).length[1];

  for(i=0; i<n; i++) 
    x[i] = exp(ll((*info).nu, 1, x0, (*info).max_conv, (*info).integr_info) +
	       ll((*info).nu, 0, x[i], (*info).max_conv, (*info).integr_info) +
	       ll((*info).nu, 1, L-x0-x[i], (*info).max_conv, (*info).integr_info));
}

/* calculates un-scaled joint density of locations of two XO given that there are two */
void xoprob_subsub_b(double *x, int n, struct gamma_data *info)
{
  /* (*info).length[0] is chr length */
  /* (*info).length[1] is distance between the first and second XOs */
  int i;
  double x1, L;

  L = (*info).length[0];
  x1 = (*info).length[1];

  for(i=0; i<n; i++) 
    x[i] = exp(ll((*info).nu, 1, x[i], (*info).max_conv, (*info).integr_info) +
	       ll((*info).nu, 0, x1, (*info).max_conv, (*info).integr_info) +
	       ll((*info).nu, 1, L-x1-x[i], (*info).max_conv, (*info).integr_info));
}


/* calculates distribution of location of first XO given that there are two */
void xoprob_sub(double *x, int n, struct gamma_data *info)
{
  int i;

  for(i=0; i<n; i++) {
    (*info).length[1] = x[i];
    x[i] = my_integrate((void (*)(double *, int, void *))xoprob_subsub,
			(void *)info, 0.0, (*info).length[0]-x[i], (*info).integr_info);
  }
}

/* calculates distribution of location of first XO given that there are two */
void xoprob_sub_b(double *x, int n, struct gamma_data *info)
{
  int i;

  for(i=0; i<n; i++) {
    (*info).length[1] = x[i];
    x[i] = my_integrate((void (*)(double *, int, void *))xoprob_subsub_b,
			(void *)info, 0.0, (*info).length[0]-x[i], (*info).integr_info);
  }
}

/* inter-crossover density */
void ioden(double *nu, double *x, double *y, int *n, int *max_conv)
{
  int i;
  struct integr_data temp;

  for(i=0; i< *n; i++)
    y[i] = exp(ll(*nu, 0, x[i], *max_conv, temp));
}

/* inter-crossover density */
void firstden(double *nu, double *x, double *y, int *n, int *max_conv)
{
  int i;
  struct integr_data temp;

  for(i=0; i< *n; i++)
    y[i] = exp(ll(*nu, 1, x[i], *max_conv, temp));
}


void GammaCoincidence(double *nu, double *x, double *y, int *n,
		      int *max_conv)
{
  int i, j;

  for(i=0; i<*n; i++) {
    for(j=1; j<*max_conv; j++) 
      y[i] += mydgamma(x[i], (double)j * *nu, 2.0 * *nu);
    y[i] /= 2.0;
  }
}


void StahlCoincidence(double *nu, double *p, double *x, double *y, int *n,
		      int *max_conv)
{
  int i, j;

  for(i=0; i<*n; i++) {
    y[i] = 0.0;
    for(j=1; j<*max_conv; j++) 
      y[i] += mydgamma(x[i], (double)j * *nu, 2.0 * (1.0 - *p) * *nu);
    y[i] = (y[i]/2.0) + *p;
  }
}

/* end of GammaDensities.c */
