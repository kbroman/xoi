/**********************************************************************
 *
 * GammaS.c
 *
 * copyright (c) 1998-2012, Karl W Broman
 * 
 * last modified Oct, 2012
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
 * Contains: GammaS, GammaMax, GammaMax2 
 *
 * Functions to calculate and maximize the likelihood for the Gamma
 * model, given data on the inter-crossover distances
 *
 **********************************************************************/

#include <math.h>
#include "GammaS.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "zeroin.h"

/**********************************************
 * GammaS
 *
 * This function is to be called from R
 *
 * It calculates the likelihood for the gamma
 * renewal model for each of a set of different
 * parameter values.
 *
 * n_length = number of inter-XO lengths
 * length   = the lengths (in Morgans)
 * type     = censor type (0=uncensored; 1=right-censored; 2=initial point
 *                         3=whole chromosome)
 * n_nu     = number of values at which the log likelihood will be calculated
 * nu       = vector of values at which the log likelihood will be calculated
 * loglik   = on exit, the log likelihood at each nu
 *
 * max.conv = maximum number of convolutions
 * center   = if 1, re-scale the log likelihoods to have 0 at the maximum
 *
 **********************************************/
void GammaS(int *n_length, double *length, int *type, int *n_nu,
	    double *nu, double *loglik, int *max_conv, int *center,
	    double *integr_tol, int *maxsubd, int *minsubd)
{
  int i;
  double mx=0.0;
  struct gamma_data info;

  info.n_length = *n_length;
  info.length = length;
  info.type = type;
  info.max_conv = *max_conv;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  for(i=0; i < *n_nu; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    loglik[i] = -calcLL(nu[i], &info);

    if(i==0) mx = loglik[i];
    else { if(mx < loglik[i]) mx=loglik[i]; }
  }

  if(*center) for(i=0; i < *n_nu; i++) { loglik[i] -= mx; }
}


/**********************************************************************
 *
 * GammaMax
 *
 * n_length = number of inter-XO lengths
 * length   = the lengths (in Morgans)
 * type     = censor type (0=uncensored; 1=right-censored; 2=initial point
 *                         3=whole chromosome)
 * low      = lower value for nu to look at
 * high     = upper value for nu to look at
 * nu       = on entry, estimated value.  on exit, the MLE
 * loglik   = on exit, the log likelihood at the MLE
 * max.conv = maximum number of convolutions
 * tol      = tolerance for convergence
 **********************************************************************/

void GammaMax(int *n_length, double *length, int *type,
	      double *low, double *high, double *nu,
	      double *loglik, int *max_conv, double *tol, 
	      double *integr_tol, int *maxsubd, int *minsubd)
{
  struct gamma_data info;

  info.max_conv = *max_conv;
  info.n_length = *n_length;
  info.type = type;
  info.length=length;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  *nu = Brent_fmin(*low, *high, (double (*)(double, void *))calcLL, 
		   (void *)(&info), *tol);

  *loglik = -calcLL(*nu, &info);
}


/**********************************************************************
 *
 * GammaSE
 *
 * n_length, length, type, max_conv as in GammaS
 * nu = MLE
 * se = on exit, the estimated SE
 * secderiv = on exit, the second derivative at the MLE
 * h = initial nu difference for estimating second derivative
 * hstep = factor by which to reduce h
 * tol = tolerance for convergence
 **********************************************************************/

void GammaSE(int *n_length, double *length, int *type,
	     double *nu, double *se, double *secderiv,
	     int *max_conv, double *h, double *hstep, 
	     double *tol, int *maxit, 
	     double *integr_tol, int *maxsubd, int *minsubd)
{
  double f, fph, fmh;
  double cur_secderiv;
  int i, j;
  struct gamma_data info;

  info.max_conv = *max_conv;
  info.n_length = *n_length;
  info.type = type;
  info.length = length;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  f = -calcLL(*nu, &info);
  fph = -calcLL(*nu + *h, &info);
  fmh = -calcLL(*nu - *h, &info);

  cur_secderiv = (fph - 2.0*f + fmh)/(*h * *h);

  for(i=0; i< *maxit; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    *h /= *hstep;

    fph = -calcLL(*nu + *h, &info);
    fmh = -calcLL(*nu - *h, &info);
    for(j=0; j < *n_length; j++)

    *secderiv = (fph - 2.0*f + fmh)/(*h * *h);

    if(fabs(*secderiv-cur_secderiv) < *tol) break;
    cur_secderiv = *secderiv;
  }

  *se = sqrt( -1.0 / *secderiv);
}

/**********************************************************************
 *
 * GammaInterval
 *
 * Calculate a likelihood support interval 
 * (drop is the amount to drop in terms of natural logarithm)
 *
 *
 **********************************************************************/

void GammaInterval(int *n_length, double *length, int *type,
		   double *low, double *high, 
		   double *nu, double *interval,
		   double *interval_level, double *drop,
		   int *max_conv, double *tol, int *maxit, 
		   double *integr_tol, int *maxsubd, int *minsubd)
{
  double temptol;
  int tempmaxit;
  struct gamma_data info;

  /* maximum */
  info.max_conv = *max_conv;
  info.n_length = *n_length;
  info.type = type;
  info.length = length;
  info.drop = *drop;

  setup_integr_par(*integr_tol, *maxsubd, *minsubd, &(info.integr_info));

  R_CheckUserInterrupt(); /* check for ^C */

  info.maxloglik = -calcLL(*nu, &info);


  R_CheckUserInterrupt(); /* check for ^C */

  /* lower limit */
  temptol = *tol;
  tempmaxit = *maxit;
  interval[0] = Rxoi_zeroin(*low, *nu, (double (*)(double, void *))calcLLmdrop, 
                            (void *)(&info), &temptol, &tempmaxit);
  interval_level[0] = -calcLL(interval[0], &info);

  R_CheckUserInterrupt(); /* check for ^C */

  /* upper limit */
  temptol = *tol;
  tempmaxit = *maxit;
  interval[1] = Rxoi_zeroin(*nu, *high, (double (*)(double, void *))calcLLmdrop, 
                            (void *)(&info), &temptol, &tempmaxit);
  interval_level[1] = -calcLL(interval[1], &info);
}

/* end of GammaS.c */

