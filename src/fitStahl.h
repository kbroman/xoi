/**********************************************************************
 *
 * fitStahl.h
 *
 * copyright (c) 2009, Karl W Broman
 * 
 * last modified Jun, 2009
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

void R_stahl_loglik(int *n_ind, int *n_xo, double *xoloc, double *chrlen,
		    int *n_nu, double *nu, double *p, double *loglik,
		    int *max_conv, double *intgr_tol, int *max_subd, 
		    int *min_subd);

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
		  int max_conv, double intgr_tol, int max_subd, int min_subd);

/* end of fitStahl.h */

