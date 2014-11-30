/**********************************************************************
 *
 * kfunc.h
 *
 * Code to calculate the 1-d version of Ripley's K function
 *
 * Karl W Broman
 * First written 7 April 2006
 * Last modified 28 April 2006
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
             double *area, double *rate, double *tol);

void kfunc(int ngrp, int *n, double **Loc, double *maxl,
           int n_d, double *d, double exclude, double *k,
           double *area, double *rate, int tot, double tol);

/* end of kfunc.h */
