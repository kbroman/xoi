/**********************************************************************
 *
 * recrate.h
 *
 * copyright (c) 2008, Karl W Broman
 *
 * last modified Aug, 2008
 * first written Aug, 2008
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
 * Contains: R_est_recrate, recrate
 *
 **********************************************************************/

/* R wrapper */
void R_est_recrate(int *n_mar, double *gen, double *phy,
		   int *n_pos, double *pos, double *recrate,
		   double *window, double *work);

/**********************************************************************
 * estimate smoothed recombination rate 
 *
 * n_mar = number of markers
 * gen   = cM positions of markers (non-decreasing)
 * phy   = Mbp positions of markers (strictly increasing)
 * n_pos = number of positions at which to calculate the rec rate
 * pos   = Mbp positions at which to calculate rec rate (non-decreasing)
 * recrate = on exit, to contain the estimated rec rates (same length as pos)
 * window = length of sliding window (in Mbp)
 * work  = workspace of length n_mar-1
 **********************************************************************/
void est_recrate(int n_mar, double *gen, double *phy,
		 int n_pos, double *pos, double *recrate,
		 double window, double *work);

/* end of recrate.h */
