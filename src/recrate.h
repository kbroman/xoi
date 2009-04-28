/**********************************************************************
 *
 * recrate.h
 *
 * copyright (c) 2008, Karl W Broman
 *
 * last modified Aug, 2008
 * first written Aug, 2008
 * Licensed under the GNU General Public License version 2 (June, 1991)
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
