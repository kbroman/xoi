/**********************************************************************
 *
 * coincidence.h
 *
 * copyright (c) 2006-7, Karl W Broman
 *
 * last modified Apr, 2007
 * first written Dec, 2006
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * Contains: R_est_coi, est_coi, runningmean 
 *
 **********************************************************************/

/* R wrapper */
void R_est_coi(int *n_ind, int *n_mar, int *n_pair,
	       double *map, int *geno, double *d, 
	       double *coi1, double *coi2, 
	       int *n_keep, double *window);

/* estimate coincidence function */
void est_coi(int n_ind, int n_mar, int n_pair,
	     double *map, int **Geno, double *d,
	     double *coi1, double *coi2, 
	     int *n_keep, double window);

/**********************************************************************
 * runningmean
 *
 * Get running mean or sum within a specified bp-width window
 *
 * method = 1 -> sum
 *        = 2 -> mean
 *        = 3 -> median
 *
 **********************************************************************/
void runningmean(int n, double *pos, double *value, double *result, 
		 double window, int method, double *work);

/* end of coincidence.h */
