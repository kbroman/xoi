/**
 * coi_um.c
 *
 * Estimate the coincidence as a function of micron distance, with
 * data on XO locations in microns plus SC length in microns.
 **/

#include "coi_um.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h> /* for Rprintf() */

/*
 * n = number of cells
 * xoloc = XO locations in each cell (in microns)
 * n_xo = number of XO in each cell (length n)
 * sclength = SC lengths (in microns, length n)
 * group = vector of groups with common intensity function, {1, ..., n_group}
 * n_group = number of groups
 * intwindow = window for smoothing intensity function
 * coiwindow = window for smoothing coincidence
 * intloc = positions at which to calculate intensity
 * n_intloc = length of intloc
 * coiloc = values at which to calculate coincidence
 * n_coiloc = length of coiloc
 * intensity = vector of length n_intloc * n_group, to contain estimated intensities
 * coincidence = vector of length n_coiloc, to contain estimated coincidence

 */
void est_coi_um(int n, double **XOLoc, int *n_xo, double *sclength,
                int *group, int n_group, double intwindow,
                double coiwindow, double *intloc, int n_intloc,
                double *coiloc, int n_coiloc,
                double **Intensity, double *coincidence)
{
    int i;

    /* estimate the intensity functions */
    for(i=0; i<n_group; i++)
        est_coi_um_intensity(n, XOLoc, n_xo, sclength,
                             group, i+1, intwindow,
                             intloc, n_intloc, Intensity[i]);
}

/* to be called from R */
void R_est_coi_um(int *n, double *xoloc, int *n_xo, double *sclength,
                  int *group, int *n_group, double *intwindow,
                  double *coiwindow, double *intloc, int *n_intloc,
                  double *coiloc, int *n_coiloc,
                  double *intensity, double *coincidence)
{
    double **XOLoc, **Intensity;
    int i;

    /* set up ragged array for XO locations */
    XOLoc = (double **)R_alloc(*n, sizeof(double *));
    XOLoc[0] = xoloc;
    for(i=1; i<*n; i++)
        XOLoc[i] = XOLoc[i-1] + n_xo[i-1];
  
    /* set up matrix for intensity values */
    Intensity = (double **)R_alloc(*n_group, sizeof(double *));
    Intensity[0] = intensity;
    for(i=1; i<*n_group; i++)
        Intensity[i] = Intensity[i-1] + *n_intloc;
  
    est_coi_um(*n, XOLoc, n_xo, sclength, group, *n_group,
               *intwindow, *coiwindow, intloc, *n_intloc,
               coiloc, *n_coiloc, Intensity, coincidence);
}

/* estimate intensity function for one group */
void est_coi_um_intensity(int n, double **XOLoc, int *n_xo,
                          double *sclength, int *group,
                          int which_group, double intwindow,
                          double *intloc, int n_intloc,
                          double *intensity)
{
    int i, j, k, count;

    for(i=0; i<n_intloc; i++) {
        intensity[i] = 0.0;

        count=0;
        for(j=0; j<n; j++) {
            if(group[j] == which_group) {
                for(k=0; k<n_xo[j]; k++) {
                    if(XOLoc[j][k]/sclength[j] >= intloc[i]-intwindow/2.0 &&
                       XOLoc[j][k]/sclength[j] <= intloc[i]+intwindow/2.0)
                        intensity[i] += 1.0;
                }
                count++;
            }
        }
        
        intensity[i] /= (double)count;
        if(intloc[i] < intwindow/2.0)
            intensity[i] /= (intloc[i] + intwindow/2.0);
        else if(intloc[i] > 1.0-intwindow/2.0)
            intensity[i] /= (1.0-intloc[i] + intwindow/2.0);
        else
            intensity[i] /= intwindow;
    }

}
