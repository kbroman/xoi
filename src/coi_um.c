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
 * centromeres = positions of centromeres (in microns, length n)
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
                double *centromeres, int *group, int n_group,
                double intwindow, double coiwindow,
                double *intloc, int n_intloc,
                double *coiloc, int n_coiloc,
                double **Intensity, double *coincidence)
{
    int i, total_xo;
    double **IntensityVals, *intensityvals;
    double **AdjustedXOPos, *adjustedxopos;

    /* count total number of crossovers */
    total_xo = 0;
    for(i=0; i<n; i++)
        total_xo += n_xo[i];
    /* space for intensity values; same structure as XOLoc */
    intensityvals = (double *)R_alloc(total_xo, sizeof(double));
    IntensityVals = (double **)R_alloc(n, sizeof(double *));
    IntensityVals[0] = intensityvals;
    for(i=1; i<n; i++)
        IntensityVals[i] = IntensityVals[i-1] + n_xo[i-1];
    /* space for adjusted XO positions; same structure as XOLoc */
    adjustedxopos = (double *)R_alloc(total_xo, sizeof(double));
    AdjustedXOPos = (double **)R_alloc(n, sizeof(double *));
    AdjustedXOPos[0] = adjustedxopos;
    for(i=1; i<n; i++)
        AdjustedXOPos[i] = AdjustedXOPos[i-1] + n_xo[i-1];

    /* get adjusted XO positions: p-arm in (0, 0.5) and q-arm in (0.5, 1) */
    calc_adjusted_xo_pos(n, XOLoc, n_xo, sclength, centromeres, AdjustedXOPos);

    /* estimate the intensity functions */
    for(i=0; i<n_group; i++)
        est_coi_um_intensity(n, AdjustedXOPos, n_xo, sclength, centromeres,
                             group, i+1, intwindow,
                             intloc, n_intloc, Intensity[i]);

    /* for each XO, find intensity at nearest calculated position */
    grab_intensities(n, AdjustedXOPos, n_xo, group, intloc, n_intloc,
                     Intensity, IntensityVals);

    /* estimate coincidence */
    est_coi_um_coincidence(n, XOLoc, IntensityVals, n_xo,
                           sclength, centromeres, intwindow, coiwindow,
                           coiloc, n_coiloc, coincidence);

}

/* to be called from R */
void R_est_coi_um(int *n, double *xoloc, int *n_xo, double *sclength,
                  double *centromeres, int *group, int *n_group,
                  double *intwindow, double *coiwindow,
                  double *intloc, int *n_intloc,
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

    est_coi_um(*n, XOLoc, n_xo, sclength, centromeres, group, *n_group,
               *intwindow, *coiwindow, intloc, *n_intloc,
               coiloc, *n_coiloc, Intensity, coincidence);
}

/* estimate intensity function for one group */
void est_coi_um_intensity(int n, double **AdjustedXOPos, int *n_xo,
                          double *sclength, double *centromeres,
                          int *group, int which_group,
                          double intwindow,
                          double *intloc, int n_intloc,
                          double *intensity)
{
    int i, j, k, count;

    /* this is definitely not the most efficient way to do this */
    /* sufficient for small data sets, but should be re-worked */
    for(i=0; i<n_intloc; i++) {
        intensity[i] = 0.0;

        count=0;
        for(j=0; j<n; j++) {
            if(group[j] == which_group) {
                for(k=0; k<n_xo[j]; k++) {
                    if(AdjustedXOPos[j][k] >= intloc[i]-intwindow/2.0 &&
                       AdjustedXOPos[j][k] <= intloc[i]+intwindow/2.0)
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

/* grab the intensities that correspond to each XOLoc position */
void grab_intensities(int n, double **XOLoc, int *n_xo,
                      int *group, double *intloc, int n_intloc,
                      double **Intensity, double **IntensityVal)
{
    int i, j, wh;

    for(i=0; i<n; i++) {
        for(j=0; j<n_xo[i]; j++) {
            wh = find_index_of_closest_value(XOLoc[i][j], n_intloc, intloc);
            IntensityVal[i][j] = Intensity[group[i]-1][wh];
        }
    }
}

/* find index of element in vec that is closest to x */
/* with ties, we just pick the first one */
int find_index_of_closest_value(double x, int n, double *vec)
{
    int i, index;
    double minimum, cur;

    index=0;
    minimum=fabs(vec[0]-x);

    for(i=1; i<n; i++) {
        cur = fabs(vec[i]-x);
        if(cur < minimum) {
            index=i;
            minimum = cur;
        }
    }
    return(index);
}

/* calculate the adjusted XO positions */
/* p-arm in (0,0.5); q-arm in (0.5, 1) */
void calc_adjusted_xo_pos(int n, double **XOLoc, int *n_xo,
                          double *sclength, double *centromeres,
                          double **AdjustedXOPos)
{
    int j, k;

    for(j=0; j<n; j++) {
        for(k=0; k<n_xo[j]; k++) {
            /* position -> (0,0.5) for p-arm and (0.5,1) for q-arm */
            if(XOLoc[j][k] <= centromeres[j])
                AdjustedXOPos[j][k] = XOLoc[j][k]/centromeres[j]/2.0;
            else
                AdjustedXOPos[j][k] = (XOLoc[j][k]-centromeres[j])/(sclength[j]-centromeres[j])/2.0 + 0.5;
        }
    }
}

/* estimate coincidence */
void est_coi_um_coincidence(int n, double **XOLoc, double **IntensityVals,
                            int *n_xo, double *sclength, double *centromeres,
                            double intwindow, double coiwindow, double *coiloc,
                            int n_coiloc, double *coincidence)
{
    int i, j1, j2, k;
    double *denom, d;
    double factor1, factor2;

    /* space for denominator; zero it out */
    denom = (double *)R_alloc(n_coiloc, sizeof(double));
    for(k=0; k<n_coiloc; k++)
        coincidence[k] = denom[k] = 0.0;

    for(i=0; i<n; i++) { /* loop over cells */

        for(k=0; k<n_coiloc; k++)
            denom[k] += (sclength[i] - coiloc[k]); /* totally not sure about this */

        /* loop over pairs of XO locations */
        for(j1=0; j1<n_xo[i]-1; j1++) {
            for(j2=j1+1; j2<n_xo[i]; j2++) {

                d = fabs(XOLoc[i][j1] - XOLoc[i][j2]); /* distance between XOs */
                for(k = 0; k<n_coiloc; k++) {
                    if(fabs(d - coiloc[k]) < coiwindow/2.0) {

                        /* scale intensities differently on the p-arm and the q-arm */
                        factor1 = XOLoc[i][j1] < centromeres[i] ? centromeres[i]*2.0 : (sclength[i]-centromeres[i])*2.0;
                        factor2 = XOLoc[i][j2] < centromeres[i] ? centromeres[i]*2.0 : (sclength[i]-centromeres[i])*2.0;

                        coincidence[k] += 1.0/(IntensityVals[i][j1]/factor1 *
                                               IntensityVals[i][j2]/factor2 *
                                               coiwindow);
                    }
                }

            }
        }
    }

    for(k=0; k<n_coiloc; k++) {
        coincidence[k] /= (denom[k]);
        if(coiloc[k] < coiwindow/2.0) /* adjust for edge effect */
            coincidence[k] *= coiwindow/(coiwindow/2.0 + coiloc[k]);
    }
}
