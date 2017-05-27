/**********************************************************************
 *
 * coincidence.c
 *
 * copyright (c) 2006-7, Karl W Broman
 *
 * last modified Apr, 2007
 * first written Dec, 2006
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
 *     at https://www.r-project.org/Licenses/GPL-3
 *
 * Contains: R_est_coi, est_coi, runningmean
 *
 **********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "coincidence.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>

/* R wrapper */
void R_est_coi(int *n_ind, int *n_mar, int *n_pair,
               double *map, int *geno, double *d,
               double *coi1, double *coi2,
               int *n_keep, double *window)
{
    int **Geno, j;

    /* make genotype data doubly-indexed */
    Geno = (int **)R_alloc(*n_mar, sizeof(int *));
    Geno[0] = geno;
    for(j=1; j < *n_mar; j++)
        Geno[j] = Geno[j-1] + *n_ind;

    est_coi(*n_ind, *n_mar, *n_pair, map, Geno,
            d, coi1, coi2, n_keep, *window);
}

/* estimate coincidence function */
void est_coi(int n_ind, int n_mar, int n_pair,
             double *map, int **Geno, double *d,
             double *coi1, double *coi2,
             int *n_keep, double window)
{
    double *rf, *mapd, *top, *bottom, *temp, *work;
    int i, j, k, s, nmarm1, *temp_idx;

    nmarm1 = n_mar - 1;

    /* allocate space */
    rf = (double *)R_alloc(nmarm1, sizeof(double));
    mapd = (double *)R_alloc(nmarm1, sizeof(double));
    top = (double *)R_alloc(n_pair, sizeof(double));
    bottom = (double *)R_alloc(n_pair, sizeof(double));
    temp = (double *)R_alloc(n_pair, sizeof(double));
    temp_idx = (int *)R_alloc(n_pair, sizeof(int));
    work = (double *)R_alloc(n_pair, sizeof(double));

    R_CheckUserInterrupt(); /* check for ^C */

    /* midpoints of intervals */
    for(i=0; i<nmarm1; i++)
        mapd[i] = (map[i]+map[i+1])/2.0;

    R_CheckUserInterrupt(); /* check for ^C */

    /* inter-interval distances */
    for(j=0, s=0; j<nmarm1-1; j++)
        for(k=(j+1); k<nmarm1; k++, s++)
            d[s] = mapd[k] - mapd[j];

    R_CheckUserInterrupt(); /* check for ^C */

    /* recombination fractions */
    for(j=0; j<nmarm1; j++) {
        rf[j] = 0.0;
        for(i=0; i<n_ind; i++) {
            if(Geno[j][i] != Geno[j+1][i])
                rf[j] += 1.0;
        }
        rf[j] /= (double)n_ind;

        R_CheckUserInterrupt(); /* check for ^C */
    }

    /* top and bottom of the coincidence function */
    for(j=0, s=0; j<nmarm1-1; j++) {
        for(k=(j+1); k<nmarm1; k++, s++) {

            top[s] = 0.0;
            bottom[s] = rf[j]*rf[k];

            for(i=0; i<n_ind; i++) {
                if(Geno[j][i] != Geno[j+1][i] &&
                   Geno[k][i] != Geno[k+1][i])
                    top[s] += 1.0;
            }
            top[s] /= (double)n_ind;

            R_CheckUserInterrupt(); /* check for ^C */
        }
    }

    /* ratio, then smooth */
    for(i=0; i<n_pair; i++) {
        if(fabs(bottom[i]) < 1e-12) coi2[i] = NA_REAL; /* to be ignored */
        else coi2[i] = top[i]/bottom[i];
    }

    R_CheckUserInterrupt(); /* check for ^C */

    /* sort d, and also top and bottom to match */
    /* first, create an index */
    for(i=0; i<n_pair; i++)
        temp_idx[i] = i;
    rsort_with_index(d, temp_idx, n_pair);

    R_CheckUserInterrupt(); /* check for ^C */

    /* sort then running means on coi2 */
    for(i=0; i<n_pair; i++)
        temp[i] = coi2[temp_idx[i]];
    runningmean(n_pair, d, temp, coi2, window, 2, work);

    R_CheckUserInterrupt(); /* check for ^C */

    /* sort top and then do running mean */
    for(i=0; i<n_pair; i++)
        temp[i] = top[temp_idx[i]];
    runningmean(n_pair, d, temp, top, window, 2, work);

    R_CheckUserInterrupt(); /* check for ^C */

    /* sort bottom and then do running mean */
    for(i=0; i<n_pair; i++)
        temp[i] = bottom[temp_idx[i]];
    runningmean(n_pair, d, temp, bottom, window, 2, work);

    R_CheckUserInterrupt(); /* check for ^C */

    for(i=0; i<n_pair; i++)
        coi1[i] = top[i]/bottom[i];

    R_CheckUserInterrupt(); /* check for ^C */

    /* now just save the unique values */
    for(j=0, i=1, *n_keep=1; i<n_pair; i++) {
        if(d[i] > d[j]) {
            coi1[*n_keep] = coi1[i];
            coi2[*n_keep] = coi2[i];
            d[*n_keep] = d[i];
            (*n_keep)++;
            j = i;
        }
    }
}


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
                 double window, int method, double *work)
{
    int lo, ns;
    int i, j;

    window /= 2.0;

    lo=0;
    for(i=0; i<n; i++) {
        result[i] = 0.0; ns=0;
        for(j=lo; j<n; j++) {
            if(pos[j] < pos[i]-window) lo = j+1;
            else if(pos[j] > pos[i]+window) break;
            else {

                if(!ISNAN(value[j])) {
                    if(method==1 || method==2)
                        result[i] += value[j];
                    else
                        work[ns] = value[j];

                    ns++;
                }
            }
        }
        if(method==2) result[i] /= (double)ns;
        if(method==3) {
            R_rsort(work, ns);
            if(ns % 2) /* odd */
                result[i] = work[(ns-1)/2];
            else /* even */
                result[i] = (work[ns/2-1]+work[ns/2])/2.0;
        }
    }

}



/* end of coincidence.c */
