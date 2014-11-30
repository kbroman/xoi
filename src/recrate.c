/**********************************************************************
 *
 * recrate.c
 *
 * copyright (c) 2008-2010, Karl W Broman
 *
 * last modified May, 2010
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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "recrate.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>

/* R wrapper */
void R_est_recrate(int *n_mar, double *gen, double *phy,
                   int *n_pos, double *pos, double *recrate,
                   double *window, double *work)
{
    est_recrate(*n_mar, gen, phy, *n_pos, pos, recrate, *window, work);
}

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
                 double window, double *work)
{
    int i, j;
    double denom=0.0, half=window/2.0, start, end;

    /* get estimated rate in each interval */
    for(i=0; i<n_mar-1; i++)
        work[i] = (gen[i+1]-gen[i])/(phy[i+1]-phy[i]);

    for(i=0; i<n_pos; i++) {
        start = pos[i]-half;
        end = pos[i]+half;

        if(start < phy[0]) {
            if(end < phy[1]) {
                recrate[i] = work[0];
            }
            else {
                recrate[i] += work[0]*(phy[1]-phy[0])/(end-phy[0]);
                for(j=1; j<n_mar-1; j++) {
                    if(phy[j+1] > end) {
                        recrate[i] += work[j]*(end-phy[j])/(end-phy[0]);
                        break;
                    }
                    else {
                        recrate[i] += work[j]*(phy[j+1]-phy[j])/(end-phy[0]);
                    }
                }
            }
        } /* end overlap p-ter */

        else if(end > phy[n_mar-1]) {
            denom = phy[n_mar-1]-start;
            if(start > phy[n_mar-2])
                recrate[i] = work[n_mar-2];
            else {
                for(j=0; j<n_mar-1; j++) {
                    if(phy[j+1] > start) {
                        if(phy[j] < start)
                            recrate[i] += work[j]*(phy[j+1]-start)/denom;
                        else
                            recrate[i] += work[j]*(phy[j+1]-phy[j])/denom;
                    }
                }
            }
        } /* end overlap q-ter */

        else {
            for(j=0; j<n_mar-1; j++) {
                if(phy[j] > end) break;
                if(phy[j+1] > start) {
                    if(phy[j+1] > end) {
                        if(phy[j] < start) {
                            recrate[i] = work[j];
                        }
                        else {
                            recrate[i] += work[j]*(end-phy[j])/window;
                        }
                    }
                    else {
                        if(phy[j] < start)
                            recrate[i] += work[j]*(phy[j+1]-start)/window;
                        else
                            recrate[i] += work[j]*(phy[j+1]-phy[j])/window;
                    }
                }
            }
        } /* end internal case */



    } /* end loop over positions */
}

/* end of recrate.c */
