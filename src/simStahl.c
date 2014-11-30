/**********************************************************************
 *
 * simStahl.c
 *
 * copyright (c) 2006-2014, Karl W Broman
 *
 * last modified Sep, 2014
 * first written Nov, 2006
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
 * Contains: simGamma
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GammaS.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "simStahl.h"

void simStahl(int *n_sim, double *nu, double *p, double *L,
              int *nxo, double *loc, int *max_nxo,
              int *n_bins4start)
{
    double **Loc, scale;
    double curloc=0.0, u;
    double *startprob, step;
    int i, j, n_nixo;

    /* re-organize loc as a doubly index array */
    Loc = (double **)R_alloc(*n_sim, sizeof(double *));
    Loc[0] = loc;
    for(i=1; i < *n_sim; i++)
        Loc[i] = Loc[i-1] + *max_nxo;

    GetRNGstate();

    if(fabs(*nu - 1.0) < 1e-8) { /* looks like a Poisson model */
        for(i=0; i< *n_sim; i++) {
            R_CheckUserInterrupt(); /* check for ^C */

            nxo[i] = rpois(*L);
            if(nxo[i] > *max_nxo)
                error("Exceeded maximum number of crossovers.");

            for(j=0; j < nxo[i]; j++)
                Loc[i][j] = runif(0.0, *L);
        }
    }
    else {
        scale = 1.0 / (2.0 * *nu * (1.0 - *p));

        /* set up starting distribution */
        startprob = (double *)R_alloc(*n_bins4start, sizeof(double));
        step = *L/(double)*n_bins4start;

        startprob[0] = 2.0*(1.0 - *p)*pgamma(0.5*step, *nu, scale, 0, 0)*step;
        for(i=1; i< *n_bins4start; i++) {
            R_CheckUserInterrupt(); /* check for ^C */

            startprob[i] = startprob[i-1] +
                2.0*(1.0 - *p)*pgamma(((double)i+0.5)*step, *nu, scale, 0, 0)*step;
        }

        for(i=0; i< *n_sim; i++) {
            R_CheckUserInterrupt(); /* check for ^C */

            nxo[i] = 0;

            /* locations of chiasmata from the gamma model */
            /* shape = nu, rate = 2*nu*(1-p) [scale = 1/{2*nu*(1-p)}] */

            u = unif_rand();
            if( u > startprob[*n_bins4start-1] )
                curloc = *L+1;
            else {
                for(j=0; j< *n_bins4start; j++) {
                    if(u <= startprob[j]) {
                        curloc = ((double)j+0.5)*step;
                        if(unif_rand() < 0.5) {
                            nxo[i] = 1;
                            Loc[i][0] = curloc;
                        }
                        break;
                    }
                }
            }

            if(curloc < *L) {
                while(curloc < *L) {
                    curloc += rgamma(*nu, scale);
                    if(curloc < *L && unif_rand() < 0.5) {
                        if(nxo[i] > *max_nxo)
                            error("Exceeded maximum number of crossovers.");

                        Loc[i][nxo[i]] = curloc;
                        (nxo[i])++;
                    }
                }
            }

            /* locations of crossovers from the no interference mechanism */
            if(*p > 0) {
                n_nixo = rpois(*L * *p);
                if(n_nixo > *max_nxo)
                    error("Exceeded maximum number of crossovers.");

                for(j=0; j < n_nixo; j++)
                    Loc[i][nxo[i]+j] = runif(0.0, *L);
                nxo[i] += n_nixo;
            }
        }
    }

    /* sort the results */
    for(i=0; i< *n_sim; i++)
        R_rsort(Loc[i], nxo[i]);

    PutRNGstate();
}


/* version when nu = m+1 is an integer
 *
 * m = interference parameter (m=0 gives no interference)
 * p = proportion of chiasmata from no interference process
 * L = length of chromosome (in cM)
 * Lstar = revised length for simulating numbers of chiasmata, for case of obligate chiasma
 *         on same scale as L
 * nxo = on output, the number of crossovers
 * Loc = on output, the locations of the crossovers
 * max_nxo = maximum no. crossovers allowed (length of loc)
 * obligate_chiasma = 1 if require at least one chiasma (0 otherwise)
 *
 */
void simStahl_int(int n_sim, int m, double p, double L,
                  double Lstar, int *nxo, double **Loc,
                  int max_nxo, int obligate_chiasma)
{
    int i, j, k, n_nichi, n_pts, n_ichi, first, max_pts;
    double *ptloc;
    double lambda1, lambda2;

    /* space for locations of chiasmata and intermediate pts */
    max_pts = 2*max_nxo*(m+1);
    ptloc = (double *)R_alloc(max_pts, sizeof(double));

    GetRNGstate();

    if(m==0) { /* looks like a Poisson model */
        for(i=0; i< n_sim; i++) {
            R_CheckUserInterrupt(); /* check for ^C */

            if(obligate_chiasma) {
                /* no. chiasmata, required >= 1 */
                while((n_ichi = rpois(Lstar/50.0)) == 0);
                /* no crossovers by thinning 1/2 */
                nxo[i] = rbinom((double)n_ichi, 0.5);
            }
            else
                nxo[i] = rpois(Lstar/100.0);

            if(nxo[i] > max_nxo)
                error("Exceeded maximum number of crossovers.");

            for(j=0; j < nxo[i]; j++)
                Loc[i][j] = runif(0.0, L);
        }
    }
    else {
        lambda1 = Lstar/50.0 * (m+1) * (1.0 - p);
        lambda2 = Lstar/50.0 * p;

        for(i=0; i< n_sim; i++) {
            while(1) {
                R_CheckUserInterrupt(); /* check for ^C */

                /* simulate no. chiasmata + intermediate pts from interference process */
                n_pts = rpois(lambda1);

                /* simulate location of the first */
                first = random_int(0, m);

                if(first > n_pts) n_ichi = 0;
                else n_ichi = n_pts/(m+1) + (int)(first < (n_pts % (m+1)));

                /* simulate no. chiamata from the no-interference model */
                n_nichi = rpois(lambda2);

                if(!obligate_chiasma || n_ichi + n_nichi > 0) break;
            }

            /* simulate no. chiasmta + intermediate points */
            /* first check if we have space */
            if(n_pts > max_pts) {
                ptloc = (double *)S_realloc((char *)ptloc, n_pts*2, max_pts, sizeof(double));
                max_pts = n_pts*2;
            }

            for(j=0; j<n_pts; j++)
                ptloc[j] = runif(0.0, L);

            /* sort them */
            R_rsort(ptloc, n_pts);

            /* take every (m+1)st */
            for(j=first, k=0; j<n_pts; j += (m+1), k++)
                ptloc[k] = ptloc[j];
            n_ichi = k;

            /* simulate chiasmata from no-interference model */
            for(j=0; j<n_nichi; j++)
                ptloc[n_ichi + j] = runif(0.0, L);

            /* sort the combined ones */
            R_rsort(ptloc, n_ichi + n_nichi);

            /* thin by 1/2 */
            nxo[i] = 0;
            for(j=0; j<n_ichi + n_nichi; j++) {
                if(unif_rand() < 0.5) {
                    Loc[i][nxo[i]] = ptloc[j];
                    (nxo[i])++;
                }
            }

        } /* loop over no. simulations */
    } /* m > 0 */


    PutRNGstate();
}

/* simStahl_int: called from R */
void R_simStahl_int(int *n_sim, int *m, double *p, double *L,
                    double *Lstar, int *nxo, double *loc,
                    int *max_nxo, int *obligate_chiasma)
{
    double **Loc;
    int i;

    /* re-organize loc as a doubly index array */
    Loc = (double **)R_alloc(*n_sim, sizeof(double *));
    Loc[0] = loc;
    for(i=1; i < *n_sim; i++)
        Loc[i] = Loc[i-1] + *max_nxo;

    simStahl_int(*n_sim, *m, *p, *L, *Lstar, nxo, Loc, *max_nxo, *obligate_chiasma);
}


/* random integer from {low, low+1, ..., high} */
int random_int(int low, int high)
{
    if(high < low)
        error("Must have high >= low");

    if(high == low) return(low); /* just one possible value */

    return((int)(unif_rand()*(double)(high-low+1)) + low);
}

/* end of simStahl.c */
