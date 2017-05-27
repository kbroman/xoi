/**********************************************************************
 *
 * fitStahl.c
 *
 * copyright (c) 2009-2012, Karl W Broman
 *
 * last modified Nov, 2012
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
 *     at https://www.r-project.org/Licenses/GPL-3
 *
 * Part of the R/xoi package
 * Contains: R_stahl_loglik, stahl_loglik, stahl_loglik_byind
 *           R_stahl_loglik_F2, stahl_loglik_F2
 *           addlog
 *
 * Functions to calculate and maximize the likelihood for the Stahl
 * model, given data on the crossover locations
 *
 **********************************************************************/

#include <math.h>
#include "GammaS.h"
#include "fitStahl.h"
#include "stahl_util.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>

void R_stahl_loglik(int *n_ind, int *n_xo, double *xoloc, double *chrlen,
                    int *n_nu, double *nu, double *p, double *loglik,
                    int *max_conv, double *intgr_tol, int *max_subd,
                    int *min_subd, int *constant_chrlen)
{
    int i;
    double **XOloc;

    /* reorganize XOloc as 2-d array */
    XOloc = (double **)R_alloc(*n_ind, sizeof(double *));
    XOloc[0] = xoloc;
    for(i=1; i<*n_ind; i++) {
        XOloc[i] = XOloc[i-1] + n_xo[i-1];
    }

    stahl_loglik(*n_ind, n_xo, XOloc, chrlen, *n_nu, nu, p,
                 loglik, *max_conv, *intgr_tol, *max_subd, *min_subd,
                 *constant_chrlen);
}

void R_stahl_loglik_F2(int *n_ind, int *n_alternatives, int *n_products,
                       int *n_xo_per,
                       double *xoloc, double *chrlen,
                       int *n_nu, double *nu, double *p, double *loglik,
                       int *max_conv, double *intgr_tol, int *max_subd,
                       int *min_subd, int *constant_chrlen)
{
    int i;
    double **XOloc;

    /* reorganize XOloc as 2-d array */
    XOloc = (double **)R_alloc(*n_products, sizeof(double *));
    XOloc[0] = xoloc;
    for(i=1; i<*n_products; i++) {
        XOloc[i] = XOloc[i-1] + n_xo_per[i-1];
    }

    stahl_loglik_F2(*n_ind, n_alternatives, *n_products, n_xo_per,
                    XOloc, chrlen, *n_nu, nu, p,
                    loglik, *max_conv, *intgr_tol, *max_subd, *min_subd,
                    *constant_chrlen);
}



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
 * constant_chrlen = 1 if all chrlen's are the same; otherwise 0
 **********************************************************************/

void stahl_loglik(int n_ind, int *n_xo, double **XOloc, double *chrlen,
                  int n_nu, double *nu, double *p, double *loglik,
                  int max_conv, double intgr_tol, int max_subd, int min_subd,
                  int constant_chrlen)
{
    int i, j;
    double *indll;

    /* allocate space for individual log likelihoods */
    indll = (double *)R_alloc(n_ind, sizeof(double));

    for(i=0; i<n_nu; i++) {
        stahl_loglik_byind(n_ind, n_xo, XOloc, chrlen, nu[i], p[i], indll,
                           max_conv, intgr_tol, max_subd, min_subd, constant_chrlen);

        loglik[i] = 0.0;
        for(j=0; j<n_ind; j++)
            loglik[i] += indll[j];
    }

}

/* stahl_loglik_F2: for intercross */
void stahl_loglik_F2(int n_ind, int *n_alternatives, int n_products, int *n_xo_per,
                     double **XOloc, double *chrlen, int n_nu, double *nu, double *p,
                     double *loglik, int max_conv, double intgr_tol, int max_subd, int min_subd,
                     int constant_chrlen)
{
    int i, j, k, curspot;
    double *indll, ll=0.0;

    /* allocate space for individual log likelihoods */
    indll = (double *)R_alloc(n_products, sizeof(double));

    for(i=0; i<n_nu; i++) {
        /* calculate log likelihoods for each individual product */
        stahl_loglik_byind(n_products, n_xo_per, XOloc, chrlen, nu[i], p[i], indll,
                           max_conv, intgr_tol, max_subd, min_subd, constant_chrlen);

        loglik[i] = 0.0;
        for(j=0, curspot=0; j<n_ind; j++) { /* loop over individuals */
            for(k=0; k<n_alternatives[j]; k++, curspot += 2) { /* loop over alternatives for each ind'l */

                if(k==0) ll=indll[curspot] + indll[curspot+1];
                else ll = addlog(ll, indll[curspot] + indll[curspot+1]);

            } /* end loop over alternativs within ind'l */
            loglik[i] += ll;
        } /* end loop over ind'l */

    } /* end loop over parameter values */
}


/**********************************************************************
 * stahl_loglik_byind: calculate log likelihood for the Stahl model
 *       This calculates the log likelihood for each of a vector of individuals
 *
 * n_ind  = no. individuals
 * n_xo   = vector with no. crossovers for each individual
 * xoloc  = crossover locations xoloc[i][j] = XO j in individual i
 * chrlen = vector of chromosome lengths
 * nu     = interference parameter at which loglik is to be calculated
 * p      = proportion of crossovers from the no-interference pathway
 * loglik = (on output) the log likelihood for each individual (length n_ind)
 * max_conv  = maximum number of convolutions
 * intgr_tol = tolerance for integration
 * max_subd  = maximum number of subdivisions in integration
 * min_subd  = minimum number of subdivisions in integration
 * constant_chrlen = 1 if all chrlen's are the same; otherwise 0
 **********************************************************************/

void stahl_loglik_byind(int n_ind, int *n_xo, double **xoloc, double *chrlen,
                        double nu, double p, double *loglik,
                        int max_conv, double intgr_tol, int max_subd, int min_subd,
                        int constant_chrlen)
{
    int i, j, k, s, n_pat, max_n_xo;
    double *y, *z, cury, curz, patll;
    double noXOll=0.0, noXOll_ni=0.0, noXOll_i=0.0;
    int n_y, n_z;
    struct gamma_data info;

    max_n_xo=1;
    for(i=0; i<n_ind; i++)
        if(max_n_xo < n_xo[i]) max_n_xo = n_xo[i];
    y = (double *)R_alloc(max_n_xo, sizeof(double));
    z = (double *)R_alloc(max_n_xo, sizeof(double));

    info.max_conv = max_conv;

    setup_integr_par(intgr_tol, max_subd, min_subd, &(info.integr_info));

    if(constant_chrlen) { /* calculate the no XO logliklihoods just once */
        noXOll_ni = log(oneminus_Gstar_stahl(chrlen[0], 1.0, p, max_conv, info.integr_info));
        noXOll_i =  log(oneminus_Gstar_stahl(chrlen[0], nu, 1.0-p, max_conv, info.integr_info));
        noXOll = noXOll_ni + noXOll_i;
    }

    for(j=0; j<n_ind; j++) {
        R_CheckUserInterrupt(); /* check for ^C */

        if(n_xo[j] == 0) { /* no crossovers */
            if(constant_chrlen)
                loglik[j] = noXOll;
            else
                loglik[j] = log(oneminus_Gstar_stahl(chrlen[j], 1.0, p, max_conv, info.integr_info)) +
                    log(oneminus_Gstar_stahl(chrlen[j], nu, 1.0-p, max_conv, info.integr_info));
        }
        else {
            n_pat = 1 << n_xo[j]; /* no. patterns */
            for(k=0; k<n_pat; k++) { /* loop over patterns */

                cury = curz = 0;
                n_y = n_z = 0;
                for(s=0; s<n_xo[j]; s++) {
                    if(k & (1 << s)) { /* is XO from NI pathway? */
                        y[n_y] = xoloc[j][s] - cury;
                        cury = xoloc[j][s];
                        n_y++;
                    }
                    else {
                        z[n_z] = xoloc[j][s] - curz;
                        curz = xoloc[j][s];
                        n_z++;
                    }
                }

                /* no interference pathway */
                if(n_y==0) { /* no crossovers */
                    if(constant_chrlen)
                        patll = noXOll_ni;
                    else
                        patll = log(oneminus_Gstar_stahl(chrlen[j], 1.0, p, max_conv, info.integr_info));
                }
                else if(n_y==1) { /* exactly one crossover */
                    patll = log(gstar_stahl(y[0], 1.0, p, max_conv)) +
                        log(oneminus_Fstar_stahl(chrlen[j]-y[0], 1.0, p, max_conv));
                }
                else { /* more than one crossover */
                    patll = log(gstar_stahl(y[0], 1.0, p, max_conv));
                    for(s=1; s<n_y; s++)
                        patll += log(fstar_stahl(y[s], 1.0, p, max_conv));
                    patll += log(oneminus_Fstar_stahl(chrlen[j]-cury, 1.0, p, max_conv));
                }

                /* interference pathway */
                if(n_z==0) { /* no crossovers */
                    if(constant_chrlen)
                        patll += noXOll_i;
                    else
                        patll += log(oneminus_Gstar_stahl(chrlen[j], nu, 1.0-p, max_conv, info.integr_info));
                }
                else if(n_z==1) { /* exactly one crossover */
                    patll += log(gstar_stahl(z[0], nu, 1.0 - p, max_conv)) +
                        log(oneminus_Fstar_stahl(chrlen[j]-z[0], nu, 1.0 - p, max_conv));
                }
                else { /* more than one crossover */
                    patll += log(gstar_stahl(z[0], nu, 1.0 - p, max_conv));
                    for(s=1; s<n_z; s++)
                        patll += log(fstar_stahl(z[s], nu, 1.0 - p, max_conv));
                    patll += log(oneminus_Fstar_stahl(chrlen[j]-curz, nu, 1.0 - p, max_conv));
                }

                if(k==0)
                    loglik[j] = patll;
                else
                    loglik[j] = addlog_stahl(loglik[j], patll);

            } /* loop over patterns */

        } /* there are crossovers */
    } /* loop over individuals */

}

/**********************************************************************
 *
 * addlog
 *
 * Calculate addlog(a,b) = log[exp(a) + exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.
 *
 **********************************************************************/
#define THRESH 200.0

double addlog(double a, double b)
{
    if(b > a + THRESH) return(b);
    else if(a > b + THRESH) return(a);
    else return(a + log1p(exp(b-a)));
}

/* end of fitStahl.c */
