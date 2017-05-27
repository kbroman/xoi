/**********************************************************************
 *
 * chiasma.c
 *
 * copyright (c) 1998-2006, Karl W Broman
 *
 * last modified Dec, 2006
 * first written ~Oct, 1998
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
 * Contains: chiasma, chiasma_estep, chiasma_mstep
 *
 *     This program is to be called from R.  It takes data for
 * the number of crossovers (recombination events) observed per
 * chromosome for a set of meioses, and estimates the underlying
 * chiasma distribution and the corresponding crossover distribution
 * under several models, all assuming no chromatid interference:
 *
 *     (1)   Poisson dist'n with Pr(0 chiasmata) = 0
 *     (2)   Poisson dist'n
 *     (3)   Chiasma dist'n to vary freely, but with Pr(0 chiasmata) = 0
 *     (4)   Chiasma dist'n to vary freely
 *
 * Input:  [all pointers]
 *     xo             = vector containing the numbers of crossovers
 *     n_xo           = length of xo
 *     max_ch         = max number of chiasmata allowed
 *     p_ch           = dist'n of # of chiasmata (max_ch+1)*5
 *     p_xo           = dist'n of # of crossovers (max_ch+1)*5
 *     lambda         = mean # XO's for models 1, 2
 *     work           = workspace of length (max_ch+1)*n_xo+2+2*(max_ch+1)
 *     n_iter         = number of EM iterations (length 5: max + 5 x # done)
 *     tol            = tolerance value
 *
 * p_ch and lambda should be given starting values to begin
 * if they're zero, the program chooses the starting values
 *
 * in the program, lambda = mean # of chiasmata
 * in the input/output, lambda = mean # of crossovers = mean chiasmata/2
 *
 **********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "chiasma.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

/* main chiasma function */
void chiasma(int *xo, int *n_xo, int *max_ch,
             double *p_ch, double *p_xo, double *lambda,
             double *work, int *n_iter, double *tol)
{
    int i, j, k;
    int maxchp1, flag;
    double *w, *oldlambda, *oldp, temp=0.0;

    maxchp1 = *max_ch + 1;
    w = work;
    oldlambda = w + maxchp1 * *n_xo;
    oldp = oldlambda + 2;

    for(i=0; i < 4; i++) {  /* loop over the four models */
        R_CheckUserInterrupt(); /* check for ^C */

        if(i==1) { /* model 2 requires no EM */
            lambda[i] = 0.0;
            for(j=0; j < *n_xo; j++) lambda[i] += (double)xo[j];
            lambda[i] /= (double)(*n_xo);
            lambda[i] *= 2.0;
            n_iter[i+1] = 0;
        }
        else {
            /* starting values */
            if(i==0 && lambda[i] < 1e-14) {
                lambda[i] = 0.0;
                for(j=0; j < *n_xo; j++) {
                    lambda[i] += (double)xo[j];
                }
                lambda[i] /= (double)(*n_xo);
            }
            else {
                for(j=0; j < maxchp1; j++) {
                    oldp[j+(i-2)*maxchp1] = p_ch[j+i*maxchp1] = p_ch[j+(i-2)*maxchp1];
                }
                if(i==2) { oldp[(i-2)*maxchp1] = p_ch[i*maxchp1] = 0.0; }
            }
            /* lambda from crossover mean to chiasma mean */
            lambda[i] *= 2.0;
            oldlambda[i] = lambda[i];

            /* begin EM algorithm */
            for(j=0; j < *n_iter; j++) { /* loop over EM steps */
                R_CheckUserInterrupt(); /* check for ^C */

                /* E step */
                chiasma_estep(xo, *n_xo, w, maxchp1, p_ch+i*maxchp1,
                              lambda, i);

                /* M step */
                chiasma_mstep(xo, *n_xo, w, maxchp1, p_ch+i*maxchp1,
                              lambda, i, *n_iter, *tol);

                /* converged? */
                flag = 0;
                if(i==0) {
                    if(fabs(lambda[i] - oldlambda[i]) > *tol/100.0) {
                        flag = 1;
                    }
                    oldlambda[i] = lambda[i];
                }
                else {
                    for(k=0; k < maxchp1; k++) {
                        if(fabs(p_ch[k+i*maxchp1] - oldp[k+(i-2)*maxchp1]) > *tol) {
                            flag = 1;
                        }
                        oldp[k+(i-2)*maxchp1] = p_ch[k+i*maxchp1];
                    }
                }
                if(flag==0) { break; }
            }
            n_iter[i+1] = j+1;
        }

        /* determine p_ch for i=0,1 */
        if(i==0 || i==1) {  /* poisson dist'n */
            for(j=0; j<maxchp1; j++)
                p_ch[j+i*maxchp1] = dpois((double)j, lambda[i], 0);
        }
        if(i==0) {  /* truncated poisson dist'n */
            p_ch[i*maxchp1] = 0.0;
            for(j=1; j<maxchp1; j++)
                p_ch[j+i*maxchp1] /= (1.0-exp(-lambda[i]));
        }

        /* determin p_xo */
        for(j=0; j<maxchp1; j++) { /* mixture dist'n */
            p_xo[j+i*maxchp1] = 0.0;
            for(k=j; k<maxchp1; k++) {
                if(i==0 || i==1) {
                    temp = dpois((double)k, lambda[i], 0);
                }
                if(i==0) {
                    if(k==0) temp = 0.0;
                    else temp /= (1.0-exp(-lambda[i]));
                }
                if(i==2 || i==3) {
                    if(k==0 && i==2) temp = 0.0;
                    else temp = p_ch[k+i*maxchp1];
                }
                p_xo[j+i*maxchp1] += dbinom((double)j,(double)k,0.5,0)*temp;
            }
        }

        /* lambda from chiasma mean to crossover mean */
        lambda[i] /= 2.0;

    } /* loop over models */

} /* end of chiasma() */





/* E step */
void chiasma_estep(int *xo, int n_xo, double *w, int maxchp1,
                   double *p_ch, double *lambda, int model)
{
    int i, j;
    double temp;

    for(j=0; j<maxchp1; j++) {

        if(model== 0 || model==1) {
            temp = dpois((double)j, lambda[model], 0);
            if(model==0) {
                if(j==0) temp = 0.0;
                else temp /= (1.0-exp(-lambda[model]));
            }
        }
        else {
            if(model==2 && j==0) { temp = 0.0; }
            else temp = p_ch[j];
        }

        for(i=0; i<n_xo; i++) {
            if(j<xo[i]) w[i+j*n_xo] = 0.0;
            else w[i+j*n_xo] += temp*dbinom((double)xo[i],(double)j,0.5,0);
        }
    }
    for(i=0; i<n_xo; i++) {
        temp = 0.0;
        for(j=0; j<maxchp1; j++) temp += w[i+j*n_xo];
        for(j=0; j<maxchp1; j++) w[i+j*n_xo] /= temp;
    }
}



/* M step */
/* n_iter and tol are provided for the trucated poisson, where  *
 * I need to use Newton's method to do the maximization         */
void chiasma_mstep(int *xo, int n_xo, double *w, int maxchp1,
                   double *p_ch, double *lambda, int model,
                   int n_iter, double tol)
{
    int i, j;
    double sum, oldlambda;
    double f, fp;

    if(model==0 || model==1) {

        /* expected sum */
        sum = 0.0;
        for(i=0; i<n_xo; i++) {
            for(j=1; j<maxchp1; j++) {
                sum += (double)j * w[i+j*n_xo];
            }
        }

        if(model==1) { /* lambda = ave */
            lambda[model] = sum / (double)n_xo;
        }
        else { /* have to do newton's method */

            oldlambda = lambda[model];

            for(j=0; j<n_iter; j++) {

                f = -(double)n_xo/(1.0-exp(-lambda[model])) + sum/lambda[model];
                fp = -sum/lambda[model]/lambda[model] +
                    (double)n_xo*exp(-lambda[model])/
                    (1.0-exp(-lambda[model]))/(1.0-exp(-lambda[model]));
                lambda[model] -= f/fp;

                if(fabs(oldlambda - lambda[model]) < tol/100.0) { break; }
                oldlambda = lambda[model];
            }

        }

    }

    else {
        for(j=0; j<maxchp1; j++) {
            p_ch[j] = 0.0;
            for(i=0; i<n_xo; i++) {
                p_ch[j] += w[i+j*n_xo];
            }
            p_ch[j] /= (double)n_xo;
        }
        if(model==2) { p_ch[0] = 0.0; }

    }

}


/* end of chiasma.c */
