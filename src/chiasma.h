/**********************************************************************
 *
 * chiasma.h
 *
 * copyright (c) 1998-2006, Karl W Broman
 *
 * last modified Nov, 2006
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
 *     at http://www.r-project.org/Licenses/GPL-3
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

/* main chiasma function */
void chiasma(int *xo, int *n_xo, int *max_ch,
             double *p_ch, double *p_xo, double *lambda,
             double *work, int *n_iter, double *tol);



/* E step */
void chiasma_estep(int *xo, int n_xo, double *w, int maxchp1,
                   double *p, double *lambda, int model);



/* M step */
void chiasma_mstep(int *xo, int n_xo, double *w, int maxchp1,
                   double *p, double *lambda, int model,
                   int n_iter, double tol);


/* end of chiasma.h */
