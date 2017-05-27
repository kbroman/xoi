/**********************************************************************
 *
 * gammaS.h
 *
 * copyright (c) 1998-2013, Karl W Broman
 *
 * last modified Jan, 2013
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
 * Part of the R/xoi package
 *
 * Functions to calculate and maximize the likelihood for the Gamma
 * model, given data on the inter-crossover distances
 *
 **********************************************************************/

#define USE_RINTEGRATE
#include <R_ext/PrtUtil.h>

/* global variables */
struct integr_data {
    int maxsubd, lenw, *iwork;
    double *dwork;
    double reltol, abstol;
    double result, abserr;
    int last, neval, ier;
    int minsubd;
};

struct gamma_data {
    int max_conv, n_length, *type;
    double *length, drop, maxloglik, nu, lenmx, p;
    struct integr_data integr_info;
};

/* from GammaS.c */
void GammaS(int *n_length, double *length, int *type, int *n_nu,
            double *nu, double *loglik, int *max_conv, int *center,
            double *integr_tol, int *maxsubd, int *minsubd);
void GammaMax(int *n_length, double *length, int *type,
              double *low, double *high, double *nu,
              double *loglik, int *max_conv, double *tol,
              double *integr_tol, int *maxsubd, int *minsubd);
void GammaSE(int *n_length, double *length, int *type,
             double *nu, double *se, double *secderiv,
             int *max_conv, double *h, double *hstep,
             double *tol, int *maxit, double *integr_tol, int *maxsubd, int *minsubd);
void GammaInterval(int *n_length, double *length, int *type,
                   double *low, double *high,
                   double *nu, double *interval,
                   double *interval_level, double *drop,
                   int *max_conv, double *tol, int *maxit,
                   double *integr_tol, int *maxsubd, int *minsubd);

/* Optimize function used by R's optimize
   Taken from R ver 2.15.2, in src/appl/fmin.c
*/
double rxoi_Brent_fmin(double ax, double bx, double (*f)(double, void *),
                       void *info, double tol);


/* from GammaUtil.c */
double mydgamma(double x, double shape, double rate);
double mypgamma(double x, double shape, double rate);
double mysgamma(double x, double shape, double rate);
double offenddist(double nu, double length, int max_conv,
                  struct integr_data theintegrdata);
void offenddist_sub(double *x, int n, struct gamma_data *info);
double sumconv(double length, double nu, int max_conv,
               double (*f)(double x, double shape, double rate));
void setup_integr_par(double integr_tol, int maxsubd, int minsubd,
                      struct integr_data *info);
double my_integrate(void f(double *x, int n, void *ex),
                    void *ex, double lo, double hi,
                    struct integr_data info);


/* from GammaLoglik.c */
double ll(double nu, int type, double length, int max_conv,
          struct integr_data theintegrdata);
double calcLL(double nu, struct gamma_data *info);
double calcLLmdrop(double nu, struct gamma_data *info);

/* from GammaDensities.c */
void location_given_one(double *nu, double *x, double *y, int *n, double *L,
                        int *max_conv, double *integr_tol, int *maxsubd, int *minsubd);
void lg1_sub(double *x, int n, struct gamma_data *info);
void first_given_two(double *nu, double *L, double *x,
                     double *y, int *n, int *max_conv,
                     double *integr_tol, int *maxsubd, int *minsubd);
void joint_given_two(double *nu, double *L, double *x, double *y,
                     double *z, int *n, int *max_conv,
                     double *integr_tol, int *maxsubd, int *minsubd);
void distance_given_two(double *nu, double *L, double *x,
                        double *y, int *n, int *max_conv,
                        double *integr_tol, int *maxsubd, int *minsubd);
void distance_given_two_sub(double *x, int n, struct gamma_data *info);
void xoprob(double *nu, double *L, double *pr, int *max_conv,
            double *integr_tol, int *maxsubd, int *minsubd);
void xoprob_subsub(double *x, int n, struct gamma_data *info);
void xoprob_sub(double *x, int n, struct gamma_data *info);
void xoprob_subsub_b(double *x, int n, struct gamma_data *info);
void xoprob_sub_b(double *x, int n, struct gamma_data *info);
void ioden(double *nu, double *x, double *y, int *n, int *max_conv);
void firstden(double *nu, double *x, double *y, int *n, int *max_conv);
void GammaCoincidence(double *nu, double *x, double *y, int *n,
                      int *max_conv);
void StahlCoincidence(double *nu, double *p, double *x, double *y, int *n,
                      int *max_conv);


/* from GammaIntegral.c */
#ifndef USE_RINTEGRATE
double my_simp(void func(double *x, int n, void *info), void *info,
               double lo, double hi, double reltol, double abstol,
               int max_subd, int min_subd);
double my_trapzd(double a, double b, int nn,
                 void func(double *x, int n, void *info),
                 void *info);
#endif

/* from simStahl.c */
void simStahl(int *n_sim, double *nu, double *p, double *L,
              int *nxo, double *loc, int *max_nxo,
              int *n_bins4start);

/* end of gammaS.h */
