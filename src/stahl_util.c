/**********************************************************************
 *
 * stahl_util.c
 *
 * copyright (c) 2009, Karl W Broman
 *
 * last modified Jun, 2009
 * first written Jun, 2009
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * Part of the R/xoi package
 * Contains:
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GammaS.h"
#include "stahl_util.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>

#define THRESH 200.0

/**********************************************************************
 *
 * addlog_stahl
 *
 * Calculate addlog(a,b) = log[exp(a) + exp(b)]
 *
 * This makes use of the function log1p(x) = log(1+x) provided
 * in R's math library.
 *
 **********************************************************************/
double addlog_stahl(double a, double b)
{
    if(a < -THRESH) return(b);
    if(b < -THRESH) return(a);
    if(b > a + THRESH) return(b);
    else if(a > b + THRESH) return(a);
    else return(a + log1p(exp(b-a)));
}


/**********************************************************************
 * offenddist_stahl
 *
 * calculates Pr(no crossovers on chromosome)
 **********************************************************************/
double offenddist_stahl(double nu, double p, double length, int max_conv,
                        struct integr_data theintegrdata)
{
    struct gamma_data info;

    info.p = p;
    info.nu = nu;
    info.max_conv = max_conv;

    return(1.0 - my_integrate((void (*)(double *, int, void *))offenddist_stahl_sub,
                              (void *)(&info), 0.0, length, theintegrdata));
}


/**********************************************************************
 * offenddist_stahl_sub
 *
 * calculates sum from k=1 to max_conv of (1/2)^k f(length,k*v,2*p*v)
 **********************************************************************/
void offenddist_stahl_sub(double *x, int n, struct gamma_data *info)
{
    int i;

    for(i=0; i<n; i++)
        x[i] = (*info).p * sumconv_stahl(x[i], (*info).nu, (*info).p, (*info).max_conv,
                                         mysgamma);
}



/**********************************************************************
 * sumconv_stahl
 *
 * calculates sum from k=1 to max_conv of (1/2)^k f(length,k*nu,2*p*nu)
 **********************************************************************/
double sumconv_stahl(double length, double nu, double p, int max_conv,
                     double (*f)(double x, double shape, double rate))
{
    int i;
    double a=0.0, b;

    b = log(0.5);
    for(i=1; i<=max_conv; i++)
        a += f(length, (double)i*nu, 2.0*p*nu) * exp((double)i*b);

    return(a);
}


double fstar_stahl(double length, double nu, double p, int max_conv)
{
    return(sumconv_stahl(length, nu, p, max_conv, mydgamma));
}


double oneminus_Fstar_stahl(double length, double nu, double p, int max_conv)
{
    if(p < 1e-12) return(1.0);
    return(sumconv_stahl(length, nu, p, max_conv, mysgamma));
}

double gstar_stahl(double length, double nu, double p, int max_conv)
{
    return(p*oneminus_Fstar_stahl(length, nu, p, max_conv));
}

double oneminus_Gstar_stahl(double length, double nu, double p, int max_conv,
                            struct integr_data theintegrdata)
{
    if(p < 1e-12) return(1.0);
    return(offenddist_stahl(nu, p, length, max_conv, theintegrdata));
}




/* end of stahl_util.c */
