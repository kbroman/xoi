/**********************************************************************
 *
 * GammaIntegral.c
 *
 * copyright (c) 1998-2006, Karl W Broman
 *
 * last modified Nov, 2006
 * first written ~Feb, 1998
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
 * Contains: my_simp, my_trapzd
 *
 **********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "GammaS.h"
#include <R.h>

/*
  integration functions from Numerical Recipes in C
*/


#ifndef USE_RINTEGRATE

/* qsimp */
double my_simp(void func(double *x, int n, void *info), void *info,
               double lo, double hi, double reltol, double abstol,
               int max_subd, int min_subd)
{
    int j;
    double s, st, ost, os;

    ost = os = -1.0e30;

    for(j=0; j < max_subd; j++) {

        st = my_trapzd(lo, hi, j, func, info);

        s = (4.0*st-ost)/3.0;

        if(j >= min_subd && fabs(s-os) < reltol*(fabs(os)+abstol))
            return s;

        os = s;
        ost = st;
    }

    warning("Integration didn't converge.\n");
    return -99.9;
}


double my_trapzd(double a, double b, int nn,
                 void func(double *x, int n, void *info),
                 void *info)
{
    double x, tnm, sum, del, temp[2];
    static double s;
    int it, j;

    if(nn==1) {
        temp[0] = a;
        temp[1] = b;
        func(temp, 2, info);
        return(s=0.5*(b-a)*(temp[0]+temp[1]));
    }
    else {
        for(it=1, j=1; j<nn-1; j++) it <<= 1;
        tnm=it;
        del = (b-a)/tnm;
        x = a+0.5*del;
        for(sum=0.0, j=1; j<=it; j++, x+= del) {
            temp[0] = x;
            func(temp, 1, info);
            sum+=temp[0];
        }
        s = 0.5*(s+(b-a)*sum/tnm);
        return s;
    }
}

#endif

/* end of GammaIntegral.c */
