/**********************************************************************
 *
 * GammaUtil.c
 *
 * copyright (c) 1998-2006, Karl W Broman
 * 
 * last modified Dec, 2006
 * first written ~Feb, 1998
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * Part of the R/xoi package
 * Contains: lgamma, dgamma, pgamma, sgamma, gser, gcf
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GammaS.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>

/***************************************************
 * mydgamma
 *
 * returns density of gamma dist at x
 ***************************************************/
double mydgamma(double x, double shape, double rate)
{
  return(dgamma(x, shape, 1/rate, 0));
}


/***************************************************
 * mypgamma
 *
 * returns cdf of gamma dist at x
 ***************************************************/
double mypgamma(double x, double shape, double rate)
{
  return(pgamma(x, shape, 1/rate, 1, 0));
}

/***************************************************
 * mysgamma
 *
 * returns 1-cdf of gamma dist at x
 ***************************************************/
double mysgamma(double x, double shape, double rate)
{
  return(pgamma(x, shape, 1/rate, 0, 0));
}

/**********************************************************************
 * offenddist
 * 
 * calculates Pr(no crossovers on chromosome)
 **********************************************************************/
double offenddist(double nu, double length, int max_conv, 
		  struct integr_data theintegrdata)
{
  struct gamma_data info;

  info.nu = nu;
  info.max_conv = max_conv;

  return(1.0 - my_integrate((void (*)(double *, int, void *))offenddist_sub,
			    (void *)(&info), 0.0, length, theintegrdata));
}


/**********************************************************************
 * offenddist_sub
 *
 * calculates sum from k=1 to max_conv of (1/2)^k f(length,k*v,2*v)
 **********************************************************************/
void offenddist_sub(double *x, int n, struct gamma_data *info)
{
  int i;

  for(i=0; i<n; i++) 
    x[i] = sumconv(x[i], (*info).nu, (*info).max_conv, mysgamma);
}



/**********************************************************************
 * sumconv
 *
 * calculates sum from k=1 to max_conv of (1/2)^k f(length,k*v,2*v)
 **********************************************************************/
double sumconv(double length, double nu, int max_conv,
	       double (*f)(double x, double shape, double rate))
{
  int i;
  /*  double a=0.0, b=0.5, c=nu, d=2.0*nu; */

  /*  for(i=0; i<max_conv; i++, b /= 2.0, c += nu) */
  /*    a += b*f(length,c,d); */
  double a=0.0, b;

  b = log(0.5);
  for(i=1; i<=max_conv; i++) 
    a += f(length, (double)i*nu, 2.0*nu) * exp((double)i*b);

  return(a);
}


/**********************************************************************
 * setup_integr_par
 * 
 * set up the parameters for intergration by Rdqags
 * or my integration function 
 **********************************************************************/
void setup_integr_par(double integr_tol, int maxsubd, int minsubd,
		      struct integr_data *info)
{
  (*info).reltol = integr_tol;
  (*info).abstol = integr_tol*10.0;
  (*info).maxsubd = maxsubd;
  (*info).minsubd = minsubd;
#ifdef USE_RINTEGRATE
  (*info).lenw = 4* maxsubd;
  (*info).iwork = (int *)R_alloc(maxsubd, sizeof(int));
  (*info).dwork = (double *)R_alloc(4* maxsubd, sizeof(double));
#endif
}

/**********************************************************************
 * my_integrate
 * 
 * A function to make it easier to use R's Rdqags
 **********************************************************************/
double my_integrate(void f(double *x, int n, void *ex),
		    void *ex, double lo, double hi, 
		    struct integr_data info)
{
  int i;
  double result;

#ifdef USE_RINTEGRATE
  double a, b, step, temp;

  
  step = (hi-lo)/(double)info.minsubd;

  result = 0.0;

  for(i=0; i<info.minsubd; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    a = lo + (double)i*step;
    b = a + step;

    Rdqags(f, ex, &a, &b, &(info.abstol), &(info.reltol), 
	   &temp, &(info.abserr), &(info.neval), &(info.ier),
	   &(info.maxsubd), &(info.lenw), &(info.last), 
	   info.iwork, info.dwork);
    result += temp;
    if(info.ier) 
      warning("Integration error: %d\n", info.ier);
  }
#else
  result = my_simp(f, ex, lo, hi, info.reltol, info.abstol,
		   info.maxsubd, info.minsubd);
#endif

  return(result);

}

/* end of GammaUtil.c */
