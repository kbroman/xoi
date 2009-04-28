/**********************************************************************
 * 
 * GammaLoglik.c
 * 
 * copyright (c) 1999-2006, Karl W Broman
 * 
 * last modified Dec, 2006
 * first written ~Jan, 2001
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * Part of the R/xoi package
 * Contains: ll, sumconv, integrand
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GammaS.h"
#include <R_ext/Utils.h>

/* log likelihood function for a single point */
double ll(double nu, int type, double length, int max_conv, 
	  struct integr_data theintegrdata)
{
  if(type==0) /* interior segment on chromosome */
    return(log(sumconv(length, nu, max_conv, mydgamma)));

  else if(type==1 || type==2) /* first or last segment on chromosome */
    return(log(sumconv(length, nu, max_conv, mysgamma)));

  else /* whole chromosome */ {
    return(log(offenddist(nu, length, max_conv, theintegrdata)));
  }
}

/**********************************************************************
 * calcLL
 * 
 * calculate negative log likelihood 
 **********************************************************************/
double calcLL(double nu, struct gamma_data *info)
{
  int i;
  double a=0.0;

  for(i=0; i < (*info).n_length; i++) {
    R_CheckUserInterrupt(); /* check for ^C */

    a -= ll(nu, (*info).type[i], (*info).length[i], (*info).max_conv,
	    (*info).integr_info);
  }

  return(a);
}

/**********************************************************************
 * calcLLmdrop
 *
 * calculates difference between neg lik and maximum - drop
 **********************************************************************/
double calcLLmdrop(double nu, struct gamma_data *info)
{
  return(((*info).maxloglik - (*info).drop) + calcLL(nu, info));
}

/* end of GammaLoglik.c */
