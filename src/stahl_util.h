/**********************************************************************
 *
 * stahl_util.h
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
double addlog_stahl(double a, double b);


/**********************************************************************
 * offenddist_stahl
 * 
 * calculates Pr(no crossovers on chromosome)
 **********************************************************************/
double offenddist_stahl(double nu, double p, double length, int max_conv, 
			struct integr_data theintegrdata);

/**********************************************************************
 * offenddist_stahl_sub
 *
 * calculates sum from k=1 to max_conv of (1/2)^k f(length,k*v,2*p*v)
 **********************************************************************/
void offenddist_stahl_sub(double *x, int n, struct gamma_data *info);

/**********************************************************************
 * sumconv_stahl
 *
 * calculates sum from k=1 to max_conv of (1/2)^k f(length,k*nu,2*p*nu)
 **********************************************************************/
double sumconv_stahl(double length, double nu, double p, int max_conv,
		     double (*f)(double x, double shape, double rate));


double fstar_stahl(double length, double nu, double p, int max_conv);

double oneminus_Fstar_stahl(double length, double nu, double p, int max_conv);

double gstar_stahl(double length, double nu, double p, int max_conv);

double oneminus_Gstar_stahl(double length, double nu, double p, int max_conv,
			    struct integr_data theintegrdata);





/* end of stahl_util.h */
