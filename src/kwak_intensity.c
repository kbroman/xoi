#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "kwak.h"


/**********************************************************************
 *
 * get_intensity
 *
 * get intensity probability for each divided marker interval.
 *
 **********************************************************************/
void R_get_intensity(int *xovec, double *window, double *center,
                     int *n_pos, int *n_xo,  int *n_center,
                     double *marker, double *intensity, int *n_ind)
{
    get_intensity(xovec, *window, center, *n_pos, *n_xo, *n_center, marker, intensity, *n_ind);
}

void get_intensity(int *xovec, double window, double *center,
                   int n_pos, int n_xo, int n_center, double *marker,
                   double *intensity, int n_ind)
{
    int i, j, ob;
    double weight, prob;

    for(i = 0 ; i < n_center ; i++)
        {
            ob = xovec[0];
            prob = 0;
            for(j = 0 ; j < n_xo ; j++)
                {
                    if(ob == xovec[j*3])
                        {
                            if( ( marker[xovec[j*3+1]-1] <= center[i] + window/2 &&
                                  marker[xovec[j*3+1]-1] >= center[i] - window/2 ) ||
                                ( marker[xovec[j*3+2]-1] <= center[i] + window/2 &&
                                  marker[xovec[j*3+2]-1] >= center[i] - window/2 ) ||
                                ( marker[xovec[j*3+1]-1] <= center[i] - window/2 &&
                                  marker[xovec[j*3+2]-1] >= center[i] + window/2 )  )
                                {
                                    weight = MIN(center[i]+window/2, marker[xovec[j*3+2]-1]) - MAX(center[i]-window/2, marker[xovec[j*3+1]-1]) ;
                                    prob += weight / ( marker[xovec[j*3+2]-1] - marker[xovec[j*3+1]-1] );
                                    //    Rprintf("j : %d, (%f,%f) ob : %d,  weight : %f,  prob : %f \n", j, center[i] - window/2, center[i] + window/2, ob, weight, prob);
                                }
                        } else {
                        if (prob > 1) prob = 1;
                        intensity[i] = intensity[i] +  prob;
                        //    Rprintf("(%f,%f) ob : %d,  intensity[%i] : %f \n", center[i] - window/2, center[i] + window/2, ob, i, intensity[i]);
                        ob = xovec[j*3];
                        prob = 0;
                        if( ( marker[xovec[j*3+1]-1] <= center[i] + window/2 &&
                              marker[xovec[j*3+1]-1] >= center[i] - window/2 ) ||
                            ( marker[xovec[j*3+2]-1] <= center[i] + window/2 &&
                              marker[xovec[j*3+2]-1] >= center[i] - window/2 ) ||
                            ( marker[xovec[j*3+1]-1] <= center[i] - window/2 &&
                              marker[xovec[j*3+2]-1] >= center[i] + window/2 )        )
                            {
                                weight = MIN(center[i]+window/2, marker[xovec[j*3+2]-1]) - MAX(center[i]-window/2, marker[xovec[j*3+1]-1]) ;
                                prob += weight / ( marker[xovec[j*3+2]-1] - marker[xovec[j*3+1]-1] );
                                //    Rprintf("(%f,%f) ob : %d,  weight : %f,  prob : %f \n", center[i] - window/2, center[i] + window/2, ob, weight, prob);
                                if (j==(n_xo-1))
                                    {
                                        if (prob > 1) prob = 1;
                                        intensity[i] = intensity[i] +  prob;
                                        //        Rprintf("j : %d (%f,%f) ob : %d,  intensity[%i] : %f \n", j, center[i] - window/2, center[i] + window/2, ob, i, intensity[i]);
                                    }
                            }
                    }
                }

            /* contained window */
            weight = MIN(center[i]+window/2, marker[n_pos-1]) - MAX(center[i]-window/2, marker[0]);

            intensity[i] /= (weight*(double)n_ind/100.0);

        } /* loop over center */
}
