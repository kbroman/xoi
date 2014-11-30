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
 * get_coincidence
 *
 * get coincidence function
 *
 **********************************************************************/

void R_get_coincidence(int *xovec, double *int_dat, double *window,
                       double *center, int *n_xo, int *n_pos,
                       int *n_center, int *start_d, double *marker,
                       double *coincidence)
{
    get_coincidence(xovec, int_dat, *window, center, *n_xo, *n_pos, *n_center, *start_d, marker, coincidence);
}

void get_coincidence(int *xovec, double *int_dat, double window,
                     double *center, int n_xo, int n_pos, int n_center,
                     int start_d, double *marker, double *coincidence)
{
    int h, i, j, ob, n_coinc, dino;
    double *coinc1, *coinc2, *coincidence2;
    double left_prob, right_prob;
    double weight;

    n_coinc = n_center - start_d - 1;
    coinc1 = (double *) R_alloc ( (n_coinc * n_center), sizeof(double) );
    coinc2 = (double *) R_alloc ( (n_coinc * n_center), sizeof(double) );
    coincidence2 = (double *) R_alloc ( n_center, sizeof(double) );

    for(h = 0 ; h < n_coinc * n_center; h++ )
        {
            coinc1[h] = 0;
            coinc2[h] = 0;
        }

    for(h = 0 ; h < n_center ; h++ )
        coincidence2[h] = 0 ;


    for(h = 0 ; h < n_coinc ; h++)
        {
            //      Rprintf("%d ", h);
            for(i = h + start_d ; i < n_center ; i++)
                {
                    ob = xovec[0];
                    left_prob = 0;
                    right_prob = 0;
                    for(j = 0 ; j < n_xo ; j++ )
                        {
                            if(ob == xovec[j*3])
                                {
                                    if(  ( marker[xovec[j*3+1]-1] < center[h] + window/2 &&
                                           marker[xovec[j*3+1]-1] > center[h] - window/2 ) ||
                                         ( marker[xovec[j*3+2]-1] < center[h] + window/2 &&
                                           marker[xovec[j*3+2]-1] > center[h] - window/2 ) ||
                                         ( marker[xovec[j*3+1]-1] < center[h] - window/2 &&
                                           marker[xovec[j*3+2]-1] > center[h] + window/2 )        )
                                        {
                                            weight = MIN(center[h]+window/2, marker[xovec[j*3+2]-1]) - MAX(center[h]-window/2, marker[xovec[j*3+1]-1]) ;
                                            left_prob = left_prob + weight / ( marker[xovec[j*3+2]-1] - marker[xovec[j*3+1]-1] );
                                            if ( marker[xovec[j*3+2]-1] > center[i] - window/2 )
                                                left_prob = 0; // in this case, xo interval is too large,
                                        }

                                    if( ( marker[xovec[j*3+1]-1] < center[i] + window/2 &&
                                          marker[xovec[j*3+1]-1] > center[i] - window/2 ) ||
                                        ( marker[xovec[j*3+2]-1] < center[i] + window/2 &&
                                          marker[xovec[j*3+2]-1] > center[i] - window/2 ) ||
                                        ( marker[xovec[j*3+1]-1] < center[i] - window/2 &&
                                          marker[xovec[j*3+2]-1] > center[i] + window/2 )     )
                                        {
                                            weight = MIN(center[i]+window/2, marker[xovec[j*3+2]-1]) - MAX(center[i]-window/2, marker[xovec[j*3+1]-1]) ;
                                            right_prob = right_prob + weight / ( marker[xovec[j*3+2]-1] - marker[xovec[j*3+1]-1] ) ;
                                            //                Rprintf("right_prob : %f  ", right_prob);
                                            if ( marker[xovec[j*3+1]-1] < center[h] + window/2 )
                                                right_prob = 0; // in this case, xo interval is too large,
                                        }

                                } else {
                                if( left_prob != 0 && right_prob != 0 )
                                    {
                                        if (left_prob > 1) left_prob = 1;
                                        if (right_prob > 1) right_prob = 1;
                                        coinc1[h*n_center + i] += ( left_prob * right_prob ) ;
                                        coinc2[h*n_center + i] = ( int_dat[h] * int_dat[i] ) ;
                                        //            Rprintf("%f ", coinc[h*n_center + i]);
                                    }
                                ob = xovec[j*3];
                                left_prob = 0;
                                right_prob = 0;

                                if( ( marker[xovec[j*3+1]-1] < center[h] + window/2 &&
                                      marker[xovec[j*3+1]-1] > center[h] - window/2 ) ||
                                    ( marker[xovec[j*3+2]-1] < center[h] + window/2 &&
                                      marker[xovec[j*3+2]-1] > center[h] - window/2 ) ||
                                    ( marker[xovec[j*3+1]-1] < center[h] - window/2 &&
                                      marker[xovec[j*3+2]-1] > center[h] + window/2 )  )
                                    {
                                        weight = MIN(center[h]+window/2, marker[xovec[j*3+2]-1]) - MAX(center[h]-window/2, marker[xovec[j*3+1]-1]) ;
                                        left_prob = left_prob + weight / ( marker[xovec[j*3+2]-1] - marker[xovec[j*3+1]-1] );
                                        if ( marker[xovec[j*3+2]-1] > center[i] - window/2 )
                                            left_prob = 0; // in this case, xo interval is too large,
                                    }

                                if( ( marker[xovec[j*3+1]-1] < center[i] + window/2 &&
                                      marker[xovec[j*3+1]-1] > center[i] - window/2 ) ||
                                    ( marker[xovec[j*3+2]-1] < center[i] + window/2 &&
                                      marker[xovec[j*3+2]-1] > center[i] - window/2 ) ||
                                    ( marker[xovec[j*3+1]-1] < center[i] - window/2 &&
                                      marker[xovec[j*3+2]-1] > center[i] + window/2 )         )
                                    {
                                        weight = MIN(center[i]+window/2, marker[xovec[j*3+2]-1]) - MAX(center[i]-window/2, marker[xovec[j*3+1]-1]) ;
                                        right_prob = right_prob + weight / ( marker[xovec[j*3+2]-1] - marker[xovec[j*3+1]-1] ) ;
                                        if ( marker[xovec[j*3+1]-1] < center[h] + window/2 )
                                            right_prob = 0; // in this case, xo interval is too large,

                                    }
                            }
                        }
                }
        }

    for(h = 0 ; h < n_coinc ; h++)
        {
            for(i = h + start_d ; i < n_center ; i++)
                {
                    coinc2[h*n_center + i] = ( int_dat[h] * int_dat[i] ) ;
                }
        }

    for(i = 0 ; i < n_center ; i++)
        {
            dino = MIN(n_center - i, n_coinc);
            for(j = 0 ; j < dino; j++)
                {
                    coincidence[i] += coinc1[i+j*n_center + j];
                    coincidence2[i] += coinc2[i+j*n_center + j];
                }

            coincidence[i] = coincidence[i] / coincidence2[i];
        }
}
