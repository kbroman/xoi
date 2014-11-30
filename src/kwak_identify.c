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
 * identify_xo
 *
 * get xo vector and turn this into 3 * (N_xo) matrix with reorg_geno function.
 * After that,
 * xomat[k][0] : observation number that kth xo happened
 * xomat[k][1] : xo starting marker position
 * xomat[k][2] : xo ending marker position.
 *
 **********************************************************************/
void R_identify_xo(int *sdat, int *n_ind, int *n_pos, int *n_xo,
                   int *left, int *right, int *ind_id, int *ob_ind)
{
    identify_xo(sdat, *n_ind, *n_pos, *n_xo, left, right, ind_id, ob_ind);
}

void identify_xo(int *sdat, int n_ind, int n_pos, int n_xo, int *left,
                 int *right, int *ind_id, int *ob_ind)
{
    int i,j;
    int bef = 0;
    int befj = 0;
    int count = 0;
    //  reorg_geno(n_ind, n_pos, sdat, &Geno);

    //  for(k = 0 ; k < n_xo ; k++)
    for( i = 0 ; i < n_ind ; i++ )
        {
            ob_ind[i] = count ;
            for ( j = 0 ; j < n_pos ; j++ )
                {
                    if ( ( bef == 1 && sdat[n_ind * j + i] == 2 ) || ( bef == 2 && sdat[n_ind * j + i] == 1 ) ) // xo happened -> save xo happened position.
                        {
                            //  if(count < n_xo)
                            {
                                ind_id[count] = (i + 1) ;
                                left[count] = befj + 1 ;
                                right[count] = j + 1 ;

                            }
                            count++ ;

                            bef = sdat[n_ind * j + i];
                            befj = j;
                        }

                    if ( bef != 0 && sdat[n_ind * j + i] == bef )  // xo not happened -> former position update
                        befj = j;

                    if ( bef == 0 && ( sdat[n_ind * j + i] == 1 || sdat[n_ind * j + i] == 2 ) ) // save 1st position that have not 0 value.
                        {

                            bef = sdat[n_ind * j + i];
                            befj = j;
                        }

                }
            bef = 0;
            befj = 0;
        }
}
