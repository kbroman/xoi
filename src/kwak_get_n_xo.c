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
 * reorg_geno
 *
 * Reorganize the genotype data so that it is a doubly indexed array
 * rather than a single long vector
 *
 * Afterwards, geno indexed like Geno[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
    int i;

    *Geno = (int **)R_alloc(n_pos, sizeof(int *));

    (*Geno)[0] = geno;
    for(i=1; i< n_pos; i++)
        (*Geno)[i] = (*Geno)[i-1] + n_ind;

}


/**********************************************************************
 *
 * get_N_xo
 *
 * Count the # of xo happened.
 * By counting this number, we assign memory to xo_vec and xo_matrix.
 *
 **********************************************************************/
void R_get_N_xo(int *n_ind, int *n_pos, int *dat, int *n_xo)
{
    int **Geno;
    reorg_geno(*n_ind, *n_pos, dat, &Geno);
    *n_xo = get_N_xo(*n_ind, *n_pos, Geno);

}

int get_N_xo(int n_ind, int n_pos, int **Geno)
{
    int i,j;
    int bef = 0;
    int n_xo = 0;

    for( i = 0 ; i < n_ind ; i++ )
        {
            for ( j = 0 ; j < n_pos ; j++ )
                {
                    if ( ( bef == 1 && Geno[j][i] == 2 ) || ( bef == 2 && Geno[j][i] == 1 ) )
                        {
                            n_xo++;
                            bef = Geno[j][i];
                        }

                    if ( bef == 0 && ( Geno[j][i] == 1 || Geno[j][i] == 2 ) )
                        bef = Geno[j][i];

                }
            bef = 0;
        }
    return( n_xo );
}
