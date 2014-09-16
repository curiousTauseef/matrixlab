#include "matrix.h"


/** \brief Computes 2-D convolution
 *
 * \param[in] A Input matrix
 * \param[in] mask Input kernel/mask
 * \param[in] scratch Scratch matrix for temporary calculations
 * \param[in] result Matrix to store the result
 * \return Convolved output matrix
 *
 */

MATRIX mat_conv2(MATRIX A, MATRIX mask, MATRIX scratch, MATRIX result)
{
    int i, j, k, l, m, n, o, p, ii, jj, mm, nn, flag = 0;
    mtype acc = 0.0;
    m = MatCol(A);
    n = MatRow(A);
    o = MatCol(mask);
    p = MatRow(mask);
    if((o%2)!=1 ||(p%2)!=1) gen_error(GEN_SIZE_ERROR);
    ii = (o - 1)/2;
    jj = (p - 1)/2;
    mm = ii+ii+m;
    nn = jj+jj+n;
    k = ii+m;
    l = jj+n;
    if(scratch == NULL)
    {
        if((scratch = mat_creat(nn, mm, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        flag = 1;
        for(i=0; i<mm; ++i)
        {
            for(j=0; j<nn; ++j)
            {
                if(i<ii || j<jj || i>=k || j>=l ) scratch[j][i] = 0.0;
                else scratch[j][i] = A[j-jj][i-ii];
            }
        }
    }

    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(i=0; i<m; ++i)
    {
        for(j=0; j<n; ++j)
        {
            acc = 0.0;
            for(k = -ii; k<=ii; ++k)
                for(l = -jj; l<=jj; ++l)
                    acc += scratch[j+jj+l][i+ii+k]*mask[jj-l][ii-k];
            result[j][i] = acc;
        }
    }
    if(flag==1) mat_free(scratch);
    return result;
}

