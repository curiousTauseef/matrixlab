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
    m = MatCol(A);
    n = MatRow(A);
    o = MatCol(mask);
    p = MatRow(mask);
    if((o%2)!=1 ||(p%2)!=1) gen_error(GEN_SIZE_ERROR);
    ii = (p-1)/2;
    jj = (o-1)/2;
    mm = jj+jj+m;
    nn = ii+ii+n;
    l = jj+m;
    k = ii+n;
    if(scratch==NULL)
    {
        if((scratch = mat_creat(nn, mm, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        flag = 1;
    }
#pragma omp parallel for private(j) firstprivate(mm, ii, jj, k, l)
    for(i=0; i<nn; ++i)
    {
        for(j=0; j<mm; ++j)
        {
            if(i<ii || j<jj || i>=k || j>=l ) scratch[i][j] = 0.0;
            else scratch[i][j] = A[i-ii][j-jj];
        }
    }

    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
#pragma omp parallel for private(j) firstprivate(m, ii, jj, k, l)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            mtype acc = 0.0;
            for(k = -ii; k<=ii; ++k)
            {
                for(l = -jj; l<=jj; ++l)
                {
                    acc += scratch[i+ii+k][j+jj+l]*mask[ii-k][jj-l];
                }
            }
            result[i][j] = acc;
        }
    }
    if(flag==1) mat_free(scratch);
    return result;
}

