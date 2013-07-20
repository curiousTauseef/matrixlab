#include <stdio.h>
#include <math.h>
#include "matrix.h"

mtype mat_innerprod(MATRIX A, MATRIX B)
{
    int	i, j, m, n;
    mtype p = 0.0;
    MATRIX	C;
    m = MatCol(A);
    n = MatRow(A);
    C = mat_mul_dot(A, B, NULL);
    if (C == NULL) gen_error(GEN_MATH_ERROR);
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            p+= C[i][j];
        }
    mat_free(C);
    return (p);
}

mtype mat_norm_inf(MATRIX A)
{
    int	i, j, m, n;
    mtype max = 0.0, p ;
    m = MatCol(A);
    n = MatRow(A);
    for (i=0; i<n; i++)
    {
        p = 0.0;
        for (j=0; j<m; ++j)
        {
            p+= (mtype)fabs(A[i][j]);
        }
        if(p>max) max = p;
    }
    return max;
}

mtype mat_norm_one(MATRIX A)
{
    int	i, j, m, n;
    mtype max = 0.0, p ;
    m = MatCol(A);
    n = MatRow(A);
    for (i=0; i<m; ++i)
    {
        p = 0.0;
        for (j=0; j<n; ++j)
        {
            p+= (mtype)fabs(A[j][i]);
        }
        if(p>max) max = p;
    }
    return max;
}

mtype mat_norm_p(MATRIX A, mtype p)
{
    int	i, j, m, n;
    mtype norm_ = 0.0;
    m = MatCol(A);
    n = MatRow(A);
    for (i=0; i<m; ++i)
    {
        for (j=0; j<n; ++j)
        {
            norm_+= (mtype)pow(fabs(A[j][i]), p);
        }
    }
    return (mtype)pow(norm_, (1/p));
}

