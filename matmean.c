#include <stdio.h>
#include "matrix.h"


mtype mat_mean( MATRIX A)
{
    int m, n;
    mtype mn;
    m = MatCol(A);
    n = MatRow(A);

    mn = mat_sum(A)/(m*n);
    return (mn);
}

MATRIX mat_mean_row( MATRIX A )
{
    int	i, j, m, n;
    MATRIX B;
    m = MatCol(A);
    n = MatRow(A);
    B = mat_creat(n, 1, ZERO_MATRIX);
    if (B==NULL) return NULL;
    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
    {
        for (j=0; j<m; ++j) B[i][0] += A[i][j];
        B[i][0] /= (mtype)m;
    }
    return (B);
}

MATRIX mat_mean_col( MATRIX A )
{
    int	i, j, m, n;
    MATRIX B;
    m = MatCol(A);
    n = MatRow(A);
    B = mat_creat(1, m, ZERO_MATRIX);
    if (B==NULL) return NULL;
    #pragma omp parallel for private(j)
    for (i=0; i<m; ++i)
    {
        for (j=0; j<n; ++j) B[0][i] += A[j][i];
        B[0][i] /= (mtype)n;
    }
    return (B);
}


