#include <stdio.h>
#include "matrix.h"


mtype mat_sum( MATRIX A )
{
    int	i, j, m, n;
    mtype mn = 0.0;
    m = MatCol(A);
    n = MatRow(A);

    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j) mn += A[i][j];
    return (mn);
}

MATRIX mat_sum_row( MATRIX A )
{
    int	i, j, m, n;
    MATRIX B;
    m = MatCol(A);
    n = MatRow(A);
    B = mat_creat(n, 1, ZERO_MATRIX);
    if (B==NULL) return NULL;
    #pragma omp parallel for private(j)
    for (j=0; j<m; ++j)
        for (i=0; i<n; ++i) B[i][0] += A[i][j];
    return (B);
}

MATRIX mat_sum_col( MATRIX A )
{
    int	i, j, m, n;
    MATRIX B;
    m = MatCol(A);
    n = MatRow(A);
    B = mat_creat(1, m, ZERO_MATRIX);
    if (B==NULL) return NULL;
    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j) B[0][j] += A[i][j];
    return (B);
}


