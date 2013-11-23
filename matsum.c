#include "matrix.h"


mtype mat_sum(MATRIX A)
{
    int i, j, m, n;
    mtype mn = 0.0;
    m = MatCol(A);
    n = MatRow(A);
    for (i=0; i<n; ++i)
        for(j=0; j<m; ++j) mn += A[i][j];
    return mn;
}

MATRIX mat_sum_row(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, 1, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(j=0; j<m; ++j)
        for(i=0; i<n; ++i) result[i][0] += A[i][j];
    return result;
}

MATRIX mat_sum_col(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(1, m, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j) result[0][j] += A[i][j];
    return result;
}

