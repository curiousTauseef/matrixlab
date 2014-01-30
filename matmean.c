#include "matrix.h"


mtype mat_mean(MATRIX A)
{
    int m, n;
    mtype mn;
    m = MatCol(A);
    n = MatRow(A);
    mn = mat_sum(A)/(m*n);
    return mn;
}

MATRIX mat_mean_row(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, 1, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        result[i][0] = 0.0;
        for(j=0; j<m; ++j) result[i][0] += A[i][j];
        result[i][0]/=(mtype)m;
    }
    return result;
}

MATRIX mat_mean_col(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(1, m, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<m; ++i)
    {
        result[0][i] = 0.0;
        for(j=0; j<n; ++j) result[0][i] += A[j][i];
        result[0][i]/=(mtype)n;
    }
    return result;
}

