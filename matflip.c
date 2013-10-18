#include <stdio.h>
#include "matrix.h"

MATRIX mat_fliplr(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result == NULL) if((result = mat_creat( MatRow(A), MatCol(A), UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            result[i][j] = A[i][m-j-1];
        }
    return (result);
}

MATRIX mat_flipud(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result == NULL) if((result = mat_creat( MatRow(A), MatCol(A), UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            result[i][j] = A[n-i-1][j];
        }
    return result;
}

