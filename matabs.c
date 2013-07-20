#include <stdio.h>
#include <math.h>
#include "matrix.h"

MATRIX mat_abs(MATRIX A)
{
    int	i, j, m, n;
    MATRIX B;
    m = MatCol(A);
    n = MatRow(A);

    if ((B = mat_creat( n, m, UNDEFINED )) == NULL)
        return (NULL);

    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            B[i][j] = (mtype)fabs(A[i][j]);
        }
    return (B);
}


