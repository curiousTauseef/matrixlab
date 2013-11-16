#include "matrix.h"


MATRIX mat_abs(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result==NULL)if((result = mat_creat( n, m, UNDEFINED )) == NULL)
        return (NULL);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            result[i][j] = (mtype)fabs(A[i][j]);
        }
    return (result);
}

