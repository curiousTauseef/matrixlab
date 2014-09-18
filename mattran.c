#include "matrix.h"


MATRIX mat_tran( MATRIX A, MATRIX result)
{
    int	i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result== NULL)if((result = mat_creat( m,n, UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<m; ++i)
        for (j=0; j<n; ++j)
        {
            result[i][j] = A[j][i];
        }
    return (result);
}

