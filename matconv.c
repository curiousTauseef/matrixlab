#include <stdio.h>
#include "matrix.h"


INT_VECTOR mat_2int_vec(MATRIX a)
{
    int m, n, i, j, l;
    INT_VECTOR v;
    m = MatCol(a);
    n = MatRow(a);
    l = m*n;
    if((v = int_vec_creat(l, UNDEFINED))==NULL) int_vec_error(INT_VECSTACK_MALLOC);
    l = 0;
    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            v[l] = (int)a[i][j];
            ++l;
        }
    return v;
}

MATRIX int_vec2_mat(INT_VECTOR a, int direction)
{
    int n, i;
    MATRIX M;
    n = Int_VecLen(a);
    if(direction ==1)
    {
        if((M = mat_creat(n, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        #pragma omp parallel for
        for (i=0; i<n; ++i) M[i][0] = a[i];
    }
    else
    {
        if((M = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        #pragma omp parallel for
        for (i=0; i<n; ++i) M[0][i] =a[i];
    }
    return M;
}

MATRIX mat_vectorize(MATRIX A, MATRIX result)
{
    int	i, j, m, n, flag = 0, k=0;
    m = MatCol(A);
    n = MatRow(A);
    if(result==A) flag =1;
    if(result== NULL|| result ==A)if ((result = mat_creat( MatNumel(A), 1, UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            result[k][0]= A[i][j];
            ++k;
        }
    if(flag) mat_free(A);
    return (result);
}

