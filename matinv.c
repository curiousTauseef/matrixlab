#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "matrix.h"

MATRIX mat_inv(MATRIX a, MATRIX result)
{
    MATRIX A, B, P;
    int	i, n;

    n = MatCol(a);
    if(MatRow(a)!=n)
    {
        return mat_error(MAT_INVERSE_NOT_SQUARE);
    }
    A = mat_copy(a, NULL);
    B = mat_creat(n, 1, UNDEFINED);
    if(result == NULL)result = mat_creat(n, n, UNDEFINED);
    P = mat_creat( n, 1, UNDEFINED );

    if (mat_lu(A, P) == -1)
    {
        mat_free(A);
        mat_free(B);
        mat_free(result);
        mat_free(P);
        return mat_error(MAT_INVERSE_ILL_COND);
    }
    for (i=0; i<n; ++i)
    {
        mat_fill(B, ZERO_MATRIX);
        B[i][0] = 1.0;
        mat_backsubs1(A, B, result, P, i);
    }
    mat_free(A);
    mat_free(B);
    mat_free(P);
    return (result);
}

MATRIX mat_reg_inv(MATRIX a, mtype r_constant, MATRIX result)
{
    int m, n, i;
    MATRIX A = NULL, B = NULL, P = NULL;
    m =MatCol(a);
    n = MatRow(a);
    A = mat_copy(a, NULL);
    if(m!=n) return mat_error(MAT_SIZEMISMATCH);
    for(i=0; i<m; ++i) A[i][i]+= r_constant;


    B = mat_creat(n, 1, UNDEFINED);
    if(result == NULL)result = mat_creat(n, n, UNDEFINED);
    P = mat_creat( n, 1, UNDEFINED );

    if (mat_lu(A, P) == -1)
    {
        mat_free(A);
        mat_free(B);
        mat_free(result);
        mat_free(P);
        return mat_error(MAT_INVERSE_ILL_COND);
    }
    for (i=0; i<n; ++i)
    {
        mat_fill(B, ZERO_MATRIX);
        B[i][0] = 1.0;
        mat_backsubs1(A, B, result, P, i);
    }
    mat_free(A);
    mat_free(B);
    mat_free(P);
    return (result);
}


MATRIX mat_chol(MATRIX a, MATRIX result)
{
    int	i, j, k, n;
    mtype sum = 0.0, p = 0.0;

    n = MatCol(a);
    if(MatRow(a)!=n)
    {
        return mat_error(MAT_INVERSE_NOT_SQUARE);
    }
    result = mat_copy(a, result);

    for (i=0; i<n; ++i)
    {
        for (j=0; j<n; ++i)
        {
            for (sum=result[i][j], k=i-1; k>=0; k--) sum -= result[i][k]*result[j][k];
            {
                if(i==j)
                {
                    if(sum<=0.0) mat_error(MAT_CHOLESKY_FAILED);
                    p = sqrt(sum);
                }
                else result[j][i] =sum/p;
            }
        }
    }
    return (result);
}







