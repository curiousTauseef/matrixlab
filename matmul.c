#include <stdio.h>
#include "matrix.h"

#define __min(x,y) (((x)>(y))?(y):(x))

MATRIX mat_mul(MATRIX A, MATRIX B, MATRIX result)
{
    int i, j, k, m, n, o, p;
    m = MatCol(A);
    n = MatRow(A);
    o = MatCol(B);
    p = MatRow(B);
    if(m!=p) return mat_error(MAT_SIZEMISMATCH);
    if(result==NULL) if ((result = mat_creat( n, o, UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j, k)
    for (i=0; i<n; ++i)
        for (j=0; j<o; ++j)
            for (k=0, result[i][j]=0.0; k<m; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
    return (result);
}

MATRIX mat_mul_fast(MATRIX A, MATRIX B, MATRIX result)
{
    int m, n, o, p, hm, hn, ho;
    MATRIX A11 = NULL, A12 = NULL, A21 = NULL, A22 = NULL;
    MATRIX B11 = NULL, B12 = NULL, B21 = NULL, B22 = NULL;
    MATRIX C11 = NULL, C12 = NULL, C21 = NULL, C22 = NULL;
    MATRIX M1 = NULL, M2 = NULL, M3 = NULL, M4 = NULL, M5 = NULL, M6 = NULL, M7 = NULL;
    MATRIX tmp1 = NULL, tmp2 = NULL;
    m = MatCol(A);
    n = MatRow(A);
    o = MatCol(B);
    p = MatRow(B);
    if(m!=p) return mat_error(MAT_SIZEMISMATCH);
    if(result==NULL) if ((result = mat_creat( n, o, ZERO_MATRIX )) == NULL)
            return mat_error(MAT_MALLOC);

    if(m<128)
    {
        return mat_mul(A, B, result);
    }
    hm = m/2;
    hn = n/2;
    ho = o/2;

    A11 = mat_xcopy(A, 0, hn, 0, hm, NULL);
    A12 = mat_xcopy(A, 0, hn, hm, m, NULL);
    A21 = mat_xcopy(A, hn, n, 0, hm, NULL);
    A22 = mat_xcopy(A, hn, n, hm, m, NULL);

    B11 = mat_xcopy(B, 0, hn, 0, ho, NULL);
    B12 = mat_xcopy(B, 0, hn, ho, m, NULL);
    B21 = mat_xcopy(B, hn, n, 0, ho, NULL);
    B22 = mat_xcopy(B, hn, n, ho, m, NULL);

 /*   C11 = mat_mul_fast(A11, B11, NULL);
    tmp1 = mat_mul_fast(A12, B21, NULL);
    C11 = mat_add(C11, tmp1, C11);
    mat_free(tmp1);

    C12 = mat_mul_fast(A11, B12, NULL);
    tmp2 = mat_mul_fast(A12, B22, NULL);
    C12 = mat_add(C12, tmp2, C12);
    mat_free(tmp2);

    C21 = mat_mul_fast(A21, B11, NULL);
    tmp1 = mat_mul_fast(A22, B21, NULL);
    C21 = mat_add(C21, tmp1, C21);
    mat_free(tmp1);

    C22 = mat_mul_fast(A21, B12, NULL);
    tmp2 = mat_mul_fast(A22, B22, NULL);
    C22 = mat_add(C22, tmp2, C22);
    mat_free(tmp2);*/

    tmp1 = mat_add(A11, A22, NULL);
    tmp2 = mat_add(B11, B22, NULL);
    M1 = mat_mul(tmp1, tmp2, NULL);

    tmp1 = mat_add(A21, A22, tmp1);
    M2 = mat_mul(tmp1, B11, NULL);

    tmp1 = mat_sub(B12, B22, tmp1);
    M3 = mat_mul(A11, tmp1, NULL);

    tmp1 = mat_sub(B21, B11, tmp1);
    M4 = mat_mul(A22, tmp1, NULL);

    tmp1 = mat_add(A11, A12, tmp1);
    M5 = mat_mul(tmp1, B22, NULL);

    tmp1 = mat_sub(A21, A11, tmp1);
    tmp2 = mat_add(B11, B12, tmp2);
    M6 = mat_mul(tmp1, tmp2, NULL);

    tmp1 = mat_sub(A12, A22, tmp1);
    tmp2 = mat_add(B21, B22, tmp2);
    M7 = mat_mul(tmp1, tmp2, NULL);

    mat_free(tmp1);
    mat_free(tmp2);

    C11 = mat_add(M1, M4, NULL);
    C11 = mat_sub(C11, M5, C11);
    C11 = mat_add(C11, M7, C11);

    C12 = mat_add(M3, M5, NULL);

    C21 = mat_add(M2, M4, NULL);

    C22 = mat_sub(M1, M2, NULL);
    C22 = mat_add(C11, M3, C22);
    C22 = mat_add(C11, M6, C22);

    result = mat_xjoin(C11, C12, C21, C22, result);

    mat_free(A11);
    mat_free(A12);
    mat_free(A21);
    mat_free(A22);

    mat_free(B11);
    mat_free(B12);
    mat_free(B21);
    mat_free(B22);

    mat_free(C11);
    mat_free(C12);
    mat_free(C21);
    mat_free(C22);

    mat_free(M1);
    mat_free(M2);
    mat_free(M3);
    mat_free(M4);
    mat_free(M5);
    mat_free(M6);
    mat_free(M7);

    return result;
}

MATRIX mat_muls(MATRIX A, mtype s, MATRIX B)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(B== NULL)if ((B = mat_creat( n, m, UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            B[i][j] = A[i][j] * s;
        }
    return (B);
}

MATRIX mat_mul_dot(MATRIX A, MATRIX B, MATRIX C)
{
    int	i, j, m, n, o, p;
    m = MatCol(A);
    n = MatRow(A);
    o = MatCol(B);
    p = MatRow(B);

    if(C== NULL)if ((C = mat_creat( MatRow(A), MatCol(A), UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            if(o==m &&p==n) C[i][j] = A[i][j] * B[i][j];
            else if(o==1 && p!=1) C[i][j] = A[i][j] * B[i][0];
            else if(p==1 && o!=1) C[i][j] = A[i][j] * B[0][j];
            else gen_error(GEN_SIZEMISMATCH);
        }
    return (C);
}

mtype mat_diagmul(MATRIX A)
{
    int i, n;
    mtype result = 1.0;
    n = MatRow(A);
    for (i=0; i<n; ++i)
    {
        result *= A[i][i];
    }
    return (result);
}


INT_VECTOR int_vec_mul(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result == NULL)if ((result = int_vec_creat( m, UNDEFINED )) == NULL)
            int_vec_error(INT_VEC_MALLOC);
    if(m!=Int_VecLen(B))gen_error(GEN_SIZEMISMATCH);
    for (i=0; i<m; ++i) result[i] = A[i] * B[i];
    return (result);
}

INT_VECTOR int_vec_muls(INT_VECTOR A, int x, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result == NULL)if ((result = int_vec_creat( m, UNDEFINED )) == NULL)
            int_vec_error(INT_VEC_MALLOC);
    #pragma omp parallel for
    for (i=0; i<m; ++i) result[i] = A[i] * x;
    return (result);
}


