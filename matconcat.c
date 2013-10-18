#include <stdio.h>
#include "matrix.h"


MATRIX mat_concat(MATRIX A, MATRIX B, int dim)
{
    int i, j, m, n, o, p;
    MATRIX C;
    if(A == NULL)
    {
        return mat_copy(B, NULL);
    }
    else
    {
        m = MatCol(A);
        n = MatRow(A);
    }
    if(B==NULL)
    {
        return mat_copy(A, NULL);
    }
    else
    {
        o = MatCol(B);
        p = MatRow(B);
    }
    if((dim==1)&&((m==o) ||!((m==0)&&(o==0))))
    {
        if((C = mat_creat(n+p, m, UNDEFINED))==NULL) return NULL;
        #pragma omp parallel for private(j)
        for(i=0; i<m; i++)
        {
            for (j = 0; j<n; j++) C[j][i] = A[j][i];
            for (j = 0; j<p; j++) C[j+n][i] = B[j][i];
        }
        return C;
    }
    if((dim==2)&&((n==p) ||!((n==0)&&(p==0))))
    {
        if((C = mat_creat(n, m+o, UNDEFINED))==NULL) return NULL;
        #pragma omp parallel for private(j)
        for(i=0; i<n; i++)
        {
            for (j = 0; j<m; j++) C[i][j] = A[i][j];
            for (j = 0; j<o; j++) C[i][j+m] = B[i][j];
        }
        return C;
    }
    return mat_error(MAT_SIZEMISMATCH);
}


INT_VECTOR int_vec_concat(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result)
{
    int i, m, n;
    m = Int_VecLen(A);
    n = Int_VecLen(B);
    if(result==NULL)
    {
        if((result = int_vec_creat(m+n, UNDEFINED))==NULL) return NULL;
    }
    else
    {
        if(Int_VecLen(result)!=(m+n)) return int_vec_error(INT_VEC_SIZEMISMATCH);
    }
    #pragma omp parallel for
    for(i=0; i<m; ++i) result[i] = A[i];
    for(i=0; i<n; ++i) result[i+m] = B[i];
    return result;
}


