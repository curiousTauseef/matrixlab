#include <stdio.h>
#include "matrix.h"



MATVEC_DPOINTER mat_max(MATRIX A, int dim)
{
    int m, n, i, j;
    MATVEC_DPOINTER p = NULL;
    MATRIX B;
    INT_VECTOR indices;
    p = matvec_creat();
    m = MatCol(A);
    n = MatRow(A);
    if (dim==1 && n>1)
    {
        if((B = mat_creat(1, m, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        if((indices = int_vec_creat(m, UNDEFINED))==NULL) int_vec_error(INT_VEC_MALLOC);
        for (i = 0; i<m; ++i)
        {
            B[0][i] = A[0][i];
            indices[i] = 0;
            for (j = 1; j<n; ++j)
                if(A[j][i]>B[0][i])
                {
                    B[0][i] = A[j][i];
                    indices[i] = j;
                }
        }
        p[0] = B;
        p[1] = indices;
        return p;
    }
    if (dim==2 && m>1)
    {
        if((B = mat_creat(n, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        if((indices = int_vec_creat(n, UNDEFINED))==NULL) int_vec_error(INT_VEC_MALLOC);
        for (i = 0; i<n; ++i)
        {
            B[i][0] = A[i][0];
            indices[i] = 0;
            for (j = 1; j<m; ++j)
                if(A[i][j]>B[i][0])
                {
                    B[i][0] = A[i][j];
                    indices[i] = j;
                }
        }
        p[0] = B;
        p[1] = indices;
        return p;
    }
    mat_error (MAT_SIZEMISMATCH);
    return  NULL;
}

MATVEC_DPOINTER mat_min(MATRIX A, int dim)
{
    int m, n, i, j;
    MATVEC_DPOINTER p = NULL;
    MATRIX B;
    INT_VECTOR indices;
    p = matvec_creat();
    m = MatCol(A);
    n = MatRow(A);
    if (dim==1 && n>1)
    {
        if((B = mat_creat(1, m, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        if((indices = int_vec_creat( m, UNDEFINED))==NULL) int_vec_error(INT_VEC_MALLOC);
        for (i = 0; i<m; ++i)
        {
            B[0][i] = A[0][i];
            indices[i] = 0;
            for (j = 1; j<n; ++j)
                if(A[j][i]<B[0][i])
                {
                    B[0][i] = A[j][i];
                    indices[i] = j;
                }
        }
        p[0] = B;
        p[1] = indices;
        return p;
    }
    if (dim==2 && m>1)
    {
        if((B = mat_creat(n, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        if((indices = int_vec_creat(n, UNDEFINED))==NULL) int_vec_error(INT_VEC_MALLOC);
        for (i = 0; i<n; ++i)
        {
            B[i][0] = A[i][0];
            indices[i] = 0;
            for (j = 1; j<m; ++j)
                if(A[i][j]<B[i][0])
                {
                    B[i][0] = A[i][j];
                    indices[i] = j;
                }
        }
        p[0] = B;
        p[1] = indices;
        return p;
    }
    mat_error (MAT_SIZEMISMATCH);
    return NULL;
}


