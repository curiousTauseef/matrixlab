#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "matrix.h"


MATRIX mat_pinv(MATRIX a, MATRIX result)
{
    MATRIX	D, T;
    T = mat_tran(a, NULL);
    D = mat_mul(T, a, NULL);
    D = mat_inv(D, D);
    if(D==NULL) return mat_error(MAT_INVERSE_ILL_COND);
    result = mat_mul(D, T, result);
    mat_free(T);
    mat_free(D);
    return result;
}

MATRIX mat_wpinv(MATRIX a, MATRIX w, MATRIX result)
{
    MATRIX	D, B, T;
    T = mat_tran(a, NULL);
    D = mat_mul(T, w, NULL);
    B = mat_mul(D, a, NULL);
    mat_free(T);
    B = mat_inv(B, B);
    if(B==NULL) return mat_error(MAT_INVERSE_ILL_COND);
    result = mat_mul(B, D, result);
    mat_free(D);
    mat_free(B);
    return result;
}

