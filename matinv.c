#include "matrix.h"
#include <stdlib.h>


/** \brief Computes the inverse of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ A^{-1} \f$
 *
 */

MATRIX mat_inv(MATRIX A, MATRIX result)
{
    MATRIX a, B, P;
    int i, n;
    n = MatCol(A);
    if(MatRow(A)!=n) return mat_error(MAT_INVERSE_NOT_SQUARE);
    a = mat_copy(A, NULL);
    B = mat_creat(n, 1, UNDEFINED);
    if(result==NULL) if((result = mat_creat(n, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    P = mat_creat(n, 1, UNDEFINED);
    if(mat_lu(a, P)==-1)
    {
        mat_free(a);
        mat_free(B);
        mat_free(result);
        mat_free(P);
        return mat_error(MAT_INVERSE_ILL_COND);
    }
    for(i=0; i<n; ++i)
    {
        mat_fill_type(B, ZERO_MATRIX);
        B[i][0] = 1.0;
        mat_backsubs1(a, B, result, P, i);
    }
    mat_free(a);
    mat_free(B);
    mat_free(P);
    return result;
}

/** \brief Computes the regularized inverse of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] r Regularizing constant
 * \param[in] result Matrix to store the result
 * \return \f$ \left(A+rI\right)^{-1} \f$
 *
 */

MATRIX mat_reg_inv(MATRIX A, mtype r, MATRIX result)
{
    int m, n, i;
    MATRIX a = NULL, B = NULL, P = NULL;
    m = MatCol(A);
    n = MatRow(A);
    a = mat_copy(A, NULL);
    if(m!=n) return mat_error(MAT_SIZEMISMATCH);
    for(i=0; i<m; ++i) A[i][i]+= r;
    B = mat_creat(n, 1, UNDEFINED);
    if(result==NULL) if((result = mat_creat(n, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    P = mat_creat(n, 1, UNDEFINED);
    if(mat_lu(a, P)==-1)
    {
        mat_free(a);
        mat_free(B);
        mat_free(result);
        mat_free(P);
        return mat_error(MAT_INVERSE_ILL_COND);
    }
    for(i=0; i<n; ++i)
    {
        mat_fill_type(B, ZERO_MATRIX);
        B[i][0] = 1.0;
        mat_backsubs1(a, B, result, P, i);
    }
    mat_free(a);
    mat_free(B);
    mat_free(P);
    return result;
}


