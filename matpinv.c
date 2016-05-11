#include "matrix.h"
#include <stdlib.h>


/** \brief Computes pseudo-inverse of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \left( A^TA\right)^{-1}A^T \f$
 *
 */

MATRIX mat_pinv(MATRIX A, MATRIX result)
{
    MATRIX D, T;
    T = mat_tran(A, NULL);
    D = mat_mul(T, A, NULL);
    D = mat_inv(D, D);
    if(D==NULL) return mat_error(MAT_INVERSE_ILL_COND);
    result = mat_mul(D, T, result);
    mat_free(T);
    mat_free(D);
    return result;
}

/** \brief Computes weighted pseudo-inverse of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] w Weight matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \left( A^TWA\right)^{-1}A^TW \f$
 *
 */

MATRIX mat_wpinv(MATRIX A, MATRIX w, MATRIX result)
{
    MATRIX D, B, T;
    T = mat_tran(A, NULL);
    D = mat_mul(T, w, NULL);
    B = mat_mul(D, A, NULL);
    mat_free(T);
    B = mat_inv(B, B);
    if(B==NULL) return mat_error(MAT_INVERSE_ILL_COND);
    result = mat_mul(B, D, result);
    mat_free(D);
    mat_free(B);
    return result;
}

