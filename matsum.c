#include "matrix.h"

/** \brief Computes element-sum of a matrix
 *
 * \param[in] A Input matrix
 * \return \f$ \textrm{sum}( \mathbf{A} ) \f$
 *
 */

mtype mat_sum(MATRIX A)
{
    int i, j, m, n;
    mtype mn = 0.0;
    m = MatCol(A);
    n = MatRow(A);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j) mn += A[i][j];
    }
    return mn;
}

/** \brief Computes row-sum of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \mathbf{A} \mathbf{1} \f$
 *
 */

MATRIX mat_sum_row(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        result[i][0] = 0.0;
        for(j=0; j<m; ++j) result[i][0] += A[i][j];
    }
    return result;
}

/** \brief Computes column-sum of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \mathbf{1}^T \mathbf{A} \f$
 *
 */

MATRIX mat_sum_col(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(1, m, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(j=0; j<m; ++j) result[0][j] = 0.0;
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j) result[0][j] += A[i][j];
    }
    return result;
}

/** \brief Computes element-sum of an integer vector
 *
 * \param[in] A Input integer vector
 * \return \f$ \textrm{sum}( A ) \f$
 *
 */

int int_vec_sum(INT_VECTOR A)
{
    int i, m, mn = 0;
    m = Int_VecLen(A);
    for(i=0; i<m; ++i) mn += A[i];
    return mn;
}

