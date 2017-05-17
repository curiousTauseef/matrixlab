#include "matrix.h"


/** \brief Computes the mean of a matrix
 *
 * \param A Input matrix
 * \return \f$ \textrm{mean}(\mathbf{A}) \f$
 *
 */

mtype mat_mean(MATRIX A)
{
    int m, n;
    mtype mn;
    m = MatCol(A);
    n = MatRow(A);
    mn = mat_sum(A)/(m*n);
    return mn;
}

/** \brief Computes row-mean of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \mathbf{A} \mathbf{1}/\textrm{\#cols} \f$
 *
 */

MATRIX mat_mean_row(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, 1, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        result[i][0] = 0.0;
        for(j=0; j<m; ++j) result[i][0] += A[i][j];
        result[i][0]/=(mtype)m;
    }
    return result;
}

/** \brief Computes column-mean of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \mathbf{1}^T \mathbf{A}/\textrm{\#rows} \f$
 *
 */

MATRIX mat_mean_col(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(1, m, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    for(j=0; j<m; ++j) result[0][j] = 0.0;
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j) result[0][j] += A[i][j];
    }
    for(j=0; j<m; ++j) result[0][j]/=(mtype)n;
    return result;
}

/** \brief Computes element-mean of an integer vector
 *
 * \param[in] A Input integer vector
 * \return \f$ \textrm{mean}( A ) \f$
 *
 */

mtype int_vec_mean(INT_VECTOR A)
{
    int i, m, mn = 0;
    m = Int_VecLen(A);
    for(i=0; i<m; ++i) mn += A[i];
    return mn/(mtype)m;
}

