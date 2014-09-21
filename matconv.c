#include "matrix.h"


/** \brief Converts a matrix to an integer vector
 *
 * \param[in] A Input matrix
 * \param[out] v Output vector
 * \return Output vector
 *
 */

INT_VECTOR mat_2int_vec(MATRIX A)
{
    int m, n, i, j, l;
    INT_VECTOR v;
    m = MatCol(A);
    n = MatRow(A);
    l = m*n;
    if((v = int_vec_creat(l, UNDEFINED))==NULL) int_vec_error(INT_VECSTACK_MALLOC);
    l = 0;
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            v[l] = (int)A[i][j];
            ++l;
        }
    return v;
}

/** \brief Converts an integer vector to a matrix
 *
 * \param[in] a Input vector
 * \param[in] dir Conversion direction
 * \return Output matrix
 *
 */

MATRIX int_vec2_mat(INT_VECTOR a, int dir)
{
    int n, i;
    MATRIX M;
    n = Int_VecLen(a);
    if(dir==0)
    {
        if((M = mat_creat(n, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        #pragma omp parallel for
        for(i=0; i<n; ++i) M[i][0] = a[i];
    }
    else
    {
        if((M = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        #pragma omp parallel for
        for(i=0; i<n; ++i) M[0][i] =a[i];
    }
    return M;
}

/** \brief Reshapes a matrix to a vector
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ vec(\mathbf{A}) \f$
 *
 */

MATRIX mat_vectorize(MATRIX A, MATRIX result)
{
    int i, j, m, n, flag = 0, k=0;
    m = MatCol(A);
    n = MatRow(A);
    if(result==A) flag =1;
    if(result==NULL || result==A)if((result = mat_creat( MatNumel(A), 1, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            result[k][0]= A[i][j];
            ++k;
        }
    if(flag) mat_free(A);
    return result;
}

/** \brief Reshapes transpose of a matrix to a vector
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ vec(\mathbf{A}^T) \f$
 *
 */

MATRIX mat_vectorize_tr(MATRIX A, MATRIX result)
{
    int i, j, m, n, flag = 0, k=0;
    m = MatCol(A);
    n = MatRow(A);
    if(result==A) flag =1;
    if(result==NULL || result==A)if((result = mat_creat(1, MatNumel(A), UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            result[0][k] = A[i][j];
            ++k;
        }
    if(flag) mat_free(A);
    return result;
}

