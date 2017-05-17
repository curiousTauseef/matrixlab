#include "matrix.h"


/** \brief Computes absolute value of matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \textrm{abs}(\mathbf{A}) \f$
 *
 */

MATRIX mat_abs(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
        return (NULL);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            result[i][j] = (mtype)fabs(A[i][j]);
        }
    return (result);
}

/** \brief Computes absolute value of an integer vector
 *
 * \param[in] A Input integer vector
 * \param[in] result Vector to store the result
 * \return \f$ \textrm{abs}(A) \f$
 *
 */

INT_VECTOR int_vec_abs(INT_VECTOR A, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result==NULL) if((result = int_vec_creat(m, UNDEFINED))==NULL)
            int_vec_error(INT_VEC_MALLOC);
    for(i=0; i<m; ++i) result[i] = abs(A[i]);
    return result;
}

