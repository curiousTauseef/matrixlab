#include "matrix.h"


/** \brief Computes element-wise matrix division
 *
 * \param[in] A First input matrix
 * \param[in] B Second input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ A./B \f$
 *
 */

MATRIX mat_div_dot(MATRIX A, MATRIX B, MATRIX result)
{
    int i, j, m, n, o, p;
    m = MatCol(A);
    n = MatRow(A);
    o = MatCol(B);
    p = MatRow(B);
    if(result==NULL) if((result = mat_creat(MatRow(A), MatCol(A), UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            if(o==m &&p==n) result[i][j] = A[i][j]/B[i][j];
            else if(o==1 && p!=1) result[i][j] = A[i][j]/B[i][0];
            else if(p==1 && o!=1) result[i][j] = A[i][j]/B[0][j];
            else gen_error(GEN_SIZEMISMATCH);
        }
    return result;
}

/** \brief Divides a matrix by a scalar
 *
 * \param[in] A Input matrix
 * \param[in] s Scalar
 * \param[in] result Matrix to store the result
 * \return \f$ \dfrac{A}{s} \f$
 *
 */

MATRIX mat_divs(MATRIX A, mtype s, MATRIX result)
{
    return mat_muls(A, (1.0f/s), result);
}

/** \brief Computes element-wise integer vector division
 *
 * \param[in] A First input vector
 * \param[in] B Second input vector
 * \param[in] result Vector to store the result
 * \return \f$ A./B \f$
 *
 */

INT_VECTOR int_vec_div(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result==NULL) if((result = int_vec_creat(m, UNDEFINED))==NULL)
            int_vec_error(INT_VEC_MALLOC);
    if(m!=Int_VecLen(B)) gen_error(GEN_SIZEMISMATCH);
    #pragma omp parallel for
    for(i=0; i<m; ++i) result[i] = A[i]/B[i];
    return result;
}

/** \brief Divides an integer vector by a scalar
 *
 * \param[in] A Input vector
 * \param[in] s Scalar
 * \param[in] result Vector to store the result
 * \return \f$ \dfrac{A}{s} \f$
 *
 */

INT_VECTOR int_vec_divs(INT_VECTOR A, int x, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result==NULL) if((result = int_vec_creat(m, UNDEFINED))==NULL)
            int_vec_error(INT_VEC_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m; ++i) result[i] = A[i]/x;
    return result;
}

