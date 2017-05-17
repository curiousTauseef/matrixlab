#include "matrix.h"

/** \brief Computes element-wise matrix inverse
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \mathbf{11}^T./\mathbf{A} \f$
 *
 */

MATRIX mat_inv_dot(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = 1.0/A[i][j];
        }
    }
    return (result);
}

/** \brief Computes element-wise matrix division
 *
 * \param[in] A First input matrix
 * \param[in] B Second input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \mathbf{A}./\mathbf{B} \f$
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
    if(o==m &&p==n)
    {
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j)
            {
                result[i][j] = A[i][j]/B[i][j];
            }
        }
    }
    else if(o==1 && p!=1)
    {
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j)
            {
                result[i][j] = A[i][j]/B[i][0];
            }
        }
    }
    else if(p==1 && o!=1)
    {
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j)
            {
                result[i][j] = A[i][j]/B[0][j];
            }
        }
    }
    else gen_error(GEN_SIZEMISMATCH);
    return result;
}

/** \brief Divides a matrix by a scalar
 *
 * \param[in] A Input matrix
 * \param[in] s Scalar
 * \param[in] result Matrix to store the result
 * \return \f$ \dfrac{\mathbf{A}}{s} \f$
 *
 */

MATRIX mat_divs(MATRIX A, mtype s, MATRIX result)
{
    return mat_muls(A, (1.0f/s), result);
}


/** \brief Divides a scalar by a matrix
 *
 * \param[in] A Input matrix
 * \param[in] s Scalar
 * \param[in] result Matrix to store the result
 * \return \f$ s\mathbf{11}^T./\mathbf{A} \f$
 *
 */

MATRIX mat_divs_inv(MATRIX A, mtype s, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = s/A[i][j];
        }
    }
    return (result);
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

INT_VECTOR int_vec_divs(INT_VECTOR A, int s, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result==NULL) if((result = int_vec_creat(m, UNDEFINED))==NULL)
            int_vec_error(INT_VEC_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m; ++i) result[i] = A[i]/s;
    return result;
}


/** \brief Computes element-wise integer vector inverse
 *
 * \param[in] A Input vector
 * \param[in] result Vector to store the result
 * \return \f$ 1./A \f$
 *
 */

INT_VECTOR int_vec_inv(INT_VECTOR A, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result==NULL) if((result = int_vec_creat(m, UNDEFINED))==NULL)
            int_vec_error(INT_VEC_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m; ++i) result[i] = 1.0/A[i];
    return result;
}



/** \brief Divides a scalar by an integer vector
 *
 * \param[in] A Input vector
 * \param[in] s Scalar
 * \param[in] result Vector to store the result
 * \return \f$ s1./A \f$
 *
 */

INT_VECTOR int_vec_divs_inv(INT_VECTOR A, int s, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(A);
    if(result==NULL) if((result = int_vec_creat(m, UNDEFINED))==NULL)
            int_vec_error(INT_VEC_MALLOC);
    #pragma omp parallel for
    for(i=0; i<m; ++i) result[i] = s/A[i];
    return result;
}


