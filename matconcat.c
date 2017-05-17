#include "matrix.h"


/** \brief Concatenates two matrices
 *
 * \param[in] A Input first matrix
 * \param[in] B Input second matrix
 * \param[in] dim Concatenation direction (ROWS/COLS)
 * \return \f$ \left[\begin{array}{cc} A & B \end{array}\right] \f$ or \f$ \left[\begin{array}{c} A \\ B \end{array}\right] \f$
 *
 */

MATRIX mat_concat(MATRIX A, MATRIX B, int dim)
{
    int i, j, m, n, o, p;
    MATRIX result;
    if(A==NULL)
    {
        return mat_copy(B, NULL);
    }
    else
    {
        m = MatCol(A);
        n = MatRow(A);
    }
    if(B==NULL)
    {
        return mat_copy(A, NULL);
    }
    else
    {
        o = MatCol(B);
        p = MatRow(B);
    }
    if((dim==ROWS)&&((m==o) ||!((m==0)&&(o==0))))
    {
        if((result = mat_creat(n+p, m, UNDEFINED))==NULL) return NULL;
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j)
            {
                result[i][j] = A[i][j];
            }
        }
        #pragma omp parallel for private(j)
        for(i=0; i<p; ++i)
        {
            for(j=0; j<m; ++j)
            {
                result[i+n][j] = B[i][j];
            }
        }
        return result;
    }
    if((dim==COLS)&&((n==p) ||!((n==0)&&(p==0))))
    {
        if((result = mat_creat(n, m+o, UNDEFINED))==NULL) return NULL;
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j) result[i][j] = A[i][j];
            for(j=0; j<o; ++j) result[i][j+m] = B[i][j];
        }
        return result;
    }
    return mat_error(MAT_SIZEMISMATCH);
}

/** \brief Concatenates two integer vectors
 *
 * \param[in] a Input first vector
 * \param[in] b Input second vector
 * \param[in] result Vector to store the result
 * \return \f$ \left[\begin{array}{cc} a & b \end{array}\right] \f$ or \f$ \left[\begin{array}{c} a \\ b \end{array}\right] \f$
 *
 */

INT_VECTOR int_vec_concat(INT_VECTOR a, INT_VECTOR b, INT_VECTOR result)
{
    int i, m, n;
    m = Int_VecLen(a);
    n = Int_VecLen(b);
    if(result==NULL)
    {
        if((result = int_vec_creat(m+n, UNDEFINED))==NULL) return NULL;
    }
    else
    {
        if(Int_VecLen(result)!=(m+n)) return int_vec_error(INT_VEC_SIZEMISMATCH);
    }
    #pragma omp parallel for
    for(i=0; i<m; ++i) result[i] = a[i];
    #pragma omp parallel for
    for(i=0; i<n; ++i) result[i+m] = b[i];
    return result;
}

