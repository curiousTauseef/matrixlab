#include "matrix.h"


/** \brief Computes the transpose of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \mathbf{A}^T \f$
 *
 */

MATRIX mat_tran(MATRIX A, MATRIX result)
{
    int	i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result==NULL) if((result = mat_creat(m,n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<m; ++i)
    {
        for (j=0; j<n; ++j)
        {
            result[i][j] = A[j][i];
        }
    }
    return result;
}

