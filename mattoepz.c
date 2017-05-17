#include "matrix.h"
#include <stdlib.h>

/** \brief Computes the symmetric Toeplitz matrix from a co-efficient matrix
 *
 * \param[in] R Input coefficient matrix
 * \param[in] result Matrix to store the result
 * \return \f$ \textrm{symtoep}(\mathbf{R}) \f$
 *
 */

MATRIX mat_symtoeplz(MATRIX R, MATRIX result)
{
    int i, j, n;
    n = MatRow(R);
    if(result==NULL) if((result = mat_creat(n, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<n; ++j)
        {
            result[i][j] = R[abs(i-j)][0];
        }
    }
    return result;
}

