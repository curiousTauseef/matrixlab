#include "matrix.h"


static mtype signa[2] = {1.0, -1.0};

/** \brief Computes a minor of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] i Row index
 * \param[in] j Column index
 * \return Minor \f$ M_{ij} \f$
 *
 */

mtype mat_minor(MATRIX A, int i, int j)
{
    MATRIX S = NULL;
    mtype result;
    S = mat_submat(A, i, j, S);
    result = mat_det(S);
    mat_free(S);
    return result;
}

/** \brief Computes a cofactor of a matrix
 *
 * \param[in] A Input matrix
 * \param[in] i Row index
 * \param[in] j Column index
 * \return Cofactor \f$ C_{ij} \f$
 *
 */

mtype mat_cofact(MATRIX A, int i, int j)
{
    mtype result;
    result = signa[(i+j)%2]*mat_minor(A, i, j);
    return result;
}

/** \brief Computes the determinant of a matrix
 *
 * \param[in] A Input matrix
 * \return \f$ \textrm{det} \left( A \right) \f$
 *
 */

mtype mat_det(MATRIX A)
{
    MATRIX b, P;
    int i, j, n;
    mtype result;
    n = MatRow(A);
    b = mat_copy(A, NULL);
    P = mat_creat(n, 1, UNDEFINED);
    i = mat_lu(b, P);
    switch(i)
    {
    case -1:
        result = 0.0;
        break;

    default:
        result = 1.0;
        for(j=0; j<n; ++j)
        {
            result *= A[(int)P[j][0]][j];
        }
        result *= signa[i%2];
        break;
    }

    mat_free(b);
    mat_free(P);
    return result;
}

