#include "matrix.h"


/** \brief Dumps a matrix in the stdout
 *
 * \param[in] A Input matrix
 *
 */

void mat_dump(MATRIX A)
{
    char s[] = "%.16g ";
    mat_fdumpf(A, s, stdout);
}

/** \brief Dumps a matrix using a given format specifier in the stdout
 *
 * \param[in] A Input matrix
 * \param[in] s Format specifier
 *
 */

void mat_dumpf(MATRIX A, const char *s)
{
    mat_fdumpf(A, s, stdout);
}

/** \brief Dumps a matrix in an opened file
 *
 * \param[in] A Input matrix
 * \param[in] fp Pointer to an opened file
 *
 */

void mat_fdump(MATRIX A, MAT_FILEPOINTER fp)
{
    char s[] = "%.16g ";
    mat_fdumpf(A, s, fp);
}

/** \brief Dumps a matrix using a given format specifier in an opened file
 *
 * \param[in] A Input matrix
 * \param[in] s Format specifier
 * \param[in] fp Pointer to an opened file
 *
 */

void mat_fdumpf(MATRIX A, const char *s, MAT_FILEPOINTER fp)
{
    int i, j, m, n;
    if(A==NULL) gen_error(GEN_NOT_FOUND);
    m = MatCol(A);
    n = MatRow(A);
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            fprintf(fp, s, A[i][j]);
        }
        fprintf(fp, "\n");
    }
    fflush(fp);
}

/** \brief Dumps an integer vector in the stdout
 *
 * \param[in] A Input vector
 *
 */

void int_vec_dump(INT_VECTOR A)
{
    char s[] = "%d ";
    int_vec_fdumpf(A, s, stdout);
}

/** \brief Dumps an integer vector using a given format specifier in the stdout
 *
 * \param[in] A Input vector
 * \param[in] s Format specifier
 *
 */

void int_vec_dumpf(INT_VECTOR A, const char *s)
{
    int_vec_fdumpf(A, s, stdout);
}

/** \brief Dumps an integer vector in an opened file
 *
 * \param[in] A Input vector
 * \param[in] fp Pointer to an opened file
 *
 */

void int_vec_fdump(INT_VECTOR A, MAT_FILEPOINTER fp)
{
    char s[] = "%d ";
    int_vec_fdumpf(A, s, fp);
}

/** \brief Dumps an integer vector using a given format specifier in an opened file
 *
 * \param[in] A Input vector
 * \param[in] s Format specifier
 * \param[in] fp Pointer to an opened file
 *
 */

void int_vec_fdumpf(INT_VECTOR A, const char *s, MAT_FILEPOINTER fp)
{
    int i, n;
    if(A==NULL) gen_error(GEN_NOT_FOUND);
    n = Int_VecLen(A);
    for(i=0; i<n; ++i)
    {
        fprintf(fp, s, A[i]);
    }
    fflush(fp);
}

