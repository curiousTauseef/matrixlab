#include "matrix.h"


MATRIX mat_dump(MATRIX A)
{
    char s[] = "%g ";
    return mat_fdumpf(A, s, stdout);
}

MATRIX mat_dumpf(MATRIX A, const char *s)
{
    return mat_fdumpf(A, s, stdout);
}

MATRIX mat_fdump(MATRIX A, FILE *fp)
{
    char s[] = "%g ";
    return mat_fdumpf(A, s, fp);
}

MATRIX mat_fdumpf(MATRIX A, const char *s, FILE *fp)
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
    return A;
}

INT_VECTOR int_vec_dump(INT_VECTOR A)
{
    char s[] = "%d ";
    return int_vec_fdumpf(A, s, stdout);
}

INT_VECTOR int_vec_dumpf(INT_VECTOR A, const char *s)
{
    return int_vec_fdumpf(A, s, stdout);
}

INT_VECTOR int_vec_fdump(INT_VECTOR A, FILE *fp)
{
    char s[] = "%d ";
    return int_vec_fdumpf(A, s, fp);
}

INT_VECTOR int_vec_fdumpf(INT_VECTOR A, const char *s, FILE *fp)
{
    int i, n;
    if(A==NULL) gen_error(GEN_NOT_FOUND);
    n = Int_VecLen(A);
    for(i=0; i<n; ++i)
    {
        fprintf(fp, s, A[i]);
    }
    fflush(fp);
    return A;
}

