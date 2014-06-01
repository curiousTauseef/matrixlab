#include <float.h>
#include "matrix.h"


MATSTACK mat_omp(MATRIX A, MATRIX b, int k, mtype tol, MATSTACK result)
{
    int i, j, m, n, nc, cs, l;
    mtype max_val;
    MATRIX At = NULL, As = NULL, r = NULL, x = NULL, xs = NULL, ip = NULL, an = NULL;
    INT_VECTOR supp = NULL, sv = NULL;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL)
    {
        if((result = matstack_creat(2))==NULL) matstack_error(MATSTACK_MALLOC);
        result[0] = mat_creat(m, 1, ZERO_MATRIX);
        result[1] = mat_creat(1, 1, UNDEFINED);
    }
    x = result[0];
    if(k>n) k = n;
    At = mat_tran(A, NULL);
    an = mat_creat(1, m, UNDEFINED);
    supp = int_vec_creat(m, ZERO_INT_VECTOR);
    sv = int_vec_creat(0, UNDEFINED);

    for(j=0; j<m; ++j)
    {
        nc = eps;
        for(i=0; i<n; ++i) nc += A[i][j]*A[i][j];
        an[0][j] = 1.0/nc;
    }
    r = mat_copy(b, NULL);
    for(l=0; l<k; ++l)
    {
        ip = mat_mul(At, r, ip);
        ip = mat_mul_dot(ip, an, ip);
        max_val = 0;
        cs = -1;
        for(i=0; i<m; ++i)
        {
            if(supp[i]==0 && fabs(ip[i][0])>max_val)
            {
                max_val = fabs(ip[i][0]);
                cs = i;
            }
        }
        supp[cs] = 1;
        sv = int_vec_append(sv, cs);
        As = mat_get_sub_matrix_from_cols(A, sv, NULL);
        xs = mat_least_squares(As, b, NULL);
        for(i=0; i<Int_VecLen(sv); ++i) x[sv[i]][0] = xs[i][0];
        r = mat_mul(As, xs, r);
        r = mat_sub(b, r, r);
        result[1][0][0] = mat_norm_p(r,2);
        mat_free(As);
        mat_free(xs);
        if(result[1][0][0]<tol) break;
    }
    mat_free(At);
    mat_free(an);
    mat_free(r);
    mat_free(ip);
    int_vec_free(supp);
    int_vec_free(sv);
    return result;
}


