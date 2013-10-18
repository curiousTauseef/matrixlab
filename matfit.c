#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

#define MAX_ITERS_RB 20


MATRIX mat_linear_ls_fit(MATRIX data, MATRIX Y, int degree)
{
    int i, j, n;
    MATRIX A, X;
    n = MatRow(data);
    A = mat_creat(n, degree + 1, ONES_MATRIX);
    for (i = 0; i<n; ++i)
    {
        for (j = degree-1; j>=0; --j) A[i][j] = data[i][0]*A[i][j+1];
    }
    X = mat_least_squares(A, Y, NULL);
    mat_free(A);
    return X;
}

MATRIX mat_least_squares(MATRIX A, MATRIX Y, MATRIX result)
{
    int m, n, o, i, j, k;
    MATRIX Apinv;
    if(MatRow(A)!=MatRow(Y)) return mat_error(MAT_SIZEMISMATCH);
    Apinv = mat_pinv(A, NULL);
    if (Apinv ==NULL) return mat_error(MAT_INVERSE_ILL_COND);
    m = MatCol(Apinv);
    n = MatRow(Apinv);
    o = MatCol(Y);
    if (result == NULL) if ((result= mat_creat( n, o, UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);
    for (i=0; i<n; ++i)
        for (j=0; j<o; ++j)
            for (k=0, result[i][j]=0.0; k<m; ++k)
            {
                result[i][j] += Apinv[i][k] * Y[k][j];
            }

    mat_free(Apinv);
    return result;
}

MATRIX mat_w_least_squares(MATRIX A, MATRIX Y, MATRIX w, MATRIX result)
{
    int m, n, o, i, j, k;
    MATRIX Awpinv, W;
    if(MatRow(w)!= MatRow(Y)) return  mat_error(MAT_SIZEMISMATCH);
    W = mat_creat_diag(w);
    if(MatRow(A)!= MatRow(Y)) return mat_error(MAT_SIZEMISMATCH);
    Awpinv = mat_wpinv(A, W, NULL);
    if (Awpinv == NULL) return mat_error(MAT_INVERSE_ILL_COND);

    m = MatCol(Awpinv);
    n = MatRow(Awpinv);
    o = MatCol(Y);
    if (result == NULL) if ((result= mat_creat( n, o, UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);
    for (i=0; i<n; ++i)
        for (j=0; j<o; ++j)
            for (k=0, result[i][j]=0.0; k<m; ++k)
            {
                result[i][j] += Awpinv[i][k] * Y[k][j];
            }
    mat_free(Awpinv);
    mat_free(W);
    return result;
}

MATRIX mat_rob_least_squares(MATRIX A, MATRIX Y, int lossfunc)
{
    int n, k;
    int flag = 0;
    mtype med = 0, madn_ = 0, norm_th= 0;
    MATRIX X = 0, res = 0, res_ = 0, W = 0, tmp1 = 0, tmp2 = 0;
    n = MatRow(A);
    W = mat_creat(n, 1, ONES_MATRIX);
    tmp1 = mat_abs(Y);
    norm_th = 0.0001f*mat_sum(tmp1);
    mat_free(tmp1);

    for (k = 0; k<MAX_ITERS_RB && flag == 0; ++k)
    {
        X = mat_w_least_squares(A, Y, W, X);
        mat_free(W);

        tmp1 = mat_mul(A, X, NULL);
        res_ = mat_sub(tmp1, Y, NULL);
        mat_free(tmp1);
        if (k == 0)
        {
            med = mat_median(res_);
            tmp1 = mat_subs(res_, med,NULL);
            tmp2 = mat_abs(tmp1);
            mat_free(tmp1);
            madn_ = mat_median(tmp2)+(float)eps;/* *6.9414 */
            mat_free(tmp2);
        }
        res = mat_abs(res_);
        if(mat_sum(res)<norm_th) flag =1;
        if(k!=(MAX_ITERS_RB -1))
        {
            switch(lossfunc)
            {
            case MAT_LOSS_HUBER:
                W = mat_huber_wt(res, madn_, madn_);
                break;
            case MAT_LOSS_BISQUARE:
                W = mat_bisquare_wt(res, madn_, madn_);
                break;
            default:
                W = mat_bisquare_wt(res, madn_, madn_);
            }
        }
        mat_free(res);
        mat_free(res_);
    }
    if(k!=MAX_ITERS_RB)mat_free(W);
    return X;
}

MATRIX mat_robust_fit(MATRIX data, MATRIX Y, int degree, int lossfunc)
{
    int i, j, n;
    MATRIX A = 0, X = 0;
    n = MatRow(data);
    A = mat_creat(n, degree+1, ONES_MATRIX);
    for (i = 0; i<n; ++i)
    {
        for (j = degree-1; j>=0; --j) A[i][j] = data[i][0]*A[i][j+1];
    }
    X = mat_rob_least_squares(A, Y, lossfunc);
    mat_free(A);
    return X;
}

