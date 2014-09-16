#include <string.h>
#include "matrix.h"
#define MAX_ITERS_RB 20


/** \brief Polynomial model using least squares
 *
 * \param[in] A Data matrix \f$ N \times 1 \f$
 * \param[in] Y Observation matrix \f$ N \times 1 \f$
 * \return
 *
 */

MATRIX mat_linear_ls_fit(MATRIX A, MATRIX Y, int deg, MATRIX result)
{
    int i, j, n;
    MATRIX B;
    n = MatRow(A);
    B = mat_creat(n, deg+1, ONES_MATRIX);
    for(i=0; i<n; ++i)
    {
        for(j=deg-1; j>=0; --j) B[i][j] = A[i][0]*B[i][j+1];
    }
    result = mat_least_squares(B, Y, result);
    mat_free(B);
    return result;
}

MATRIX mat_least_squares(MATRIX A, MATRIX Y, MATRIX result)
{
    int m, n, o, i, j, k;
    MATRIX Apinv;
    if(MatRow(A)!=MatRow(Y)) return mat_error(MAT_SIZEMISMATCH);
    Apinv = mat_pinv(A, NULL);
    if(Apinv ==NULL) return mat_error(MAT_INVERSE_ILL_COND);
    m = MatCol(Apinv);
    n = MatRow(Apinv);
    o = MatCol(Y);
    if(result==NULL) if((result= mat_creat(n, o, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i)
        for(j=0; j<o; ++j)
            for(k=0, result[i][j]=0.0; k<m; ++k)
            {
                result[i][j] += Apinv[i][k]*Y[k][j];
            }

    mat_free(Apinv);
    return result;
}

MATRIX mat_w_least_squares(MATRIX A, MATRIX Y, MATRIX w, MATRIX result)
{
    int m, n, o, i, j, k;
    MATRIX Awpinv, W;
    if(MatRow(w)!= MatRow(Y)) return  mat_error(MAT_SIZEMISMATCH);
    W = mat_creat_diag(w, NULL);
    if(MatRow(A)!= MatRow(Y)) return mat_error(MAT_SIZEMISMATCH);
    Awpinv = mat_wpinv(A, W, NULL);
    if(Awpinv==NULL) return mat_error(MAT_INVERSE_ILL_COND);

    m = MatCol(Awpinv);
    n = MatRow(Awpinv);
    o = MatCol(Y);
    if(result==NULL) if((result = mat_creat(n, o, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i)
        for(j=0; j<o; ++j)
            for(k=0, result[i][j]=0.0; k<m; ++k)
            {
                result[i][j] += Awpinv[i][k]*Y[k][j];
            }
    mat_free(Awpinv);
    mat_free(W);
    return result;
}

MATRIX mat_rob_least_squares(MATRIX A, MATRIX Y, int lossfunc, MATRIX result)
{
    int n, k;
    int flag = 0;
    mtype med = 0, madn_ = 0, norm_th= 0;
    MATRIX res = NULL, res_ = NULL, W = NULL, tmp1 = NULL, tmp2 = NULL;
    n = MatRow(A);
    W = mat_creat(n, 1, ONES_MATRIX);
    tmp1 = mat_abs(Y, NULL);
    norm_th = 0.0001f*mat_sum(tmp1);
    mat_free(tmp1);
    tmp1 = NULL;
    for(k=0; k<MAX_ITERS_RB && flag==0; ++k)
    {
        result = mat_w_least_squares(A, Y, W, result);
        tmp1 = mat_mul(A, result, tmp1);
        res_ = mat_sub(tmp1, Y, res_);
        if(k==0)
        {
            med = mat_median(res_);
            tmp1 = mat_subs(res_, med, tmp1);
            tmp2 = mat_abs(tmp1, tmp2);
            madn_ = mat_median(tmp2)* 1.4826+(mtype)eps;/* *6.9414 */
            mat_free(tmp2);
        }
        res = mat_abs(res_, res);
        if(mat_sum(res)<norm_th) flag = 1;
        if(k!=(MAX_ITERS_RB-1))
        {
            switch(lossfunc)
            {
            case MAT_LOSS_HUBER:
                W = mat_huber_wt(res, 1.0, madn_*1.345, W);
                break;
            case MAT_LOSS_BISQUARE:
                W = mat_bisquare_wt(res, 1.0, madn_*4.685, W);
                break;
            default:
                W = mat_bisquare_wt(res, 1.0, madn_*4.685, W);
            }
        }
    }
    mat_free(W);
    mat_free(res);
    mat_free(res_);
    mat_free(tmp1);
    return result;
}

MATRIX mat_robust_fit(MATRIX A, MATRIX Y, int deg, int lossfunc, MATRIX result)
{
    int i, j, n;
    MATRIX B = NULL;
    n = MatRow(A);
    B = mat_creat(n, deg+1, ONES_MATRIX);
    for(i=0; i<n; ++i)
    {
        for(j=deg-1; j>=0; --j) B[i][j] = A[i][0]*B[i][j+1];
    }
    result = mat_rob_least_squares(B, Y, lossfunc, result);
    mat_free(B);
    return result;
}

