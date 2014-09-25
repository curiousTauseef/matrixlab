#include "matrix.h"


MATRIX mat_mds(MATRIX d, int dims, int type, MATRIX result)
{
    switch(type)
    {
    case MAT_MDS_METRIC:
        result = __mat_mds_metric(d, dims, result);
        break;
    case MAT_MDS_NONMETRIC:
        result = __mat_mds_nonmetric(d, dims, result);
        break;
    default:
        gen_error(GEN_BAD_TYPE);
    }
    return result;
}

mtype __mat_mulfunc(mtype x, mtype y)
{
    return (x*y);
}

MATRIX __mat_mds_metric(MATRIX d, int dims, MATRIX result)
{
    MATRIX P = NULL, B = NULL, J = NULL, tmp = NULL;
    INT_VECTOR tmp2 = NULL;
    MATSTACK E = NULL;
    if(MatRow(d)!=MatCol(d)) mat_error(MAT_SIZEMISMATCH);
    if(result== NULL) if((result = mat_creat(MatRow(d), dims, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    P = mat_gfunc(d, __mat_sqrfunc, P);
    J = mat_creat(MatRow(d), MatCol(d), UNIT_MATRIX);
    J = mat_subs(J, 1/((mtype)MatCol(d)), J);
    tmp = mat_mul(J, P, tmp);
    B = mat_mul(tmp, J, B);
    mat_free(tmp);
    mat_free(J);
    B = mat_muls(B, -0.5, B);

    E = mat_eig_sym(B, E);
    tmp2 = int_vec_creat(dims, SERIES_INT_VECTOR);
    tmp2 = int_vec_subs(tmp2, MatRow(E[0])-1, tmp2);
    tmp2 = int_vec_muls(tmp2, -1, tmp2);

    J = mat_get_sub_matrix_from_cols(E[1], tmp2, NULL);
    tmp = mat_get_sub_matrix_from_rows(E[0], tmp2, NULL);
    tmp = mat_gfunc(tmp, __mat_sqrtfunc, tmp);

    B = mat_tran(tmp, NULL);
    result = mat_bsxfun(J, B, __mat_mulfunc, result);
    mat_free(J);
    mat_free(tmp);
    mat_free(B);
    int_vec_free(tmp2);
    matstack_free(E);
    return result;
}

MATRIX __mat_mds_nonmetric(MATRIX d, int dims, MATRIX result)
{
    if(MatRow(d)!=MatCol(d)) mat_error(MAT_SIZEMISMATCH);
    if(result== NULL) if((result = mat_creat(MatRow(d), dims, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    return result;
}


