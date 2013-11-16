#include <string.h>
#include <stdlib.h>
#include "matrix.h"


int mats_isnan(mtype x)
{
    volatile mtype temp = x;
    return temp != x;
}

int mats_isinf(mtype x)
{
    volatile mtype temp = x;
    if ((temp == x) && ((temp - x) != 0.0))
        return (x < 0.0 ? -1 : 1);
    else return 0;
}

void mat_nextline(void)
{
    mat_fnextline(stdout);
}

void mat_fnextline(FILE *fp)
{
    fprintf(fp, "\n");
}

MATRIX mat_bsxfun(MATRIX a, MATRIX b, mtype (*pt2func)(mtype, mtype), MATRIX result)
{
    int m, n, o, p, i, j;
    m = MatRow(a);
    n = MatCol(a);

    o = MatRow(b);
    p = MatCol(b);
    if(m<o && n==p && m==1)
    {
        if(result== NULL) if((result = mat_creat(o, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<o; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*pt2func)(a[0][j], b[i][j]);
            }
        }
    }
    else if(m>o && n==p && o==1)
    {
        if(result== NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*pt2func)(a[i][j], b[0][j]);
            }
        }
    }
    else if(m==o && n<p && n==1)
    {
        if(result== NULL) if((result = mat_creat(m, p, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<p; ++j)
            {
                result[i][j] = (*pt2func)(a[i][0], b[i][j]);
            }
        }
    }
    else if(m==o && n>p && p==1)
    {
        if(result== NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*pt2func)(a[i][j], b[i][0]);
            }
        }
    }
    else if(m==1 && p==1)
    {
        if(result== NULL) if((result = mat_creat(o, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<o; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*pt2func)(a[0][j], b[i][0]);
            }
        }
    }
    else if(n==1 && o==1)
    {
        if(result== NULL) if((result = mat_creat(m, p, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<p; ++j)
            {
                result[i][j] = (*pt2func)(a[i][0], b[0][j]);
            }
        }
    }
    else mat_error(MAT_SIZEMISMATCH);
    return result;
}

INT_VECTOR int_vec_permute_vect(int n, int k, INT_VECTOR result)
{
    int i, j, d;
    if(result== NULL) if((result = int_vec_creat(k, SERIES_INT_VECTOR))==NULL)
            return int_vec_error(INT_VEC_MALLOC);
    for(i=0; i<n-1; i++)
    {
        j = result[i];
        d = i + (rand()%(n-i));
        result[i] = result[d];
        result[d] = j;
    }
    return result;
}

INT_VECTOR mat_get_sub_vector(INT_VECTOR data, INT_VECTOR indices)
{
    int i, n;
    INT_VECTOR subvec;
    n = Int_VecLen(indices);
    subvec = int_vec_creat(n, UNDEFINED);
    for(i=0; i<n; ++i)subvec[i] = data[indices[i]];
    return subvec;
}

MATRIX mat_get_sub_matrix_from_rows(MATRIX data, INT_VECTOR indices, MATRIX result)
{
    int i, j, k, n;
    k = MatCol(data);
    n = Int_VecLen(indices);
    if(result==NULL) if((result = mat_creat(n, k, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i)
    {
        for(j=0; j<k; ++j) result[i][j] = data[indices[i]][j];
    }
    return result;
}

MATRIX mat_get_sub_matrix_from_cols(MATRIX data, INT_VECTOR indices, MATRIX result)
{
    int i, j, k, n;
    k = MatRow(data);
    n = Int_VecLen(indices);
    if(result==NULL) if((result = mat_creat(k, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i)
    {
        for(j=0; j<k; ++j) result[j][i] = data[j][indices[i]];
    }
    return result;
}

MATRIX mat_calc_dist_sq(MATRIX data, MATRIX curr_data, MATRIX result)
{
    int i, j, m, n;
    mtype dist;
    n = MatRow(data);
    m = MatCol(data);
    if(result==NULL) if((result = mat_creat(n, 1, ZERO_MATRIX))==NULL) return mat_error(MAT_MALLOC);;
    for(i=0; i<n; ++i)
    {
        dist = 0.0;
        for(j=0; j<m; ++j)
        {
            dist += (data[i][j]-curr_data[i][0])*(data[i][j]-curr_data[i][0]);
        }
        result[i][0] = dist;
    }
    return result;
}

INT_VECTOR mat_find_within_dist(MATRIX data, MATRIX curr_data, mtype range)
{
    int i, count = 0, n;
    INT_VECTOR indices = 0;
    MATRIX dist = NULL;
    mtype tmp0;
    n = MatRow(data);
    tmp0 = range*range;
    dist = mat_calc_dist_sq(data, curr_data, dist);
    for(i=0; i<n; ++i)
    {
        if(dist[i][0]<tmp0)
        {
            ++count;
            indices = (int*) realloc(indices, (count+1)*sizeof(int));
            indices[count] = i;
        }
    }
    indices[0] = count;
    indices = (int *)indices+1;
    mat_free(dist);
    return indices;
}

MATRIX mat_pick_row(MATRIX data, int r, MATRIX result)
{
    int m, i;
    m = MatCol(data);
    if(result==NULL) if((result = mat_creat(1, m, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(i=0; i<m; ++i) result[0][i] = data[r][i];
    return result;
}

MATRIX mat_pick_col(MATRIX data, int c, MATRIX result)
{
    int n, i;
    n = MatRow(data);
    if(result==NULL) if((result = mat_creat(n, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i) result[i][0] = data[i][c];
    return result;
}

void __mat_cart2pol(mtype x, mtype y, mtype *rho, mtype *th)
{
    *rho = (mtype)sqrt(x*x + y*y);
    *th = (mtype)atan2(y, x);
}

#define __cart2pol(x, y, rho, th)\
{\
    rho = (mtype)sqrt(x*x + y*y);\
    th = (mtype)atan2(y, x);\
}

void __mat_pol2cart(mtype rho, mtype th, mtype *x, mtype *y)
{
    *x = (mtype)(rho*cos(th));
    *y = (mtype)(rho*sin(th));
}

#define __pol2cart(rho, th, x, y)\
{\
    x = (mtype)(rho*cos(th));\
    y = (mtype)(rho*sin(th));\
}

MATRIX mat_cart2pol(MATRIX A, int dim, MATRIX result)
{
    int m, n, i;
    m = MatCol(A);
    n = MatRow(A);
    if(dim==0 && n>1)
    {
        if(result==NULL)if((result = mat_creat(2, m, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<m; ++i)
            __cart2pol(A[0][i], A[1][i], result[0][i], result[1][i]);
        return result;
    }
    if(dim==1 && m>1)
    {
        if(result==NULL)if((result = mat_creat(n, 2, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<n; ++i)
            __cart2pol(A[i][0], A[i][1], result[i][0], result[i][1]);
        return result;
    }
    return mat_error (MAT_SIZEMISMATCH);
}

MATRIX mat_pol2cart(MATRIX A, int dim, MATRIX result)
{
    int m, n, i;
    m = MatCol(A);
    n = MatRow(A);
    if(dim==0 && n>1)
    {
        if(result==NULL)if((result = mat_creat(2, m, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<m; ++i)
            __pol2cart(A[0][i], A[1][i], result[0][i], result[1][i]);
        return result;
    }
    if(dim==1 && m>1)
    {
        if(result==NULL)if((result = mat_creat(n, 2, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<n; ++i)
            __pol2cart(A[i][0], A[i][1], result[i][0], result[i][1]);
        return result;
    }
    return mat_error (MAT_SIZEMISMATCH);
}

int __mat_powerof2(int width, int *m, int *twopm)
{
    *m = (int)ceil(log((mtype)width)/log((mtype)2.0));
    *twopm = (int)pow((mtype)2.0, (mtype)*m);

    return 1;
}

INT_VECTOR int_vec_unique(INT_VECTOR a)
{
    int i, l, uni_l = 0;
    mtype *ordered = NULL, **p_ordered;
    INT_VECTOR int_ordered = NULL;
    SEARCH_TREE s;
    p_ordered = &ordered;
    l = Int_VecLen(a);
    s = mat_bs_make_null();
    for (i = 0; i<l; ++i)
    {
        s = mat_bs_insert((mtype)a[i],s);
    }
    uni_l = mat_bs_inorder(s, uni_l, p_ordered);
    int_ordered = int_vec_creat(uni_l, UNDEFINED);
    for(i = 0; i< uni_l; ++i) int_ordered[i] = (int)ordered[i];
    free(ordered);
    s = mat_bs_free(s);
    return int_ordered;
}

