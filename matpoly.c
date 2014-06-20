#include "matrix.h"


MATSTACK mat_cheby_series_table, mat_binom_series_table;

MATRIX mat_poly_eval(MATRIX a, mtype x, int dir, MATRIX result)
{
    int m, n, i, j;
    mtype r;
    m = MatRow(a);
    n = MatCol(a);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for(j=m-1; j>=0; --j) r=r*x+a[j][i];
            result[0][i] = r;
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for(j=n-1; j>=0; --j) r=r*x+a[i][j];
            result[i][0] = r;
        }
    }
    return result;
}

MATRIX mat_poly_diff(MATRIX a, int dir, MATRIX result)
{
    int m, n, i, j;
    m = MatRow(a);
    n = MatCol(a);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(m-1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            for(j=m-2; j>=0; --j) result[j][i] = (j+1)*a[j+1][i];
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, n-1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=n-2; j>=0; --j) result[i][j]= (j+1)*a[i][j+1];
        }
    }
    return result;
}

MATRIX mat_poly_diff_eval(MATRIX a, mtype x, int dir, MATRIX result)
{
    int m, n, i, j;
    mtype r;
    m = MatRow(a);
    n = MatCol(a);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for(j=m-2; j>=0; --j) r=r*x+a[j+1][i]*(j+1);
            result[0][i] = r;
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for(j=n-2; j>=0; --j) r=r*x+a[i][j+1]*(j+1);
            result[i][0] = r;
        }
    }
    return result;
}

MATRIX mat_poly_add(MATRIX a, MATRIX b, MATRIX result)
{
    int i, na, nb;
    na = MatCol(a);
    nb = MatCol(b);
    if(na==nb) return mat_add(a, b, result);
    if(na<nb)
    {
        if(result==NULL) if((result = mat_creat(1, nb, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<na; ++i) result[0][i] = a[0][i]+b[0][i];
        for(i=na; i<nb; ++i) result[0][i] = b[0][i];
    }
    else
    {
        if(result==NULL) if((result = mat_creat(1, na, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<nb; ++i) result[0][i] = a[0][i]+b[0][i];
        for(i=nb; i<na; ++i) result[0][i] = a[0][i];
    }
    return result;
}

MATRIX mat_poly_scale(MATRIX a, mtype s, MATRIX result)
{
    return mat_muls(a, s, result);
}

MATRIX mat_poly_shift(MATRIX a, int s, MATRIX result)
{
    int i, n;
    n = MatCol(a)+s;
    if(result==NULL) if((result = mat_creat(1, (n>0)?n:0, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    mat_fill(result, 0.0);
    for(i=(s>0)?s:0; i<n; ++i) result[0][i] = a[0][i-s];
    return result;
}

void mat_cheby_init()
{
    mat_cheby_series_table = matstack_creat(512);
    mat_cheby_series_table[0] = mat_creat(1, 1, ONES_MATRIX);
    mat_cheby_series_table[1] = mat_creat(1, 2, ZERO_MATRIX);
    mat_cheby_series_table[1][0][1] = 1;
}

void mat_binom_init()
{
    mat_binom_series_table = matstack_creat(512);
    mat_binom_series_table[0] = mat_creat(1, 1, ONES_MATRIX);
    mat_binom_series_table[1] = mat_creat(1, 2, ONES_MATRIX);
}

MATRIX mat_cheby(int n)
{
    if(mat_cheby_series_table[n]==NULL)
    {
        MATRIX p = mat_cheby(n-1);
        MATRIX pp = mat_cheby(n-2);
        MATRIX npp = mat_poly_scale(pp, -1.0, NULL);
        MATRIX px2 = mat_poly_shift(p, 1, NULL);
        px2 = mat_poly_scale(px2, 2, px2);
        mat_cheby_series_table[n] = mat_poly_add(px2, npp, px2);
        mat_free(npp);
    }
    return mat_cheby_series_table[n];
}

mtype mat_binom(int n, int k)
{
    int i;
    if(mat_binom_series_table[n]==NULL)
    {
        mat_binom_series_table[n] = mat_creat(1, n+1, UNDEFINED);
        mat_binom_series_table[n][0][0] = 1;
        mat_binom_series_table[n][0][n] = 1;
        for(i=1; i<n; ++i)
        {
            mat_binom_series_table[n][0][i] = mat_binom(n-1, i-1)+mat_binom(n-1, i);
        }
    }
    return mat_binom_series_table[n][0][k];
}

MATRIX mat_cheby_coeffs_to_poly(MATRIX coeffs, MATRIX result)
{
    int i, n;
    n = MatCol(coeffs);
    if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    mat_fill(result, 0.0);
    result[0][0] = coeffs[0][0]/2;
    for(i=1; i<n; ++i)
    {
        MATRIX sp = mat_poly_scale(mat_cheby(i), coeffs[0][i], NULL);
        result = mat_poly_add(result, sp, result);
        mat_free(sp);
    }
    return result;
}

MATRIX mat_cheby_approx(mtype (*f)(mtype), mtype a, mtype b, int n, MATRIX result)
{
    int i, j;
    mtype r, c, xs, sum, m;
    MATRIX ys = NULL, coeffs = NULL, tmp = NULL, ctmp = NULL;
    r = (a-b)/2.0;
    c = -(a+b)/2.0;
    ys = mat_creat(1, n+1, UNDEFINED);
    coeffs = mat_creat(1, n+1, UNDEFINED);
    for(i=0; i<=n; ++i)
    {
        xs = r*cos(MAT_PI*(i+0.5)/(n+1));
        ys[0][i] = f(xs);
    }
    for(j=0; j<=n; ++j)
    {
        sum = 0.0;
        for(i=0; i<=n; ++i)
        {
            sum += cos(MAT_PI*j*(i+0.5)/(n+1))*ys[0][i];
        }
        coeffs[0][j] = sum*2.0/(n+1);
    }
    mat_free(ys);
    result = mat_cheby_coeffs_to_poly(coeffs, result);
    mat_free(coeffs);
    m = r;
    for(i=1; i<=n; ++i)
    {
        result[0][i] /= m;
        m *= r;
    }
    if(fabs(c)>eps)
    {
        tmp = result;
        result = mat_creat(1, n+1, ZERO_MATRIX);
        ctmp = mat_creat(1, n+1, UNDEFINED);
        ctmp[0][0] = 1;
        for(i=1; i<=n; ++i) ctmp[0][i] = ctmp[0][i-1]*c;
        mat_binom_init();
        for(j=0; j<=n; ++j)
        {
            for(i=j; i<=n; ++i)
            {
                result[0][j] += tmp[0][i]*mat_binom(i,j)*ctmp[0][i];
            }
            result[0][j] /= ctmp[0][j];
        }
        mat_free(tmp);
        mat_free(ctmp);
    }
    return result;
}

