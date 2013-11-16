#include "matrix.h"
#define SIGN(a, b) ((b) < 0 ? -fabs(a) : fabs(a) )


MATSTACK mat_pca(MATRIX data, int pca_type)
{
    int i, j, k, k2, m, n;
    MATRIX evals, interm;
    MATSTACK tmmps0 = NULL;
    MATRIX symmat, symmat2;
    m = MatCol(data);
    n = MatRow(data);

    switch(pca_type)
    {
    case MAT_PCA_CORRELATION:
        tmmps0 = mat_corcol(data);
        symmat = tmmps0[1];
        break;
    case MAT_PCA_COVARIANCE:
        tmmps0 = mat_covcol(data);
        symmat = tmmps0[1];
        break;
    case MAT_PCA_SUMOFSQUARES:
        symmat = mat_scpcol(data);
        break;
    default:
        tmmps0 = mat_covcol(data);
        symmat = tmmps0[1];
        break;
    }
    evals = mat_creat(m, 1, UNDEFINED);
    interm = mat_creat(m, 1, UNDEFINED);
    symmat2 = mat_copy(symmat, NULL);
    mat_tred2(symmat, evals, interm);
    mat_tqli(evals, interm, symmat);

    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            interm[j][0] = tmmps0[2][i][j];
        }
        for(k=0; k<3; ++k)
        {
            tmmps0[2][i][k] = 0.0;
            for(k2=0; k2<m; ++k2)
            {
                tmmps0[2][i][k] += interm[k2][0] * symmat[k2][m-k-1];
            }
        }
    }

    for(j=0; j<m; ++j)
    {
        for(k=0; k<m; ++k)
        {
            interm[k][0] = symmat2[j][k];
        }
        for(i=0; i<3; ++i)
        {
            symmat2[j][i] = 0.0;
            for (k2=0; k2<m; ++k2)
            {
                symmat2[j][i] += interm[k2][0] * symmat[k2][m-i-1];
            }
            if(evals[m-i-1][0]>0.0005)
                symmat2[j][i] /= (mtype)sqrt(evals[m-i-1][0]);
            else
                symmat2[j][i] = 0.0;
        }
    }
    mat_free(evals);
    mat_free(interm);
    return tmmps0;
}

MATSTACK mat_eig_sym(MATRIX symmat, MATSTACK result)
{
    int m, n;
    MATRIX interm, tmp_result0 = NULL, tmp_result1 = NULL;
    INT_VECTOR indcs = NULL;
    MATSTACK tmp = NULL;
    m = MatCol(symmat);
    n = MatRow(symmat);
    if(m!=n) mat_error(MAT_SIZEMISMATCH);
    if(result==NULL)
    {
        if ((result = matstack_creat(2)) == NULL)
            return matstack_error(MATSTACK_MALLOC);
        result[0] = NULL;
        result[1] = NULL;
    }
    interm = mat_creat(m, 1, UNDEFINED);
    tmp_result0 = mat_creat(m, 1, UNDEFINED);
    tmp_result1 = mat_copy(symmat, tmp_result1);
    mat_tred2(tmp_result1, tmp_result0, interm);
    mat_tqli(tmp_result0, interm, tmp_result1);

    tmp = mat_qsort(tmp_result0, ROWS, tmp);
    result[0] = mat_copy(tmp[0], result[0]);
    indcs = mat_2int_vec(tmp[1]);
    result[1] = mat_get_sub_matrix_from_cols(tmp_result1, indcs, result[1]);
    int_vec_free(indcs);
    mat_free(interm);
    mat_free(tmp_result0);
    mat_free(tmp_result1);

    return result;
}

MATSTACK mat_corcol(MATRIX data)
{
    mtype x;
    MATRIX stddev;
    int i, j, j1, j2, m, n;
    MATSTACK corr = matstack_creat(3);
    m = MatCol(data);
    n = MatRow(data);
    corr[0] = mat_creat(1, m, ZERO_MATRIX);
    stddev = mat_creat(1, m, ZERO_MATRIX);
    corr[1] = mat_creat(m, m, ZERO_MATRIX);
    corr[2] = mat_creat(n, m, ZERO_MATRIX);
    for(j=0; j<m; ++j)
    {
        corr[0][0][j] = 0.0;
        for(i=0; i<n; ++i)
        {
            corr[0][0][j] += data[i][j];
        }
        corr[0][0][j] /= (mtype)n;
    }
    for(j=0; j<m; ++j)
    {
        stddev[0][j] = 0.0;
        for(i=0; i<n; ++i)
        {
            stddev[0][j] += ((data[i][j] - corr[0][0][j])*(data[i][j] - corr[0][0][j]));
        }
        stddev[0][j] /= (mtype)n;
        stddev[0][j] = (mtype)sqrt(stddev[0][j]);
        if(stddev[0][j]<=Eps) stddev[0][j] = 1.0;
    }
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            corr[2][i][j] = data[i][j]-corr[0][0][j];
            x = (mtype)sqrt((mtype)n);
            x *= stddev[0][j];
            corr[2][i][j] = corr[2][i][j]/x;
        }
    }
    for(j1=0; j1<m-1; ++j1)
    {
        corr[1][j1][j1] = 1.0;
        for(j2=j1+1; j2<m; ++j2)
        {
            corr[1][j1][j2] = 0.0;
            for(i=0; i<n; ++i)
            {
                corr[1][j1][j2] += (corr[2][i][j1]*corr[2][i][j2]);
            }
            corr[1][j2][j1] = corr[1][j1][j2];
        }
    }
    corr[1][m-1][m-1] = 1.0;
    mat_free(stddev);
    return corr;
}

MATSTACK mat_covcol(MATRIX data)
{
    int i, j, j1, j2, m, n;
    MATSTACK covar = matstack_creat(3);
    m = MatCol(data);
    n = MatRow(data);
    covar[0] = mat_creat(1, m, ZERO_MATRIX);
    covar[1] = mat_creat(m, m, ZERO_MATRIX);
    covar[2] = mat_creat(n, m, ZERO_MATRIX);
    for(j=0; j<m; ++j)
    {
        covar[0][0][j] = 0.0;
        for(i=0; i<n; ++i)
        {
            covar[0][0][j] += data[i][j];
        }
        covar[0][0][j] /= (mtype)n;
    }

    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            covar[2][i][j] = data[i][j]-covar[0][0][j];
        }
    }
    for(j1=0; j1<m; ++j1)
    {
        for(j2=j1; j2<m; ++j2)
        {
            covar[1][j1][j2] = 0.0;
            for(i=0; i<n; ++i)
            {
                covar[1][j1][j2] += covar[2][i][j1]*covar[2][i][j2];
            }
            covar[1][j1][j2]/= (mtype)n;
            covar[1][j2][j1] = covar[1][j1][j2];
        }
    }
    return covar;
}

MATRIX mat_scpcol(MATRIX data)
{
    int i, j1, j2, m, n;
    MATRIX spco = NULL;
    m = MatCol(data);
    n = MatRow(data);
    spco = mat_creat(m, m, ZERO_MATRIX);
    for(j1=0; j1<m; ++j1)
    {
        for(j2=j1; j2<m; ++j2)
        {
            spco[j1][j2] = 0.0;
            for(i=0; i<n; ++i)
            {
                spco[j1][j2] += data[i][j1]*data[i][j2];
            }
            spco[j2][j1] = spco[j1][j2];
        }
    }
    return spco;
}

void mat_tred2(MATRIX a, MATRIX d, MATRIX e)
{
    int l, k, j, i, n;
    mtype scale, hh, h, g, f;
    n = MatRow(a);
    for(i=n-1; i>0; --i)
    {
        l = i - 1;
        h = scale = 0.0;
        if(l>0)
        {
            for(k = 0; k <= l; ++k)
                scale += (mtype)fabs(a[i][k]);
            if(scale==0.0)
                e[i][0] = a[i][l];
            else
            {
                for(k=0; k<=l; ++k)
                {
                    a[i][k] /= scale;
                    h += a[i][k]*a[i][k];
                }
                f = a[i][l];
                g = (mtype)(f>0 ? -sqrt(h) : sqrt(h));
                e[i][0] = scale*g;
                h -= f * g;
                a[i][l] = f-g;
                f = 0.0;
                for(j=0; j<=l; ++j)
                {
                    a[j][i] = a[i][j]/h;
                    g = 0.0;
                    for(k=0; k<=j; ++k)
                        g += a[j][k] * a[i][k];
                    for(k=j+1; k<=l; ++k)
                        g += a[k][j]*a[i][k];
                    e[j][0] = g/h;
                    f += e[j][0]*a[i][j];
                }
                hh = f/(h+h);
                for(j=0; j<=l; ++j)
                {
                    f = a[i][j];
                    e[j][0] = g = e[j][0]-hh*f;
                    for(k=0; k<=j; ++k)
                        a[j][k] -= (f*e[k][0]+g*a[i][k]);
                }
            }
        }
        else
            e[i][0] = a[i][l];
        d[i][0] = h;
    }
    d[0][0] = 0.0;
    e[0][0] = 0.0;
    for(i=0; i<n; ++i)
    {
        l = i-1;
        if(d[i][0])
        {
            for(j=0; j<=l; ++j)
            {
                g = 0.0;
                for(k=0; k<=l; ++k)
                    g += a[i][k]*a[k][j];
                for(k=0; k<=l; ++k)
                    a[k][j] -= g*a[k][i];
            }
        }
        d[i][0] = a[i][i];
        a[i][i] = 1.0;
        for(j=0; j<=l; ++j)
            a[j][i] = a[i][j] = 0.0;
    }
}

void mat_tqli(MATRIX d, MATRIX e, MATRIX z)
{
    int m, l, iter, i, k, n;
    mtype s, r, p, g, f, dd, c, b;
    n = MatRow(z);

    for(i=1; i<n; ++i) e[i-1][0] = e[i][0];
    e[n-1][0] = 0.0;
    for(l=0; l<n; ++l)
    {
        iter = 0;
        do
        {
            for(m=l; m<n-1; ++m)
            {
                dd = (mtype)(fabs(d[m][0]) + fabs(d[m+1][0]));
                if(fabs(e[m][0])+dd == dd) break;
            }
            if(m!=l)
            {
                if(iter++ == 50) gen_error(GEN_NOT_CONVERGED);
                g = (d[l+1][0] - d[l][0]) / (2.0f * e[l][0]);
                r = (mtype)sqrt((g * g) + 1.0f);
                g = d[m][0] - d[l][0] + e[l][0] / (mtype)(g + SIGN(r, g));
                s = c = 1.0;
                p = 0.0;
                for(i=m-1; i>=l; --i)
                {
                    f = s*e[i][0];
                    b = c*e[i][0];
                    if(fabs(f) >= fabs(g))
                    {
                        c = g/f;
                        r = (mtype)sqrt((c*c)+1.0f);
                        e[i+1][0] = f*r;
                        c *= (s = 1.0f/r);
                    }
                    else
                    {
                        s = f/g;
                        r = (mtype)sqrt((s*s)+1.0f);
                        e[i+1][0] = g*r;
                        s *= (c = 1.0f/r);
                    }
                    g = d[i+1][0]-p;
                    r = (d[i][0]-g)*s+2.0f*c*b;
                    p = s*r;
                    d[i+1][0] = g+p;
                    g = c*r-b;
                    for(k=0; k<n; ++k)
                    {
                        f = z[k][i+1];
                        z[k][i+1] = s*z[k][i]+c*f;
                        z[k][i] = c*z[k][i]-s*f;
                    }
                }
                d[l][0] = d[l][0]-p;
                e[l][0] = g;
                e[m][0] = 0.0;
            }
        }
        while(m != l);
    }
}

