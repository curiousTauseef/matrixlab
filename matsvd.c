#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )
#define MAX(a,b) ((a)>(b) ? (a) : (b))

static mtype PYTHAG(mtype a, mtype b)
{
    mtype at = fabs(a), bt = fabs(b), ct, result;

    if(at>bt){ct = bt / at; result = at*sqrt(1.0+ct*ct);}
    else if(bt>0.0){ct = at/bt; result = bt*sqrt(1.0+ct*ct);}
    else result = 0.0;
    return(result);
}

MATSTACK mat_svd(mtype **a, int m, int n, mtype *w, mtype **v)
{
    int flag, i, its, j, jj, k, l, nm;
    mtype c, f, h, s, x, y, z;
    mtype anorm = 0.0, g = 0.0, scale = 0.0;
    mtype *rv1;
    MATSTACK ms = matstack_creat(3);
    ms[0] = mat_creat(MatRow(a),MatCol(a), UNDEFINED);
    ms[1] = mat_creat(1, MatCol(a), UNDEFINED);
    ms[2] = mat_creat(MatRow(a),MatCol(a), UNDEFINED);


    if (m < n)
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }

    rv1 = (mtype *)malloc((unsigned int) n*sizeof(mtype));

/* Householder reduction to bidiagonal form */
    for(i=0; i<n; ++i)
    {
        /* left-hand reduction */
        l = i+1;
        rv1[i] = scale*g;
        g = s = scale = 0.0;
        if(i<m)
        {
            for(k=i; k<m; ++k)
                scale += fabs((mtype)a[k][i]);
            if(scale)
            {
                for(k=i; k<m; ++k)
                {
                    a[k][i] = (mtype)((mtype)a[k][i]/scale);
                    s += ((mtype)a[k][i]*(mtype)a[k][i]);
                }
                f = (mtype)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f*g-s;
                a[i][i] = (mtype)(f-g);
                if(i!=n-1)
                {
                    for(j=l; j<n; ++j)
                    {
                        for(s=0.0, k=i; k<m; ++k)
                            s += ((mtype)a[k][i]*(mtype)a[k][j]);
                        f = s/h;
                        for(k=i; k<m; ++k)
                            a[k][j] += (mtype)(f*(mtype)a[k][i]);
                    }
                }
                for(k=i; k<m; ++k)
                    a[k][i] = (mtype)((mtype)a[k][i]*scale);
            }
        }
        w[i] = (mtype)(scale*g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if(i<m && i!=n-1)
        {
            for(k=l; k<n; ++k)
                scale += fabs((mtype)a[i][k]);
            if(scale)
            {
                for(k=l; k<n; ++k)
                {
                    a[i][k] = (mtype)((mtype)a[i][k]/scale);
                    s += ((mtype)a[i][k]*(mtype)a[i][k]);
                }
                f = (mtype)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f*g-s;
                a[i][l] = (mtype)(f-g);
                for(k=l; k<n; ++k)
                    rv1[k] = (mtype)a[i][k]/h;
                if(i!= m-1)
                {
                    for(j=l; j<m; ++j)
                    {
                        for(s=0.0, k=l; k<n; ++k)
                            s += ((mtype)a[j][k]*(mtype)a[i][k]);
                        for(k=l; k<n; ++k)
                            a[j][k] += (mtype)(s*rv1[k]);
                    }
                }
                for(k=l; k<n; ++k)
                    a[i][k] = (mtype)((mtype)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((mtype)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for(i=n-1; i>=0; --i)
    {
        if(i<n-1)
        {
            if(g)
            {
                for(j=l; j<n; ++j)
                    v[j][i] = (mtype)(((mtype)a[i][j] / (mtype)a[i][l]) / g);
                    /* mtype division to avoid underflow */
                for(j=l; j<n; ++j)
                {
                    for(s=0.0, k=l; k<n; ++k)
                        s += ((mtype)a[i][k]*(mtype)v[k][j]);
                    for(k=l; k<n; ++k)
                        v[k][j] += (mtype)(s*(mtype)v[k][i]);
                }
            }
            for(j=l; j<n; ++j)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for(i=n-1; i>=0; --i)
    {
        l = i+1;
        g = (mtype)w[i];
        if(i<n-1)
            for(j=l; j<n; ++j)
                a[i][j] = 0.0;
        if(g)
        {
            g = 1.0/g;
            if(i != n-1)
            {
                for(j=l; j<n; ++j)
                {
                    for(s=0.0, k=l; k<m; ++k)
                        s += ((mtype)a[k][i]*(mtype)a[k][j]);
                    f = (s/(mtype)a[i][i])*g;
                    for(k=i; k<m; ++k)
                        a[k][j] += (mtype)(f*(mtype)a[k][i]);
                }
            }
            for(j=i; j<m; ++j)
                a[j][i] = (mtype)((mtype)a[j][i]*g);
        }
        else
        {
            for(j=i; j<m; ++j)
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for(k=n-1; k>=0; --k)
    {                             /* loop over singular values */
        for(its=0; its<30; ++its)
        {                         /* loop over allowed iterations */
            flag = 1;
            for(l=k; l>=0; --l)
            {                     /* test for splitting */
                nm = l-1;
                if(fabs(rv1[l])+anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if(fabs((mtype)w[nm])+anorm == anorm)
                    break;
            }
            if(flag)
            {
                c = 0.0;
                s = 1.0;
                for(i=l; i<=k; ++i)
                {
                    f = s*rv1[i];
                    if(fabs(f)+anorm != anorm)
                    {
                        g = (mtype)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (mtype)h;
                        h = 1.0/h;
                        c = g*h;
                        s = (-f*h);
                        for(j=0; j<m; ++j)
                        {
                            y = (mtype)a[j][nm];
                            z = (mtype)a[j][i];
                            a[j][nm] = (mtype)(y*c+z*s);
                            a[j][i] = (mtype)(z*c-y*s);
                        }
                    }
                }
            }
            z = (mtype)w[k];
            if(l==k)
            {                  /* convergence */
                if(z<0.0)
                {              /* make singular value nonnegative */
                    w[k] = (mtype)(-z);
                    for(j=0; j<n; ++j)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if(its>=30)
            {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = (mtype)w[l];
            nm = k-1;
            y = (mtype)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = PYTHAG(f, 1.0);
            f = ((x-z)*(x+z)+h*((y/(f+SIGN(g, f)))-h))/x;

            /* next QR transformation */
            c = s = 1.0;
            for(j=l; j<=nm; ++j)
            {
                i = j+1;
                g = rv1[i];
                y = (mtype)w[i];
                h = s*g;
                g = c*g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c+g*s;
                g = g*c-x*s;
                h = y*s;
                y = y*c;
                for(jj=0; jj<n; ++jj)
                {
                    x = (mtype)v[jj][j];
                    z = (mtype)v[jj][i];
                    v[jj][j] = (mtype)(x*c+z*s);
                    v[jj][i] = (mtype)(z*c-x*s);
                }
                z = PYTHAG(f, h);
                w[j] = (mtype)z;
                if(z)
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = (c*g)+(s*y);
                x = (c*y)-(s*g);
                for(jj=0; jj<m; ++jj)
                {
                    y = (mtype)a[jj][j];
                    z = (mtype)a[jj][i];
                    a[jj][j] = (mtype)(y*c+z*s);
                    a[jj][i] = (mtype)(z*c-y*s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (mtype)x;
        }
    }
    free(rv1);
    return(NULL);
}


