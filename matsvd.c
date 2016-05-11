#include "matrix.h"

MATSTACK mat_svd(MATRIX a, int niters, MATSTACK result)
{
    int i, j, k, l, m, n, iter;
    mtype rmax = 0.0, c = 0.0, f = 0.0, p = 0.0, q= 0.0, r1 = 0.0, r2 = 0.0, rtol = eps, s = 0.0;
    mtype x, y, z, scale = 0.0, tol;
    MATRIX ac = NULL, sigma = NULL, e = NULL, v = NULL;
    m = MatRow(a);
    n = MatCol(a);
    if(result==NULL)
    {
        result = matstack_creat(3);
        result[0] = mat_copy(a, NULL);
        result[1] = mat_creat(n, 1, UNDEFINED);
        result[2] = mat_creat(n, n, UNDEFINED);
    }
    else
    {
        if(MatRow(result[0])!=m ||MatCol(result[0])!=n)
        {
            mat_free(result[0]);
            result[0] = mat_creat(m, n, UNDEFINED);
        }
        result[0] = mat_copy(a, result[0]);
        if(MatRow(result[1])!=n ||MatCol(result[1])!=1)
        {
            mat_free(result[1]);
            result[1] = mat_creat(n, 1, UNDEFINED);
        }
        if(MatRow(result[2])!=n ||MatCol(result[2])!=n)
        {
            mat_free(result[2]);
            result[2] = mat_creat(n, n, UNDEFINED);
        }
    }
    ac = result[0];
    sigma = result[1];
    v = result[2];
    e = mat_creat(1, n, UNDEFINED);
    for(i=0; i<n; ++i)
    {
        e[0][i] = scale*p;
        s = 0.0, scale = 0.0;
        if(i<m)
        {
            for(j=i; j<m; ++j) scale += fabs(ac[j][i]);
            if(scale>0.0)
            {
                for(j=i; j<m; ++j)
				{
					ac[j][i] /= scale;
					s += ac[j][i]*ac[j][i];
				}
                f = ac[i][i];
                p = sqrt(s);
                if(f>=0.0) p = -p;
                q = f*p-s;
                ac[i][i] = f-p;
                if(i!=n-1)
                {
                    for(j=i+1; j<n; ++j)
                    {
                        s = 0;
                        for(k=i; k<m; ++k) s += ac[k][i]*ac[k][j];
                        f = s/q;
                        for(k=i; k<m; ++k) ac[k][j] += f*ac[k][i];
                    }
                }
                for (j=i; j<m; ++j) ac[j][i] *= scale;
            }

        }
        sigma[i][0] = scale*p;
        s = 0.0;
        p = 0.0;
        scale = 0.0;
        if(i<m && i!=n-1)
        {
            for(j=i+1; j<n; ++j) scale += fabs(ac[i][j]);
            if(scale>0.0)
            {
                for(j=i+1; j<n; ++j)
                {
                    ac[i][j] /= scale;
                    s += ac[i][j]*ac[i][j];
                }
                f = ac[i][i+1];
                p = sqrt(s);
                if(f>=0.0) p =- p;
                q = f*p-s;
                ac[i][i+1] = f-p;
                for(j=i+1; j<n; ++j) e[0][j] = ac[i][j]/q;
                if(i!=m-1)
                {
                    for(j=i+1; j<m; ++j)
                    {
                        s = 0.0;
                        for(k=i+1; k<n; ++k) s += ac[j][k]*ac[i][k];
                        for(k=i+1; k<n; ++k) ac[j][k] += s*e[0][k];
                    }
                }
                for(j=i+1; j<n; ++j) ac[i][j] *= scale;
            }
        }
        r1 = fabs(sigma[i][0])+fabs(e[0][i]);
        if(r1>rmax) rmax = r1;
    }

    for(i=n-1; i>=0; --i)
    {
        if(i<n-1)
        {
            if(fabs(p)>0.0)
            {
                q = ac[i][i+1]*p;
                for(j=i+1; j<n; ++j) v[j][i] = ac[i][j]/q;
                for(j=i+1; j<n; ++j)
                {
                    s = 0.0;
                    for(k=i+1; k<n; ++k) s += ac[i][k]*v[k][j];
                    for(k=i+1; k<n; ++k) v[k][j] += s*v[k][i];
                }
            }
            for(j=i+1; j<n; ++j)
            {
                v[i][j] = 0.0;
                v[j][i] = 0.0;
            }
        }
        v[i][i] = 1.0;
        p = e[0][i];
    }
    for(i=n-1; i>=0; --i)
    {
        p = sigma[i][0];
        if(i<n-1) for(j=i+1; j<n; ++j) ac[i][j] = 0.0;
        if(fabs(p)>0.0)
        {
            q = ac[i][i]*p;
            if(i<n-1)
            {
                for(j=i+1; j<n; ++j)
                {
                    s = 0.0;
                    for(k=i+1; k<m; ++k) s += ac[k][i]*ac[k][j];
                    f = s/q;
                    for(k=i; k<m; ++k) ac[k][j] += f*ac[k][i];
                }
            }
            for(j=i; j<m; ++j) ac[j][i] /= p;
        }
        else for(j=i; j<m; ++j) ac[j][i] = 0.0;
        ++ac[i][i];
    }

    tol = rtol*rmax;
    for(k=n-1; k>=0; --k)
    {
        for(iter=0; iter<niters; ++iter)
        {
            char flg = 0;
            for(l=k; l>=0; --l)
            {
                if(fabs(e[0][l])<=tol)
                {
                    flg = 1;
                    break;
                }
                if(fabs(sigma[l-1][0])<=tol) break;
            }
            if(flg==0)
            {
                c = 0.0;
                s = 1.0;
                for(i=l; i<=k; ++i)
                {
                    f = s*e[0][i];
                    e[0][i] = c*e[0][i];
                    if(fabs(f)<=tol) break;
                    p = sigma[i][0];
                    sigma[i][0] = sqrt(f*f+p*p);
                    c = p/sigma[i][0];
                    s = -f/sigma[i][0];
                    for(j=0; j<m; ++j)
                    {
                        r1 = ac[j][l-1];
                        r2 = ac[j][i];
                        ac[j][l-1] = r1*c+r2*s;
                        ac[j][i] = c*r2-s*r1;
                    }
                }
            }
            z = sigma[k][0];
            if(l==k)
            {
                if(z<0.0)
                {
                    sigma[k][0] = -z;
                    for(j=0; j<n; ++j) v[j][k] = -v[j][k];
                }
                break;
            }
            x = sigma[l][0];
            y = sigma[k-1][0];
            p = e[0][k-1];
            q = e[0][k];
            f = ((y-z)*(y+z)+(p-q)*(p+q))/(2.0*q*y);
            p = sqrt(1.0+f*f);
            if(f<0.0) p=-p;
            f = ((x-z)*(x+z)+q*(y/(f+p)-q))/x;
            c = 1.0;
            s = 1.0;
            for(i=l+1; i<=k; ++i)
            {
                p = e[0][i];
                y = sigma[i][0];
                q = s*p;
                p = c*p;
                e[0][i-1] = sqrt(f*f+q*q);
                c = f/e[0][i-1];
                s = q/e[0][i-1];
                f = c*x+s*p;
                p = c*p-s*x;
                q = s*y;
                y = c*y;
                for(j=0; j<n; ++j)
                {
                    x = v[j][i-1];
                    z = v[j][i];
                    v[j][i-1] = c*x+s*z;
                    v[j][i] = c*z-s*x;
                }
                sigma[i-1][0] = sqrt(f*f+q*q);
                if(fabs(sigma[i-1][0])>0.0)
                {
                    c = f/sigma[i-1][0];
                    s = q/sigma[i-1][0];
                }
                f = c*p+s*y;
                x = c*y-s*p;
                for(j=0; j<m; ++j)
                {
                    y = ac[j][i-1];
                    z = ac[j][i];
                    ac[j][i-1] = c*y+s*z;
                    ac[j][i] = c*z-s*y;
                }
            }
            e[0][l] = 0.0;
            e[0][k] = f;
            sigma[k][0] = x;
        }
    }
    mat_free(e);
    return result;
}
