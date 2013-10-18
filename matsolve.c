#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "matrix.h"



int mat_lu(MATRIX A, MATRIX P)
{
    int i, j, k, n;
    int maxi;
    mtype c, c1, tmp;
    int p;

    n = MatCol(A);
    for(p=0,i=0; i<n; ++i)
    {
        P[i][0] = (mtype)i;
    }

    for(k=0; k<n; ++k)
    {
        for(i=k, maxi=k, c=0.0; i<n; ++i)
        {
            c1 = (mtype)fabs( A[(int)P[i][0]][k] );
            if(c1>c)
            {
                c = c1;
                maxi = i;
            }
        }
        if(k != maxi)
        {
            ++p;
            tmp = P[k][0];
            P[k][0] = P[maxi][0];
            P[maxi][0] = tmp;
        }

        if( A[(int)P[k][0]][k] == 0.0 )
            return (-1);

        for(i=k+1; i<n; ++i)
        {
            A[(int)P[i][0]][k] = A[(int)P[i][0]][k] / A[(int)P[k][0]][k];
            for(j=k+1; j<n; ++j)
            {
                A[(int)P[i][0]][j] -= A[(int)P[i][0]][k] * A[(int)P[k][0]][j];
            }
        }
    }
    return (p);
}

void mat_backsubs1(MATRIX A, MATRIX B, MATRIX X, MATRIX P, int xcol)
{
    int i, j, k, n;
    mtype sum;

    n = MatCol(A);

    for(k=0; k<n; ++k)
    {
        for(i=k+1; i<n; ++i)
            B[(int)P[i][0]][0] -= A[(int)P[i][0]][k] * B[(int)P[k][0]][0];
    }

    X[n-1][xcol] = B[(int)P[n-1][0]][0] / A[(int)P[n-1][0]][n-1];
    for(k=n-2; k>=0; --k)
    {
        sum = 0.0;
        for(j=k+1; j<n; ++j)
        {
            sum += A[(int)P[k][0]][j] * X[j][xcol];
        }
        X[k][xcol] = (B[(int)P[k][0]][0]-sum)/A[(int)P[k][0]][k];
    }
}

MATRIX mat_lsolve(MATRIX a, MATRIX b)
{
    MATRIX A, B, X, P;
    int n;

    n = MatCol(a);
    A = mat_copy(a, NULL);
    B = mat_copy(b, NULL);
    X = mat_creat(n, 1, ZERO_MATRIX);
    P = mat_creat(n, 1, UNDEFINED);

    mat_lu(A, P);
    mat_backsubs1(A, B, X, P, 0);

    mat_free(A);
    mat_free(B);
    mat_free(P);
    return (X);
}

MATRIX mat_cholesky(MATRIX a, MATRIX result)
{
    int i, k, j, m, n;
    mtype r = 0.0, epsnorm, tmp0;
    m = MatRow(a);
    n = MatCol(a);
    if(m!=n) gen_error(GEN_SIZEMISMATCH);
    if(result== NULL)if ((result = mat_copy(a, NULL)) == NULL)
            return mat_error(MAT_MALLOC);

    for(k=0; k<n; ++k) if(a[k][k]>r) r = a[k][k];
    epsnorm = (mtype)eps*r;

    for(k=0; k<n; ++k)
    {
        tmp0 = 0.0;
        for(j=0; j<k; ++j) tmp0 += result[j][k]*result[j][k];
        for(j=k+1; j<n; ++j) result[j][k] = 0.0;
        r = result[k][k]-tmp0;
        if(r<=epsnorm)
        {
            mat_free(result);
            return NULL;
        };

        result[k][k] = r = (float)sqrt(r);

        for(j=k+1; j<n; ++j)
        {
            tmp0 = 0.0;
            for(i=0; i<k; ++i)tmp0 += result[i][j]*result[i][k];
            result[k][j] = (result[k][j] - tmp0)/r;
        }
    }
    return result;
}


